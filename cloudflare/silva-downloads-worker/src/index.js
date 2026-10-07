import { AwsClient } from "aws4fetch";
import { objectKey, publishedReleases } from "./keys.js";

// Large database objects need enough time for a slow client's initial request.
const URL_TTL_SECONDS = 3600;
const noStore = { "Cache-Control": "no-store" };

function error(message, status, extraHeaders = {}) {
  return Response.json({ error: message }, {
    status,
    headers: { ...noStore, ...extraHeaders },
  });
}

function configured(env) {
  return Boolean(env.SILVA_BUCKET && env.R2_ACCOUNT_ID && env.R2_BUCKET_NAME &&
    env.R2_ACCESS_KEY_ID && env.R2_SECRET_ACCESS_KEY);
}

function authorized(request, token) {
  return request.headers.get("Authorization") === `Bearer ${token}`;
}

export default {
  async fetch(request, env) {
    // Only allowlisted objects receive GET-only URLs.
    if (request.method !== "GET") {
      return error("Method not allowed", 405, { Allow: "GET" });
    }

    const requestUrl = new URL(request.url);
    if (requestUrl.pathname !== "/v1/releases" && requestUrl.pathname !== "/v1/download-url") {
      return error("Not found", 404);
    }
    if (!env.SILVA_APP_TOKEN) return error("Worker is not configured", 503);
    if (!authorized(request, env.SILVA_APP_TOKEN)) {
      return error("Unauthorized", 401, { "WWW-Authenticate": "Bearer" });
    }

    if (requestUrl.pathname === "/v1/releases") {
      if (!configured(env)) return error("Worker is not configured", 503);
      try {
        const releases = [];
        for (const release of publishedReleases) {
          const groups = [
            ["manifest.json"],
            ["SILVA.arb.zst", "SILVA.arb"],
            ["SILVA.sidx.zst", "SILVA.sidx"],
          ];
          const available = await Promise.all(groups.map(async (files) => {
            const keys = files.map((file) => objectKey(release, file)).filter(Boolean);
            const objects = await Promise.all(keys.map((key) => env.SILVA_BUCKET.head(key)));
            return objects.some((object) => object && object.size > 0);
          }));
          if (available.every(Boolean)) {
            releases.push(release);
          }
        }
        return Response.json({ releases }, { headers: noStore });
      } catch {
        return error("Storage lookup failed", 502);
      }
    }
    const releaseValues = requestUrl.searchParams.getAll("release");
    const fileValues = requestUrl.searchParams.getAll("file");
    if (releaseValues.length !== 1 || fileValues.length !== 1) {
      return error("Specify one release and one file", 400);
    }

    const key = objectKey(releaseValues[0], fileValues[0]);
    if (!key) return error("Release or file is not available", 404);

    if (!configured(env)) {
      return error("Worker is not configured", 503);
    }

    let object;
    try {
      object = await env.SILVA_BUCKET.head(key);
    } catch {
      return error("Storage lookup failed", 502);
    }
    if (!object) return error("Object not found in R2", 404);

    const objectUrl = new URL(`https://${env.R2_ACCOUNT_ID}.r2.cloudflarestorage.com`);
    objectUrl.pathname = `/${env.R2_BUCKET_NAME}/${key}`;
    objectUrl.searchParams.set("X-Amz-Expires", String(URL_TTL_SECONDS));

    try {
      const signer = new AwsClient({
        service: "s3",
        region: "auto",
        accessKeyId: env.R2_ACCESS_KEY_ID,
        secretAccessKey: env.R2_SECRET_ACCESS_KEY,
      });
      const signed = await signer.sign(new Request(objectUrl, { method: "GET" }), {
        aws: { signQuery: true },
      });
      return Response.json({
        release: releaseValues[0],
        file: fileValues[0],
        key,
        size_bytes: object.size,
        url: signed.url.toString(),
        expires_in_seconds: URL_TTL_SECONDS,
      }, { headers: noStore });
    } catch {
      return error("Could not sign download URL", 502);
    }
  },
};
