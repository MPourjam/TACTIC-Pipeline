import { AwsClient } from "aws4fetch";
import { objectKey } from "./keys.js";

// Large database objects need enough time for a slow client's initial request.
const URL_TTL_SECONDS = 3600;
const noStore = { "Cache-Control": "no-store" };

function error(message, status, extraHeaders = {}) {
  return Response.json({ error: message }, {
    status,
    headers: { ...noStore, ...extraHeaders },
  });
}

export default {
  async fetch(request, env, ctx) {
    // Configure Cloudflare Access on the entire Worker (All traffic).
    // Access sets this context only after authenticating the request.
    if (!ctx.access) return error("Cloudflare Access authentication required", 403);

    if (request.method !== "GET") {
      return error("Method not allowed", 405, { Allow: "GET" });
    }

    const requestUrl = new URL(request.url);
    if (requestUrl.pathname !== "/v1/download-url") {
      return error("Not found", 404);
    }

    const releaseValues = requestUrl.searchParams.getAll("release");
    const fileValues = requestUrl.searchParams.getAll("file");
    if (releaseValues.length !== 1 || fileValues.length !== 1) {
      return error("Specify one release and one file", 400);
    }

    const key = objectKey(releaseValues[0], fileValues[0]);
    if (!key) return error("Release or file is not available", 404);

    if (!env.SILVA_BUCKET || !env.R2_ACCOUNT_ID || !env.R2_BUCKET_NAME ||
        !env.R2_ACCESS_KEY_ID || !env.R2_SECRET_ACCESS_KEY) {
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
