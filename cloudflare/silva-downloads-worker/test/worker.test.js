import assert from "node:assert/strict";
import test from "node:test";
import worker from "../src/index.js";

const endpoint = "https://example.workers.dev/v1/download-url";
const appToken = "test-app-token";

function authenticatedRequest(url, options = {}) {
  const headers = new Headers(options.headers);
  headers.set("Authorization", `Bearer ${appToken}`);
  return new Request(url, { ...options, headers });
}

test("issues a GET-only download URL with the application token", async () => {
  const response = await worker.fetch(
    authenticatedRequest(`${endpoint}?release=144&file=SILVA.sidx.zst`),
    {
      SILVA_APP_TOKEN: appToken,
      R2_ACCOUNT_ID: "test-account",
      R2_BUCKET_NAME: "silva-dbs",
      R2_ACCESS_KEY_ID: "test-key",
      R2_SECRET_ACCESS_KEY: "test-secret",
      SILVA_BUCKET: { async head() { return { size: 42 }; } },
    },
    {},
  );
  assert.equal(response.status, 200);
  const result = await response.json();
  const url = new URL(result.url);
  assert.equal(result.file, "SILVA.sidx.zst");
  assert.equal(result.size_bytes, 42);
  assert.equal(url.hostname, "test-account.r2.cloudflarestorage.com");
  assert.equal(url.pathname, "/silva-dbs/silva/144/sina-1.6.0/k10-fast/linux-amd64/SILVA.sidx.zst");
  assert.equal(url.searchParams.get("X-Amz-Expires"), "3600");
  assert.ok(url.searchParams.get("X-Amz-Signature"));
  assert.ok(!JSON.stringify(result).includes("test-secret"));
});

test("lists only complete, nonempty published releases with authentication", async () => {
  const response = await worker.fetch(
    authenticatedRequest("https://example.workers.dev/v1/releases"),
    {
      SILVA_APP_TOKEN: appToken,
      R2_ACCOUNT_ID: "test-account",
      R2_BUCKET_NAME: "silva-dbs",
      R2_ACCESS_KEY_ID: "test-key",
      R2_SECRET_ACCESS_KEY: "test-secret",
      SILVA_BUCKET: {
        async head(key) {
          if (key.includes("/144/") && key.includes("/SILVA.sidx")) return null;
          return { size: 42 };
        },
      },
    },
    {},
  );
  assert.equal(response.status, 200);
  assert.deepEqual(await response.json(), { releases: ["138.2"] });
  assert.equal(response.headers.get("Cache-Control"), "no-store");
});

test("reports missing server configuration for release discovery", async () => {
  const response = await worker.fetch(
    authenticatedRequest("https://example.workers.dev/v1/releases"),
    { SILVA_APP_TOKEN: appToken },
    {},
  );
  assert.equal(response.status, 503);
});

test("rejects object paths outside the allowlist", async () => {
  const response = await worker.fetch(
    authenticatedRequest(`${endpoint}?release=144&file=..%2Fprivate.txt`),
    { SILVA_APP_TOKEN: appToken },
    {},
  );
  assert.equal(response.status, 404);
});

test("does not sign an R2 object that does not exist", async () => {
  let lookedUpKey;
  const response = await worker.fetch(
    authenticatedRequest(`${endpoint}?release=144&file=SILVA.sidx.zst`),
    {
      SILVA_APP_TOKEN: appToken,
      R2_ACCOUNT_ID: "test-account",
      R2_BUCKET_NAME: "silva-dbs",
      R2_ACCESS_KEY_ID: "test-key",
      R2_SECRET_ACCESS_KEY: "test-secret",
      SILVA_BUCKET: {
        async head(key) {
          lookedUpKey = key;
          return null;
        },
      },
    },
    {},
  );
  assert.equal(response.status, 404);
  assert.equal(
    lookedUpKey,
    "silva/144/sina-1.6.0/k10-fast/linux-amd64/SILVA.sidx.zst",
  );
});

test("download endpoint rejects writes", async () => {
  const response = await worker.fetch(
    new Request(`${endpoint}?release=144&file=SILVA.sidx.zst`, { method: "PUT" }),
    {},
    {},
  );
  assert.equal(response.status, 405);
  assert.equal(response.headers.get("Allow"), "GET");
});

test("rejects missing and incorrect tokens before any R2 lookup", async () => {
  let lookups = 0;
  const env = {
    SILVA_APP_TOKEN: appToken,
    SILVA_BUCKET: { async head() { lookups++; return { size: 42 }; } },
  };
  for (const url of ["https://example.workers.dev/v1/releases",
    `${endpoint}?release=144&file=SILVA.sidx`]) {
    for (const headers of [{}, { Authorization: "Bearer wrong-token" }]) {
      const response = await worker.fetch(new Request(url, { headers }), env, {});
      assert.equal(response.status, 401);
      assert.equal(response.headers.get("WWW-Authenticate"), "Bearer");
      assert.equal(response.headers.get("Cache-Control"), "no-store");
    }
  }
  assert.equal(lookups, 0);
});

test("fails closed when the application token secret is absent", async () => {
  const response = await worker.fetch(
    authenticatedRequest("https://example.workers.dev/v1/releases"),
    {},
    {},
  );
  assert.equal(response.status, 503);
});
