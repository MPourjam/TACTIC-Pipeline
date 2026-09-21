import assert from "node:assert/strict";
import test from "node:test";
import worker from "../src/index.js";

const endpoint = "https://example.workers.dev/v1/download-url";

test("refuses requests without Cloudflare Access", async () => {
  const response = await worker.fetch(
    new Request(`${endpoint}?release=144&file=SILVA.sidx`),
    {},
    {},
  );
  assert.equal(response.status, 403);
});

test("rejects object paths outside the allowlist", async () => {
  const response = await worker.fetch(
    new Request(`${endpoint}?release=144&file=..%2Fprivate.txt`),
    {},
    { access: {} },
  );
  assert.equal(response.status, 404);
});

test("does not sign an R2 object that does not exist", async () => {
  let lookedUpKey;
  const response = await worker.fetch(
    new Request(`${endpoint}?release=144&file=SILVA.sidx`),
    {
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
    { access: {} },
  );
  assert.equal(response.status, 404);
  assert.equal(
    lookedUpKey,
    "silva/144/sina-1.6.0/k10-fast/linux-amd64/SILVA.sidx",
  );
});
