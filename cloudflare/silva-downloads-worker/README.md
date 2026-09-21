# SILVA download URL Worker

This Worker authorizes users with Cloudflare Access and issues short-lived,
GET-only URLs for a fixed set of objects in the private `silva-dbs` R2 bucket.
The large files are downloaded directly from R2, not through the Worker.
It does not run inside the pipeline Docker image.

## Cloudflare setup

1. Create an **account-owned** R2 API token with **Object Read only** access,
   scoped to the `silva-dbs` bucket. Keep its Access Key ID and Secret Access
   Key private. The Worker needs these credentials to sign S3 API URLs; the R2
   binding alone cannot create S3 presigned URLs.
2. In Workers & Pages, deploy this project with the build root set to
   `cloudflare/silva-downloads-worker`. The default deploy command is
   `npx wrangler deploy`. The Worker name must match `name` in `wrangler.jsonc`.
3. Add `R2_ACCESS_KEY_ID` and `R2_SECRET_ACCESS_KEY` as **Worker secrets** under
   Settings > Variables and Secrets (or use `npx wrangler secret put NAME`).
   Never put their values in this repository or the pipeline Docker image.
   The first deployment can succeed without secrets, but the endpoint returns
   503 until they are configured.
4. In the Worker's **Access** tab, select **Protect this Worker behind Access**,
   **All traffic**, and configure a policy for your users. This must protect
   production and preview URLs. The Worker refuses requests that have no
   Access context even if this dashboard setting is accidentally removed.
5. Confirm the `R2_ACCOUNT_ID`, bucket name, and uploaded object keys in
   `wrangler.jsonc` and `src/keys.js`. If the bucket uses a jurisdiction-specific
   endpoint, change the S3 hostname in `src/index.js` accordingly.

`package.json` and npm are used to build the Worker on Cloudflare or on a
developer machine. They are not needed in `Dockerfile.TACTIC`.

## API

After a user authenticates through Cloudflare Access:

```text
GET /v1/download-url?release=144&file=SILVA.sidx
```

The JSON response includes `url`, `key`, `size_bytes`, and
`expires_in_seconds`. The URL is valid for one hour and authorizes only a GET
of the requested object. Treat it as a secret until it expires; do not log it.
The `file` parameter accepts only the exact names in `src/keys.js`: the
release manifest and SILVA ARB/SIDX artifacts. SortMeRNA files are hosted
separately and are not available through this Worker. A missing R2 object
returns 404.

The laptop client should fetch the release manifest first, then request the
required artifact. It must verify the downloaded size and SHA-256 against the
manifest before atomically installing it. If a large download fails after the
URL expires, request a new URL and retry. R2 object keys and hashes in the
manifest must correspond to the uploaded objects.

## Local development

Run `npm install`, `npm test`, and `npm run dev` from this directory. For local
development only, place `R2_ACCESS_KEY_ID` and `R2_SECRET_ACCESS_KEY` in an
untracked `.dev.vars` file. Do not use real production secrets in tests. The
local R2 binding is simulated by default, so it will not see production bucket
objects unless you explicitly configure a remote binding.
