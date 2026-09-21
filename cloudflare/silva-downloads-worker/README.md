# SILVA download URL Worker

This Worker provides public release discovery and short-lived, GET-only URLs
for a fixed set of SILVA objects in the private `silva-dbs` R2 bucket. Anyone
can request the allowed files without signing in or supplying credentials.
The large files are downloaded directly from R2. The Worker runs on Cloudflare;
its R2 signing credentials remain server-side.

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
4. In the Worker's **Access** tab, disable **Protect this Worker behind Access**
   for the production URL used by the pipeline. Deploying this code does not
   remove an existing dashboard Access policy. If the endpoint still redirects
   to a login page, check for a matching Access application in Zero Trust.
   The R2 bucket can stay private: the public Worker grants GET access only to
   the release and file names allowed in `src/keys.js`.
5. Confirm the `R2_ACCOUNT_ID`, bucket name, and uploaded object keys in
   `wrangler.jsonc` and `src/keys.js`. If the bucket uses a jurisdiction-specific
   endpoint, change the S3 hostname in `src/index.js` accordingly.

`package.json` and npm are used to build the Worker on Cloudflare or on a
developer machine. They are not needed in `Dockerfile.TACTIC`.

## API

Both endpoints accept anonymous HTTPS GET requests:

```text
GET /v1/releases
```

This returns a newest-first release list, including only allowlisted releases
for which a nonempty manifest and an allowed ARB and SIDX variant exist in R2.
Compressed `.zst` artifacts are supported; uncompressed artifacts are considered
only when enabled in `src/keys.js`. Add new releases there after uploading the
manifest and both database artifacts, then deploy the Worker.

```text
GET /v1/download-url?release=144&file=SILVA.sidx
```

The JSON response includes `url`, `key`, `size_bytes`, and
`expires_in_seconds`. The URL is valid for one hour and authorizes only a GET
of the requested object. Clients request a fresh URL when needed; no client
Access token, account, or API key is required.
The `file` parameter accepts only the exact names in `src/keys.js`: the
release manifest and SILVA ARB/SIDX artifacts. SortMeRNA files are hosted
separately and are not available through this Worker. A missing R2 object
returns 404.

The laptop client should fetch the release manifest first, then request the
required artifact. It must verify the downloaded size and SHA-256 against the
manifest before atomically installing it. If a large download fails after the
URL expires, request a new URL and retry. R2 object keys and hashes in the
manifest must correspond to the uploaded objects.

## Pipeline client

The current pipeline uses this Worker by default. It first calls `/v1/releases`,
then downloads candidate manifests until it finds one matching the machine
platform, byte order, SINA version, and SHA-256 of its local SINA binary.
For each artifact, the pipeline first downloads the uncompressed `SILVA.arb`
or `SILVA.sidx`. It requests the corresponding `.zst` file only when the
uncompressed object returns HTTP 404 and the manifest includes its
`arb_zst_sha256` or `index_zst_sha256`. Downloads are checked for size and
SHA-256; compressed files are decompressed in chunks using the Python
`zstandard` dependency and checked against the uncompressed hash. Verified
files are installed as `/databases/SILVA/<release>/SILVA.arb` and `SILVA.sidx`
and reused on later runs. Other download failures use the SILVA FTP fallback.

Windows, WSL, and Docker application instances download without a sign-in
step. The pipeline does not read `SILVA_ACCESS_TOKEN` or send authentication
headers. No npm, cloudflared, or R2 keys are needed on client machines.
`SILVA_WORKER_URL` overrides the Worker origin. To explicitly use the old FTP
download and local index build, set `SILVA_SOURCE=ftp`.
By default, unavailable Worker/R2 services, missing remote artifacts, or an
unavailable compatible release trigger the original SILVA FTP download and
local SINA index build. The pipeline logs this fallback and its RAM requirement;
index construction defaults to one thread. The requested release is preserved
(`latest` is resolved by SILVA FTP on fallback). Corrupt files, malformed
manifests, and local filesystem failures remain errors.

After deploying the Worker and disabling Access, check the production endpoint
without authentication (or open it in a private browser window):

```bash
curl --fail --silent --show-error https://silva-downloads-worker.mohsenpm50.workers.dev/v1/releases
```

Expect a JSON release list, not a sign-in page. Rebuild the pipeline image to
include the updated anonymous client.

## Local development

Run `npm install`, `npm test`, and `npm run dev` from this directory. For local
development only, place `R2_ACCESS_KEY_ID` and `R2_SECRET_ACCESS_KEY` in an
untracked `.dev.vars` file. Do not use real production secrets in tests. The
local R2 binding is simulated by default, so it will not see production bucket
objects unless you explicitly configure a remote binding.
