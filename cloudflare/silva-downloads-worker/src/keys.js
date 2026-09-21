// Only these published releases and artifact names can receive download URLs.
const releases = new Set(["138.2", "144"]);
const databaseFiles = new Set([
  "manifest.json",
  "SILVA.arb",
  "SILVA.arb.zst",
  "SILVA.sidx",
  "SILVA.sidx.zst",
]);

export function objectKey(release, file) {
  if (!releases.has(release)) return null;
  if (databaseFiles.has(file)) {
    return `silva/${release}/sina-1.6.0/k10-fast/linux-amd64/${file}`;
  }
  return null;
}
