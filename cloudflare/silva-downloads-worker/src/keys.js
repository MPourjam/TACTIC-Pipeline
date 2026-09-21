// Only these published releases and artifact names can receive download URLs.
const releases = new Set(["138.2", "144"]);
const databaseFiles = new Set([
  "manifest.json",
  "SILVA.arb",
  "SILVA.arb.zst",
  "SILVA.sidx",
  "SILVA.sidx.zst",
]);
const sortmernaFiles = new Set([
  "silva-arc-16s-id95.fasta",
  "silva-bac-16s-id90.fasta",
]);

export function objectKey(release, file) {
  if (!releases.has(release)) return null;
  if (databaseFiles.has(file)) {
    return `silva/${release}/sina-1.6.0/k10-fast/linux-amd64/${file}`;
  }
  if (sortmernaFiles.has(file)) {
    return `silva/${release}/sortmerna/${file}`;
  }
  return null;
}
