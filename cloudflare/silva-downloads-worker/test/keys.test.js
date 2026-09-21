import assert from "node:assert/strict";
import test from "node:test";
import { objectKey } from "../src/keys.js";

test("maps the uploaded SILVA release artifacts", () => {
  assert.equal(
    objectKey("144", "SILVA.sidx"),
    "silva/144/sina-1.6.0/k10-fast/linux-amd64/SILVA.sidx",
  );
  assert.equal(
    objectKey("138.2", "manifest.json"),
    "silva/138.2/sina-1.6.0/k10-fast/linux-amd64/manifest.json",
  );
  assert.equal(
    objectKey("144", "SILVA.sidx.zst"),
    "silva/144/sina-1.6.0/k10-fast/linux-amd64/SILVA.sidx.zst",
  );
  assert.equal(
    objectKey("144", "silva-arc-16s-id95.fasta"),
    "silva/144/sortmerna/silva-arc-16s-id95.fasta",
  );
});

test("rejects unknown releases, paths, and sensitive object names", () => {
  assert.equal(objectKey("145", "SILVA.sidx"), null);
  assert.equal(objectKey("144/../138.2", "SILVA.sidx"), null);
  assert.equal(objectKey("144", "../SILVA.sidx"), null);
  assert.equal(objectKey("144", "fake1.fasta"), null);
  assert.equal(objectKey("144", "SILVA_DB.tar.gz"), null);
});
