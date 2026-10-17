# USEARCH 11

The supplied `usearch11.0.667_i86linux64` is the exact executable in the source
analysis directory. It reports `usearch v11.0.667_i86linux64`.

MD5: `fa050a3029d33b7b25036a2cc8c6da97`.

This matches the author's published checksum for the same binary at
https://github.com/rcedgar/usearch_old_binaries/ . That repository states that
USEARCH versions 5 through 11 have been donated to the public domain under CC0-1.0.
The upstream license is at
https://github.com/rcedgar/usearch_old_binaries/blob/main/LICENSE .
The SHA-256 is recorded in `metadata/source_files.tsv` and each run's `software.json`.

The binary targets Linux x86-64. The launcher resolves it relative to the bundle,
so it does not require installation or the original working directory.
If executable permissions are lost during transfer, run:

```bash
chmod +x tools/usearch11.0.667_i86linux64
```
