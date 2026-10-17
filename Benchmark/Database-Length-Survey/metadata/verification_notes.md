# Verification Record

Packaging and checks were performed on 2026-10-16 (UTC).

| Check | Result |
|---|---|
| Original survey scripts | All eight bundled copies match their source SHA-256 values |
| Small zOTU FASTAs | All 62 distinct payload checksums match the input manifest |
| Large archive | All 276 members pass ZIP CRC / inventory checks; 275 payload/manifest members match SHA-256, plus the checksum list itself |
| Available measurements | 541 file/script observations reference available inputs; 537 have saved numerical baselines and all match |
| Missing historical measurements | 94 explicitly recorded; complete reproduction of every historical report is not claimed |
| Python 3.10 unit tests | Seven passed, including legacy statistic equivalence, empty reads/files, malformed FASTQ rejection, checksum rejection and ZIP/extracted runner modes |
| zOTU comparison task | All 66 recovered observations reproduced, zero mismatches |
| Global zOTU statistics task | All 64 recovered observations reproduced, zero mismatches |
| FASTQ task checks | All six FASTQ task selections exercised against the real ZIP, one or two source paths per task, zero mismatches |
| Conda specifications | Both YAML files parse; survey needs only Python, while TACTICBench explicitly includes its declared Seaborn dependency |
| Existing TACTICBench data ZIP | Outer SHA-256 rechecked successfully, unchanged |

The real-data runner outputs and their `validation.json` summaries are preserved
under `results/validation_checks/`. The packaging scanner measured every included
payload; an additional complete second FASTQ survey over the entire ZIP was not
run. Fresh Conda environments were not installed. Tests used the existing
Python 3.10 TACTICBench environment; no upstream TACTIC/DADA2 pipeline was run.

The four full-read ZIEL-I observations lack a saved numerical baseline and are
fresh measurements from the recovered raw inputs. Do not interpret them as four
independently verified historical numeric results.
