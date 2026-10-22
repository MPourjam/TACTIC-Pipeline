# Minimal Reproduction Environment

The primary alignment-cut workflow imports NumPy and Matplotlib, calls SeqKit,
and uses the supplied USEARCH 11 executable. All remaining imports in that
workflow are Python standard-library modules.

`stepwise-diversity-reproduction.yaml` therefore declares only Python, NumPy,
Matplotlib's non-GUI package, and SeqKit. Conda resolves their required libraries.
USEARCH 11 is bundled separately under `tools/` and is recorded by checksum.

`tested-linux-64.explicit.txt` records the full solved package set used for the
isolated reproduction checks. To recreate that exact Linux environment:

```bash
conda create -n stepwise-diversity-reproduction --file environment/tested-linux-64.explicit.txt
conda activate stepwise-diversity-reproduction
```

Use either this explicit specification or the shorter YAML specification when
creating the environment. The explicit file is platform-specific; the YAML pins
the direct analysis tool versions while allowing dependency resolution.

The historical primer-sweep and sequence-length plotting helpers additionally
import SciPy or call Cutadapt. The old Nextflow file also requires Nextflow/Java
and refers to historical script paths. Those branches are preserved for context
and are not invoked by the portable reproduction commands.
