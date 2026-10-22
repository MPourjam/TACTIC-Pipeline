
nextflow.enable.dsl=2

/*
 * Parameter definitions
 * ---------------------
 * Adjust these to control the sweep.
 */

params.ref_fasta          = 'data/silva_138.2_SSURef_NR99.fasta'
params.primer_table       = 'primers.tsv'

// scripts moved into ./scripts
params.stepwise_script    = 'scripts/stepwise_16S_diversity.py'
params.fasta_stats_script = 'scripts/fasta_stats.py'

params.outdir_root        = 'results'          // top-level output directory

// Cutadapt-related parameter ranges
params.error_rates   = [0.1, 0.2]             // --error-rate
params.min_overlaps  = [10, 15]               // --overlap
params.min_lengths   = [100, 200]             // --minimum-length
params.max_lengths   = [1000, 1500]           // --maximum-length

// Primer sweeps
params.forward_primers = ['27F']              // names in primers.tsv with direction=F
// Each entry is a string exactly as you would pass to --reverse-targets
params.reverse_sets     = [
    '338R,519R,785R,926R,1492R',
    '519R,785R'
]

// Misc options
params.rc_reverse = true                      // pass --rc-reverse or not
params.threads    = 20                         // threads for cutadapt

/*
 * Normalize scalar/list params to lists (for safety)
 */

def asList(x) {
    (x instanceof List) ? x : [ x ]
}

def ERROR_RATES     = asList(params.error_rates)
def MIN_OVERLAPS    = asList(params.min_overlaps)
def MIN_LENGTHS     = asList(params.min_lengths)
def MAX_LENGTHS     = asList(params.max_lengths)
def FORWARD_PRIMERS = asList(params.forward_primers)
/*
 * Main process:
 * For each parameter combination:
 *   1. Run stepwise_amplicons.py simulate
 *   2. For each resulting *.fasta, run fasta_stats.py
 *   3. Write *.fasta_stats.txt alongside the fasta in the same run directory
 */

process RUN_STEPWISE_AND_STATS {

    tag { run_id }

    cpus params.threads
    publishDir "${params.outdir_root}", mode: 'copy', saveAs: { it }

    input:
    tuple val(run_id), val(p)
    path ref_fasta
    path primer_table
    path stepwise_script
    path fasta_stats_script

    output:
    // The whole run directory with FASTAs, summary table and stats
    path "${run_id}"

    script:
    // Build optional flag for rc_reverse
    def rcFlag = params.rc_reverse ? '--rc-reverse' : ''
    """
    set -euo pipefail

    mkdir -p ${run_id}

    echo "[INFO] Running stepwise_amplicons for: ${run_id}"

    python ${stepwise_script} simulate \\
      --ref ${ref_fasta} \\
      --primers ${primer_table} \\
      --forward ${p.forward} \\
      --reverse-targets '${p.reverse_targets}' \\
      ${rcFlag} \\
      --min-length ${p.min_length} \\
      --max-length ${p.max_length} \\
      --error-rate ${p.error_rate} \\
      --min-overlap ${p.min_overlap} \\
      --threads ${task.cpus} \\
      --outdir ${run_id}

    echo "[INFO] Running fasta_stats.py on all FASTA outputs for: ${run_id}"

    shopt -s nullglob
    for f in ${run_id}/*.fasta; do
        base=\$(basename "\$f" .fasta)
        python ${fasta_stats_script} "\$f" > "${run_id}/\${base}.fasta_stats.txt"
    done

    echo "[INFO] Finished run: ${run_id}"
    """
}

workflow {
    /*
     * Channels: reference files / scripts
     */
    REF_FASTA          = Channel.value(file(params.ref_fasta))
    PRIMER_TABLE       = Channel.value(file(params.primer_table))
    STEPWISE_SCRIPT    = Channel.value(file(params.stepwise_script))
    FASTA_STATS_SCRIPT = Channel.value(file(params.fasta_stats_script))

    /*
     * Build all parameter combinations using Nextflow channels
     */
    PARAM_COMBOS = Channel
        .fromList(ERROR_RATES)
        .combine(Channel.fromList(MIN_OVERLAPS))
        .combine(Channel.fromList(MIN_LENGTHS))
        .combine(Channel.fromList(MAX_LENGTHS))
        .combine(Channel.fromList(FORWARD_PRIMERS))
        .combine(Channel.fromList(REVERSE_SETS))
        .map { er, ov, minL, maxL, fwd, revset ->
            def clean_revset = revset.replace(',', '-')
            def run_id = "er${er}_ov${ov}_min${minL}_max${maxL}_${fwd}_${clean_revset}"
            
            tuple(
                run_id,
                [
                    error_rate     : er,
                    min_overlap    : ov,
                    min_length     : minL,
                    max_length     : maxL,
                    forward        : fwd,
                    reverse_targets: revset
                ]
            )
        }

    RUN_STEPWISE_AND_STATS(
        PARAM_COMBOS,
        REF_FASTA,
        PRIMER_TABLE,
        STEPWISE_SCRIPT,
        FASTA_STATS_SCRIPT
    )
}
    """
}
