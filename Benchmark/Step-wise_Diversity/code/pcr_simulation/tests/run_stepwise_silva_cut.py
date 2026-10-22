from ..silva_algn_cut_and_count import generate_amplicon_fastas, SILVA_ALGN_PATHS, _count_unique_taxonomies, plot_from_stats_tsv


print("Generating amplicon FASTAs from complete SILVA DB...")

# print(_count_unique_taxonomies(SILVA_ALGN_PATHS[2]))
# exit()

generate_amplicon_fastas(
    SILVA_ALGN_PATHS[4],  # GG2 Full Species Taxonomy
    start_pos=1150,
    end_pos=43026,
    step_size=2000,
    min_amplicon_length=1200,
    max_amplicon_length=1600,
    output_dir="./outputs_GG2_cut_clust",
    clust_ids=[0.995, 0.987, 0.97]
)
