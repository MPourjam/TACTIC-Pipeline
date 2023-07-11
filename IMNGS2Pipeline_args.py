# The config yaml file for IMNGS2 Pipeline The arguemnts should be placed with the
# name of pipeline step (CamelCase) and the arguments added as indentaion to.
# If an argument parser exists with the same name of pipeline step (e.g: MergePairsArgs),
# then the arguments would be implemented in processing.
#####################
### Preprocessing ###
#####################
MergePairsArgs:
  fastq_maxdiffs: 50
  fastq_pctid: 50
  fastq_minmergelen: 200
  fastq_maxmergelen: 600
TrimBothSidesArgs:
  stripleft: 5
  stripright: 5
TrimOneSideArgs:
  stripleft: 5
FilterBothSidesArgs:
  fastq_maxee_rate: 0.002
FilterOneSideArgs:
  fastq_maxee_rate: 0.002
  fastq_truncqual: 10
DereplicationArgs:
  sizein: true
  sizeout: true
# sort_sequences:
ClusterZOTUsArgs:
  minsize: 2
Filter16SArgs:
  e: 0.1
  num_alignments: 1
# prepare_zotus:
BuildZOTUTableArgs:
  id: 0.97
FilterZOTUAbundanceArgs:
  abundance_cutoff: 0.0
# select_zotu_seqs:
AddTaxArgs:
  turn: all
################
### Analysis ###
################
TrimSidesArgs:
  stripleft: 5
  stripright: 5
ComplexTICArgs:
  family_sim: 0.90
  genus_sim: 0.95
  species_sim: 0.97