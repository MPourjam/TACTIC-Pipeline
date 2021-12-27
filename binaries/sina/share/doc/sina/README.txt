SINA - reference based multiple sequence alignment
**************************************************

[image: latest][image] [image: Bioconda][image] [image:
downloads][image] [image: TravisCI][image] [image: CircleCI][image]
[image: Read the Docs][image] [image: Codecov][image]

SINA aligns nucleotide sequences to match a pre-existing MSA using a
graph based alignment algorithm similar to PoA. The graph approach
allows SINA to incorporate information from many reference sequences
building without blurring highly variable regions. While pure NAST
implementations depend highly on finding a good match in the reference
database, SINA is able to align sequences relatively distant to
references with good quality and will yield a robust result for query
sequences with many close reference.


Features
========

* Speed. Aligning 100,000 full length rRNA against the SILVA NR
  takes 40 minutes on a mid-sized 2018 desktop computer. Aligning
  1,000,000 V4 amplicons takes about 60 minutes.

* Accuracy. SINA is used to build the SILVA SSU and LSU rRNA
  databases.

* Classification. SINA includes an LCA based classification module.

* ARB. SINA is able to directly read and write ARB format files such
  as distributed by the SILVA project.


Online Version
==============

An online version for submitting small batches of sequences is made
available by the SILVA project as part of their ACT: Alignment,
Classification and Tree Service. In addition to SINA's alignment and
classification stages, ACT allows directly building phylogenetic trees
with RAxML or FastTree from your sequences and (optionally) additional
sequences chosen using SINA's add-neighbors feature.


Installing SINA
===============

The preferred way to install SINA locally is via Bioconda. If you have
a working Bioconda installation, just run:

   conda create -n sina sina
   conda activate sina

Alternatively, self-contained images are available at
https://github.com/epruesse/SINA/releases. Choose the most recent
"tar.gz" appropriate for your operating system and unpack:

   tar xf sina-1.6.0-linux.tar.gz
   cd sina-1.6.0
   ./sina


Documentation
=============

The full documentation is available at https://sina.readthedocs.io.

The algorithm is explained in the paper:

   Elmar Pruesse, Jörg Peplies, Frank Oliver Glöckner; *SINA: Accurate
   high-throughput multiple sequence alignment of ribosomal RNA
   genes.* Bioinformatics 2012; 28 (14): 1823-1829.
   doi:10.1093/bioinformatics/bts252
