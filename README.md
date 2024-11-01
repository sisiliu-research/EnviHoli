# EnviHoli
Code repository for sedimentary ancient DNA (sedaDNA) shotgun sequencing (metagenomics) data analysis  
# General content
## Scripts in This Repository

### Bash Scripts
This repository contains 10 bash scripts located in the `bash_scripts` directory. These scripts are:

01-fastqc-clumpify-fastp-dedupe.sh
02-bowtie2-a1.sh
03-bowtie2-a2.sh
04-bowtie2-a3.sh
05-merge-split-Siberian.sh
06-metaDMG.sh
07-post-metadmg-lca.sh
08-combine_lca_il.sh
09-mismatch-ili.sh
10-combine_dmgout.sh

### External Scripts
In addition, there are 5 external scripts in the `external_scripts` directory, which are written in various languages (e.g., Python, Ruby) and are not part of the bash scripts.
combine_dmgout.R
combine_lca.R
dedup_sam.py
mismatch.R
post-metadmg-lca.R
