# EnviHoli
Code repository for sedimentary ancient DNA (sedaDNA) shotgun sequencing (metagenomics) data analysis.  
Contributed by: [Sisi Liu](mailto:sisi.liu@awi.de)/[Sisi Liu](mailto:sisi.liu.research@gmail.com), [Lars Harms](mailto:lars.harms@awi.de), [Christiane Böckel](mailto:christiane.boeckel@awi.de), [Kathleen R. Stoof-Leichsenring](Kathleen.Stoof-Leichsenring@awi.de)

# General content
## Scripts and Files in This Repository

### Bash Scripts
This repository contains 10 bash scripts located in the `bash_scripts` directory. These scripts are:

00-bowtie2-build-bac0.sh  
01-fastqc-clumpify-fastp-dedupe.sh  
02-bowtie2-a1.sh  
03-bowtie2-a2.sh  
04-bowtie2-a3.sh  
05-merge-sort.sh  
06-metaDMG.sh  
07-post-metadmg-lca.sh  
08-combine_lca.sh  
09-mismatch.sh  
10-combine_dmgout.sh  

### External Scripts
In addition, there are 5 external scripts in the `external_scripts` directory, which are written in various languages (Python and R) and are not part of the bash scripts.  

1. combine_dmgout.R  
2. combine_lca.R  
3. dedup_sam.py  
4. mismatch.R  
5. post-metadmg-lca.R  

### Files
In addition, there are 10 files in the `external_files` directory, which are descriptions of sources of raw shotgun sequencing data (raw_shotgun_data_sources.txt), taxonomic reference data (taxonomic_reference_database.xlsx), and age-depth models of 8 lake cores (e.g., Age-depth_*_shotgun.csv)

# Detailed Description
This section provides an in-depth look at the data analysis's features and functionality.

## Installation

Before running the scripts, make sure to install the required dependencies.
1. [Fastqc](https://anaconda.org/bioconda/fastqc)  
2. [bbmap](https://anaconda.org/bioconda/bbmap)  
3. [fastp](https://anaconda.org/bioconda/fastp)
4. [bowtie2](https://anaconda.org/bioconda/bowtie2)
5. [samtools](https://anaconda.org/bioconda/samtools)
6. [gz-sort](https://github.com/keenerd/gz-sort)
7. [metaDMG](https://github.com/miwipe/metaDMG_installation)  
8. [r-base](https://anaconda.org/conda-forge/r-base)
9. [python](https://anaconda.org/conda-forge/python)
10. [fasta-splitter](https://kirill-kryukov.com/study/tools/fasta-splitter/)

## Dependencies' manuals

Instructions on how to use the dependencies.
1. [Fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  
2. [bbmap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/)  
3. [fastp](https://github.com/OpenGene/fastp)  
4. [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
5. [samtools](https://www.htslib.org/doc/samtools.html)
6. [gz-sort](https://github.com/keenerd/gz-sort)
7. [metaDMG](https://metadmg-dev.github.io/metaDMG-core/)
8. [fasta-splitter](https://kirill-kryukov.com/study/tools/fasta-splitter/)

## Usage

### I. Shotgun sequencing data quality check -> deduplication -> adapter trimming and merging of paired-end reads in parallel -> deduplication -> quality check

1. Input raw sequencing paired end fastq files: there are two files, ${FILEBASE}.R1.fastq.gz and ${FILEBASE}.R2.fastq.gz (or ${FILEBASE}_R1.fastq.gz and ${FILEBASE}_R2.fastq.gz, depending on sequencing company), for each sequencing id ${FILEBASE}.
2. Script: bash_scripts/01-fastqc-clumpify-fastp-dedupe.sh
3. Output for next step (alignment): *fastp_dedupe_merged.fq.gz

### II. Taxonomic reference database establishment and end-to-end alignment in Bowtie2

1. Source data for taxonomic reference database establishment: external_files/taxonomic_reference_sources.txt
2. Script for taxonomic reference database establishment: bash_scripts/00_0-fasta-splitter-bacteriaRefseq.sh and 00_1-bowtie2-build-bac0.sh (bowtie2-build for Bacteria refseq database establishmen. Other database using the same script with different path-to-db and splited size)
3. Input merged shotgun sequencing data: *fastp_dedupe_merged.fq.gz
4. Script for alignment against taxonomic reference database: bash_scripts/02-bowtie2-a1.sh, bash_scripts/03-bowtie2-a2.sh, bash_scripts/04-bowtie2-a3.sh
5. Output for next step (merge and sort alignments): ${FILEBASE}.$(basename $DB).bam (${FILEBASE} is fastq file id; $(basename $DB) is taxonomic reference database name. In total, there are 147 alignment bam files per seqencing file.)

### III. Merge and sort alignments
Motivation: To make sure alignments have been sorted by readID; sort the sam file instead of bam file due to size of headers of merged bam file > 2GB.
1. Input all alignments:${FILEBASE}.$(basename $DB).bam
2. Script for merge and sort: bash_scripts/05-merge-sort.sh
3. Outoup for next step (taxonomic classification and ancient damage pattern analysis):${FILEBASE}_L30.sorted.sam.gz

### IV. Taxonomic classification and ancient damage pattern analysis
Taxonomic profile: [Wang et al., 2022](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.14006); 
Ancient pattern: [Michelsen et al., 2022](https://www.biorxiv.org/content/10.1101/2022.12.06.519264v1)
1. Input sorted alignments:${FILEBASE}_L30.sorted.sam.gz
2. Script: bash_scripts/06-metaDMG.sh
3. Output structure: see [Michelsen et al., 2022](https://www.biorxiv.org/content/10.1101/2022.12.06.519264v1)
4. Taxonomic reference dump files used for classification are archived on Zenodo: 10.5281/zenodo.17974857

### V. Post-processing of MetaDMG
1. Attach full lineage and key ranks based on tax_id: bash_scripts/07-post-metadmg-lca.sh and external_scripts/post-metadmg-lca.R
2. Combine taxonomic classification results: bash_scripts/08-combine_lca.sh and external_scripts/combine_lca.R
3. Calculate ATCG substitutions frequency and attach lineage information: bash_scripts/09-mismatch.sh and external_scripts/mismatch.R
4. Combine C>T rate of metaDMGout.csv and attach lineage information: bash_scripts/10-combine_dmgout.sh and external_scripts/combine_dmgout.R





