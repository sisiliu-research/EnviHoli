#!/bin/bash

#=====
# slurm to split refseq into 40 parts by considering computational capacity of HPC
# "how many parts", denpending on file size of *.fa and HPC 
# Taking bacteria refseq as an example
#
# by Sisi Liu
# 
# contact: sisi.liu@awi.de
#
#=====

#SBATCH --account=
#SBATCH --job-name=b0in
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

module load perl/5.35.0-gcc12.1.0
module load bowtie2/2.5.1


export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cd /path/to/bowtie2db/bac

srun cat /path/to/refseq/bacteria/*.fna > /path/to/bowtie2db/bac/bac.fa
srun /path/to/programmes/fasta-splitter.pl --n-parts 40 bac.fa
srun rm bac.fa

i=1
for file in bac.part*; do\
  bname=$(basename "$file" | cut -d. -f1)
  mv $file $bname.$i
  i=$(expr ${i} + 1)
done
#---- END to split whole bacteria_refseq database into 40 parts ----



