#!/bin/bash

#========
# slurm to build taxonomic reference database
# taking /bac00 as an example, including /bac00/bac.{1..9}
# by Sisi Liu
# 
# contact: sisi.liu@awi.de
#========

#SBATCH --account=
#SBATCH --job-name=b0in
#SBATCH --partition=fat
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

module load bowtie2/2.5.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cd /path/to/bowtie2db/bac00

for file in bac.?; do\
  srun bowtie2-build --threads ${SLURM_CPUS_PER_TASK} $file $file
done
