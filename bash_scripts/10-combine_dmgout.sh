#!/bin/bash

#===========================================================================
# Combine C>T rate of metaDMGout.csv and attach taxlineage
# on several sample using arrays
# Version 0.4 
#
# by Sisi Liu
# 
# contact: sisi.liu@awi.de
#
# slurm options and variables under >set required variables< 
# have to be modified by the user
#=============================================================================

#SBATCH --account=
#SBATCH --job-name=c_out
#SBATCH --partition=smp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

#=========================== CPU =================
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

#=========================== Modules =================
RTOOL="r/4.2.2"

#===================================================================
SCRIPT="/path/to/scripts"

WORK="/path/to/output/directory"
IN_FOLDER="out.metaDMG"
IN_FILE_END="_metaDMGout.csv"
IN_FILE_SAMPLE_END="_L30"

OUT_FOLDER="out.metaDMG.post"
OUT_FILE="lele_metaDMGout"

SAMPLE_INFOR="/path/to/sample/metadata/file.csv"
TAXLINEAGE="/path/to/taxonomy/Siberian_taxlineage.RData"

#====  processing ====
cd ${WORK}
mkdir -p ${WORK}/${OUT_FOLDER}

module load ${RTOOL}
# run R with args
srun Rscript ${SCRIPT}/combine_dmgout.R $IN_FOLDER $IN_FILE_END $IN_FILE_SAMPLE_END $SAMPLE_INFOR $TAXLINEAGE $OUT_FOLDER $OUT_FILE 

module unload ${RTOOL}
