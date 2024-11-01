#!/bin/bash

#===========================================================================
# Taxonomic classification and ancient damage pattern analysis
# on several sample using arrays
# Version 0.5 
#
# by Sisi Liu
# 
# contact: sisi.liu@awi.de
#
# slurm options and variables under >set required variables< 
# have to be modified by the user
#=============================================================================

#SBATCH --account=
#SBATCH --job-name=idmg1
#SBATCH --partition=smp
#SBATCH --time=24:00:00
#SBATCH --qos=48h
#SBATCH --array=1-20%20
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

# given variables (please do not change)
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# modules
VIEW="samtools/1.16.1"
DMG_CORE="metaDMG/0.37"
DMG_CPP_PATH="/path/to/metadmg_test/0.38/metaDMG-cpp/metaDMG-cpp"

#===================================================================
WORK="/path/to/output/directory"
ALIGN="/path/to/bowtie2/output/directory"
OUT_DMG="out.metaDMG"
# file name
cd ${ALIGN}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))

#== variables
# sorted sam file
OUT_SORT="${FILEBASE}_L30.sorted.sam.gz"
# path to sorted sam
PATH_SORT="${ALIGN}/${FILEBASE}/${OUT_SORT}"

echo ${PATH_SORT}
echo ${PATH_SORT_MOD}

#== load
module load ${VIEW}
module load ${DMG_CORE}

#== Processing
# make folder for output
mkdir -p ${WORK}/${OUT_DMG}/${FILEBASE}
cd ${WORK}/${OUT_DMG}/${FILEBASE}
# create config file
nam="/albedo/work/projects/p_biodiv_dbs/Siberian_taxonomies/Siberian_names.dmp"
nod="/albedo/work/projects/p_biodiv_dbs/Siberian_taxonomies/Siberian_nodes.dmp"
acc="/albedo/work/projects/p_biodiv_dbs/Siberian_taxonomies/Siberian_modify_combined_acc2taxid_v1.gz"

echo "create yaml"
srun metaDMG config ${PATH_SORT} --names $nam --nodes $nod --acc2tax $acc --parallel-samples 1 --cores-per-sample ${CPU} --custom-database --bayesian --config-file config_${FILEBASE}.yaml -m ${DMG_CPP_PATH}

echo "compute"
srun metaDMG compute ${WORK}/${OUT_DMG}/${FILEBASE}/config_${FILEBASE}.yaml

echo "csv"
srun metaDMG convert ${WORK}/${OUT_DMG}/${FILEBASE}/config_${FILEBASE}.yaml --output ${FILEBASE}_metaDMGout.csv

echo "pdf"
srun metaDMG plot ${WORK}/${OUT_DMG}/${FILEBASE}/config_${FILEBASE}.yaml --output ${FILEBASE}.pdf

#== unload
module unload ${DMG_CORE}
module unload ${VIEW}
