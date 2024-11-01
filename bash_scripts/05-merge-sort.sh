#!/bin/bash

#===========================================================================
# Merge and sort alignments on several sample using arrays
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
#SBATCH --job-name=L1
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


# set required variables (adapt according to your own requirements)
#===================================================================
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

# given variables (please do not change)
#===================================================================

VIEW="samtools/1.16.1"
module load ${VIEW}

#== find the file name
INDIR="/path/to/bowtie2/output/directory"
cd ${INDIR}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
echo ${FILEBASE}
WORKFOLDER="${INDIR}/${FILEBASE}"
echo Working in ${WORKFOLDER}

#== Processing
# Merge BAM Files to SAM.GZ
cd ${WORKFOLDER}
srun samtools merge --verbosity 5 ${WORKFOLDER}/${FILEBASE}.L30.merged.sam.gz ${WORKFOLDER}/*.bam -@ ${CPU}

# Sort
srun samtools view -@ ${CPU} -H ${WORKFOLDER}/${FILEBASE}.L30.merged.sam.gz | sed '0,/^@HD/s/SO:unsorted/SO:queryname/' | pigz > ${WORKFOLDER}/${FILEBASE}.L30.Header.sam.gz
srun samtools view -@ ${CPU} ${WORKFOLDER}/${FILEBASE}.L30.merged.sam.gz | pigz > ${WORKFOLDER}/${FILEBASE}.L30.alignment.sam.gz
srun /path/to/tools/gz-sort/gz-sort -S 30G -P 30 ${WORKFOLDER}/${FILEBASE}.L30.alignment.sam.gz ${WORKFOLDER}/${FILEBASE}.L30.alignment.sort.sam.gz
srun zcat ${WORKFOLDER}/${FILEBASE}.L30.Header.sam.gz ${WORKFOLDER}/${FILEBASE}.L30.alignment.sort.sam.gz | samtools view -h -o ${WORKFOLDER}/${FILEBASE}_L30.sorted.sam.gz
rm ${WORKFOLDER}/${FILEBASE}.L30.Header.sam.gz ${WORKFOLDER}/${FILEBASE}.L30.alignment.sam.gz ${WORKFOLDER}/${FILEBASE}.L30.alignment.sort.sam.gz

#== unload
module unload ${VIEW}

