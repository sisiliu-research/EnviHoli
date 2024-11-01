#!/bin/bash

#===========================================================================
# Alignment against taxonomic reference database package 3 on several sample using arrays
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
#SBATCH --job-name=l2a3
#SBATCH --partition=smp
#SBATCH --time=20:00:00
#SBATCH --qos=48h
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=64
#SBATCH --mem=150G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


# set required variables (adapt according to your own requirements)
#===================================================================
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

ALIGN="bowtie2/2.5.1"
VIEW="samtools/1.16.1"

# given variables
#===================================================================
SCRIPTSP="/path/to/scripts"
WORK="/path/to/project/output/directory"
OUTDIR="output-01-18"
OUT_FASTP="out.dedupe"
OUT_ALIGN="out.bowtie2"

END_R1="_fastp_dedupe_R1.fq.gz"
END_R2="_fastp_dedupe_R1.fq.gz"
END_MERGED="_fastp_dedupe_merged.fq.gz"

INDIR="${WORK}/${OUT_FASTP}"

cd ${INDIR}
FILE=$(ls *${END_MERGED} | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILEBASE=${FILE%${END_MERGED}}
OUT_MERGED="${FILEBASE}_fastp_dedupe_merged.fq.gz"
SAMF="${WORK}/${OUT_ALIGN}/${FILEBASE}"

# prepare environment
#===================================================================
module load ${ALIGN}
module load ${VIEW}

#mkdir -p ${WORK}/${OUT_ALIGN}
#mkdir -p ${WORK}/${OUT_ALIGN}/${FILEBASE}

#== Other database
# which path 
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/bowtie2_db_add"

#== Phylonorway Contigs
for DB in ${DBP}/PhylonorwayContigs/PhylonorwayContigs.{1..10}; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== Selected mammalian genomes
for DB in ${DBP}/selected_mammals/selected_mammals.{1..10}; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== Selected plant mitochondrial genomes
DBi="selected_plants_mito/selected_plants_mito"
# full path to DB
DB="${DBP}/${DBi}"

echo Mapping ${OUT_MERGED} against $DB
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
rm ${SAMF}/${FILEBASE}.$(basename $DB).sam

#== Selected plant chloroplast genomes
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/build"
DBi="selected_plants/selected_plants_chloro"
# full path to DB
DB="${DBP}/${DBi}"

echo Mapping ${OUT_MERGED} against $DB
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
rm ${SAMF}/${FILEBASE}.$(basename $DB).sam

#== IMG/VR v4: an expanded database of uncultivated virus genomes
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/bowtie2_db_23_08"
for DB in ${DBP}/imgvr_clean/imgvr.clean.{1..5}; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done


#==
module unload ${ALIGN}
module unload ${VIEW}


