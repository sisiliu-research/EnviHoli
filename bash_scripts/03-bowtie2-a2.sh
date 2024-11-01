#!/bin/bash

#===========================================================================
# Alignment against taxonomic reference database package 2 on several sample using arrays
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
#SBATCH --job-name=l2a2
#SBATCH --partition=smp
#SBATCH --time=48:00:00
#SBATCH --qos=48h
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=64
#SBATCH --mem=220G
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

mkdir -p ${WORK}/${OUT_ALIGN}
mkdir -p ${WORK}/${OUT_ALIGN}/${FILEBASE}

#== Invertebrate refseq
# which path 
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/bowtie2_db_23_08"

# whih DB?
for DB in ${DBP}/invert00/invert.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/invert01/invert.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== NCBI NT database, 2021 January
for DB in ${DBP}/nt2100/library.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2101/library.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2102/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2103/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2104/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2105/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2106/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

for DB in ${DBP}/nt2107/library.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# remove duplicates
cd ${SAMF}
for FILE_IN in ${FILEBASE}.library.{1..40}.sam; do\
	FILE_END="${FILE_IN%.sam}"
	FILE_OUT="${FILE_END}.dedup.sam"
	echo Mapping ${FILE_IN} ${FILE_OUT}
	srun ${SCRIPTSP}/dedup_sam.py $FILE_IN $FILE_OUT
	echo "dedup sam to bam"
	srun samtools view -bS ${SAMF}/$FILE_OUT > ${SAMF}/${FILE_END}.bam
done

rm ${SAMF}/*.library.{1..40}.sam
rm ${SAMF}/*.library.{1..40}.dedup.sam

#==
module unload ${ALIGN}
module unload ${VIEW}


