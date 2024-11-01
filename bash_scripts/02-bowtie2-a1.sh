#!/bin/bash

#===========================================================================
# Alignment against taxonomic reference database 1 on several sample using arrays
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
#SBATCH --job-name=l2a1
#SBATCH --partition=smp
#SBATCH --time=48:00:00
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
SCRIPTSP="/path/to/scripts" # execute external python script: dedup_sam.py
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

#== archaea, fungi, viral refseq
# which path 
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/bowtie2_db_23_08"
# wich DB
DBi="afv/afv"
# full path to DB
DB="${DBP}/${DBi}"
#== mapping
echo Mapping ${OUT_MERGED} against $DB
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
rm ${SAMF}/${FILEBASE}.$(basename $DB).sam

#== Bacteria refseq
# whih DB?
for DB in ${DBP}/bac00/bac.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac01/bac.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac02/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac03/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac04/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac05/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac06/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/bac07/bac.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done
#== plant, protozoa, plastid, mitochondrion refseq
# whih DB?
for DB in ${DBP}/pppm/pppm.{1..4}; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# remove duplicates: share refseq between plant and mitochondrion package
FILE_IN="${FILEBASE}.pppm.5.sam"
FILE_OUT="${FILEBASE}.pppm.5.dedup.sam"
for DB in ${DBP}/pppm/pppm.5; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun ${SCRIPTSP}/dedup_sam.py ${SAMF}/$FILE_IN ${SAMF}/$FILE_OUT
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).dedup.sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== Vertebrate mammalian
# whih DB?
for DB in ${DBP}/vert_mam00/vert_mam.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/vert_mam01/vert_mam.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/vert_mam02/vert_mam.??; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== Vertebrate other
# whih DB?
for DB in ${DBP}/vert_other00/vert_other.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/vert_other01/vert_other.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

# whih DB?
for DB in ${DBP}/vert_other02/vert_other.?; do\
	echo Mapping ${OUT_MERGED} against $DB
	srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
	srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
	rm ${SAMF}/${FILEBASE}.$(basename $DB).sam
done

#== Ancient mammalian
# which path 
DBP="/path/to/p_biodiv_dbs/bowtie2_ngsLCA/build"
# wich DB
DBi="ancient_mammals/ancient_mammals"
# full path to DB
DB="${DBP}/${DBi}"

# mapping
echo Mapping ${OUT_MERGED} against $DB
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
rm ${SAMF}/${FILEBASE}.$(basename $DB).sam

#== Chloroplast genome from Stefano Meucci paper
# wich DB
DBi="chloroplast_stefano/chloroplast_stefano"
# full path to DB
DB="${DBP}/${DBi}"

# mapping
echo Mapping ${OUT_MERGED} against $DB
srun bowtie2 --threads ${CPU} -k 1000 -x $DB -U ${INDIR}/${OUT_MERGED} --no-unal -S ${SAMF}/${FILEBASE}.$(basename $DB).sam
srun samtools view -bS ${SAMF}/${FILEBASE}.$(basename $DB).sam > ${SAMF}/${FILEBASE}.$(basename $DB).bam
rm ${SAMF}/${FILEBASE}.$(basename $DB).sam


#==
module unload ${ALIGN}
module unload ${VIEW}


