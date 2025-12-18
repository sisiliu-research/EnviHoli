#!/bin/bash

#===========================================================================
# slurm batch script to run the shotgun pipeline step 1
# on several sample using arrays
# Version 0.4 
#
# by Lars Harms
# 
# contact: lars.harms@awi.de
#
# slurm options and variables under >set required variables< 
# have to be modified by the user
#=============================================================================

#SBATCH --account=
#SBATCH --job-name=s1il1
#SBATCH --partition=smp
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --mem=100G
#SBATCH --array=1-20%4
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=


# set required variables (adapt according to your own requirements)
#===================================================================
WORK="/path/to/working/directory"
OUTDIR="/path/to/output/dir"
INDIR="/path/to/rawdata"
R1_ENDING=".R1.fastq.gz"
R2_ENDING=".R2.fastq.gz"
CLUMPIFY="TRUE"
DEDUPE="TRUE"

FASTQC="fastqc/0.11.9"
BBTOOLS="bbmap/39.01"
FASTP="fastp/0.23.2"

# given variables (please do not change)
#===================================================================
OUT_FQC="out.fastqc"
PRE="pre"
POST="post"
OUT_CLUMPIFY="out.clumpify"
OUT_FASTP="out.fastp"
OUT_DEDUPE="out.dedupe"

CPU=${SLURM_CPUS_PER_TASK}
MEM=100
# prepare environment
#===================================================================
mkdir -p ${OUTDIR}/${OUT_FQC}/${PRE}
mkdir -p ${OUTDIR}/${OUT_FQC}/${POST}
mkdir -p ${OUTDIR}/${OUT_CLUMPIFY}
mkdir -p ${OUTDIR}/${OUT_FASTP}
mkdir -p ${OUTDIR}/${OUT_DEDUPE}


cd ${INDIR}
FILE_R1=$(ls *${R1_ENDING} | sed -n ${SLURM_ARRAY_TASK_ID}p)
FILE_R2=$(ls *${R2_ENDING} | sed -n ${SLURM_ARRAY_TASK_ID}p)

FILEBASE=${FILE_R1%${R1_ENDING}}

OUT_R1_CL="${FILEBASE}_clumpify_R1.fq.gz"
OUT_R2_CL="${FILEBASE}_clumpify_R2.fq.gz"

OUT_R1="${FILEBASE}_fastp_R1.fq.gz"
OUT_R2="${FILEBASE}_fastp_R2.fq.gz"

OUT_MERGED="${FILEBASE}_fastp_merged.fq.gz"

OUT_R1_DD="${FILEBASE}_fastp_dedupe_R1.fq.gz"
OUT_R2_DD="${FILEBASE}_fastp_dedupe_R2.fq.gz"

OUT_MERGED_DD="${FILEBASE}_fastp_dedupe_merged.fq.gz"

cd ${WORK}

# tasks to be performed
#===================================================================

# FASTQC PRE
#----------
module load ${FASTQC}
srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${PRE} -t 2 ${INDIR}/${FILE_R1} ${INDIR}/${FILE_R2}
module unload ${FASTQC}

# CLUMPIFY
#----------
if [ ${CLUMPIFY} == "FALSE" ]; then
	echo "Removal of read duplications using clumpify is turned off."
else
	module load ${BBTOOLS}
	srun clumpify.sh in=${INDIR}/${FILE_R1} in2=${INDIR}/${FILE_R2} out=${OUTDIR}/${OUT_CLUMPIFY}/${OUT_R1_CL} out2=${OUTDIR}/${OUT_CLUMPIFY}/${OUT_R2_CL} dedupe=t optical=f
	module unload ${BBTOOLS}
fi

# FASTP
#----------
if [ ${CLUMPIFY} == "FALSE" ]; then
	module load ${FASTP}
	srun fastp --in1 ${INDIR}/${FILE_R1} --in2 ${INDIR}/${FILE_R2} --out1 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} --out2 ${OUTDIR}/${OUT_FASTP}/${OUT_R2} -m --merged_out ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED} -a auto --poly_g_min_len 10 --poly_x_min_len 10 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 30 --low_complexity_filter --complexity_threshold 30 --correction --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 -w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
	module unload ${FASTP}
else
	module load ${FASTP}
	srun fastp --in1 ${OUTDIR}/${OUT_CLUMPIFY}/${OUT_R1_CL} --in2 ${OUTDIR}/${OUT_CLUMPIFY}/${OUT_R2_CL} --out1 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} --out2 ${OUTDIR}/${OUT_FASTP}/${OUT_R2} -m --merged_out ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED} -a auto --poly_g_min_len 10 --poly_x_min_len 10 --qualified_quality_phred 15 --unqualified_percent_limit 40 --n_base_limit 5 --length_required 30 --low_complexity_filter --complexity_threshold 30 --correction --overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 -w ${CPU} --verbose --json=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.json --html=${OUTDIR}/${OUT_FASTP}/${FILEBASE}.html
	module unload ${FASTP}
fi

# DEDUPE
#----------
if [ ${DEDUPE} == "FALSE" ]; then
	echo "Removal of read duplications using dedupe is turned off."
else
	module load ${BBTOOLS}
	srun dedupe.sh -Xmx${MEM}g in=${OUTDIR}/${OUT_FASTP}/${OUT_R1} out=${OUTDIR}/${OUT_DEDUPE}/${OUT_R1_DD} ac=f
	srun dedupe.sh -Xmx${MEM}g in=${OUTDIR}/${OUT_FASTP}/${OUT_R2} out=${OUTDIR}/${OUT_DEDUPE}/${OUT_R2_DD} ac=f
	srun dedupe.sh -Xmx${MEM}g in=${OUTDIR}/${OUT_FASTP}/${OUT_MERGED} out=${OUTDIR}/${OUT_DEDUPE}/${OUT_MERGED_DD} ac=f
	module unload ${BBTOOLS}
fi

# FASTQC POST
#----------
if [ ${DEDUPE} == "FALSE" ]; then
	module load ${FASTQC}
	srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${POST} -t 3 ${OUTDIR}/${OUT_FASTP}/${OUT_R1} 	${OUTDIR}/${OUT_FASTP}/${OUT_R2} ${OUTDIR}/${OUT_FASTP}/${OUT_MERGED}
	module unload ${FASTQC}
else
	module load ${FASTQC}
	srun fastqc -q -o ${OUTDIR}/${OUT_FQC}/${POST} -t 3 ${OUTDIR}/${OUT_DEDUPE}/${OUT_R1_DD} 	${OUTDIR}/${OUT_DEDUPE}/${OUT_R2_DD} ${OUTDIR}/${OUT_DEDUPE}/${OUT_MERGED_DD}
	module unload ${FASTQC}
fi

