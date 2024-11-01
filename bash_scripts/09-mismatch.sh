#!/bin/bash 

#===== bash + Call external R to do ======
# 1. load *.mismatches.txt.gz (a tab-separated ( sep = “\t”) file)
# 2. calculate ATCG substitutions frequency
# 3. add taxlineage
# Contact: Sisi Liu, sliu@awi.de

#=========================== Envi Albedo: case-dependent ================= 

#SBATCH --account=
#SBATCH --job-name=il_mis
#SBATCH --partition=smp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=

#=========================== CPU =================
CPU=${SLURM_CPUS_PER_TASK}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

#=========================== Modules =================
RTOOL="r/4.2.2"

#=========================== Envi working: case-dependent =================
SCRIPT="/path/to/scripts"

WORK="/path/to/output/directory"
IN_FOLDER="out.metaDMG"
IN_FILE_END="_L30.mismatches.txt.gz"

OUT_FOLDER="out.metaDMG.post"
OUT_FILE="lele_mismatch" # <lake_name>_mismatch

SAMPLE_INFOR="/path/to/sample/metadata/file.csv"
TAXLINEAGE="/path/to/taxonomy/Siberian_taxlineage.RData"

#=========================== R Processing =================
module load ${RTOOL}

# run R with args
cd ${WORK}
srun Rscript ${SCRIPT}/mismatch.R $IN_FOLDER $IN_FILE_END $SAMPLE_INFOR $TAXLINEAGE $OUT_FOLDER $OUT_FILE

module unload ${RTOOL}



