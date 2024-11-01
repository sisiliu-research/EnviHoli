#!/bin/bash 

#===== bash + Call external R to do ======
# 1. load *.lca.txt.gz (a tab-separated ( sep = “\t”) file)
# 2. R, aggragate reads by taxid
# 3. R, attach info. from rankedlineage.dmp and taxidlineage.dmp, 
# 4. R, attach ngs IDs (e.g., 191011_SND405_A_L006_APMG-5-1..)

# Contact: Sisi Liu, sliu@awi.de

#=========================== Envi Albedo: case-dependent ================= 

#SBATCH --account=
#SBATCH --job-name=l1_lca
#SBATCH --partition=smp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --array=1-22%22
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
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
LCAPATH="/path/to/metaDMG/output"

# taxonomy files (NCBI tree), upto date
FULL="/path/to/p_biodiv_dbs/Siberian_taxonomies/Siberian_fullnamelineage.dmp"
RANK="/path/to/p_biodiv_dbs/Siberian_taxonomies/Siberian_rankedlineage.dmp"
NODES="/path/to/p_biodiv_dbs/Siberian_taxonomies/Siberian_nodes.dmp"

# lca gz file
LCA_GZ="_L30.lca.txt.gz"

#=========================== Which LCA_GZ ? =================
# pass the list of folder names to FILEBASE, orders to slurm jobs (assume folder name == ngs id, e.g., 191011_SND405_A_L006_APMG-5-1)
cd ${LCAPATH}
FILEBASE=$(basename $(find . -mindepth 1 -maxdepth 1 -type d | sed -n ${SLURM_ARRAY_TASK_ID}p))
# location of LCA_GZ
WORKPATH="${LCAPATH}/${FILEBASE}/data/lca"

#=========================== R Processing =================
module load ${RTOOL}
# run R with args
srun Rscript ${SCRIPT}/post-metadmg-lca.R $WORKPATH $FILEBASE $LCA_GZ $FULL $RANK $NODES

#=========================== unload modules =================
module unload ${RTOOL}
