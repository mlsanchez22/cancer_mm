#!/bin/bash
#SBATCH -J dr-schizo-launcher
#SBATCH -t 0
#SBATCH --mem 8G
#SBATCH -n 1
#SBATCH -o pipeline_custom_circuits.log

###################################################
# SET THE FOLLOWING VARIABLES TO MATCH YOUR NEEDS #
###################################################
SINGULARITY_MOD='Singularity/3.9.7'
MAMBA_MOD='Mamba/4.14.0-0'
PARAMS_FILE='/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/code/09_drug_repository/09_4_drexml_gtex/pancreas_tcga_paad/params_custom_circuits_male.yml'
CONDA_INIT_SCRIPT='/mnt/lustre/expanse/software/easybuild/Rocky/8.5/Skylake/software/Mamba/4.14.0-0/etc/profile.d/conda.sh'
CONDA_ENV='nextflow'
UPDATE_NEXTFLOW='TRUE'

# The PIPELINE_REPO variable must point to the github repository (user/project) or the main.nf full path (in your local -cloned- repository)
#PIPELINE_REPO='babelomics/dr-pipeline'
PIPELINE_REPO='/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/code/09_drug_repository/dr-pipeline/main.nf'

########
# INIT #
########
# check params file
[[ -f "${PARAMS_FILE}" ]] || { echo "Params file ${PARAMS_FILE} not found!" && exit 1; }

# check if singularity module is loaded
[[ $(module is-loaded ${SINGULARITY_MOD}) ]] || module load ${SINGULARITY_MOD}

# check if mamba module is loaded
[[ $(module is-loaded ${MAMBA_MOD}) ]] || module load ${MAMBA_MOD}

# initialize conda
source ${CONDA_INIT_SCRIPT} || { echo "The conda environment couldn't be initialized!" && exit 1; }

# check conda env
if [[ ! $(conda info --envs |grep ${CONDA_ENV}) ]]
then
  echo "${CONDA_ENV} environment not found. Trying to create a new conda environment..."
  conda create -y -n ${CONDA_ENV} bioconda::nextflow
fi

# activate the conda environment
conda activate ${CONDA_ENV} || { echo "The conda environment couldn't be activated!" && exit 1; }

# update nextflow
[[ "${UPDATE_NEXTFLOW}" ]] && nextflow self-update

########
# MAIN #
########
# run the workflow from the github repository
nextflow run ${PIPELINE_REPO} -profile slurm,singularity -resume -params-file ${PARAMS_FILE}

