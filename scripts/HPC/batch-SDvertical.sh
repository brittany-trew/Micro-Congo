#!/bin/bash
#SBATCH -J MBHt       # job name
#SBATCH -c 1                       # number of cores
#SBATCH -t 0-23:00                 # runtime
#SBATCH -p davies,sapphire,shared  # partitions
#SBATCH --mem=190G              # memory in MB
#SBATCH --array=8-12              # no jobs: x years × 25 tiles
#SBATCH -o R/outputs/MBHet_%A.out        # ONE stdout file for whole array
#SBATCH -e R/outputs/MBHet_%A.err        # ONE stderr file for whole array

# Paths
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_19
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_19.sif"

# Batch offset
export ARRAY_OFFSET=0

# Create a unique temporary directory for this array task
scratch_lab="/n/netscratch/davies_lab"
scratch_root="${scratch_lab}/Lab/btrew/tmp"
mkdir -p "$scratch_root" || exit 1

job_tmp=$(mktemp -d \
    "${scratch_root}/job_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_XXXXXX") \
    || exit 1

export TMPDIR="$job_tmp"

# Clean up only this task's temporary directory on exit
trap 'rm -rf -- "$job_tmp"' EXIT

echo "Temporary directory: $TMPDIR"

singularity exec \
  --bind "${my_packages}:/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2" \
  --bind "${job_tmp}:${job_tmp}" \
  --env "TMPDIR=${job_tmp}" \
  --env R_LIBS_USER=/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2 \
  "$rstudio_singularity_image" \
  Rscript --vanilla scripts/mondobai/ColSDVertical.R \
  "$((ARRAY_OFFSET + SLURM_ARRAY_TASK_ID))"