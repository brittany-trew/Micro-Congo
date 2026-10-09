#!/bin/bash
#SBATCH -J imbS       # job name
#SBATCH -c 1                       # number of cores
#SBATCH -t 0-23:00                 # runtime
#SBATCH -p davies,sapphire,shared  # partitions
#SBATCH --mem=15000               # memory in MB
#SBATCH --array=1-4              # no jobs: x years × 25 tiles
#SBATCH -o R/outputs/IMHetS_%A.out        # ONE stdout file for whole array
#SBATCH -e R/outputs/IMHetS_%A.err        # ONE stderr file for whole array

# Paths
my_packages=${HOME}/R/ifxrstudio/RELEASE_3_19
rstudio_singularity_image="/n/singularity_images/informatics/ifxrstudio/ifxrstudio:RELEASE_3_19.sif"

# Batch offset
export ARRAY_OFFSET=0

# Run historical model
singularity exec \
  --bind ${my_packages}:/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2 \
  --env R_LIBS_USER=/home/rstudio/R/x86_64-pc-linux-gnu-library/4.2 \
  $rstudio_singularity_image \
  Rscript --vanilla scripts/imbalanga/ColSDVerticalStrata.R $((ARRAY_OFFSET + SLURM_ARRAY_TASK_ID))
