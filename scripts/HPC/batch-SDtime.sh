#!/bin/bash
#SBATCH -J MBtime       # job name
#SBATCH -c 1                       # number of cores
#SBATCH -t 0-23:00                 # runtime
#SBATCH -p davies,sapphire,shared,serial_requeue  # partitions
#SBATCH --mem=200G               # memory in MB
#SBATCH --array=1-55              # no jobs: x years × 25 tiles
#SBATCH -o R/outputs/MBtime_%A.out        # ONE stdout file for whole array
#SBATCH -e R/outputs/MBtime_%A.err        # ONE stderr file for whole array

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
  Rscript --vanilla scripts/mondobai/ColSDTime.R $((ARRAY_OFFSET + SLURM_ARRAY_TASK_ID))
