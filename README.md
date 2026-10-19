## Requirements

The workflow uses R and Google Earth Engine. Microclimate modelling and
subsequent processing were run on a high-performance computing (HPC)
cluster using the SLURM scheduler.

The full workflow is computationally intensive and was designed for HPC
execution. SLURM submission scripts are provided in `scripts/HPC/`.
Before use, adapt file paths, account names, partitions and resource
requests to your computing environment.

## Workflow

Run R scripts from the repository root so that relative paths resolve
correctly. First, edit `scripts/parameters.R` to specify file paths and
site settings.

### 1. Process LiDAR data

Run the following R scripts in order:

- `scripts/ProcessLiDAR.R`
- `scripts/StructuralClasses.R`
- `scripts/StructuralStrata.R`

### 2. Prepare microclimate model inputs

Download land-cover and PAI data from Google Earth Engine, eg:

- `scripts/gee/DownloadSentinel.js`
- `scripts/gee/DownloadMODIS.js`

Run the site-specific land-cover processing script, eg: 

- `scripts/land_cover/process_landcover_Imbalanga.R`

followed by:

- `scripts/CalcPAI.R`
- `scripts/GroundMetrics.R`
- `scripts/processMacroclimate.R`

### 3. Run the microclimate model and process outputs

Submit the shell scripts from a terminal on the HPC cluster using
`sbatch`. For example:

```bash
sbatch scripts/HPC/batch-imbalanga.sh
```

The submission stages are listed below. Allow each stage to complete
before submitting any stage that requires its outputs.

1. `batch-imbalanga.sh` — run the model for Imbalanga; use the
   corresponding submission scripts for other sites.
2. `batch-output.sh` — process model outputs into daily summaries.
3. `batch-dailymax.sh`
4. `batch-SDvertical.sh`
5. `batch-SDtime.sh`
6. `batch-SDverticalZones.sh`
7. `batch-SDtemporalZones.sh`
8. `batch-yearVars.sh`

All submission scripts are located in `scripts/HPC/`.

### 4. Analyse results

- `scripts/strata_means.R`
- `scripts/scale_variables.R`
- `scripts/tmax_warming.R`
- `scripts/models_PAI.R`
- `scripts/modelsCD.R`

### 5. Generate Main Figures

- `scripts/Fig1.R`
- `scripts/Fig2.R`
- `scripts/Fig3.R`



Run the statistical analysis and figure-generation scripts using the
outputs produced above.