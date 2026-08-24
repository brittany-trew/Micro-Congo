# Micro-Congo

Code associated with analyses of vertical forest structure and microclimate
in Congolese tropical forests.

Paper Title: Vertical microclimate diversity shapes thermal refugia in Congolese rainforests.
Authors: Brittany T. Trew, Rebecca Senior, Evan G. Hockridge, Tom R. Bishop, Charlene Janion-Scheepers, Gwili Gibbon, Guillaume Baltus, Andrew B. Davies

Code Author: Brittany Trew, contact brittany.trew@fas.harvard.edu.

## Overview

This repository contains R scripts used to process LiDAR-derived forest
structure data, model below-canopy microclimates, and analyse vertical
microclimatic variation across study sites for LiDAR sites in Odzala-Kokoua NP, Congo.

## Repository structure

- `scripts/` — R scripts used for data processing and analysis
- `scripts/lidar/` — LiDAR processing
- `scripts/microclimate/` — microclimate modelling
- `scripts/analysis/` — statistical analyses
- `scripts/figures/` — figure generation

## Data

Raw and processed datasets are not included in this repository because of
their size.

## Requirements

Analyses were conducted in R and Google Earth Engine.

## Workflow

1. Process LiDAR data and derive forest structural metrics.
2. Generate below-canopy temperature estimates.
3. Calculate vertical microclimate metrics.
4. Analyse relationships between forest structure and microclimate.
5. Generate manuscript figures.

## Citation

Citation information will be added following publication.