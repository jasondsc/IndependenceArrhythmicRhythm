# Recommendations for Quantifying Rhythmic and Arrhythmic Components of Brain Activity

## Background

Brain activity is composed of both arrhythmic and rhythmic signal elements. Growing interest in the field aims to quantify the independent contribution of each signal element to behaviour and disease. Yet, little work explores whether current spectral modelling approaches can accurately recover the independent contribution of these neural signal elements. 

Here we test how three different methods for quantifying brain rhythms impact the interpretation of findings. 

We demonstrate that spectral detrending methods (i.e., subtracting the arrhythmic fit from the spectrum in either log-log or linear space) introduce spurious relationships between spectral model parameters. This ultimately challenges the robustness, reproducibility, and interpretability of findings. 

## Data

Data generated for this project (simulations) are aviable from the corresponding author upon reasonable request. Some examples of the simulated time series are on [OSF](https://osf.io/nx693/overview) while some of the organized outputs are included in the GitHub (see Folder). Empirical data are accessible from the CamCAN website.

Note, we include some examples of the simulated time series, and the organized outputs after PSD and spectral modelling.

## Description of the project

Here are the different steps we did :

- We simulate neural time series data using the NeuroDSP toolbox to generate ground truth data to parameterize (Simulate_neural_timeseries.ipynb) This relies on the sim_peak_oscillation function
- We model the simulated neural time series data using the ms-specparam function as implemented in brainstorm and organize the outputs into a csv file for statistical analyses (Specparam_timeseries_and_organize_outputs.m). Within this function we compute the corrected linear and log spectra
- We evaluate how the choice of methodology for quantifying brain rhythms impacts the interpretation of the results (Statistical_analysis_and_plotting.Rmd). Here we use the corrected linear, log, or modelled rhythmic spectra to compute alpha power and run statistical tests.
- We then test whether the divergence between methods is similarly observed in an empirical dataset (CamCan_restingstateAlpha.Rmd)

## Manuscript and Citation

This work is on [bioRxiv](https://doi.org/10.1101/2025.09.24.678322). Please cite da Silva Castanheira et al., 2025. If you have any questions please contact the authors of the paper.


