# Recommendations for Quantifying Rhythmic and Arrhythmic Components of Brain Activity

## Background

Brain activity is composed of both arrhythmic and rhythmic signal elements. Growing interest in the field aims to quantify the independent contribution of each signal element to behaviour and disease. Yet, little work explores whether current spectral modelling approaches can accurately recover the independent contribution of these neural signal elements. 

Here we test how three different methods for quantifying brain rhythms impact the interpretation of findings. 

We demonstrate that spectral detrending methods (i.e., subtracting the arrhythmic fit from the spectrum in either log-log or linear space) introduce spurious relationships between spectral model parameters. This ultimately challenges the robustness, reproducibility, and interpretability of findings. 


## Description of the project

Here are the different steps we did :

- We simulate neural time series data using the NeuroDSP toolbox to generate ground truth data to parameterize (Simulate_neural_timeseries.ipynb)
- We model the simulated neural time series data using the ms-specparam function as implemented in brainstorm and organize the outputs into a csv file for statistical analyses (Specparam_timeseries_and_organize_outputs.m)
- We evaluate how the choice of methodology for quantifying brain rhythms impacts the interpretation of the results (Statistical_analysis_and_plotting.Rmd)

## Manuscript and Citation

This work on [bioRxiv](XXXXXXX). Please cite da Silva Castanheira et al., 2025. If you have any questions please contact the authors of the paper.


