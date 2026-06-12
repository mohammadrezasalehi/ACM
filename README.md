# ACM
A framework for mapping amplitude-dependent modulatory interactions and higher-order functional connectivity in resting-state fMRI.


# Amplitude-Conditional Modulation of Functional Connectivity

This repository contains MATLAB codes of the manuscript:

"Amplitude-Conditioned Brain Dynamics Reveal a Global Inhibition and Targeted Facilitation of Network Modulation"

## Overview

Traditional functional connectivity characterizes pairwise interactions between brain regions. This project introduces an amplitude-conditioned framework for quantifying how spontaneous fluctuations in a third region modulate connectivity among all other brain regions.

The framework generates modulation spectra and extracts four interpretable modulatory traits:

* Modulatory Sensitivity (Variance)
* Directional Control (Slope)
* Volatility (Zero-Crossing Rate)
* Nonlinear Complexity (Polynomial Degree)

The repository also includes implementations of surrogate-based validation, reliability analysis, and behavioral prediction.

## Repository Structure

src/03_ACM
Core computational functions of the manuscript

data/inter_parcellated/Schaefer200_Kong17/REST1_LR
Example inputs data

## Requirements

MATLAB R2018b or later

Required toolboxes:

* Statistics and Machine Learning Toolbox
* Signal Processing Toolbox
*cifti-matlab-master
*cmap-master
*workbench

## Data Availability

Resting-state fMRI data were obtained from the Human Connectome Project (HCP).

https://db.humanconnectome.org

Users must obtain HCP data directly through the official repository.

## Citation

If you use this code, please cite.
