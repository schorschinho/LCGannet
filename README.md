# Osprey
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/schorschinho/osprey)](https://github.com/schorschinho/osprey/releases)
[![GitHub Release Date](https://img.shields.io/github/release-date/schorschinho/osprey)](https://github.com/schorschinho/osprey/releases)
[![Build Status](https://dev.azure.com/OspreyCI/Osprey/_apis/build/status/develop-pull?branchName=develop)](https://dev.azure.com/OspreyCI/Osprey/_build/latest?definitionId=1&branchName=develop)
[![GitHub commit activity](https://img.shields.io/github/commit-activity/m/schorschinho/osprey?foo=bar)](https://github.com/schorschinho/osprey/commits/develop)
[![GitHub last commit](https://img.shields.io/github/last-commit/schorschinho/osprey)](https://github.com/schorschinho/osprey/commits/develop)
[![Website](https://img.shields.io/website?down_color=lightgrey&down_message=offline&up_color=green&up_message=online&url=https%3A%2F%2Fschorschinho.github.io%2Fosprey)](https://schorschinho.github.io/osprey)
[![License](https://img.shields.io/github/license/schorschinho/osprey)](https://github.com/schorschinho/osprey/blob/develop/LICENSE.md)
[![DOI](https://zenodo.org/badge/197683508.svg)](https://doi.org/10.5281/zenodo.14226197)
[![GitHub Total Downloads](https://img.shields.io/github/downloads/schorschinho/osprey/total)](https://github.com/schorschinho/osprey/releases)

<img src="graphics/osprey.png" alt="Osprey" width="200"/>

[Osprey Documentation](https://schorschinho.github.io/osprey)

Osprey is an all-in-one software suite for state-of-the art processing and
quantitative analysis of in-vivo magnetic resonance spectroscopy (MRS) data.

### Features
- 1-file job definition system for reproducible data analysis
- Automated recognition of input file format and sequence origin
- Fully-automated loading and pre-processing pipeline for optimal SNR, linewidth, phasing, and alignment
- Integrated linear-combination modeling module
- Full density-matrix simulated basis sets (using real pulses and sequence timings; including effects of localization and spectral editing) for many metabolites and common sequence implementations
- Functions to create custom basis sets and import basis sets from LCModel or Tarquin
- Integrated voxel co-registration and segmentation module (requires SPM12)
- Quantification based on tissue fractions and (customizable) metabolite/tissue water relaxation times
- GUI to display data, quality assessment, and quantitative results at each step of the analysis
- Rich 'Overview' GUI panel for batched datasets to visualize distributions of metabolite estimates and mean +/- SD spectra

### Supported methods
- Conventional MRS (STEAM, PRESS, semi-LASER, LASER)
- MEGA editing
- Hadamard-encoded editing (HERMES, HERCULES)
- Dual voxel (PRIAM)

### Supported sequence implementations
- Philips (Philips product sequences; Johns Hopkins patches)
- Siemens (Siemens product and WIP sequences; Johns Hopkins patches; CMRR sequences; Jamie Near sequence)
- GE (GE product sequences; Ralph Noeske sequence)

### Supported file formats
- Philips: SDAT/SPAR, DATA/LIST, SIN/LAB/RAW (coming soon)
- Siemens: TWIX/DAT, RDA (single- & multi-file), DICOM (DCM/IMA, single- & multi-file)
- GE: P

## Getting started

### Prerequisites

Osprey requires [MATLAB](https://www.mathworks.com/products/matlab.html) and
has been tested on version 2017a and newer. The following toolboxes are
required for full functionality:

- Optimization
- Statistics and Machine Learning

### Installation

Download the latest **Osprey** code from its [GitHub
repository](https://github.com/schorschinho/osprey), then extract and add the
entire folder (with subfolders) to your MATLAB path. Make sure to regularly
check for updates, as we frequently commit new features, bug fixes, and improved
functions.

To perform voxel co-registration and tissue segmentation, download **SPM12**
[from the UCL website](http://www.fil.ion.ucl.ac.uk/spm/software/spm12/), then
extract and add to your MATLAB path. If you run an Apple Silicon processor 
(M1 and later), please download the [SPM development version from GitHub](https://github.com/spm/spm).

If you want to use the `Osprey` Graphical User Interface (GUI),
please download the following toolboxes from the MATLAB File Exchange:

- [GUI Layout Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)
      (David Sampson)

- [Widget Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/66235-widgets-toolbox)
      (Robin Jackey)

Download both toolboxes in the MATLAB toolbox format (.mltbx). You can
double-click to install. MATLAB will automatically add the toolboxes to its
path.

Make sure to remove FID-A and Gannet from your MATLAB path.

## Contact, Feedback, Suggestions

For any sort of questions, feedback, suggestions, or critique, please visit the [Osprey support forum](https://forum.mrshub.org/c/mrs-software/osprey/10) on the [MRSHub](https://www.mrshub.org).

We also welcome your direct contributions to Osprey here in the GitHub repository.

## Developers

- [Georg Oeltzschner](mailto:goeltzs1@jhu.edu)
- [Helge J. Zöllner](mailto:hzoelln2@jhu.edu)
- [Richard A. E. Edden](mailto:raee2@jhu.edu)

Should you publish material that made use of Osprey, please cite the following publication:

[G Oeltzschner, HJ Zöllner, SCN Hui, M Mikkelsen, MG Saleh, S Tapper, RAE Edden. Osprey: Open-Source Processing, Reconstruction  & Estimation of Magnetic Resonance Spectroscopy Data. J Neurosci Meth 343:108827 (2020).](https://doi.org/10.1016/j.jneumeth.2020.108827)

Please also refer to the version you have used during your analysis using the Zenodo DOI:

[![DOI](https://zenodo.org/badge/197683508.svg)](https://doi.org/10.5281/zenodo.14226197)

## Acknowledgements

This work has been supported by NIH grants R01 EB016089, P41 EB15909, P41 EB031771, R01 EB023963, and K99 AG062230.

We wish to thank collaborators and partners for providing LCModel basis sets and control files. If you use these resources for your analysis of the following data types, please mention the respective individuals in your acknowledgements:
- Siemens 7T STEAM (TE = 5 ms): Dr. Dinesh Deelchand (University of Minnesota)
- Siemens 3T and 7T SPECIAL (TE = 8.5/9 ms): Dr. Ariane Fillmer (PTB Berlin)

We also wish to thank the following individuals for their contributions to the
development of Osprey and shared processing code:

- Jamie Near (McGill University, Montreal)
- Ralph Noeske (GE Healthcare, Berlin)
- Peter Barker (Johns Hopkins University, Baltimore, MD)
- Robin de Graaf (Yale School of Medicine, New Haven, CT)
- Philipp Ehses (German Center for Neurodegenerative Diseases, Bonn)
- Wouter Potters (UMC Amsterdam)
- Xiangrui Li (Ohio State University, Columbus, OH)
- Peter Van Schuerbeek (UZ Brussel)

We are particularly grateful for the incredible [raincloud plot tools](https://github.com/RainCloudPlots/RainCloudPlots) developed by Micah Allen, Davide Poggiali, Kirstie Whitaker, Tom Rhys Marshall, and Rogier Kievit. Should you make use of the OspreyOverview raincloud plots, please consider citing their original publications:

>    * Allen M, Poggiali D, Whitaker K et al. Raincloud plots: a multi-platform tool for robust data
> visualization [version 1; peer review: 2 approved].
> Wellcome Open Res 2019, 4:63. DOI: 10.12688/wellcomeopenres.15191.1
>    * Allen M, Poggiali D, Whitaker K, Marshall TR, Kievit R. (2018) RainCloudPlots tutorials and codebase (Version v1.1). Zenodo. http://doi.org/10.5281/zenodo.3368186
