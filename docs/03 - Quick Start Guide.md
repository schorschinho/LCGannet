## Quick Start Guide

This chapter is intended to serve as a quick rundown of the `Osprey` workflow,
as well as syntax cheat sheet for more experienced users. The full documentation for
each of the following commands can be found in their respective sections and
will be released shortly.

Please familiarize yourself with the basic `Osprey` commands. The number of
commands is deliberately kept to a minimum, as is the number of possible input
arguments to these commands.

### OspreyJob

At the beginning of every analysis, you initialize a fresh `MRSCont` data
container with a job (defined in the file `jobFile`) using the following
command:

```
[MRSCont] = OspreyJob(jobFile);
```

`OspreyJob` will add all necessary information required for data analysis to
the data container, and perform a few basic sanity checks.

A valid Osprey job file contains:
- a list of MRS (and, optionally, structural imaging) data files to be loaded
- basic information on the MRS sequence used
- settings for data modeling
- an output folder to store the results and exported files

Osprey job files are described in `.m` format. Please refer to
the `OspreyJob` chapter for detailed instructions on how to create a job. A good way to understand the structure of a job file is to study the example jobs in the `Osprey/exampledata` directory, and to create your own job starting from there.

Before you create your own job, consider the recommendations on how to organize
your data in the previous chapter.

### OspreyLoad

After loading a job file, you can start the data loading by invoking the
function

```
[MRSCont] = OspreyLoad(MRSCont);
```

`OspreyLoad` will automatically determine the file format, load the raw data,
and, if multi-array receiver coils have been used for acquisition, perform the
receiver-channel combination.

Please refer to the `OspreyLoad` chapter for details on the supported file
formats and the coil combination techniques.

### OspreyProcess

`OspreyLoad` only loads the raw data and organizes them in a uniform way, in
order to streamline further processing. The next command is a cornerstone of the
entire Osprey pipeline, as it will initiate a series of automatic processing
steps:

```
[MRSCont] = OspreyProcess(MRSCont);
```

`OspreyProcess` performs all state-of-the-art MRS post-processing techniques
to ensure optimal spectral quality, including eddy-current correction,
frequency-and-phase correction, removal of residual water, baseline correction,
and automated phasing and frequency referencing. In addition, several basic
quality metrics (linewidth, SNR) are determined.

`OspreyProcess` is also responsible for handling all aspects of data
processing that are specific to *J*-difference-edited ("spectral editing")
experiments, such as MEGA, HERMES, and HERCULES. These additional steps include
alignment of sub-spectra to reduce subtraction artefacts, as well as the
calculation of (Hadamard-encoded) difference spectra.

Finally, `OspreyProcess` can write the processed spectra into files formats that are readable by spectral fitting software packages like LCModel, TARQUIN, or jMRUI:
- To write LCModel-readable `.RAW` files, you need to set a switch `opts.saveLCM = 1` in your job file. You can then load the `.RAW` file via the `Other` option in the LCModel file type selection menu.
- To write jMRUI/TARQUIN-readable `.TXT` files, you need to set a switch `opts.saveJMRUI = 1` in your job file.
- Please be aware that LCModel may prompt you to enter the number of FID data points (a positive integer number, e.g. 2048), dwell time (in seconds, e.g. 0.0005), and the static magnetic field strength in Hz (eg 123.26 for a Siemens 2.89 T magnet). In TARQUIN, you will probably need to enter the echo time in seconds (e.g. 0.03 for a 30-ms acquisition). You can extract this information from the Osprey `MRSCont` data container, or ...
- ... by setting the switch `opts.saveVendor = 1` in your job file, `Osprey` will create single vendor-specific files, i.e. `.SDAT/.SPAR` files for Philips data, and `.RDA` files for Siemens data, regardless of the raw data format. This function has been introduced to sidestep the above-mentioned problems with `.RAW` and `.TXT` files. LCModel, TARQUIN and jMRUI accept the `.SDAT/.SPAR` and `.RDA` file formats, and should be able to read the header information correctly.

All exported third-party format files are stored in sub-directories of the output folder specified in the job file.

Please refer to the `OspreyProcess` chapter for details on the implementation
of the various post-processing steps.

### OspreyCoreg

You can provide a T1-weighted structural image in NIfTI format (`*.nii`) in
order to perform voxel co-registration and segmentation. The function
`OspreyCoreg` will extract voxel geometry information from the MRS data
headers, and leverage SPM12 NIfTI handling functions to create an SPM image
volume containing a voxel mask. The voxel masks will be saved in a sub-directory (`/VoxelMasks/`) of the output folder specified in the job file.

The path to the structural image is provided in the job definition. After
running `OspreyLoad`, you can call the co-registration routine using the
command

```
[MRSCont] = OspreyCoreg(MRSCont);
```

Please refer to the `OspreyCoreg` chapter for details.

### OspreySeg

`OspreySeg` calls the SPM12 tissue segmentation algorithm ("New
Segment") to segment the T1-weighted structural image into grey matter (GM),
white matter (WM), and cerebro-spinal fluid (CSF). The relative volume fractions
for each tissue class is subsequently used to perform water-scaled concentration
estimation, including CSF correction and tissue-specific relaxation correction.

Please refer to the `OspreySeg` chapter for details.

### OspreyFit

`OspreyFit` will be released shortly.

This function will perform linear-combination modeling of the processed spectra.
The Osprey fitting algorithm is using metabolite basis functions for the
specific MRS sequence, which are calculated from 2-D spatially resolved full
density-matrix simulations, using the real sequence timings and pulse waveforms
wherever possible.

The linear-combination modeling performed by `OspreyFit` is similar to the
fitting carried out by the algorithms in LCModel, TARQUIN, or QUEST/AQSES in
jMRUI. In all brevity, the algorithm performs a non-linear least-squares
optimization (Levenberg-Marquardt) of several amplitude, phase, lineshape, and
baseline parameters in order to minimize the difference between the modeled
spectra and the data.

`OspreyFit` will save all relevant modeling parameters into the `MRSCont` data
container. These values will be subsequently used in `OspreyQuant` to produce
quantitative estimates for all metabolites specified in the basis set.

Please refer to the `OspreyFit` chapter for details on the algorithmic
implementation, basis set creation, and available modeling settings.

### OspreyQuant

`OspreyQuant` will be released shortly.
