# Supporting Code and Data for BGCMIP mansucript

This repository holds the supporting code for the following mansucript:

Kearney, K.A., Cheng, W, and Hermann, A.J. A biogeochemical model intercomparison for the eastern Bering Sea shelf.  Submitted for review to JGR: Oceans.

The contents include the scripts used to compile several variants of the ROMS model for the Bering10K domain, parameterize and run each simulation, analyze output, and create figures for the paper.

Due to storage space limitations, output netCDF files (37 TB) could not be included here.  Instead, we preserve the necessary scripts to reproduce this output, as well as selected post-processed datasets required to reproduce the figures in our paper.

## Folder locations

In its original setup, the folders included in this repository were located on the UW HPC mox.hyak.uw.edu machine.  During the final analysis phase, all data was then moved to klone.hyak.uw.edu.  Residual path names reflecting these locations may be found in the following code.  This table presents the equivalent path names across each machine.  If attempting to reproduce the simulations or analysis in this paper, machine-specific references to various folders will need to be updated accordingly.

| Folder | mox | klone |
| ------ | --- | ----- |
| ``simulations/`` | ``/gscratch/bumblereem/kearney/BGC_hindcasts_workdir/`` | ``/gscratch/cicoes/GR011377_bgcmip/BGC_hindcasts_workdir/`` |
| ``ROMS_Data/`` | ``/gscratch/bumblereem/ROMS_Data/`` | ``/gscratch/cicoes/GR011377_bgcmip/ROMS_Data/`` |
| ``analysis/`` | ``/gscratch/bumblereem/kearney/BGC_hindcasts_workdir/analysis/`` | ``/gscratch/cicoes/GR011377_bgcmip/BGC_hindcasts_workdir/analysis/``


## ROMS\_Data/

ROMS input data.  Due to storage considerations, the full dataset (507 GB) is unable to be made fully available through this platform.  However, we include .cdl printouts corresponding to the original .nc files; these can be perused to see the structure of input required, including variable names and dimensions.  Data available upon request: kelly dot kearney at noaa dot gov.

## simulations/

This folder is the ROMS application folder for all simulations run for this manuscript.  It includes the following scripts:

- ``compilebgc.sh``: Calls build script to compile ROMS in the different variants.
- ``build_roms.sh``: Modified version of standard ROMS build script that sets the proper CPP flags for the 3 BGC model variants.
- ``bgc_hindcasts.py``: python scripts that configures standard input files for all simulations, call the ROMS executable, and restart simulations as necessary to complete all simulations.
- ``bgc_hindcasts.sub``: SLURM submission script that shows the specific calls to ``bgc_hindcasts.py`` that were used to complete all simulations.

The ``bgcmip_*/In/`` subfolders hold the ROMS input parameter files for each simulation: physical ocean (ocean), biological (bpar), ice (ipar), and station (spos).  The ``bgcmip_*/Log/`` subfolders include standard error and standard output for each simulation.  Note that time-stepping in the standard output files has been trimmed for storage space consideration, but all relevant configuration information has been retained.

## analysis/

This folder includes the primary Matlab analysis script (bgcmip_paper2.m) used to calculate metrics and generate figures for this paper.  It also includes all code dependencies of that script and the temporary output datasets (in .mat format) required to reproduce those figures.

## roms-kate-ice/

This submodule is a snapshot of the ROMS source code compiled for this application.  It reflects the [beringnpz/roms](https://github.com/beringnpz/roms) code, main branch, commit c861bc98bfb071906809f01aa974886cba1745c1.

## bering-Apps/

This submodule is a snapshot of the Bering Sea (and related) ROMS application and sub-application data (i.e. files associated with the domain, biological modules, etc.) that are the basis for the input files in the ``simulations/`` folder.  It reflects the [beringnpz/bering-Apps](https://github.com/beringnpz/bering-Apps) code, master branch, commit commit f959da05007b8c04cced507802bdbdf9cbbbade4.

## romscom/

This submodule is a snapshot of the ROMS Communication Toolbox code used to configure and run the simulations.  It reflects the [beringnpz/romscom](https://github.com/beringnpz/romscom) code, master branch, commit 492a3d634ec0d278444bbd8a55778fae5599d42d.




