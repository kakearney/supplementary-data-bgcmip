# Supporting Code and Data for BGCMIP mansucript

## Summary

This repository holds the supporting code for the following mansucript:

Kearney, K.A., Cheng, W, and Hermann, A.J. A biogeochemical model intercomparison for the eastern Bering Sea shelf.  Submitted for review to JGR: Oceans.

The contents include the scripts used to compile several variants of the ROMS model for the Bering10K domain, parameterize and run each simulation, analyze output, and create figures for the paper.

Due to storage space limitations, output netCDF files (37 TB) could not be included here.  Instead, we preserve the necessary scripts to reproduce this output, as well as selected post-processed datasets required to reproduce the figures in our paper.

In its original setup, the folders included in this repository were located on the UW HPC mox.hyak.uw.edu machine.  During the final analysis phase, all data was then moved to klone.hyak.uw.edu.  Residual path names reflecting these locations may be found in the following code.  This table presents the equivalent path names across each machine.  If attempting to reproduce the simulations or analysis in this paper, machine-specific references to various folders will need to be updated accordingly.

| Folder | mox | klone |
| ------ | --- | ----- |
| ``simulations/`` | ``/gscratch/bumblereem/kearney/BGC_hindcasts_workdir/`` | ``/gscratch/cicoes/GR011377_bgcmip/BGC_hindcasts_workdir/`` |
| ``ROMS_Datasets/`` | ``/gscratch/bumblereem/ROMS_Datasets/`` | ``/gscratch/cicoes/GR011377_bgcmip/ROMS_Datasets/`` |
| ``analysis/`` | ``/gscratch/bumblereem/kearney/BGC_hindcasts_workdir/analysis/`` | ``/gscratch/cicoes/GR011377_bgcmip/BGC_hindcasts_workdir/analysis/``
|``bering-Apps/``|``/gscratch/bumblereem/kearney/bering-Apps/``| ``/gscratch/cicoes/GR011377_bgcmip/bering-Apps``|
|``roms-kate-ice/``|``/gscratch/bumblereem/roms-kate-ice``| ``/gscratch/cicoes/GR011377_bgcmip/roms-kate-ice``|

## Contents

### ROMS\_Datasets/

ROMS input data.  Due to storage considerations, the full dataset (507 GB) is unable to be made available through this platform.  However, we include .cdl printouts corresponding to the original .nc files; these can be perused to see the structure of input required, including variable names, units, and dimensions.  Symbolic links to the hard copy of the files are also included (primarily to remind me where that data copy is located).  Data available upon request from kelly dot kearney at noaa dot gov.

### simulations/

This folder is the ROMS application folder for all simulations run for this manuscript.  It includes the following scripts:

- ``compilebgc.sh``: Calls build script to compile ROMS in the different variants.
- ``build_roms.sh``: Modified version of standard ROMS build script that sets the proper CPP flags for the 3 BGC model variants.
- ``bgc_hindcasts.py``: python scripts that configures standard input files for all simulations, call the ROMS executable, and restart simulations as necessary to complete all simulations.
- ``bgc_hindcasts.sub``: SLURM submission script that shows the specific calls to ``bgc_hindcasts.py`` that were used to complete all simulations.

The ``bgcmip_*/In/`` subfolders hold the ROMS input parameter files for each simulation: physical ocean (ocean), biological (bpar), ice (ipar), and station (spos).  The ``bgcmip_*/Log/`` subfolders include standard error and standard output for each simulation.  Note that time-stepping in the standard output files has been trimmed for storage space consideration, but all relevant configuration information has been retained.

### analysis/

This folder includes the primary Matlab analysis script (bgcmip_paper_final.m) used to calculate metrics and generate figures for this paper.  It also includes some additional .m files called directly by that script and the temporary output datasets (in .mat format) required to reproduce those figures.

Note that if you wish to run this code locally, this script uses large number of third party functions, primarily small plotting or data analysis utilities.  These can be downloaded from GitHub
( [barareaneg    ](https://github.com/kakearney/barareaneg-pkg  )
, [Bering10KInput](https://github.com/beringnpz/Bering10KInput  )
, [labelaxes     ](https://github.com/kakearney/labelaxes-pkg   )
, [legendflex    ](https://github.com/kakearney/legendflex-pkg  )
, [mergeaxes     ](https://github.com/kakearney/mergeaxes-pkg   )
, [minmax        ](https://github.com/kakearney/minmax-pkg      )
, [offsetaxis    ](https://github.com/kakearney/offsetaxis-pkg  )
, [pcolorpad     ](https://github.com/kakearney/pcolorpad-pkg   )
, [plotboxpos    ](https://github.com/kakearney/plotboxpos-pkg  )
, [plotgrid      ](https://github.com/kakearney/plotgrid-pkg    )
, [roms          ](https://github.com/kakearney/roms            )
, [shapeprjread  ](https://github.com/kakearney/shapeprjread-pkg)
, [singlepatch   ](https://github.com/kakearney/singlepatch-pkg )
, [suplabel      ](https://github.com/kakearney/suplabel-pkg    )
, [uniquecell    ](https://github.com/kakearney/uniquecell-pkg  )
, [wraptext      ](https://github.com/kakearney/wraptext-pkg    )
) or from the MatlabCentral File Exchange
( [Climate Data Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/70338 )
, [export_fig          ](https://www.mathworks.com/matlabcentral/fileexchange/23629 )
, [interparc           ](https://www.mathworks.com/matlabcentral/fileexchange/34874 )
, [subaxis             ](https://www.mathworks.com/matlabcentral/fileexchange/3696  )
, [wprctile            ](https://www.mathworks.com/matlabcentral/fileexchange/16920 )
, [RunLength           ](https://www.mathworks.com/matlabcentral/fileexchange/41813 )
, [crameri             ](https://www.mathworks.com/matlabcentral/fileexchange/68546 )
, [hex2rgb             ](https://www.mathworks.com/matlabcentral/fileexchange/46289 )
, [treemap             ](https://www.mathworks.com/matlabcentral/fileexchange/17192 )
, [yaml                ](https://www.mathworks.com/matlabcentral/fileexchange/106765)
) and must be added to your Matlab path for full functionality.


### roms-kate-ice/

This submodule is a snapshot of the ROMS source code compiled for this application.  It reflects the [beringnpz/roms](https://github.com/beringnpz/roms) code, main branch, commit c861bc98bfb071906809f01aa974886cba1745c1.  This is a fork of Kate Hedstrom's version of ROMS with sea ice.

### bering-Apps/

This submodule is a snapshot of the Bering Sea (and related) ROMS application and sub-application data (i.e. files associated with the domain, biological modules, etc.) that are the basis for the input files in the ``simulations/`` folder.  It reflects the [beringnpz/bering-Apps](https://github.com/beringnpz/bering-Apps) code, master branch, commit commit f959da05007b8c04cced507802bdbdf9cbbbade4.

### romscom/

This submodule is a snapshot of the ROMS Communication Toolbox code used to configure and run the simulations.  It reflects the [beringnpz/romscom](https://github.com/beringnpz/romscom) code, master branch, commit 492a3d634ec0d278444bbd8a55778fae5599d42d.




