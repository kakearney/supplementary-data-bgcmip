#!/bin/bash
#
# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Copyright (c) 2002-2020 The ROMS/TOMS Group                           :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#                                                                       :::
# ROMS/TOMS Compiling BASH Script                                       :::
#                                                                       :::
# Script to compile an user application where the application-specific  :::
# files are kept separate from the ROMS source code.                    :::
#                                                                       :::
# Q: How/why does this script work?                                     :::
#                                                                       :::
# A: The ROMS makefile configures user-defined options with a set of    :::
#    flags such as ROMS_APPLICATION. Browse the makefile to see these.  :::
#    If an option in the makefile uses the syntax ?= in setting the     :::
#    default, this means that make will check whether an environment    :::
#    variable by that name is set in the shell that calls make. If so   :::
#    the environment variable value overrides the default (and the      :::
#    user need not maintain separate makefiles, or frequently edit      :::
#    the makefile, to run separate applications).                       :::
#                                                                       :::
# Usage:                                                                :::
#                                                                       :::
#    ./build_roms.bash [options]                                        :::
#                                                                       :::
# Options:                                                              :::
#                                                                       :::
#    -j [N]       Compile in parallel using N CPUs                      :::
#                  omit argument for all available CPUs                 :::
#                                                                       :::
#    -p macro     Prints any Makefile macro value. For example,         :::
#                                                                       :::
#                  build.bash -p FFLAGS                                 :::
#                                                                       :::
#    -noclean     Do not clean already compiled objects                 :::
#                                                                       :::
#    -db          Compile in debug mode (USE_DEBUG)                     :::
#                                                                       :::
#    -variant [x] CPP-combo variant to compile.  Can be:                :::
#                 phys:    no bio, default ice [default]                :::
#                 cobalt:  BIO_COBALT bio, default ice                  :::
#                 bestnpz: BEST_NPZ bio, default ice                    :::
#                 banas:   BIO_BANAS bio, default ice                   :::
#                                                                       :::
#                                                                       :::
# Notice that sometimes the parallel compilation fail to find MPI       :::
# include file "mpif.h".                                                :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

which_MPI=openmpi                             # default, overwritten below

parallel=0
clean=1
dprint=0
compiledebug=0
compileserial=0
compilevar="phys"
compdate=$(date "+%Y%m%d%H%M")


MY_CPP_FLAGS=

while [ $# -gt 0 ]
do
  case "$1" in
    -j )
      shift
      parallel=1
      test=`echo $1 | grep '^[0-9]\+$'`
      if [ "$test" != "" ]; then
        NCPUS="-j $1"
        shift
      else
        NCPUS="-j"
      fi
      ;;

    -p )
      shift
      clean=0
      dprint=1
      debug="print-$1"
      shift
      ;;

    -noclean )
      shift
      clean=0
      ;;
      
    -variant )
      shift
      compilevar="$1"
      shift
      ;;
      
    -db )
      shift
      compiledebug=1
      ;;
      
    -serial )
      shift
      compileserial=1
      ;;
      
    -datestr )
      shift
      compdate="$1"
      shift
      ;;
      
    * )
      echo ""
      echo "$0 : Unknown option [ $1 ]"
      echo ""
      echo "Available Options:"
      echo ""
      echo "-j [N]      Compile in parallel using N CPUs"
      echo "              omit argument for all avaliable CPUs"
      echo ""
      echo "-p macro    Prints any Makefile macro value"
      echo "              For example:  build.bash -p FFLAGS"
      echo ""
      echo "-noclean    Do not clean already compiled objects"
      echo ""
      exit 1
      ;;
  esac
done

# Set the CPP option defining the particular application. This will
# determine the name of the ".h" header file with the application
# CPP definitions.

export   ROMS_APPLICATION=BERING_10K

# Set a local environmental variable to define the path to the directories
# where all this project's files are kept.

export        MY_ROOT_DIR=${PWD}
export     MY_PROJECT_DIR=${MY_PROJECT_DIR:=/gscratch/bumblereem/kearney/testBeringApp}

# The path to the user's local current ROMS source code.
#
# If using svn locally, this would be the user's Working Copy Path (WCPATH).
# Note that one advantage of maintaining your source code locally with svn
# is that when working simultaneously on multiple machines (e.g. a local
# workstation, a local cluster and a remote supercomputer) you can checkout
# the latest release and always get an up-to-date customized source on each
# machine. This script is designed to more easily allow for differing paths
# to the code and inputs on differing machines.

 export       MY_ROMS_SRC=/gscratch/bumblereem/kearney/roms-kate-ice

# Set path of the directory containing makefile configuration (*.mk) files.
# The user has the option to specify a customized version of these files
# in a different directory than the one distributed with the source code,
# ${MY_ROMS_SRC}/Compilers. If this is the case, you need to keep these
# configurations files up-to-date.

 export         COMPILERS=${MY_ROMS_SRC}/Compilers
#export         COMPILERS=${HOME}/Compilers/ROMS

#--------------------------------------------------------------------------
# Set tunable CPP options.
#--------------------------------------------------------------------------
#
# Sometimes it is desirable to activate one or more CPP options to run
# different variants of the same application without modifying its header
# file. If this is the case, specify each options here using the -D syntax.
# Notice also that you need to use shell's quoting syntax to enclose the
# definition.  Both single or double quotes work. For example,
#
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DAVERAGES"
#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -DDEBUGGING"
#
# can be used to write time-averaged fields. Notice that you can have as
# many definitions as you want by appending values.

#export      MY_CPP_FLAGS="${MY_CPP_FLAGS} -D"

#--------------------------------------------------------------------------
# Compiler options.
#--------------------------------------------------------------------------
#
# Other user defined environmental variables. See the ROMS makefile for
# details on other options the user might want to set here. Be sure to
# leave the switches meant to be off set to an empty string or commented
# out. Any string value (including off) will evaluate to TRUE in
# conditional if-statements.

 export           USE_MPI=on            # distributed-memory parallelism
 export        USE_MPIF90=on            # compile with mpif90 script
#export         which_MPI=mpich         # compile with MPICH library
#export         which_MPI=mpich2        # compile with MPICH2 library
#export         which_MPI=mvapich2      # compile with MVAPICH2 library
 export         which_MPI=openmpi       # compile with OpenMPI library

#export        USE_OpenMP=on            # shared-memory parallelism

 export              FORT=ifort
#export              FORT=gfortran
#export              FORT=pgi

#export         USE_DEBUG=on            # use Fortran debugging flags
 export         USE_LARGE=on            # activate 64-bit compilation
#export       USE_NETCDF4=on            # compile with NetCDF-4 library
#export   USE_PARALLEL_IO=on            # Parallel I/O with NetCDF-4/HDF5


#--------------------------------------------------------------------------
# If Earth System Model (ESM) coupling, set location of ESM component
# libraries and modules. The strategy is to compile and link each ESM
# component separately first and ROMS last since it is driving the
# coupled system. Only the ESM components activated are considered and
# the rest ignored.
#--------------------------------------------------------------------------

# export        WRF_SRC_DIR=${HOME}/ocean/repository/WRF
#
# if [ -n "${USE_DEBUG:+1}" ]; then
#   export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_ciceG
#   export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coampsG
#   export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcmG
#   export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wamG
# # export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrfG
#   export      WRF_LIB_DIR=${WRF_SRC_DIR}
# else
#   export     CICE_LIB_DIR=${MY_PROJECT_DIR}/Build_cice
#   export   COAMPS_LIB_DIR=${MY_PROJECT_DIR}/Build_coamps
#   export    REGCM_LIB_DIR=${MY_PROJECT_DIR}/Build_regcm
#   export      WAM_LIB_DIR=${MY_PROJECT_DIR}/Build_wam
#   export      WRF_LIB_DIR=${MY_PROJECT_DIR}/Build_wrf
# # export      WRF_LIB_DIR=${WRF_SRC_DIR}
# fi

#--------------------------------------------------------------------------
# If applicable, use my specified library paths.
#--------------------------------------------------------------------------

 export USE_MY_LIBS=              # use system default library paths
#export USE_MY_LIBS=yes           # use my customized library paths

MY_PATHS=${COMPILERS}/my_build_paths.bash

if [ "${USE_MY_LIBS}" = "yes" ]; then
  source ${MY_PATHS} ${MY_PATHS}
fi

# KAK: For initial testing, I'm just setting the sepcific paths that we used
# in our older build script, and otherwise relying on defaults.

export   PATH=/gscratch/sw/intel-201703/compilers_and_libraries_2017.2.174/linux/mpi/intel64/bin/mpif90:$PATH
export   NETCDF_INCDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/include # netcdf include
export   NETCDF_LIBDIR=/sw/netcdf-fortran+c-4.4.1.1_icc-17/lib     # netcdf lib


#--------------------------------------------------------------------------
# The rest of this script sets the path to the users header file and
# analytical source files, if any. See the templates in User/Functionals.
#--------------------------------------------------------------------------
#
# If applicable, use the MY_ANALYTICAL_DIR directory to place your
# customized biology model header file (like fennel.h, nemuro.h, ecosim.h,
# etc).

 export     MY_HEADER_DIR=${MY_PROJECT_DIR}

 export MY_ANALYTICAL_DIR=${MY_PROJECT_DIR}

# Put the binary to execute in the following directory.

 export            BINDIR=${MY_PROJECT_DIR}

# Put the f90 files in a project specific Build directory to avoid conflict
# with other projects.

if [ -n "${USE_DEBUG:+1}" ]; then
 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_romsG
else
 export       SCRATCH_DIR=${MY_PROJECT_DIR}/Build_roms
fi

# Go to the users source directory to compile. The options set above will
# pick up the application-specific code from the appropriate place.

 cd ${MY_ROMS_SRC}

# Stop if activating both MPI and OpenMP at the same time.

if [ -n "${USE_MPI:+1}" ] && [ -n "${USE_OpenMP:+1}" ]; then
  echo "You cannot activate USE_MPI and USE_OpenMP at the same time!"
  exit 1
fi

#--------------------------------------------------------------------------
# Compile: Physics only
#--------------------------------------------------------------------------


# Turn on debugging if necessary

if [ "$compiledebug" -eq 1 ]; then
  export USE_DEBUG=on
  rbase="romsG"
  debugstr="in debug mode"
else
  export  USE_DEBUG=
  rbase="romsM"
  debugstr=""
fi

# Switch to serial

if [ "$compileserial" -eq 1 ]; then
  export USE_MPI=
  export USE_MPIF90=
  if [ "$compiledebug" -eq 1 ]; then
    rbase="romsG"
  else
    rbase="romsS"
  fi
fi

# Variant flags and details

case "$compilevar" in
  phys )
    longname="physics-only"
    ;;
  cobalt )
    longname="COBALT"
    export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBIO_COBALT"
    ;;
  bestnpz )
    longname="BEST_NPZ"
    export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBEST_NPZ -DCARBON"
    ;;
  banas )
    longname="BIO_BANAS"
    export MY_CPP_FLAGS="${MY_CPP_FLAGS} -DBIO_BANAS"
    ;;
  * )
    echo "Unknown variant"
    echo "  (options: phys, cobalt, bestnpz)"
    echo ""
    exit 1
    ;;
esac


# Scratch directory
  
export  SCRATCH_DIR=${MY_PROJECT_DIR}/Build_${compilevar}_${compdate}

# Remove build directory.

if [ $clean -eq 1 ]; then
  echo "Cleaning..."
  make clean
fi

# Compile (the binary will go to BINDIR set above).

if [ $dprint -eq 1 ]; then
  make $debug
else
  if [ $parallel -eq 1 ]; then
    make $NCPUS
  else
    echo ""
    echo "Compiling $longname variant $debugstr"
    make &> buildroms_log.txt
    if [ $? -ne 0 ]; then
      mv buildroms_log.txt ${SCRATCH_DIR}/buildroms_log.txt
      echo "  Compilation failed: see ${SCRATCH_DIR}/buildroms_log.txt for details"
    else
      mv buildroms_log.txt ${SCRATCH_DIR}/buildroms_log.txt
      mv ${MY_PROJECT_DIR}/${rbase} ${MY_PROJECT_DIR}/${rbase}_${compilevar}_${compdate}
      echo "  Success: ${rbase}_${compilevar}_${compdate} created"
    fi
  fi
fi