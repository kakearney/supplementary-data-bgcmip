# Input arguments:
# -d, --dryrun: Run in dry run mode (set up but don't call ROMS), default false
# --bio [name1 name2 ...]: bio model to run, default ["bestnpz", "cobalt", "banas", "phys"]
# --nbudget: Switch to N-budget archiving options 
# --loop: Repeat run, starting from end of previous (changes initial conditions)
# --bury: For BEST_NPZ or COBALT, switch to using bury-everything variants
# --noinfauna: For BEST_NPZ, turns off infauna grazing

# BGC Hindcasts: testing

import os
import sys
import argparse
from datetime import datetime, timedelta
import netCDF4 as nc
import subprocess

import romscom.romscom as rc
import romscom.rcutils as rcu

# Parse input

ap = argparse.ArgumentParser()
ap.add_argument("--bio", nargs="+", choices=["bestnpz", "cobalt", "banas", "phys"], default=["bestnpz", "cobalt", "banas", "phys"])
ap.add_argument("-d", "--dryrun", action="store_true")
ap.add_argument("--nbudget", action="store_true")
ap.add_argument("--loop", action="store_true")
ap.add_argument("--bury", action="store_true")
ap.add_argument("--noinfauna", action="store_true")
args = ap.parse_args()

# BGC setup

biocpp   = {"bestnpz": "BEST_NPZ", 
            "cobalt":  "BIO_COBALT", 
            "banas":   "BIO_BANAS", 
            "phys":    ""}
          
romsexec = {"bestnpz":     "./romsM_bestnpz_202305311611",
            "cobalt":       "./romsM_cobalt_202306231457_oneill",
            "banas":         "./romsM_banas_202305241654",
            "phys":           "./romsM_phys_202305171514"}
            
if args.bury:
    romsexec["bestnpz"] = "./romsM_bestnpz_202310311628"  # version with fracBurial input parameter 
    romsexec["cobalt"] =   "./romsM_cobalt_202310311639"  # version with NLIM_ONEILL and COBALT_ALLBURIAL on

mpiexec = "mpirun"

# Start and end times

if args.nbudget:
    maxdate = datetime(1994,1,1)
else:
    maxdate = datetime(2020,1,1)
    
dstart = datetime(1990,1,15,0,0)

# Read default parameters

appfol = "../bering-Apps/Apps/Bering_BGC_variants/"

ocean   = rc.readparamfile(os.path.join(appfol, "bering_ocean.yaml"), tconvert=True)
station = rc.readparamfile(os.path.join(appfol, "bering_spos.yaml"))
ice     = rc.readparamfile(os.path.join(appfol, "bering_ipar.yaml"))

bio = {};
for btmp in args.bio:
    if btmp != "phys":
        bio[btmp] = rc.readparamfile(os.path.join(appfol, f"bering_bpar_{btmp}.yaml"))
        
# Bio parameter adjustments

if args.bury and "bestnpz" in args.bio:
    bio["bestnpz"]['fracUnburied'] = 0.0
    if args.noinfauna:
        bio["bestnpz"]["prefD" ] = 0.0 # turn off benthos feeding
        bio["bestnpz"]["prefPL"] = 0.0
        bio["bestnpz"]["prefPS"] = 0.0

# Set output flags

if args.nbudget:
    
    # Turn on all source/sink flux diagnostics
    
    if "banas" in args.bio:
        for k in bio["banas"]['Dout']: 
            bio["banas"]['Dout'][k] = True # BIO_BANAS: turn on all diagnostics
            
    if "bestnpz" in args.bio:
        for k in bio["bestnpz"]['Dout']:
            if k.startswith("iflx"):
                bio["bestnpz"]['Dout'][k] = True
                
    if "cobalt" in args.bio:
        for k in bio["cobalt"]['Dout']:
            if k.startswith("iflx"):
                bio["cobalt"]['Dout'][k] = True
                
    # Also turn on advection/diffusion related diagnostics
    
    for s in args.bio:
        bio[s]['Dout']['iThadv'] = True
        bio[s]['Dout']['iTvadv'] = True
        bio[s]['Dout']['iThdif'] = True
        bio[s]['Dout']['iTvdif'] = True
    
    # Keep only tracer state variables plus zeta, aice, and hice in averages, 
    # and remove diagostics unrelated to bio tracers
    
    for k in ocean['Aout']:
        ocean['Aout'][k] = False
    ocean['Aout']['idFsur'] = True
    ocean['Aout']['idTvar'] = True
    ocean['Aout']['idIceTvar'] = True
    ocean['Aout']['idBeTvar'] = True
    
    for k in ocean['Dout']:
        ocean['Dout'][k] = False
        
    for k in ice['Aout']:
        ice['Aout'][k] = False
    ice['Aout']['idAice'] = True
    ice['Aout']['idHice'] = True
        
    # Not doing anything with history files, so turn it all off
    
    for k in ocean['Hout']:
        ocean['Hout'][k] = False
    for k in ice['Hout']:
        ice['Hout'][k] = False
        
    for s in args.bio:
        for k in bio[s]['Hout']:
            bio[s]['Hout'][k] = False
else:
    # Turn on diagnostics related to primary model-comparison metrics
    
    if "banas" in args.bio:
        for k in bio["banas"]['Dout']: 
            bio["banas"]['Dout'][k] = True # BIO_BANAS: turn on all diagnostics
        
    if "bestnpz" in args.bio:
        bio["bestnpz"]["Dout"]["iprod_PhS"]        = True
        bio["bestnpz"]["Dout"]["iprod_PhL"]        = True
        bio["bestnpz"]["Dout"]["iprod_MZL"]        = True
        bio["bestnpz"]["Dout"]["iprod_Cop"]        = True
        bio["bestnpz"]["Dout"]["iprod_NCaS"]       = True
        bio["bestnpz"]["Dout"]["iprod_EupS"]       = True
        bio["bestnpz"]["Dout"]["iprod_NCaO"]       = True
        bio["bestnpz"]["Dout"]["iprod_EupO"]       = True
        bio["bestnpz"]["Dout"]["iprod_Jel"]        = True
        bio["bestnpz"]["Dout"]["iprod_Ben"]        = True
        bio["bestnpz"]["Dout"]["iprod_IcePhL"]     = True
        bio["bestnpz"]["Dout"]["ipar"]             = True
        bio["bestnpz"]["Dout"]["inolims"]          = True
        bio["bestnpz"]["Dout"]["inoliml"]          = True
        bio["bestnpz"]["Dout"]["inhlims"]          = True
        bio["bestnpz"]["Dout"]["inhliml"]          = True

    if "cobalt" in args.bio:
        bio["cobalt"]["Dout"]["inpp_sm"]   = True
        bio["cobalt"]["Dout"]["inpp_md"]   = True
        bio["cobalt"]["Dout"]["inpp_lg"]   = True
        bio["cobalt"]["Dout"]["inpp_di"]   = True
        bio["cobalt"]["Dout"]["ifratio"]   = True
        bio["cobalt"]["Dout"]["iprod_smz"] = True
        bio["cobalt"]["Dout"]["iprod_mdz"] = True
        bio["cobalt"]["Dout"]["iprod_lgz"] = True
    
# Time variables 

ocean['DSTART'] = dstart
# ocean['NTIMES'] = enddate - ocean['DSTART']

# Set archiving time steps and file size

ocean['NRST'] = timedelta(days=4)
ocean['NSTA'] = timedelta(days=1)
ocean['NHIS'] = timedelta(days=4)
ocean['NDEFHIS'] = timedelta(days=40)
if args.nbudget:
    ocean['NDIA'] = timedelta(days=1)
    ocean['NAVG'] = timedelta(days=1)
    ocean['NDEFAVG'] = timedelta(days=7)
    ocean['NDEFDIA'] = timedelta(days=7)
else:
    ocean['NDIA'] = timedelta(days=4)
    ocean['NAVG'] = timedelta(days=4)
    ocean['NDEFAVG'] = timedelta(days=40)
    ocean['NDEFDIA'] = timedelta(days=40)

# Change a few input files to reflect data folder location

datafol = "../ROMS_Datasets"

ocean["GRDNAME"]  = f"{datafol}/grids/AlaskaGrids_Bering10K.nc"
ocean["TIDENAME"] = f"{datafol}/OTPS/tides_OTPS_Bering10K.nc"
ocean["NUDNAME"]  = f"{datafol}/initial/nudgingcoeff_Bering10K.nc"

# Adjust input based on years to run 
# (template file only holds a few years worth)
# Note: Ideally we'd just run the whole thing at once.  But ROMS does *not* like
# it when you supply more than ~100 forcing multi-files (there's nothing in the
# documentation about this... but somewhere under the hood a variable must be
# allocated to only hold a limited number of characters, and if you surpass that
# things go sideways in a messy-crash-with-very-unhelpful-errors sort of way!
# So, to accomodate, we run in blocks).

blocksz = 10 # years per call, to keep file numbers under the ROMS limit

for bioname in args.bio:
        
    # Start with time-invariant inputs
    
    pflag = bioname == "phys"

    if args.loop:
        if args.nbudget:
            simdir  = f"bgcmip_loop_nbudget_{bioname}"
            simname = f"bgcmip_loop_nbudget_{bioname}"
        else:
            simdir  = f"bgcmip_loop_{bioname}"
            simname = f"bgcmip_loop_{bioname}"
    else:
        if args.nbudget:
            simdir  = f"bgcmip_nbudget_{bioname}"
            simname = f"bgcmip_nbudget_{bioname}"
        else:
            simdir  = f"bgcmip_{bioname}"
            simname = f"bgcmip_{bioname}"
            
    if args.bury:
        simdir  = simdir  + "bury"
        simname = simname + "bury"
        
    if args.noinfauna:
        simdir  = simdir  + "noinfauna"
        simname = simname + "noinfauna"

    # Status message
    
    print("*****************************************************")
    print(f"Preparing simulation: {simname}")
    print("-----------------------------------------------------")

    # Set up sim folders

    fol = rc.simfolders(simdir, create=True)

    # Write accessory files

    bpar = os.path.join(fol['in'], f"{simname}_bpar.in") # bio
    ipar = os.path.join(fol['in'], f"{simname}_ipar.in") # ice
    spos = os.path.join(fol['in'], f"{simname}_spos.in") # stations

    if not pflag:
        rc.dict2standardin(bio[bioname], compress=False, file=bpar)
    rc.dict2standardin(ice, compress=False, file=ipar)
    rc.dict2standardin(station, compress=False, file=spos)

    if not pflag:
        ocean['BPARNAM'] = bpar
    ocean['IPARNAM'] = ipar
    ocean['SPOSNAM'] = spos

    # Point to correct initialization file
    # For loop 1: bio-specific version, same WOA data but units vary between bio models
    # For loop 2: shift last restart slice from to loop 1 to Jan 15, 1990
    
    if pflag:
        ocean['ININAME'] = os.path.join(datafol, 'initial', 'ini_hindcast_unnested_Bering10K_BEST_NPZ.nc')
        ocean["CLMNAME"] = os.path.join(datafol, 'initial', 'ini_hindcast_unnested_Bering10K_BEST_NPZ.nc')
    else:
        if args.loop:
            
            loop1simname = simname.replace("_loop", "")
            loop1simname = loop1simname.replace("_nbudget", "")
            
            ini1 = os.path.join('ininame_second_loop', f"{loop1simname}_lastrst.nc")
            ini2 = os.path.join('ininame_second_loop', f"{loop1simname}_lastrst_shiftedTo19900115.nc")
            
            ocean['ININAME'] = ini2
            
            if not os.path.isfile(ocean['ININAME']):
                print("  Creating ININAME file")
                # Info from previous loop
                loop1fol = rc.simfolders(loop1simname, create=False)
                loop1rst = rcu.parserst(os.path.join(loop1fol['out'], loop1simname))
                # Find last restart file and last index
                f = nc.Dataset(loop1rst['lastfile'], 'r')
                tval = f.variables['ocean_time'][:]
                idx = tval.argmax()
                # Slice last index into new file
                subprocess.run(["ncks", "-d", f"ocean_time,{idx},{idx}", loop1rst['lastfile'], ini1])
                # Shift date
                subprocess.run(["cdo", f"settaxis,1990-01-15,00:00:00", ini1, ini2])
        else:
            ocean['ININAME'] = os.path.join(datafol, 'initial', f"ini_hindcast_unnested_Bering10K_{biocpp[bioname]}.nc")
           
        # Clima file is the same regardless of looping, uses WOA data
        ocean["CLMNAME"] = os.path.join(datafol, 'initial', f"ini_hindcast_unnested_Bering10K_{biocpp[bioname]}.nc")

    # varinfo setup: point to correct bio-specific varinfo (allowing use of identical bry files)

    if pflag:
        ocean['VARNAME'] = os.path.join(appfol, f"varinfo_{bioname}.dat")
    else:
        ocean['VARNAME'] = os.path.join(appfol, f"varinfo_{bioname}_scaledbry.dat")
    
    # Break simulation into blocks to avoid the (undocumented) limit on ROMS 
    # input multi-file number.
    
    for yr1 in range(dstart.year,maxdate.year+1,blocksz):
        
        endblock = min(maxdate, datetime(yr1+blocksz,1,1))
        yrs = range(max(yr1-1,dstart.year), min(endblock.year+2, maxdate.year+1)) # Adds a little overlap for restarting from one block to the next
        
        ocean["FRCNAME"] = [f"{datafol}/BarrowCO2/atmo_co2_barrow_1970_2020.nc", 
                            f"{datafol}/Iron/ESM4_Bering10K_iron_dust_clim.nc",
                            f"{datafol}/salinity/sss.clim.nc",
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-Pair-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-Qair-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-Tair-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-Uwind-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-Vwind-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-rain-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-swrad-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-atmos-northPacific-lwrad-{x}.nc" for x in yrs],
                            [f"{datafol}/GloFAS/GloFAS_runoff_Bering10K_{x}.nc" for x in yrs],
                            [f"{datafol}/GloFAS/GloFAS-based_nutrientflux_Bering10K_{x}.nc" for x in yrs]
                           ]
        ocean["NFFILES"] = len(ocean["FRCNAME"])

        ocean["BRYNAME"] = [f"{datafol}/WOA2018/WOA2018_Bering10K_N30_brybgc.nc",
                            [f"{datafol}/CFS/{x}/CFS-ocean-Bering10K-N30-bryocn-{x}.nc" for x in yrs],
                            [f"{datafol}/CFS/{x}/CFS-ocean-ESPER-Bering10K-N30-brycarbon-{x}.nc"  for x in yrs]
                           ]
        ocean["NBCFILES"] = len(ocean["BRYNAME"])

        # ocean['NTIMES'] = min(maxdate, datetime(yr1+blocksz,1,1)) - ocean['DSTART']

        # Reset restart/initialize marker, in case changed by one of the other bio 
        # sims (runtodate will adjust as necessary)
    
        if args.loop:
            ocean['NRREC'] = -1
        else:
            ocean['NRREC'] = 0    
    
        # Run
    
        romscmd = [mpiexec, romsexec[bioname]]
    
        status = rc.runtodate(ocean, simdir, simname, enddate=endblock, romscmd=romscmd, 
                     dryrunflag=args.dryrun)
                     
        if status != "success":
            break
                 
                 
