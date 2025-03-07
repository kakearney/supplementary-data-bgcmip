/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2017 The ROMS/TOMS Group
**
**   Licensed under a MIT/X style license
**
**   See License_ROMS.txt
**
*******************************************************************************
**
**  Options for BERING simulation
*/

#undef NETCDF4              // use classic netCDF 
#undef PARALLEL_IO          // no parallel input/output
#undef OFFLINE_FLOATS       // could be used for floats in offline version... but not now
#undef USE_CICE             // no CICE ice... should be undefined by default but hitting issues building

/* general */

#define CURVGRID            // use curvilinear coordinates
#define MASKING             // use land/sea masking
#define NONLIN_EOS          // nonlinear equation of state
#define SOLVE3D             // 3D primitive equations
#define SALINITY            // have salinity
#ifdef SOLVE3D               
# undef SPLINES             // turn off option for parabolic splines reconstruction of vertical derivatives
#endif                      
#undef FLOATS               // toggle on/off floats
#define STATIONS            // toggle on/off stations output
#undef WET_DRY              // no wetting/drying of grid cells
                          
#undef T_PASSIVE            // no passive tracers
#ifdef T_PASSIVE            
# define ANA_PASSIVE        // ... but if on, use analytical initial conditions for them
#endif                 
                
/* salinity nudging */ 
                
#define SCORRECTION         // freshwater flux correction

/* ice */

#ifdef SOLVE3D
# define  ICE_MODEL         // Turn on default ice model with...
# ifdef ICE_MODEL           
#  define  ICE_THERMO       // ... ice thermodynamic component
#  define  ICE_MK           // ... Mellor-Kantha thermodynamics
#  undef   ICE_ALB_EC92
#  undef   ICE_SMOOTH
#  define  ICE_MOMENTUM     // ... ice momentum component
#  define  ICE_MOM_BULK     // ... something related to ice-water stress computation
#  define  ICE_EVP          // ... elastic-viscous-plastic rheology
#  define  ICE_ADVECT       // ... advect ice tracers
#  define  ICE_SMOLAR       // ... MPDATA advection scheme
#  define  ICE_UPWIND       // ... upwind advection scheme
#  define  ICE_BULK_FLUXES  // ... ice in bulk flux computation
#  define  ANA_AIOBC        // ... analytical aice boundary conditions (defaults to 0)
#  define  ANA_HIOBC        // ... analytical hice boundary conditions (defaults to 0)
#  define  ANA_HSNOBC       // ... analytical snow boundary conditions (defaults to 0)
# endif
#endif

/* output stuff */
 
#define NO_WRITE_GRID       // Don't write grid arrays
#undef OUT_DOUBLE           // Don't force double precision
#define RST_SINGLE          // Use single precision for restart files
#define AVERAGES            // Write out averages output
#undef AVERAGES2            // No secondary averages output
 
/* advection, dissipation, pressure grad, etc. */
 
#ifdef SOLVE3D
# define DJ_GRADPS          // use splines density Jacobian (Shchepetkin, 2000) in pressure graident term
#endif
 
#define UV_ADV              // turn on advection terms
#define UV_COR              // turn on Coriolis terms
#define UV_SADVECTION       // turn on splines vertical advection
 
#define UV_VIS2             // turn on harmonic horizontal mixing, momentum 
#define UV_SMAGORINSKY      // turn on Smagorinky-like viscosity 
#define VISC_3DCOEF         // turn on time-invarant horizontal viscosity at rho-points
#define MIX_S_UV            // mixing along constant S-surfaces 
#define VISC_GRID           // scale viscosity coefficient by grid size

#ifdef SOLVE3D
# define TS_DIF2            // turn on harmonic horizontal mixing, tracers 
# define MIX_GEO_TS         // mix along geopotential (constant z) surfaces
# define DIFF_GRID          // scales diffusion coefficients by grid size
#endif
 
 
/* vertical mixing */
 
#ifdef SOLVE3D
# define SOLAR_SOURCE       // solar radiation source term
 
# define LMD_MIXING         // Use Large et al. (1994) interior closure with ...
# ifdef LMD_MIXING
#  define LMD_RIMIX         // ... diffusivity due to shear instability
#  define LMD_CONVEC        // ... convective mixing due to shear instability
#  define LMD_SKPP          // ... surface boundary layer KPP mixing
#  undef LMD_BKPP           // ... no bottom boundary KPP mixing
#  define LMD_NONLOCAL      // ... nonlocal transport
#  define LMD_SHAPIRO       // ... Shapiro filtering boundary layer depth
#  undef LMD_DDMIX          // ... no double-diffusive mixing
# endif
 
# undef GLS_MIXING          // Don't use alternative mixing schemes
# undef MY25_MIXING

# if defined GLS_MIXING || defined MY25_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
# endif
#endif
 
/* surface forcing */
 
#ifdef SOLVE3D
# define CORE_FORCING       // input humidity is specific humidity, not relative
# define BULK_FLUXES        // use bulk fluxes computation...
# define CCSM_FLUXES        // ... specifically, the CCSM version of bulk fluxes computation
# if defined BULK_FLUXES || defined CCSM_FLUXES
#  define LONGWAVE_OUT      // compute outgoing longwave radiation (with downward provided as input)
#  define DIURNAL_SRFLUX    // add diurnal cycle to daily-averaged shortwave input
#  define EMINUSP           // compute evap minus precip
#  undef ANA_SRFLUX         // no analytical surface fluxes
#  undef ALBEDO             // Don't calculate shortwave using global albedo equation (use input plus diurnal instead)
#  define ALBEDO_CURVE      // albedo function of lat from Large and Yeager
#  undef LONGWAVE           // Not using net longwave
# endif
#endif
 
/* surface and side corrections */
 
#ifdef SOLVE3D
# undef SRELAXATION         // No salinity relaxation
# undef QCORRECTION         // No net heat flux correction
#endif
 
#ifdef SOLVE3D
# undef TCLIMATOLOGY        // No tracer climatology
# undef TCLM_NUDGING        // No tracer nudging
#endif
 
/* point sources (rivers, line sources) */
/* Using Runoff instead now             */

#ifdef SOLVE3D
# define RUNOFF             // Add runoff as an additional rain field
# define ANA_PSOURCE        // Use analytical point sources
#endif
 
/* tides */
 
#define LTIDES              // Turn on tides (Not a ROMS CPP option, just used here to turn some stuff on/off in bulk)
#ifdef LTIDES
# undef FILTERED            // don't turn on filters... may need on eventually  KAK: what filters, exactly?
# define SSH_TIDES          // impose tidal elevation
# define UV_TIDES           // impose tidal currents
# define ADD_FSOBC          // add tidal elevation to processed OBC data
# define ADD_M2OBC          // add tidal currents  to processed OBC data
# undef RAMP_TIDES          // don't ramp tidal forcing over a day
# define TIDES_ASTRO        // calculate astronomical phase argument
# define POT_TIDES          // impose potential tides
# define UV_LDRAG           // turn on linear bottom friction
# define UV_DRAG_GRID       // use spatially-varying linear coefficient of bottom drag
# define LIMIT_BSTRESS      // limit bottom stress to not change direction of momentum
# undef UV_QDRAG
#else
# define UV_QDRAG           // quadratic bottom stress
#endif
 
/* Boundary conditions...careful with grid orientation */
/* BERING_10K: north = Russia (northwest-ish), closed
               south = Gulf of Alaska (southeast-ish), open
               east = Bering Strait, closed but with momentum point source (see ana_psource)
               west = North Pacific (south of Aleutians), open */
 
#define RADIATION_2D
 
/* roms quirks */
 
#ifdef SOLVE3D
# define ANA_BSFLUX
# define ANA_BTFLUX
#else
# define ANA_SMFLUX
#endif

/* MPI stuff (see https://www.myroms.org/projects/src/ticket/747)      */
/* not strictly necessary, but enforces same behavior as older version */

# define BOUNDARY_ALLREDUCE /* use mpi_allreduce in mp_boundary */
# undef  COLLECT_ALLGATHER  /* use mpi_allgather in mp_collect  */
# define COLLECT_ALLREDUCE  /* use mpi_allreduce in mp_collect  */
# define REDUCE_ALLGATHER   /* use mpi_allgather in mp_reduce   */
# undef  REDUCE_ALLREDUCE   /* use mpi_allreduce in mp_reduce   */

/* Diagnostics */

# define DIAGNOSTICS_TS



/* COBALT module options */

#ifdef BIO_COBALT
# undef TS_MPDATA
# define TS_HSIMT
#elif defined SOLVE3D
# define TS_U3HADVECTION
# define TS_C4VADVECTION
#endif

/*#define DEBUG_COBALT */
/*#define COBALT_CONSERVATION_TEST */
/*#define COBALT_NOSOURCE */
/*#define COBALT_DO_NOTHING  */

#if defined BIO_COBALT
# undef FILTERED
# undef AVERAGES2
# define OPTIC_MANIZZA     /* Manizza light attenuation... */
# define COASTAL_ATTEN     /* ... with additional coastal attenuation */
# define COBALT_MINERALS
# undef COBALT_PHOSPHORUS
# define COBALT_IRON
# define NO_IRON_COAST
# define COBALT_CARBON
# define DIAGNOSTICS_BIO
/*# define BENTHIC  */
/*# define TIMESERIES */
# undef ANA_ALK_DIC
# undef ANA_BIOLOGY        /* analytical biology initial conditions */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# undef COASTDIAT
#endif

#ifdef BEST_NPZ
# define JELLY             /* Add jellyfish */
# define IRON_LIMIT        /* Add iron  */
# define BENTHIC           /* Add benthos (infauna and detritus) */
# define ICE_BIO           /* Add ice bio (PhL, NO3, NH4) */
# undef CLIM_ICE_1D
# define DIAPAUSE          /* turn on seasonal vertical migration for large copepods */
# define OPTIC_MANIZZA     /* Manizza light attenuation... */
# define COASTAL_ATTEN     /* ... with additional coastal attenuation */
# if defined CARBON
#  define CARBON_FLUX      /* For river fluxes of DIC,TA */
#  define OXYGEN           /* For oxygen cycling */
# endif
# define DIAGNOSTICS_BIO   /* diagnostics on */
# define DIAGBIOAVG        /* averages-like diagnostics for select diags */
# define ANA_ICEBIOBC      /* Analytical ice bio boundary conditions */
# undef ANA_BIOLOGY
# define TCLM_NUDGING      /* Nudging of tracer climatology for iron */
# undef ANA_TCLIMA         /* analytical tracers climatology for iron */
# undef TCLIMATOLOGY       /* Processing of tracer climatology for iron (now deprecated) */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
#endif

#ifdef BIO_BANAS
# define OPTIC_MANIZZA     /* Manizza light attenuation... */
# define COASTAL_ATTEN     /* ... with additional coastal attenuation */
# define ANA_BPFLUX        /* analytical bottom passive tracers fluxes */
# define ANA_SPFLUX        /* analytical surface passive tracers fluxes */
# define DIAGNOSTICS_BIO   /* turn on diagnostics */
# define DIAGBIOAVG        /* averages-like diagnostics */
#endif