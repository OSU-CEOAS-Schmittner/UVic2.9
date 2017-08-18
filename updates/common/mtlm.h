!======================== include file "mtlm.h" ========================

!**********************************************************************
! this file is based on code that may have had the following copyright:
! (c) CROWN COPYRIGHT 1997, U.K. METEOROLOGICAL OFFICE.

! Permission has been granted by the authors to the public to copy
! and use this software without charge, provided that this Notice and
! any statement of authorship are reproduced on all copies. Neither the
! Crown nor the U.K. Meteorological Office makes any warranty, express
! or implied, or assumes any liability or responsibility for the use of
! this software.
!**********************************************************************

!-----------------------------------------------------------------------
! Atmospheric CO2 variables
!-----------------------------------------------------------------------
! CO2   = CO2 mass mixing ratio (kg/kg).
! EPCO2 = Ratio of molecular weights of CO2

      real EPCO2, CO2
      parameter (EPCO2 = 1.5194)

      common /land_r/ CO2(POINTS)

!-----------------------------------------------------------------------
! Driving variables
!-----------------------------------------------------------------------
! LYING_SNOW = Snow mass (kg/m2).
! LW         = Surface longwave radiation (W/m**2).
! SW         = Surface shortwave radiation (W/m**2).
! SWN        = Net shortwave radiation (W/m**2).
! PSTAR      = Surface pressure (Pa).
! Q          = Specific humidity (kg/kg).
! SAT_D      = Diurnal surface air temperature (K).
! TS1        = Sub-surface temperature (K).
! TIMEDAY    = Time of day (s).
! W          = Water incident at the soil surface (mm/day).
! WIND       = Wind speed (m/s).
! LW_OUT     = total net longwave radiation (W/m2)
! LHC        = Latent heat of condensation (J/kg).
! LHF        = Latent heat of fusion (J/kg).
! SIGMA      = Stefan-Boltzman constant (W/m2/K4).

      real LYING_SNOW, LW, SW, SWN, PSTAR, Q, SAT_D, TS1, TIMEDAY
      real W, WIND, LW_OUT, LHC, LHF, SIGMA

      common /land_r/ LYING_SNOW(POINTS), LW(POINTS), SW(POINTS)
      common /land_r/ SWN(POINTS), PSTAR(POINTS), Q(POINTS)
      common /land_r/ SAT_D(POINTS), TS1(POINTS), TIMEDAY
      common /land_r/ WIND(POINTS), LW_OUT(POINTS), LHC, LHF, SIGMA

!-----------------------------------------------------------------------
! Variables for forcing with climatological means
!-----------------------------------------------------------------------
! DTEMP_DAY    = Diurnal temperature range (K).
! LYING_SNOW_C = Snow mass (kg/m2).
! QS           = Saturated specific humidity (kg/kg).
! RAIN         = Rainfall rate (kg/m2/s).
! SNOW         = Snowfall rate (kg/m2/s).
! SNOWMELT     = Snow melt (mm/day).
! SW_C         = Surface shortwave radiation (W/m**2).
! SUN          = Normalized solar radiation.
! SURF_ROFF    = Surface runoff (mm/day).
! TIME_MAX     = GMT at which maximum temperature occurs (s).
! T_C          = Air temperature (K).
! W_C          = Water incident at the soil surface (mm/day).
! RH_C         = relative humidity from atmospheric model

      real DTEMP_DAY, QS, RAIN, SNOW, SNOWMELT, SW_C, SUN, SURF_ROFF
      real TIME_MAX, T_C, W_C, RH_C

      common /land_r/ DTEMP_DAY(POINTS), QS(POINTS), RAIN(POINTS)
      common /land_r/ SNOW(POINTS), SNOWMELT(POINTS), SW_C(POINTS)
      common /land_r/ SUN(POINTS,STEPSM), SURF_ROFF(POINTS)
      common /land_r/ TIME_MAX(POINTS), T_C(POINTS), W_C(POINTS)
      common /land_r/ RH_C(POINTS)

!-----------------------------------------------------------------------
! Model parameters
!-----------------------------------------------------------------------
! TIMESTEP    = Timestep for daily calculations (s).

      real TIMESTEP
      common /land_r/ TIMESTEP

!-----------------------------------------------------------------------
! Soil parameters
!-----------------------------------------------------------------------
! ALBSOIL   = Soil albedo.
! ALBSNOW   = Snow albedo.
! ALBLAND   = Surface albedo.
! Z0S       = Roughness length for bare soil (m).
! ROOTDEP   = Rootdepth (m).
! HCAP_SOIL = Soil heat capacity (W/m3/K).
! HCON_SOIL = Soil heat conductivity (W/m/K).

      real ALBSOIL, ALBSNOW, ALBLAND, Z0S, ROOTDEP, HCAP_SOIL, HCON_SOIL

!      parameter (ROOTDEP=1.0, HCAP_SOIL=3.3E5, HCON_SOIL=0.23)
      parameter (ROOTDEP=1.0, HCAP_SOIL=3.3E5, HCON_SOIL=0.75)
      common /land_r/ ALBSOIL(POINTS), ALBSNOW(POINTS)
      common /land_r/ ALBLAND(POINTS), Z0S(POINTS)

!-----------------------------------------------------------------------
! Vegetation parameters
!-----------------------------------------------------------------------
! ALBSNC    = Cold deep snow albedo.
! ALBSNF    = Snow free albedo.
! CATCH     = Canopy capacity (kg/m2).
! Z0        = Vegetative roughness length (m).
! VEG_FRAC  = Vegetated fraction.
! FRAC_VS   = Total fraction of gridbox covered by vegetation and soil.
! FRAC_MIN  = Minimum areal fraction for PFTs.
! FRAC_SEED = "Seed" fraction for PFTs.

      real ALBSNC, ALBSNF, CATCH, Z0, VEG_FRAC, FRAC_VS, FRAC_MIN
      real FRAC_SEED
      parameter (FRAC_MIN=1.0E-6, FRAC_SEED=0.01)

      common /land_r/ ALBSNC(POINTS,NPFT), ALBSNF(POINTS,NPFT)
      common /land_r/ CATCH(POINTS,NPFT)
      common /land_r/ Z0(POINTS,NPFT), VEG_FRAC(POINTS)
      common /land_r/ FRAC_VS(POINTS)

!-----------------------------------------------------------------------
! IVM Prognostics
!-----------------------------------------------------------------------
! CS   = Soil carbon (kg C/m2).
! FRAC = Areal coverage.
! LAI  = Leaf area index.

      real CS, FRAC, LAI

      common /land_r/ CS(POINTS), FRAC(POINTS,NTYPE)
      common /land_r/ LAI(POINTS,NPFT)

!-----------------------------------------------------------------------
! IVM Fluxes and diagnostics
!-----------------------------------------------------------------------
! C_VEG       = Vegetation carbon (kg C/m2).
! CV          = Gridbox mean vegetation carbon (kg C/m2).
! FTIME       = Weighting factor for accumulations.
! HT          = Canopy height (m).
! GPP         = Gross Primary Productivity (kg C/m2/s).
! G_LEAF      = Leaf turnover rate (/yr).
! G_LEAF_DAY  = Daily mean leaf turnover rate (/yr).
! G_LEAF_PHEN = Daily leaf turnover rate including phenology (/yr).
! NPP         = Net Primary Productivity (kg C/m2/s).
! RESP_S      = Soil respiration rate (kg C/m2/s).
! LIT_C_T     = Gridbox mean carbon litter (kg C/m2/yr).
! BF          = Burn fraction

      real C_VEG, CV, FTIME, HT, GPP, G_LEAF, G_LEAF_DAY
      real G_LEAF_PHEN, NPP, RESP_S, LIT_C_T, BF

      common /land_r/ C_VEG(POINTS,NPFT), CV(POINTS), FTIME
      common /land_r/ HT(POINTS,NPFT), GPP(POINTS,NPFT)
      common /land_r/ G_LEAF(POINTS,NPFT), G_LEAF_DAY(POINTS,NPFT)
      common /land_r/ G_LEAF_PHEN(POINTS,NPFT), NPP(POINTS,NPFT)
      common /land_r/ RESP_S(POINTS), LIT_C_T(POINTS), BF

!-----------------------------------------------------------------------
! Driving variables for TRIFFID
!-----------------------------------------------------------------------
! G_LEAF_DR  = Accumulated leaf turnover rate (/yr).
! NPP_DR     = Accumulated Net Primary Productivity (kg C/m2/yr).
! RESP_S_DR  = Accumulated soil respiration rate (kg C/m2/yr).
! RESP_W_DR  = Accumulated wood respiration rate (kg C/m2/yr).

      real G_LEAF_DR, NPP_DR, RESP_S_DR, RESP_W_DR

      common /land_r/ G_LEAF_DR(POINTS,NPFT), NPP_DR(POINTS,NPFT)
      common /land_r/ RESP_S_DR(POINTS), RESP_W_DR(POINTS,NPFT)

!-----------------------------------------------------------------------
! Hydrology variables
!-----------------------------------------------------------------------
! ESUB       = Sublimation (kg/m2/s).
! FSMC       = Moisture availability factor.
! M          = Total soil moisture (kg/m2).
! MNEG       = Negative soil moisture (kg/m2).

      real ESUB, FSMC, M, MNEG

      common /land_r/ ESUB(POINTS), FSMC(POINTS), M(POINTS)
      common /land_r/ MNEG(POINTS)

!-----------------------------------------------------------------------
! Temperatures
!-----------------------------------------------------------------------
! TSOIL = Temperature of bare soil (K).
! TSTAR = Surface temperature (K).

      real TSOIL, TSTAR

      common /land_r/ TSOIL(POINTS), TSTAR(POINTS,NPFT)

!-----------------------------------------------------------------------
! Inputs defining locations
!-----------------------------------------------------------------------

! LAND_PTS   = Number of land points.
! LAND_INDEX = Indices of land points.
! VEG_PTS    = Number of land points which include the nth vegetation
!              type.
! VEG_INDEX  = Indices of land points which include the nth vegetation
!              type.

      integer LAND_PTS, LAND_INDEX, VEG_PTS, VEG_INDEX

      common /land_i/ LAND_PTS, LAND_INDEX(POINTS), VEG_PTS(NPFT)
      common /land_i/ VEG_INDEX(POINTS,NPFT)

! LAT             = Latitude (degrees)
! LATMIN,LATMAX   = Latitudinal limits of the area (degrees).
! GAREA           = grid area (m2)
! LONG            = Longitude (degrees)
! LONGMIN,LONGMAX = Longitudinal limits of the area (degrees).

      real LAT, LATMIN, LATMAX, LONG, LONGMIN, LONGMAX, GAREA

      common /land_r/ LAT(POINTS),LATMIN, LATMAX, GAREA(POINTS)
      common /land_r/ LONG(POINTS), LONGMIN, LONGMAX

!-----------------------------------------------------------------------
! Time parameters
!-----------------------------------------------------------------------
! ISTEP        = step counter.
! STEP_DAY     = Number of steps in a day.
! LAND_COUNTER = number time steps for the model.

      integer ISTEP, STEP_DAY, LAND_COUNTER

      common /land_i/ ISTEP, STEP_DAY, LAND_COUNTER

! SEC_DAY  = Number of seconds in a day (s).
! DAY_YEAR = Number of days in a year (days).
! SEC_YEAR = Number of seconds in a year (s).

      real DAY_YEAR, SEC_DAY, SEC_YEAR
      common /land_r/ DAY_YEAR, SEC_DAY, SEC_YEAR

! INT_VEG   = .T. for interactive vegetation
! VEG_EQUIL = .T. if the vegetation equilibrium

      logical INT_VEG, VEG_EQUIL

      common /land_l/ INT_VEG, VEG_EQUIL

!-----------------------------------------------------------------------
! Variables defining anthropogenic disturbance
!-----------------------------------------------------------------------
! FRACA     = Areal fraction of agriculture.
! G_ANTH    = Anthropogenic disturbance rate (/yr).

      real FRACA, G_ANTH

      common /land_r/ FRACA(POINTS), G_ANTH(POINTS)

!-----------------------------------------------------------------------
! New variables required for LAND module
!-----------------------------------------------------------------------
! ET       = Evapotranspiration (kg/m2/s).
! LE       = Latent heat flux (W/m2).
! TSTAR_GB = Grid box average surface temperature (K).
! SH       = Sensible heat flux (W/m2).

      real ET, LE, TSTAR_GB, SH

      common /land_r/ ET(POINTS), LE(POINTS), TSTAR_GB(POINTS)
      common /land_r/ SH(POINTS)

! DAY_PHEN = IN Number of days between phenology.
! DAY_TRIF = IN Number of days between TRIFFID.

      integer  DAY_PHEN, DAY_TRIF

      common /land_i/ DAY_PHEN, DAY_TRIF

! L_PHEN = IN .T. if phenology to be updated.
! L_TRIF = IN .T. if vegetation to be updated.

      logical L_PHEN, L_TRIF

      common /land_l/ L_PHEN, L_TRIF

# if defined O_time_averages
!-----------------------------------------------------------------------
! time averaged arrays
!-----------------------------------------------------------------------
! ntatsl        = time step counter for time averaging
! ta_TS1        = time averaged TS1
! ta_TSTAR      = time averaged TSTAR
! ta_ALBLAND    = time averaged ALBLAND
! ta_ET         = time averaged ET
! ta_M          = time averaged M
! ta_CS         = time averaged CS
! ta_RESP_S     = time averaged RESP_S
! ta_LIT_C_T    = time averaged LIT_C_T
! ta_BURN       = time averaged BURN
! ta_GPP        = time averaged GPP
! ta_NPP        = time averaged NPP
! ta_HT         = time averaged HT
! ta_LAI        = time averaged LAI
! ta_C_VEG      = time averaged C_VEG
! ta_LYING_SNOW = time averaged LYING_SNOW
! ta_SURF_ROFF  = time averaged SURF_ROFF
! ta_FRAC       = time averaged FRAC

      integer ntatsl

      common /land_i/ ntatsl

      real ta_TS1, ta_TSTAR_GB, ta_ALBLAND, ta_ET, ta_M, ta_CS
      real ta_RESP_S, ta_LIT_C_T, ta_BURN, ta_GPP, ta_NPP, ta_HT
      real ta_LAI, ta_C_VEG, ta_LYING_SNOW, ta_SURF_ROFF, ta_FRAC

      common /land_r/ ta_TS1(POINTS), ta_TSTAR_GB(POINTS)
      common /land_r/ ta_ALBLAND(POINTS), ta_ET(POINTS), ta_M(POINTS)
      common /land_r/ ta_CS(POINTS), ta_RESP_S(POINTS)
      common /land_r/ ta_LIT_C_T(POINTS), ta_BURN(POINTS)
      common /land_r/ ta_GPP(POINTS,NPFT), ta_NPP(POINTS,NPFT)
      common /land_r/ ta_HT(POINTS,NPFT), ta_LAI(POINTS,NPFT)
      common /land_r/ ta_C_VEG(POINTS,NPFT), ta_LYING_SNOW(POINTS)
      common /land_r/ ta_SURF_ROFF(POINTS), ta_FRAC(POINTS,NTYPE)

# endif
# if defined O_time_step_monitor
!-----------------------------------------------------------------------
! time average integrals
!-----------------------------------------------------------------------
! ntatil              = time step counter for time average integrals
! tai_CS              = time averaged integral CS
! tai_RESP_S          = time averaged integral RESP_S
! tai_LIT_C_T         = time averaged integral LIT_C_T
! tai_BURN            = time averaged integral BURN
! tai_CV              = time averaged integral CV
! tai_NPP             = time averaged integral NPP
! tai_GPP             = time averaged integral GPP
! tai_HT              = time averaged integral HT
! tai_LAI             = time averaged integral LAI
! tai_LYING_SNOW      = time averaged integral LYING_SNOW
! tai_TSOIL           = time averaged integral TSOIL
! tai_TSTAR           = time averaged integral TSTAR
! tai_M               = time averaged integral M
! tai_ET              = time averaged integral ET
! tai_clnd            = average total carbon in land
! tai_cfa2l           = average total flux atmosphere to land

      integer ntatil

      common /land_i/ ntatil

      real tai_CS, tai_RESP_S, tai_LIT_C_T, tai_BURN, tai_CV, tai_NPP
      real tai_GPP, tai_HT, tai_LAI, tai_LYING_SNOW, tai_TSOIL
      real tai_TSTAR, tai_M, tai_ET, tai_clnd, tai_cfa2l

      common /land_r/ tai_CS, tai_RESP_S, tai_LIT_C_T, tai_BURN, tai_CV
      common /land_r/ tai_NPP, tai_GPP, tai_HT, tai_LAI, tai_LYING_SNOW
      common /land_r/ tai_TSOIL, tai_TSTAR, tai_M, tai_ET, tai_clnd
      common /land_r/ tai_cfa2l
# endif
