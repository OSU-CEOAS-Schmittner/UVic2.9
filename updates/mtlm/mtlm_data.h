!===================== include file "mtlm_data.h" ======================

!     data statements for mtlm

!-----------------------------------------------------------------------
! Functional Type dependent parameters

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
! C3       = 1 for C3 Plants, 0 for C4 Plants.
! AGR      = 1 for agricultural type, 0 for non-agricultural.

      integer C3(NPFT), AGR(NPFT)

! ALPHA    = Quantum efficiency (mol CO2/mol PAR photons).
! A_WL     = Allometric coefficient relating the target woody biomass to
!            the leaf area index (kg C/m2).
! A_WS     = Woody biomass as a multiple of live stem biomass.
! B_WL     = Allometric exponent relating the target woody biomass to
!            the leaf area index.
! DGL_DM   = Rate of change of leaf turnover rate with moisture
!            availability.
! DGL_DT   = Rate of change of leaf turnover rate with temperature (/K)
! DQCRIT   = Critical humidity deficit (kg H2O/kg air).
! ETA_SL   = Live stemwood coefficient (kg C/m/LAI).
! FSMC_OF  = Moisture availability below which leaves are dropped.
! F0       = CI/CA for DQ = 0.
! GLMIN    = Minimum leaf conductance for H2O
! G_AREA   = Disturbance rate (/360days).
! G_GROW   = Rate of leaf growth (/360days).
! G_LEAF_0 = Minimum turnover rate for leaves (/360days).
! G_ROOT   = Turnover rate for root biomass (/360days).
! G_WOOD   = Turnover rate for woody biomass (/360days).
! KPAR     = PAR Extinction coefficient (m2 leaf/m2 ground).
! LAI_MAX  = Maximum projected LAI.
! LAI_MIN  = Minimum projected LAI.
! NL0      = Top leaf nitrogen concentration (kg N/kg C).
! NR_NL    = Ratio of root nitrogen concentration to leaf nitrogen
!            concentration.
! NS_NL    = Ratio of stem nitrogen concentration to leaf nitrogen
!            concentration.
! OMEGA    = Leaf scattering coefficient for PAR.
! R_GROW   = Growth respiration fraction.
! SIGL     = Specific density of leaf carbon (kg C/m2 leaf).
! TLEAF_OF = Temperature below which leaves are dropped.
! TLOW     = Lower temperature for photosynthesis (deg C)
! TUPP     = Upper temperature for photosynthesis (deg C)
! Mar 16, 2016 Andreas next line moved here from mtlm.h
! MAF      = Moisture availability factor PFT/soil dependence

      real ALPHA(NPFT), A_WL(NPFT), A_WS(NPFT), B_WL(NPFT)
      real DGL_DM(NPFT), DGL_DT(NPFT), DQCRIT(NPFT)
      real ETA_SL(NPFT), FSMC_OF(NPFT), F0(NPFT)
      real GLMIN(NPFT), G_AREA(NPFT), G_GROW(NPFT)
      real G_LEAF_0(NPFT), G_ROOT(NPFT),G_WOOD(NPFT)
      real KPAR(NPFT), LAI_MAX(NPFT), LAI_MIN(NPFT)
      real NL0(NPFT), NR_NL(NPFT), NS_NL(NPFT), OMEGA(NPFT)
      real R_GROW(NPFT), SIGL(NPFT), TLEAF_OF(NPFT)
      real TLOW(NPFT), TUPP(NPFT)
      real MAF(NTYPE)

!-----------------------------------------------------------------------
!                        BT     NT    C3G    C4G     S
!-----------------------------------------------------------------------
      data C3      /      1,     1,     1,     0,     1 /
      data AGR     /      0,     0,     1,     1,     0 /
      data ALPHA   /   0.06,  0.06,  0.06, 0.040,  0.06 /
      data A_WL    /   0.65,  0.65, 0.005, 0.005,  0.10 /
      data A_WS    /  10.00, 10.00,  1.00,  1.00, 10.00 /
      data B_WL    /  1.667, 1.667, 1.667, 1.667, 1.667 /
      data DGL_DM  /  100.0, 100.0, 100.0, 100.0, 100.0 /
      data DGL_DT  /    9.0,   9.0,   0.0,   0.0,   9.0 /
      data DQCRIT  /  0.090, 0.060, 0.100, 0.075, 0.100 /
      data ETA_SL  /   0.01,  0.01,  0.01,  0.01,  0.01 /
      data F0      /  0.875, 0.875, 0.900, 0.800, 0.900 /
      data FSMC_OF /   0.85,  0.60,  0.05,  0.00,  0.50 /
      data GLMIN   / 1.0E-6,1.0E-6,1.0E-6,1.0E-6,1.0E-6 /
      data G_AREA  /  0.004, 0.004,  0.10,  0.10,  0.05 /
      data G_GROW  /  20.00, 20.00, 20.00, 20.00, 20.00 /
      data G_LEAF_0/   0.25,  0.25,  0.25,  0.25,  0.25 /
      data G_ROOT  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      data G_WOOD  /   0.01,  0.01,  0.20,  0.20,  0.05 /
      data KPAR    /   0.50,  0.50,  0.50,  0.50,  0.50 /
      data LAI_MAX /   8.00,  8.00,  3.50,  3.50,  3.50 /
      data LAI_MIN /   3.00,  3.00,  1.00,  1.00,  1.00 /
      data NL0     /  0.036, 0.030, 0.054, 0.027, 0.027 /
      data NR_NL   /   2.00,  2.00,  2.00,  2.00,  2.00 /
      data NS_NL   /   0.10,  0.10,  1.00,  1.00,  0.10 /
      data OMEGA   /   0.15,  0.15,  0.15,  0.17,  0.15 /
      data R_GROW  /   0.25,  0.25,  0.25,  0.25,  0.25 /
      data SIGL    / 0.0375,0.1000,0.0250,0.0500,0.0500 /
      data TLEAF_OF/ 273.15,243.15,258.15,258.15,243.15 /
      data TLOW    /  -10.0, -15.0,  -5.0,   8.0, -10.0 /
      data TUPP    /   33.0,  25.0,  33.0,  42.0,  33.0 /

      data MAF / 1.0, 1.0, 0.95, 0.95, 0.97, 0.95 /

