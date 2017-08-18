      subroutine SF_STOM  (LAND_PTS, LAND_INDEX, FT, CO2, FSMC, HT, IPAR
     &,                    LAI, PSTAR, Q1, RA, TSTAR, ZERODEGC, EPCO2
     &,                    EPSILON, GPP, NPP, RESP_P, RESP_W, GC)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Routine to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes

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

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of vegetated points.
! LAND_INDEX = IN Index of vegetated points.
! FT        = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), FT, J, K, L

! CO2      = IN Atmospheric CO2 concentration (kg CO2/kg air).
! FSMC     = IN Soil water factor.
! HT       = IN Canopy height (m).
! IPAR     = IN Incident PAR (W/m2).
! LAI      = IN Leaf area index.
! PSTAR    = IN Surface pressure (Pa).
! Q1       = IN Specific humidity of level 1 (kg H2O/kg air).
! RA       = IN Aerodynamic resistance (s/m).
! TSTAR    = IN Surface temperature (K).
! ZERODEGC = IN Zero Celsius (K).
! EPCO2    = IN Ratio of molecular weights of CO2 and dry air.
! EPSILON  = Ratio of molecular weights of water and dry air.
! GPP      = OUT Gross Primary Productivity (kg C/m2/s).
! NPP      = OUT Net Primary Productivity (kg C/m2/s).
! RESP_P   = OUT Plant respiration rate (kg C/m2/sec).
! RESP_W   = OUT Wood respiration rate (kg C/m2/sec).
! GC       = INOUT Canopy resistance to H2O (m/s).
! ANETC    = WORK Net canopy photosynthesis (mol CO2/m2/s).
! CI       = WORK Internal CO2 pressure (Pa).
! DQ       = WORK Specific humidity deficit (kg H2O/kg air).
! DQC      = WORK Canopy level specific humidity deficit
!            (kg H2O/kg air).
! FPAR     = WORK PAR absorption factor.
! LAI_BAL  = WORK Leaf area index in balanced growth state.
! NL       = WORK Mean leaf nitrogen concentration (kg N/kg C).
! NL_BAL   = WORK Mean leaf nitrogen concentration in balanced growth
!            state (kg N/kg C).
! N_LEAF   = WORK Nitrogen contents of the leaf and root (kg N/m2).
! N_STEM   = WORK Nitrogen contents of the stem (kg N/m2).
! QS       = WORK Saturated specific humidity(kg H2O/kg air).
! RA_RC    = WORK Ratio of aerodynamic resistance to canopy resistance.
! RDC      = WORK Canopy dark respiration, without soil water dependence
!            (mol CO2/m2/s).
! RESP_P_G = WORK Plant growth respiration rate (kg C/m2/sec).
! RESP_P_M = WORK Plant maintenance respiration rate (kg C/m2/sec).
! ROOT     = WORK Root carbon (kg C/m2).

      real CO2(POINTS), FSMC(POINTS), HT(POINTS), IPAR(POINTS)
      real LAI(POINTS), PSTAR(POINTS), Q1(POINTS), RA(POINTS)
      real TSTAR(POINTS), ZERODEGC, EPCO2, EPSILON, GPP(POINTS)
      real RESP_P(POINTS), RESP_W(POINTS), GC(POINTS), ANETC(POINTS)
      real NPP(POINTS), CI(POINTS), DQ(POINTS), DQC(POINTS)
      real FPAR(POINTS), LAI_BAL(POINTS), NL(POINTS), NL_BAL(POINTS)
      real N_LEAF(POINTS), N_ROOT(POINTS), N_STEM(POINTS), QS(POINTS)
      real RA_RC(POINTS), RDC(POINTS), RESP_P_G(POINTS)
      real RESP_P_M(POINTS), ROOT(POINTS)

! Local parameters
! ITER = Number of iterations to determine the canopy climate.

      integer ITER
      parameter (ITER=1)

!-----------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit and the surface
! density of the air
!-----------------------------------------------------------------------
      call QSAT (POINTS, LAND_PTS, LAND_INDEX, EPSILON, ZERODEGC
     &,          QS, TSTAR, PSTAR)
      do J=1,LAND_PTS
        L = LAND_INDEX(J)
        DQ(L) = MAX(0.0,(QS(L) - Q1(L)))
      enddo

!-----------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------
      do J=1,LAND_PTS
        L = LAND_INDEX(J)
        FPAR(L) = (1 - EXP(-KPAR(FT)*LAI(L))) / KPAR(FT)
      enddo

!-----------------------------------------------------------------------
! Iterate to ensure that the canopy humidity deficit is consistent with
! the H2O flux. Ignore the (small) difference between the canopy and
! reference level CO2 concentration. Initially set the canopy humidity
! deficit using the previous value of GC.
!-----------------------------------------------------------------------
      do K=1,ITER

!-----------------------------------------------------------------------
! Diagnose the canopy level humidity deficit
!-----------------------------------------------------------------------
        do J=1,LAND_PTS
          L = LAND_INDEX(J)
          RA_RC(L) = RA(L) * GC(L)
          DQC(L) = DQ(L) / (1 + RA_RC(L))
        enddo

!-----------------------------------------------------------------------
! Call CANOPY to calculate the canopy resistance and photosynthesis
!-----------------------------------------------------------------------
        call CANOPY (LAND_PTS, LAND_INDEX, FT, DQC, IPAR, TSTAR, CO2
     &,              PSTAR, FPAR, FSMC, LAI, ZERODEGC, EPCO2, GC
     &,              ANETC, CI, RDC)

      enddo

      do J=1,LAND_PTS
        L = LAND_INDEX(J)

!-----------------------------------------------------------------------
! Assume that root biomass is equal to balanced growth leaf biomass
!-----------------------------------------------------------------------
        LAI_BAL(L) = (A_WS(FT)*ETA_SL(FT)*HT(L)/A_WL(FT))
     &             **(1.0/(B_WL(FT)-1))
        ROOT(L) = SIGL(FT) * LAI_BAL(L)

!-----------------------------------------------------------------------
! Calculate the actual and balanced mean leaf nitrogen concentration
! assuming perfect light acclimation
!-----------------------------------------------------------------------
        NL(L) = (FPAR(L) / LAI(L)) * NL0(FT)
        NL_BAL(L) = (1 - EXP(-KPAR(FT)*LAI_BAL(L)))
     &            / (KPAR(FT)*LAI_BAL(L)) * NL0(FT)

!-----------------------------------------------------------------------
! Calculate the total nitrogen content of the leaf, root and stem
!-----------------------------------------------------------------------
        N_LEAF(L) = NL(L) * SIGL(FT) * LAI(L)
        N_ROOT(L) = NR_NL(FT) * NL_BAL(L) * ROOT(L)
        N_STEM(L) = NS_NL(FT) * NL_BAL(L) * ETA_SL(FT) * HT(L) * LAI(L)

!-----------------------------------------------------------------------
! Calculate the Gross Primary Productivity, the plant maintenance
! respiration rate, and the wood maintenance respiration rate
! in kg C/m2/sec
!-----------------------------------------------------------------------
        GPP(L) = 12.0E-3 * (ANETC(L) + RDC(L)*FSMC(L))
        RESP_P_M(L) = 12.0E-3 * RDC(L)
     &     * (N_LEAF(L)*FSMC(L) + N_STEM(L) + N_ROOT(L)) / N_LEAF(L)
        RESP_W(L) = 12.0E-3 * RDC(L) * N_STEM(L) / N_LEAF(L)

!-----------------------------------------------------------------------
! Calculate the total plant respiration and the Net Primary Productivity
!-----------------------------------------------------------------------
        RESP_P_G(L) = R_GROW(FT) * (GPP(L) - RESP_P_M(L))
        RESP_P(L) = RESP_P_M(L) + RESP_P_G(L)
        NPP(L) = GPP(L) - RESP_P(L)

      enddo
#endif

      return
      end
