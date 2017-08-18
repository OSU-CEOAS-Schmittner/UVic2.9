      subroutine CANOPY (LAND_PTS, LAND_INDEX, FT, DQC, IPAR, TSTAR
     &,                  CO2C, PSTAR, FPAR, FSMC, LAI, ZERODEGC
     &,                  EPCO2, GC, ANETC, CI, RDC)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Calculates the canopy resistance, net photosynthesis and transpiration
! by scaling-up the leaf level response using the "Big-Leaf" approach
! of Sellers et al. (1994)

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
! FT         = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), FT, I, L

! DQC      = IN Canopy level specific humidity deficit (kg H2O/kg air).
! PSTAR    = IN Surface pressure (Pa).
! CO2C     = IN Canopy level CO2 concentration (kg CO2/kg air).
! IPAR     = IN Incident PAR (W/m2).
! TSTAR    = IN Surface temperature (K).
! FPAR     = IN PAR absorption factor.
! FSMC     = IN Soil water factor.
! LAI      = IN Leaf area index (m2 leaf/m2 ground).
! ZERODEGC = IN Zero Celsius (K).
! EPCO2    = IN Ratio of molecular weights of CO2 and dry air.
! ANETC    = OUT Net canopy photosynthesis (mol CO2/m2/s).
! CI       = OUT Internal CO2 concentration (mol CO2/m3).
! GC       = OUT Canopy conductance for H2O (m/s).
! RDC      = OUT Canopy dark respiration (mol CO2/m2/s).
! ANETL    = WORK Net leaf photosynthesis (mol CO2/m2/s/LAI).
! APAR     = WORK PAR absorbed by the top leaf (W/m2).
! CA       = WORK Canopy level CO2 pressure (Pa).
! GL       = WORK Leaf conductance for H2O (m/s).
! OA       = WORK Atmospheric O2 pressure (Pa).
! RD       = WORK Dark respiration of top leaf (mol CO2/m2/s).

      real CO2C(POINTS), DQC(POINTS), PSTAR(POINTS), IPAR(POINTS)
      real TSTAR(POINTS), FPAR(POINTS), FSMC(POINTS), LAI(POINTS)
      real ZERODEGC, EPCO2, ANETC(POINTS), CI(POINTS), GC(POINTS)
      real RDC(POINTS), ANETL(POINTS), APAR(POINTS), CA(POINTS)
      real GL(POINTS), OA(POINTS), RD(POINTS)

! Local parameters
! EPO2  = Ratio of molecular weights of O2 and dry air.
! O2    = Atmospheric concentration of oxygen (kg O2/kg air).

      real EPO2, O2
      parameter (EPO2=1.106, O2=0.23)

!-----------------------------------------------------------------------
! Calculate the atmospheric pressures of CO2 and O2 and the PAR absorbed
! by the top leaf
!-----------------------------------------------------------------------
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        CA(L) = CO2C(L) / EPCO2 * PSTAR(L)
        OA(L) = O2 / EPO2 * PSTAR(L)
        APAR(L) = (1 - OMEGA(FT)) * IPAR(L)
      enddo

!-----------------------------------------------------------------------
! Call the leaf level model for the top leaf of the C3 and C4 plants
!-----------------------------------------------------------------------
      call LEAF (LAND_PTS, LAND_INDEX, FT, ZERODEGC, DQC, APAR, TSTAR
     &,          CA, OA, PSTAR, FSMC, GL, ANETL, CI, RD)

!-----------------------------------------------------------------------
! Scale-up to the canopy level
!-----------------------------------------------------------------------
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        ANETC(L) = ANETL(L) * FPAR(L)
        GC(L) = FPAR(L) * GL(L)
        RDC(L) = RD(L) * FPAR(L)
      enddo

      return
      end

      subroutine LEAF (LAND_PTS, LAND_INDEX, FT, ZERODEGC, DQ, APAR, TL
     &,                CA, OA, PSTAR, FSMC, GL, AL, CI, RD)

!-----------------------------------------------------------------------
! Calculates the leaf resistance and net photosynthesis using:
!  (i) Collatz et al. (1992) C3/C4 photosynthesis model
! (ii) Jacobs (1994) CI/CA closure.

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
! FT         = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), FT, J, L

! ZERODEGC = IN Zero Celsius (K).
! DQ       = IN Canopy level specific humidity deficit (kg H2O/kg air).
! APAR     = IN Absorbed PAR (W/m2)
! TL       = IN Leaf temperature (K).
! CA       = IN Canopy CO2 pressure (Pa).
! OA       = IN Atmospheric O2 pressure (Pa).
! PSTAR    = IN Atmospheric pressure (Pa).
! FSMC     = IN Soil water factor.
! GL       = OUT Leaf conductance for H2O (m/s).
! AL       = OUT Net Leaf photosynthesis (mol CO2/m2/s).
! RD       = OUT Dark respiration (mol CO2/m2/s).
! CI       = OUT Internal CO2 pressure (Pa).
! ACR      = WORK Absorbed PAR (mol photons/m2/s).
! B1,B2,B3 = WORK Coefficients of the quadratic.
! CCP      = WORK Photorespiratory compensatory point (mol/m3).
! CONV     = WORK Factor for converting mol/m3 into Pa (J/mol).
! DENOM    = WORK Denominator in equation for VCM
! GLCO2    = WORK Leaf conductance for CO2 (m/s).
! KC       = WORK C3 Michaelis constant for CO2 (Pa)
! KO       = WORK C3 Michaelis constant for O2 (Pa).
! QTENF    = WORK Q10 function.
! TAU      = WORK CO2/O2 specificity ratio.
! TDEGC    = WORK Leaf temperature (deg C).
! VCM      = WORK Maximum rate of carboxylation of Rubisco
!            (mol CO2/m2/s).
! VCMAX    = WORK Maximum rate of carboxylation of Rubisco - without
!            the temperature factor (mol CO2/m2/s).
! WL       = WORK Gross leaf phtosynthesis (mol CO2/m2/s).
! WCARB    = WORK Carboxylation limited gross photosynthetic rate
!            (mol CO2/m2/s).
! WLITE    = WORK Light limited gross photosynthetic rate
!            (mol CO2/m2/s).
! WEXPT    = WORK export limited gross photosynthetic rate
!            (mol CO2/m2/s).
! WP       = WORK Smoothed minimum of Carboxylation and Light limited
!            gross photosynthesis (mol CO2/m2/s).

      real ZERODEGC, DQ(POINTS), APAR(POINTS), TL(POINTS), CA(POINTS)
      real OA(POINTS), PSTAR(POINTS), FSMC(POINTS), GL(POINTS)
      real AL(POINTS), RD(POINTS), CI(POINTS), ACR(POINTS), B1(POINTS)
      real B2(POINTS), B3(POINTS), CCP(POINTS), CONV(POINTS)
      real DENOM(POINTS), GLCO2(POINTS), KC(POINTS), KO(POINTS)
      real TAU(POINTS), TDEGC(POINTS), VCM(POINTS), VCMAX(POINTS)
      real QTENF(POINTS), WL(POINTS), WCARB(POINTS), WLITE(POINTS)
      real WEXPT(POINTS), WP(POINTS)

! CLOS_INDEX = WORK Index of land points with closed stomata.
! CLOS_PTS   = WORK Number of land points with closed stomata.
! OPEN_INDEX = WORK Index of land points with open stomata.
! OPEN_PTS   = WORK Number of land points with open stomata.

      integer CLOS_INDEX(POINTS), CLOS_PTS, OPEN_INDEX(POINTS), OPEN_PTS

! Local parameters
! R                     = Gas constant (J/K/mol).
! RATIO                 = Ratio of leaf resistance for CO2 to leaf
! BETA1                 = Coupling coefficient for co-limitation.
! BETA2                 = Coupling coefficient for co-limitation.
!                         resistance for H2O.
! FDC, FDC3, FDC4       = Dark respiration coefficient
! NEFFC, NEFFC3, NEFFC4 = Constant relating VCMAX and leaf N from
!                         Schulze et al. 1994 (AMAX = 0.4E-3 * NL,
!                         assuming dry matter is 40% carbon by mass)
!                         and Jacobs 1994: C3 : VCMAX = 2 * AMAX
!                                          C4 : VCMAX = AMAX (mol/m2/s)

      real R, RATIO, BETA1, BETA2, FDC, FDC3, FDC4, NEFFC, NEFFC3
      real NEFFC4
      parameter (R=8.3144, RATIO=1.6, BETA1=0.83, BETA2=0.93)
      parameter (FDC3=0.015, FDC4=0.025, NEFFC3=0.64E-3, NEFFC4=0.32E-3)

      if (C3(FT) .eq. 1) then
        FDC = FDC3
        NEFFC = NEFFC3
      else
        FDC = FDC4
        NEFFC = NEFFC4
      endif

!----------------------------------------------------------------------
! Initialize counters
!----------------------------------------------------------------------
      CLOS_PTS = 0
      OPEN_PTS = 0

      do J=1,LAND_PTS
        L = LAND_INDEX(J)

!----------------------------------------------------------------------
! Calculate the points with closed stomata
!----------------------------------------------------------------------
        if (FSMC(L).eq.0.0 .or. DQ(L).ge.DQCRIT(FT)
     &                     .or. APAR(L).eq.0.0) then
          CLOS_PTS = CLOS_PTS + 1
          CLOS_INDEX(CLOS_PTS) = J
        else
          OPEN_PTS = OPEN_PTS + 1
          OPEN_INDEX(OPEN_PTS) = J
        endif

!----------------------------------------------------------------------
! Calculate the factor for converting mol/m3 into Pa (J/m3).
!----------------------------------------------------------------------
        CONV(L) = R * TL(L)
      enddo

!----------------------------------------------------------------------
! Calculate the photosynthetic parameters
!----------------------------------------------------------------------
      if (C3(FT) .eq. 1) then
        do J=1,OPEN_PTS
          L = LAND_INDEX(OPEN_INDEX(J))
          VCMAX(L) = NEFFC * NL0(FT)
          TDEGC(L) = TL(L) - ZERODEGC
          TAU(L) = 2600.0 * (0.57 ** (0.1 * (TDEGC(L) - 25.0)))
          CCP(L) = 0.5 * OA(L) / TAU(L)
          QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
          DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &             * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
          VCM(L) = QTENF(L) / DENOM(L)
          RD(L) = FDC * QTENF(L)
        enddo
      else
        do J=1,OPEN_PTS
          L = LAND_INDEX(OPEN_INDEX(J))
          VCMAX(L) = NEFFC * NL0(FT)
          TDEGC(L) = TL(L) - ZERODEGC
          CCP(L) = 0.0
          QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
          DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &             * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
          VCM(L) = QTENF(L) / DENOM(L)
          RD(L) = FDC * QTENF(L)
        enddo
      endif

      do J=1,CLOS_PTS
        L = LAND_INDEX(CLOS_INDEX(J))
        VCMAX(L) = NEFFC * NL0(FT)
        TDEGC(L) = TL(L) - ZERODEGC
        QTENF(L) = VCMAX(L) * (2.0 ** (0.1 * (TDEGC(L) - 25.0)))
        DENOM(L) = (1 + EXP (0.3 * (TDEGC(L) - TUPP(FT))))
     &           * (1 + EXP (0.3 * (TLOW(FT) - TDEGC(L))))
        VCM(L) = QTENF(L) / DENOM(L)
        RD(L) = FDC * QTENF(L)
      enddo

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      do J=1,OPEN_PTS
        L = LAND_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the internal CO2 pressure (Jacobs, 1994).
!----------------------------------------------------------------------
        CI(L) = (CA(L) - CCP(L)) * F0(FT)
     &        * (1 - DQ(L) / DQCRIT(FT)) + CCP(L)

!----------------------------------------------------------------------
! Convert absorbed PAR into mol PAR photons/m2/s
!----------------------------------------------------------------------
        ACR(L) = APAR(L) / 2.19E5

      enddo

!----------------------------------------------------------------------
! Calculate the gross photosynthesis for RuBP-Carboxylase, Light and
! Export limited photosynthesis (Collatz et al., 1992).
!----------------------------------------------------------------------
      if (C3(FT) .eq. 1) then
        do J=1,OPEN_PTS
          L = LAND_INDEX(OPEN_INDEX(J))
          KC(L) = 30.0 * (2.1 ** (0.1 * (TDEGC(L) - 25.0)))
          KO(L) = 30000.0 * (1.2 ** (0.1 * (TDEGC(L) - 25.0)))
          WCARB(L) = VCM(L) * (CI(L) - CCP(L))
     &             / (CI(L) + KC(L) * (1. + OA(L) / KO(L)))
          WLITE(L) = ALPHA(FT) * ACR(L) * (CI(L) - CCP(L))
     &             / (CI(L) + 2 * CCP(L))
          WEXPT(L) = 0.5 * VCM(L)
        enddo
      else
        do J=1,OPEN_PTS
          L = LAND_INDEX(OPEN_INDEX(J))
          WCARB(L) = VCM(L)
          WLITE(L) = ALPHA(FT) * ACR(L)
          WEXPT(L) = 20000.0 * VCM(L) * CI(L) / PSTAR(L)
        enddo
      endif

!----------------------------------------------------------------------
! Calculate the co-limited rate of gross photosynthesis
!----------------------------------------------------------------------
      do J=1,OPEN_PTS
        L = LAND_INDEX(OPEN_INDEX(J))
        B1(L) = BETA1
        B2(L) = - (WCARB(L) + WLITE(L))
        B3(L) = WCARB(L) * WLITE(L)
        WP(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))
      enddo

      do J=1,OPEN_PTS
        L = LAND_INDEX(OPEN_INDEX(J))
        B1(L) = BETA2
        B2(L) = - (WP(L) + WEXPT(L))
        B3(L) = WP(L) * WEXPT(L)
        WL(L) = -B2(L)/(2*B1(L))
     &         - SQRT(B2(L)*B2(L)/(4*B1(L)*B1(L)) - B3(L)/B1(L))
      enddo

!----------------------------------------------------------------------
! Carry out calculations for points with open stomata
!----------------------------------------------------------------------
      do J=1,OPEN_PTS
        L = LAND_INDEX(OPEN_INDEX(J))

!----------------------------------------------------------------------
! Calculate the net rate of photosynthesis
!----------------------------------------------------------------------
        AL(L) = (WL(L) - RD(L)) * FSMC(L)

!----------------------------------------------------------------------
! Diagnose the leaf conductance
!----------------------------------------------------------------------
        GLCO2(L) = (AL(L) * CONV(L)) / (CA(L) - CI(L))
        GL(L) = RATIO * GLCO2(L)
      enddo

!----------------------------------------------------------------------
! Close stomata at points with negative or zero net photosynthesis
! or where the leaf resistance exceeds its maximum value.
!----------------------------------------------------------------------
      do J=1,OPEN_PTS
        L = LAND_INDEX(OPEN_INDEX(J))
        if (GL(L).le.GLMIN(FT) .or. AL(L).le.0.0) then
          GL(L) = GLMIN(FT)
          GLCO2(L) = GL(L) / RATIO
          AL(L) = -RD(L) * FSMC(L)
          CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
        endif
      enddo

!----------------------------------------------------------------------
! Define fluxes and conductances for points with closed stomata
!----------------------------------------------------------------------
      do J=1,CLOS_PTS
        L = LAND_INDEX(CLOS_INDEX(J))
        GL(L) = GLMIN(FT)
        GLCO2(L) = GL(L) / RATIO
        AL(L) = -RD(L) * FSMC(L)
        CI(L) = CA(L) - AL(L) * CONV(L) / GLCO2(L)
      enddo
#endif

      return
      end
