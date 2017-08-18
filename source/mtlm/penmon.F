      subroutine PENMON (POINTS, LAND_PTS, LAND_INDEX, DZ_SOIL
     &,                  HCON_SOIL, RS, Z0, LW, SWN, PSTAR, Q1, T1
     &,                  TS1, WIND, Z1, LC, LF, EPSILON, SIGMA, TM
     &,                  ZERODEGC, E, LE, SH, G, RADNET, TSTAR, MNEG
     &,                  LW_OUT, LYING_SNOW)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Routine to calculate the evaporation using an extended version of
! the Penman-Monteith equation.

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

! POINTS     = IN Total number of land points.
! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.

      integer POINTS, LAND_PTS, LAND_INDEX(POINTS), I, L

! Surface variables
! DZ_SOIL     = IN Soil layer thickness (m).
! HCON_SOIL   = IN Soil heat capacity (W/m/K).
! RS          = IN Surface resistance (s/m).
! Z0          = IN Roughness length (m).
! MNEG        = IN Negative soil moisture (kg/m2).
! LYING_SNOW  = IN Snow (kg/m2).
! Driving variables
! LW          = IN Downward LW (W/m2).
! SWN         = IN Net downward SW (W/m2).
! PSTAR       = IN Surface pressure (Pa).
! Q1          = IN Specific humidity (kg/kg).
! T1          = IN Atmospheric temperature (K).
! TS1         = IN Sub-surface temperature (K).
! WIND        = IN Windspeed (m/s).
! Parameters
! Z1          = IN Reference height (m).
! LC          = IN Latent heat of condensation (J/kg).
! LF          = IN Latent heat of fusion (J/kg).
! EPSILON     = IN Ratio of molecular weights of water and dry air.
! SIGMA       = IN Stefan-Boltzman constant (W/m2/K4).
! TM          = IN Melting point of fresh water (K).
! ZERODEGC    = IN Zero Celsius (K).
! Outputs
! E           = OUT Evapotranspiration (kg/m2/s).
! LE          = OUT Latent heat flux (W/m2).
! SH          = OUT Sensible heat flux (W/m2).
! G           = OUT Ground heat flux (W/M2).
! RADNET      = OUT Surface net radiation (W/m2).
! TSTAR       = OUT Surface temperature (K).
! LW_OUT      = OUT net longwave radiation (W/m2).
! Work Variables
! AHAT        = WORK "Available energy" (W/m2).
! AS1         = WORK 2*HCON_SOIL/DZ_SOIL (W/m2/K).
! CHN         = WORK Neutral transfer coefficients.
! DENOM       = WORK Denominator of PM equation.
! DQ1         = WORK Humidity deficit (kg/kg).
! DQS_DT      = WORK Rate of change of saturated specific humidity with
!               temperature (kg/kg/K).
! DUM         = WORK Workspace variable.
! LAT         = WORK Latent heat constant (J/kg).
! NUMER       = WORK Numerator of PM equation.
! RESF        = WORK 1/(1+RS/RA).
! QS1         = WORK Saturated specific humidity at (T1,PSTAR) (kg/kg).
! RA          = WORK Aerodynamic resistance (s/m).
! RHOSTAR     = WORK Surface air density (kg/m3).
! ZETAH,ZETAM = WORK Temporaries in calculation of CHN.

      real DZ_SOIL ,HCON_SOIL, RS(POINTS), Z0(POINTS)
      real MNEG(POINTS), LYING_SNOW(POINTS), LW(POINTS)
      real SWN(POINTS), PSTAR(POINTS), Q1(POINTS), T1(POINTS)
      real TS1(POINTS), WIND(POINTS), Z1, LC, LF, EPSILON, SIGMA
      real TM, ZERODEGC, E(POINTS), LE(POINTS), SH(POINTS)
      real G(POINTS), RADNET(POINTS), TSTAR(POINTS)
      real LW_OUT(POINTS), AHAT(POINTS), AS1, CHN, DENOM
      real DQ1(POINTS), DQS_DT(POINTS), DUM, LAT(POINTS), NUMER
      real RESF, QS1(POINTS), RA(POINTS), RHOSTAR(POINTS), ZETAH
      real ZETAM

! Local parameters
! R        = Gas constant (J/kg/K).
! CP       = Specific heat of dry air at constant pressure (J/kg/K).
! KARMAN   = Von Karman's constant.

      real R, CP, KARMAN,KARMAN_SQ
      parameter (R=287.05, CP=1005.0, KARMAN=0.4, KARMAN_SQ=0.16 )

      AS1 = 2*HCON_SOIL/DZ_SOIL

!----------------------------------------------------------------------
! Calculate the surface air density (use atmospheric rather than
! surface temperature).
!----------------------------------------------------------------------
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        RHOSTAR(L) = PSTAR(L)/(R*T1(L))
      enddo

!----------------------------------------------------------------------
! Calculate the saturated specific humidity, its gradient w.r.t.
! temperature, and the humidity deficit.
!----------------------------------------------------------------------
      call QSAT (POINTS, LAND_PTS, LAND_INDEX, EPSILON, ZERODEGC
     &,          QS1, T1, PSTAR)

!CDIR NODEP
      do I=1,LAND_PTS
       L = LAND_INDEX(I)
       if (LYING_SNOW(L) .le. 50.) then
         LAT(L) = LC
       else
         LAT(L) = LC + LF
       endif
       DQS_DT(L) = (EPSILON*LAT(L)*QS1(L))/(R*T1(L)*T1(L))
       DQ1(L) = QS1(L) - Q1(L)

!-----------------------------------------------------------------------
! Calculate available energy when the surface temperature is equal to
! the atmospheric temperature (AHAT).
!-----------------------------------------------------------------------
        AHAT(L) = SWN(L) + LW(L) - SIGMA*T1(L)*T1(L)*T1(L)*T1(L)
     &                           - AS1*(T1(L)-TS1(L))

!-----------------------------------------------------------------------
! Calculate the neutral bulk transfer coefficient and aerodynamic
! resistance.
!-----------------------------------------------------------------------
        ZETAM = LOG((Z1 + Z0(L)) / Z0(L))
        ZETAH = LOG((Z1 + Z0(L)) / (0.1*Z0(L)))
        CHN = KARMAN_SQ / (ZETAH * ZETAM)
        RA(L) = 1.0 / (CHN * WIND(L))

!----------------------------------------------------------------------
! Calculate the evaporation rate and diagnose the surface temperature.
!----------------------------------------------------------------------
        RESF = 1.0/(1.0+RS(L)/RA(L))
        DUM = RHOSTAR(L)*CP/RA(L) + 4*SIGMA*T1(L)**3 + AS1
        NUMER = (DQS_DT(L)*AHAT(L) + DUM*DQ1(L)) * RESF
        DENOM = RESF*LAT(L)*DQS_DT(L)  + RA(L)*DUM/RHOSTAR(L)
        E(L) = NUMER / DENOM
        if (MNEG(L).lt.0.) E(L) = amin1(E(L),0.)
        LE(L) = LAT(L) * E(L)
        TSTAR(L) = T1(L) +
     &             (AHAT(L) - LAT(L)*RHOSTAR(L)*DQ1(L)*RESF/RA(L))
     &             /(DUM + DQS_DT(L)*LAT(L)*RHOSTAR(L)*RESF/RA(L))
        SH(L) = RHOSTAR(L)*CP/RA(L)*(TSTAR(L)-T1(L))
        LW_OUT(L) = LW(L) - SIGMA*TSTAR(L)*TSTAR(L)*TSTAR(L)*TSTAR(L)
        RADNET(L) = SWN(L) + LW_OUT(L)
        G(L) = RADNET(L) - LE(L) - SH(L)
      enddo
#endif

      return
      end
