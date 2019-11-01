      subroutine MTLM_STATE (POINTS, LAND_PTS, LAND_INDEX, DZ_SOIL
     &,                      HCAP_SOIL, KS, THETA_SAT, LF, TM, TIMESTEP
     &,                      G, RAIN, SNOW, E, ESUB, M, LYING_SNOW, TS1
     &,                      RUNOFF, SNOWMELT, MNEG)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Routine to update land surface prognostic variables
! (soil moisture, soil temperature, lying snow mass).
! Also diagnoses runoff and snowmelt

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

! POINTS     = IN Maximum number of land points.
! LAND_PTS   = IN Number of land points.
! LAND_INDEX = IN Indices of landpoints

      integer POINTS, LAND_PTS, LAND_INDEX(POINTS)
      integer I, L

! Surface parameters
! DZ_SOIL   = IN Soil layer thickness (m).
! HCAP_SOIL = IN Soil heat capacity (W/m3/K).
! TIMESTEP  = IN Timestep (s).
! KS        = IN Saturated hydraulic conductivity
! THETA_SAT = IN Saturated volumetric soil moisture concentration
!             (m3/m3).
! LF        = Latent heat of fusion (J/kg).
! TM        = Melting point of fresh water (K).
! G         = IN Ground heat flux (W/m2).

      real DZ_SOIL, HCAP_SOIL, TIMESTEP, KS, THETA_SAT, LF, TM
      real G(POINTS)

! Driving variables
! RAIN       = IN Rainfall (kg/m2/s).
! SNOW       = IN Snowfall (kg/m2/s).
! E          = IN Evapotranspiration (kg/m2/s).
! ESUB       = IN Sublimation (kg/m2/s).
! M          = INOUT Soil moisture (kg/m2).
! NEG        = INOUT Negative Soil moisture (kg/m2).
! LYING_SNOW = INOUT Lying snow (kg/m2).
! TS1        = INOUT Soil temperature (K).
! RUNOFF     = OUT Runoff (kg/m2/s).
! SNOWMELT   = OUT Snow melt (kg/m2/s).

      real RAIN(POINTS), SNOW(POINTS), E(POINTS), ESUB(POINTS)
      real M(POINTS), MNEG(POINTS), LYING_SNOW(POINTS)
      real TS1(POINTS), RUNOFF(POINTS), SNOWMELT(POINTS)

! Local parameters
! B          = Clapp-Hornberger exponent.

      real B
      parameter (B=6.6)

!----------------------------------------------------------------------
! Update the soil temperature and diagnose snowmelt.
!----------------------------------------------------------------------
!CDIR NODEP
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        TS1(L) = TS1(L) + TIMESTEP* G(L) / (DZ_SOIL*HCAP_SOIL)
        SNOWMELT(L) = 0.0
        if ((LYING_SNOW(L).gt.0.0) .and. (TS1(L).gt.TM)) then
          SNOWMELT(L) = HCAP_SOIL*DZ_SOIL*(TS1(L)-TM)
     &                / (LF*TIMESTEP)
          if ((SNOWMELT(L)).gt.(SNOW(L)-ESUB(L)+LYING_SNOW(L)/TIMESTEP))
     &      then
!           limit snowmelt to amount of snow and fix soil temperature
            SNOWMELT(L) = SNOW(L)-ESUB(L)+LYING_SNOW(L)/TIMESTEP
            TS1(L) = TS1(L) - SNOWMELT(L)*LF*TIMESTEP/HCAP_SOIL*DZ_SOIL
          else
            TS1(L) = TM
          endif
        endif

!----------------------------------------------------------------------
! Update the lying snow.
!----------------------------------------------------------------------
        LYING_SNOW(L) = LYING_SNOW(L)
     &                + TIMESTEP*(SNOW(L)-ESUB(L)-SNOWMELT(L))
        if (LYING_SNOW(L) .lt. 0.) then
!          if negative snow, change extra sublimation to evaporation
           ESUB(L) = ESUB(L) + LYING_SNOW(L)/TIMESTEP
           E(L) = E(L) - LYING_SNOW(L)/TIMESTEP
           TS1(L) = TS1(L) + LF*LYING_SNOW(L)/(DZ_SOIL*HCAP_SOIL)
           LYING_SNOW(L) = 0.
         endif

!----------------------------------------------------------------------
! Update the soil moisture.
!----------------------------------------------------------------------
        RUNOFF(L) = KS*(M(L)/(1000.0*DZ_SOIL*THETA_SAT))**(2*B+3)
        M(L) = M(L) + TIMESTEP*(RAIN(L)+SNOWMELT(L)-E(L)-RUNOFF(L))
!       keep track of any negative moisture for conservation
        if (M(L) + MNEG(L) .lt. 0.) then
          MNEG(L) = MNEG(L) + M(L)
          M(L) = 0.
        else
          M(L) = MNEG(L) + M(L)
          MNEG(L) = 0.
        endif
      enddo
#endif

      return
      end
