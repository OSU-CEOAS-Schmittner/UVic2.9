      subroutine SWRAD (POINTS, LAND_PTS, LAND_INDEX, ALBSNF, ALBSNC
     &,                LYING_SNOW, SW, TSTAR, TM, SWN, ALBLAND)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Calculate net shortwave radiation.

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

! ALBSNF      = IN Snow free albedo.
! ALBSNC      = IN Cold deep snow albedo.
! LYING_SNOW  = IN Lying snow (kg/m**2).
! SW          = IN Downward SW radiation at surface (W/m2).
! TSTAR       = IN Surface temperature (K).
! TM          = IN Melting point of fresh water (K).
! SWN         = OUT Net SW radiation at surface (W/m2).
! ALBLAND     = OUT Albedo
! ALBEDO      = Albedo.
! ALBSNOW     = Deep snow albedo.
! ALBSNOW_MIN = Minimum deep snow albedo.

      real ALBSNF(POINTS), ALBSNC(POINTS), LYING_SNOW(POINTS)
      real SW(POINTS), TSTAR(POINTS), TM, SWN(POINTS), ALBLAND(POINTS)
      real ALBEDO, ALBSNOW, ALBSNOW_MIN

! Local parameters
! DELT  = Temperature range below TM over which snow albedo is reduced.
! MASKD = Snow masking depth (kg/m**2).

      real DELT, MASKD
      parameter (DELT=2.0, MASKD=5.0)

      do L=1,LAND_PTS
        I = LAND_INDEX(L)
        ALBEDO = ALBSNF(I)
        ALBSNOW = ALBSNC(I)
        ALBSNOW_MIN = 0.7*ALBSNC(I) + 0.3*ALBSNF(I)
        if (LYING_SNOW(I) .gt. 0.0) then
          if (TSTAR(I).gt.(TM - DELT)) then
            ALBSNOW = ALBSNC(I) + 0.3*(ALBSNF(I) - ALBSNC(I))
     &                               *(TSTAR(I) - TM + DELT)/DELT
            ALBSNOW = MAX(ALBSNOW,ALBSNOW_MIN)
          endif
          ALBEDO = ALBSNF(I) + (ALBSNOW - ALBSNF(I))
     &                          *(1.0 - EXP(-LYING_SNOW(I)/MASKD))
        endif
        ALBLAND(I) = ALBEDO
!       incoming shortwave is already reduced by the albedo in the embm
        SWN(I) = SW(I)
      enddo
#endif

      return
      end
