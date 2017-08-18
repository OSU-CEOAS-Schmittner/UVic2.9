      subroutine MICROBE (POINTS, LAND_PTS, LAND_INDEX
     &,                   CS, STH_SOIL, V_SAT, V_WILT, TSOIL, RESP_S)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Calculates the soil respiration based on a simplified version of the
! model of Raich et al. (1991).

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

! CS         = IN Soil carbon (kg C/m2).
! STH_SOIL   = IN Top layer soil moisture as a fraction of saturation
!              (m3/m3).
! V_SAT      = IN Volumetric soil moisture concentration at saturation
!              (m3 H2O/m3 soil).
! V_WILT     = IN Volumetric soil moisture concentration below which
!              stomata close as a fraction of saturation
!              (m3 H2O/m3 soil).
! TSOIL      = IN Soil temperature (K).
! RESP_S     = OUT Soil respiration (kg C/m2/s).
! FSTH,FTEMP = WORK Factors describing the influence of soil moisture
!              and soil temperature respectively on soil respiration.
! STH_OPT    = WORK Fractional soil moisture at which respiration is
!              maximum.
! STH_WILT   = WORK Wilting soil moisture as a fraction of saturation.

      real CS(POINTS), STH_SOIL(POINTS), V_SAT(POINTS), V_WILT(POINTS)
      real TSOIL(POINTS), RESP_S(POINTS), FSTH, FTEMP, STH_OPT, STH_WILT

! Local parameters
! KAPS = Specific soil respiration rate at 25 deg and optimum soil
!        moisture (/s).
! Q10  = Q10 factor for soil respiration.

      real KAPS, Q10
      parameter (KAPS=0.35E-8, Q10=2.0)

      do I=1,LAND_PTS

        L = LAND_INDEX(I)
        if (V_SAT(L) .gt. 0.0) then

          STH_WILT = V_WILT(L) / V_SAT(L)
          STH_OPT = 0.5 * (1 + STH_WILT)

          if (STH_SOIL(L) .le. STH_WILT) then
            FSTH = 0.2
          elseif (STH_SOIL(L) .gt. STH_WILT .and.
     &            STH_SOIL(L) .le. STH_OPT) then
            FSTH = 0.2 + 0.8 * ((STH_SOIL(L) - STH_WILT)
     &                        / (STH_OPT - STH_WILT))
          elseif (STH_SOIL(L) .gt. STH_OPT) then
            FSTH = 1 - 0.8 * (STH_SOIL(L) - STH_OPT)
          endif

          FTEMP = Q10 ** (0.1 * (TSOIL(L) - 298.15))

          RESP_S(L) = KAPS * CS(L) * FSTH * FTEMP

        else

          RESP_S(L) = 0.0

        endif

      enddo
#endif

      return
      end
