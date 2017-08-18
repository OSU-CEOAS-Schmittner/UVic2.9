      subroutine LEAF_LIT (LAND_PTS, LAND_INDEX, N, FSMC, TSTAR, G_LEAF)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Calculates the leaf turnover rate as a function of temperature and
! soil water availability
!-----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of vegetated points.
! LAND_INDEX = IN Index of vegetated points.
! N        = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), N, J, L

! FSMC   = IN Soil moisture availability factor.
! TSTAR  = IN Surface temperature (K).
! G_LEAF = OUT Rate of leaf turnover (/360days).
! FM,FT  = WORK Soil moisture and leaf temperature amplifiers of leaf
!          turnover.

      real FSMC(POINTS), TSTAR(POINTS), G_LEAF(POINTS), FM, FT

!-----------------------------------------------------------------------
! Calculate the leaf turnover rate
!-----------------------------------------------------------------------
      do J=1,LAND_PTS
        L = LAND_INDEX(J)
        FT = 1.0
        FM = 1.0
        if (TSTAR(L) .lt. TLEAF_OF(N)) then
          FT = 1.0 + DGL_DT(N)*(TLEAF_OF(N)-TSTAR(L))
        elseif (FSMC(L) .lt. FSMC_OF(N)) then
          FM = 1.0 + DGL_DM(N)*(FSMC_OF(N)-FSMC(L))
        endif
        G_LEAF(L) = G_LEAF_0(N)*FT*FM
      enddo
#endif

      return
      end
