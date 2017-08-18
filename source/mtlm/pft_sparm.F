      subroutine PFT_SPARM  (LAND_PTS, LAND_INDEX, N, ALBSOIL, HT, LAI
     &,                      ALBSNC_T,ALBSNF_T,CATCH_T, Z0_T)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Routine to calculate the land surface parameters of a given PFT from
! its areal fraction and structural properties.

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

! LAND_PTS   = IN Number of vegetated points.
! LAND_INDEX = IN Index of vegetated points.
! N         = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), N, J, L

! ALBSOIL  = IN Soil albedo.
! HT       = IN Vegetation height (m).
! LAI      = IN Leaf area index.
! ALBSNC_T = OUT Snow-covered albedo.
! ALBSNF_T = OUT Snow-free albedo.
! CATCH_T  = OUT Canopy capacity (kg/m2).
! Z0_T     = OUT Roughness length (m).
! FLIT     = WORK Weighting factor for albedo.

      real ALBSOIL(POINTS), HT(POINTS), LAI(POINTS), ALBSNC_T(POINTS)
      real ALBSNF_T(POINTS), CATCH_T(POINTS), Z0_T(POINTS), FLIT

!-----------------------------------------------------------------------
! Surface parameters for each Plant Functional Type
!-----------------------------------------------------------------------
! ALBSNC_MAX  = Snow-covered albedo for large LAI.
! ALBSNC_MIN  = Snow-covered albedo for zero LAI.
! ALBSNF_MAX  = Snow-free albedo for large LAI.
! DZ0V_DH     = Rate of change of vegetation roughness length with
!               height.
! CATCH0      = Minimum canopy capacity (kg/m2).
! DCATCH_DLAI = Rate of change of canopy capacity with LAI.
! INFIL_F     = Infiltration enhancement factor.
! KEXT        = Light extinction coefficient.
! ROOTD_FT    = Rootdepth (m).

      real ALBSNC_MAX(NPFT), ALBSNC_MIN(NPFT), ALBSNF_MAX(NPFT)
      real DZ0V_DH(NPFT), CATCH0(NPFT), DCATCH_DLAI(NPFT)
      real INFIL_F(NPFT), KEXT(NPFT), ROOTD_FT(NPFT)

!----------------------------------------------------------------------
!                           BT    NT   C3G   C4G    S
!----------------------------------------------------------------------
      DATA ALBSNC_MAX  /  0.35, 0.35, 0.60, 0.60, 0.40 /
      DATA ALBSNC_MIN  /  0.65, 0.65, 0.80, 0.80, 0.80 /
      DATA ALBSNF_MAX  /  0.16, 0.16, 0.20, 0.20, 0.20 /
      DATA DZ0V_DH     /  0.05, 0.05, 0.10, 0.10, 0.10 /
      DATA CATCH0      /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA DCATCH_DLAI /  0.05, 0.05, 0.05, 0.05, 0.05 /
      DATA INFIL_F     /  4.00, 4.00, 2.00, 2.00, 2.00 /
      DATA KEXT        /  0.50, 0.50, 0.50, 0.50, 0.50 /
      DATA ROOTD_FT    /  3.00, 1.00, 0.50, 0.50, 0.50 /

      do J=1,LAND_PTS
        L = LAND_INDEX(J)
        FLIT = 1.0 - EXP(-KEXT(N) * LAI(L))
        ALBSNC_T(L) = ALBSNC_MIN(N) * (1 - FLIT)
     &              + ALBSNC_MAX(N) * FLIT
        ALBSNF_T(L) = ALBSOIL(L) * (1 - FLIT) + ALBSNF_MAX(N) * FLIT
        Z0_T(L) = DZ0V_DH(N) * HT(L)
        CATCH_T(L) = CATCH0(N) + DCATCH_DLAI(N) * LAI(L)
      enddo
#endif

      return
      end
