      subroutine PHENOL (LAND_PTS, LAND_INDEX, N, G_LEAF, HT
     &,                  DTIME_PHEN, G_LEAF_PHEN, LAI)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Parametrizes leaf phenological changes and updates the leaf area
! index and the leaf turnover rate.

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
! N        = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), N, J, L

! G_LEAF      = IN Rate of leaf turnover (/360days).
! HT          = IN Canopy height (m).
! DTIME_PHEN  = IN Timestep (years).
! G_LEAF_PHEN = OUT Rate of leaf turnover including leaf phenology
!               (/360days).
! LAI         = INOUT Leaf area index.
! DPHEN       = WORK Increment to phenological state.
! LAI_BAL     = WORK Balanced growth LAI.
! PHEN        = WORK Phenological state.

      real G_LEAF(POINTS), HT(POINTS), DTIME_PHEN, G_LEAF_PHEN(POINTS)
      real LAI(POINTS), DPHEN, LAI_BAL(POINTS), PHEN(POINTS)

!-----------------------------------------------------------------------
! Diagnose the phenological state
!-----------------------------------------------------------------------
      do J=1,LAND_PTS
        L = LAND_INDEX(J)
        LAI_BAL(L) = (A_WS(N)*ETA_SL(N)*HT(L)
     &               /A_WL(N))**(1.0/(B_WL(N)-1))
        PHEN(L) = LAI(L)/LAI_BAL(L)
      enddo

!-----------------------------------------------------------------------
! Update the phenological state and output the leaf turnover rate in
! terms of the balanced growth LAI
!-----------------------------------------------------------------------
      do J=1,LAND_PTS
        L = LAND_INDEX(J)

        if (G_LEAF(L).gt.2*G_LEAF_0(N)) then
          DPHEN = -DTIME_PHEN*G_GROW(N)
          DPHEN = MAX(DPHEN,(0.01-PHEN(L)))
          G_LEAF_PHEN(L) = -DPHEN/DTIME_PHEN
        else
          DPHEN = DTIME_PHEN*G_GROW(N)*(1.0-PHEN(L))
          DPHEN = MIN(DPHEN,(1.0-PHEN(L)))
          G_LEAF_PHEN(L) = PHEN(L)*G_LEAF(L)
        endif

!-----------------------------------------------------------------------
! Update the leaf area index
!-----------------------------------------------------------------------
        PHEN(L) = PHEN(L) + DPHEN
        LAI(L) = PHEN(L)*LAI_BAL(L)

      enddo
#endif

      return
      end
