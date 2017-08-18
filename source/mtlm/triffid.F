      subroutine TRIFFID (LAND_PTS, LAND_INDEX, FORW, GAMMA, FRAC_VS
     &,                   FRAC_AGR, FRAC_MIN, FRAC_SEED, DENOM_MIN
     &,                   BF, G_LEAF, NPP, RESP_S, RESP_W, CS, FRAC
     &,                   HT, LAI, C_VEG, CV, LIT_C, LIT_C_T)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Simulates changes in vegetation structure, areal
! coverage and the carbon contents of vegetation and soil.
! can be used to advance these variables dynamically
! (GAMMA=1/TIMESTEP) or to iterate toward  equilibrium
! (GAMMA --> 0.0, FORW=1.0).

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

! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.

      integer LAND_PTS, LAND_INDEX(POINTS)
      integer L, N, T

! FORW      = IN Forward timestep weighting.
! FRAC_VS   = IN Total fraction of gridbox covered by veg or soil.
! GAMMA     = IN Inverse timestep (/360days).
! FRAC_AGR  = IN Fraction of agriculture.
! FRAC_MIN  = IN Minimum areal fraction for PFTs.
! FRAC_SEED = IN "Seed" fraction for PFTs.
! DENOM_MIN = IN Minimum value for the denominator of the update
!             equation. Ensures that gradient descent does not lead
!             to an unstable solution.
! BF        = Burn fraction
! G_LEAF    = IN Turnover rate for leaf and fine root biomass
!             (/360days).
! NPP       = INOUT Net primary productivity (kg C/m2/360days).
! RESP_S    = INOUT Soil respiration (kg C/m2/360days).
! RESP_W    = INOUT Wood maintenance respiration (kg C/m2/360days).
! CS        = INOUT Soil carbon (kg C/m2).
! FRAC      = INOUT Fractional cover of each Functional Type.
! HT        = INOUT Vegetation height (m).
! LAI       = INOUT Leaf area index.
! C_VEG     = OUT Total carbon content of the vegetation (kg C/m2).
! CV        = OUT Gridbox mean vegetation carbon (kg C/m2).
! LIT_C     = OUT Carbon Litter (kg C/m2/360days).
! LIT_C_T   = OUT Gridbox mean carbon litter (kg C/m2/360days).
! DCVEG     = WORK Change in vegetation carbon during the timestep
!             (kg C/m2/timestep).
! DFRAC     = WORK Change in areal fraction during the timestep
!             (/timestep).
! FRAC_FLUX = WORK PFT fraction to be used in the calculation of
!             the gridbox mean fluxes.
! LAI_BAL   = WORK Leaf area index in balanced growth state.
! LEAF      = WORK Leaf biomass (kg C/m2).
! PC_S      = WORK Net carbon flux available for spreading
!             (kg C/m2/yr).
! PHEN      = WORK Phenological state.
! ROOT      = WORK Root biomass (kg C/m2).
! WOOD      = WORK Woody biomass (kg C/m2).
! DFRAC_AGR = WORK Increment to areal fraction from agriculture.
!             (/timestep).
! DFA       = WORK Increment to areal fraction from agriculture that
!             is not burnt. (/timestep).

      real FORW, FRAC_VS(POINTS), GAMMA, FRAC_AGR(POINTS), FRAC_MIN
      real FRAC_SEED, DENOM_MIN, BF, G_LEAF(POINTS,NPFT)
      real NPP(POINTS,NPFT), RESP_S(POINTS), RESP_W(POINTS,NPFT)
      real CS(POINTS), FRAC(POINTS,NTYPE), HT(POINTS,NPFT)
      real LAI(POINTS,NPFT), C_VEG(POINTS,NPFT), CV(POINTS)
      real LIT_C(POINTS,NPFT), LIT_C_T(POINTS), DCVEG(POINTS,NPFT)
      real DFRAC(POINTS,NPFT),  FRAC_FLUX, LAI_BAL(POINTS,NPFT)
      real LEAF(POINTS,NPFT), PC_S(POINTS,NPFT), PHEN(POINTS,NPFT)
      real ROOT(POINTS,NPFT), WOOD(POINTS,NPFT)
      real DFRAC_AGR(POINTS,NPFT), DFA

!----------------------------------------------------------------------
! Loop through Functional Types
!----------------------------------------------------------------------
      do N=1,NPFT

!----------------------------------------------------------------------
! Loop through TRIFFID points
!----------------------------------------------------------------------
        do T=1,LAND_PTS
          L=LAND_INDEX(T)

!----------------------------------------------------------------------
! Diagnose the balanced-growth leaf area index and the associated leaf,
! wood, root and total vegetation carbon
!----------------------------------------------------------------------
          LAI_BAL(L,N) = (A_WS(N)*ETA_SL(N)*HT(L,N)
     &                   /A_WL(N))**(1.0/(B_WL(N)-1))
          LEAF(L,N) = SIGL(N)*LAI_BAL(L,N)
          ROOT(L,N) = LEAF(L,N)
          WOOD(L,N) = A_WL(N)*(LAI_BAL(L,N)**B_WL(N))
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)

!----------------------------------------------------------------------
! Diagnose the phenological state
!----------------------------------------------------------------------
          PHEN(L,N) = LAI(L,N)/LAI_BAL(L,N)
        enddo

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
        call VEGCARB (LAND_PTS, LAND_INDEX, N, FORW, GAMMA, DENOM_MIN
     &,               G_LEAF(1,N), NPP(1,N), RESP_W(1,N), LEAF(1,N)
     &,               ROOT(1,N), WOOD(1,N), DCVEG(1,N), PC_S(1,N))

      enddo

!-----------------------------------------------------------------------
! Diagnose the new value of Canopy Height, Leaf Area Index and Total
! Vegetation Carbon
!-----------------------------------------------------------------------
      do N=1,NPFT

        do T=1,LAND_PTS
          L=LAND_INDEX(T)

          HT(L,N) = WOOD(L,N)/(A_WS(N)*ETA_SL(N))
     &              *(A_WL(N)/WOOD(L,N))**(1.0/B_WL(N))
          LAI_BAL(L,N) = LEAF(L,N)/SIGL(N)
          LAI(L,N) = PHEN(L,N)*LAI_BAL(L,N)
          C_VEG(L,N) = LEAF(L,N) + ROOT(L,N) + WOOD(L,N)

        enddo

      enddo

!----------------------------------------------------------------------
! Update the areal coverage of each functional type
!----------------------------------------------------------------------
      call LOTKA (LAND_PTS, LAND_INDEX, C_VEG, FORW, FRAC_VS
     &,           FRAC_AGR, FRAC_MIN, FRAC_SEED, DENOM_MIN, GAMMA
     &,           LAI_BAL, PC_S, FRAC, DFRAC, DFRAC_AGR)

!----------------------------------------------------------------------
! Diagnose the litter fall from the carbon balance of each vegetation
! type
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        LIT_C_T(L) = 0.
        do N=1,NPFT
# if defined O_crop_data
          DFA = DFRAC_AGR(L,N)*(1. - BF)
          FRAC_FLUX = FRAC(L,N) - (1. - FORW)*(DFRAC(L,N) + DFA)
          LIT_C(L,N) = NPP(L,N) - GAMMA/FRAC_FLUX*(C_VEG(L,N)*FRAC(L,N)
     &               - (C_VEG(L,N) - DCVEG(L,N))*(FRAC(L,N)
     &               - DFRAC(L,N) - DFA))
# else
          FRAC_FLUX = FRAC(L,N) - (1. - FORW)*(DFRAC(L,N))
          LIT_C(L,N) = NPP(L,N) - GAMMA/FRAC_FLUX*(C_VEG(L,N)*FRAC(L,N)
     &               - (C_VEG(L,N) - DCVEG(L,N))*(FRAC(L,N)
     &               - DFRAC(L,N)))
# endif
          LIT_C_T(L) = LIT_C_T(L) + FRAC_FLUX*LIT_C(L,N)
        enddo
      enddo

!----------------------------------------------------------------------
! Call SOILCARB to update the soil carbon content
!----------------------------------------------------------------------
      call SOILCARB (POINTS, LAND_PTS, LAND_INDEX, FORW, GAMMA
     &,              DENOM_MIN, LIT_C_T, RESP_S, CS)

!----------------------------------------------------------------------
! Diagnose the gridbox mean vegetation carbon
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        CV(L) = 0.0
        do N=1,NPFT
          CV(L) = CV(L) + FRAC(L,N)*C_VEG(L,N)
        enddo
      enddo
#endif

      return
      end
