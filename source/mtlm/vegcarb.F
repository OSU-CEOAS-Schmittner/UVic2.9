       subroutine VEGCARB (LAND_PTS, LAND_INDEX, N, FORW, GAMMA
     &,                    DENOM_MIN, G_LEAF, NPP, RESP_W, LEAF, ROOT
     &,                    WOOD, DCVEG, PC_S)

#if defined O_mtlm
!----------------------------------------------------------------------
! Updates carbon contents of the vegetation.

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
!----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.

      integer LAND_PTS, LAND_INDEX(POINTS)
      integer L, N, T

! FORW        = IN Forward timestep weighting.
! GAMMA       = IN Inverse timestep (/360days).
! DENOM_MIN   = IN Minimum value for the denominator of the update
!               equation. Ensures that gradient descent does not lead
!               to an unstable solution.
! G_LEAF      = IN Turnover rate for leaf and fine root biomass
!               (/360days).
! NPP         = INOUT Net primary productivity (kg C/m2/360days).
! RESP_W      = INOUT Wood maintenance respiration (kg C/m2/360days).
! LEAF        = INOUT Leaf biomass (kg C/m2).
! ROOT        = INOUT Root biomass (kg C/m2).
! WOOD        = INOUT Woody biomass (kg C/m2).
! DCVEG       = OUT Change in vegetation carbon during the timestep
!               (kg C/m2/timestep).
! PC_S        = OUT Net carbon flux available for spreading
!               (kg C/m2/360days).
! DFPAR_DLAI  = WORK Rate of change of FPAR with leaf area index.
! DLAI        = WORK Increment to the leaf area index.
! DLAMG_DLAI  = WORK Required for calculation of equilibrium increments.
! DLIT_DLAI   = WORK Required for calculation of equilibrium increments.
! DNPP_DLAI   = WORK Rate of change of NPP with leaf area index
!               (kg C/m2/360days/LAI).
! DPC_DLAI    = WORK Rate of change of PC with leaf area index
!               (kg C/m2/360days/LAI).
! DPCG_DLAI   = WORK Rate of change of PC_G with leaf area index
!                (kg C/m2/360days/LAI).
! DRESPW_DLAI = WORK Rate of change of RESP_W with leaf area index
! FPAR        = WORK PAR interception factor.
! LAI         = WORK Leaf area index.
! LAMBDA_G    = WORK Fraction of NPP available for spreading.
! LIT_C_L     = WORK Local rate of Carbon Litter production
!               (kg C/m2/360days).
! PC          = WORK Net carbon flux available to vegetation
!               (kg C/m2/360days)
! PC_G        = WORK Net carbon flux available for growth
!               (kg C/m2/360days).

      real FORW, GAMMA, DENOM_MIN, G_LEAF(POINTS), NPP(POINTS)
      real RESP_W(POINTS), LEAF(POINTS), ROOT(POINTS), WOOD(POINTS)
      real DCVEG(POINTS), PC_S(POINTS), DFPAR_DLAI, DLAI, DLAMG_DLAI
      real DLIT_DLAI, DNPP_DLAI(POINTS), DPC_DLAI(POINTS)
      real DPCG_DLAI(POINTS), DRESPW_DLAI, FPAR, LAI(POINTS), LAMBDA_G
      real LIT_C_L(POINTS), PC(POINTS), PC_G(POINTS)

      do T=1,LAND_PTS
        L=LAND_INDEX(T)

        LAI(L) = LEAF(L)/SIGL(N)
!----------------------------------------------------------------------
! Calculate the local production rate for carbon litter
!----------------------------------------------------------------------
        LIT_C_L(L) = G_LEAF(L)*LEAF(L)+G_ROOT(N)*ROOT(L)
     &               + G_WOOD(N)*WOOD(L)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the vegetation
!----------------------------------------------------------------------
        PC(L) = NPP(L) - LIT_C_L(L)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
        DLIT_DLAI = (G_LEAF(L)*LEAF(L)+G_ROOT(N)*ROOT(L))/LAI(L)
     &            + B_WL(N)*G_WOOD(N)*WOOD(L)/LAI(L)

        FPAR = (1 - EXP(-KPAR(N)*LAI(L)))/KPAR(N)
        DFPAR_DLAI = EXP(-KPAR(N)*LAI(L))

        DNPP_DLAI(L) = NPP(L)*DFPAR_DLAI/FPAR
     &               + (1-R_GROW(N))*RESP_W(L)
     &               *(DFPAR_DLAI/FPAR-B_WL(N)/LAI(L))

        LAMBDA_G = 1 - (LAI(L) - LAI_MIN(N))
     &                /(LAI_MAX(N) - LAI_MIN(N))
        DLAMG_DLAI = -1.0/(LAI_MAX(N) - LAI_MIN(N))

        PC_G(L) = LAMBDA_G * NPP(L) - LIT_C_L(L)
        DPCG_DLAI(L) = LAMBDA_G*DNPP_DLAI(L)
     &               + DLAMG_DLAI*NPP(L)
     &               - DLIT_DLAI
        DPC_DLAI(L) = DNPP_DLAI(L) - DLIT_DLAI

      enddo

!----------------------------------------------------------------------
! Update vegetation carbon contents
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        DCVEG(L) = LEAF(L)+ROOT(L)+WOOD(L)
      enddo

      call GROWTH (LAND_PTS, LAND_INDEX, N, DPCG_DLAI, FORW
     &,            GAMMA, DENOM_MIN, PC_G, LEAF, ROOT, WOOD)

      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        DCVEG(L) = LEAF(L)+ROOT(L)+WOOD(L)-DCVEG(L)
      enddo

!----------------------------------------------------------------------
! Diagnose the carbon available for spreading and apply implicit
! corrections to the driving fluxes.
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        DLAI = LEAF(L)/SIGL(N) - LAI(L)
        PC_S(L) = PC(L) + FORW*DPC_DLAI(L)*DLAI - DCVEG(L)*GAMMA

        FPAR = (1 - EXP(-KPAR(N)*LAI(L)))/KPAR(N)
        DFPAR_DLAI = EXP(-KPAR(N)*LAI(L))
        DRESPW_DLAI = RESP_W(L)*B_WL(N)/LAI(L)

        NPP(L) = NPP(L) + FORW*DNPP_DLAI(L)*DLAI
        RESP_W(L) = RESP_W(L) + FORW*DRESPW_DLAI*DLAI
      enddo

      return
      end

      subroutine GROWTH (LAND_PTS, LAND_INDEX, N, DPCG_DLAI, FORW
     &,                  GAMMA, DENOM_MIN, PC_G, LEAF, ROOT,WOOD)

!----------------------------------------------------------------------
! Increments leaf, root and wood carbon.

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
!----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.
! N          = IN Plant functional type.

      integer LAND_PTS, LAND_INDEX(POINTS), N, I, L

! DPCG_DLAI  = IN Rate of change of PC_G with leaf area index
!              (kg C/m2/360days/LAI).
! FORW       = IN Forward timestep weighting.
! GAMMA      = IN Inverse timestep (/360days).
! DENOM_MIN  = IN Minimum value for the denominator of the update
!              equation. Ensures that gradient descent does not lead
!              to an unstable solution.
! PC_G       = IN Net carbon flux available for growth
!              (kg C/m2/360days).
! LEAF       = INOUT Leaf biomass (kg C/m2).
! ROOT       = INOUT Root biomass (kg C/m2).
! WOOD       = INOUT Woody biomass (kg C/m2).
! DENOM      = WORK Denominator of update equation.
! DLEAF      = WORK Increments to leaf biomass (kg C/m2).
! DROOT      = WORK Increments to root biomass (kg C/m2).
! DWOOD      = WORK Increments to woody biomass (kg C/m2).
! DL_DW      = WORK Rate of change of leaf carbon with wood carbon.
! DLAI_DW    = WORK Rate of change of leaf area index with wood carbon
!              (LAI m2/kg C).
! DR_DW      = WORK Rate of change of root carbon with wood carbon.
! NUMER      = WORK Numerator of the update equation.
! WOOD_MAX   = WORK Maximum wood carbon (kg C/m2).
! WOOD_MIN   = WORK Minimum wood carbon (kg C/m2).

      real DPCG_DLAI(POINTS), FORW, GAMMA, DENOM_MIN, PC_G(POINTS)
      real LEAF(POINTS), ROOT(POINTS), WOOD(POINTS), DENOM, DLEAF
      real DROOT, DWOOD, DL_DW, DLAI_DW, DR_DW, NUMER, WOOD_MAX
      real WOOD_MIN

      do I=1,LAND_PTS
        L=LAND_INDEX(I)

!----------------------------------------------------------------------
! Calculate the increment to the wood carbon
!----------------------------------------------------------------------
        DL_DW = LEAF(L)/(B_WL(N)*WOOD(L))
        DR_DW = DL_DW
        DLAI_DW = DL_DW/SIGL(N)

        NUMER = PC_G(L)
        DENOM = (1+DL_DW+DR_DW)*GAMMA-FORW*DLAI_DW*DPCG_DLAI(L)
        DENOM = MAX(DENOM,DENOM_MIN)

        DWOOD = NUMER/DENOM

!----------------------------------------------------------------------
! Ensure that the local leaf area index does not drop below its
! minimum value or exceed its maximum value.
!----------------------------------------------------------------------
        WOOD_MIN = A_WL(N)*LAI_MIN(N)**B_WL(N)
        WOOD_MAX = A_WL(N)*LAI_MAX(N)**B_WL(N)
        DWOOD = MAX((WOOD_MIN-WOOD(L)),DWOOD)
        DWOOD = MIN((WOOD_MAX-WOOD(L)),DWOOD)

!----------------------------------------------------------------------
! Diagnose the increments to leaf and root carbon
!----------------------------------------------------------------------
        DLEAF = SIGL(N)*((WOOD(L)+DWOOD)/A_WL(N))**(1.0/B_WL(N))
     &         -LEAF(L)
        DROOT = DLEAF

!----------------------------------------------------------------------
! Update carbon contents
!----------------------------------------------------------------------
        LEAF(L) = LEAF(L)+DLEAF
        ROOT(L) = ROOT(L)+DROOT
        WOOD(L) = WOOD(L)+DWOOD

      enddo
#endif

      return
      end
