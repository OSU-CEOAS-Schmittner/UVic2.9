      subroutine LOTKA (LAND_PTS, LAND_INDEX, C_VEG, FORW, FRAC_VS
     &,                 FRAC_AGR, FRAC_MIN, FRAC_SEED, DENOM_MIN
     &,                 GAMMA, LAI, PC_S, FRAC, DFRAC, DFRAC_AGR)

#if defined O_mtlm
!-----------------------------------------------------------------------
! Updates fractional coverage of each functional type.
! Based on the Lotka-Volterra equations of interspecies
! competition.
!-----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.
! DOM        = IN Dominance hierachy.

      integer LAND_PTS, LAND_INDEX(POINTS), DOM(POINTS,NPFT)
      integer K, L, M, N, T

! C_VEG      = IN Carbon content of vegetation (kg C/m2).
! FORW       = IN Forward timestep weighting.
! FRAC_VS    = IN Total fractional cover of vegetation and soil.
! FRAC_AGR   = IN Fraction of agriculture.
! FRAC_MIN   = IN Minimum areal fraction for PFTs.
! FRAC_SEED  = IN "Seed" fraction for PFTs.
! DENOM_MIN  = IN Minimum value for the denominator of the update
!              equation. Ensures that gradient descent does not lead
!              to an unstable solution.
! GAMMA      = IN Inverse timestep (/360days).
! LAI        = IN Leaf area index.
! PC_S       = IN Net carbon flux available for spreading
!              (kg C/m2/360days).
! FRAC       = INOUT Fractional cover of each Functional Type.
! DFRAC      = OUT Increment to areal fraction during the timestep.
!              (/timestep).
! DFRAC_AGR  = OUT Increment to areal fraction from agriculture.
!              (/timestep).
! B          = WORK Mean rate of change of vegetation fraction over the
!              timestep (/360days).
! DB_DFRAC   = WORK Rate of change of B with vegetation fraction.
! COM        = WORK Coefficients representing the influence of one type
!              (second argument) on another (first argument).
! DIFF_SUM   = WORK Difference divided by sum for competing canopy
!              heights.
! HC1, HC2, HC3, HC4 = WORK Competing canopy heights (m).
! NOSOIL     = WORK Fractional area not available to vegetation.
! SPACE      = WORK Space available for invasion.

      real C_VEG(POINTS,NPFT), FORW, FRAC_VS(POINTS), FRAC_AGR(POINTS)
      real FRAC_MIN, FRAC_SEED, DENOM_MIN, GAMMA, LAI(POINTS,NPFT)
      real PC_S(POINTS,NPFT), FRAC(POINTS,NTYPE), DFRAC(POINTS,NPFT)
      real DFRAC_AGR(POINTS,NPFT), B(POINTS,NPFT)
      real DB_DFRAC(POINTS,NPFT,NPFT), COM(POINTS,NPFT,NPFT), DIFF_SUM
      real HC1, HC2, HC3, HC4, NOSOIL(POINTS), SPACE(POINTS,NPFT)

! Local parameters
! POW       = Power in sigmoidal function.

      real POW
      parameter (POW=20.)

!----------------------------------------------------------------------
! Define competition coefficients and the dominance hierachy
!----------------------------------------------------------------------
      do N=1,NPFT
        do M=1,NPFT
          do T=1,LAND_PTS
            L = LAND_INDEX(T)
            COM(L,N,M) = 1.
          enddo
        enddo
      enddo

      do T=1,LAND_PTS
        L = LAND_INDEX(T)

        HC1 = A_WL(1)/(A_WS(1)*ETA_SL(1))*(LAI(L,1)**(B_WL(1) - 1))
        HC2 = A_WL(2)/(A_WS(2)*ETA_SL(2))*(LAI(L,2)**(B_WL(2) - 1))
        DIFF_SUM = (HC1-HC2)/(HC1+HC2)

        COM(L,1,2) = 1./(1+EXP(POW*DIFF_SUM))    ! BT vs NT
        COM(L,1,3) = 0.                          ! BT vs C3G
        COM(L,1,4) = 0.                          ! BT vs C4G
        COM(L,1,5) = 0.                          ! BT vs S

        COM(L,2,1) = 1.-COM(L,1,2)               ! NT vs BT
        COM(L,2,3) = 0.                          ! NT vs C3G
        COM(L,2,4) = 0.                          ! NT vs C4G
        COM(L,2,5) = 0.                          ! NT vs S

        HC3 = A_WL(3)/(A_WS(3)*ETA_SL(3))*(LAI(L,3)**(B_WL(3) - 1))
        HC4 = A_WL(4)/(A_WS(4)*ETA_SL(4))*(LAI(L,4)**(B_WL(4) - 1))
        DIFF_SUM = (HC3-HC4)/(HC3+HC4)

        COM(L,3,4) = 1./(1+EXP(POW*DIFF_SUM))    ! C3G vs C4G
        COM(L,4,3) = 1.-COM(L,3,4)               ! C4G vs C3G

        COM(L,5,3) = 0.                          ! S vs C3G
        COM(L,5,4) = 0.                          ! S vs C4G

        if (HC1 .ge. HC2) then
          DOM(L,1) = 1
          DOM(L,2) = 2
        elseif (HC1 .lt. HC2) then
          DOM(L,1) = 2
          DOM(L,2) = 1
        endif

        DOM(L,3) = 5

        if (HC3 .ge. HC4) then
          DOM(L,4) = 3
          DOM(L,5) = 4
        elseif (HC3 .lt. HC4) then
          DOM(L,4) = 4
          DOM(L,5) = 3
        endif

      enddo

!----------------------------------------------------------------------
! Calculate the space available for the expansion of each FT
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L = LAND_INDEX(T)
        NOSOIL(L) = 1. - FRAC_VS(L)
      enddo
      do K=1,NPFT
        do T=1,LAND_PTS
          L = LAND_INDEX(T)
          N = DOM(L,K)
          SPACE(L,N) = 1. - NOSOIL(L) - FRAC_MIN*(NPFT-K)
        enddo
      enddo
      do N=1,NPFT
        do M=1,NPFT
          do T=1,LAND_PTS
            L = LAND_INDEX(T)
            SPACE(L,N) = SPACE(L,N) - COM(L,N,M)*FRAC(L,M)
          enddo
        enddo
      enddo

!----------------------------------------------------------------------
! Calculate the variables required for the implicit calculation.
! Divide the update equation by FRAC to eliminate the (unstable)
! bare soil solution.
!----------------------------------------------------------------------
      do N=1,NPFT
        do T=1,LAND_PTS
          L = LAND_INDEX(T)
          B(L,N) = PC_S(L,N)*SPACE(L,N)/C_VEG(L,N) - G_AREA(N)
          do M=1,NPFT
            DB_DFRAC(L,N,M) = -COM(L,N,M)*PC_S(L,N)/C_VEG(L,N)
          enddo
        enddo
      enddo

!----------------------------------------------------------------------
! Update the areal fractions
!----------------------------------------------------------------------
      call COMPETE (LAND_PTS, LAND_INDEX, DOM, C_VEG, B, DB_DFRAC
     &,             FORW, GAMMA, NOSOIL, FRAC_MIN, FRAC_SEED
     &,             DENOM_MIN, FRAC_AGR, FRAC, DFRAC, DFRAC_AGR)

      return
      end

      subroutine COMPETE (LAND_PTS, LAND_INDEX, DOM, C_VEG, B, DB_DFRAC
     &,                   FORW, GAMMA, NOSOIL, FRAC_MIN, FRAC_SEED
     &,                   DENOM_MIN, FRAC_AGR, FRAC, DFRAC, DFRAC_AGR)

!-----------------------------------------------------------------------
! Updates fractional coverage of each functional type.
! Requires a dominance hierachy as input.
!-----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "mtlm_data.h"

! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.
! DOM        = IN Dominance hierachy.

      integer LAND_PTS, LAND_INDEX(POINTS), DOM(POINTS,NPFT)
      integer K, L, M, N, T

! C_VEG             = IN Carbon content of vegetation (kg C/m2).
! B                 = IN Mean rate of change of vegetation fraction
!                     over the timestep (/360days).
! DB_DFRAC          = IN Rate of change of B with vegetation fraction.
! FORW              = IN Forward weighting factor.
! GAMMA             = IN Inverse timestep (/360days).
! NOSOIL            = IN Fractional area not available to vegetation.
! FRAC_MIN          = IN Minimum areal fraction for PFTs.
! FRAC_SEED         = IN "Seed" fraction for PFTs.
! DENOM_MIN         = IN Minimum value for the denominator of the update
!                     equation. Ensures that gradient descent does not
!                     lead to an unstable solution.
! FRAC_AGR          = IN Fraction of agriculture.
! FRAC              = INOUT Updated areal fraction.
! DFRAC             = OUT Increment to areal fraction.
! DFRAC_AGR         = OUT Increment to areal fraction from agriculture.
! DENOM             = WORK Denominator of update equation.
! FRACN,FRACM       = WORK Fractions used in the spreading calculation.
! NUMER             = WORK Numerator of the update equation.
! SPACE             = WORK Available space.
! P1,P2,Q1,Q2,R1,R2 = WORK Coefficients in simultaneous equations.
! ARF               = WORK variable for imposing agriculture fraction.

      real C_VEG(POINTS,NPFT), B(POINTS,NPFT)
      real DB_DFRAC(POINTS,NPFT,NPFT), FORW, GAMMA, NOSOIL(POINTS)
      real FRAC_MIN, FRAC_SEED, DENOM_MIN, FRAC_AGR(POINTS)
      real FRAC(POINTS,NTYPE), DFRAC(POINTS,NPFT)
      real DFRAC_AGR(POINTS,NPFT), DENOM, FRACN, FRACM, NUMER
      real SPACE(POINTS), P1, P2, Q1, Q2, R1, R2, ARF

!----------------------------------------------------------------------
! Initializations. Set increments to zero and define the space
! available to the dominant type leaving space for the seeds of others.
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        do N=1,NPFT
          DFRAC(L,N) = 0.0
        enddo
        SPACE(L) = 1. - NOSOIL(L) - FRAC_MIN*(NPFT - 1)

      enddo

# if defined O_crop_data || defined O_pasture_data || defined O_agric_data
!----------------------------------------------------------------------
! Set specified agriculture area (by reducing non-agriculture fraction)
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        ARF = 1.e-20
         do N=1,NPFT
          if (AGR(N) .eq. 0) ARF = ARF + FRAC(L,N) - FRAC_MIN
        enddo

        ARF = 1.-(1.-NOSOIL(L)-FRAC_MIN*(NPFT-1))*(1.-FRAC_AGR(L))/ARF
        if (ARF .gt. 1.) then
          ARF = 1.
        elseif (ARF .lt. 0.) then
          ARF = 0.
        endif

        do N=1,NPFT
          if (AGR(N) .eq. 0) then
            DFRAC_AGR(L,N) = (FRAC_MIN - FRAC(L,N))*ARF
            FRAC(L,N) = FRAC(L,N) + DFRAC_AGR(L,N)
          else
            DFRAC_AGR(L,N) = 0.
          endif
        enddo

      enddo

# endif
!----------------------------------------------------------------------
! Calculate the increments to the tree fractions
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        N = DOM(L,1)
        M = DOM(L,2)
        FRACN = FRAC(L,N)
        FRACN = MAX(FRACN,FRAC_SEED)
        FRACM = FRAC(L,M)
        FRACM = MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN - FORW*DB_DFRAC(L,N,N)
        P2 = GAMMA/FRACM - FORW*DB_DFRAC(L,M,M)
        Q1 = -FORW*DB_DFRAC(L,N,M)
        Q2 = -FORW*DB_DFRAC(L,M,N)
        R1 = B(L,N)
        R2 = B(L,M)
        do K=1,NPFT
          R1 = R1 + FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
          R2 = R2 + FORW*(DB_DFRAC(L,M,K)*DFRAC(L,K))
        enddo

        NUMER = R1 - (Q1/P2)*R2
        DENOM = P1 - (Q1/P2)*Q2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N) + DFRAC(L,N)
        if (FRAC(L,N) .lt. FRAC_MIN) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + FRAC_MIN
          FRAC(L,N) = FRAC_MIN
        elseif (FRAC(L,N) .gt. SPACE(L)) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + SPACE(L)
          FRAC(L,N) = SPACE(L)
        endif

        SPACE(L) = SPACE(L) - FRAC(L,N) + FRAC_MIN

        NUMER = R2 - Q2*DFRAC(L,N)
        DENOM = P2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,M) = NUMER/DENOM
        FRAC(L,M) = FRAC(L,M) + DFRAC(L,M)
        if (FRAC(L,M) .lt. FRAC_MIN) then
          DFRAC(L,M) = DFRAC(L,M) - FRAC(L,M) + FRAC_MIN
          FRAC(L,M) = FRAC_MIN
        elseif (FRAC(L,M) .gt .SPACE(L)) then
          DFRAC(L,M) = DFRAC(L,M) - FRAC(L,M) + SPACE(L)
          FRAC(L,M) = SPACE(L)
        endif

        SPACE(L) = SPACE(L) - FRAC(L,M) + FRAC_MIN

      enddo

!----------------------------------------------------------------------
! Calculate the increment to the shrub fraction
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        N = DOM(L,3)
        FRACN = FRAC(L,N)
        FRACN = MAX(FRACN,FRAC_SEED)

        DENOM = GAMMA/FRACN - FORW*DB_DFRAC(L,N,N)
        DENOM = MAX(DENOM,DENOM_MIN)

        NUMER = B(L,N)
        do K=1,NPFT
          NUMER = NUMER + FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
        enddo

        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N) + DFRAC(L,N)
        if (FRAC(L,N) .lt. FRAC_MIN) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + FRAC_MIN
          FRAC(L,N) = FRAC_MIN
        elseif (FRAC(L,N) .gt. SPACE(L)) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + SPACE(L)
          FRAC(L,N) = SPACE(L)
        endif

        SPACE(L) = SPACE(L) - FRAC(L,N) + FRAC_MIN

      enddo

!----------------------------------------------------------------------
! Calculate the increments to the grass fractions
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        N = DOM(L,4)
        M = DOM(L,5)
        FRACN = FRAC(L,N)
        FRACN = MAX(FRACN,FRAC_SEED)
        FRACM = FRAC(L,M)
        FRACM = MAX(FRACM,FRAC_SEED)

        P1 = GAMMA/FRACN-FORW*DB_DFRAC(L,N,N)
        P2 = GAMMA/FRACM-FORW*DB_DFRAC(L,M,M)
        Q1 = -FORW*DB_DFRAC(L,N,M)
        Q2 = -FORW*DB_DFRAC(L,M,N)
        R1 = B(L,N)
        R2 = B(L,M)
        do K=1,NPFT
          R1 = R1 + FORW*(DB_DFRAC(L,N,K)*DFRAC(L,K))
          R2 = R2 + FORW*(DB_DFRAC(L,M,K)*DFRAC(L,K))
        enddo

        NUMER = R1 - (Q1/P2)*R2
        DENOM = P1 - (Q1/P2)*Q2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,N) = NUMER/DENOM
        FRAC(L,N) = FRAC(L,N) + DFRAC(L,N)
        if (FRAC(L,N) .lt. FRAC_MIN) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + FRAC_MIN
          FRAC(L,N) = FRAC_MIN
        elseif (FRAC(L,N) .gt. SPACE(L)) then
          DFRAC(L,N) = DFRAC(L,N) - FRAC(L,N) + SPACE(L)
          FRAC(L,N) = SPACE(L)
        endif

        SPACE(L) = SPACE(L) - FRAC(L,N) + FRAC_MIN

        NUMER = R2 - Q2*DFRAC(L,N)
        DENOM = P2
        DENOM = MAX(DENOM,DENOM_MIN)
        DFRAC(L,M) = NUMER/DENOM
        FRAC(L,M) = FRAC(L,M) + DFRAC(L,M)
        if (FRAC(L,M) .lt. FRAC_MIN) then
          DFRAC(L,M) = DFRAC(L,M) - FRAC(L,M) + FRAC_MIN
          FRAC(L,M) = FRAC_MIN
        elseif (FRAC(L,M) .gt. SPACE(L)) then
          DFRAC(L,M) = DFRAC(L,M) - FRAC(L,M) + SPACE(L)
          FRAC(L,M) = SPACE(L)
        endif

        SPACE(L) = SPACE(L) - FRAC(L,M) + FRAC_MIN

      enddo

# if defined O_crop_data || defined O_pasture_data || defined O_agric_data
!----------------------------------------------------------------------
! Set specified agricultural area (by removing spread)
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        ARF = 1.e-20
        SPACE(L) = 1. - NOSOIL(L) - FRAC_MIN*(NPFT - 1)
         do N=1,NPFT
          if (AGR(N) .eq. 0) ARF = ARF + FRAC(L,N) - FRAC_MIN
        enddo

        ARF = 1.-(1.-NOSOIL(L)-FRAC_MIN*(NPFT-1))*(1.-FRAC_AGR(L))/ARF
        if (ARF .gt. 1.) then
          ARF = 1.
        elseif (ARF .lt. 0.) then
          ARF = 0.
        endif

        do N=1,NPFT
          if (AGR(N) .eq. 0) then
            DFRAC_AGR(L,N) = DFRAC_AGR(L,N) + (FRAC_MIN - FRAC(L,N))*ARF
            FRAC(L,N) = FRAC(L,N) + (FRAC_MIN - FRAC(L,N))*ARF
          endif
        enddo

      enddo

# endif
!----------------------------------------------------------------------
! Diagnose the new bare soil fraction
!----------------------------------------------------------------------
      do T=1,LAND_PTS

        L = LAND_INDEX(T)
        FRAC(L,SOIL) = 1. - NOSOIL(L)
        do N=1,NPFT
          FRAC(L,SOIL) = FRAC(L,SOIL) - FRAC(L,N)
        enddo

      enddo
#endif

      return
      end
