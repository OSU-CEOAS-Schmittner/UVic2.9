      subroutine SOILCARB (POINTS, LAND_PTS, LAND_INDEX, FORW, GAMMA
     &,                    DENOM_MIN, LIT_C_T, RESP_S, CS)

#if defined O_mtlm
!----------------------------------------------------------------------
! Updates carbon contents of the soil.

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

! POINTS     = IN Total number of land points.
! LAND_PTS   = IN Number of points on which TRIFFID may operate.
! LAND_INDEX = IN Indices of land points on which TRIFFID may operate.

      integer POINTS, LAND_PTS, LAND_INDEX(POINTS), L, T

! FORW       = IN Forward timestep weighting.
! GAMMA      = IN Inverse timestep (/360days).
! DENOM_MIN  = IN Minimum value for the denominator of the update
!              equation. Ensures that gradient descent does not lead
!              to an unstable solution.
! LIT_C_T    = IN Total carbon litter (kg C/m2/360days).
! RESP_S     = INOUT Soil respiration (kg C/m2/360days).
! CS         = INOUT Soil carbon (kg C/m2).
! DCS        = WORK Increment to the soil carbon (kg C/m2).
! DPC_DCS    = WORK Rate of change of PC with soil carbon (/360days).
! PC         = WORK Net carbon accumulation in the soil
!              (kg C/m2/360days).

      real FORW, GAMMA, LIT_C_T(POINTS), DENOM_MIN, RESP_S(POINTS)
      real CS(POINTS), DCS(POINTS), DPC_DCS(POINTS), PC(POINTS)

      do T=1,LAND_PTS
        L=LAND_INDEX(T)

!----------------------------------------------------------------------
! Diagnose the net local carbon flux into the soil
!----------------------------------------------------------------------
        PC(L) = LIT_C_T(L)-RESP_S(L)

!----------------------------------------------------------------------
! Variables required for the implicit and equilibrium calculations
!----------------------------------------------------------------------
        DPC_DCS(L) = RESP_S(L)/CS(L)

!----------------------------------------------------------------------
! Save current value of soil carbon
!----------------------------------------------------------------------
        DCS(L) = CS(L)

      enddo

!----------------------------------------------------------------------
! Update soil carbon
!----------------------------------------------------------------------
      call DECAY (POINTS, LAND_PTS, LAND_INDEX, DPC_DCS, FORW, GAMMA
     &,           DENOM_MIN, PC, CS)

!----------------------------------------------------------------------
! Apply implicit correction to the soil respiration rate.
!----------------------------------------------------------------------
      do T=1,LAND_PTS
        L=LAND_INDEX(T)

        DCS(L) = CS(L) - DCS(L)
        RESP_S(L) = RESP_S(L) + FORW*DPC_DCS(L)*DCS(L)
# if defined O_mtlm_carbon_13 || O_mtlm_carbon_14
!       correct soil respiration
        RESP_S(L) = LIT_C_T(L) - DCS(L)*GAMMA
# endif

      enddo

      return
      end

      subroutine DECAY (POINTS, LAND_PTS, LAND_INDEX, DPC_DCS, FORW
     &,                 GAMMA, DENOM_MIN, PC, CS)

!-----------------------------------------------------------------------
! Updates carbon contents of the soil.

**********************************************************************
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

      integer POINTS, LAND_PTS, LAND_INDEX(POINTS), L, T

! DPC_DCS    = IN Rate of change of PC with soil carbon (yr).
! FORW       = IN Forward timestep weighting.
! GAMMA      = IN Inverse timestep (/360days).
! DENOM_MIN  = IN Minimum value for the denominator of the update
!              equation. Ensures that gradient descent does not lead
!              to an unstable solution.
! PC         = IN Net carbon flux into the soil (kg C/m2/360days).
! CS         = INOUT Soil carbon (kg C/m2).
! DENOM      = WORK Denominator of update equation.
! NUMER      = WORK Numerator of the update equation.
! CS_MIN     = Minimum soil carbon (kg C/m2).

      real DPC_DCS(POINTS), FORW, GAMMA, DENOM_MIN, PC(POINTS)
      real CS(POINTS), DENOM, NUMER, CS_MIN
      parameter (CS_MIN=1.0E-6)

      do T=1,LAND_PTS
        L=LAND_INDEX(T)
        NUMER = PC(L)
        DENOM = GAMMA+FORW*DPC_DCS(L)
        DENOM = MAX(DENOM,DENOM_MIN)
        CS(L) = CS(L)+NUMER/DENOM
        CS(L) = MAX(CS_MIN,CS(L))
      enddo
#endif

      return
      end
