      subroutine MTLM (is, ie, js, je)

#if defined O_mtlm
!-----------------------------------------------------------------------
! land-surface and vegetation model
! The land model is a modified form of the MOSES and TRIFFID

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
      include "csbc.h"
      include "mtlm.h"
# if defined O_mtlm_carbon_13
      include "mtlmc13.h"
# endif
# if defined O_mtlm_carbon_14
      include "mtlmc14.h"
# endif
      include "mtlm_data.h"
      include "switch.h"

      integer ie, is, je, js, NOTCROP(NPFT)
      integer I, J, L, N, KITER             ! Loop counters

!-----------------------------------------------------------------------
! Driving variables
!-----------------------------------------------------------------------
! LW_OUT_P = WORK Surface longwave (W/m**2)
! LW_OUT_S = WORK Surface longwave (W/m**2)
! PAR      = WORK Photosynthetically active
! SWN_P    = WORK Net shortwave for each PFT
! SWN_S    = WORK Net shortwave for bare soil
! RTSB     = WORK Inverse of the timestep for burning

      real LW_OUT_P(POINTS,NPFT), LW_OUT_S(POINTS), PAR(POINTS)
      real SWN_P(POINTS,NPFT), SWN_S(POINTS)
      real RTSB

!-----------------------------------------------------------------------
! Other variables which need to be saved between calls
!-----------------------------------------------------------------------
! ALBLAND_SOIL = WORK surface albedo
! ALBLAND_VEG  = WORK surface albedo

      real ALBLAND_SOIL(POINTS), ALBLAND_VEG(POINTS,NPFT)

!-----------------------------------------------------------------------
! Work Carbon variables.
!-----------------------------------------------------------------------
! FTIME_PHEN = Weighting factor for accumulations.
! GC         = Canopy conductance (m/s). phenology (/yr).
! LAI_BAL    = Balanced growth LAI.
! LIT_C      = Carbon litter (kg C/m2/yr).
! RESP_P     = Plant respiration rate (kg C/m2/s).
! RESP_W     = Wood maintenance respiration rate (kg C/m2/s).

      real FTIME_PHEN, GC(POINTS,NPFT), LAI_BAL(POINTS,NPFT)
      real LIT_C(POINTS,NPFT), RESP_P(POINTS,NPFT), RESP_W(POINTS,NPFT)

!-----------------------------------------------------------------------
! Work Hydrology variables
!-----------------------------------------------------------------------
! ESOIL      = Soil evaporation (kg/m2/s).
! ETRAN      = Transpiration (kg/m2/s).
! G_T        = ground heat flux on veg tiles (W/m2).
! G          = ground heat flux (W/m2).
! GSOIL      = Bare soil ground heat flux (W/m2).
! SH_T       = Sensible heat flux on veg tiles (W/m2
! HSOIL      = Bare soil sensible heat flux (W/m2).
! LESOIL     = Latent heat of soil evaporation (W/m2).
! LETRAN     = Latent heat of transpiration (W/m2).
! RADNET_T   = Net radiation on veg tiles (W/m2).
! RADNET     = Net Radiation (W/m2).
! RADNETSOIL = Bare soil net radiation (W/m2).
! STHU       = Unfrozen soil moisture as a fraction of saturation.
! V_ROOT     = Liquid volumetric soil moisture concentration in the
!              rootzone (m3 H2O/m3 soil).

      real ESOIL(POINTS), ETRAN(POINTS,NPFT), G_T(POINTS,NPFT)
      real G(POINTS), GSOIL(POINTS), SH_T(POINTS,NPFT), HSOIL(POINTS)
      real LESOIL(POINTS), LETRAN(POINTS,NPFT), RADNET_T(POINTS,NPFT)
      real RADNET(POINTS), RADNETSOIL(POINTS), STHU(POINTS)
      real V_ROOT(POINTS)

!-----------------------------------------------------------------------
! Other surface quantities
!-----------------------------------------------------------------------
! RA     = Aerodynamic resistance (s/m).
! RS     = Surface resistance (s/m).
! RSSOIL = Surface resistance for bare soil (s/m).
! V_SAT  = Volumetric soil moisture concentration at saturation
!          (m3 H2O/m3 soil).
! V_WILT = Volumetric soil moisture concentration below which stomata
!          close (m3 H2O/m3 soil).

      real RA(POINTS), RS(POINTS,NPFT), RSSOIL(POINTS), V_SAT(POINTS)
      real V_WILT(POINTS)

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
! ITER_EQ = Number of TRIFFID iterations for gradient descent to
!           equilibrium.

      integer ITER_EQ
      parameter (ITER_EQ=10)

! DENOM_MIN  = Minimum value for the denominator of the update equation.
!              Ensures that gradient descent does not lead to an
!              unstable solution.
! GAMMA_EQ   = Inverse timestep for gradient descent to equilibrium
!              (/yr).
! CP         = Specific heat of dry air at constant pressure (J/kg/K).
! EPSILON    = Ratio of molecular weights of water and dry air.
! KARMAN     = Von Karman's constant.
! TM         = Melting point of fresh water (K).
! ZERODEGC   = Zero Celsius (K).
! MSAT       = Saturated soil moisture (kg/m2).
! RSS        = Surface resistance for bare soil (s/m
! SATCON     = Saturated hydraulic conductivity (kg/m2/s).
! V_CRIT     = Volumetric soil moisture concentration above which
!              stomata are not sensitive to soil water (m3 H2O/m3 soil).
! VSAT       = Volumetric soil moisture concentration at saturation
!              (m3 H2O/m3 soil).
! VWILT      = Volumetric soil moisture concentration below which
!              stomata close (m3 H2O/m3 soil).
! Z0_SOIL    = Roughness length for bare soil (m).
! Z1_REF     = Reference level height (m).
! DTIME_PHEN = Timestep for phenology (/yr).
! FORW       = Forward timestep weighting.
! GAMMA      = Inverse timestep (/yr).
! RHO_WATER  = density of pure water (Kg/m3)
! PI         = PI

      real DENOM_MIN ,GAMMA_EQ, CP, EPSILON, KARMAN, TM, ZERODEGC
      real MSAT, RSS, SATCON, V_CRIT, VSAT, VWILT, Z0_SOIL, Z1_REF
      real DTIME_PHEN, FORW, GAMMA, RHO_WATER, PI

      parameter (DENOM_MIN=1.0E-6, GAMMA_EQ=1.0E-1, CP=1005.0)
      parameter (EPSILON=0.62198, KARMAN=0.4, TM=273.15)
      parameter (ZERODEGC=273.15, RSS=100.0, SATCON=0.0005)
      parameter (V_CRIT=0.34, VSAT=0.458, VWILT=0.13)
      parameter (MSAT=1000*ROOTDEP*VSAT, Z0_SOIL=0.0003, Z1_REF=10.0)
      parameter (RHO_WATER = 1000.0, PI=3.14159265358979323846)

# if defined O_mtlm_carbon_13
      real B13NPP, G13NPP(NPFT), B13RES, G13RES
# endif
# if defined O_mtlm_carbon_14
      real B14NPP, G14NPP(NPFT), B14RES, G14RES
# endif

!-----------------------------------------------------------------------
! Factors for accumulation variables
!-----------------------------------------------------------------------
      FTIME = TIMESTEP/REAL(SEC_DAY*DAY_TRIF)
      FTIME_PHEN = TIMESTEP/REAL(SEC_DAY*DAY_PHEN)

!-----------------------------------------------------------------------
! Calculate timestep values of the driving data.
!-----------------------------------------------------------------------
      LAND_COUNTER = LAND_COUNTER + 1
      ISTEP = mod(LAND_COUNTER-1,STEP_DAY) + 1

      L_TRIF = .false.
      L_PHEN = .false.
      if (ISTEP .eq. STEP_DAY) then
        L = LAND_COUNTER/STEP_DAY
        if (INT_VEG .and. MOD(L,DAY_TRIF) .eq. 0.) L_TRIF=.true.
        if (MOD(L,DAY_PHEN) .eq. 0.) L_PHEN=.true.
      endif

      TIMEDAY = (REAL(ISTEP) - 0.5)*SEC_DAY/STEP_DAY

      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        SAT_D(L) = T_C(L) + 0.5*DTEMP_DAY(L)
     &         *COS(2*PI*(TIMEDAY-TIME_MAX(L))/SEC_DAY)
        LW(L) = SIGMA*T_C(L)**4*(4.0*SAT_D(L)/T_C(L)-3.0)
        SW(L) = SUN(L,ISTEP)
      enddo

      call QSAT (POINTS, LAND_PTS, LAND_INDEX, EPSILON, ZERODEGC
     &,          QS, SAT_D, PSTAR)

!CDIR NODEP
      do I=1,LAND_PTS

        L = LAND_INDEX(I)
        if (RH_C(L) .gt. 0.) then
          Q(L) = QS(L)*RH_C(L)
        else
          Q(L) = 0.
        endif

!-----------------------------------------------------------------------
! Approximate PAR, assume atmosphere is neutral and setup vector fields.
!-----------------------------------------------------------------------
        PAR(L) = 0.45 * SW(L)
        Z0S(L) = Z0_SOIL
        V_SAT(L) = VSAT
        V_WILT(L) = VWILT
        RA(L) = 0.0

      enddo

!-----------------------------------------------------------------------
! Calculate the water and carbon fluxes for each FT.
!-----------------------------------------------------------------------
      do N=1,NPFT

!-----------------------------------------------------------------------
! Calculate the soil moisture stress factor.
!-----------------------------------------------------------------------
        V_ROOT(L) = M(L)*MAF(N)/(RHO_WATER*ROOTDEP)

        if (V_ROOT(L) .gt. V_CRIT) then
          FSMC(L) = 1.0
        elseif (V_ROOT(L) .le. VWILT) then
          FSMC(L) = 0.0
        else
          FSMC(L) = (V_ROOT(L) - VWILT)
     &              /(V_CRIT - VWILT)
        endif

!-----------------------------------------------------------------------
! Calculate available energy when the surface temperature is equal to
! the atmospheric temperature (AHAT).
!-----------------------------------------------------------------------
        call SWRAD (POINTS, LAND_PTS, LAND_INDEX, ALBSNF(1,N)
     &,             ALBSNC(1,N), LYING_SNOW, SW, TSTAR(1,N), TM
     &,             SWN_P(1,N), ALBLAND_VEG(1,N))

!-----------------------------------------------------------------------
! Calculate the canopy resistance and the resistance factor.
!-----------------------------------------------------------------------
        call SF_STOM (LAND_PTS, LAND_INDEX, N, CO2, FSMC
     &,               HT(1,N), PAR, LAI(1,N), PSTAR, Q, RA, SAT_D
     &,               ZERODEGC, EPCO2, EPSILON, GPP(1,N), NPP(1,N)
     &,               RESP_P(1,N), RESP_W(1,N), GC(1,N))

!CDIR NODEP
        do I=1,LAND_PTS
          L = LAND_INDEX(I)
          RS(L,N) = 1.0 / GC(L,N)
!-----------------------------------------------------------------------
! If snow assume no surface resistance
!-----------------------------------------------------------------------
          if (LYING_SNOW(L) .gt. 50.0) RS(L,N) = 0.0
        enddo

!-----------------------------------------------------------------------
! Solve the surface energy balance.
!-----------------------------------------------------------------------
        call PENMON (POINTS, LAND_PTS, LAND_INDEX, ROOTDEP
     &,              HCON_SOIL, RS(1,N), Z0(1,N), LW, SWN_P(1,N)
     &,              PSTAR, Q, SAT_D, TS1, WIND, Z1_REF, LHC, LHF
     &,              EPSILON, SIGMA, TM, ZERODEGC, ETRAN(1,N)
     &,              LETRAN(1,N), SH_T(1,N), G_T(1,N), RADNET_T(1,N)
     &,              TSTAR(1,N), MNEG, LW_OUT_P(1,N), LYING_SNOW)

!-----------------------------------------------------------------------
! Calculate the leaf turnover rate.
!-----------------------------------------------------------------------
        call LEAF_LIT (LAND_PTS, LAND_INDEX, N
     &,                FSMC, TSTAR(1,N), G_LEAF(1,N))

      enddo ! FT Loop

!-----------------------------------------------------------------------
! Calculate the bare soil evaporation.
!-----------------------------------------------------------------------
      call SWRAD (POINTS, LAND_PTS, LAND_INDEX ,ALBSOIL, ALBSNOW
     &,           LYING_SNOW, SW, TSOIL, TM, SWN_S, ALBLAND_SOIL)

!CDIR NODEP
!-----------------------------------------------------------------------
! Calculate the soil moisture stress factor.
!-----------------------------------------------------------------------
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        V_ROOT(L) = M(L)*MAF(SOIL)/(RHO_WATER*ROOTDEP)
        if (V_ROOT(L) .gt. V_CRIT) then
          FSMC(L) = 1.0
        elseif (V_ROOT(L) .le. VWILT) then
          FSMC(L) = 0.0
        else
          FSMC(L) = (V_ROOT(L) - VWILT)
     &            /(V_CRIT - VWILT)
        endif
        if (FSMC(L) .gt. (RSS/1.0E6)) then
          RSSOIL(L) = RSS/FSMC(L)
        else
          RSSOIL(L) = 1.0E6
        endif

!-----------------------------------------------------------------------
! If snow assume no surface resistance
!-----------------------------------------------------------------------
        if (LYING_SNOW(L) .gt. 50.0) RSSOIL(L) = 0.0
      enddo

!-----------------------------------------------------------------------
! Solve the soil surface energy balance.
!-----------------------------------------------------------------------
      call PENMON (POINTS, LAND_PTS, LAND_INDEX, ROOTDEP, HCON_SOIL
     &,            RSSOIL, Z0S, LW, SWN_S, PSTAR, Q, SAT_D, TS1, WIND
     &,            Z1_REF, LHC, LHF, EPSILON, SIGMA, TM, ZERODEGC
     &,            ESOIL, LESOIL, HSOIL, GSOIL, RADNETSOIL, TSOIL
     &,            MNEG, LW_OUT_S, LYING_SNOW)

!-----------------------------------------------------------------------
! Calculate gridbox mean fluxes and surface temperature.
!-----------------------------------------------------------------------

!!  !CDIR NODEP
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        ET(L) = 0.0
        LE(L) = 0.0
        SH(L) = 0.0
        RADNET(L) = 0.0
        G(L) = 0.0
        TSTAR_GB(L) = 0.0
        ALBLAND(L) = 0.0
        LW_OUT(L) = 0.0
        SWN(L) = 0.0
!!  !CDIR unroll=5
        do N=1,NPFT
          ET(L) = ET(L) + FRAC(L,N)*ETRAN(L,N)
          LE(L) = LE(L) + FRAC(L,N)*LETRAN(L,N)
          TSTAR_GB(L) = TSTAR_GB(L) + FRAC(L,N)*TSTAR(L,N)
          SH(L) = SH(L) + FRAC(L,N)*SH_T(L,N)
          G(L) = G(L) + FRAC(L,N)*G_T(L,N)
          RADNET(L) = RADNET(L) + FRAC(L,N)*RADNET_T(L,N)
          ALBLAND(L) = ALBLAND(L) + FRAC(L,N)*ALBLAND_VEG(L,N)
          LW_OUT(L) = LW_OUT(L) + FRAC(L,N)*LW_OUT_P(L,N)
          SWN(L) = SWN(L) + FRAC(L,N)*SWN_P(L,N)
        enddo
        ET(L) = ET(L) + (1 - VEG_FRAC(L))*ESOIL(L)
        ALBLAND(L) = ALBLAND(L) + (1 - VEG_FRAC(L))*ALBLAND_SOIL(L)
        LW_OUT(L) = LW_OUT(L) + (1 - VEG_FRAC(L))*LW_OUT_S(L)
        SWN(L) = SWN(L) + (1 - VEG_FRAC(L))*SWN_S(L)
!       If snow is lying assume that evapotranspiration is zero
        if (LYING_SNOW(L) .gt. 50.0) then
          ESUB(L) = ET(L)
          ET(L) = 0.
        else
          ESUB(L) = 0.
        endif
        LE(L) = LE(L) + (1 - VEG_FRAC(L)) * LESOIL(L)
        TSTAR_GB(L) = TSTAR_GB(L) + (1 - VEG_FRAC(L)) * TSOIL(L)
        SH(L) = SH(L) + (1 - VEG_FRAC(L)) * HSOIL(L)
        G(L) = G(L) + (1 - VEG_FRAC(L)) * GSOIL(L)
        RADNET(L) = RADNET(L) + (1 - VEG_FRAC(L)) * RADNETSOIL(L)
      enddo

!-----------------------------------------------------------------------
! Calculate the soil respiration
!-----------------------------------------------------------------------
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        STHU(L) = V_ROOT(L)/VSAT
      enddo
      call MICROBE (POINTS, LAND_PTS, LAND_INDEX, CS, STHU, V_SAT
     &,             V_WILT, TS1, RESP_S)

!-----------------------------------------------------------------------
! Update phenology accumulation variable.
!-----------------------------------------------------------------------

!!  !CDIR NODEP
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
!!  !CDIR unroll=5
        do N=1,NPFT
          G_LEAF_DAY(L,N) = G_LEAF_DAY(L,N) + G_LEAF(L,N)*FTIME_PHEN
        enddo
      enddo

!-----------------------------------------------------------------------
! Update leaf phenological state
!-----------------------------------------------------------------------
      if (L_PHEN) then
# if defined O_mtlm_segday
        if (DAY_PHEN .gt. 1) then
          DTIME_PHEN = real(DAY_PHEN)*segtim/DAY_YEAR
        else
          DTIME_PHEN = real(DAY_PHEN)/DAY_YEAR
        endif
# else
        DTIME_PHEN = real(DAY_PHEN)/DAY_YEAR
# endif
        do N=1,NPFT
          call PHENOL (LAND_PTS, LAND_INDEX, N, G_LEAF_DAY(1,N)
     &,                HT(1,N), DTIME_PHEN, G_LEAF_PHEN(1,N), LAI(1,N))
          do I=1,LAND_PTS
            L = LAND_INDEX(I)
            G_LEAF_DAY(L,N) = 0.0
          enddo
        enddo
      endif ! End of PHENOL call

!-----------------------------------------------------------------------
! Update the land surface state
!-----------------------------------------------------------------------
# if defined O_mtlm_segday
      call MTLM_STATE (POINTS, LAND_PTS, LAND_INDEX, ROOTDEP, HCAP_SOIL
     &,                SATCON, VSAT, LHF, TM, TIMESTEP*segtim, G, RAIN
     &,                SNOW, ET, ESUB, M, LYING_SNOW, TS1, SURF_ROFF
     &,                SNOWMELT, MNEG)
# else
      call MTLM_STATE (POINTS, LAND_PTS, LAND_INDEX, ROOTDEP, HCAP_SOIL
     &,                SATCON, VSAT, LHF, TM, TIMESTEP, G, RAIN, SNOW
     &,                ET, ESUB, M, LYING_SNOW, TS1, SURF_ROFF
     &,                SNOWMELT, MNEG)
# endif

!-----------------------------------------------------------------------
! Accumulate TRIFFID driving variables
!-----------------------------------------------------------------------
# if defined O_mtlm_carbon_13
      do N=1,NPFT
        B13NPP = AC13NPP(N)*RC13A
        B13NPP = MAX(0.5*RC13STD, MIN(2.*RC13STD,B13NPP))
        G13NPP(N) = B13NPP/(1.+B13NPP)
      enddo
# endif
# if defined O_mtlm_carbon_14
      do N=1,NPFT
        B14NPP = AC14NPP(N)*RC14A
        B14NPP = MAX(0.5*RC14STD, MIN(2.*RC14STD,B14NPP))
        G14NPP(N) = B14NPP/(1.+B14NPP)
      enddo
# endif
      do I=1,LAND_PTS
        L = LAND_INDEX(I)
        RESP_S_DR(L) = RESP_S_DR(L)+RESP_S(L)*FTIME*SEC_YEAR
# if defined O_mtlm_carbon_13
        B13RES = CS13(L)/(CS(L)-CS13(L))
        B13RES = MAX(0.5*RC13STD,MIN(2.*RC13STD,B13RES))
        G13RES = B13RES/(1.+B13RES)
        RESP_S13(L) = G13RES*RESP_S(L)
        RESP_S_DR13(L) = RESP_S_DR13(L)+RESP_S13(L)*FTIME*SEC_YEAR
# endif
# if defined O_mtlm_carbon_14
        B14RES = CS14(L)/(CS(L)-CS14(L))
        B14RES = MAX(0.5*RC14STD,MIN(2.*RC14STD,B14RES))
        G14RES = B14RES/(1.+B14RES)
        RESP_S14(L) = G14RES*RESP_S(L)
        RESP_S_DR14(L) = RESP_S_DR14(L)+RESP_S14(L)*FTIME*SEC_YEAR
# endif
        do N=1,NPFT
          NPP_DR(L,N) = NPP_DR(L,N) + NPP(L,N)*FTIME*SEC_YEAR
# if defined O_mtlm_carbon_13
          NPP13(L,N) = G13NPP(N)*NPP(L,N)
          NPP_DR13(L,N) = NPP_DR13(L,N) + NPP13(L,N)*FTIME*SEC_YEAR
# endif
# if defined O_mtlm_carbon_14
          NPP14(L,N) = G14NPP(N)*NPP(L,N)
          NPP_DR14(L,N) = NPP_DR14(L,N) + NPP14(L,N)*FTIME*SEC_YEAR
# endif
          G_LEAF_DR(L,N) = G_LEAF_DR(L,N) + G_LEAF_PHEN(L,N)*FTIME
          RESP_W_DR(L,N) = RESP_W_DR(L,N) + RESP_W(L,N)*FTIME*SEC_YEAR
        enddo
      enddo

!---------------------------------------------------------------------
! Form the mean for variables given to the atmospheric and ocean model
!---------------------------------------------------------------------
      do j=2,jmt-1
        do i=2,imt-1
          L = land_map(i,j)
          if (L .ne. 0) then
            sbc(i,j,iro) = sbc(i,j,iro) + SURF_ROFF(L)*TIMESTEP*0.1
            sbc(i,j,isca) = sbc(i,j,isca) + (1. - ALBLAND(L))*TIMESTEP
            sbc(i,j,ievap) = sbc(i,j,ievap)
     &                     + 0.1*(ET(L) + ESUB(L))*TIMESTEP
            sbc(i,j,ilwr) = sbc(i,j,ilwr) - 1000.*(LW_OUT(L) - SW_C(L)
     &                    + SWN(L))*TIMESTEP
            sbc(i,j,isens) = sbc(i,j,isens) + 1000.*SH(L)*TIMESTEP
# if defined O_carbon
            do N=1,NPFT
              sbc(i,j,inpp) = sbc(i,j,inpp)
     &                      + NPP(L,N)*FRAC(L,N)*TIMESTEP
            enddo
            sbc(i,j,isr) = sbc(i,j,isr) + RESP_S(L)*TIMESTEP
# endif
# if defined O_mtlm_carbon_13
            do N=1,NPFT
              sbc(i,j,inpp13) = sbc(i,j,inpp13)
     &                      + NPP13(L,N)*FRAC(L,N)*TIMESTEP
            enddo
            sbc(i,j,isr13) = sbc(i,j,isr13) + RESP_S13(L)*TIMESTEP
# endif
# if defined O_mtlm_carbon_14
            do N=1,NPFT
              sbc(i,j,inpp14) = sbc(i,j,inpp14)
     &                      + NPP14(L,N)*FRAC(L,N)*TIMESTEP
            enddo
            sbc(i,j,isr14) = sbc(i,j,isr14) + RESP_S14(L)*TIMESTEP
# endif
          endif
        enddo
      enddo
      atlnd = atlnd + TIMESTEP

!----------------------------------------------------------------------
! Update the vegetation areal coverages, structural parameters,
! and soil carbon.
!----------------------------------------------------------------------
      if (L_TRIF) then

        if (VEG_EQUIL) then
          FORW = 1.0
          GAMMA = GAMMA_EQ
          KITER = ITER_EQ
        else
          FORW = 0.0
# if defined O_mtlm_segday
          if (DAY_TRIF .gt. 1) then
            GAMMA = DAY_YEAR/(segtim*real(DAY_TRIF))
          else
            GAMMA = DAY_YEAR/real(DAY_TRIF)
          endif
# else
          GAMMA = DAY_YEAR/DAY_TRIF
# endif
          KITER = 1
        endif

!       diagnose the total land carbon before triffid
        do j=2,jmt-1
          do i=2,imt-1
            L = land_map(i,j)
            if (L .ne. 0) then
              sbc(i,j,iburn) = CS(L) - RESP_S_DR(L)/GAMMA
              do N=1,NPFT
                sbc(i,j,iburn) = sbc(i,j,iburn) + FRAC(L,N)*C_VEG(L,N)
     &                         + FRAC(L,N)*NPP_DR(L,N)/GAMMA
              enddo
            endif
          enddo
        enddo

        do I=1,KITER
          call TRIFFID (LAND_PTS, LAND_INDEX, FORW, GAMMA, FRAC_VS
     &,                 FRACA, FRAC_MIN, FRAC_SEED, DENOM_MIN, BF
     &,                 G_LEAF_DR, NPP_DR, RESP_S_DR, RESP_W_DR, CS
     &,                 FRAC, HT, LAI, C_VEG, CV, LIT_C, LIT_C_T
# if defined O_mtlm_carbon_13
     &,                 NPP_DR13, RESP_S_DR13
     &,                 RC13STD, CS13, C_VEG13, CV13
# endif
# if defined O_mtlm_carbon_14
     &,                 NPP_DR14, RESP_S_DR14
     &,                 RC14STD, CS14, C_VEG14, CV14
# endif
     &                  )
        enddo
        LAND_COUNTER = 0.

!----------------------------------------------------------------------
! Diagnose the amount of vegetation burnt (land use change emissions)
!----------------------------------------------------------------------
# if defined O_mtlm_segday
        RTSB = 1./(DAY_TRIF*segtim*SEC_DAY)
# else
        RTSB = 1./(DAY_TRIF*SEC_DAY)
# endif
        do j=2,jmt-1
          do i=2,imt-1
            L = land_map(i,j)
            if (L .ne. 0) then
              sbc(i,j,iburn) = sbc(i,j,iburn) - CS(L)
              do N=1,NPFT
                sbc(i,j,iburn) = sbc(i,j,iburn) - FRAC(L,N)*C_VEG(L,N)
              enddo
              sbc(i,j,iburn) =  sbc(i,j,iburn)*RTSB
            else
              sbc(i,j,iburn) = 0.
            endif
          enddo
        enddo
# if defined O_global_sums
        ntlbc = 0
# endif

!----------------------------------------------------------------------
! Zero the accumulated driving variables
!----------------------------------------------------------------------
        do I=1,LAND_PTS
          L = LAND_INDEX(I)
          RESP_S_DR(L) = 0.0
# if defined O_mtlm_carbon_13
          RESP_S_DR13(L) = 0.0
# endif
# if defined O_mtlm_carbon_14
          RESP_S_DR14(L) = 0.0
# endif
          do N=1,NPFT
            NPP_DR(L,N) = 0.0
# if defined O_mtlm_carbon_13
            NPP_DR13(L,N) = 0.0
# endif
# if defined O_mtlm_carbon_14
            NPP_DR14(L,N) = 0.0
# endif
            G_LEAF_DR(L,N) = 0.0
            RESP_W_DR(L,N) = 0.0
          enddo
        enddo

!----------------------------------------------------------------------
! Derive vegetation parameters from the new areal fractions and
! structural properties.
!----------------------------------------------------------------------
        do N=1,NPFT
          call PFT_SPARM (LAND_PTS, LAND_INDEX, N, ALBSOIL
     &,                   HT(1,N), LAI(1,N), ALBSNC(1,N), ALBSNF(1,N)
     &,                   CATCH(1,N), Z0(1,N))
        enddo

!----------------------------------------------------------------------
! Define other vegetation parameters
!----------------------------------------------------------------------
        do I=1,LAND_PTS
          L = LAND_INDEX(I)
          VEG_FRAC(L) = 0.0
          do N=1,NPFT
            VEG_FRAC(L) = VEG_FRAC(L) + FRAC(L,N)
          enddo
          FRAC_VS(L) = VEG_FRAC(L) + FRAC(L,SOIL)
        enddo

      endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!     accumulate time averages
!-----------------------------------------------------------------------
      call ta_mtlm_tavg (is, ie, js, je, 1)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     accumulate time step integrals
!-----------------------------------------------------------------------
      call ta_mtlm_tsi (is, ie, js, je, 1)
# endif
#endif

      return
      end
