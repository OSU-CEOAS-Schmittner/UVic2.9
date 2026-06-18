      subroutine setmtlm (is, ie, js, je)

#if defined O_mtlm
!-----------------------------------------------------------------------
!     Initialize the land surface and vegetation model
!-----------------------------------------------------------------------

      implicit none

      include "size.h"
      include "calendar.h"
      include "cembm.h"
      include "csbc.h"
      include "grdvar.h"
      include "coord.h"
      include "switch.h"
      include "levind.h"
      include "atm.h"
# if defined O_crop_data || defined O_pasture_data || defined O_agric_data
      include "veg.h"
# endif
# if defined O_ice && defined O_landice_data
#  if defined O_ice_cpts
      include "cpts.h"
#  endif
      include "ice.h"
# endif
      include "tmngr.h"
      include "mtlm.h"
# if defined O_mtlm_carbon_13
      include "mtlmc13.h"
# endif
# if defined O_mtlm_carbon_14
      include "mtlmc14.h"
# endif
      include "mtlm_data.h"

      character(120) :: fname, new_file_name

      integer i, ie, iem1, is, isp1, iou, j, je, jem1, js, jsp1, k, L, n

      logical exists, inqvardef

      real LAI_BAL, epsln
      parameter (epsln=1.0e-20)

      isp1 = is+1
      iem1 = ie-1
      jsp1 = js+1
      jem1 = je-1
      atlnd = 1.

!-----------------------------------------------------------------------
!     Default parameters
!-----------------------------------------------------------------------
# if defined O_embm
      LHC = vlocn*1.e-4
      LHF = flice*1.e-4
# else
      LHC = 2.501E6
      LHF = 0.334E6
# endif
      SIGMA = 5.67E-8
      DAY_YEAR = yrlen
      SEC_DAY = daylen
      SEC_YEAR = DAY_YEAR*SEC_DAY
      STEP_DAY = INT(SEC_DAY/TIMESTEP)
!     set longitude to half the STEP_DAY interval to line up with noon.
!     longitude is set to a constant to treat all latitudes equally.
      LONG(:) = 360./STEP_DAY/2.
      LAND_COUNTER = 0
      dtlnd = TIMESTEP

!----------------------------------------------------------------------
!     Initialization of arrays
!----------------------------------------------------------------------
      PSTAR(:) = 1.e5
      DTEMP_DAY(:) = 0.
      LYING_SNOW(:) = 0.
      TSTAR(:,:) = 280.
      TSOIL(:) = 280.0
      TS1(:) = 280.0
      CS(:) = 10.0
      M(:) = 242.
      MNEG(:) = 0.
      FSMC(:) = 1.0
      RESP_S_DR(:) = 0.0
      ALBSOIL(:) = 0.3
      ALBSNOW(:) = 0.6
      Z0S(:) = 0.0003
      FRACA(:) = 0.0
      FRAC(:,:) = 0.1
      LAI(:,1:2) = 6.
      LAI(:,3:5) = 2.
      HT(:,1:2) = 21.46
      HT(:,3:4) = 0.794
      HT(:,5) = 1.587
      NPP_DR(:,:) = 0.0
      G_LEAF_DR(:,:) = 0.0
      RESP_W_DR(:,:) = 0.0

!-----------------------------------------------------------------------
!     Define externally dependent arrays
!-----------------------------------------------------------------------
      L = 0
      do j=jsp1,jem1
        do i=isp1,iem1
          if (kmt(i,j) .le. klmax) then
            L = L + 1
            GAREA(L) = dxt(i)*dyt(j)*cst(j)*1e-4
# if defined O_mtlm_pressure
            PSTAR(L) = PSTAR(L)*exp((-1.)*elev(i,j)/7.3e5)
# endif
# if defined O_crop_data || defined O_pasture_data || defined O_agric_data
            FRACA(L) = agric(i,j,2)
# endif
            LAT(L) = tlat(i,j)
            sbc(i,j,isca) = 1. - ALBSOIL(L)
            sbc(i,j,ievap) = 0.
            sbc(i,j,isens) = 0.
            sbc(i,j,ilwr) = 0.
          endif
        enddo
      enddo

      if (L .gt. POINTS) then
        print*, "==> Error: Number of land points is inconsistent"
        print*, "==>        set POINTS in size.h to: ", L
        stop
      endif

!----------------------------------------------------------------------
!     Initialize the non-vegetation fractions
!----------------------------------------------------------------------
      FRAC(:,SOIL) = 1.0
      do N=1,NPFT
        FRAC(:,SOIL) = FRAC(:,SOIL) - FRAC(:,N)
      enddo

!----------------------------------------------------------------------
!     Initialize the vegetation carbon contents
!----------------------------------------------------------------------
      CV(:) = 0.
      G_LEAF_PHEN(:,:) = 0.0
      do N=1,NPFT
        LAI_BAL = (A_WS(N)*ETA_SL(N)*HT(1,N)/A_WL(N))
     &                 **(1.0/(B_WL(N)-1))
        C_VEG(:,N) = 2*SIGL(N)*LAI_BAL
     &             + A_WS(N)*ETA_SL(N)*HT(:,N)*LAI_BAL
        CV(:) = CV(:) + C_VEG(:,N)*FRAC(:,N)
      enddo

!----------------------------------------------------------------------
!     Derive vegetation parameters from the areal fractions and the
!     structural properties.
!----------------------------------------------------------------------
      L = 0
      LAND_INDEX(:) = 0
      do j=jsp1,jem1
        do i=isp1,iem1
          if (kmt(i,j) .le. klmax) then
            L = L + 1
            LAND_INDEX(L) = L
          endif
        enddo
      enddo
      do N=1,NPFT
        call PFT_SPARM (L, LAND_INDEX, N, ALBSOIL, HT(1,N), LAI(1,N)
     &,                 ALBSNC(1,N), ALBSNF(1,N), CATCH(1,N), Z0(1,N))
      enddo

!----------------------------------------------------------------------
!     Define other vegetation parameters
!----------------------------------------------------------------------
      VEG_FRAC(:) = 0.0
      do N=1,NPFT
        VEG_FRAC(:) = VEG_FRAC(:) + FRAC(:,N)
      enddo
      FRAC_VS(:) = VEG_FRAC(:) + FRAC(:,SOIL)

!----------------------------------------------------------------------
!       Reading a restart file
!----------------------------------------------------------------------
      if (.not. init) then
        fname = new_file_name ('restart_mtlm.nc')
        inquire (file=trim(fname), exist=exists)
        if (exists) call mtlm_rest_in (fname, is, ie, js, je)
# if defined O_restart_2
        fname = new_file_name ("restart_2_mtlm.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) call mtlm_rest_in (fname, is, ie, js, je)
# endif
      endif
# if defined O_mtlm_carbon_13
      do N=1,NPFT
         if (C3(N).eq.1) then
            ac13npp(N) = 0.979
         else
            ac13npp(N) = 0.993
         endif
c         ac13npp(N) = 1.
      enddo
      print*, "ac13npp:", ac13npp(:)

      if (BF.ne.0.) then
        print*, "Error in setmtlm:"
        print*, "if O_mtlm_carbon_13 option is used BF must be zero"
        stop
      endif
# endif
# if defined O_mtlm_carbon_14
      do N=1,NPFT
         if (C3(N).eq.1) then
            ac14npp(N) = 0.958
         else
            ac14npp(N) = 0.986
         endif
c         ac14npp(N) = 1.
      enddo
      print*, "ac14npp:", ac14npp(:)

      if (BF.ne.0.) then
        print*, "Error in setmtlm:"
        print*, "if O_mtlm_carbon_14 option is used BF must be zero"
        stop
      endif
# endif

!----------------------------------------------------------------------
!     Create the VEG_INDEX array of land points with each type
!----------------------------------------------------------------------
      L = 0
      land_map(:,:) = 0
      LAND_PTS = 0
      LAND_INDEX(:) = 0
      do j=jsp1,jem1
        do i=isp1,iem1
          if (kmt(i,j) .le. klmax) then
            L = L + 1
# if defined O_ice && defined O_landice_data
            if (aicel(i,j,2) .lt. 0.5 .and. tmsk(i,j) .lt. 0.5) then
              land_map(i,j) = L
              LAND_PTS = LAND_PTS + 1
              LAND_INDEX(LAND_PTS) = L
            endif
# else
            if (tmsk(i,j) .lt. 0.5) then
              land_map(i,j) = L
              LAND_PTS = LAND_PTS + 1
              LAND_INDEX(LAND_PTS) = L
            endif
# endif
          endif
        enddo
      enddo
      do N=1,NPFT
        VEG_PTS(N) = 0
        do J=1,LAND_PTS
          L = LAND_INDEX(J)
          if (FRAC(L,N) .gt. FRAC_MIN + epsln) then
            VEG_PTS(N) = VEG_PTS(N) + 1
            VEG_INDEX(VEG_PTS(N),N) = L
          endif
        enddo
      enddo

# if defined O_mtlm_segday
      if (DAY_TRIF .lt. segtim) then
        print*, ""
        print*, "==> Warning: DAY_TRIF is set to be less than segtim."
        print*, "             with option mtlm_segday, triffid will"
        print*, "             only be done once every coupling time."
      endif
      if (DAY_PHEN .lt. segtim) then
        print*, ""
        print*, "==> Warning: DAY_PHEN is set to be less than segtim."
        print*, "             with option mtlm_segday, phenology will"
        print*, "             only be done once every coupling time."
      endif
      if (segtim .lt. 1.) then
        print*, ""
        print*, "==> Error: segtim must be greater than one when using"
        print*, "           the option mtlm_segday. Turn off this"
        print*, "           option if segtim is less than one."
        stop
      endif
# endif
      if (STEP_DAY .gt. STEPSM) then
        print*, ""
        print*, "==> Error: STEPSM is too small. Increase TIMESTEP or "
        print*, "           set STEPSM in size.h to: ", STEP_DAY
        stop
      endif
# if defined O_mtlm_segday
      if (mod(SEC_DAY*segtim,TIMESTEP) .gt. 1.e-6) then
        print*, ""
        print*, "==> Error: there must be an integral number of mtlm "
        print*, "           timesteps in a coupling time."
        stop
      endif
      if (DAY_TRIF .gt. 1) then
        if (mod(FLOAT(DAY_TRIF),segtim) .gt. 1.e-6) then
          print*, '==> Error: there must be an integral number coupling'
     &,     ' timesteps within DAY_TRIF when using O_mtlm_segday.'
          stop
        else
          DAY_TRIF = int(float(DAY_TRIF)/segtim)
        endif
      endif
      if (DAY_PHEN .gt. 1) then
        if (mod(FLOAT(DAY_PHEN),segtim) .gt. 1.e-6) then
          print*, '==> Error: there must be an integral number coupling'
     &,     ' timesteps within DAY_PHEN when using O_mtlm_segday.'
          stop
        else
          DAY_PHEN = int(float(DAY_PHEN)/segtim)
        endif
      endif
# else
      if (mod(SEC_DAY,TIMESTEP) .gt. 1.e-6) then
        print*, ""
        print*, "==> Error: there must be an integral number of mtlm "
        print*, "           timesteps in a day."
        stop
      endif
# endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!     zero time average accumulators
!-----------------------------------------------------------------------

      call ta_mtlm_tavg (is, ie, js, je, 0)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_mtlm_tsi (is, ie, js, je, 0)
# endif
#endif

      return
      end
