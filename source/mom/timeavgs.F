      subroutine avgset (xt, xu, yt, yu, zt, zw, imkmax
     &,          cvxz, cvx, cvy, cvz, javgr, imav, jmav, levav)

#if defined O_mom
!-----------------------------------------------------------------------
!     setup the "averaging" grid for use with option "time_averages" in
!     MOM for accumulating and saving time averaged data.

!     Warning: whenever the definition of this time  averaging
!              grid is changed below, it must be re-installed by
!              running the "timeavgs.F" module.

!     input:
!       imt   = number of longitudes in the MOM grid
!       jmt   = number of latitudes in the MOM grid
!       km    = number of levels in the MOM grid
!       xt    = longitudinal coordinates of "t" points (deg)
!       xu    = longitudinal coordinates of "u" points (deg)
!       yt    = latitudinal coordinates of "t" points (deg)
!       yu    = latitudinal coordinates of "u" points (deg)
!       zt    = depth of "t" points (cm)
!       zw    = depth of "w" points (bottom of "t" cells) (cm)
!       imkmax= imt*km for dimensioning purposes

!     output:
!       cvxz  = indicates points in (i,k) plane for time averaging data
!       cvx   = indicates points in longitude for time averaging data
!       cvy   = indicates points in latitude for time averaging data
!       cvz   = indicates points in depth for time averaging data
!       javgr = number of the latitude row containing averaging points
!       imav  = number of longitudinal points for averaging data
!       jmav  = number of latitudinal points for averaging data
!       levav = number of depth points for averaging data
!-----------------------------------------------------------------------

      implicit none

      integer imkmax, jrow, i, k, j, imav, jmav, levav, lcvxz
      integer lcvx, lcvy, lcvz

      real cksum, checksum

      include "size.h"
      include "stdunits.h"

      integer cvxz(imkmax), cvx(imt), cvy(jmt), cvz(km), javgr(jmt)

      logical savex(imt), savey(jmt), savez(km)

      real xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km)
      real xtav(imt), xuav(imt), ytav(jmt), yuav(jmt)
      real ztav(km), zwav(km)

!     initialize all points as not being on the "averaging" grid

!     savex(i) = (t,f) implies that the model grid point longitude
!                corresponding to grid point index "i" (is, is not)
!                on the "averaging" grid
!     savey(j) = (t,f) implies that the model grid point latitude
!                corresponding to grid point index "j" (is, is not)
!                on the "averaging" grid
!     savez(k) = (t,f) implies that the model grid point depth
!                correspondingto grid point index "k" (is, is not)
!                on the "averaging" grid

!     a particular grid point with index (i,j,k) is on the "averaging"
!     grid when "savex(i)", "savey(j)", & "savez(k)" are all true.

      do jrow=1,jmt
        savey(jrow) = .true.
      enddo

      do i=1,imt
        savex(i) = .true.
      enddo

      do k=1,km
        savez(k) = .true.
      enddo
      savey(1)   = .false.
      savey(jmt) = .false.

!     save indices for jrows which define the "averaging" grid

      j = 0
      do jrow=1,jmt
        javgr(jrow) = 0
        if (savey(jrow)) then
          j     = j + 1
          javgr(jrow) = j
        endif
      enddo

!----------------------------------------------------------------------
!     calculate the size of the "averaging" grid
!----------------------------------------------------------------------

      imav = 0
      do i=1,imt
        if (savex(i)) imav = imav+1
      enddo

      jmav = 0
      do jrow=1,jmt
        if (savey(jrow)) jmav = jmav + 1
      enddo

      levav = 0
      do k=1,km
        if (savez(k)) levav = levav + 1
      enddo

!----------------------------------------------------------------------
!     establish control vectors for gathering data from the MOM grid
!     onto the "averaging" grid
!----------------------------------------------------------------------

      lcvxz = 0
      do k=1,km
        if (savez(k)) then
          do i=1,imt
            if (savex(i)) then
              lcvxz = lcvxz + 1
              cvxz(lcvxz) = (k-1)*imt + i
            endif
          enddo
        endif
      enddo

      lcvx = 0
      do i=1,imt
        if (savex(i)) then
          lcvx = lcvx + 1
          cvx(lcvx) = i
        endif
      enddo

      lcvy = 0
      do jrow=1,jmt
        if (savey(jrow)) then
          lcvy = lcvy + 1
          cvy(lcvy) = jrow
        endif
      enddo

      lcvz = 0
      do k=1,km
        if (savez(k)) then
          lcvz = lcvz + 1
          cvz(lcvz) = k
        endif
      enddo

      do i=1,imav
        xtav(i) = xt(cvx(i))
        xuav(i) = xu(cvx(i))
      enddo
      do j=1,jmav
        ytav(j) = yt(cvy(j))
        yuav(j) = yu(cvy(j))
      enddo
      do k=1,levav
        ztav(k) = zt(cvz(k))
        zwav(k) = zw(cvz(k))
      enddo

!----------------------------------------------------------------------
!     show the "averaging" grid coordinates
!----------------------------------------------------------------------

      write (stdout,'(//,20x,a)')
     & 'T I M E   A V E R A G E S   G R I D   I N I T I A L I Z E D'
      write (stdout,'(//1x,a,/)')
     & ' Time averages will be taken at the following grid points:'
      write (stdout,9900) imav,' grid "t" longitudes:'
     &,                   (xtav(i),i=1,imav)
      write (stdout,9900) imav,' grid "v" longitudes:'
     &,                   (xuav(i),i=1,imav)
      write (stdout,9900) jmav,' grid "t" latitudes:'
     &,                   (ytav(j),j=1,jmav)
      write (stdout,9900) jmav,' grid "v" latitudes:'
     &,                   (yuav(j),j=1,jmav)
      write (stdout,9900) levav,' "t" grid depths:'
     &,                   (ztav(k),k=1,levav)
      write (stdout,9900) levav,' "w" grid depths:'
     &,                   (zwav(k),k=1,levav)

!---------------------------------------------------------------------
!     compute a grid checksum
!---------------------------------------------------------------------

      cksum = 0.0
      cksum = cksum + checksum (xtav, imav, 1)
      cksum = cksum + checksum (ytav, jmav, 1)
      cksum = cksum + checksum (ztav, levav, 1)
      cksum = cksum + checksum (xuav, imav, 1)
      cksum = cksum + checksum (yuav, jmav, 1)
      cksum = cksum + checksum (zwav, levav, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Time average grid checksum = ',cksum
      write (stdout,'(/)')
      return
9900  format (1x,i4,a,/ (1x,10f10.2))
      end

      subroutine avgvar (j, jrow, w, varu, vart, stf, smf, mapt)

!-----------------------------------------------------------------------

!     "avgvar" produces a four-dimensional dataset of time averaged
!      quantities on a grid defined by the "timeavgs.F" module.
!      this grid may be the entire MOM grid, a region of it, or any
!      coarser subset of it.  For example:
!      if MOM uses a 2 deg x 2 deg grid, the "averaging" grid
!      could be a 6 deg x 8 deg subset of part or the entire domain.
!      spatial resolution of the "averaging" grid is controlled in
!      the USER INPUT section of "avgset".

!      how to use:

!      1) set up the "averaging" grid in "avgset" (in this module)
!      2) run this "timeavgs.F" module and follow its directions
!      3) enable the "time_averages" option
!-----------------------------------------------------------------------
# if defined O_time_averages

      implicit none

      character(120) :: fname
      character(32) :: nstamp

      integer ntrec, nyear, nmonth, nday, nhour, nmin, nsec
      integer n3dvar, n3dtra, n3dvam, n2dvar, imkmav, imjmav, lenrow
      integer javgr, jrow, jav, n, i, j, nvel, navgts, jj, l, k, is
      integer ie, navgp, imav, jmav, levav, cvxz, cvx, cvy, cvz

      real spbuf, spbuf2, rnavgt, rnta_conv, rnta_sscar, rnta_botcar
      real avg3d, avg2d, avgper, time, ztav, xuav, ytav, zwav, xtav
      real yuav, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "cregin.h"
      include "timeavgs.h"
      include "npzd.h"

!     the following parameters are needed to define the workspace area

!     n3dvar  = #  of 3d fields (#  of prognostic variables + one for w)
!     n2dvar  = #  of 2d fields ("stf(nt)","smf(2)", and "psi")

      parameter (n3dvar=nt+2+1, n3dtra=nt)
      parameter (n3dvam=n3dvar-1, n2dvar=nt+2+1)
      parameter (imkmav=imtav*kmav, imjmav=imtav*jmtav
     &,          lenrow=imkmav*n3dvar+imtav*n2dvar)

      common /avgblk/ xtav(imtav), xuav(imtav), ytav(jmtav+2)
      common /avgblk/ yuav(jmtav+2), ztav(kmav), zwav(kmav)
      common /avgblk/ avg3d(imkmav,n3dvar), avg2d(imtav,n2dvar)
      common /avgblk/ spbuf(imkmav,n3dvar,jmtav)
      common /avgblk/ spbuf2(imtav,n2dvar,jmtav)
      character(12) :: name3d, name2d
      common /avgblc/ name3d(n3dvar), name2d(n2dvar+1)
      character(10) :: mapt(nt), mapta(nt)
      common /avgblc/ mapta

!     imav    =  number of points in longitude on the "averaging" grid
!     jmav    =  number of points in latitude on the "averaging" grid
!     levav   =  number of points in depth on the "averaging" grid

!     javgr   =  "averaging" grid "jrow" number. it maps "jrow" on the
!                MOM grid into a "jrow" on the "averaging" grid.

!     cvxz    =  control vector used in gathering data from array
!                dimensioned (imt,km) into one dimensioned (imav,levav)
!     cvx     =  control vector used in gathering data from array
!                dimensioned (imt) into one dimensioned (imav)
!     cvy     =  control vector used in gathering data from array
!                dimensioned (jmt) into one dimensioned (jmav)
!     cvz     =  control vector used in gathering data from array
!                dimensioned (km) into one dimensioned (levav)

!     navgts  =  a counter for tracking the number of time steps
!                 within an averaging period.
!     navgp   =  a counter for tracking the number of averaging
!                 periods.

      common /avgbli/ cvxz(imtkm), cvx(imt), cvy(jmt), cvz(km)
      common /avgbli/ javgr(jmt), imav, jmav, levav, navgts, navgp

      include "coord.h"
      include "docnam.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "calendar.h"
#  if defined O_embm
      include "cembm.h"
      include "atm.h"
#  endif

      real varu(imtkm,jmw,2), w(imtkm)
      real vart(imtkm,jmw,n3dtra)
      real stf(imt,1:jmw,nt), smf(imt,1:jmw,2)
      real tmask(imt,km), umask(imt,km)

      mapta(:) =  mapt(:)

!-----------------------------------------------------------------------
!     only process those "jrows" that are on the "averaging" grid
!-----------------------------------------------------------------------

      if (javgr(jrow) .ne. 0) then

        jav = javgr(jrow)

!-----------------------------------------------------------------------
!       integrate three dimensional data only on those longitudes
!       that are included in the "averaging" grid
!-----------------------------------------------------------------------

        do n=1,n3dtra
          do i=1,imkmav
            spbuf(i,n,jav) = spbuf(i,n,jav) + vart(cvxz(i),j,n)
          enddo
        enddo
        do n=1,2
          do i=1,imkmav
            nvel = n + n3dtra
            spbuf(i,nvel,jav) = spbuf(i,nvel,jav) + varu(cvxz(i),j,n)
          enddo
        enddo
        n = n3dvar
        do i=1,imkmav
          spbuf(i,n,jav) = spbuf(i,n,jav) + w(cvxz(i))
        enddo

!-----------------------------------------------------------------------
!       integrate two dimensional fields only on those longitudes
!       that are included in the "averaging" grid
!-----------------------------------------------------------------------

        do n=1,n2dvar
          if (n .le. nt) then
#  if defined O_embm && !defined O_save_fluxes_with_virtual
            if (trim(mapt(n)).ne.'temp' .and. trim(mapt(n)).ne.'salt')
     &        then
              do i=1,imtav
                spbuf2(i,n,jav) = spbuf2(i,n,jav) + stf(cvx(i),j,n)
     &                          - vflux(i,jrow)*gaost(n)
              enddo
            else
              do i=1,imtav
                spbuf2(i,n,jav) = spbuf2(i,n,jav) + stf(cvx(i),j,n)
              enddo
            endif
#  else
            do i=1,imtav
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + stf(cvx(i),j,n)
            enddo
#  endif
          elseif (n .lt. n2dvar) then
            do i=1,imtav
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + smf(cvx(i),j,n-nt)
            enddo
          else
            do i=1,imtav
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + ps(cvx(i),jrow,1)
#  endif
#  if defined O_stream_function
              spbuf2(i,n,jav) = spbuf2(i,n,jav) + psi(cvx(i),jrow,1)
#  endif
            enddo
          endif
        enddo

!-----------------------------------------------------------------------
!       update integration counter once per time step on the last row
!-----------------------------------------------------------------------

        if (jav .eq. jmtav) navgts = navgts + 1

      endif

      return

      entry avgout

!=======================================================================

!     output all averaged quantities

!=======================================================================

!-----------------------------------------------------------------------
!     construct time mean quantities on the "averaging" grid
!-----------------------------------------------------------------------

      rnavgt = c1/navgts
#  if defined O_save_convection
      rnta_conv = 0.
      if (nta_conv .gt. 0) rnta_conv = c1/nta_conv
#  endif
#  if defined O_save_carbon_carbonate_chem
      rnta_sscar = 0.
      if (nta_sscar .gt. 0) rnta_sscar = c1/nta_sscar
#  endif
      do jj=1,jmtav

!       construct the time average

        do n=1,n3dvar
          do l=1,imkmav
            avg3d(l,n) = rnavgt*spbuf(l,n,jj)
          enddo
        enddo

        do n=1,n2dvar
          do l=1,imtav
            avg2d(l,n) = rnavgt*spbuf2(l,n,jj)
          enddo
        enddo
#  if defined O_gent_mcwilliams

        do k=1,km
          do i=2,imtav
            ta_vetiso(i,k,jj+1) = rnavgt*ta_vetiso(i,k,jj+1)
            ta_vntiso(i,k,jj+1) = rnavgt*ta_vntiso(i,k,jj+1)
            ta_vbtiso(i,k,jj+1) = rnavgt*ta_vbtiso(i,k,jj+1)
          enddo
        enddo
#  endif
#  if defined O_carbon && defined O_carbon_14

        do k=1,km
          do i=2,imtav
            ta_dc14(i,k,jj+1) = rnavgt*ta_dc14(i,k,jj+1)
          enddo
        enddo
#  endif
#  if defined O_save_kv

        do k=1,km
          do i=2,imtav
            ta_diff_cbt(i,k,jj+1) = rnavgt*ta_diff_cbt(i,k,jj+1)
          enddo
        enddo
#  endif
#  if defined O_save_npzd

        do k=1,kpzd
          do i=2,imtav
            ta_rnpp(i,k,jj+1) = rnavgt*ta_rnpp(i,k,jj+1)
            ta_rgraz(i,k,jj+1) = rnavgt*ta_rgraz(i,k,jj+1)
            ta_rmorp(i,k,jj+1) = rnavgt*ta_rmorp(i,k,jj+1)
            ta_rmorpt(i,k,jj+1) = rnavgt*ta_rmorpt(i,k,jj+1)
            ta_rmorz(i,k,jj+1) = rnavgt*ta_rmorz(i,k,jj+1)
            ta_rexcr(i,k,jj+1) = rnavgt*ta_rexcr(i,k,jj+1)
#   if defined O_npzd_nitrogen
            ta_rnpp_D(i,k,jj+1) = rnavgt*ta_rnpp_D(i,k,jj+1)
            ta_rgraz_D(i,k,jj+1) = rnavgt*ta_rgraz_D(i,k,jj+1)
            ta_rmorp_D(i,k,jj+1) = rnavgt*ta_rmorp_D(i,k,jj+1)
            ta_rnfix(i,k,jj+1) = rnavgt*ta_rnfix(i,k,jj+1)
#   endif
          enddo
        enddo

        do k=1,km
          do i=2,imtav
            ta_rremi(i,k,jj+1) = rnavgt*ta_rremi(i,k,jj+1)
            ta_rexpo(i,k,jj+1) = rnavgt*ta_rexpo(i,k,jj+1)
            ta_rexpocal(i,k,jj+1) = rnavgt*ta_rexpocal(i,k,jj+1)
#   if defined O_npzd_nitrogen
            ta_rdeni(i,k,jj+1) = rnavgt*ta_rdeni(i,k,jj+1)
#   endif
          enddo
        enddo
        do i=2,imtav
          ta_rprocal(i,jj+1) = rnavgt*ta_rprocal(i,jj+1)
        enddo
#  endif
#  if defined O_save_convection

        do i=2,imtav
          ta_totalk(i,jj+1) = rnta_conv*ta_totalk(i,jj+1)
          ta_vdepth(i,jj+1) = rnta_conv*ta_vdepth(i,jj+1)
          ta_pe(i,jj+1) = rnta_conv*ta_pe(i,jj+1)
        enddo
#  endif
#  if defined O_save_carbon_carbonate_chem

        do i=2,imtav
          ta_sspH(i,jj+1) = rnta_sscar*ta_sspH(i,jj+1)
          ta_ssCO3(i,jj+1) = rnta_sscar*ta_ssCO3(i,jj+1)
          ta_ssOc(i,jj+1) = rnta_sscar*ta_ssOc(i,jj+1)
          ta_ssOa(i,jj+1) = rnta_sscar*ta_ssOa(i,jj+1)
          ta_sspCO2(i,jj+1) = rnta_sscar*ta_sspCO2(i,jj+1)
        enddo
#  endif

        tmask(:,:) = c0
        umask(:,:) = c0
        do i=1,imtav
          if (kmt(i,jj+1) .gt. 0) tmask(i,1:kmt(i,jj+1)) = c1
          if (kmu(i,jj+1) .gt. 0) umask(i,1:kmt(i,jj+1)) = c1
        enddo

        if (jj .eq. 1) ntrec = 0
        is = 1
        ie = imt
        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        if (jj .eq. 1) call def_tavg
        call def_tavg_mom (fname)

#  if !defined O_save_time_relyear0
!       make output time relative to year 1
        time = time - 1.
#  endif
        avgper = timavgper*accel
        if (avgper .le. 1e-6) avgper = 0.
#  if defined O_save_time_endper
        tmp = 0.
#  elif defined O_save_time_startper
        tmp = 1.
#  else
        tmp = 0.5
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
        time = time - tmp*avgper/360.
#   elif defined O_calendar_gregorian
        time = time - tmp*avgper/365.25
#   else
        time = time - tmp*avgper/365.
#   endif
#  else
#   if defined O_calendar_360_day
        time = time*360. - tmp*avgper
#   elif defined O_calendar_gregorian
        time = time*365.25 - tmp*avgper
#   else
        time = time*365. - tmp*avgper
#   endif
#  endif

        call mom_tavg_out (fname, is, ie, jj+1, jj+1, imt, jmt, km, nt
     &,   xt, yt, zt, xu, yu, zw, dxt, dyt, dzt, dxu, dyu, dzw, avgper
     &,   time, nstamp, mapta, avg3d(1,1), avg3d(1,nt+1), avg3d(1,nt+2)
     &,   avg3d(1,nt+3), avg2d(is,1), avg2d(is,nt+1), avg2d(is,nt+2)
     &,   redctn
#  if defined O_gent_mcwilliams
     &,   ta_vetiso(is,1,jj+1), ta_vntiso(is,1,jj+1)
     &,   ta_vbtiso(is,1,jj+1)
#  endif
#  if defined O_carbon && defined O_carbon_14
     &,   ta_dc14(is,1,jj+1)
#  endif
#  if defined O_save_kv
     &,   ta_diff_cbt(is,1,jj+1)
#  endif
#  if defined O_save_npzd
     &,   kpzd, ta_rnpp(is,1,jj+1), ta_rgraz(is,1,jj+1)
     &,   ta_rmorp(is,1,jj+1), ta_rmorpt(is,1,jj+1)
     &,   ta_rmorz(is,1,jj+1), ta_rremi(is,1,jj+1)
     &,   ta_rexcr(is,1,jj+1), ta_rexpo(is,1,jj+1)
     &,   ta_rexpocal(is,1,jj+1), ta_rprocal(is,jj+1)
#   if defined O_npzd_nitrogen
     &,   ta_rnpp_D(is,1,jj+1), ta_rgraz_D(is,1,jj+1)
     &,   ta_rmorp_D(is,1,jj+1), ta_rnfix(is,1,jj+1)
     &,   ta_rdeni(is,1,jj+1)
#   endif
#  endif
#  if defined O_save_convection
     &,   ta_totalk(is,jj+1), ta_vdepth(is,jj+1), ta_pe(is,jj+1)
#  endif
#  if defined O_save_carbon_carbonate_chem
     &,   ta_sspH(is,jj+1), ta_ssCO3(is,jj+1), ta_ssOc(is,jj+1)
     &,   ta_ssOa(is,jj+1), ta_sspCO2(is,jj+1)
#  endif
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
     &,   avg2d(is,nt+3)
#  endif
#  if defined O_stream_function
     &,   avg2d(is,nt+3)
#  endif
     &,   kmt(is,jj+1), mskhr(is,jj+1), tmask(is,1), umask(is,1)
     &,   tlat(is,jj+1), tlon(is,jj+1), ulat(is,jj+1), ulon(is,jj+1)
     &,   tgarea(is,jj+1), ugarea(is,jj+1), ntrec)

        if (jj .eq. 1)
     &    write (stdout,'(a,i4,a,a,a,i10,a,a)') '=> Ocn time means # '
     &,     ntrec, ' written to ',trim(fname),' on ts = ',itt, ', '
     &,     stamp

!-----------------------------------------------------------------------
!       zero out the "averaging" data the for the next averaging period
!-----------------------------------------------------------------------

        do n=1,n3dvar
          do l=1,imkmav
            spbuf(l,n,jj) = c0
          enddo
        enddo

        do n=1,n2dvar
          do l=1,imtav
            spbuf2(l,n,jj) = c0
          enddo
        enddo

      enddo
#  if defined O_gent_mcwilliams

      ta_vetiso(:,:,:) = c0
      ta_vntiso(:,:,:) = c0
      ta_vbtiso(:,:,:) = c0
#  endif
#  if defined O_carbon && defined O_carbon_14
      ta_dc14(:,:,:) = c0
#  endif
#  if defined O_save_kv
      ta_diff_cbt(:,:,:) = c0
#  endif
#  if defined O_save_npzd
      ta_rnpp(:,:,:) = c0
      ta_rgraz(:,:,:) = c0
      ta_rmorp(:,:,:) = c0
      ta_rmorpt(:,:,:) = c0
      ta_rmorz(:,:,:) = c0
      ta_rexcr(:,:,:) = c0
      ta_rremi(:,:,:) = c0
      ta_rexpo(:,:,:) = c0
      ta_rexpocal(:,:,:) = c0
      ta_rprocal(:,:) = c0
#   if defined O_npzd_nitrogen
      ta_rnpp_D(:,:,:) = c0
      ta_rgraz_D(:,:,:) = c0
      ta_rmorp_D(:,:,:) = c0
      ta_rdeni(:,:,:) = c0
      ta_rnfix(:,:,:) = c0
#   endif
#  endif
#  if defined O_save_convection
      ta_totalk(:,:) = c0
      ta_vdepth(:,:) = c0
      ta_pe(:,:) = c0
#  endif
#  if defined O_save_carbon_carbonate_chem
      ta_sspH(:,:) = c0
      ta_ssCO3(:,:) = c0
      ta_ssOc(:,:) = c0
      ta_ssOa(:,:) = c0
      ta_sspCO2(:,:) = c0
#  endif

!-----------------------------------------------------------------------
!     zero out the "averaging" counter for the next averaging period
!-----------------------------------------------------------------------

      navgts = 0
#  if defined O_save_convection
      nta_conv = 0
#  endif
#  if defined O_save_carbon_carbonate_chem
      nta_sscar = 0
#  endif

      return

      entry avgi

!-----------------------------------------------------------------------
!     initialize counters for tracking number of time steps within
!     an averaging period and the number of averaging periods.
!-----------------------------------------------------------------------

      navgts = 0
      navgp  = 0
#  if defined O_save_convection
      nta_conv = 0
#  endif
#  if defined O_save_carbon_carbonate_chem
      nta_sscar = 0
#  endif

!-----------------------------------------------------------------------
!     setup the "averaging" grid (data will be time averaged only on
!     these grid cells)
!-----------------------------------------------------------------------

      call avgset (xt, xu, yt, yu, zt, zw, imtkm
     &,         cvxz, cvx, cvy, cvz, javgr, imav, jmav, levav)

!-----------------------------------------------------------------------
!     verify that the "averaging" grid size matches parameters set
!     in "timeavgs.h"
!-----------------------------------------------------------------------

      if (imav .ne. imtav .or. jmav .ne. jmtav .or. levav .ne. kmav)then
        write (stdout,*) '=> Error: number of grid points in averaging'
     &,' grid does not match the parameter setting in timeavgs.h'
        write (stdout,*) ' imtav=',imtav, ': avgset returns', imav
        write (stdout,*) ' jmtav=',jmtav, ': avgset returns', jmav
        write (stdout,*) '  kmav=',kmav,  ': avgset returns', levav
        stop '=>avgi'
      endif

!-----------------------------------------------------------------------
!     initialize time mean averages to "zero"
!-----------------------------------------------------------------------

      do jj=1,jmtav
        do n=1,n3dvar
          do l=1,imkmav
            spbuf(l,n,jj) = c0
          enddo
        enddo

        do n=1,n2dvar
          do l=1,imtav
            spbuf2(l,n,jj) = c0
          enddo
        enddo
      enddo
#  if defined O_gent_mcwilliams
      ta_vetiso(:,:,:) = c0
      ta_vntiso(:,:,:) = c0
      ta_vbtiso(:,:,:) = c0
#  endif
#  if defined O_save_npzd
      ta_rnpp(:,:,:) = c0
      ta_rgraz(:,:,:) = c0
      ta_rmorp(:,:,:) = c0
      ta_rmorpt(:,:,:) = c0
      ta_rmorz(:,:,:) = c0
      ta_rexcr(:,:,:) = c0
      ta_rremi(:,:,:) = c0
      ta_rexpo(:,:,:) = c0
      ta_rexpocal(:,:,:) = c0
      ta_rprocal(:,:) = c0
#   if defined O_npzd_nitrogen
      ta_rnpp_D(:,:,:) = c0
      ta_rgraz_D(:,:,:) = c0
      ta_rmorp_D(:,:,:) = c0
      ta_rnfix(:,:,:) = c0
      ta_rdeni(:,:,:) = c0
#   endif
#  endif
#  if defined O_save_convection
      ta_totalk(:,:) = c0
      ta_vdepth(:,:) = c0
      ta_pe(:,:) = c0
#  endif
#  if defined O_save_carbon_carbonate_chem
      ta_sspH(:,:) = c0
      ta_ssCO3(:,:) = c0
      ta_ssOc(:,:) = c0
      ta_ssOa(:,:) = c0
      ta_sspCO2(:,:) = c0
#  endif
# endif
#endif

      return
      end
