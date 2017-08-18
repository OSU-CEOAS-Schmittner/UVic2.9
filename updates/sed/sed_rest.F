      subroutine sed_rest_in (fname, ids, ide, jds, jde)

#if defined O_sed
!=======================================================================
!     input routine for sediment model restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================

      implicit none

      character(*) :: fname
      character(32) :: nstamp
      character(3) :: a3
      character(120) :: var1, var2

      integer i, iou, j, ln, ntrec, ids, ide, jds, jde, l, m, n
      integer ils, ile, jls, jle, lls, lle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical inqvardef, warning

      real tmp
      real, allocatable :: tmpij(:,:)
      real, allocatable :: tmpip(:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "tmngr.h"
      include "switch.h"
      include "csbc.h"
      include "sed.h"

      integer map_tmp(imt,jmt)

      real data(imt,jmt)

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1
      if (.not. inqvardef('atsed', iou)) then
        print*, "==> Warning: sediment model restart not found"
        return
      endif

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global read domain)
!-----------------------------------------------------------------------
      lls = 1
      lle = ipmax
      allocate ( tmpip(lls:lle) )

      ils = max(ids,1)
      ile = min(ide,imt)
      jls = max(jds,1)
      jle = min(jde,jmt)
      allocate ( tmpij(ils:ile,jls:jle) )

      map_tmp(:,:) = map_sed(:,:)

!-----------------------------------------------------------------------
!     read 1d data (t)
!-----------------------------------------------------------------------
      tmp = itt
      call getvars ('itt', iou, ntrec, tmp, c1, c0)
      itt = tmp
      tmp = irstdy
      call getvars ('irstdy', iou, ntrec, tmp, c1, c0)
      irstdy = tmp
      tmp = msrsdy
      call getvars ('msrsdy', iou, ntrec, tmp, c1, c0)
      msrsdy = tmp
      tmp = year0
      call getvars ('year', iou, 1, tmp, c1, c0)
      nyear = tmp
      tmp = month0
      call getvars ('month', iou, 1, tmp, c1, c0)
      nmonth = tmp
      tmp = day0
      call getvars ('day', iou, 1, tmp, c1, c0)
      nday = tmp
      tmp = hour0
      call getvars ('hour', iou, 1, tmp, c1, c0)
      nhour = tmp
      tmp = min0
      call getvars ('minute', iou, 1, tmp, c1, c0)
      nmin = tmp
      tmp = sec0
      call getvars ('second', iou, 1, tmp, c1, c0)
      nsec = tmp
      call mkstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = carblith
      call getvars ('carblith', iou, 1, tmp, c1, c0)
      carblith = tmp
      if (init_time_in) then
        itt = 0
        irstdy = 0
        msrsdy = 0
        relyr = 0.
        call mkstmp (stamp, year0, month0, day0, hour0, min0, sec0)
        carblith = 0.
      endif
      tmp = atsed
      call getvars ('atsed', iou, 1, tmp, c1, c0)
      atsed = tmp
      tmp = weathflx
      call getvars ('weathflx', iou, 1, tmp, c1, c0)
      weathflx = tmp
      tmp = sed_year
      call getvars ('sed_year', iou, 1, tmp, c1, c0)
      sed_year = tmp

!-----------------------------------------------------------------------
!     read 2d data (x,y)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ln = ic(1)*ic(2)

      tmpij(:,:) = kmt(ils:ile,jls:jle)
      call getvara ('kmt', iou, ln, ib, ic, tmpij, c1, c0)

      map_tmp(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call getvara ('map_sed', iou, ln, ib, ic, tmpij, c1, c0)
      map_tmp(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      warning = .false.
      do j=jls,jle
        do i=ils,ile
          if (map_tmp(i,j) .ne. map_sed(i,j))  warning = .true.
        enddo
      enddo
      if (warning) then
        print*, ' '
        print*, '==> Warning: restart sediment mask is different.'
      endif

!-----------------------------------------------------------------------
!     read 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = 1
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibtemp)
      call getvara ('sbc_btemp', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibtemp) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibsalt)
      call getvara ('sbc_bsalt', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibsalt) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ircal)
      call getvara ('sbc_rcal', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ircal) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,irorg)
      call getvara ('sbc_rorg', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,irorg) = tmpij(ils:ile,jls:jle)
# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibdic)
      call getvara ('sbc_bdic', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibdic) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibdicfx)
      call getvara ('sbc_bdicfx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibdicfx) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_npzd_alk
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibalk)
      call getvara ('sbc_balk', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibalk) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibalkfx)
      call getvara ('sbc_balkfx', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibalkfx) = tmpij(ils:ile,jls:jle)
# endif
# if defined O_npzd_o2
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibo2)
      call getvara ('sbc_bo2', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ibo2) = tmpij(ils:ile,jls:jle)
# endif

      call unloadsed (ipmax, zrct, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_zrct', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      zrct(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, water_z_p, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_water_z_p', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      water_z_p(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, k1, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_k1', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      k1(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, k2, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_k2', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      k2(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, k3, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_k3', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      k3(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, csat, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_csat', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      csat(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, rc, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_rc', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      rc(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, ttrorg, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_ttrorg', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      ttrorg(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, ttrcal, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_ttrcal', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      ttrcal(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, sed_ml_mass, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_ml_mass', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      sed_ml_mass(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, c_advect, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_c_advect', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      c_advect(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, rain_cal_p, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_rain_cal_p', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      rain_cal_p(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, rain_org_p, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_rain_org_p', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      rain_org_p(1:ipmax) = tmpip(1:ipmax)

      call unloadsed (ipmax, co3_p, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('sed_co3_p', iou, ln, ib, ic, tmpij, c1, c0)
      where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      co3_p(1:ipmax) = tmpip(1:ipmax)

!-----------------------------------------------------------------------
!     read 4d data (x,y,z,t)
!-----------------------------------------------------------------------

      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = 1
      ic(3) = 1
      ib(4) = 1
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)

      n = 0
      call getdimlen ('slev', iou, n)
      n = min(n, nzmax)
      do l=1,n
        ib(3) = l

        do m=1,3
          write(a3, '(i1)') m

          tmpip(1:ipmax) = carb(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call getvara ('sed_carb_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
          call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          carb(l,m,1:ipmax) = tmpip(1:ipmax)

          tmpip(1:ipmax) = dcpls(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call getvara ('sed_dcpls_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
          call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          dcpls(l,m,1:ipmax) = tmpip(1:ipmax)

          tmpip(1:ipmax) = dcmin(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call getvara ('sed_dcmin_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
          call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          dcmin(l,m,1:ipmax) = tmpip(1:ipmax)

        enddo

        tmpip(1:ipmax) = pore(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_pore', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        pore(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = form(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_form', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        form(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = o2(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_o2', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        o2(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = orggg(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_orggg', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        orggg(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = orgml(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_orgml', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        orgml(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = calgg(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_calgg', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        calgg(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = calml(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_calml', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        calml(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = dopls(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_dopls', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        dopls(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = domin(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_domin', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        domin(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = dbpls(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_dbpls', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        dbpls(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = dbmin(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_dbmin', iou, ln, ib, ic, tmpij, c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        dbmin(l,1:ipmax) = tmpip(1:ipmax)

      enddo

      n = 0
      call getdimlen ('blev', iou, n)
      n = min(n, ibmax)
      do l=1,n
        ib(3) = l

        tmpip(1:ipmax) = buried_mass(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_buried_mass', iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        buried_mass(l,1:ipmax) = tmpip(1:ipmax)

        tmpip(1:ipmax) = buried_calfrac(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('sed_buried_calfrac', iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (map_tmp(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        buried_calfrac(l,1:ipmax) = tmpip(1:ipmax)

      enddo

      deallocate (tmpip)

      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Sed restart read from ',trim(fname),' on ', nstamp

      deallocate (tmpij)

      return
      end

      subroutine sed_rest_def (fname)
!=======================================================================
!     definition routine for sediment model restarts

!   inputs:
!     fname = file name
!=======================================================================

      implicit none

      character(*) :: fname
      character(3) :: a3

      integer iou, ntrec, igs, ige, ig, jgs, jge, jg, m
      integer it(10), id_time, id_xt, id_yt, id_xt_e, id_yt_e
      integer id_slev, id_slev_e, id_blev, id_blev_e

      real c100, c1e3, c1e6, c1e20

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "iounit.h"
      include "mw.h"
      include "tmngr.h"
      include "sed.h"

      c100 = 100.
      c1e3 = 1.e3
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global atributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('longitude', iou, ig, id_xt)
      call defdim ('latitude', iou, jg, id_yt)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('slev', iou, nzmax, id_slev)
      call defdim ('slev_edges', iou, nzmax+1, id_slev_e)
      call defdim ('blev', iou, ibmax, id_blev)
      call defdim ('blev_edges', iou, ibmax+1, id_blev_e)

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'years since 0-1-1')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('itt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itt', ' ',' ')
      call defvar ('irstdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'irstdy', ' ',' ')
      call defvar ('msrsdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'msrsdy', ' ',' ')
      call defvar ('year', iou, 1, it, c0, c0, ' ', 'D'
     &, 'year', ' ',' ')
      call defvar ('month', iou, 1, it, c0, c0, ' ', 'D'
     &, 'month', ' ',' ')
      call defvar ('day', iou, 1, it, c0, c0, ' ', 'D'
     &, 'day', ' ',' ')
      call defvar ('hour', iou, 1, it, c0, c0, ' ', 'D'
     &, 'hour', ' ',' ')
      call defvar ('minute', iou, 1, it, c0, c0, ' ', 'D'
     &, 'minute', ' ',' ')
      call defvar ('second', iou, 1, it, c0, c0, ' ', 'D'
     &, 'second', ' ',' ')
      call defvar ('carblith', iou, 1, it, c0, c0, ' ', 'D'
     &, 'carblith', ' ',' ')
      call defvar ('atsed', iou, 1, it, c0, c0, ' ', 'D'
     &, 'atsed', ' ',' ')
      call defvar ('weathflx', iou, 1, it, c0, c0, ' ', 'D'
     &, 'weathflx', ' ',' ')
      call defvar ('sed_year', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sed_year', ' ',' ')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yt
      call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')
      it(1) = id_slev
      call defvar ('slev', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sediment level', ' ','1')
      it(1) = id_slev_e
      call defvar ('slev_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'sediment level edges', ' ','1')
      it(1) = id_blev
      call defvar ('blev', iou, 1, it, c0, c0, ' ', 'D'
     &, 'burial level', ' ','1')
      it(1) = id_blev_e
      call defvar ('blev_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'burial level edges', ' ','1')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      call defvar ('kmt', iou, 2, it, c0, c1e3, ' ', 'I'
     &,  'kmt', ' ' ,' ')
      call defvar ('map_sed', iou , 2, it, c0, c1e6, ' ', 'I'
     &, 'map_sed', ' ', ' ')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_time

      call defvar ('sbc_btemp', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc btemp', ' ', ' ')
      call defvar ('sbc_bsalt', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc bsalt', ' ', ' ')
      call defvar ('sbc_rcal', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc rcal', ' ', ' ')
      call defvar ('sbc_rorg', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc rorg', ' ', ' ')
# if defined O_carbon
      call defvar ('sbc_bdic', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc bdic', ' ', ' ')
      call defvar ('sbc_bdicfx', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc bdicfx', ' ', ' ')
# endif
# if defined O_npzd_alk
      call defvar ('sbc_balk', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc balk', ' ', ' ')
      call defvar ('sbc_balkfx', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc balkfx', ' ', ' ')
# endif
# if defined O_npzd_o2
      call defvar ('sbc_bo2', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sbc bo2', ' ', ' ')
# endif
      call defvar ('sed_zrct', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'zrct', ' ', ' ')
      call defvar ('sed_water_z_p', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'water_z_p', ' ', ' ')
      call defvar ('sed_k1', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'k1', ' ', ' ')
      call defvar ('sed_k2', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'k2', ' ', ' ')
      call defvar ('sed_k3', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'k3', ' ', ' ')
      call defvar ('sed_csat', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'csat', ' ', ' ')
      call defvar ('sed_rc', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'rc', ' ', ' ')
      call defvar ('sed_ttrorg', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'ttrorg', ' ', ' ')
      call defvar ('sed_ttrcal', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'ttrcal', ' ', ' ')
      call defvar ('sed_ml_mass', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'sed_ml_mass', ' ', ' ')
      call defvar ('sed_c_advect', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'c_advect', ' ', ' ')
      call defvar ('sed_rain_cal_p', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'rain_cal_p', ' ', ' ')
      call defvar ('sed_rain_org_p', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'rain_org_p', ' ', ' ')
      call defvar ('sed_co3_p', iou, 3, it,  -c1e20, c1e20, ' '
     &, 'D', 'co3_p', ' ', ' ')

!-----------------------------------------------------------------------
!     define 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_slev
      it(4) = id_time

      do m=1,3
        write(a3, '(i1)') m

        call defvar ('sed_carb_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,   ' ', 'D', 'sed_carb '//trim(a3), ' ', ' ')

        call defvar ('sed_dcpls_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,   ' ', 'D', 'sed_dcpls_'//trim(a3), ' ', ' ')

        call defvar ('sed_dcmin_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,   ' ', 'D', 'sed_dcmin_'//trim(a3), ' ', ' ')

      enddo

      call defvar ('sed_pore', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_pore', ' ', ' ')

      call defvar ('sed_form', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_form', ' ', ' ')

      call defvar ('sed_o2', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_o2', ' ', ' ')

      call defvar ('sed_orggg', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_orggg', ' ', ' ')

      call defvar ('sed_orgml', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_orgml', ' ', ' ')

      call defvar ('sed_calgg', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_calgg', ' ', ' ')

      call defvar ('sed_calml', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_calml', ' ', ' ')

      call defvar ('sed_dopls', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_dopls', ' ', ' ')

      call defvar ('sed_domin', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_domin', ' ', ' ')

      call defvar ('sed_dbpls', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_dbpls', ' ', ' ')

      call defvar ('sed_dbmin', iou, 4, it,  -c1e20, c1e20, ' ', 'D'
     &, 'sed_dbmin', ' ', ' ')

      it(3) = id_blev

      call defvar ('sed_buried_mass', iou, 4, it,  -c1e20, c1e20
     &, ' ', 'D', 'buried_mass', ' ', ' ')

      call defvar ('sed_buried_calfrac', iou, 4, it,  -c1e20, c1e20
     &, ' ', 'D', 'buried_calfrac', ' ', ' ')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine sed_rest_out (fname, ids, ide, jds, jde)
!=======================================================================
!     output routine for sediment model restarts

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname              = file name
!     ids, ide ...       = start and end index for data domain
!=======================================================================

      implicit none

      character(*) :: fname
      character(3) :: a3
      character(32) :: nstamp
      character(120) :: var1, var2

      integer i, iou, j, ln, l, m, n, ntrec, ids, ide, jds, jde, igs
      integer ige, ig, jgs, jge, jg, ils, ile, jls, jle, lls, lle
      integer ib(10), ic(10), id_xt, id_yt, nyear, nmonth, nday, nhour
      integer nmin, nsec

      real tmp, c100, c1e3, c1e20
      real, allocatable :: tmpij(:,:), tmpi(:), tmpj(:), tmpie(:)
      real, allocatable :: tmpje(:), tmpip(:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "switch.h"
      include "tmngr.h"
      include "csbc.h"
      include "sed.h"

      real xt_e(imt+1), yt_e(jmt+1), data(imt,jmt), slev(nzmax)
      real slev_e(nzmax+1), blev(ibmax), blev_e(ibmax+1)

      lls = 1
      lle = ipmax

      allocate ( tmpip(lls:lle) )

      c100 = 100.
      c1e3 = 1.e3
      c1e20 = 1.e20
      nstamp = stamp

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpi(igs:ige) )
      allocate ( tmpj(jgs:jge) )
      allocate ( tmpie(igs:ige+1) )
      allocate ( tmpje(jgs:jge+1) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      if (init_time_out) then
        tmp = 0.
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = 0.
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
        call mkstmp (nstamp, year0, month0, day0, hour0, min0, sec0)
        carblith = 0.
      else
        tmp = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call putvars ('time', iou, ntrec, tmp, c1, c0)
        tmp = itt
        call putvars ('itt', iou, ntrec, tmp, c1, c0)
        tmp = iday(imodeltime)
        call putvars ('irstdy', iou, ntrec, tmp, c1, c0)
        tmp = msday(imodeltime)
        call putvars ('msrsdy', iou, ntrec, tmp, c1, c0)
      endif
      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('year', iou, 1, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, 1, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, 1, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, 1, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, 1, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, 1, tmp, c1, c0)
      tmp = carblith
      call putvars ('carblith', iou, 1, tmp, c1, c0)
      tmp = atsed
      call putvars ('atsed', iou, 1, tmp, c1, c0)
      tmp = weathflx
      call putvars ('weathflx', iou, 1, tmp, c1, c0)
      tmp = sed_year
      call putvars ('sed_year', iou, 1, tmp, c1, c0)

!-----------------------------------------------------------------------
!     write 1d data (x, y or z)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ig
      tmpi(igs:ige) = xt(igs:ige)
      call putvara ('longitude', iou, ig, ib, ic, tmpi, c1, c0)
      tmpi(igs:ige) = xu(igs:ige)
      call putvara ('longitude_V', iou, ig, ib, ic, tmpi, c1, c0)

      ic(1) = jg
      tmpj(jgs:jge) = yt(jgs:jge)
      call putvara ('latitude', iou, jg, ib, ic, tmpj, c1, c0)
      tmpj(jgs:jge) = yu(jgs:jge)
      call putvara ('latitude_V', iou, jg, ib, ic, tmpj, c1, c0)

      ic(1) = ig + 1
      call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xt_e(igs:ige+1)
      call putvara ('longitude_edges', iou, ig+1, ib, ic, tmpie, c1, c0)

      ic(1) = jg + 1
      call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
      call putvara ('latitude_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

      ib(1) = 1
      do l=1,nzmax
        slev(l) = float(l)
        slev_e(l) = float(l) - 0.5
      enddo
      slev_e(nzmax+1) = slev(nzmax) + 0.5
      ic(1) = nzmax
      call putvara ('slev', iou, nzmax, ib, ic, slev, 1., 0.)
      ic(1) = nzmax+1
      call putvara ('slev_edges', iou, nzmax+1, ib, ic, slev_e, 1., 0.)
      ib(1) = 1
      do l=1,ibmax
        blev(l) = float(l)
        blev_e(l) = float(l) - 0.5
      enddo
      blev_e(ibmax+1) = blev(ibmax) + 0.5
      ic(1) = ibmax
      call putvara ('blev', iou, ibmax, ib, ic, blev, 1., 0.)
      ic(1) = ibmax+1
      call putvara ('blev_edges', iou, ibmax+1, ib, ic, blev_e
     &, 1., 0.)

!-----------------------------------------------------------------------
!     write 2d data (x,y)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ln = ic(1)*ic(2)
      tmpij(ils:ile,jls:jle) = kmt(ils:ile,jls:jle)
      call putvara ('kmt', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = map_sed(ils:ile,jls:jle)
      call putvara ('map_sed', iou, ln, ib, ic, tmpij, c1, c0)

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibtemp)
      call putvara ('sbc_btemp', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibsalt)
      call putvara ('sbc_bsalt', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ircal)
      call putvara ('sbc_rcal', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,irorg)
      call putvara ('sbc_rorg', iou, ln, ib, ic, tmpij, c1, c0)
# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibdic)
      call putvara ('sbc_bdic', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibdicfx)
      call putvara ('sbc_bdicfx', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_npzd_alk
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibalk)
      call putvara ('sbc_balk', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibalkfx)
      call putvara ('sbc_balkfx', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_npzd_o2
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ibo2)
      call putvara ('sbc_bo2', iou, ln, ib, ic, tmpij, c1, c0)
# endif

      tmpip(1:ipmax) = zrct(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_zrct', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = water_z_p(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_water_z_p', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = k1(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_k1', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = k2(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_k2', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = k3(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_k3', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = csat(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_csat', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = rc(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_rc', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = ttrorg(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_ttrorg', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = ttrcal(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_ttrcal', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = sed_ml_mass(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_ml_mass', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = c_advect(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_c_advect', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = rain_cal_p(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_rain_cal_p', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = rain_org_p(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_rain_org_p', iou, ln, ib, ic, tmpij, c1, c0)

      tmpip(1:ipmax) = co3_p(1:ipmax)
      call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('sed_co3_p', iou, ln, ib, ic, tmpij, c1, c0)

!-----------------------------------------------------------------------
!     write 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = 1
      ic(3) = 1
      ib(4) = 1
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)

      do l=1,nzmax
        ib(3) = l

        do m=1,3
          write(a3, '(i1)') m

          tmpip(1:ipmax) = carb(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call putvara ('sed_carb_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)

          tmpip(1:ipmax) = dcpls(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call putvara ('sed_dcpls_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)

          tmpip(1:ipmax) = dcmin(l,m,1:ipmax)
          call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
          tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
          call putvara ('sed_dcmin_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)

        enddo

        tmpij(ils:ile,jls:jle) = map_sed(ils:ile,jls:jle)
        call putvara ('map_sed', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = pore(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_pore', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = form(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_form', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = o2(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_o2', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = orggg(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_orggg', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = orgml(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_orgml', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = calgg(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_calgg', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = calml(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_calml', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = dopls(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_dopls', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = domin(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_domin', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = dbpls(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_dbpls', iou, ln, ib, ic, tmpij, c1, c0)

        tmpip(1:ipmax) = dbmin(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_dbmin', iou, ln, ib, ic, tmpij, c1, c0)

      enddo

      do l=1,ibmax
        ib(3) = l

        tmpip(1:ipmax) = buried_mass(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_buried_mass', iou, ln, ib, ic, tmpij
     &,   c1, c0)

        tmpip(1:ipmax) = buried_calfrac(l,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('sed_buried_calfrac', iou, ln, ib, ic, tmpij
     &,   c1, c0)

      enddo

      deallocate (tmpip)

      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Sed restart written to ',trim(fname),' on ', nstamp

      deallocate ( tmpij )
      deallocate ( tmpi )
      deallocate ( tmpj )
      deallocate ( tmpie )
      deallocate ( tmpje )
#endif

      return
      end
