      subroutine mom_tbt_def (fname, imt, jmt, km, nt, xt, yt
     &,                       calendar, expnam, runstamp, mapt)

#if defined O_mom && defined O_mom_tbt
!=======================================================================
!     definition routine for ocean time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapt         = tracer map
!=======================================================================

      implicit none

      integer iou, j, n, imt, jmt, km, nt, igs, ige, ig, jgs, jge
      integer jg, kgs, kge, kg, it(10), iu(10), id_time
      integer id_xt, id_xu, id_yt, id_yu, id_zt, id_zw, id_xt_e
      integer id_xu_e, id_yt_e, id_yu_e, id_zt_e, id_zw_e

      character(*) :: fname, calendar, expnam, runstamp
      character(3) :: a3
      character(10) :: mapt(nt)

      real xt(imt), yt(jmt)
      real c0, c1, c100, c500, c1e3, c1e4, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1
      kgs = 1
      kge = km
      kg  = kge-kgs+1

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
      call defdim ('depth', iou, kg, id_zt)
      call defdim ('longitude_V', iou, ig, id_xu)
      call defdim ('latitude_V', iou, jg, id_yu)
      call defdim ('depth_W', iou, kg, id_zw)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('depth_edges', iou, kg+1, id_zt_e)
      call defdim ('longitude_V_edges', iou, ig+1, id_xu_e)
      call defdim ('latitude_V_edges', iou, jg+1, id_yu_e)
      call defdim ('depth_W_edges', iou, kg+1, id_zw_e)

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
# if defined O_units_time_years
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'years since 1-1-1')
#  else
     &, 'time', 'time', 'years since 0-1-1')
#  endif
# else
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'days since 1-1-1')
#  else
     &, 'time', 'time', 'days since 0-1-1')
#  endif
# endif
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_avgper', iou, 1, it, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','day')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      call defvar ('G_dxt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'width t grid', ' ', 'm')
      it(1) = id_yt
      call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      call defvar ('G_dyt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'height t grid', ' ', 'm')
      it(1) = id_xu
      call defvar ('longitude_V', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      call defvar ('G_dxu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'width u grid', ' ', 'm')
      it(1) = id_yu
      call defvar ('latitude_V', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      call defvar ('G_dyu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'height u grid', ' ', 'm')
      it(1) = id_zt
      call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      call defvar ('G_dzt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'thickness t grid', ' ', 'm')
      it(1) = id_zw
      call defvar ('depth_W', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'longitude', 'degrees_east')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'longitude', 'degrees_east')
      it(1) = id_zt_e
      call defvar ('depth_edges', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth edges', ' ', 'm')
      it(1) = id_zw_e
      call defvar ('depth_W_edges', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth of w grid edges', ' ', 'm')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('G_kmt', iou, 2, it, c0, c1e6, ' ', 'I'
     &, 'ocean grid depth level', 'model_level_number' ,'1')
      call defvar ('G_mskhr', iou, 2, it, c0, c1e6, ' ', 'I'
     &, 'horizontal region mask', ' ' ,'1')
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degrees_east')
      call defvar ('G_latU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid longitude', 'longitude', 'degrees_east')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nt
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n
        call defvar ('F_flux_'//trim(a3), iou ,3, it, -c1e6, c1e6, ' '
     &,   'F', 'flux tracer '//trim(a3), ' ', '?')
      enddo

!-----------------------------------------------------------------------
!     define time dependent 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_zt
      iu(3) = id_zt
      it(4) = id_time
      iu(4) = id_time
      do n=1,nt
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n
        call defvar ('O_tracer_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,   ' ', 'F', 'tracer '//trim(a3), ' ', '?')
        call defvar ('O_dt_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F',  'change in time tracer '//trim(a3),' ', '?')
        call defvar ('O_advX_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F','-zonal advection tracer '//trim(a3), ' ', '?')
        call defvar ('O_advY_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', '-meridional advection tracer '//trim(a3), ' ', '?')
        call defvar ('O_advZ_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', '-vertical advection tracer '//trim(a3), ' ', '?')
        call defvar ('O_diffX_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'zonal diffusion tracer '//trim(a3), ' ', '?')
        call defvar ('O_diffY_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'meridional diffusion tracer '//trim(a3), ' ', '?')
        call defvar ('O_diffZ_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'vertical diffusion tracer '//trim(a3), ' ', '?')
# if defined O_source_term || defined O_npzd || defined O_carbon_14
        call defvar ('O_source_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'change due to sources'//trim(a3), ' ', '?')
# endif
        call defvar ('O_convect_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'change due to convection tracer '//trim(a3), ' '
     &,   '?')
# if defined O_fourfil || defined O_firfil
        call defvar ('O_filter_'//trim(a3), iou, 4, it, -c1e6, c1e6
     &,   ' ', 'F', 'change due to filtering tracer '//trim(a3), ' '
     &,   '?')
# endif
      enddo

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mom_tbt_out (fname, ids, ide, jds, jde, imt, jmt, km
     &,                       nt, xt, yt, zt, xu, yu, zw, dxt, dyt, dzt
     &,                       dxu, dyu, dzw, avgper, time, stamp, mapt
     &,                       tbt, stf, tlat, tlon, ulat, ulon, kmt
     &,                       mskhr, tm)

!=======================================================================
!     output routine for ocean tracer balances

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if  no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths
!     avgper       = length of averaging period
!     time         = time in years
!     tbt, ...     = data to be written

!   local variables
!     igs, ige, jgs, jge = global write domain start and end indicies
!     ig, jg             = global write domain size
!     ils, ile, jls, jle = local domain start and end indicies
!     it                 = t grid axis definitions (x,y,t default)
!     iu                 = u grid axis definitions (x,y,t default)
!     is                 = start for write on each axis (x,y,t default)
!     ic                 = count for write on each axis (x,y,t default)
!     id_...             = id's for axis (used for it, iu or defvar)
!     iou                = io unit (ncid)
!     ln                 = length of data to be written
!=======================================================================

      implicit none

      integer iou, j, ln, n, ntrec, imt, jmt, km, nt, ids, ide
      integer jds, jde, igs, ige, ig, jgs, jge, jg, kgs, kge, kg
      integer ils, ile, jls, jle, kls, kle, lls, lle, ib(10)
      integer ic(10), kmt(ids:ide,jds:jde), mskhr(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      character(*) :: fname, stamp
      character(3) :: a3
      character(10) :: mapt(nt)

      real xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km), avgper
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt), dzt(km), dzw(km)
      real tbt(ids:ide,jds:jde,km,nt,11), stf(ids:ide,jds:jde,nt)
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
      real tm(ids:ide,jds:jde,km), um(ids:ide,jds:jde,km)
      real time, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), zt_e(km+1), zw_e(km+1), c0, c1, c10, c100
      real c1e3, c1e4, c1e5, c1e6, p1, p001, p035, cal2J
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpijk(:,:,:), tmpijkm(:,:,:)
      real, allocatable :: tmpijl(:,:,:), tmpijlm(:,:,:)
      real, allocatable :: tmpi(:), tmpj(:), tmpk(:), tmpl(:)
      real, allocatable :: tmpie(:), tmpje(:), tmpke(:), tmple(:)

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e5 = 1.e5
      c1e6 = 1.e6
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      cal2J = 2.389e-05

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1
      kgs = 1
      kge = km
      kg  = kge-kgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)
      kls = max(1,kgs)
      kle = min(km,kge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )
      allocate ( tmpijk(ils:ile,jls:jle,kls:kle) )
      allocate ( tmpijkm(ils:ile,jls:jle,kls:kle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, time, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)

      if (ntrec .eq. 1) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpk(kgs:kge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )
        allocate ( tmpke(kgs:kge+1) )

        ib(1) = 1
        ic(1) = ig
        tmpi(igs:ige) = xt(igs:ige)
        call putvara ('longitude', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxt(igs:ige)
        call putvara ('G_dxt', iou, ig, ib, ic, tmpi, c100, c0)
        tmpi(igs:ige) = xu(igs:ige)
        call putvara ('longitude_V', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxu(igs:ige)
        call putvara ('G_dxu', iou, ig, ib, ic, tmpi, c100, c0)

        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('latitude', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyt(jgs:jge)
        call putvara ('G_dyt', iou, jg, ib, ic, tmpj, c100, c0)
        tmpj(jgs:jge) = yu(jgs:jge)
        call putvara ('latitude_V', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyu(jgs:jge)
        call putvara ('G_dyu', iou, jg, ib, ic, tmpj, c100, c0)

        ic(1) = kg
        tmpk(kgs:kge) = zt(kgs:kge)
        call putvara ('depth', iou, kg, ib, ic, tmpk, c100, c0)
        tmpk(kgs:kge) = dzt(kgs:kge)
        call putvara ('G_dzt', iou, kg, ib, ic, tmpk, c100, c0)
        tmpk(kgs:kge) = zw(kgs:kge)
        call putvara ('depth_W', iou, kg, ib, ic, tmpk, c100, c0)

        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('longitude_edges', iou, ig+1, ib, ic, tmpie
     &,   c1, c0)
        call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xu_e(igs:ige+1)
        call putvara ('longitude_V_edges', iou, ig+1, ib, ic, tmpie
     &,   c1, c0)

        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('latitude_edges', iou, jg+1, ib, ic, tmpje
     &,   c1, c0)
        call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
        call putvara ('latitude_V_edges', iou, jg+1, ib, ic, tmpje
     &,   c1, c0)

        ic(1) = kg + 1
        call edge_maker (1, zt_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zt_e(kgs:kge+1)
        call putvara ('depth_edges', iou, kg+1, ib, ic, tmpke
     &,   c100, c0)
        call edge_maker (2, zw_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zw_e(kgs:kge+1)
        call putvara ('depth_W_edges', iou, kg+1, ib, ic, tmpke
     &,   c100, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpk )
        deallocate ( tmpie )
        deallocate ( tmpje )
        deallocate ( tmpke )

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = kmt(ils:ile,jls:jle)
        call putvara ('G_kmt', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = mskhr(ils:ile,jls:jle)
        call putvara ('G_mskhr', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('G_latT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('G_lonT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulat(ils:ile,jls:jle)
        call putvara ('G_latU', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulon(ils:ile,jls:jle)
        call putvara ('G_lonU', iou, ln, ib, ic, tmpij, c1, c0)

      endif

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      do n=1,nt
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpij(ils:ile,jls:jle) = stf(ils:ile,jls:jle,n)
        call putvaramsk('F_flux_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   tmpijm, c1, c0)
      enddo

!-----------------------------------------------------------------------
!     write 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      do n=1,nt
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,1)
        call putvaramsk('O_dt_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,2)
        call putvaramsk('O_advX_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,3)
        call putvaramsk('O_advY_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,4)
        call putvaramsk('O_advZ_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,5)
        call putvaramsk('O_diffX_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,6)
        call putvaramsk('O_diffY_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,7)
        call putvaramsk('O_diffZ_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
# if defined O_source_term || defined O_npzd || defined O_carbon_14
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,8)
        call putvaramsk('O_source_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
# endif
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,9)
        call putvaramsk('O_convect_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
# if defined O_fourfil || defined O_firfil
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,10)
        call putvaramsk('O_filter_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
# endif
        tmpijk(ils:ile,jls:jle,kls:kle) =
     &    tbt(ils:ile,jls:jle,kls:kle,n,11)
        call putvaramsk('O_tracer_'//trim(a3), iou, ln, ib, ic, tmpijk
     &,   tmpijkm, c1, c0)
      enddo

      deallocate (tmpij)
      deallocate (tmpijm)
      deallocate (tmpijk)
      deallocate (tmpijkm)
#endif

      return
      end
