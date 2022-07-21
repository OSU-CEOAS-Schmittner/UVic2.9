      subroutine co2emitdata

#if defined O_co2emit_data || defined O_co2emit_data_transient
!=======================================================================
!     routine to read and interpolate one dimensional forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer iou, n, ln, ib(10), ic(10)

      logical inqvardef, exists

      real dat(3,2), data_time, tim(3), wt1, wt3

      real, allocatable :: data(:,:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "cembm.h"
      include "switch.h"
      include "tmngr.h"

      name = "F_co2emit.nc"

      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Error: ", trim(fname), " does not exist."
          stop 'co2emitdata'
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln,2) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          text = 'years'
          call getatttext (iou, 'time', 'units', text)
          if (trim(text) .eq. "days since 1-1-1")
     &      time(:) = time(:)/yrlen - 1.
          if (trim(text) .eq. "days since 0-1-1")
     &       time(:) = time(:)/yrlen
          if (trim(text) .eq. "years since 1-1-1")
     &      time(:) = time(:) - 1.
          exists = inqvardef('F_co2efuel', iou)
          if (.not. exists) then
            print*, "==>  Warning: F_co2efuel data does not exist."
            data(:,1) = 0.
          else
            call getvara ('F_co2efuel', iou, ln, ib, ic, data(:,1)
     &,       c1, c0)
          endif
          exists = inqvardef(trim('F_co2eland'), iou)
          if (.not. exists) then
            print*, "==>  Warning: F_co2eland data does not exist."
            data(:,2) = 0.
          else
            call getvara ('F_co2eland', iou, ln, ib, ic, data(:,2)
     &,       c1, c0)
          endif
        endif
        tim(:) = time(1)
        dat(:,1) = data(1,1)
        dat(:,2) = data(1,2)
      endif

# if defined O_co2emit_data_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      tim(2) = min(time(ln), max(time(1), data_time))
# else
      tim(2) = min(time(ln), max(time(1), co2_yr))
# endif

      if (tim(2) .le. time(1)) then
        dat(2,1) = data(1,1)
        dat(2,2) = data(1,2)
      elseif (tim(2) .ge. time(ln)) then
        dat(2,1) = data(ln,1)
        dat(2,2) = data(ln,2)
      else
        if (tim(2) .gt. tim(3)) then
          do n=2,ln
            if (time(n-1) .le. tim(2) .and. time(n) .ge. tim(2)) then
              tim(1) = time(n-1)
              dat(1,1) = data(n-1,1)
              dat(1,2) = data(n-1,2)
              tim(3) = time(n)
              dat(3,1) = data(n,1)
              dat(3,2) = data(n,2)
            endif
          enddo
        endif
        wt1 = 1.
        if (tim(3) .ne. tim(1)) wt1 = (tim(3)-tim(2))/(tim(3)-tim(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        dat(2,1) = dat(1,1)*wt1 + dat(3,1)*wt3
        dat(2,2) = dat(1,2)*wt1 + dat(3,2)*wt3
      endif

!     convert flux from kg s-1 to g cm-2 s-1
      co2emit_fuel = dat(2,1)*1.e3/atmsa
      co2emit_land = dat(2,2)*1.e3/atmsa
# if defined O_co2emit_data_fuel && !defined O_co2emit_data_land
      co2emit = co2emit_fuel
# elif defined O_co2emit_data_land && !defined O_co2emit_data_fuel
      co2emit = co2emit_land
# else
      co2emit = co2emit_fuel + co2emit_land

# endif
#endif

      return
      end

      subroutine co2ccndata
#if defined O_co2ccn_data || defined O_co2ccn_data_transient || defined O_co2emit_track_co2 || defined O_co2emit_track_co2_transient
!=======================================================================
!     routine to read and interpolate one dimensional forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer iou, n, ln, ib(10), ic(10)

      logical inqvardef, exists, track

      real avg_co2, dat(3), data_time, pk, tim(3), wt1, wt3, fa

      real, allocatable :: data(:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "cembm.h"
      include "switch.h"
      include "tmngr.h"
      include "atm.h"
# if defined O_carbon_co2_2d

      real tmp, dmsk(imt,jmt)
# endif

!     fa is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      fa = 1./(4.138e-7*rhoatm*shc)
      track = .false.

# if defined O_co2emit_track_co2 && !defined O_co2emit_track_sat
!     response timescale
      pk = 1.0/(float(ntrack_co2)/4.)

# endif
      if (.not. allocated (time)) then
# if defined O_co2emit_track_co2
        name = "A_co2track.nc"
# else
        name = "A_co2.nc"
# endif
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Error: ", trim(fname), " does not exist."
          stop '=>co2ccn'
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          text = 'years'
          call getatttext (iou, 'time', 'units', text)
          if (trim(text) .eq. "days since 1-1-1")
     &      time(:) = time(:)/yrlen - 1.
          if (trim(text) .eq. "days since 0-1-1")
     &       time(:) = time(:)/yrlen
          if (trim(text) .eq. "years since 1-1-1")
     &      time(:) = time(:) - 1.
          exists = inqvardef('A_co2', iou)
          if (.not. exists) then
            print*, "==>  Warning: A_co2 data does not exist."
          else
            call getvara ('A_co2', iou, ln, ib, ic, data, c1, c0)
          endif
         endif
        tim(:) = time(1)
        dat(:) = data(1)
      endif

# if defined O_co2ccn_data_transient || defined O_co2emit_track_co2_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      tim(2) = min(time(ln), max(time(1), data_time))
# else
      tim(2) = min(time(ln), max(time(1), co2_yr))
# endif

      if (tim(2) .le. time(1)) then
        dat(2) = data(1)
      elseif (tim(2) .ge. time(ln)) then
        dat(2) = data(ln)
      else
        if (tim(2) .gt. tim(3)) then
          do n=2,ln
            if (time(n-1) .le. tim(2) .and. time(n) .ge. tim(2)) then
              tim(1) = time(n-1)
              dat(1) = data(n-1)
              tim(3) = time(n)
              dat(3) = data(n)
            endif
          enddo
        endif
        wt1 = 1.
        if (tim(3) .ne. tim(1)) wt1 = (tim(3)-tim(2))/(tim(3)-tim(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        dat(2) = dat(1)*wt1 + dat(3)*wt3
      endif

# if defined O_co2emit_track_co2

      itrack_co2 = itrack_co2 + 1
      if (itrack_co2 .gt. ntrack_co2) itrack_co2 = 1
#  if defined O_carbon_co2_2d
      dmsk(:,:) = 1.
      call areaavg (at(1,1,2,ico2), dmsk, tmp)
      co2ccn = tai_co2ccn + tmp
#  endif
      track_co2(itrack_co2) = co2ccn
      avg_co2 = 0.
      track = .true.
      do n=1,ntrack_co2
        if (track_co2(n) .ge. 1.e20) track = .false.
        avg_co2 = avg_co2 + track_co2(n)
      enddo
      avg_co2 = avg_co2/ntrack_co2
#  if !defined O_co2emit_track_sat

      if (track) then
!       use simple proportional cotrol
        co2emit = pk*(dat(2) - avg_co2)/(segtim*daylen*fa)
      endif
#  endif
# endif
# if defined O_co2emit_track_sat

#  if defined O_co2emit_track_co2
!     no tracking if "track_co2_co2emit && track_sat_co2emit"
      track = .false.
#  else
      track = .true.
      do n=1,ntrack_sat
        if (track_sat(n) .ge. 1.e20) track = .false.
      enddo

#  endif
# endif
# if defined O_co2ccn_data || defined O_co2ccn_data_transient
#  if defined O_carbon_co2_2d
      dmsk(:,:) = 1.
      call areaavg (at(1,1,1,ico2), dmsk, tmp)
      at(:,:,1,ico2) = at(:,:,1,ico2) + dat(2) - tmp
      call areaavg (at(1,1,2,ico2), dmsk, tmp)
      at(:,:,2,ico2) = at(:,:,2,ico2) + dat(2) - tmp
#  endif
      co2ccn = dat(2)
# else
      if (.not. track) then
#  if defined O_carbon_co2_2d
        dmsk(:,:) = 1.
        call areaavg (at(1,1,1,ico2), dmsk, tmp)
        at(:,:,1,ico2) = at(:,:,1,ico2) + dat(2) - tmp
        call areaavg (at(1,1,2,ico2), dmsk, tmp)
        at(:,:,2,ico2) = at(:,:,2,ico2) + dat(2) - tmp
#  endif
        co2ccn = dat(2)
      endif
# endif
#endif

      return
      end

      subroutine satdata
#if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient || defined O_embm_vcs
!=======================================================================
!     routine to read and interpolate one dimensional forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer iou, n, ln, ib(10), ic(10)

      logical inqvardef, exists, track

      real avg_sat, dat(3), data_time, fa, pk, tim(3)
      real tmp, wt1, wt3

      real, allocatable :: data(:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "cembm.h"
      include "atm.h"
# if defined O_landice_data
#  if defined O_ice_cpts
      include "cpts.h"
#  endif
      include "ice.h"
# endif
      include "switch.h"
      include "tmngr.h"

      real dmsk(imt,jmt), sat(imt,jmt)

# if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient
!     fa is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      fa = 1./(4.138e-7*rhoatm*shc)

!     temperature to CO2 conversion divided by response time scale
!     (approximate CO2 concentration divided by climate sensitivity)
      pk = 300./3.0/(float(ntrack_sat)/4.)
      track = .true.

      if (.not. allocated (time)) then
        name = "A_sattrack.nc"
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "==> Error: ", trim(fname), " does not exist."
          stop '=>satdata'
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          text = 'years'
          call getatttext (iou, 'time', 'units', text)
          if (trim(text) .eq. "days since 1-1-1")
     &      time(:) = time(:)/yrlen - 1.
          if (trim(text) .eq. "days since 0-1-1")
     &       time(:) = time(:)/yrlen
          if (trim(text) .eq. "years since 1-1-1")
     &      time(:) = time(:) - 1.
          exists = inqvardef('A_sat', iou)
          if (.not. exists) then
            print*, "==>  Warning: A_sat data does not exist."
          else
            call getvara ('A_sat', iou, ln, ib, ic, data, c1, c0)
            text = "C"
            call getatttext (iou, 'A_sat', 'units', text)
!           convert to model units (C)
            if (trim(text) .eq. "K")
     &      where (data(:) .lt. 1.e30) data(:) = data(:) - 273.15
          endif
        endif
        tim(:) = time(1)
        dat(:) = data(1)
      endif

# if defined O_co2emit_track_sat_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      tim(2) = min(time(ln), max(time(1), data_time))
# else
      tim(2) = min(time(ln), max(time(1), co2_yr))
# endif

      if (tim(2) .le. time(1)) then
        dat(2) = data(1)
      elseif (tim(2) .ge. time(ln)) then
        dat(2) = data(ln)
      else
        if (tim(2) .gt. tim(3)) then
          do n=2,ln
            if (time(n-1) .le. tim(2) .and. time(n) .ge. tim(2)) then
              tim(1) = time(n-1)
              dat(1) = data(n-1)
              tim(3) = time(n)
              dat(3) = data(n)
            endif
          enddo
        endif
        wt1 = 1.
        if (tim(3) .ne. tim(1)) wt1 = (tim(3)-tim(2))/(tim(3)-tim(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        dat(2) = dat(1)*wt1 + dat(3)*wt3
      endif
# endif

      sat(:,:) = at(:,:,2,isat) - elev(:,:)*rlapse
# if defined O_landice_data
     &         - hicel(:,:,2)*rlapse
# endif
# if defined O_sealev || defined O_sealev_data
     &         - elev_sealev(:,:)*rlapse
# endif
      dmsk(:,:) = 1.
      call areaavg (sat, dmsk, tmp)

      itrack_sat = itrack_sat + 1
      if (itrack_sat .gt. ntrack_sat) itrack_sat = 1
      track_sat(itrack_sat) = tmp

      avg_sat = 0.
      do n=1,ntrack_sat
        if (track_sat(n) .ge. 1.e10) track = .false.
        avg_sat = avg_sat + track_sat(n)
      enddo
      avg_sat = avg_sat/ntrack_sat

# if defined O_co2emit_track_co2
!     no tracking if "track_co2_co2emit && track_sat_co2emit"
      track = .false.
# endif
# if defined O_co2emit_track_sat
      if (track) then
!       use simple proportional control
        co2emit = pk*(dat(2) - avg_sat)/(segtim*daylen*fa)
      endif
# endif
#endif

      return
      end
