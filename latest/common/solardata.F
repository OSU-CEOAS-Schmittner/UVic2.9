      subroutine solardata

#if defined O_solar_data || defined O_solar_data_transient
!=======================================================================
!     routine to read and interpolate solar forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer iou, n, ln, ib(10), ic(10)

      logical exists

      real dat(3), data_time, tim(3)
      real c1e3, wt1, wt3

      real, allocatable :: data(:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "cembm.h"
      include "tmngr.h"

      c1e3 = 1.e3
      if (.not. allocated (time)) then
        fname = "A_solarcon.nc"
        name = new_file_name (fname)
        inquire (file=trim(name), exist=exists)
        if (.not. exists) then
          print*, "==> Warning: ", trim(name), " does not exist."
          ln = 3
          allocate ( time(ln) )
          allocate ( data(ln) )
          time(:) = year0
          data(:) = solarconst
        else
          call openfile (name, iou)
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
          if (.not. exists) then
            print*, "==>  Warning: A_solarcon data does not exist."
          else
            call getvara ('A_solarcon', iou, ln, ib, ic, data
     &,       c1e3, c0)
          endif
        endif
        tim(:) = time(1)
        dat(:) = data(1)
      endif

# if defined O_solar_data_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      tim(2) = min(time(ln), max(time(1), data_time))
#  if defined O_solar_data_transient_repyr
      if (data_time .gt. solar_yr) then
        first_time = .true.
        intrp = .true.
        tim(2) = data_time - int(data_time) + solar_yr
        print*, "Warning => solardata: repeating year", solar_yr
      endif
#  endif
# else
      tim(2) = min(time(ln), max(time(1), solar_yr))
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

      solarconst = dat(2)
#endif

      return
      end
