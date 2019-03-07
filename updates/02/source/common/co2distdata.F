       subroutine co2distdata
!=======================================================================
!     read and interpolate co2 distribution data
!=======================================================================
#if defined O_carbon_co2_2d
# if defined O_co2emit_data || defined O_co2emit_data_transient

      implicit none

      character(120) :: fname, name, vname, new_file_name, text

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp, exists, inqvardef

      real data_time, wt3, wt1, yrl(3), iyr(3)

      real, allocatable :: time(:)

      save time, ln, yrl, first_time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "atm.h"
      include "cembm.h"
      include "levind.h"
      include "tmngr.h"

      real dmsk(imt,jmt), tmp, tmpij(imtm2,jmtm2)

      name = "F_co2dist.nc"
# if defined O_co2emit_data_fuel && !defined O_co2emit_data_land
      vname = "F_co2dfuel"
# elif defined O_co2emit_data_land && !defined O_co2emit_data_fuel
      vname = "F_co2dland"
# else
      vname = "F_co2dist"
# endif

      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Warning => ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          time(:) = year0
          co2dist(:,:,:) = 1.
          first_time = .false.
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          exists = inqvardef(trim(vname), iou)
          if (.not. exists) then
            print*, "==>  Warning: F_co2dist data does not exist."
          endif
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
          first_time = .true.
        endif
        iyr(:) = 0
        yrl(:) = 0.
      else
        first_time = .false.
      endif

      intrp = .false.
#  if defined O_co2emit_data_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      yrl(2) = min(time(ln), max(time(1), data_time))
#  if defined O_co2emit_data_transient_repyr
      if (data_time .gt. co2_yr) then
        first_time = .true.
        intrp = .true.
        tim(2) = data_time - int(data_time) + co2_yr
        print*, "Warning => co2distdata: repeating year", co2_yr
      endif
#  endif
#  else
      yrl(2) = min(time(ln), max(time(1), co2_yr))
#  endif

      if (yrl(2) .gt. time(1) .and. yrl(2) .lt. time(ln)) intrp = .true.

      if (first_time .or. yrl(2) .gt. yrl(3)) then
!       read data
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (exists) then
          dmsk(:,:) = 1.
          if (intrp) then
            do n=2,ln
              if (time(n-1) .le. yrl(2) .and. time(n) .ge. yrl(2)) then
                yrl(1) = time(n-1)
                iyr(1) = n-1
                yrl(3) = time(n)
                iyr(3) = n
              endif
            enddo
            call openfile (fname, iou)
            ib(3) = iyr(1)
            print*, "=> reading co2dist data for year:",yrl(1)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, c1, c0)
            co2dist(2:imtm1,2:jmtm1,1) = tmpij(1:imtm2,1:jmtm2)
            call embmbc (co2dist(:,:,1))
            call areatot (co2dist(:,:,1), dmsk, tmp)
            if (tmp .ne. 0) co2dist(:,:,1) = co2dist(:,:,1)*atmsa/tmp
            ib(3) = iyr(3)
            print*, "=> reading co2dist data for year:",yrl(3)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, c1, c0)
            co2dist(2:imtm1,2:jmtm1,3) = tmpij(1:imtm2,1:jmtm2)
            call embmbc (co2dist(:,:,3))
            call areatot (co2dist(:,:,3), dmsk, tmp)
            if (tmp .ne. 0) co2dist(:,:,3) = co2dist(:,:,3)*atmsa/tmp
          else
            if (yrl(2) .le. time(1)) then
              n = 1
              yrl(1) = time(1)
              yrl(3) = time(1)
              iyr(n) = 1
            else
              n = 3
              yrl(1) = time(ln)
              yrl(3) = time(ln)
              iyr(n) = ln
            endif
            call openfile (fname, iou)
            ib(3) = iyr(n)
            print*, "=> reading co2dist data for year:",yrl(2)
            call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,       tmpij, c1, c0)
            co2dist(2:imtm1,2:jmtm1,2) = tmpij(1:imtm2,1:jmtm2)
            call embmbc (co2dist(:,:,2))
            call areatot (co2dist(:,:,2), dmsk, tmp)
            if (tmp .ne. 0) co2dist(:,:,2) = co2dist(:,:,2)*atmsa/tmp
            co2dist(:,:,1) = co2dist(:,:,2)
            co2dist(:,:,3) = co2dist(:,:,2)
          endif
        endif
      endif

      if (intrp) then
!       interpolate data
        wt1 = 1.
        if (yrl(3) .ne. yrl(1)) wt1 = (yrl(3)-yrl(2))/(yrl(3)-yrl(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        do j=1,jmt
          do i=1,imt
            co2dist(i,j,2) = co2dist(i,j,1)*wt1 + co2dist(i,j,3)*wt3
          enddo
        enddo
      elseif (yrl(2) .le. time(1)) then
        co2dist(:,:,2) = co2dist(:,:,1)
      elseif (yrl(2) .ge. time(ln)) then
        co2dist(:,:,2) = co2dist(:,:,3)
      endif
      call embmbc (co2dist(1,1,2))
# endif
#endif

      return
      end
