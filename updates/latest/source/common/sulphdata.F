      subroutine sulphdata

#if defined O_sulphate_data || defined O_sulphate_data_transient
!=======================================================================
!     read and interpolate sulphate data
!=======================================================================

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
      include "atm.h"
      include "calendar.h"
      include "cembm.h"
      include "levind.h"
      include "tmngr.h"

      real tmpij(imtm2,jmtm2)

      name = "A_sulphod.nc"
      vname = "A_sulphod"

      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Warning => ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          time(:) = year0
          sulph(:,:,:) = 0.
          first_time = .false.
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
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
          exists = inqvardef(trim(vname), iou)
          if (.not. exists) then
            print*, "==>  Warning: A_sulpfor data does not exist."
          endif
          first_time = .true.
        endif
        iyr(:) = 0
        yrl(:) = 0.
      else
        first_time = .false.
      endif

# if defined O_sulphate_data_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      yrl(2) = min(time(ln), max(time(1), data_time))
      ice_yr = data_time
# else
      yrl(2) = min(time(ln), max(time(1), sulph_yr))
# endif

      intrp = .false.
      if (yrl(2) .gt. time(1) .and. yrl(2) .lt. time(ln)) intrp = .true.

      if (first_time .or. yrl(2) .gt. yrl(3)) then
!       read data
        ib(:) = 1
        ic(:) = imtm2
        ic(2) = jmtm2
        ic(3) = 1
        fname = new_file_name (name)
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
          print*, "=> reading sulphate data for year:",yrl(1)
          call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          sulph(2:imtm1,2:jmtm1,1) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (sulph(:,:,1))
          call openfile (fname, iou)
          ib(3) = iyr(3)
          print*, "=> reading sulphate data for year:",yrl(3)
          call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          sulph(2:imtm1,2:jmtm1,3) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (sulph(:,:,3))
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
          print*, "=> reading sulphate data for year:",yrl(2)
          call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          sulph(2:imtm1,2:jmtm1,2) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (sulph(:,:,2))
          sulph(:,:,1) = sulph(:,:,2)
          sulph(:,:,3) = sulph(:,:,2)
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
            sulph(i,j,2) = sulph(i,j,1)*wt1 + sulph(i,j,3)*wt3
          enddo
        enddo
      elseif (yrl(2) .le. time(1)) then
        sulph(:,:,2) = sulph(:,:,1)
      elseif (yrl(2) .ge. time(ln)) then
        sulph(:,:,2) = sulph(:,:,3)
      endif
      call embmbc (sulph(1,1,2))
#endif

      return
      end
