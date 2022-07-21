       subroutine icedata

#if defined O_landice_data || defined O_landice_data_transient
!=======================================================================
!     read and interpolate ice data
!     for a specific data set that gives ice thickness (relative to
!     present day) and area every 1000 years from -19000 to 2000
!     calendar years (21kbp to 0kbp)
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp, exists, inqvardef

      real data_time, wt3, wt1, c100, yrl(3), iyr(3)

      real, allocatable :: time(:)

      save time, ln, yrl, first_time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "atm.h"
      include "calendar.h"
      include "cembm.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
      include "levind.h"
      include "tmngr.h"

      real tmpij(imtm2,jmtm2)

      c100 = 100.

!     recalculate land masks to be consistent with aicel
      fname = new_file_name ("G_mskt.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('G_mskt', iou)
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        if (exists) call getvara ('G_mskt', iou, imtm2*jmtm2, ib, ic
     &,   tmpij, c1, c0)
        tmsk(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
      else
        tmsk(:,:) = 0.
        do j=2,jmtm1
          do i=2,imtm1
            if (kmt(i,j) .gt. 0.) tmsk(i,j) = 1.
          enddo
        enddo
      endif
      call embmbc (tmsk)

      name = "L_ice.nc"

      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Warning => ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          time(:) = year0
          aicel(:,:,:) = 0.
          hicel(:,:,:) = 0.
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
          first_time = .true.
        endif
        iyr(:) = 0
        yrl(:) = 0.
      else
        first_time = .false.
      endif

      intrp = .false.
# if defined O_landice_data_transient
      data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      yrl(2) = min(time(ln), max(time(1), data_time))
#  if defined O_landice_data_transient_repyr
      if (data_time .gt. ice_yr) then
        first_time = .true.
        intrp = .true.
        tim(2) = data_time - int(data_time) + ice_yr
        print*, "Warning => icedata: repeating year", ice_yr
      endif
#  endif
# else
      yrl(2) = min(time(ln), max(time(1), ice_yr))
# endif

      if (yrl(2) .gt. time(1) .and. yrl(2) .lt. time(ln)) intrp = .true.

      if (first_time .or. yrl(2) .gt. yrl(3)) then
!       read data
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
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
          print*, "=> reading ice data for year:",yrl(1)
          call getvara ('L_icethk', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c100, c0)
          hicel(2:imtm1,2:jmtm1,1) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (hicel(:,:,1))
          call getvara ('L_icefra', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          aicel(2:imtm1,2:jmtm1,1) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (aicel(:,:,1))
          call openfile (fname, iou)
          ib(3) = iyr(3)
          print*, "=> reading ice data for year:",yrl(3)
          call getvara ('L_icethk', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c100, c0)
          hicel(2:imtm1,2:jmtm1,3) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (hicel(:,:,3))
          call getvara ('L_icefra', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          aicel(2:imtm1,2:jmtm1,3) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (aicel(:,:,3))
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
          print*, "=> reading ice data for year:",yrl(2)
          call getvara ('L_icethk', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c100, c0)
          hicel(2:imtm1,2:jmtm1,2) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (hicel(:,:,2))
          call getvara ('L_icefra', iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          aicel(2:imtm1,2:jmtm1,2) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (aicel(:,:,2))
          hicel(:,:,1) = hicel(:,:,2)
          hicel(:,:,3) = hicel(:,:,2)
          aicel(:,:,1) = aicel(:,:,2)
          aicel(:,:,3) = aicel(:,:,2)
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
            hicel(i,j,2) = hicel(i,j,1)*wt1 + hicel(i,j,3)*wt3
            aicel(i,j,2) = aicel(i,j,1)*wt1 + aicel(i,j,3)*wt3
            if (aicel(i,j,2) .lt. 0.5) then
              aicel(i,j,2) = 0.
              hicel(i,j,2) = 0.
            else
              aicel(i,j,2) = 1.
            endif
          enddo
        enddo
      endif
      call embmbc (hicel(1,1,2))
      call embmbc (aicel(1,1,2))

      do j=1,jmtm1
        do i=1,imtm1
          if (aicel(i,j,2) .ge. 0.5) tmsk(i,j) = 0.
        enddo
      enddo
      call embmbc (tmsk)
#endif

      return
      end
