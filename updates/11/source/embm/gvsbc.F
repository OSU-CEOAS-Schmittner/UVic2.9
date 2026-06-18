      subroutine gvsbc

#if defined O_embm
!=======================================================================
!     calculates albedo and dalton numbers over vegetation
!     may read and interpolate agricutural land data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name, text

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp, inqvardef, exists

      real data_time, wt3, wt1, z0, yrv(3), iyr(3)

      real, allocatable :: time(:)

      save time, ln, yrv, first_time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
      include "veg.h"
      include "levind.h"
      include "tmngr.h"

# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
      real tmpij(imtm2,jmtm2)

      fname = "L_agricfra.nc"
      fname = new_file_name (fname)
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        fname = "L_cropfra.nc"
        fname = new_file_name (fname)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, 'Warning => agriculture fraction file does not exist.'
           agric(:,:,:) = 0.
        endif
      endif

      if (exists) then
        if (.not. allocated (time)) then
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
          iyr(:) = 0
          yrv(:) = 0.
          agric(:,:,:) = 0.
        else
          first_time = .false.
        endif
        intrp = .false.

#  if defined O_crop_data_transient || defined O_pasture_data_transient || defined O_agric_data_transient
        data_time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        yrv(2) = min(time(ln), max(time(1), data_time))
#  else
        yrv(2) = min(time(ln), max(time(1), agric_yr))
#  endif
        if (yrv(2).gt.time(1) .and. yrv(2).lt.time(ln)) intrp = .true.

        if (first_time .or. yrv(2) .gt. yrv(3)) then
!         read data
          fname = new_file_name (fname)
          ib(:) = 1
          ic(:) = 1
          ic(1) = imtm2
          ic(2) = jmtm2
          if (intrp) then
            do n=2,ln
              if (time(n-1) .le. yrv(2) .and. time(n) .ge. yrv(2)) then
                yrv(1) = time(n-1)
                iyr(1) = n-1
                yrv(3) = time(n)
                iyr(3) = n
              endif
            enddo

            call openfile (fname, iou)
            ib(3) = iyr(1)
            agric(2:imtm1,2:jmtm1,1) = 0.
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
            exists = inqvardef('L_cropfra', iou)
            if (exists) then
              call getvara ('L_cropfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,1) = agric(2:imtm1,2:jmtm1,1)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
# if defined O_pasture_data || defined O_agric_data
            exists = inqvardef('L_pastfra', iou)
            if (exists) then
              call getvara ('L_pastfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,1) = agric(2:imtm1,2:jmtm1,1)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
            call embmbc (agric(:,:,1))

            call openfile (fname, iou)
            ib(3) = iyr(3)
            agric(2:imtm1,2:jmtm1,3) = 0.
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
            exists = inqvardef('L_cropfra', iou)
            if (exists) then
              call getvara ('L_cropfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,3) = agric(2:imtm1,2:jmtm1,3)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
            exists = inqvardef('L_pastfra', iou)
            if (exists) then
              call getvara ('L_pastfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,3) = agric(2:imtm1,2:jmtm1,3)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
            call embmbc (agric(1,1,3))

          else
            if (yrv(2) .le. time(1)) then
              n = 1
              yrv(1) = time(1)
              yrv(3) = time(1)
              iyr(n) = 1
            else
              n = 3
              yrv(1) = time(ln)
              yrv(3) = time(ln)
              iyr(n) = ln
            endif
            call openfile (fname, iou)
            ib(3) = iyr(n)
            print*, "=> reading agriculture data for year:",yrv(n)
            agric(2:imtm1,2:jmtm1,2) = 0.
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
            exists = inqvardef('L_cropfra', iou)
            if (exists) then
              call getvara ('L_cropfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,2) = agric(2:imtm1,2:jmtm1,2)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
            exists = inqvardef('L_pastfra', iou)
            if (exists) then
              call getvara ('L_pastfra', iou, imtm2*jmtm2, ib, ic, tmpij
     &,         c1, c0)
              agric(2:imtm1,2:jmtm1,2) = agric(2:imtm1,2:jmtm1,2)
     &                                 + tmpij(1:imtm2,1:jmtm2)
            endif
# endif
            call embmbc (agric(1,1,2))
            agric(:,:,1) = agric(:,:,2)
            agric(:,:,3) = agric(:,:,2)
          endif
        endif

        if (intrp) then
!         interpolate data
          wt1 = 1.
          if (yrv(3) .ne. yrv(1)) wt1 = (yrv(3)-yrv(2))/(yrv(3)-yrv(1))
          wt1 = max(0., min(1., wt1))
          wt3 = 1. - wt1
          do j=1,jmt
            do i=1,imt
              agric(i,j,2) = agric(i,j,1)*wt1 + agric(i,j,3)*wt3
            enddo
          enddo
          call embmbc (agric(1,1,2))
        endif
      else
        first_time = .true.
        intrp = .false.

      endif
# else
      first_time = .true.
      intrp = .false.
# endif

!-----------------------------------------------------------------------
!     calculate surface coalbedo and Dalton number for land
!-----------------------------------------------------------------------

      if (intrp .or. first_time) then
        do j=1,jmt
          do i=1,imt
# if defined O_mtlm
            if (iveg(i,j) .gt. 0 .and. land_map(i,j) .eq. 0 .and.
     &          tmsk(i,j) .lt. 0.5) then
# else
            if (iveg(i,j) .gt. 0 .and. tmsk(i,j) .lt. 0.5) then
# endif
# if defined O_crop_data ||  defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
              sbc(i,j,isca) = 1. - (agric(i,j,2)*veg_alb(iagric)
     &                      + veg_alb(iveg(i,j))*(1. - agric(i,j,2)))
              z0 = veg_rl(iveg(i,j))*(1.-agric(i,j,2))
     &           + veg_rl(iagric)*agric(i,j,2)
!             Dalt = (k**2)/(log(z/z0)*log(z/z0q))
!             where k=0.4, z=10 meters and z0q=z0*e**(-2)
              veg_dalt(i,j) = 0.16/(log(10./z0)*log(73.89/z0))
            endif
          enddo
        enddo
# else
              sbc(i,j,isca) = 1.- veg_alb(iveg(i,j))
            endif
          enddo
        enddo
        do n=1,nveg
          z0 = veg_rl(n)
!         Dalt = (k**2)/(log(z/z0)*log(z/z0q))
!         where k=0.4, z=10 meters and z0q=z0*e**(-2)
          veg_dalt(n) = 0.16/(log(10./z0)*log(73.89/z0))
        enddo
# endif
      endif
#endif

      return
      end
