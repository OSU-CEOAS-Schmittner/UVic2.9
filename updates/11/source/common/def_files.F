      subroutine def_tavg
!=======================================================================
!     defines tavg files for UVic_ESCM
!=======================================================================

      implicit none

      character(120) :: fname

#if defined O_embm
      call def_tavg_embm (fname)
#endif
#if defined O_ism
      call def_tavg_ism (fname)
#endif
#if defined O_mtlm
      call def_tavg_mtlm (fname)
#endif
#if defined O_mom
      call def_tavg_mom (fname)
#endif
#if defined O_sed
      call def_tavg_sed (fname)
#endif

      return
      end

      subroutine def_rest (last)
!=======================================================================
!     defines rest files for UVic_ESCM

!     input:
!       last = last restart flag (1 = last)
!=======================================================================

      implicit none

      character(120) :: fname
      integer last

#if defined O_embm
      call def_rest_embm (last, fname)
#endif
#if defined O_ism
      call def_rest_ism (last, fname)
#endif
#if defined O_mtlm
      call def_rest_mtlm (last, fname)
#endif
#if defined O_mom
      call def_rest_mom (last, fname)
#endif
#if defined O_sed
      call def_rest_sed (last, fname)
#endif

      return
      end

      subroutine def_tsi
!=======================================================================
!     defines tsi files for UVic_ESCM
!=======================================================================

      implicit none

      character(120) :: fname

#if defined O_embm
      call def_tsi_embm (fname)
#endif
#if defined O_ism
      call def_tsi_ism (fname)
#endif
#if defined O_mtlm
      call def_tsi_mtlm (fname)
#endif
#if defined O_mom
      call def_tsi_mom (fname)
#endif
#if defined O_sed
      call def_tsi_sed (fname)
#endif
#if defined O_global_sums || defined O_co2emit_diag
      call def_tsi_gsums (fname)
#endif

      return
      end

      subroutine inqdefined (name, defined)

!=======================================================================
!     keeps track of which files have been defined. assumes the file is
!     about to be defined.

!     input:
!       name = file name

!     output:
!       defined = logical flag
!=======================================================================

      integer n, max_num_files
      parameter (max_num_files=201)

      character(120) :: name
      character(120), allocatable :: file_names(:)

      logical defined

      save file_names

      if (.not. allocated (file_names)) then
        allocate ( file_names(max_num_files) )
        file_names(1:max_num_files) = " "
      endif

      n = 1
      defined = .false.
      do while (.not. defined .and. file_names(n) .ne. " ")
        if (trim(name) .eq. trim(file_names(n))) defined = .true.
        n = n + 1
      enddo
!     always leave file_names(max_num_files) = " " to avoid more
!     testing in the "do while" loop
      if (.not. defined) then
        if (n .lt. max_num_files) then
          file_names(n) = name
        else
          stop 'maximum number of files exceeded in inqdefined'
        endif
      endif

      return
      end

#if defined O_embm
      subroutine def_tavg_embm (fname)
!=======================================================================
!     defines tavg files for the embm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "atm.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"
      include "iounit.h"
      include "tmngr.h"

      defined = .false.
      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_embm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call embm_tavg_def (name, imt, jmt, nat, ncat, xt, yt
     &,                     calendar, expnam, runstamp, mapat)
      endif
      fname = name

      return
      end

      subroutine def_rest_embm (last, fname)
!=======================================================================
!     defines restart file for the embm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_embm.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_embm',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call embm_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_embm (fname)
!=======================================================================
!     defines tsi files for the embm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_embm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call embm_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined O_ism
      subroutine def_tavg_ism (fname)
!=======================================================================
!     defines tavg files for the ism

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_ism',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
!        call ism_tavg_def (name, imt, jmt, xt, yt, calendar, expnam
!     &,                    runstamp)
      endif
      fname = name

      return
      end

      subroutine def_rest_ism (last, fname)
!=======================================================================
!     defines restart file for the ism

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_ism.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_ism',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call ism_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_ism (fname)
!=======================================================================
!     defines tsi files for the ism

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_ism',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
!        call ism_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined O_mtlm
      subroutine def_tavg_mtlm (fname)
!=======================================================================
!     defines tavg files for the mtlm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_mtlm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mtlm_tavg_def (name, imt, jmt, NPFT, NTYPE, xt, yt
     &,                     calendar, expnam, runstamp)
      endif
      fname = name

      return
      end

      subroutine def_rest_mtlm (last, fname)
!=======================================================================
!     defines restart file for the mtlm

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_mtlm.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_mtlm',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call mtlm_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_mtlm (fname)
!=======================================================================
!     defines tsi files for the mtlm

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_mtlm',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mtlm_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined O_mom
      subroutine def_tavg_mom (fname)
!=======================================================================
!     defines tavg files for the mom

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "iounit.h"
      include "mw.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_mom',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mom_tavg_def (name, imt, jmt, km, nt, kpzd, xt, yt
     &,                    calendar, expnam, runstamp, mapt)
      endif
      fname = name

      return
      end

      subroutine def_rest_mom (last, fname)
!=======================================================================
!     defines restart file for the mom

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

       implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_mom.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_mom',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call mom_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_mom (fname)
!=======================================================================
!     defines tsi files for the mom

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_mom',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call mom_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif

#if defined O_sed
      subroutine def_tavg_sed (fname)
!=======================================================================
!     defines tavg files for sed

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "size.h"
      include "coord.h"
      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tavg_sed',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call sed_tavg_def (name, imt, jmt, xt, yt, calendar, expnam
     &,                    runstamp)
      endif
      fname = name

      return
      end

      subroutine def_rest_sed (last, fname)
!=======================================================================
!     defines restart file for the sed

!     input:
!       last = last rest flag (1 = last)

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, name, file_stamp, new_file_name
      character(32) :: nstamp

      integer iou, last, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "tmngr.h"

      if (last .eq. 1) then
        name = 'restart_sed.nc'
      else
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('rest_sed',nstamp,'.nc')
      endif
      name = new_file_name (name)
      call inqdefined (name, defined)
      if (defined) then
        call opennext (name, relyr, ntrec, iou)
      else
        call opennew (name, iou)
      endif
      call sed_rest_def (name)
      fname = name

      return
      end

      subroutine def_tsi_sed (fname)
!=======================================================================
!     defines tsi files for sed

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_sed',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call sed_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif
#if defined O_global_sums || defined O_co2emit_diag

      subroutine def_tsi_gsums (fname)
!=======================================================================
!     defines tsi files for gsums

!     output:
!       fname = file name
!=======================================================================

      implicit none

      character(120) :: fname, file_stamp, name, new_file_name
      character(32) :: nstamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      logical defined

      save name
      data name /' '/

      include "iounit.h"
      include "tmngr.h"

      if (name .eq. ' ') then
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        name = file_stamp ('tsi_gsums',nstamp,'.nc')
        name = new_file_name (name)
        call inqdefined (name, defined)
        if (defined) then
          call opennext (name, relyr, ntrec, iou)
        else
          call opennew (name, iou)
        endif
        call gsums_tsi_def (name, calendar, expnam, runstamp)
      endif
      fname = name

      return
      end
#endif
