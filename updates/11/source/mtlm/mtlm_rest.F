      subroutine mtlm_rest_in (fname, ids, ide, jds, jde)

#if defined O_mtlm
!=======================================================================
!     input routine for land restarts

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

      include "size.h"
      include "coord.h"
      include "grdvar.h"
      include "csbc.h"
      include "levind.h"
      include "tmngr.h"
      include "switch.h"
      include "mtlm.h"
# if defined O_mtlm_carbon_13
      include "mtlmc13.h"
# endif
# if defined O_mtlm_carbon_14
      include "mtlmc14.h"
# endif

      character(*) :: fname
      character(32) :: nstamp
      character(3) :: a3

      integer i, iou, j, k, ln, n, ntrec, ids, ide, jds, jde, ils, ile
      integer jls, jle, ib(10), ic(10), undef, L, mskl(imt,jmt)
      integer mskt(imt,jmt)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical exists, inqvardef

      real data(imt,jmt), tmp, c0, c1
      real, allocatable :: tmpij(:,:)

      c0 = 0.
      c1 = 1.
      undef = 0.

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)
      ntrec = 1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global read domain)
!-----------------------------------------------------------------------
      ils = max(ids,1)
      ile = min(ide,imt)
      jls = max(jds,1)
      jle = min(jde,jmt)

      allocate ( tmpij(ils:ile,jls:jle) )
      L = 0
!     mskl is the map of current maximum possible land points
      mskl(:,:) = 0
!     mskt is the map of current maximum possible land points in restart
      mskt(:,:) = 0
      do j=2,jmt-1
        do i=2,imt-1
          if (kmt(i,j) .le. klmax) then
            L = L + 1
            mskl(i,j) = L
          endif
          if (kmt(i,j) .le. 0) mskt(i,j) = L
        enddo
      enddo

!-----------------------------------------------------------------------
!     read 1d data (t)
!-----------------------------------------------------------------------
      tmp = undef
      call getvars ('POINTS', iou, ntrec, tmp, c1, c0)
      if (tmp .gt. POINTS .or. tmp .le. 0. .or. tmp .eq. undef) then
        i = tmp
        print*, "==> Warning: Points in restart is inconsistent"
        print*, "             restart: ", i, "model: ",POINTS
        print*, "             retarting land model from initial values"
        return
      endif
      tmp = atlnd
      call getvars ('atlnd', iou, ntrec, tmp, c1, c0)
      atlnd = tmp
      tmp = LAND_COUNTER
      call getvars ('LAND_COUNTER', iou, ntrec, tmp, c1, c0)
      LAND_COUNTER = tmp
      tmp = dayoyr
      call getvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      dayoyr = tmp
      tmp = ntlbc
      call getvars ('ntlbc', iou, ntrec, tmp, c1, c0)
      ntlbc = tmp
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
      call getvars ('year', iou, ntrec, tmp, c1, c0)
      nyear = tmp
      tmp = month0
      call getvars ('month', iou, ntrec, tmp, c1, c0)
      nmonth = tmp
      tmp = day0
      call getvars ('day', iou, ntrec, tmp, c1, c0)
      nday = tmp
      tmp = hour0
      call getvars ('hour', iou, ntrec, tmp, c1, c0)
      nhour = tmp
      tmp = min0
      call getvars ('minute', iou, ntrec, tmp, c1, c0)
      nmin = tmp
      tmp = sec0
      call getvars ('second', iou, ntrec, tmp, c1, c0)
      nsec = tmp
      call mkstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      if (init_time_in) then
        itt = 0
        irstdy = 0
        msrsdy = 0
        relyr = 0.0
        call mkstmp (stamp, year0, month0, day0, hour0, min0, sec0)
      endif

!-----------------------------------------------------------------------
!     read 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)

      tmpij(ils:ile,jls:jle) = mskt(ils:ile,jls:jle)
      call getvara ('mskl', iou, ln, ib, ic, tmpij, c1, c0)
      mskt(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = land_map(ils:ile,jls:jle)
      call getvara ('land_map', iou, ln, ib, ic, tmpij, c1, c0)
      land_map(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)

      call unloadland (POINTS, TSOIL, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('TSOIL', iou, ln, ib, ic, tmpij, c1, c0)
      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, TSOIL, imt, jmt, mskl, data)

      call unloadland (POINTS, LYING_SNOW, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('LYING_SNOW', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, LYING_SNOW, imt, jmt, mskl, data)

      call unloadland (POINTS, TS1, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('TS1', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, TS1, imt, jmt, mskl, data)

      call unloadland (POINTS, CS, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('CS', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, CS, imt, jmt, mskl, data)

      call unloadland (POINTS, CV, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('CV', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, CV, imt, jmt, mskl, data)

      call unloadland (POINTS, VEG_FRAC, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('VEG_FRAC', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, VEG_FRAC, imt, jmt, mskl, data)

      call unloadland (POINTS, FRAC_VS, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('FRAC_VS', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, FRAC_VS, imt, jmt, mskl, data)

      call unloadland (POINTS, M, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('M', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, M, imt, jmt, mskl, data)

      call unloadland (POINTS, MNEG, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('MNEG', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, MNEG, imt, jmt, mskl, data)

      call unloadland (POINTS, LIT_C_T, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('LIT_C_T', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, LIT_C_T, imt, jmt, mskl, data)

      call unloadland (POINTS, FSMC, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('FSMC', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, FSMC, imt, jmt, mskl, data)

      call unloadland (POINTS, RESP_S_DR, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('RESP_S_DR', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, RESP_S_DR, imt, jmt, mskl, data)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isca)
      call getvara ('sbc_sca', iou, ln, ib, ic, tmpij, c1, c0)
!     only use land points so other values can be changed
      do j=jls,jle
        do i=ils,ile
          if (mskl(i,j) .gt. 0) sbc(i,j,isca) = tmpij(i,j)
        enddo
      enddo

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ilwr)
      call getvara ('sbc_lwr', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ilwr) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isens)
      call getvara ('sbc_sens', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isens) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ievap)
      call getvara ('sbc_evap', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,ievap) = tmpij(ils:ile,jls:jle)

# if defined O_mtlm_carbon_13
      call unloadland (POINTS, CS13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      if (inqvardef('CS13', iou)) then
        call getvara ('CS13', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      else
        call getvara ('CS', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                           *rc13std/(1.+rc13std)
        print*,"mtlm_rest: CS13 not found in restart file"
        print*,"           CS13 initialized with CS*r/(1+r)"
      endif
      call loadland (POINTS, CS13, imt, jmt, mskl, data)

      call unloadland (POINTS, CV13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      if (inqvardef('CV13', iou)) then
        call getvara ('CV13', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      else
        call getvara ('CV', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                           *rc13std/(1.+rc13std)
        print*,"mtlm_rest: CV13 not found in restart file"
        print*,"           CV13 initialized with CV*r/(1+r)"
      endif
      call loadland (POINTS, CV13, imt, jmt, mskl, data)

      call unloadland (POINTS, RESP_S_DR13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('RESP_S_DR13', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, RESP_S_DR13, imt, jmt, mskl, data)

# endif

# if defined O_mtlm_carbon_14
      call unloadland (POINTS, CS14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      if (inqvardef('CS14', iou)) then
        call getvara ('CS14', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      else
        call getvara ('CS', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                           *rc14std/(1.+rc14std)
        print*,"mtlm_rest: CS14 not found in restart file"
        print*,"           CS14 initialized with CS*r/(1+r)"
      endif
      call loadland (POINTS, CS14, imt, jmt, mskl, data)

      call unloadland (POINTS, CV14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      if (inqvardef('CV14', iou)) then
        call getvara ('CV14', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      else
        call getvara ('CV', iou, ln, ib, ic, tmpij, c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                           *rc14std/(1.+rc14std)
        print*,"mtlm_rest: CV14 not found in restart file"
        print*,"           CV14 initialized with CV*r/(1+r)"
      endif
      call loadland (POINTS, CV14, imt, jmt, mskl, data)

      call unloadland (POINTS, RESP_S_DR14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call getvara ('RESP_S_DR14', iou, ln, ib, ic, tmpij, c1, c0)
      where (mskt(ils:ile,jls:jle) .ne. 0)
     &  data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call loadland (POINTS, RESP_S_DR14, imt, jmt, mskl, data)

# endif

# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp)
      call getvara ('sbc_npp', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,inpp) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr)
      call getvara ('sbc_sr', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isr) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn)
      call getvara ('sbc_burn', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,iburn) = tmpij(ils:ile,jls:jle)

#  if defined O_mtlm_carbon_13
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp13)
      call getvara ('sbc_npp13', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,inpp13) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr13)
      call getvara ('sbc_sr13', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isr13) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn13)
      call getvara ('sbc_burn13', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,iburn13) = tmpij(ils:ile,jls:jle)
#  endif

#  if defined O_mtlm_carbon_14
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp14)
      call getvara ('sbc_npp14', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,inpp14) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr14)
      call getvara ('sbc_sr14', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,isr14) = tmpij(ils:ile,jls:jle)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn14)
      call getvara ('sbc_burn14', iou, ln, ib, ic, tmpij, c1, c0)
      sbc(ils:ile,jls:jle,iburn14) = tmpij(ils:ile,jls:jle)
#  endif
# endif


      do n=1,npft
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n

        call unloadland (POINTS, TSTAR(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('TSTAR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, TSTAR(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, ALBSNF(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('ALBSNF_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, ALBSNF(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, ALBSNC(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('ALBSNC_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, ALBSNC(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, HT(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('HT_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, HT(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, LAI(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('LAI_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, LAI(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, C_VEG(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('CVEG_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, C_VEG(1,n), imt, jmt, mskl, data)

# if defined O_mtlm_carbon_13
        call unloadland (POINTS, C_VEG13(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        if (inqvardef('CVEG13_'//trim(a3), iou)) then
          call getvara ('CVEG13_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (mskt(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        else
!         if no restart data initialize with standard ratio
          call getvara ('CVEG_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (mskt(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                             *rc13std/(1.+rc13std)
          print*,"mtlm_rest: "
          print*,"CVEG13_"//trim(a3), " not found in restart file"
          print*,"                  initialized with CVEG*r/(1+r)"
        endif
        call loadland (POINTS, C_VEG13(1,n), imt, jmt, mskl, data)
# endif

# if defined O_mtlm_carbon_14
        call unloadland (POINTS, C_VEG14(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        if (inqvardef('CVEG14_'//trim(a3), iou)) then
          call getvara ('CVEG14_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (mskt(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        else
!         if no restart data initialize with standard ratio
          call getvara ('CVEG_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
          where (mskt(ils:ile,jls:jle) .ne. 0)
     &      data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
     &                             *rc14std/(1.+rc14std)
          print*,"mtlm_rest: "
          print*,"CVEG14_"//trim(a3), " not found in restart file"
          print*,"                  initialized with CVEG*r/(1+r)"
        endif
        call loadland (POINTS, C_VEG14(1,n), imt, jmt, mskl, data)
# endif

        call unloadland (POINTS, G_LEAF_PHEN(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('G_LEAF_PHEN_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, G_LEAF_PHEN(1,n), imt, jmt, mskl
     &,   data)

        call unloadland (POINTS, G_LEAF_DR(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('G_LEAF_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, G_LEAF_DR(1,n), imt, jmt, mskl
     &,   data)

        call unloadland (POINTS, NPP_DR(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('NPP_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, NPP_DR(1,n), imt, jmt, mskl, data)

# if defined O_mtlm_carbon_13
        call unloadland (POINTS, NPP_DR13(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('NPP_DR13_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, NPP_DR13(1,n), imt, jmt, mskl, data)
# endif

# if defined O_mtlm_carbon_14
        call unloadland (POINTS, NPP_DR14(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('NPP_DR14_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, NPP_DR14(1,n), imt, jmt, mskl, data)
# endif

        call unloadland (POINTS, RESP_W_DR(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('RESP_W_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, RESP_W_DR(1,n), imt, jmt, mskl
     &,   data)

        call unloadland (POINTS, CATCH(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('CATCH_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, CATCH(1,n), imt, jmt, mskl, data)

        call unloadland (POINTS, Z0(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('Z0_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, Z0(1,n), imt, jmt, mskl, data)

       enddo

      do n=1,ntype
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n

        call unloadland (POINTS, FRAC(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call getvara ('FRAC_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
        where (mskt(ils:ile,jls:jle) .ne. 0)
     &    data(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
        call loadland (POINTS, FRAC(1,n), imt, jmt, mskl, data)

       enddo
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Lnd restart read from ',trim(fname),' on ', nstamp

      deallocate ( tmpij )

      return
      end

      subroutine mtlm_rest_def (fname)
!=======================================================================
!     definition routine for land restarts

!   inputs:
!     fname = file name
!=======================================================================

      implicit none

      include "size.h"
      include "iounit.h"
      include "tmngr.h"

      character(*) :: fname
      character(3) :: a3

      integer i, iou, j, ln, n, ntrec, igs, ige, ig, jgs, jge, jg
      integer it(10), iu(10), id_time, id_xt, id_xu, id_yt, id_yu
      integer id_xt_e, id_xu_e, id_yt_e, id_yu_e

      real c0, c1, c1e5, c1e20

      c0 = 0.
      c1 = 1.
      c1e5 = 1.e5
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
!     write global attributes
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
      call defdim ('longitude_V', iou, ig, id_xu)
      call defdim ('latitude_V', iou, jg, id_yu)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('longitude_V_edges', iou, ig+1, id_xu_e)
      call defdim ('latitude_V_edges', iou, jg+1, id_yu_e)

!-----------------------------------------------------------------------
!       define 1d data (t)
!-----------------------------------------------------------------------
      it(:) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'years since 0-1-1')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('POINTS', iou, 1, it, c0, c0, ' ', 'D'
     &, 'POINTS', ' ',' ')
      call defvar ('atlnd', iou, 1, it, c0, c0, ' ', 'D'
     &, 'atlnd', ' ',' ')
      call defvar ('LAND_COUNTER', iou, 1, it, c0, c0, ' ', 'D'
     &, 'LAND_COUNTER', ' ',' ')
      call defvar ('dayoyr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'dayoyr', ' ',' ')
      call defvar ('ntlbc', iou, 1, it, c0, c0, ' ', 'D'
     &, 'ntlbc', ' ',' ')
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

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yt
      call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_xu
      call defvar ('longitude_V', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yu
      call defvar ('latitude_V', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time

      call defvar ('mskl', iou , 3, it, c0, c1e5, ' ', 'I'
     &, 'mskl', ' ', ' ')

      call defvar ('land_map', iou , 3, it, c0, c1e5, ' ', 'I'
     &, 'land_map', ' ', ' ')

      call defvar ('TSOIL', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'TSOIL', ' ', ' ')

      call defvar ('LYING_SNOW', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'LYING_SNOW', ' ', ' ')

      call defvar ('TS1', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'TS1', ' ', ' ')

      call defvar ('CS', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CS', ' ', ' ')

      call defvar ('CV', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CV', ' ', ' ')

      call defvar ('VEG_FRAC', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'VEG_FRAC', ' ', ' ')

      call defvar ('FRAC_VS', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'FRAC_VS', ' ', ' ')

      call defvar ('M', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'M', ' ', ' ')

      call defvar ('FSMC', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'FSMC', ' ', ' ')

      call defvar ('RESP_S_DR', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'RESP_S_DR', ' ', ' ')

      call defvar ('MNEG', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'MNEG', ' ', ' ')

      call defvar ('LIT_C_T', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'LIT_C_T', ' ', ' ')

      call defvar ('sbc_sca', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_sca', ' ', ' ')

      call defvar ('sbc_lwr', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_lwr', ' ', ' ')

      call defvar ('sbc_sens', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_sens', ' ', ' ')

      call defvar ('sbc_evap', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_evap', ' ', ' ')

#  if defined O_mtlm_carbon_13
      call defvar ('CS13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CS13', ' ', ' ')

      call defvar ('CV13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CV13', ' ', ' ')

      call defvar ('RESP_S_DR13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'RESP_S_DR13', ' ', ' ')

#  endif
#  if defined O_mtlm_carbon_14
      call defvar ('CS14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CS14', ' ', ' ')

      call defvar ('CV14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'CV14', ' ', ' ')

      call defvar ('RESP_S_DR14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'RESP_S_DR14', ' ', ' ')

#  endif
#  if defined O_carbon
      call defvar ('sbc_npp', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_npp', ' ', ' ')

      call defvar ('sbc_sr', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_sr', ' ', ' ')

      call defvar ('sbc_burn', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_burn', ' ', ' ')

#   if defined O_mtlm_carbon_13
      call defvar ('sbc_npp13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_npp13', ' ', ' ')

      call defvar ('sbc_sr13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_sr13', ' ', ' ')

      call defvar ('sbc_burn13', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_burn13', ' ', ' ')
#   endif
#   if defined O_mtlm_carbon_14
      call defvar ('sbc_npp14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_npp14', ' ', ' ')

      call defvar ('sbc_sr14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_sr14', ' ', ' ')

      call defvar ('sbc_burn14', iou , 3, it, -c1e20, c1e20, ' ', 'D'
     &, 'sbc_burn14', ' ', ' ')
#   endif
#  endif
      do n=1,nPFT
        if (n .lt. 1000) write(a3,'(i3)') n
        if (n .lt. 100) write(a3,'(i2)') n
        if (n .lt. 10) write(a3,'(i1)') n

        call defvar ('TSTAR_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'TSTAR_'//trim(a3), ' ', ' ')

        call defvar ('ALBSNF_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'ALBSNF_'//trim(a3), ' ', ' ')

        call defvar ('ALBSNC_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'ALBSNC_'//trim(a3), ' ', ' ')

        call defvar ('HT_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'HT_'//trim(a3), ' ', ' ')

        call defvar ('LAI_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'LAI_'//trim(a3) , ' ', ' ')

        call defvar ('CVEG_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'CVEG_'//trim(a3) , ' ', ' ')

        call defvar ('G_LEAF_PHEN_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'G_LEAF_PHEN_'//trim(a3), ' ', ' ')

        call defvar ('G_LEAF_DR_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'G_LEAF_DR_'//trim(a3), ' ', ' ')

        call defvar ('NPP_DR_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'NPP_DR_'//trim(a3), ' ', ' ')

        call defvar ('RESP_W_DR_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'RESP_W_DR_'//trim(a3), ' ', ' ')

        call defvar ('CATCH_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'CATCH_'//trim(a3), ' ', ' ')

        call defvar ('Z0_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'Z0_'//trim(a3), ' ', ' ')

#  if defined O_mtlm_carbon_13
        call defvar ('CVEG13_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'CVEG13_'//trim(a3) , ' ', ' ')

        call defvar ('NPP_DR13_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'NPP_DR13_'//trim(a3), ' ', ' ')
#  endif

#  if defined O_mtlm_carbon_14
        call defvar ('CVEG14_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'CVEG14_'//trim(a3) , ' ', ' ')

        call defvar ('NPP_DR14_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'NPP_DR14_'//trim(a3), ' ', ' ')
#  endif
      enddo

      do n=1,ntype
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n

        call defvar ('FRAC_'//trim(a3), iou , 3, it, -c1e20
     &,   c1e20, ' ', 'D', 'FRAC_'//trim(a3), ' ', ' ')

      enddo

      call enddef (iou)
#endif

      return
      end

      subroutine mtlm_rest_out (fname, ids, ide, jds, jde)
!=======================================================================
!     output routine for land restarts

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
#if defined O_mtlm

      implicit none

      include "size.h"
      include "coord.h"
      include "grdvar.h"
      include "csbc.h"
      include "levind.h"
      include "tmngr.h"
      include "switch.h"
      include "mtlm.h"
# if defined O_mtlm_carbon_13
      include "mtlmc13.h"
# endif
# if defined O_mtlm_carbon_14
      include "mtlmc14.h"
# endif
      include "iounit.h"

      character(*) :: fname
      character(32) :: nstamp
      character(3) :: a3

      integer i, iou, j, k, ln, n, ntrec, ids, ide, jds, jde, igs, ige
      integer ig, jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10), L
      integer mskl(imt,jmt)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      real xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1), yu_e(jmt+1)
      real data(imt,jmt), tmp, c0, c1
      real, allocatable :: tmpij(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
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
      allocate ( tmpi(ils:ile) )
      allocate ( tmpj(jls:jle) )
      allocate ( tmpie(ils:ile+1) )
      allocate ( tmpje(jls:jle+1) )
      L = 0
      mskl(:,:) = 0
      do j=2,jmt-1
        do i=2,imt-1
          if (kmt(i,j) .le. klmax) then
            L = L + 1
            mskl(i,j) = L
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      tmp = POINTS
      call putvars ('POINTS', iou, ntrec, tmp, c1, c0)
      tmp = LAND_COUNTER
      call putvars ('LAND_COUNTER', iou, ntrec, tmp, c1, c0)
      call putvars ('atlnd', iou, ntrec, atlnd, c1, c0)
      tmp = dayoyr
      call putvars ('dayoyr', iou, ntrec, tmp, c1, c0)
      tmp = ntlbc
      call putvars ('ntlbc', iou, ntrec, tmp, c1, c0)
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
      call putvars ('year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('second', iou, ntrec, tmp, c1, c0)

!-----------------------------------------------------------------------
!     write 1d data (x or y)
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
      call putvara ('longitude_edges', iou, ig+1, ib, ic, tmpie
     &, c1, c0)
      call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
      tmpie(igs:ige+1) = xu_e(igs:ige+1)
      call putvara ('longitude_V_edges', iou, ig+1, ib, ic, tmpie
     &, c1, c0)

      ic(1) = jg + 1
      call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
      call putvara ('latitude_edges', iou, jg+1, ib, ic, tmpje
     &, c1, c0)
      call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
      tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
      call putvara ('latitude_V_edges', iou, jg+1, ib, ic, tmpje
     &, c1, c0)

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

      tmpij(ils:ile,jls:jle) = mskl(ils:ile,jls:jle)
      call putvara ('mskl', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = land_map(ils:ile,jls:jle)
      call putvara ('land_map', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, TSOIL, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('TSOIL', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, LYING_SNOW, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('LYING_SNOW', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, TS1, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('TS1', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, CS, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CS', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, CV, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CV', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, VEG_FRAC, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('VEG_FRAC', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, FRAC_VS, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('FRAC_VS', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, M, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('M', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, FSMC, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('FSMC', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, RESP_S_DR, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('RESP_S_DR', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, MNEG, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('MNEG', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, LIT_C_T, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('LIT_C_T', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isca)
      call putvara ('sbc_sca', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ilwr)
      call putvara ('sbc_lwr', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isens)
      call putvara ('sbc_sens', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,ievap)
      call putvara ('sbc_evap', iou, ln, ib, ic, tmpij, c1, c0)

# if defined O_carbon
      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp)
      call putvara ('sbc_npp', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr)
      call putvara ('sbc_sr', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn)
      call putvara ('sbc_burn', iou, ln, ib, ic, tmpij, c1, c0)
# endif

# if defined O_mtlm_carbon_13
      call unloadland (POINTS, CS13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CS13', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, CV13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CV13', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, RESP_S_DR13, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('RESP_S_DR13', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp13)
      call putvara ('sbc_npp13', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr13)
      call putvara ('sbc_sr13', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn13)
      call putvara ('sbc_burn13', iou, ln, ib, ic, tmpij, c1, c0)
# endif

# if defined O_mtlm_carbon_14
      call unloadland (POINTS, CS14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CS14', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, CV14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('CV14', iou, ln, ib, ic, tmpij, c1, c0)

      call unloadland (POINTS, RESP_S_DR14, imt, jmt, mskl, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvara ('RESP_S_DR14', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,inpp14)
      call putvara ('sbc_npp14', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,isr14)
      call putvara ('sbc_sr14', iou, ln, ib, ic, tmpij, c1, c0)

      tmpij(ils:ile,jls:jle) = sbc(ils:ile,jls:jle,iburn14)
      call putvara ('sbc_burn14', iou, ln, ib, ic, tmpij, c1, c0)
# endif

      do n=1,npft
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n

        call unloadland (POINTS, TSTAR(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('TSTAR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, ALBSNF(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('ALBSNF_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, ALBSNC(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('ALBSNC_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, HT(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('HT_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, LAI(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('LAI_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, C_VEG(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('CVEG_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, G_LEAF_PHEN(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('G_LEAF_PHEN_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, G_LEAF_DR(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('G_LEAF_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, NPP_DR(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('NPP_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, RESP_W_DR(1,n), imt, jmt, mskl
     &,   data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('RESP_W_DR_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, CATCH(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('CATCH_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, Z0(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('Z0_'//trim(a3), iou, ln, ib, ic, tmpij, c1, c0)

# if defined O_mtlm_carbon_13
        call unloadland (POINTS, C_VEG13(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('CVEG13_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, NPP_DR13(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('NPP_DR13_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
# endif

# if defined O_mtlm_carbon_14
        call unloadland (POINTS, C_VEG14(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('CVEG14_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

        call unloadland (POINTS, NPP_DR14(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('NPP_DR14_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)
# endif

      enddo

      do n=1,ntype
        if (n .lt. 1000) write(a3, '(i3)') n
        if (n .lt. 100) write(a3, '(i2)') n
        if (n .lt. 10) write(a3, '(i1)') n

        call unloadland (POINTS, FRAC(1,n), imt, jmt, mskl, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvara ('FRAC_'//trim(a3), iou, ln, ib, ic, tmpij
     &,   c1, c0)

       enddo

      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Lnd restart written to ',trim(fname),' on ', nstamp

      deallocate ( tmpij )
      deallocate ( tmpi )
      deallocate ( tmpj )
      deallocate ( tmpie )
      deallocate ( tmpje )
#endif

      return
      end
