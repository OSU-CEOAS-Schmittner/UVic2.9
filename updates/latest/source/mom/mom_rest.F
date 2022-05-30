      subroutine mom_rest_in (fname, ids, ide, jds, jde)

#if defined O_mom
!=======================================================================
!     input routine for ocean restarts

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

      integer i, iou, j, ln, n, ntrec, ids, ide, jds, jde, ils
      integer ile, jls, jle, kls, kle, ib(10), ic(10), jrow, k
      integer nyear, nmonth, nday, nhour, nmin, nsec, salt_on_off

      logical inqvardef, exists, warning

      real tmp
      real, allocatable :: tmpij(:,:), tmpik(:,:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "mw.h"
      include "tmngr.h"
      include "switch.h"
#if defined O_mobi_silicon
      include "atm.h"
      include "mobi.h"
#endif

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
      kls = 1
      kle = km

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpik(ils:ile,kls:kle) )

!-----------------------------------------------------------------------
!     read 1d data (t)
!-----------------------------------------------------------------------
      tmp = nots
      call getvars ('nots', iou, ntrec, tmp, c1, c0)
      nots = tmp
      tmp = dayoyr
      tmp = itt
      call getvars ('itt', iou, 1, tmp, c1, c0)
      itt = tmp
      tmp = irstdy
      call getvars ('irstdy', iou, 1, tmp, c1, c0)
      irstdy = tmp
      tmp = msrsdy
      call getvars ('msrsdy', iou, 1, tmp, c1, c0)
      msrsdy = tmp
      tmp = relyr
      call getvars ('relyr', iou, 1, tmp, c1, c0)
      relyr = tmp
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
      if (init_time_in) then
        itt = 0
        irstdy = 0
        msrsdy = 0
        relyr = 0.
        call mkstmp (stamp, year0, month0, day0, hour0, min0, sec0)
      endif

#if defined O_mobi_silicon
      tmp = globalsilwflx
      call getvars ('globalsilwflx', iou, 1, tmp, c1, c0)
      globalsilwflx = tmp
#endif

!-----------------------------------------------------------------------
!     read 2d data (x,y)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = jle-jls+1
      ln = ic(1)*ic(2)

!     check if kmt is consistent
      tmpij(:,:) = kmt(ils:ile,jls:jle)
      call getvara ('kmt', iou, ln, ib, ic, tmpij, c1, c0)
      warning = .false.
      do j=jls,jle
        do i=ils,ile
          if (kmt(i,j) .eq. 0 .and. int(tmpij(i,j)) .ne. 0) then
            print*, ' '
            print*, '==> Error: restart kmt is different. ocean point'
            print*, '           has been changed to land. this usually'
            print*, '           causes problems in the ocean solver.'
            print*, ' '
            stop '=>mom_rest_in'
          endif
          if (kmt(i,j) .ne. int(tmpij(i,j))) warning = .true.
        enddo
      enddo
      if (warning) then
        print*, ' '
        print*, '==> Warning: restart kmt is different. this often '
        print*, '             causes problems in the ocean solver.'
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
      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk

!     update pointers to tau-1, tau, & tau+1 data in the MW based on itt
      if (wide_open_mw) then
!       rotate time levels instead of moving data
        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif

# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      tmpij(ils:ile,jls:jle) = ps(ils:ile,jls:jle,1)
      call getvara ('ps1', iou, ln, ib, ic, tmpij, c1, c0)
      ps(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = ps(ils:ile,jls:jle,2)
      call getvara ('ps2', iou, ln, ib, ic, tmpij, c1, c0)
      ps(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
!     (ubarm, ubar) is at (tau, tau+1)
      tmpij(ils:ile,jls:jle) = ubarm1(ils:ile,jls:jle,1)
      call getvara ('ubarm1', iou, ln, ib, ic, tmpij, c1, c0)
      ubarm1(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = ubarm1(ils:ile,jls:jle,2)
      call getvara ('ubarm2', iou, ln, ib, ic, tmpij, c1, c0)
      ubarm1(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = ubar(ils:ile,jls:jle,1)
      call getvara ('ubar1', iou, ln, ib, ic, tmpij, c1, c0)
      ubar(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = ubar(ils:ile,jls:jle,2)
      call getvara ('ubar2', iou, ln, ib, ic, tmpij, c1, c0)
      ubar(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)

!     initialize the kflds disk which is used for surface pressure.
      pguess(:,:) = 0.
      ptd(:,:) = ps(:,:,1)
      call oput (kflds, nwds, nkflds-1, ptd)
      call oput (kflds, nwds, nkflds, ptd)

# endif
# if defined O_stream_function
!     first do psi at "tau" then at "tau+1"
      tmpij(ils:ile,jls:jle) = psi(ils:ile,jls:jle,1)
      call getvara ('psi1', iou, ln, ib, ic, tmpij, c1, c0)
      psi(ils:ile,jls:jle,1) = tmpij(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = psi(ils:ile,jls:jle,2)
      call getvara ('psi2', iou, ln, ib, ic, tmpij, c1, c0)
      psi(ils:ile,jls:jle,2) = tmpij(ils:ile,jls:jle)
!     guess fields
      tmpij(ils:ile,jls:jle) = ptd(ils:ile,jls:jle)
      call getvara ('ptd1', iou, ln, ib, ic, tmpij, c1, c0)
      ptd(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call oput (kflds, nwds, nkflds-1, ptd)
      tmpij(ils:ile,jls:jle) = ptd(ils:ile,jls:jle)
      call getvara ('ptd2', iou, ln, ib, ic, tmpij, c1, c0)
      ptd(ils:ile,jls:jle) = tmpij(ils:ile,jls:jle)
      call oput (kflds, nwds, nkflds, ptd)
# endif

# if defined O_mom
!     construct depth arrays associated with "u" cells
      call depth_u (kmt, imt, jmt, zw, km, kmu, h, hr)

# endif
!-----------------------------------------------------------------------
!     read 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = 1
      ib(3) = 1
      ic(3) = kle-kls+1
      ib(4) = 1
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)

!     read the "tau" latitude rows
      do jrow=1,jmt
        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
          call getrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
        endif
        ib(2) = jrow
        do n=1,nt
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          var1 = 'tracer1_'//trim(a3)
          if (trim(mapt(n)) .eq. 'temp') then
            if (inqvardef('temp1', iou)) var1 = 'temp1'
          elseif (trim(mapt(n)) .eq. 'salt') then
            if (inqvardef('salt1', iou)) var1 = 'salt1'
          elseif (trim(mapt(n)) .eq. 'dic') then
            if (inqvardef('dic1', iou)) var1 = 'dic1'
          elseif (trim(mapt(n)) .eq. 'dic13') then
            if (inqvardef('dic131', iou)) var1 = 'dic131'
          elseif (trim(mapt(n)) .eq. 'alk') then
            if (inqvardef('alk1', iou)) var1 = 'alk1'
          elseif (trim(mapt(n)) .eq. 'sil') then
            if (inqvardef('sil1', iou)) var1 = 'sil1'
          elseif (trim(mapt(n)) .eq. 'o2') then
            if (inqvardef('o21', iou)) var1 = 'o21'
          elseif (trim(mapt(n)) .eq. 'po4') then
            if (inqvardef('po41', iou)) var1 = 'po41'
          elseif (trim(mapt(n)) .eq. 'dfe') then
            if (inqvardef('dfe1', iou)) var1 = 'dfe1'
          elseif (trim(mapt(n)) .eq. 'dop') then
            if (inqvardef('dop1', iou)) var1 = 'dop1'
          elseif (trim(mapt(n)) .eq. 'phyt') then
            if (inqvardef('phyt1', iou)) var1 = 'phyt1'
          elseif (trim(mapt(n)) .eq. 'phyt_phos') then
            if (inqvardef('phyt_phos1', iou)) var1 = 'phyt_phos1'
          elseif (trim(mapt(n)) .eq. 'zoop') then
            if (inqvardef('zoop1', iou)) var1 = 'zoop1'
          elseif (trim(mapt(n)) .eq. 'detr') then
            if (inqvardef('detr1', iou)) var1 = 'detr1'
          elseif (trim(mapt(n)) .eq. 'detr_phos') then
            if (inqvardef('detr_phos1', iou)) var1 = 'detr_phos1'
          elseif (trim(mapt(n)) .eq. 'detrfe') then
            if (inqvardef('detrfe1', iou)) var1 = 'detrfe1'
          elseif (trim(mapt(n)) .eq. 'no3') then
            if (inqvardef('no31', iou)) var1 = 'no31'
          elseif (trim(mapt(n)) .eq. 'don') then
            if (inqvardef('don1', iou)) var1 = 'don1'
          elseif (trim(mapt(n)) .eq. 'diaz') then
            if (inqvardef('diaz1', iou)) var1 = 'diaz1'
          elseif (trim(mapt(n)) .eq. 'din15') then
            if (inqvardef('din151', iou)) var1 = 'din151'
          elseif (trim(mapt(n)) .eq. 'don15') then
            if (inqvardef('don151', iou)) var1 = 'don151'
          elseif (trim(mapt(n)) .eq. 'phytn15') then
            if (inqvardef('phytn151', iou)) var1 = 'phytn151'
          elseif (trim(mapt(n)) .eq. 'zoopn15') then
            if (inqvardef('zoopn151', iou)) var1 = 'zoopn151'
          elseif (trim(mapt(n)) .eq. 'detrn15') then
            if (inqvardef('detrn151', iou)) var1 = 'detrn151' 
          elseif (trim(mapt(n)) .eq. 'diazn15') then
            if (inqvardef('diazn151', iou)) var1 = 'diazn151'
          elseif (trim(mapt(n)) .eq. 'doc13') then
            if (inqvardef('doc131', iou)) var1 = 'doc131'
          elseif (trim(mapt(n)) .eq. 'phytc13') then
            if (inqvardef('phytc131', iou)) var1 = 'phytc131'
          elseif (trim(mapt(n)) .eq. 'zoopc13') then
            if (inqvardef('zoopc131', iou)) var1 = 'zoopc131'
          elseif (trim(mapt(n)) .eq. 'detrc13') then
            if (inqvardef('detrc131', iou)) var1 = 'detrc131'
          elseif (trim(mapt(n)) .eq. 'diazc13') then
            if (inqvardef('diazc131', iou)) var1 = 'diazc131'
          elseif (trim(mapt(n)) .eq. 'c14') then
            if (inqvardef('c141', iou)) var1 = 'c141'
          elseif (trim(mapt(n)) .eq. 'cfc11') then
            if (inqvardef('cfc111', iou)) var1 = 'cfc111'
          elseif (trim(mapt(n)) .eq. 'cfc12') then
            if (inqvardef('cfc121', iou)) var1 = 'cfc121'
          elseif (trim(mapt(n)) .eq. 'diatn15') then
            if (inqvardef('diatn151', iou)) var1 = 'diatn151'
          elseif (trim(mapt(n)) .eq. 'caco3') then
            if (inqvardef('caco31', iou)) var1 = 'caco31'
          elseif (trim(mapt(n)) .eq. 'diatc13') then
            if (inqvardef('diatc131', iou)) var1 = 'diatc131'
          elseif (trim(mapt(n)) .eq. 'caco3c13') then
            if (inqvardef('caco3c131', iou)) var1 = 'caco3c131'
          elseif (trim(mapt(n)) .eq. 'detr_B') then
            if (inqvardef('detr1_B', iou)) var1 = 'detr1_B'
          elseif (trim(mapt(n)) .eq. 'diat') then
            if (inqvardef('diat1', iou)) var1 = 'diat1'
          endif
          exists = inqvardef(trim(var1), iou)
          if (exists) then
            tmpik(ils:ile,kls:kle) = t(ils:ile,kls:kle,j,n,tau)
            call getvara(trim(var1), iou, ln, ib, ic, tmpik, c1, c0)
            t(ils:ile,kls:kle,j,n,tau) = tmpik(ils:ile,kls:kle)
          else
            if (jrow .eq. 1) print*, '==> Warning: ', trim(var1)
     &,       ' not found'
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       The section below is used to alter the salinity for the purpose
!       of simulating changes in sea level rise. To increase by 1 PSU,
!       0.001 is added to the data read from the restart file. In
!       control.in, init = .false. and rest = .true.. This causes the
!       restart file to be read, and a new one to be written. After one
!       very very short run, this salinity addition needs to be turned
!       off or you will be continuously adding salt with every
!       subsequent run. Later in the code taum1 (tau-1.0) is set to tau
!       and this salt_on_off variable will also control the alteration
!       need to be applied to taup1 (tau=1.0).  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          salt_on_off = 0
          if (salt_on_off .eq. 1) then                       
          if (n .eq. 2) then
          t(:,:,j,2,tau) = t(:,:,j,2,tau)+0.001 
          endif
          endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
        exists = inqvardef('u1', iou)
        if (exists) then
          tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,1,tau)
          call getvara('u1', iou, ln, ib, ic, tmpik, c1, c0)
          u(ils:ile,kls:kle,j,1,tau) = tmpik(ils:ile,kls:kle)
        else
          if (jrow .eq. 1) print*, '==> Warning: u1 not found'
        endif
        exists = inqvardef('v1', iou)
        if (exists) then
          tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,2,tau)
          call getvara('v1', iou, ln, ib, ic, tmpik, c1, c0)
          u(ils:ile,kls:kle,j,2,tau) = tmpik(ils:ile,kls:kle)
        else
          if (jrow .eq. 1) print*, '==> Warning: v1 not found'
        endif
!       initialize every latitude
        if (wide_open_mw) then
!         do nothing since "tau" data is in place in the MW
        else
          call putrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
        endif
      enddo
      if (wide_open_mw) then
!       Initialize 1st and last latitude row for tau-1 to prevent
!       use of uninitialized values on boundary row.
        do j=1,jmt,jmt-1
          do k=1,km
            do i=1,imt
              u(i,k,j,1,taum1) = u(i,k,j,1,tau)
              u(i,k,j,2,taum1) = u(i,k,j,2,tau)
              do n=1,nvar-2
                t(i,k,j,n,taum1) = t(i,k,j,n,tau)
              enddo
            enddo
          enddo
        enddo
      endif

!     read the "tau+1" latitude rows
      do jrow=1,jmt
        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
        endif
        ib(2) = jrow
        do n=1,nt
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          var2 = 'tracer2_'//trim(a3)
          if (trim(mapt(n)) .eq. 'temp') then
            if (inqvardef('temp2', iou)) var2 = 'temp2'
          elseif (trim(mapt(n)) .eq. 'salt') then
            if (inqvardef('salt2', iou)) var2 = 'salt2'
          elseif (trim(mapt(n)) .eq. 'dic') then
            if (inqvardef('dic2', iou)) var2 = 'dic2'
          elseif (trim(mapt(n)) .eq. 'dic13') then
            if (inqvardef('dic132', iou)) var2 = 'dic132'
          elseif (trim(mapt(n)) .eq. 'alk') then
            if (inqvardef('alk2', iou)) var2 = 'alk2'
          elseif (trim(mapt(n)) .eq. 'sil') then
            if (inqvardef('sil2', iou)) var2 = 'sil2'
          elseif (trim(mapt(n)) .eq. 'o2') then
            if (inqvardef('o22', iou)) var2 = 'o22'
          elseif (trim(mapt(n)) .eq. 'po4') then
            if (inqvardef('po42', iou)) var2 = 'po42'
          elseif (trim(mapt(n)) .eq. 'dfe') then
            if (inqvardef('dfe2', iou)) var2 = 'dfe2'
          elseif (trim(mapt(n)) .eq. 'dop') then
            if (inqvardef('dop2', iou)) var2 = 'dop2'
          elseif (trim(mapt(n)) .eq. 'phyt') then
            if (inqvardef('phyt2', iou)) var2 = 'phyt2'
          elseif (trim(mapt(n)) .eq. 'phyt_phos') then
            if (inqvardef('phyt_phos2', iou)) var2 = 'phyt_phos2'
          elseif (trim(mapt(n)) .eq. 'zoop') then
            if (inqvardef('zoop2', iou)) var2 = 'zoop2'
          elseif (trim(mapt(n)) .eq. 'detr') then
            if (inqvardef('detr2', iou)) var2 = 'detr2'
          elseif (trim(mapt(n)) .eq. 'detr_phos') then
            if (inqvardef('detr_phos2', iou)) var2 = 'detr_phos2'
          elseif (trim(mapt(n)) .eq. 'detrfe') then
            if (inqvardef('detrfe2', iou)) var2 = 'detrfe2'
          elseif (trim(mapt(n)) .eq. 'no3') then
            if (inqvardef('no32', iou)) var2 = 'no32'
          elseif (trim(mapt(n)) .eq. 'don') then
            if (inqvardef('don2', iou)) var2 = 'don2'
          elseif (trim(mapt(n)) .eq. 'diaz') then
            if (inqvardef('diaz2', iou)) var2 = 'diaz2'
          elseif (trim(mapt(n)) .eq. 'din15') then
            if (inqvardef('din152', iou)) var2 = 'din152'
          elseif (trim(mapt(n)) .eq. 'don15') then
            if (inqvardef('don152', iou)) var2 = 'don152'
          elseif (trim(mapt(n)) .eq. 'phytn15') then
            if (inqvardef('phytn152', iou)) var2 = 'phytn152'
          elseif (trim(mapt(n)) .eq. 'zoopn15') then
            if (inqvardef('zoopn152', iou)) var2 = 'zoopn152'
          elseif (trim(mapt(n)) .eq. 'detrn15') then
            if (inqvardef('detrn152', iou)) var2 = 'detrn152'
          elseif (trim(mapt(n)) .eq. 'diazn15') then
            if (inqvardef('diazn152', iou)) var2 = 'diazn152'
          elseif (trim(mapt(n)) .eq. 'doc13') then
            if (inqvardef('doc132', iou)) var2 = 'doc132'
          elseif (trim(mapt(n)) .eq. 'phytc13') then
            if (inqvardef('phytc132', iou)) var2 = 'phytc132'
          elseif (trim(mapt(n)) .eq. 'zoopc13') then
            if (inqvardef('zoopc132', iou)) var2 = 'zoopc132'
          elseif (trim(mapt(n)) .eq. 'detrc13') then
            if (inqvardef('detrc132', iou)) var2 = 'detrc132'
          elseif (trim(mapt(n)) .eq. 'diazc13') then
            if (inqvardef('diazc132', iou)) var2 = 'diazc132'
          elseif (trim(mapt(n)) .eq. 'c14') then
            if (inqvardef('c142', iou)) var2 = 'c142'
          elseif (trim(mapt(n)) .eq. 'cfc11') then
            if (inqvardef('cfc112', iou)) var2 = 'cfc112'
          elseif (trim(mapt(n)) .eq. 'cfc12') then
            if (inqvardef('cfc122', iou)) var2 = 'cfc122'
          elseif (trim(mapt(n)) .eq. 'diatn15') then
            if (inqvardef('diatn152', iou)) var2 = 'diatn152'
          elseif (trim(mapt(n)) .eq. 'caco3') then
            if (inqvardef('caco32', iou)) var2 = 'caco32'
          elseif (trim(mapt(n)) .eq. 'diatc13') then
            if (inqvardef('diatc132', iou)) var2 = 'diatc132'
          elseif (trim(mapt(n)) .eq. 'caco3c13') then
            if (inqvardef('caco3c132', iou)) var2 = 'caco3c132'
          elseif (trim(mapt(n)) .eq. 'detr_B') then
            if (inqvardef('detr2_B', iou)) var2 = 'detr2_B'
          elseif (trim(mapt(n)) .eq. 'diat') then
            if (inqvardef('diat2', iou)) var2 = 'diat2'
          endif
          exists = inqvardef(trim(var2), iou)
          if (exists) then
            tmpik(ils:ile,kls:kle) = t(ils:ile,kls:kle,j,n,taup1)
            call getvara(trim(var2), iou, ln, ib, ic, tmpik, c1, c0)
            t(ils:ile,kls:kle,j,n,taup1) = tmpik(ils:ile,kls:kle)
          else
            if (jrow .eq. 1) print*, '==> Warning: ',trim(var2)
     &,       ' not found'
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       The section below is used to alter the salinity for the purpose
!       of simulating changes in sea level rise. See associated section
!       above for tau
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         if (salt_on_off .eq. 1) then
          if (n .eq. 2) then
          t(:,:,j,2,taup1) = t(:,:,j,2,taup1)+0.001 
          endif
         endif      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        enddo
        exists = inqvardef('u2', iou)
        if (exists) then
          tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,1,taup1)
          call getvara('u2', iou, ln, ib, ic, tmpik, c1, c0)
          u(ils:ile,kls:kle,j,1,taup1) = tmpik(ils:ile,kls:kle)
        else
          if (jrow .eq. 1) print*, '==> Warning: u2 not found'
        endif
        exists = inqvardef('v2', iou)
        if (exists) then
          tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,2,taup1)
          call getvara('v2', iou, ln, ib, ic, tmpik, c1, c0)
          u(ils:ile,kls:kle,j,2,taup1) = tmpik(ils:ile,kls:kle)
        else
          if (jrow .eq. 1) print*, '==> Warning: v2 not found'
        endif
!       initialize every latitude
        if (wide_open_mw) then
!         do nothing since "tau+1" data is in place in the MW
        else
          call putrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
        endif
      enddo

      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Ocn restart read from ',trim(fname),' on ', nstamp

      deallocate (tmpij)
      deallocate (tmpik)

      return
      end

      subroutine mom_rest_def (fname)
!=======================================================================
!     definition routine for ocean restarts

!   inputs:
!     fname = file name
!=======================================================================

      implicit none

      character(*) :: fname
      character(3) :: a3

      integer iou, n, ntrec, igs, ige, ig, jgs, jge, jg, kgs
      integer kge, kg, it(10), iu(10), id_time, id_xt, id_xu
      integer id_yt, id_yu, id_zt, id_zw, id_xt_e, id_xu_e
      integer id_yt_e, id_yu_e, id_zt_e, id_zw_e

      real c100, c1e3, c1e20, c1e6

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "iounit.h"
      include "csbc.h"
      include "mw.h"
      include "tmngr.h"

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
     &, 'time', 'time', 'years since 0-1-1')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('nots', iou, 1, it, c0, c0, ' ', 'D'
     &, 'nots', ' ',' ')
      call defvar ('itt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'itt', ' ',' ')
      call defvar ('irstdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'irstdy', ' ',' ')
      call defvar ('msrsdy', iou, 1, it, c0, c0, ' ', 'D'
     &, 'msrsdy', ' ',' ')
      call defvar ('relyr', iou, 1, it, c0, c0, ' ', 'D'
     &, 'relyr', ' ',' ')
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
      it(1) = id_zt
      call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      it(1) = id_xu
      call defvar ('longitude_V', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      it(1) = id_yu
      call defvar ('latitude_V', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      it(1) = id_zw
      call defvar ('depth_W', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')
      it(1) = id_zt_e
      call defvar ('depth_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'depth edges', ' ', 'm')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'longitude edges', ' ', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'latitude edges', ' ', 'degrees_north')
      it(1) = id_zw_e
      call defvar ('depth_W_edges', iou, 1, it, c0, c0, ' ', 'D'
     &, 'depth edges', ' ', 'm')

#if defined O_mobi_silicon
      call defvar ('globalsilwflx', iou, 1, it, c0, c0, ' ', 'D'
     &, 'globalsilwflx', ' ',' ')
#endif

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('kmt', iou, 2, it, c0, c1e3, ' ', 'I'
     &,  'kmt', ' ' ,' ')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      call defvar ('ps1', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'ps1', ' ', ' ')
      call defvar ('ps2', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'ps2', ' ', ' ')
      call defvar ('ubarm1', iou, 3, iu,  -c1e20, c1e20, ' ', 'D'
     &, 'ubarm1', ' ', ' ')
      call defvar ('ubarm2', iou, 3, iu,  -c1e20, c1e20, ' ', 'D'
     &, 'ubarm2', ' ', ' ')
      call defvar ('ubar1', iou, 3, iu,  -c1e20, c1e20, ' ', 'D'
     &, 'ubar1', ' ', ' ')
      call defvar ('ubar2', iou, 3, iu,  -c1e20, c1e20, ' ', 'D'
     &, 'ubar2', ' ', ' ')
# endif
# if defined O_stream_function
      call defvar ('psi1', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'psi1', ' ', ' ')
      call defvar ('psi2', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'psi2', ' ', ' ')
      call defvar ('ptd1', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'ptd1', ' ', ' ')
      call defvar ('ptd2', iou, 3, it,  -c1e20, c1e20, ' ', 'D'
     &, 'ptd2', ' ', ' ')
# endif

!-----------------------------------------------------------------------
!     define 4d data (x,y,z,t)
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
        if (trim(mapt(n)) .eq. 'temp') then
          call defvar ('temp1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean potential temperature at tau', ' ', ' ')
          call defvar ('temp2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean potential temperature at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'salt') then
          call defvar ('salt1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean salinity at tau', '  ', ' ')
          call defvar ('salt2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean salinity at tau+1', '  ', ' ')
        elseif (trim(mapt(n)) .eq. 'dic') then
          call defvar ('dic1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean carbon at tau', ' ', ' ')
          call defvar ('dic2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean carbon at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'dic13') then
          call defvar ('dic131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean carbon 13 at tau', ' ', ' ')
          call defvar ('dic132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean carbon 13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'alk') then
          call defvar ('alk1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'alkalinity at tau',' ', ' ')
          call defvar ('alk2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'alkalinity at tau+1',' ', ' ')
        elseif (trim(mapt(n)) .eq. 'sil') then
          call defvar ('sil1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean silicon at tau', ' ', ' ')
          call defvar ('sil2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ocean silicon at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'o2') then
          call defvar ('o21', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'oxygen at tau',' ', ' ')
          call defvar ('o22', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'oxygen at tau+1',' ', ' ')
        elseif (trim(mapt(n)) .eq. 'po4') then
          call defvar ('po41', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phosphate at tau', ' ', ' ')
          call defvar ('po42', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phosphate at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'dfe') then
          call defvar ('dfe1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'iron at tau', ' ', ' ')
          call defvar ('dfe2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'iron at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'dop') then
          call defvar ('dop1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DOP at tau', ' ', ' ')
          call defvar ('dop2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DOP at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'phyt') then
          call defvar ('phyt1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton at tau', ' ', ' ')
          call defvar ('phyt2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'phyt_phos') then
          call defvar ('phyt_phos1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton phosphorous at tau', ' ', ' ')
          call defvar ('phyt_phos2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton phosphorous at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'zoop') then
          call defvar ('zoop1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton at tau', ' ', ' ')
          call defvar ('zoop2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'detr') then
          call defvar ('detr1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus at tau', ' ', ' ')
          call defvar ('detr2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'detr_phos') then
          call defvar ('detr_phos1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus phosphorous at tau', ' ', ' ')
          call defvar ('detr_phos2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus phosphorous at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'detrfe') then
          call defvar ('detrfe1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'particulate iron at tau', ' ', ' ')
          call defvar ('detrfe2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'particulate iron at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'no3') then
          call defvar ('no31', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'nitrate at tau', ' ', ' ')
          call defvar ('no32', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'nitrate at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'don') then
          call defvar ('don1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DON at tau', ' ', ' ')
          call defvar ('don2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DON at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'diaz') then
          call defvar ('diaz1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs at tau', ' ', ' ')
          call defvar ('diaz2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'din15') then
          call defvar ('din151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'nitrate 15 at tau', ' ', ' ')
          call defvar ('din152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'nitrate 15 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'don15') then
          call defvar ('don151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'nitrate 15 at tau', ' ', ' ')
          call defvar ('don152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DON15 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'phytn15') then
          call defvar ('phytn151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton n15 at tau', ' ', ' ')
          call defvar ('phytn152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton n15 at tau+1', ' ', ' ')  
        elseif (trim(mapt(n)) .eq. 'zoopn15') then
          call defvar ('zoopn151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton n15 at tau', ' ', ' ')
          call defvar ('zoopn152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton n15 at tau+1', ' ', ' ')  
        elseif (trim(mapt(n)) .eq. 'detrn15') then
          call defvar ('detrn151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus n15 at tau', ' ', ' ')
          call defvar ('detrn152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus n15 at tau+1', ' ', ' ')  
       elseif (trim(mapt(n)) .eq. 'diazn15') then
          call defvar ('diazn151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs n15 at tau', ' ', ' ')
          call defvar ('diazn152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs n15 at tau+1', ' ', ' ')
       elseif (trim(mapt(n)) .eq. 'doc13') then
          call defvar ('doc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DOC13 at tau', ' ', ' ')
          call defvar ('doc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'DOC13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'phytc13') then
          call defvar ('phytc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton c13 at tau', ' ', ' ')
          call defvar ('phytc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'phytoplankton c13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'zoopc13') then
          call defvar ('zoopc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton c13 at tau', ' ', ' ')
          call defvar ('zoopc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'zooplankton c13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'detrc13') then
          call defvar ('detrc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus c13 at tau', ' ', ' ')
          call defvar ('detrc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'detritus c13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'diazc13') then
          call defvar ('diazc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs c13 at tau', ' ', ' ')
          call defvar ('diazc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diazotrophs c13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'c14') then
          call defvar ('c141', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'carbon 14 at tau', ' ', ' ')
          call defvar ('c142', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'carbon 14 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call defvar ('cfc111', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'CFC11 at tau', ' ', ' ')
          call defvar ('cfc112', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'CFC11 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call defvar ('cfc121', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'CFC12 at tau', ' ', ' ')
          call defvar ('cfc122', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'CFC12 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'diatn15') then
          call defvar ('diatn151', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatom n15 at tau', ' ', ' ')
          call defvar ('diatn152', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatom n15 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'caco3') then
          call defvar ('caco31', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'calcite at tau', ' ', ' ')
          call defvar ('caco32', iou, 4, it, -c1e20, c1e20, ' ','D'
     &,     'calcite at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'diatc13') then
          call defvar ('diatc131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatom c13 at tau', ' ', ' ')
          call defvar ('diatc132', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatom c13 at tau+1', ' ', ' ')
        elseif (trim(mapt(n)) .eq. 'caco3c13') then
          call defvar ('caco3c131', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'calcite c13 at tau', ' ', ' ')
          call defvar ('caco3c132', iou, 4, it, -c1e20, c1e20, ' ','D'
     &,     'calcite c13 at tau+1', ' ', ' ')
# if defined O_kk_ballast
        elseif (trim(mapt(n)) .eq. 'detr_B') then
          call defvar ('detr1_B', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ballasted detritus at tau', ' ', ' ')
          call defvar ('detr2_B', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'ballasted detritus at tau+1', ' ', ' ')
# endif
        elseif (trim(mapt(n)) .eq. 'diat') then
          call defvar ('diat1', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatoms at tau', ' ', ' ')
          call defvar ('diat2', iou, 4, it, -c1e20, c1e20, ' ', 'D'
     &,     'diatoms at tau+1', ' ', ' ')
        else
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call defvar ('tracer1_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,     ' ', 'D', 'tracer '//trim(a3)//' at tau', ' ', ' ')
          call defvar ('tracer2_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,     ' ', 'D', 'tracer '//trim(a3)//' at tau+1', ' '
     &,     ' ')
        endif
      enddo
      call defvar ('u1', iou , 4, iu, -c1e20, c1e20, ' '
     &, 'D', 'u1', ' ', ' ')
      call defvar ('u2', iou , 4, iu, -c1e20, c1e20, ' '
     &, 'D', 'u2', ' ', ' ')
      call defvar ('v1', iou , 4, iu, -c1e20, c1e20, ' '
     &, 'D', 'v1', ' ', ' ')
      call defvar ('v2', iou , 4, iu, -c1e20, c1e20, ' '
     &, 'D', 'v2', ' ', ' ')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mom_rest_out (fname, ids, ide, jds, jde)
!=======================================================================
!     output routine for ocean restarts

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

      integer i, iou, j, jrow, k, ln, n, ntrec, ids, ide, jds, jde
      integer igs, ige, ig, jgs, jge, jg, kgs, kge, kg, ils, ile
      integer jls, jle, kls, kle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      real c100, c1e3, c1e20, diag1, diag0
      real, allocatable :: tmpij(:,:), tmpik(:,:)
      real, allocatable :: tmpi(:), tmpj(:), tmpk(:)
      real, allocatable :: tmpie(:), tmpje(:), tmpke(:)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
      include "tmngr.h"
#if defined O_mobi_silicon
      include "atm.h"
      include "mobi.h"
#endif

      real xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1), yu_e(jmt+1)
      real zt_e(km+1), zw_e(km+1), tmp, bufsl(imt,km,2)
      real ext(imt,2)

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
      allocate ( tmpik(ils:ile,kls:kle) )
      allocate ( tmpi(igs:ige) )
      allocate ( tmpj(jgs:jge) )
      allocate ( tmpk(kgs:kge) )
      allocate ( tmpie(igs:ige+1) )
      allocate ( tmpje(jgs:jge+1) )
      allocate ( tmpke(kgs:kge+1) )

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
        tmp = 0.
        call putvars ('relyr', iou, ntrec, tmp, c1, c0)
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
        tmp = relyr
        call putvars ('relyr', iou, ntrec, tmp, c1, c0)
      endif
      tmp = nots
      call putvars ('nots', iou, ntrec, tmp, c1, c0)
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

      ic(1) = kg
      tmpk(kgs:kge) = zt(kgs:kge)
      call putvara ('depth', iou, kg, ib, ic, tmpk, c100, c0)
      tmpk(kgs:kge) = zw(kgs:kge)
      call putvara ('depth_W', iou, kg, ib, ic, tmpk, c100, c0)

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

      ic(1) = kg + 1
      call edge_maker (1, zt_e, zt, dzt, zw, dzw, km)
      tmpke(kgs:kge+1) = zt_e(kgs:kge+1)
      call putvara ('depth_edges', iou, kg+1, ib, ic, tmpke, c100, c0)
      call edge_maker (2, zw_e, zt, dzt, zw, dzw, km)
      tmpke(kgs:kge+1) = zw_e(kgs:kge+1)
      call putvara ('depth_W_edges', iou, kg+1, ib, ic, tmpke, c100, c0)

#if defined O_mobi_silicon
# if defined O_kk_si_compensating_sources
      tmp = rivsil*globaldisch !ratio times water flux
# else
      tmp = globalsilwflx
!      print*,'KK rest', globalsilwflx
# endif
      call putvars ('globalsilwflx', iou, 1, tmp, c1, c0)
#endif

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
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      tmpij(ils:ile,jls:jle) = ps(ils:ile,jls:jle,1)
      call putvara ('ps1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = ps(ils:ile,jls:jle,2)
      call putvara ('ps2', iou, ln, ib, ic, tmpij, c1, c0)
!     (ubarm, ubar) is at (tau, tau+1)
      tmpij(ils:ile,jls:jle) = ubarm1(ils:ile,jls:jle,1)
      call putvara ('ubarm1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = ubarm1(ils:ile,jls:jle,2)
      call putvara ('ubarm2', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = ubar(ils:ile,jls:jle,1)
      call putvara ('ubar1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = ubar(ils:ile,jls:jle,2)
      call putvara ('ubar2', iou, ln, ib, ic, tmpij, c1, c0)
# endif
# if defined O_stream_function
!     first do psi at "tau" then at "tau+1"
      tmpij(ils:ile,jls:jle) = psi(ils:ile,jls:jle,1)
      call putvara ('psi1', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = psi(ils:ile,jls:jle,2)
      call putvara ('psi2', iou, ln, ib, ic, tmpij, c1, c0)
!     guess fields
      call oget (kflds, nwds, nkflds-1, ptd)
      tmpij(ils:ile,jls:jle) = ptd(ils:ile,jls:jle)
      call putvara ('ptd1', iou, ln, ib, ic, tmpij, c1, c0)
      call oget (kflds, nwds, nkflds, ptd)
      tmpij(ils:ile,jls:jle) =  ptd(ils:ile,jls:jle)
      call putvara ('ptd2', iou, ln, ib, ic, tmpij, c1, c0)
# endif

!-----------------------------------------------------------------------
!     write 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      ib(1) = 1
      ic(1) = ile-ils+1
      ib(2) = 1
      ic(2) = 1
      ib(3) = 1
      ic(3) = kle-kls+1
      ib(4) = 1
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)

!     save the "tau" latitude rows
      do jrow=1,jmt
        if (wide_open_mw) then
          j = jrow
!         remove external mode from "tau". since psi has been updated
!         psi(,,2) is at "tau"
          if (jrow .lt. jmt) then
# if defined O_stream_function
            do i=2,imt-1
              diag1 = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
              diag0 = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
              ext(i,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
              ext(i,2) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
             enddo
            do k=1,km
              do i=2,imt-1
                if (k .le. kmu(i,jrow)) then
                  bufsl(i,k,1) = (u(i,k,j,1,tau) - ext(i,1))
                  bufsl(i,k,2) = (u(i,k,j,2,tau) - ext(i,2))
                else
                  bufsl(i,k,1) = c0
                  bufsl(i,k,2) = c0
                endif
              enddo
            enddo
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
!           since external mode has been updated, ubarm1 is at "tau"
            do k=1,km
              do i=2,imt-1
                if (k .le. kmu(i,jrow)) then
                  bufsl(i,k,1) = (u(i,k,j,1,tau) - ubarm1(i,jrow,1))
                  bufsl(i,k,2) = (u(i,k,j,2,tau) - ubarm1(i,jrow,2))
                else
                  bufsl(i,k,1) = c0
                  bufsl(i,k,2) = c0
                endif
              enddo
            enddo
# endif
            call setbcx (bufsl(1,1,1), imt, km)
            call setbcx (bufsl(1,1,2), imt, km)
          else
            do k=1,km
              do i=1,imt
                bufsl(i,k,1) = c0
                bufsl(i,k,2) = c0
              enddo
            enddo
          endif
        else
          j = jmw
          call getrow (latdisk(taudisk), nslab, jrow, u(1,1,j,1,tau)
     &,                                               t(1,1,j,1,tau))
          do k=1,km
            do i=1,imt
              bufsl(i,k,1) = u(i,k,j,1,tau)
              bufsl(i,k,2) = u(i,k,j,2,tau)
            enddo
          enddo
        endif

        ib(2) = jrow
        do n=1,nt
          if (trim(mapt(n)) .eq. 'temp') then
            var1 = 'temp1'
          elseif (trim(mapt(n)) .eq. 'salt') then
            var1 = 'salt1'
          elseif (trim(mapt(n)) .eq. 'dic') then
            var1 = 'dic1'
          elseif (trim(mapt(n)) .eq. 'dic13') then
            var1 = 'dic131'
          elseif (trim(mapt(n)) .eq. 'alk') then
            var1 = 'alk1'
          elseif (trim(mapt(n)) .eq. 'sil') then
            var1 = 'sil1'
          elseif (trim(mapt(n)) .eq. 'o2') then
            var1 = 'o21'
          elseif (trim(mapt(n)) .eq. 'po4') then
            var1 = 'po41'
          elseif (trim(mapt(n)) .eq. 'dfe') then
            var1 = 'dfe1'
          elseif (trim(mapt(n)) .eq. 'dop') then
            var1 = 'dop1'
          elseif (trim(mapt(n)) .eq. 'phyt') then
            var1 = 'phyt1'
          elseif (trim(mapt(n)) .eq. 'phyt_phos') then
            var1 = 'phyt_phos1'
          elseif (trim(mapt(n)) .eq. 'zoop') then
            var1 = 'zoop1'
          elseif (trim(mapt(n)) .eq. 'detr') then
            var1 = 'detr1'
          elseif (trim(mapt(n)) .eq. 'detr_phos') then
            var1 = 'detr_phos1'
          elseif (trim(mapt(n)) .eq. 'detrfe') then
            var1 = 'detrfe1'
          elseif (trim(mapt(n)) .eq. 'no3') then
            var1 = 'no31'
          elseif (trim(mapt(n)) .eq. 'don') then
            var1 = 'don1'
          elseif (trim(mapt(n)) .eq. 'diaz') then
            var1 = 'diaz1'
          elseif (trim(mapt(n)) .eq. 'din15') then
            var1 = 'din151'
          elseif (trim(mapt(n)) .eq. 'don15') then
            var1 = 'don151' 
          elseif (trim(mapt(n)) .eq. 'phytn15') then
            var1 = 'phytn151'  
          elseif (trim(mapt(n)) .eq. 'zoopn15') then
            var1 = 'zoopn151'
          elseif (trim(mapt(n)) .eq. 'detrn15') then
            var1 = 'detrn151'   
          elseif (trim(mapt(n)) .eq. 'diazn15') then
            var1 = 'diazn151'
          elseif (trim(mapt(n)) .eq. 'doc13') then
            var1 = 'doc131' 
          elseif (trim(mapt(n)) .eq. 'phytc13') then
            var1 = 'phytc131'
          elseif (trim(mapt(n)) .eq. 'zoopc13') then
            var1 = 'zoopc131'
          elseif (trim(mapt(n)) .eq. 'detrc13') then
            var1 = 'detrc131'
          elseif (trim(mapt(n)) .eq. 'diazc13') then
            var1 = 'diazc131'
          elseif (trim(mapt(n)) .eq. 'c14') then
            var1 = 'c141'
          elseif (trim(mapt(n)) .eq. 'cfc11') then
            var1 = 'cfc111'
          elseif (trim(mapt(n)) .eq. 'cfc12') then
            var1 = 'cfc121'
          elseif (trim(mapt(n)) .eq. 'diatn15') then
            var1 = 'diatn151'
          elseif (trim(mapt(n)) .eq. 'caco3') then
            var1 = 'caco31'
          elseif (trim(mapt(n)) .eq. 'diatc13') then
            var1 = 'diatc131'
          elseif (trim(mapt(n)) .eq. 'caco3c13') then
            var1 = 'caco3c131'
# if defined O_kk_ballast
          elseif (trim(mapt(n)) .eq. 'detr_B') then
            var1 = 'detr1_B'
# endif
          elseif (trim(mapt(n)) .eq. 'diat') then
            var1 = 'diat1'
          else
            if (n .lt. 1000) write(a3, '(i3)') n
            if (n .lt. 100) write(a3, '(i2)') n
            if (n .lt. 10) write(a3, '(i1)') n
            var1 = 'tracer1_'//trim(a3)
          endif
          tmpik(ils:ile,kls:kle) = t(ils:ile,kls:kle,j,n,tau)
          call putvara(trim(var1), iou, ln, ib, ic, tmpik, c1, c0)
        enddo
        tmpik(ils:ile,kls:kle) = bufsl(ils:ile,kls:kle,1)
        call putvara('u1', iou, ln, ib, ic, tmpik, c1, c0)
        tmpik(ils:ile,kls:kle) = bufsl(ils:ile,kls:kle,2)
        call putvara('v1', iou, ln, ib, ic, tmpik, c1, c0)
      enddo

!     save the "tau+1" latitude rows
      do jrow=1,jmt
        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
        endif 
        ib(2) = jrow
        do n=1,nt
          if (trim(mapt(n)) .eq. 'temp') then
            var2 = 'temp2'
          elseif (trim(mapt(n)) .eq. 'salt') then
            var2 = 'salt2'
          elseif (trim(mapt(n)) .eq. 'dic') then
            var2 = 'dic2'
          elseif (trim(mapt(n)) .eq. 'dic13') then
            var2 = 'dic132'
          elseif (trim(mapt(n)) .eq. 'alk') then
            var2 = 'alk2'
          elseif (trim(mapt(n)) .eq. 'sil') then
            var2 = 'sil2'
          elseif (trim(mapt(n)) .eq. 'o2') then
            var2 = 'o22'
          elseif (trim(mapt(n)) .eq. 'po4') then
            var2 = 'po42'
          elseif (trim(mapt(n)) .eq. 'dfe') then
            var2 = 'dfe2'
          elseif (trim(mapt(n)) .eq. 'dop') then
            var2 = 'dop2'
          elseif (trim(mapt(n)) .eq. 'phyt') then
            var2 = 'phyt2'
          elseif (trim(mapt(n)) .eq. 'phyt_phos') then
            var2 = 'phyt_phos2'
          elseif (trim(mapt(n)) .eq. 'zoop') then
            var2 = 'zoop2'
          elseif (trim(mapt(n)) .eq. 'detr') then
            var2 = 'detr2'
          elseif (trim(mapt(n)) .eq. 'detr_phos') then
            var2 = 'detr_phos2'
          elseif (trim(mapt(n)) .eq. 'detrfe') then
            var2 = 'detrfe2'
          elseif (trim(mapt(n)) .eq. 'no3') then
            var2 = 'no32'
          elseif (trim(mapt(n)) .eq. 'don') then
            var2 = 'don2'
          elseif (trim(mapt(n)) .eq. 'diaz') then
            var2 = 'diaz2'
          elseif (trim(mapt(n)) .eq. 'din15') then
            var2 = 'din152'
          elseif (trim(mapt(n)) .eq. 'don15') then
            var2 = 'don152' 
          elseif (trim(mapt(n)) .eq. 'phytn15') then
            var2 = 'phytn152'  
          elseif (trim(mapt(n)) .eq. 'zoopn15') then
            var2 = 'zoopn152'
          elseif (trim(mapt(n)) .eq. 'detrn15') then
            var2 = 'detrn152'   
          elseif (trim(mapt(n)) .eq. 'diazn15') then
            var2 = 'diazn152'
          elseif (trim(mapt(n)) .eq. 'doc13') then
            var2 = 'doc132' 
          elseif (trim(mapt(n)) .eq. 'phytc13') then
            var2 = 'phytc132'
          elseif (trim(mapt(n)) .eq. 'zoopc13') then
            var2 = 'zoopc132'
          elseif (trim(mapt(n)) .eq. 'detrc13') then
            var2 = 'detrc132'
          elseif (trim(mapt(n)) .eq. 'diazc13') then
            var2 = 'diazc132'
          elseif (trim(mapt(n)) .eq. 'c14') then
            var2 = 'c142'
          elseif (trim(mapt(n)) .eq. 'cfc11') then
            var2 = 'cfc112'
          elseif (trim(mapt(n)) .eq. 'cfc12') then
            var2 = 'cfc122'
          elseif (trim(mapt(n)) .eq. 'diatn15') then
            var2 = 'diatn152'
          elseif (trim(mapt(n)) .eq. 'caco3') then
            var2 = 'caco32'
          elseif (trim(mapt(n)) .eq. 'diatc13') then
            var2 = 'diatc132'
          elseif (trim(mapt(n)) .eq. 'caco3c13') then
            var2 = 'caco3c132'
# if defined O_kk_ballast
          elseif (trim(mapt(n)) .eq. 'detr_B') then
            var2 = 'detr2_B'
# endif
          elseif (trim(mapt(n)) .eq. 'diat') then
            var2 = 'diat2'
          else
            if (n .lt. 1000) write(a3, '(i3)') n
            if (n .lt. 100) write(a3, '(i2)') n
            if (n .lt. 10) write(a3, '(i1)') n
            var2 = 'tracer2_'//trim(a3)
          endif
          tmpik(ils:ile,kls:kle) = t(ils:ile,kls:kle,j,n,taup1)
          call putvara(trim(var2), iou, ln, ib, ic, tmpik, c1, c0)
        enddo
        tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,1,taup1)
        call putvara('u2', iou, ln, ib, ic, tmpik, c1, c0)
        tmpik(ils:ile,kls:kle) = u(ils:ile,kls:kle,j,2,taup1)
        call putvara('v2', iou, ln, ib, ic, tmpik, c1, c0)
      enddo

      call rdstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
      print*, '=> Ocn restart written to ',trim(fname),' on ', nstamp

      deallocate ( tmpij )
      deallocate ( tmpik )
      deallocate ( tmpi )
      deallocate ( tmpj )
      deallocate ( tmpk )
      deallocate ( tmpie )
      deallocate ( tmpje )
      deallocate ( tmpke )
#endif

      return
      end
