      subroutine getrow (lu, nwrs, nrec, u, t)

#if defined O_mom
!=======================================================================
!     get velocity and tracer data from ramdisk latitude row "nrec"

!         lu    = i/o unit number associated with "fname"
!         nwrs  = number of words to read
!         nrec  = latitude row
!         u     = velocities
!         t     = tracers
!         ntrac = number of tracers
!=======================================================================

      implicit none

      integer ntrac, lu, nrec, ntl, nwrs

      real sdisk

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "iounit.h"

      parameter (ntrac=nvar-2)

      real u(imt,km,jmw,2), t(imt,km,jmw,ntrac)
      common /distrib/ sdisk(nslab,jmt,2)

      if (lu .ne. latdisk(1) .and. lu .ne. latdisk(2)) then
        write (stdout,*) '=>Error: unit ',lu,' is not allowed. nrec='
     &, nrec
        stop '=>getrow'
      endif

      if (lu .eq. latdisk(1)) then
        ntl = 1
      elseif (lu .eq. latdisk(2)) then
        ntl = 2
      endif

      call getlat (sdisk(1,nrec,ntl), u, t)

      return

      entry putrow (lu, nwrs, nrec, u, t)

!=======================================================================
!     put velocity and tracer data to ramdisk latitude row "nrec"

!         lu    = i/o unit number associated with "fname"
!         nwrs  = number of words to read
!         nrec  = latitude row
!         u     = velocities
!         t     = tracers
!         ntrac = number of tracers
!=======================================================================

      if (lu .ne. latdisk(1) .and. lu .ne. latdisk(2)) then
        write (stdout,*) '=>Error: unit ',lu,' is not allowed. nrec='
     &, nrec
        stop '=>putrow'
      endif

      if (lu .eq. latdisk(1)) then
        ntl = 1
      elseif (lu .eq. latdisk(2)) then
        ntl = 2
      endif

      call putlat (sdisk(1,nrec,ntl), u, t)
      return
      end

      subroutine getlat (disk, u, t)

!-----------------------------------------------------------------------
!     copy contiguous portions of virtual disk to non-contiguous
!     portions of memory.
!-----------------------------------------------------------------------

      implicit none

      integer ntrac, n, k, i

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"

      parameter (ntrac=nvar-2)

      real u(imt,km,jmw,2), t(imt,km,jmw,ntrac), disk(imt,km,nvar)

      do n=1,2
        do k=1,km
          do i=1,imt
            u(i,k,1,n) = disk(i,k,n)
          enddo
        enddo
      enddo
      do n=1,ntrac
        do k=1,km
          do i=1,imt
            t(i,k,1,n) = disk(i,k,n+2)
          enddo
        enddo
      enddo
      return
      end

      subroutine putlat (disk, u, t)

!-----------------------------------------------------------------------
!     copy non-contiguous portions of memory to contiguous portions
!     of virtual disk.
!-----------------------------------------------------------------------

      implicit none

      integer ntrac, n, k, i

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"

      parameter (ntrac=nvar-2)

      real u(imt,km,jmw,2), t(imt,km,jmw,ntrac), disk(imt,km,nvar)

      do n=1,2
        do k=1,km
          do i=1,imt
            disk(i,k,n) = u(i,k,1,n)
          enddo
        enddo
      enddo
      do n=1,ntrac
        do k=1,km
          do i=1,imt
            disk(i,k,n+2) = t(i,k,1,n)
          enddo
        enddo
      enddo
      return
      end

      subroutine oget (lu, nwrs, nrec, a)

!=======================================================================
!     get 2D field data from ramdisk record "nrec"

!         lu    = i/o unit number associated with "fname"
!         nwrs  = number of words to read
!         nrec  = record
!         a     = two dimensional field
!=======================================================================

      implicit none

      integer n12, lu, nrec, ns, nwrs, ne, n

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "iounit.h"

      parameter (n12=nwds*nkflds)

      real sdisk2d(n12), a(nwds)
      save sdisk2d

      if (lu .ne. kflds) then
        write (stdout,*) '=>Error: unit # ',lu,' not allowed. rec=',nrec
        stop '=>oget'
      else
        ns = (nrec-1)*nwrs + 1
        ne = ns + nwrs - 1
        do n=ns,ne
          a(n-ns+1) = sdisk2d(n)
        enddo
      endif
      return

      entry oput (lu, nwrs, nrec, a)

!=======================================================================
!     put 2D field data to ramdisk record "nrec"

!         lu    = i/o unit number associated with "fname"
!         nwrs  = number of words to read
!         nrec  = record
!         a     = two dimensional field
!=======================================================================

      if (lu .ne. kflds) then
        write (stdout,*) '=>Error: unit # ',lu,' not allowed. rec=',nrec
        stop '=>oget'
      else
        ns = (nrec-1)*nwrs + 1
        ne = ns + nwrs - 1
        do n=ns,ne
          sdisk2d(n) = a(n-ns+1)
        enddo
      endif
      return
      end

      subroutine copy_all_rows (ifrom, ito)

!-----------------------------------------------------------------------
!     copy all latitude rows from time level "ifrom" to time level "ito"
!     for prognostic variables in the MW
!-----------------------------------------------------------------------

      implicit none

      integer n, j, k, i, ito, ifrom

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      do n=1,2
        do j=1,jmw
          do k=1,km
            do i=1,imt
              u(i,k,j,n,ito) = u(i,k,j,n,ifrom)
            enddo
          enddo
        enddo
      enddo
      do n=1,nvar-2
        do j=1,jmw
          do k=1,km
            do i=1,imt
              t(i,k,j,n,ito) = t(i,k,j,n,ifrom)
            enddo
          enddo
        enddo
      enddo
      return
      end

      subroutine euler_shuffle

!-----------------------------------------------------------------------
!     after the second pass of an euler backward step, exchange "tau"
!     and "tau+1" data, after shuffling, data will be in proper position
!     for the next time step.
!-----------------------------------------------------------------------

      implicit none

      integer n, j, k, i

      real temp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      do n=1,2
        do j=1,jmw
          do k=1,km
            do i=1,imt
              temp             = u(i,k,j,n,tau)
              u(i,k,j,n,tau)   = u(i,k,j,n,taup1)
              u(i,k,j,n,taup1) = temp
            enddo
          enddo
        enddo
      enddo
      do n=1,nt
        do j=1,jmw
          do k=1,km
            do i=1,imt
              temp             = t(i,k,j,n,tau)
              t(i,k,j,n,tau)   = t(i,k,j,n,taup1)
              t(i,k,j,n,taup1) = temp
            enddo
          enddo
        enddo
      enddo
#endif

      return
      end
