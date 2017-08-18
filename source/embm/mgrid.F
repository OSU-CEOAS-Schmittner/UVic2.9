      subroutine mgrid (phi, ap, an, as, ae, aw, b, ib, ie, jb, je,
     &                  id, jd, sweepin, levelin, epsin, sweepout,
     &                  levelout, epsout)

#if defined O_embm && defined O_embm_mgrid && defined O_cyclic
!=======================================================================
!                    uvic multigrid solver

!     implementation of a structured multigrid solver.

!     this file contains 11 routines, mgrid, ptr, sumcf, sumrsd,
!     sumdel, pgs, ctdma, tdma, resid, new_ctdma and new_tdma,
!=======================================================================
!     subroutine to update the interior phi field by making ntimes
!     sweeps of the v cycle additive correction multigrid algorithm.
!   input:
!     phi: initial guess for phi
!     ap,aw,ae,as,an: active coefficients for p,w,e,s,n nodes
!     b: accumulated fixed source term
!     sweepin: maximum number of mgrid sweeps
!     levelin: desired number of mgrid levels (>= levelout)
!     epsin: desired epsilon
!     ib,ie,jb,je: first and last interior indices in i and j
!     id,jd:  array dimensions
!   output:
!     phi: updated estimate of phi
!     sweepout: actual number of mgrid sweeps
!     levelout: actual number of mgrid levels
!     epsout: actual epsilon
!=======================================================================

      implicit none

      integer ib, ie, jb, je, id, idl, jd, jdl, ld, l, lstart
      integer sweepin, sweepout, levelin, levelout

      real(kind=8) ap(id,jd), aw(id,jd), ae(id,jd), as(id,jd), an(id,jd)
      real(kind=8) b(id,jd), phi(id,jd), rsd(id,jd), avrsd, avrsdl
      real(kind=8) epsin, epsout, l2rsd0, l2rsd, l2rsdl

!     setup the work arrays and pointers
      integer ibl(levelin), iel(levelin), jbl(levelin), jel(levelin)

      real(kind=8), allocatable, dimension(:,:) :: apl, awl, ael, asl
      real(kind=8), allocatable, dimension(:,:) :: anl, bl, rsdl, phil

!     generate indices of all grid levels
      levelout = levelin
      call ptr (ibl, iel, jbl, jel, levelout, ib, ie, jb, je, levelin)

      jdl = jd + levelout
      idl = (id/2)*2 + 1

!     allocate work arrays
      allocate (apl(idl,jdl))
      allocate (awl(idl,jdl))
      allocate (ael(idl,jdl))
      allocate (asl(idl,jdl))
      allocate (anl(idl,jdl))
      allocate (bl(idl,jdl))
      allocate (rsdl(idl,jdl))
      allocate (phil(idl,jdl))

      apl(:,:) = 0.
      awl(:,:) = 0.
      ael(:,:) = 0.
      asl(:,:) = 0.
      anl(:,:) = 0.
      bl(:,:) = 0.

!     get initial residuals
      call resid (rsd,l2rsd0,phi,ap,aw,ae,as,an,b,ib,ie,jb,je,id,jd)

!     calculate the coefficients of the coarse grid equations
      do l=2,levelout
        if (l .eq. 2) then
          call sumcf (apl,awl,ael,asl,anl,ap,aw,ae,as,an,ibl(l),iel(l),
     &                jbl(l),jel(l),ib,ie,jb,je,idl,jdl,id,jd)
        else
          call sumcf (apl,awl,ael,asl,anl,apl,awl,ael,asl,anl,ibl(l),
     &                iel(l),jbl(l),jel(l),ibl(l-1),iel(l-1),jbl(l-1),
     &                jel(l-1),idl,jdl,idl,jdl)
        endif
      enddo

      do sweepout=1,sweepin

        bl(:,:) = 0.0
        rsdl(:,:) = 0.0
        phil(:,:) = 0.0

!       sweep from fine grid to coarse grid
        call pgs (phi, ap, aw, ae, as, an, b, ib, ie, jb, je, id, jd)
        call resid (rsd,l2rsd,phi,ap,aw,ae,as,an,b,ib,ie,jb,je,id,jd)
        do l=2,levelout
          if (l .eq. 2) then
            call sumrsd (bl,rsd,ibl(l),iel(l),jbl(l),jel(l),ib,ie,jb,
     &                   je,idl,jdl,id,jd)
          else
            call sumrsd (bl,rsdl,ibl(l),iel(l),jbl(l),jel(l),ibl(l-1),
     &                   iel(l-1),jbl(l-1),jel(l-1),idl,jdl,idl,jdl)
          endif
          call pgs (phil,apl,awl,ael,asl,anl,bl,ibl(l),iel(l),jbl(l),
     &              jel(l),idl,jdl)
          call resid (rsdl,l2rsdl,phil,apl,awl,ael,asl,anl,bl,ibl(l),
     &                iel(l),jbl(l),jel(l),idl,jdl)
        enddo

!       sweep from coarse grid to finest grid
        lstart = max0 (1, levelout-1)
        do l=lstart,1,-1
          if (l .eq. 1) then
            call sumdel (phi,phil,ib,ie,jb,je,ibl(l+1),iel(l+1),
     &                   jbl(l+1),jel(l+1),id,jd,idl,jdl)
            call pgs (phi,ap,aw,ae,as,an,b,ib,ie,jb,je,id,jd)
          else
            call sumdel (phil,phil,ibl(l),iel(l),jbl(l),jel(l),ibl(l+1),
     &                   iel(l+1),jbl(l+1),jel(l+1),idl,jdl,idl,jdl)
            call pgs (phil,apl,awl,ael,asl,anl,bl,ibl(l),iel(l),jbl(l),
     &                jel(l),idl,jdl)
          endif
        enddo

!       exit if converged
        call resid (rsd,l2rsd,phi,ap,aw,ae,as,an,b,ib,ie,jb,je,id,jd)
        epsout = l2rsd/l2rsd0
        if (epsout .le. epsin) then
!         deallocate work arrays
          deallocate (apl)
          deallocate (awl)
          deallocate (ael)
          deallocate (asl)
          deallocate (anl)
          deallocate (bl)
          deallocate (rsdl)
          deallocate (phil)
          return
        endif
      enddo

!     deallocate work arrays
      deallocate (apl)
      deallocate (awl)
      deallocate (ael)
      deallocate (asl)
      deallocate (anl)
      deallocate (bl)
      deallocate (rsdl)
      deallocate (phil)

      return
      end

      subroutine ptr (ibl, iel, jbl, jel, lmax, ib, ie, jb, je, ld)
!=======================================================================
!     routine to calculate the array pointers to appropriate space in
!     the work arrays for each level of iteration.
!   input:
!     lmax: guess for lmax
!     ib,ie,jb,je: first and last interior indices in i and j
!     ld:  array dimensions
!   output:
!     lmax: number of levels from finest grid to single grid
!     ibl,iel,jbl,jel: first and last index in i and j for l level
!=======================================================================

      implicit none

      integer ib, ie, jb, je, ld, lmaxin
      integer lmax, l, ncv
      integer ibl(ld), iel(ld), jbl(ld), jel(ld)

      if (lmax .lt. 2) return

      lmaxin = lmax
      lmax = 1
      if (((ie-ib) .lt. 2) .and. ((je-jb) .lt. 2)) return

      lmax = 2
      ibl(lmax) = ib
      iel(lmax) = ie
      jbl(lmax) = jb
      jel(lmax) = jbl(lmax) + (je - jb + 2)/2 - 1
      ncv = (iel(lmax) - ibl(lmax) + 1)*(jel(lmax) - jbl(lmax) + 1)
      if (ncv .eq. (ie - ib + 1)) return

      do l=3,lmaxin
        ibl(l) = ib
        iel(l) = ie
        jbl(l) = jel(l-1) + 1
        jel(l) = jbl(l) + (jel(l-1) - jbl(l-1) + 2)/2 - 1
        ncv = (iel(l) - ibl(l) + 1)*(jel(l) - jbl(l) + 1)
        lmax = l
        if (ncv .eq. (ie - ib + 1)) return
      enddo

      return
      end

      subroutine sumcf (ap2, aw2, ae2, as2, an2, ap1, aw1, ae1, as1,
     &                  an1, ib2, ie2, jb2, je2, ib1, ie1, jb1, je1,
     &                  id2, jd2, id1, jd1)
!=======================================================================
!     routine to sum the fine grid (level 1) coefficients onto the
!     coarse grid (level 2).
!   input:
!     ap1,aw1,ae1,as1,an1: level 1 coefficients for p,w,e,s,n nodes
!     ib1,ie1,jb1,je1,ib2,ie2,jb2,je2: indices for level 1 and 2 grid
!     id1,jd1,id2,jd2: dimensions for level 1 and 2 arrays
!   output
!     ap2,aw2,ae2,as2,an2: level 2 coefficients for p,w,e,s,n nodes
!=======================================================================

      implicit none

      integer ib2, ie2, jb2, je2, ib1, ie1,jb1, je1, id1, id2, jd1, jd2
      integer i1, j1, i2, j2, iwait, jwait
      real(kind=8) ap2(id2,jd2), aw2(id2,jd2), ae2(id2,jd2)
      real(kind=8) as2(id2,jd2), an2(id2,jd2)
      real(kind=8) ap1(id1,jd1), aw1(id1,jd1), ae1(id1,jd1)
      real(kind=8) as1(id1,jd1), an1(id1,jd1)

      do j1=jb1,je1
        if ((j1-jb1+1)/2*2 .ne. (j1-jb1+1)) then
           jwait = 0
        else
           jwait = 1
        endif
        j2 = jb2+(j1-jb1+2)/2-1
        do i1=ib1,ie1
          i2 = i1
          aw2(i2,j2) = aw2(i2,j2) + aw1(i1,j1)
          ae2(i2,j2) = ae2(i2,j2) + ae1(i1,j1)
          as2(i2,j2) = as2(i2,j2) + (1 - jwait)*as1(i1,j1)
          an2(i2,j2) = an2(i2,j2) + (0 + jwait)*an1(i1,j1)
          ap2(i2,j2) = ap2(i2,j2) + ap1(i1,j1) - (0 + jwait)*as1(i1,j1)
     &               - (1 - jwait)*an1(i1,j1)
        enddo
      enddo

      return
      end

      subroutine sumrsd (bl2, rsd1, ib2, ie2, jb2, je2, ib1, ie1, jb1,
     &                   je1, id2, jd2, id1, jd1)
!=======================================================================
!     routine to sum the fine grid (level 1) residuals onto the
!     coarse grid (level 2) source coefficients.
!   input:
!     rsd1: level 1 residual for p node
!     ib1,ie1,jb1,je1,ib2,ie2,jb2,je2: indices for level 1 and 2 grid
!     id1,jd1,id2,jd2: dimensions for level 1 and 2 arrays
!   output:
!     bl2 level 2 source coefficients
!=======================================================================

      implicit none

      integer ib2, ie2, jb2, je2, ib1, ie1, jb1, je1, id1, id2
      integer jd1, jd2, i1, j1, i2, j2
      real(kind=8) bl2(id2,jd2), rsd1(id1,jd1)

      do j1=jb1,je1
        j2 = jb2 + (j1 - jb1 + 2)/2 - 1
        do i1=ib1,ie1
          bl2(i1,j2) = bl2(i1,j2) + rsd1(i1,j1)
        enddo
      enddo

      return
      end

      subroutine sumdel (phi1, phi2, ib1, ie1, jb1, je1, ib2, ie2,
     &                   jb2, je2, id1, jd1, id2, jd2)
!=======================================================================
!     routine to sum the coarse grid (level 2) corrections onto the
!     fine grid (level 1) phi solution.
!   input:
!     phi2: level 2 corrections
!     ib1,ie1,jb1,je1,ib2,ie2,jb2,je2: indices for level 1 and 2 grid
!     id1,jd1,id2,jd2: dimensions for level 1 and 2 arrays
!   output
!     phi1: level 1 solution
!=======================================================================

      implicit none

      integer ib2, ie2, jb2, je2, ib1, ie1, jb1, je1, id1, id2
      integer jd1, jd2, i1, i2, j1, j2
      real(kind=8) phi1(id1,jd1), phi2(id2,jd2)

      do j1=jb1,je1
        j2 = jb2 + (j1 - jb1 + 2)/2 - 1
        do i1=ib1,ie1
          phi1(i1,j1) = phi1(i1,j1) + phi2(i1,j2)
        enddo
      enddo

      return
      end

      subroutine pgs (phi, ap, aw, ae, as, an, b, ib, ie, jb, je,
     &                id, jd)
!=======================================================================
!     subroutine to update the interior phi field by making ntimes
!     sweeps using a point gauss-seidel algorithm with a "cyclic" tdma
!     line solver.
!   input:
!     phi: initial guess for phi
!     ap,aw,ae,as,an: active coefficients for p,w,e,s,n nodes
!     b: accumulated fixed source term
!     ib,ie,jb,je: first and last interior indices in i and j
!     id,jd:  array dimensions
!   output:
!     phi: updated estimate of phi
!=======================================================================

      implicit none

      integer, intent(in) :: ib, ie, jb, je, id, jd
      real(kind=8), intent(in), dimension(id,jd) :: ap
      real(kind=8), intent(in), dimension(id,jd) :: aw
      real(kind=8), intent(in), dimension(id,jd) :: ae
      real(kind=8), intent(in), dimension(id,jd) :: as
      real(kind=8), intent(in), dimension(id,jd) :: an
      real(kind=8), intent(in), dimension(id,jd) :: b
      real(kind=8), intent(inout), dimension(id,jd) :: phi

      integer :: i, j
      real(kind=8), dimension(id,jd) :: bplus

!     south boundary
      j = jb
      do i=ib,ie
        bplus(i,1) = b(i,j) + an(i,j)*phi(i,j+1)
      enddo
      call ctdma (phi(1,j),aw(1,j),ap(1,j),ae(1,j),bplus,ib,ie,id)

!     interior
      do j=jb+2,je-1,2
        do i=ib,ie
          bplus(i,j) = b(i,j) + an(i,j)*phi(i,j+1) + as(i,j)*phi(i,j-1)
        enddo
      enddo
      call new_ctdma (phi,aw,ap,ae,bplus,ib,ie,id,jb+2,je-1,jd)

      do j=jb+1,je-1,2
        do i=ib,ie
          bplus(i,j) = b(i,j) + an(i,j)*phi(i,j+1) + as(i,j)*phi(i,j-1)
        enddo
      enddo
      call new_ctdma (phi,aw,ap,ae,bplus,ib,ie,id,jb+1,je-1,jd)

!     north boundary
      j = je
      do i=ib,ie
        bplus(i,1) = b(i,j) + as(i,j)*phi(i,j-1)
      enddo
      call ctdma (phi(1,j),aw(1,j),ap(1,j),ae(1,j),bplus,ib,ie,id)

      return
      end subroutine pgs

      subroutine ctdma (phi, aw, ap, ae, b, ib, ie, id)
!=======================================================================
!   subroutine to do a cyclic tridiagonal matrix solve
!   input:
!     phi: initial guess for phi
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     ib,ie: first and last interior indices in i
!     id:  array dimensions
!   output:
!     phi: updated estimate of phi
!=======================================================================

      implicit none

      integer ib, ie, id, i
      real(kind=8) phi(id), aw(id), ap(id), ae(id), b(id)
      real(kind=8) factor, alpha(id), beta(id), theta(id)

      alpha(ib) = 2*ap(ib)
      do i=ib+1,ie-1
        alpha(i) = ap(i)
      enddo
      alpha(ie) = ap(ie) + ae(ie)*aw(ib)/ap(ib)

      call tdma (phi, aw, alpha, ae, b, ib, ie, id)

      beta(ib) = -ap(ib)
      do i=ib+1,ie-1
        beta(i) = 0.0
      enddo
      beta(ie) = -ae(ie)

      call tdma (theta, aw, alpha, ae, beta, ib, ie, id)

      factor = (phi(ib) + aw(ib)/ap(ib)*phi(ie))/
     &         (1.+theta(ib) + aw(ib)/ap(ib)*theta(ie))
      do i=ib,ie
        phi(i) = phi(i) - factor*theta(i)
      enddo

      return
      end

      subroutine tdma (phi, aw, ap, ae, b, ib, ie, id)
!=======================================================================
!   subroutine to do a tridiagonal matrix solve for an east-west line
!   input:
!     phi: initial guess for phi
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     ib,ie: first and last interior indices in i
!     id:  array dimensions
!   output:
!     phi: updated estimate of phi
!=======================================================================

      implicit none

      integer ib, ie, id, i
      real(kind=8) phi(id), aw(id), ap(id), ae(id), b(id), alpha(id)
      real(kind=8) beta

      beta = ap(ib)
      phi(ib) = b(ib)/beta
      do i=ib+1,ie
        alpha(i) = -ae(i-1)/beta
        beta = ap(i) + aw(i)*alpha(i)
        phi(i) = (b(i) + aw(i)*phi(i-1))/beta
      enddo
      do i=ie-1,ib,-1
        phi(i) = phi(i) - alpha(i+1)*phi(i+1)
      enddo

      return
      end

      subroutine resid (rsd, l2rsd, phi, ap, aw, ae, as, an, b,
     &                  ib, ie, jb, je, id, jd)
!=======================================================================
!     subroutine to calculate the residual at each interior c.v and
!     the average of the absolute residuals over all interior c.v.
!   input:
!     phi: updated estimate of phi field
!     ap,aw,ae,as,an: active coefficients for p,w,e,s,n nodes
!     b: accumulated fixed source term
!     ntimes: number of sweeps
!     ib,ie,jb,je: first and last interior indices in i and j
!     id,jd:  array dimensions
!   output:
!     rsd: residual array for each interior control volume
!     l2rsd: residual l2 norm for all interior control volume
!=======================================================================

      implicit none

      integer ib, ie, jb, je, id, jd, i, j
      real(kind=8) rsd(id,jd), phi(id,jd), ap(id,jd), aw(id,jd)
      real(kind=8) ae(id,jd), as(id,jd), an(id,jd), b(id,jd)
      real(kind=8) l2rsd

      l2rsd = 0.0

      j = jb
!     south west boundary
      i = ib
      rsd(i,j) = aw(i,j)*phi(ie,j) + ae(i,j)*phi(i+1,j)
     &         + an(i,j)*phi(i,j+1) + b(i,j) - ap(i,j)*phi(i,j)
      l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
!     south boundary interior
      do i=ib+1,ie-1
        rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j)
     &           + an(i,j)*phi(i,j+1) + b(i,j) - ap(i,j)*phi(i,j)
        l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
      enddo
!     south east boundary
      i = ie
      rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(ib,j)
     &         + an(i,j)*phi(i,j+1) + b(i,j) - ap(i,j)*phi(i,j)
      l2rsd = l2rsd + rsd(i,j)*rsd(i,j)

      do j=jb+1,je-1
!       west boundary interior
        i = ib
        rsd(i,j) = aw(i,j)*phi(ie,j) + ae(i,j)*phi(i+1,j)
     &           + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)
     &           + b(i,j) - ap(i,j)*phi(i,j)
        l2rsd = l2rsd +rsd(i,j)*rsd(i,j)
!       interior
        do i=ib+1,ie-1
          rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j)
     &             + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)
     &             + b(i,j) - ap(i,j)*phi(i,j)
          l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
        enddo
!       east boundary interior
        i = ie
        rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(ib,j)
     &           + as(i,j)*phi(i,j-1) + an(i,j)*phi(i,j+1)
     &           + b(i,j) - ap(i,j)*phi(i,j)
        l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
      enddo

      j = je
!     north west boundary
      i = ib
      rsd(i,j) = aw(i,j)*phi(ie,j) + ae(i,j)*phi(i+1,j)
     &         + as(i,j)*phi(i,j-1) + b(i,j) - ap(i,j)*phi(i,j)
      l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
      do i=ib+1,ie-1
!       north boundary interior
        rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(i+1,j)
     &           + as(i,j)*phi(i,j-1) + b(i,j) - ap(i,j)*phi(i,j)
        l2rsd = l2rsd + rsd(i,j)*rsd(i,j)
      enddo
!     north east boundary
      i = ie
      rsd(i,j) = aw(i,j)*phi(i-1,j) + ae(i,j)*phi(ib,j)
     &         + as(i,j)*phi(i,j-1) + b(i,j) - ap(i,j)*phi(i,j)
      l2rsd = l2rsd + rsd(i,j)*rsd(i,j)

      l2rsd = sqrt(l2rsd)

      return
      end

      subroutine new_ctdma (phi, aw, ap, ae, b, ib, ie, id, jb, je, jd)
!=======================================================================
!   subroutine to do a 1 tridiagonal matrix solve
!   input:
!     phi: initial guess for phi
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     ib,ie: first and last interior indices in i
!     id:  array dimensions
!   output:
!     phi: updated estimate of phi
!=======================================================================

      implicit none

      integer, intent(in) :: ib, ie, id
      integer, intent(in) :: jb, je, jd
      real(kind=8), intent(inout), dimension(id,jd) :: phi
      real(kind=8), intent(in), dimension(id,jd) :: aw
      real(kind=8), intent(in), dimension(id,jd) :: ap
      real(kind=8), intent(in), dimension(id,jd) :: ae
      real(kind=8), intent(in), dimension(id,jd) :: b

      integer :: i, j
      real(kind=8) :: factor
      real(kind=8), dimension(id,jd):: alpha
      real(kind=8), dimension(id,jd):: beta
      real(kind=8), dimension(id,jd):: theta

      if ( je .lt. jb ) return

      alpha(ib,jb:je:2) = 2*ap(ib,jb:je:2)

      alpha(ib+1:ie-1,jb:je:2) = ap(ib+1:ie-1,jb:je:2)

      alpha(ie,jb:je:2) = ap(ie,jb:je:2) +
     &                    ae(ie,jb:je:2)*aw(ib,jb:je:2)/ap(ib,jb:je:2)

      call new_tdma (phi, aw, alpha, ae, b, ib, ie, id, jb, je, jd)

      beta(ib,jb:je:2) = -ap(ib,jb:je:2)

      beta(ib+1:ie-1,jb:je:2) = 0.0

      beta(ie,jb:je:2) = -ae(ie,jb:je:2)

      call new_tdma (theta, aw, alpha, ae, beta, ib, ie, id, jb, je, jd)

      do j=jb,je,2
        factor = (phi(ib,j) + aw(ib,j)/ap(ib,j)*phi(ie,j))/
     &           (1.+theta(ib,j) + aw(ib,j)/ap(ib,j)*theta(ie,j))

        phi(ib:ie,j) = phi(ib:ie,j) - factor*theta(ib:ie,j)

      enddo

      return
      end subroutine new_ctdma

      subroutine new_tdma (phi, aw, ap, ae, b, ib, ie, id, jb, je, jd)
!=======================================================================
!   subroutine to do a tridiagonal matrix solve for an east-west line
!   input:
!     phi: initial guess for phi
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     ib,ie: first and last interior indices in i
!     id:  array dimensions
!   output:
!     phi: updated estimate of phi
!=======================================================================

      implicit none

      integer, intent(in) :: ib, ie, id
      integer, intent(in) :: jb, je, jd
      real(kind=8), intent(inout), dimension(id,jd) :: phi
      real(kind=8), intent(in), dimension(id,jd) :: aw
      real(kind=8), intent(in), dimension(id,jd) :: ap
      real(kind=8), intent(in), dimension(id,jd) :: ae
      real(kind=8), intent(in), dimension(id,jd) :: b

      integer :: i, j
      real(kind=8), dimension(jd) :: beta
      real(kind=8), dimension(id,jd) :: alpha

      do j = jb,je,2
        beta(j) = ap(ib,j)
        phi(ib,j) = b(ib,j)/beta(j)
      enddo

      do i=ib+1,ie
!cdir altcode=loopcnt
        do j = jb,je,2
          alpha(i,j) = -ae(i-1,j)/beta(j)
          beta(j) = ap(i,j) + aw(i,j)*alpha(i,j)
          phi(i,j) = (b(i,j) + aw(i,j)*phi(i-1,j))/beta(j)
        enddo
      enddo

      do j = jb,je,2
        do i=ie-1,ib,-1
          phi(i,j) = phi(i,j) - alpha(i+1,j)*phi(i+1,j)
        enddo
      enddo
#endif

      return
      end
