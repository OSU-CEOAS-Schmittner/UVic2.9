      subroutine invtri (z, topbc, botbc, dcb, tdt, kmz, mask, is, ie
     &,                  joff, js, je)

#if defined O_mom
# if defined O_implicitvmix || defined O_isopycmix
!=======================================================================
!     solve the vertical diffusion equation implicitly using the
!     method of inverting a tridiagonal matrix as described in
!     Numerical Recipies in Fortran, The art of Scientific Computing,
!     Second Edition, Press, Teukolsky, Vetterling, Flannery, 1992
!     pages 42,43.
!     this routine assums that the variables are defined at grid points
!     and the top and bottom b.c. are flux conditions.

!     inputs:
!     z         = right hand side terms
!     topbc     = top boundary condition
!     botbc     = bottom boundary condition
!     dcb       = vertical mixing coeff at base of cell
!     tdt       = 2 * timestep
!     kmz       = level indicator
!     mask      = land/sea mask
!     is        = index of starting longitude
!     ie        = index of ending longitude
!     js        = starting latitude row in MW
!     je        = ending latitude row in MW
!     joff      = offset between jrow on disk and j in the MW

!     outputs:
!     z         = returned solution
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, k, km1, kp1, i, is, ie

      real factu, factl, eps

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "vmixc.h"

      integer kmz(imt,jmt)

      real z(imt,km,jmw), topbc(imt,1:jmw), botbc(imt,1:jmw)
      real dcb(imt,km,jsmw:jemw), mask(imt,km,1:jmw)
      real a(imt,km,jsmw:jemw), b(imt,km,jsmw:jemw)
      real c(imt,0:km,jsmw:jemw), d(imt,km,jsmw:jemw)
      real f(imt,0:km,jsmw:jemw), e(imt,km,jsmw:jemw)
      real bet(imt,jsmw:jemw), tdt(km)

      do j=js,je
        jrow = j + joff
        do k=1,km
          km1 = max(1,k-1)
          kp1 = min(k+1,km)
          factu = dztur(k)*tdt(k)*aidif
          factl = dztlr(k)*tdt(k)*aidif
          do i=is,ie
            a(i,k,j) = -dcb(i,km1,j)*factu*mask(i,k,j)
            c(i,k,j) = -dcb(i,k,j)*factl*mask(i,kp1,j)
            f(i,k,j) = z(i,k,j)*mask(i,k,j)
            b(i,k,j) = c1 - a(i,k,j) - c(i,k,j)
          enddo
        enddo
        do i=is,ie
          a(i,1,j)  = c0
          c(i,km,j) = c0
          b(i,1,j)  = c1 - a(i,1,j) - c(i,1,j)
          b(i,km,j) = c1 - a(i,km,j) - c(i,km,j)

!         top and bottom b.c.

          f(i,1,j)  = z(i,1,j) + topbc(i,j)*tdt(1)*dztr(1)*aidif
     &                           *mask(i,1,j)
          k=max(2,kmz(i,jrow))
          f(i,k,j)   = z(i,k,j) - botbc(i,j)*tdt(k)*dztr(k)*aidif
     &                           *mask(i,k,j)
        enddo
      enddo

!     decomposition and forward substitution

      eps = 1.e-30
      do j=js,je
        do i=is,ie
            bet(i,j) = mask(i,1,j)/(b(i,1,j) + eps)
            z(i,1,j) = f(i,1,j)*bet(i,j)
        enddo
        do k=2,km
          do i=is,ie
            e(i,k,j) = c(i,k-1,j)*bet(i,j)
            bet(i,j) = mask(i,k,j)/(b(i,k,j) - a(i,k,j)*e(i,k,j) + eps)
            z(i,k,j) = (f(i,k,j) - a(i,k,j)*z(i,k-1,j))*bet(i,j)
          enddo
        enddo
      enddo

!     back substitution

      do j=js,je
        do k=km-1,1,-1
          do i=is,ie
            z(i,k,j) = z(i,k,j) - e(i,k+1,j)*z(i,k+1,j)
          enddo
        enddo
      enddo
# endif
#endif

      return
      end
