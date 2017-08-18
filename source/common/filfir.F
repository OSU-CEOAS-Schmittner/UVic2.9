      subroutine filfir (t, f, s, kl, kind, jtof, jsf, jef)

#if defined O_firfil
!=======================================================================
!     simple finite impulse response filter with [.25, .5, .25] weights
!     modified for assymetric and symmetric boundary conditions

!     input:
!             t    = array of quantity to be filtered along
!                    the first dimension.
!                    note: t(i,k) must be zero where f(i,k) = zero
!                    for this filter to work.
!             f    = mask of zeroes & ones to indicate land
!                    and ocean. zero indicates a land point
!             s    = scratch array
!             kl   = number of vertical levels to filter
!             kind = (0,1) = (symmetric, asymmetric) boundary condition
!                    symmetric is appropriate for tracers & vorticity
!                    asymmetric is appropriate for velocities
!             jtof = number of filter passes per row
!             jsf  = starting row
!             jef  = ending row

!     output:
!             t    = (imt,km) array of filtered quantities
!=======================================================================

      implicit none

      integer kl, j, jsf, jef, kind, num, n, k, i

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"

      integer jtof(jmw)

      real t(imt,kl,jsmw:jemw), f(imt,kl,jsmw:jemw), s(imt,kl,jsmw:jemw)

      do j=jsf,jef
        call setbcx (t(1,1,j), imt, kl)
      enddo
      if (kind .eq. 0) then

!-----------------------------------------------------------------------
!       apply the filter "num" times using a symmetric (no flux)
!       boundary condition
!-----------------------------------------------------------------------

        do j=jsf,jef
          num = jtof(j)
          do n=1,num
            do k=1,kl
              do i=2,imtm1
                s(i,k,j) = f(i,k,j)*(p25*(t(i-1,k,j) + t(i+1,k,j)) +
     &                   t(i,k,j)*(c1 - p25*(f(i-1,k,j) + f(i+1,k,j))))
              enddo
            enddo
            call setbcx (s(1,1,j), imt, kl)
            do k=1,kl
              do i=2,imtm1
                t(i,k,j) = f(i,k,j)*(p25*(s(i-1,k,j) + s(i+1,k,j)) +
     &                   s(i,k,j)*(c1 - p25*(f(i-1,k,j) + f(i+1,k,j))))
              enddo
            enddo
            call setbcx (t(1,1,j), imt, kl)
          enddo
        enddo
      elseif (kind .eq. 1) then

!----------------------------------------------------------------------
!       apply the filter "num" times using an asymmetric (flux)
!       boundary condition
!----------------------------------------------------------------------

        do j=jsf,jef
          num = jtof(j)
          do n=1,num
            do k=1,kl
              do i=2,imtm1
                s(i,k,j) = f(i,k,j)*(p25*t(i-1,k,j) + p5*t(i,k,j) +
     &                               p25*t(i+1,k,j))
              enddo
            enddo
            call setbcx (s(1,1,j), imt, kl)
            do k=1,kl
              do i=2,imtm1
                t(i,k,j) = f(i,k,j)*(p25*s(i-1,k,j) + p5*s(i,k,j)
     &                              +p25*s(i+1,k,j))
              enddo
            enddo
            call setbcx (t(1,1,j), imt, kl)
          enddo
        enddo
      else
        write (stdout,'(/a,i10,a)') ' error=> kind =', kind,' in filfir'
        stop '=>filfir'
      endif
#endif

      return
      end
