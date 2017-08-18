      subroutine convect_brine (joff, js, je, is, ie)

#if defined O_mom && defined O_ice && defined O_convect_brine
!-----------------------------------------------------------------------
!     solve for brine rejection for each ice category and mix S & T.
!     the brine is first rejected to some depth, defined by
!     potential density difference from the previous time step
!     between the surface and the depth. then convective mixing
!     is applied under each category to both T & S.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!-----------------------------------------------------------------------

      implicit none

      integer k, j, js, je, jrow, joff, i, is, ie, nc, kl, n

      real dens, tq, sq, drodt, drods, drhodt, drhods, ddensdtdt
      real ddensdtds, ddensdsds, cont, tr, trs, r1, dr, dtcbfdzw

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "mw.h"
      include "cembm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "coord.h"
      include "accel.h"
      include "scalar.h"
      include "state.h"
      include "dens.h"

      real tr(imt,km,jmw,nt), trs(imt,km,jmw,nt)

      cont = 0.  ! density contrast to mix brine over (g/cm^3)

      do j=js,je
        jrow = j + joff
        do i=is,ie

          if (cba0(i,jrow) .eq. 1.0) then
!           if no brine flux at all, convect entire cell
            call convct2 (t(1,1,1,1,taup1), joff, j, j, i, i, kmt)
          else

            if (cba0(i,jrow) .ne. 0.0) then
!             convect any unaccounted for area (cba0 should always be 1 or 0)
              tr(i,1:km,j,1:nt) = t(i,1:km,j,1:nt,taup1)
              call convct2 (tr, joff, j, j, i, i, kmt)
              trs(i,1:km,j,1:nt) = tr(i,1:km,j,1:nt)*cba0(i,jrow)
            else
              trs(i,1:km,j,1:nt) = 0.0
            endif

            do nc=0,ncat
              if (cba(i,jrow,nc) .ne. 0.) then
!               convect areas with brine flux (open water is nc=0)
                tr(i,1:km,j,1:nt) = t(i,1:km,j,1:nt,taup1)
                kl = 1
                if (cont .gt. 0) then
!                 find level over which to spread the flux
                  tq = tr(i,kl,j,1) - to(1)
                  sq = tr(i,kl,j,2) - so(1)
                  r1 = dens (tr(i,kl,j,1)-to(1),tr(i,kl,j,2)-so(1), 1)
                  dr = 0.
                  do while (dr .le. cont .and. kl .lt. kmt(i,jrow))
                    kl = kl + 1
                    dr = dens (tr(i,kl,j,1)-to(1),tr(i,kl,j,2)-so(1), 1)-r1
                  enddo
                endif
                dtcbfdzw = c2dtts*cbf(i,jrow,nc)/zw(kl)
                tr(i,1:kl,j,2) = tr(i,1:kl,j,2) + dtxcel(1:kl)*dtcbfdzw

                call convct2 (tr, joff, j, j, i, i, kmt)
                trs(i,1:km,j,1:nt) = trs(i,1:km,j,1:nt)
     &                             + tr(i,1:km,j,1:nt)*cba(i,jrow,nc)
              endif
            enddo
            t(i,1:km,j,1:nt,taup1) = trs(i,1:km,j,1:nt)

          endif

        enddo
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#endif

      return
      end
