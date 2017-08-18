      subroutine advupb (u, v, t, is, ie, js, je)

#if defined O_ice && defined O_ice_evp && defined O_embm
!=======================================================================
!     upstream advection of tracers with B-grid velocities

!     input:
!       u    = B-grid u component of ice velocity
!       v    = B-grid v component of ice velocity
!       t    = tracer to be advected
!     output:
!       t = advected tracer
!=====================================================================

      implicit none

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1, n

      real afe, afn, afw, dt, ue, uw, vn

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "grdvar.h"

      real afs(is:ie), tmp(is:ie,js:je), u(is:ie,js:je), v(is:ie,js:je)
      real t(is:ie,js:je)

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      call embmbc (t)
      dt = dts/float(niats)
      do n=1,niats
        tmp(:,:) = t(:,:)
        afs(:) = 0.0
        do j=jsp1,jem1
!         advection velocity on the western face of "T" cells
          uw = (u(1,j-1)*dyu(j-1) + u(1,j)*dyu(j)) * dyt2r(j)
!         advective flux on the western face of "T" cells
          afw = uw*(tmp(1,j) + tmp(2,j))
     &        + abs(uw)*(tmp(1,j) - tmp(2,j))
          do i=isp1,iem1
!           advection velocity on the eastern face of "T" cells
            ue  = (u(i,j-1)*dyu(j-1) + u(i,j)*dyu(j)) * dyt2r(j)
!           advective flux on the eastern face of "T" cells
            afe = ue*(tmp(i,j) + tmp(i+1,j))
     &          + abs(ue)*(tmp(i,j) - tmp(i+1,j))
!           advection velocity on the northern face of "T" cells
            vn  = (v(i-1,j)*dxu(i-1) + v(i,j)*dxu(i))*dxt2r(i)
!           advective flux on the northern face of "T" cells
            afn = vn*(tmp(i,j) + tmp(i,j+1))
     &          + abs(vn)*(tmp(i,j) - tmp(i,j+1))
            t(i,j) = tmp(i,j) - dt*cstr(j)*((afe - afw)*dxt2r(i)
     &               + (afn*csu(j) - afs(i)*csu(j-1))*dyt2r(j))
            afw = afe
            afs(i) = afn
          enddo
        enddo

        call embmbc (t)
      enddo
#endif

      return
      end
