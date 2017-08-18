      subroutine calc_psgrad (psgrad, uext, vext, js, je, is, ie)

# if defined O_mom
!-----------------------------------------------------------------------
!     compute the surface pressure gradients

!     inputs:

!     js   = index of starting row
!     je   = index of ending row
!     is   = index of starting longitude
!     ie   = index of ending longitude

!     outputs:

!     psgrad  = grad(surf press)
!     uext = external mode u (tau+1) for point (ie,je) only
!     vext = external mode v (tau+1) for point (ie,je) only
!-----------------------------------------------------------------------

      implicit none

      integer is, ie, js, je, js1, je1, is1, ie1, jrow, i, kz

      real fxa, r2dtuv, f3, atosp, f2, uext, vext, d1, d2, diag1, diag0
      real diag3, diag4, dubdt, dvbdt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"

      real psgrad(is:ie,js:je,2)

      js1 = max(js,2)
      je1 = min(je,jmt-1)
      is1 = max(is,2)
      ie1 = min(ie,imtm1)

!     on mixing time steps "ptd" has been multiplied by a factor of
!     two and the time step has to be adjusted also.

      if (mod(itt,nmix) .eq. 1) then
        fxa = p5
      else
        fxa = c1
      endif
      r2dtuv = c1/c2dtuv
!     f3    = c2dtuv/c2dtsf
      f3    = c2dtuv
      do jrow=js1,je1
        do i=is1,ie1
          atosp = acor*cori(i,jrow,1)
          f2    = atosp*c2dtuv
          kz = kmu(i,jrow)
          if (kz .ne. 0) then
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
            uext         = ubar(i,jrow,1)
            vext         = ubar(i,jrow,2)
            d1           = ps(i+1,jrow+1,1) - ps(i,jrow,1)
            d2           = ps(i+1,jrow,1) - ps(i,jrow+1,1)
            psgrad(i,jrow,1)  = (d1 + d2)*dxu2r(i)*csur(jrow)
            psgrad(i,jrow,2)  = (d1 - d2)*dyu2r(jrow)
# endif
# if defined O_stream_function
            diag1        = psi(i+1,jrow+1,1)-psi(i  ,jrow,1)
            diag0        = psi(i  ,jrow+1,1)-psi(i+1,jrow,1)
            uext         = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
            vext         =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            diag3        = fxa*(ptd(i+1,jrow+1)-ptd(i  ,jrow))
            diag4        = fxa*(ptd(i  ,jrow+1)-ptd(i+1,jrow))
            dubdt        = (diag3+diag4)*dyu2r(jrow)*hr(i,jrow)
            dvbdt        = (diag3-diag4)*dxu2r(i)*hr(i,jrow)*csur(jrow)
            psgrad(i,jrow,1)=r2dtuv*(dubdt + f3*zu(i,jrow,1) + f2*dvbdt)
            psgrad(i,jrow,2)=r2dtuv*(-dvbdt+ f3*zu(i,jrow,2) + f2*dubdt)
# endif
          else
            psgrad(i,jrow,1) = c0
            psgrad(i,jrow,2) = c0
          endif
        enddo
      enddo

#endif
      return
      end
