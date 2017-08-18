      subroutine neptune

#if defined O_mom && defined O_neptune
!=======================================================================
!     calculate neptune (maximum entropy) velocities

!     the option "neptune" provides a subgridscale parameterization
!     for the interaction of eddies and topography

!     reference:
!       Holloway, G., 1992: representing topographic stress for large
!     scale ocean models, J. Phys. Oceanogr., 22, 1033-1046

!     neptune is calculated as an equilibrium streamfunction given by
!     pnep=-f*snep*snep*hnep and is applied through eddy viscosity

!     hnep = model streamfunction depth
!     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))

!     the neptune length scale snep has a value of senep at the
!     equator and smoothly changes to spnep at the poles
!=======================================================================

      implicit none

      integer i, jrow

      real tl, f, snep, hnep, diag1, diag0

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cnep.h"
      include "emode.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"

      integer kmz(imt,jmt)

!     compute an array to indicate "interior" streamfunction grid boxes

      do i=1,imt
        kmz(i,1) = 0
      enddo

      do jrow=1,jmt
        kmz(1,jrow) = 0
      enddo

      do jrow=2,jmt
        do i=2,imt
          kmz(i,jrow) = min(kmu(i-1,jrow-1), kmu(i,jrow-1),
     &                      kmu(i-1,jrow), kmu(i,jrow))
        enddo
      enddo

!     calculate the topographic stress equilibrium streamfunction

!     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
!     pnep  = -f*snep*snep*hnep

      do jrow=2,jmtm1
        do i=2,imtm1
          tl = tlat(i,jrow)/radian
          f    = c2*omega*sin(tl)
          snep = spnep + (senep - spnep)*
     &           (p5 + p5*cos(c2*tl))

!         find depth on streamfunction grid

          hnep = 0
          if (kmz(i,jrow) .ne. 0) then
            hnep = zw(kmz(i,jrow))
          endif

          pnep(i,jrow) = -f*snep*snep*hnep

        enddo
# if defined O_cyclic
          pnep(1,jrow) = pnep(imtm1,jrow)
          pnep(imt,jrow) = pnep(2,jrow)
# endif
      enddo

!     calculate depth independent velocity components from pnep

      do jrow=2,jmtm1
        do i=2,imtm1
          diag1 = pnep(i+1,jrow+1) - pnep(i  ,jrow)
          diag0 = pnep(i  ,jrow+1) - pnep(i+1,jrow)
          unep(i,jrow,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
          unep(i,jrow,2) =  (diag1-diag0)*dxu2r(jrow)*hr(i,jrow)
     &                      *csur(jrow)
        enddo
# if defined O_cyclic
          unep(1,jrow,1) = unep(imtm1,jrow,1)
          unep(1,jrow,2) = unep(imtm1,jrow,2)
          unep(imt,jrow,1) = unep(2,jrow,1)
          unep(imt,jrow,2) = unep(2,jrow,2)
# endif
      enddo

#endif
      return
      end
