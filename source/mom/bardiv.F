      subroutine bardiv
#if defined O_mom
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface

!=======================================================================

!     compute uncorrected barotropic velocities and their divergence
!     for r.h.s. of surface pressure equation

!     Parallel Ocean General Circulation Modeling
!     R.D. Smith, J.K. Dukowicz and R.C. Malone
!     Physica D 60 (1992) 38-61
!     North-Holland
!     Elsevier Science Publishers B.V.

!     Implicit Free-Surface Method for the Bryan-Cox-Semtner Ocean
!     Model.
!     J.K. Dukowicz, R.D. Smith
!     Submitted to J. Geophysical Research   June 1993
!=======================================================================

      implicit none

      integer i, jrow, npt, nislsp

      real factu, factv, fx, fy, d1, d2, utwid, vtwid, tolr, tempu
      real tempv, pnew

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cprnts.h"
      include "diag.h"
      include "emode.h"
      include "grdvar.h"
      include "index.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      character(8) :: bc_symm

!-----------------------------------------------------------------------
!      construct uncorrected barotropic velocities
!      based on gradients of surface pressure
!-----------------------------------------------------------------------

      do i=1,imt
        ubarm1(i,1,1)   = c0
        ubarm1(i,jmt,1) = c0
        ubarm1(i,1,2)   = c0
        ubarm1(i,jmt,2) = c0
      enddo

      do jrow=2,jmtm1
        factu = p5*c2dtsf*csur(jrow)
        factv = p5*c2dtsf*dyur(jrow)

!       leapfrog

        if (leapfrog) then
          do i=2,imtm1
            fx = acor*c2dtsf*cori(i,jrow,1)
            fy = c1/(c1 + fx**2)
            d1 = gam*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
     &       + (c1-gam)*(ps(i+1,jrow+1,2) - ps(i,jrow,2))
            d2 = gam*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
     &       + (c1-gam)*(ps(i+1,jrow,2) - ps(i,jrow+1,2))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)

!           calculate uncorrected velocities at time tau+1

            uhat(i,jrow,1) = fy*(utwid + fx*vtwid) + ubarm1(i,jrow,1)
#  if defined O_implicit_free_surface
     &                    + ubar(i,jrow,1)
#  endif
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid) + ubarm1(i,jrow,2)
#  if defined O_implicit_free_surface
     &                   + ubar(i,jrow,2)
#  endif
          enddo
        elseif (euler1 .or. forward) then

!         forward or 1st pass of euler

          do i=2,imtm1
            fx = acor*c2dtsf*cori(i,jrow,1)
            fy = c1/(c1 + fx**2)
            d1 = ps(i+1,jrow+1,1) - ps(i,jrow,1)
            d2 = ps(i+1,jrow,1) - ps(i,jrow+1,1)
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)

!           calculate uncorrected velocities at time tau+1

            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &             + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &       + ubarm1(i,jrow,2)
          enddo
        elseif (euler2) then

!         euler 2nd pass

          do i=2,imtm1
            fx = acor*c2dtsf*cori(i,jrow,1)
            fy = c1/(c1 + fx**2)
            d1 = theta*(pguess(i+1,jrow+1) - pguess(i,jrow))
     &         + (c1-theta)*(ps(i+1,jrow+1,1) - ps(i,jrow,1))
            d2 = theta*(pguess(i+1,jrow) - pguess(i,jrow+1))
     &       + (c1-theta)*(ps(i+1,jrow,1) - ps(i,jrow+1,1))
            utwid = zu(i,jrow,1)*c2dtsf - factu*(d1 + d2)*dxur(i)
            vtwid = zu(i,jrow,2)*c2dtsf - factv*(d1 - d2)

!           calculate uncorrected velocities at time tau+1

            uhat(i,jrow,1) = fy*(utwid + fx*vtwid)
     &                    + ubarm1(i,jrow,1)
            uhat(i,jrow,2) = fy*(vtwid - fx*utwid)
     &                    + ubarm1(i,jrow,2)
          enddo
        else
          write (stdout,*) '=>Error: leapfrog, euler1, forward, euler2='
     &,   leapfrog, euler1, forward, euler2
          stop '=>bardiv'
        endif

        do i=2,imtm1
          if (kmu(i,jrow) .eq. 0 )then
            uhat(i,jrow,1) = c0
            uhat(i,jrow,2) = c0
          endif
        enddo
      enddo

      call border (uhat(1,1,1), 'u even')
      call border (uhat(1,1,2), 'u odd')

!-----------------------------------------------------------------------
!     filtering of uhat and vhat was removed as suggested by Rick Smith.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     form divergence of uncorrected barotropic velocities for
!     r.h.s. of surface pressure eqn
!-----------------------------------------------------------------------

!     compute the divergence of the sea surface pressure gradients

      call spforc (uhat, dxu, dyu, csu, h, divf)

      fx = c1/(apgr*c2dtsf)
      do jrow=2,jmtm1
        do i=2,imtm1
          divf(i,jrow) = fx*divf(i,jrow)
        enddo
      enddo

#  if defined O_implicit_free_surface
      if (euler2) then

!       correct r.h.s on 2nd pass euler backward

        fx = c1/(grav*dtsf*c2dtsf*apgr)
        do jrow=2,jmtm1
          fy = dyt(jrow)*cst(jrow)*fx
          do i=1,imtm1
            divf(i,jrow) = divf(i,jrow)
     &                   + fy*dxt(i)*(pguess(i,jrow) - ps(i,jrow,1))
          enddo
        enddo
      endif
#  endif

        call border (divf, 't even')

!-----------------------------------------------------------------------
!     calculate guess for change in pressure from guess for pressure
!-----------------------------------------------------------------------

      if (leapfrog) then
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,2)
          enddo
        enddo
      elseif (euler1 .or. forward) then
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = pguess(i,jrow)-ps(i,jrow,1)
          enddo
        enddo
      else
        do jrow=1,jmt
          do i=1,imt
            ptd(i,jrow) = c0
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     use conjugate gradient 9 point laplacian to solve directly for
!     the change in surface pressure.
!-----------------------------------------------------------------------

!     initialize coefficients for the conjugate gradient solver

      call spc9pt (dxu, dyu, csu, h, cf)

!-----------------------------------------------------------------------
!     add diag term to central coeff in the implicit free surface
!-----------------------------------------------------------------------

#  if defined O_implicit_free_surface
      do jrow=2,jmtm1
        fx = cst(jrow)*dyt(jrow)/(apgr*c2dtsf*dtsf*grav)
        do i=2,imtm1
          if (map(i,jrow) .le. 0) then
            cf(i,jrow,0,0) = cf(i,jrow,0,0) - fx*dxt(i)
          endif
        enddo
      enddo
#  endif
      npt    = 9
      variable   = 'd(surf press)'
      nislsp = 0
      bc_symm = 't even'
#  if defined O_implicit_free_surface
      tolr   = tolrfs
#  else
      tolr   = tolrsp
#  endif

      call congr (npt, variable, bc_symm, ptd, ptd, divf, res
     &,           cf
     &,           mxscan, mscan, tolr
     &,           imask, iperm, jperm, iofs, nislsp, nippts
     &,           converged, esterr)

#  if !defined O_implicit_free_surface

!     remove null space from ptd (rigid lid only)

      call checkerboard (ptd, map)
      call border (ptd, bc_symm)

!     remove mean

      call zero_level (ptd, 'surf press', map, dxt, dyt, cst)
      call border (ptd, bc_symm)
#  endif

!-----------------------------------------------------------------------
!     correct barotropic velocities with change in surface pressure
!     gradient. this completes the calculation of barotropic
!     velocity at tau+1
!     also update barotropic velocities for next time step
!-----------------------------------------------------------------------

      do jrow=2,jmtm1
        factu = p5*apgr*c2dtsf*csur(jrow)
        factv = p5*apgr*c2dtsf*dyur(jrow)
        do i=2,imtm1
          d1 = ptd(i+1,jrow+1) - ptd(i,jrow)
          d2 = ptd(i+1,jrow) - ptd(i,jrow+1)

!         tau + 1

          tempu = uhat(i,jrow,1) - factu*(d1 + d2)*dxur(i)
          tempv = uhat(i,jrow,2) - factv*(d1 - d2)
          if (leapfrog) then
#  if defined O_implicit_free_surface
            tempu = tempu - ubar(i,jrow,1)
            tempv = tempv - ubar(i,jrow,2)
#  endif

!           tau - 1  <= tau

            ubarm1(i,jrow,1) = ubar(i,jrow,1)
            ubarm1(i,jrow,2) = ubar(i,jrow,2)
          endif

!         tau  <= tau + 1

          ubar(i,jrow,1) = tempu
          ubar(i,jrow,2) = tempv
        enddo
        do i=2,imtm1
          if (kmu(i,jrow) .eq. 0 ) then
            ubar(i,jrow,1)   = c0
            ubar(i,jrow,2)   = c0
            ubarm1(i,jrow,1) = c0
            ubarm1(i,jrow,2) = c0
          endif
        enddo
      enddo

      call border (ubar(1,1,1), 'u even')
      call border (ubar(1,1,2), 'u odd')
      call border (ubarm1(1,1,1), 'u even')
      call border (ubarm1(1,1,2), 'u odd')

!-----------------------------------------------------------------------
!     update the surface pressure based upon the relaxation solution
!-----------------------------------------------------------------------

      if (leapfrog) then

!       leapfrog

        do jrow=1,jmt
          do i=1,imt
            pnew           = ptd(i,jrow) + ps(i,jrow,2)
            pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
        enddo
      elseif (euler1) then

!       eb 1st pass

        do jrow=1,jmt
          do i=1,imt
#  if defined O_implicit_free_surface
            pnew = ptd(i,jrow) + ps(i,jrow,1)
            pguess(i,jrow) = pnew
            ps(i,jrow,1)   = pnew
#  else
            pnew = ptd(i,jrow) + ps(i,jrow,1)
            pguess(i,jrow) = pnew
#  endif
          enddo
        enddo
      elseif (forward) then

!       forward

        do jrow=1,jmt
          do i=1,imt
            pnew           = ptd(i,jrow) + ps(i,jrow,1)
            pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
        enddo
      elseif (euler2) then

!       eb 2nd pass

        do jrow=1,jmt
          do i=1,imt
            pnew           = ptd(i,jrow) + pguess(i,jrow)
            pguess(i,jrow) = c3*(pnew - ps(i,jrow,1)) + ps(i,jrow,2)
            ps(i,jrow,2)   = ps(i,jrow,1)
            ps(i,jrow,1)   = pnew
          enddo
        enddo
      endif

#  if defined O_remove_ps_checkerboard
#   if !defined O_implicit_free_surface

!     test accumulation of residual checkerboard

      call checkerboard(ps(1,1,1), map)
      call border (ps(1,1,1), bc_symm)
      call zero_level (ps(1,1,1), 'surf press', map, dxt, dyt, cst)
      call border (ps(1,1,1), bc_symm)
#   endif
#  endif
# endif
#endif

      return
      end
