      subroutine tropic (c2dtsf, acor, f, itt, dtts)

#if defined O_mom
!=======================================================================

!        S O L V E   T H E   B A R O T R O P I C   E Q U A T I O N

!   There are several significant changes made in MOM_2 in the
!   calculation of the vertically averaged velocities.

!   1. Based on a 1994 rederivation of the finite difference equations
!     for the stream function by Charles Goldberg (GFDL/Princeton
!     University/Trenton State College), the coefficients in
!     the Poisson equations for dpsi differ slighty from those
!     used in MOM_1.  Designed to reduce the residuals in the
!     finite difference momentum equations, these new equations
!     for dpsi seem to be less stiff.  Tests on a variety of
!     geometries and topographies indicate that the new equations
!     converge more rapidly to the same tolerances in all solvers
!     at all resolutions.

!   2. In all three Poisson solvers for stream function, the MOM_2
!     default mode is that the values of psi and dpsi are not held
!     fixed at zero on the boundary of the "main continent".  They
!     are allowed to float in the same way that other "island
!     boundary" values float.  Although this requires an island
!     integral to be computed on the boundary of every land mass, it
!     also makes the iterative system less stiff, and again results
!     in fewer iterations to converge to the same tolerances
!     in all solvers at all resolutions.

!     Tests indicate that except on architectures where island
!     integrals are prohibitively expensive and the stream function
!     formulation therefore computationally infeasible, the cost
!     of the extra island integrals is small and the savings due
!     to reduced numbers of iterations significant.

!     The user may specify either a land mass on whose boundary the
!     stream function is later normalized to zero or that no post
!     solver normalization is to take place.  Surface pressures are
!     always normalized to have mean zero.  Options are also provided
!     for turning off the island integrals on the boundaries of
!     selected land masses; however, this practice is not recommended.

!   3. The convergence criterion in all three solvers has been
!     changed to the following:  "stop when the predicted maximum error
!     in the solved variable (dpsi or surface pressure) is less than
!     the user specified tolerance."  Convergence tolerances are now in
!     the same units as the variable solved for, and tell the user how
!     many digits of the answer are correct.Thus, if one expects a
!     maximum dpsi of 1.0e12 and desires convergence to 5 significants
!     digits, one chooses a tolerance of 1.0e7. Note that the tolerance
!     used in MOM 2 is NOT the same as "crit" in MOM 1.

!     The maximum error in the solution is predicted as follows:
!     First, convergence of the solver is assumed to be "geometric" in
!     the sense that the maximum absolute correction added to the
!     solved variable in iteration k is modeled as

!           step(k) = step(1)*(convergence_rate)**(k-1)

!     The estimated maximum error in the solution after k iterations
!     is then bounded by the sum of the missing terms in the geometric
!     series truncated after k terms:

!     sum {step(i)} = step(k)*convergence_rate/(1.0 - convergence_rate)
!       i=k+1,infinity

!     Experimental evidence indicates that the convergence rate
!     of an iterative solver remains essentially stable over many
!     iterations, and that the maximum errors when the solvers are
!     stopped by this criterion as compared to solutions obtained
!     by allowing the solvers to run to machine precision are indeed
!     less than the stated tolerances.
!=======================================================================

      implicit none

      integer i, jrow, itt, nconv, luptdb, luptd, npt

      real c2dtsf, dtts, acor, fxa, dpsi1, absmax

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "mw.h"
      include "switch.h"

      character(8) :: bc_symm

      real f(imt,jmt)

      data nconv /0/
      save nconv

      if (eoyear) nconv = 0

      call setbcx (zu(1,1,1), imt, jmt)
      call setbcx (zu(1,1,2), imt, jmt)

# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface

!-----------------------------------------------------------------------
!     solve for the "tau+1" surface pressure and barotropic velocities
!-----------------------------------------------------------------------

      call bardiv

!-----------------------------------------------------------------------
!     test accuracy of solving for change in surface pressure
!-----------------------------------------------------------------------

      if (.not.converged) then
        write (stdout,'(a,i5,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans. max(ps)='
     &, absmax(ps(1,1,2)), ' max(ptd)=',absmax(ptd)
     &, ' estimated max(err)=', esterr
        nconv = nconv + 1
        if (nconv .gt. 50) stop 'nconv > 50 in tropic.f'
      endif
# endif

# if defined O_stream_function

!-----------------------------------------------------------------------
!     construct the forcing for the stream function equation
!-----------------------------------------------------------------------

      call sfforc (zu, dxu, dyu, csu, ztd)

#  if defined O_fourfil || defined O_firfil

!     filter forcing at high latitudes

      call filz (ztd, cf)
#  endif

!-----------------------------------------------------------------------
!     read in solution from 2 previous timesteps for the
!     purpose of computing an initial guess for the present solution.
!-----------------------------------------------------------------------

      luptdb = nkflds - mod(itt,2)
      luptd  = nkflds - 1 + mod(itt,2)
      call oget (kflds, nwds, luptdb, res)
      call oget (kflds, nwds, luptd, ptd)

      fxa=c1
      if (.not. leapfrog) fxa=p5
      do jrow=1,jmt
        do i=1,imt
          ptd(i,jrow) = fxa*(c2*ptd(i,jrow)-res(i,jrow))
        enddo
      enddo
#  if defined O_cyclic
      do jrow=2,jmtm1
        ptd(1,jrow)   = ptd(imtm1,jrow)
        ptd(imt,jrow) = ptd(2,jrow)
      enddo
#  endif

!-----------------------------------------------------------------------
!     choose 5 or 9 point numerics
!-----------------------------------------------------------------------

#  if defined O_sf_5_point

!     initialize coefficients using 5 point numerics

      call sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 5
#  endif
#  if defined O_sf_9_point

!     initialize coefficients using 9 point numerics

      call sfc9pt (acor, f, c2dtsf, dxu, dyu, csu, hr, cf)
      npt = 9
#  endif

!-----------------------------------------------------------------------
!     choose a method for solving for the "tau+1" stream function change
!-----------------------------------------------------------------------

      variable   = 'dpsi'
      bc_symm    = 't odd'

#  if defined O_conjugate_gradient
      call congr (npt, variable, bc_symm, ptd, ptd, ztd, res
     &,           cf
     &,           mxscan, mscan, tolrsf
     &,           imask, iperm, jperm, iofs, nisle, nippts
     &,           converged, esterr)
#  endif

#  if defined O_oldrelax

!     use sequential over-relaxation to solve the 5 pt laplacian
!     as in the codes of Cox (1984) and Semtner (1974).

      call relax1 (npt, variable, bc_symm, ptd, ptd, ztd, res
     &,            cf
     &,            sor, mxscan, mscan, tolrsf
     &,            imask, iperm, jperm, iofs, nisle, nippts
     &,            map
     &,            converged
     &,            esterr
     &            )
#  endif

#  if defined O_hypergrid

!     use sequential over-relaxation to solve the 5 or 9 pt laplacian
!     along diagonals.

      call hyper3 (npt, variable, bc_symm, ptd, ptd, ztd, res
     &,            cf
     &,            sor, mxscan, mscan, tolrsf
     &,            imask, iperm, jperm, iofs, nisle, nippts
     &,            map
     &,            converged
     &,            esterr
     &            )
#  endif

!     correct for drifting dpsi on land mass "imain"

      if (imain .gt. 0) then
        dpsi1 = ptd(iperm(iofs(imain)+1), jperm(iofs(imain)+1))
        call con_adjust (ptd, dpsi1, map)
      endif

!-----------------------------------------------------------------------
!     test accuracy of solving for change in stream function
!-----------------------------------------------------------------------

      if (.not.converged) then
        write (stdout,'(a,i5,3(a,1pe10.3))')
     &  ' WARNING: SOLVER DID NOT CONVERGE in ',mscan
     &, ' scans. max(psi)='
     &, absmax(psi(1,1,2)), ' max(dpsi)=',absmax(ptd)
     &, ' estimated max(err)=', esterr
        nconv = nconv + 1
        if (nconv .gt. 50) stop 'nconv > 50 in tropic.f'
      endif

!-----------------------------------------------------------------------
!     update the stream function based upon the solution
!-----------------------------------------------------------------------

      if (euler2) then
        do jrow=1,jmt
          do i=1,imt
            psi(i,jrow,1) = psi(i,jrow,2) + ptd(i,jrow)
          enddo
        enddo
      else
        do jrow=1,jmt
          do i=1,imt
            res(i,jrow)    = psi(i,jrow,2) + ptd(i,jrow)
            psi(i,jrow,2)  = psi(i,jrow,1)
            psi(i,jrow,1)  = res(i,jrow)
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     save ptd to compute 1st guess for relaxation next timestep
!     (..note.. on 1st pass of euler backward timestep, bypass this
!            save, since it will be done on the 2nd pass)
!     (..note.. on a mixing timestep, alter ptd to be consistent with
!            normal, leap-frog stepping)
!-----------------------------------------------------------------------

      if (.not. euler1) then

        if (.not. leapfrog) then
          do jrow=1,jmt
            do i=1,imt
              ptd(i,jrow)=c2*ptd(i,jrow)
            enddo
          enddo
        endif

        call oput (kflds, nwds, luptdb, ptd)

      endif
# endif

      return
      end

      subroutine sfforc (zu, dxu, dyu, csu, forc)

!=======================================================================

!           S T R E A M   F U N C T I O N   F O R C I N G
!=======================================================================

      implicit none

      real cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0), p5, c0
      parameter (p5=0.5, c0=0.0)

      include "size.h"

      integer i, jrow, i1, j1

      real zu(imt,jmt,2), dxu(imt), dyu(jmt), csu(jmt)
      real forc(imt,jmt), ustuff(imt,jmt), vstuff(imt,jmt)

!-----------------------------------------------------------------------
!     initialize the forcing
!-----------------------------------------------------------------------

      do i=1,imt
        do jrow=1,jmt
          forc(i,jrow) = c0
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     multiply the u eqn by dx*cos, the v eqn by dy, then subtract their
!     partial differences to eliminate the unknown surface pressure from
!     the resulting equation
!-----------------------------------------------------------------------

      do i=1,imt-1
        do jrow=1,jmt-1
          ustuff(i,jrow) = zu(i,jrow,1)*dxu(i)*csu(jrow)
          vstuff(i,jrow) = zu(i,jrow,2)*dyu(jrow)
        enddo
      enddo

      do i1=-1,0
        do j1=-1,0
          do jrow=2,jmt-1
            do i=2,imt-1
              forc(i,jrow) = forc(i,jrow)
     &                     - cddyt(i1,j1)*ustuff(i+i1,jrow+j1)
     &                     + cddxt(i1,j1)*vstuff(i+i1,jrow+j1)
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine sfc5pt (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)

!=======================================================================

!     5  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N

!     coefficient initialization for 5 point elliptic solvers

!     inputs:

!     acor   = implicit coriolis factor (0.0 => 1.0)
!     f      = 2*omega*sin(phi(j))
!     c2dtsf = twice the time step (seconds)
!     dxu    = width of "u" grid cell (cm)
!     dyu    = height of "u" grid cell (cm)
!     csu    = cosine of "u" grid cell
!     hr     = 1/depth at "u" cells (cm)

!     outputs:

!     coeff   = 3 x 3 array of coefficients at each (i,j) point
!=======================================================================

      implicit none

      integer jj, ii, j, i, i1, j1, i2, j2

      real p5, c0, c2dtsf, acor
      parameter (p5=0.5, c0=0.0)
      real cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
      real f(imt,jmt), coef(imt,jmt,-1:1,-1:1)
      real ustuff(imt,jmt), vstuff(imt,jmt)

!-----------------------------------------------------------------------
!     initialize the coefficients
!-----------------------------------------------------------------------

      do jj=-1,1
        do ii=-1,1
          do j=1,jmt
            do i=1,imt
              coef(i,j,ii,jj) = c0
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     compute coefficients for all points
!-----------------------------------------------------------------------

      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 5 point coefficients

!     note that ne [and nw] coefficient adds to n coefficient in
!     ustuff term, but ne [and se] coefficient adds to e coefficient in
!     vstuff term for the 5 point operator.
!-----------------------------------------------------------------------

      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do j=2,jmt-1
                do  i=2,imt-1
                  coef(i,j,0,j1+j2) = coef(i,j,0,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                  coef(i,j,i1+i2,0) = coef(i,j,i1+i2,0) +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     augment coefficients for implicit treatment of coriolis term
!     all coefficients are calculated, but corner ones are zero.
!-----------------------------------------------------------------------

      if (acor .ne. 0.0) then
        do i=1,imt-1
          do j=1,jmt-1
            ustuff(i,j) = acor*hr(i,j)*(-f(i,j))
            vstuff(i,j) = acor*hr(i,j)*( f(i,j))
          enddo
        enddo
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do j=2,jmt-1
                  do  i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      return
      end

      subroutine sfc9pt  (acor, f, c2dtsf, dxu, dyu, csu, hr, coef)

!=======================================================================

!     9  P T    C O E F F I C I E N T   I N I T I A L I A Z A T I O N

!     coefficient initialization for 9 point elliptic solvers

!     inputs:

!     acor   = implicit coriolis factor (0.0 => 1.0)
!     f      = 2*omega*sin(phi(j))
!     c2dtsf = twice the time step (seconds)
!     dxu    = width of "u" grid cell (cm)
!     dyu    = height of "u" grid cell (cm)
!     csu    = cosine of "u" grid cell
!     hr     = 1/depth at "u" cells (cm)

!     outputs:

!     coeff   = 3 x 3 array of coefficients at each (i,j) point
!=======================================================================

      implicit none

      integer jj, ii, j, i, i1, j1, i2, j2

      real c0, p5, c2dtsf, acor
      parameter (c0=0.0, p5=0.5)
      real cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real csu(jmt), dxu(imt), dyu(jmt), hr(imt,jmt)
      real f(imt,jmt)
      real coef(imt,jmt,-1:1,-1:1)
      real ustuff(imt,jmt), vstuff(imt,jmt)

!-----------------------------------------------------------------------
!     initialize the work area
!-----------------------------------------------------------------------

      do jj=-1,1
        do ii=-1,1
          do j=1,jmt
            do i=1,imt
              coef(i,j,ii,jj) = c0
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     generate arrays of coefficients

!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     compute coefficients for all points
!-----------------------------------------------------------------------

      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = dxu(i)*csu(j)*hr(i,j) / (c2dtsf*dyu(j))
          vstuff(i,j) = dyu(j)*hr(i,j) / (c2dtsf*dxu(i)*csu(j))
        enddo
      enddo

!---------------------------------------------------------------------
!     calculate 9 point coefficients
!---------------------------------------------------------------------

      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do j=2,jmt-1
                do  i=2,imt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2) +
     &               cddyu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)  +
     &               cddxu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     augment coefficients for implicit treatment of coriolis term
!     all coefficients are calculated, but corner ones are zero.
!-----------------------------------------------------------------------

      if (acor .ne. 0.0) then
        do i=1,imt-1
          do j=1,jmt-1
            ustuff(i,j) = acor*hr(i,j)*(-f(i,j))
            vstuff(i,j) = acor*hr(i,j)*( f(i,j))
          enddo
        enddo
        do i1=0,1
          do j1=0,1
            do i2=-1,0
              do j2=-1,0
                do j=2,jmt-1
                  do  i=2,imt-1
                    coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &               - cddxu(i1,j1)*cddyt(i2,j2)*ustuff(i+i2,j+j2)
     &               - cddyu(i1,j1)*cddxt(i2,j2)*vstuff(i+i2,j+j2)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      endif

      return
      end

      subroutine spforc (zu, dxu, dyu, csu, h, forc)

!=======================================================================
!           S U R F A C E   P R E S S U R E   F O R C I N G
!=======================================================================

      implicit none

      integer i, j, i1, j1

      real p5
      parameter (p5=0.5)
      real cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real zu(imt,jmt,2), dxu(imt), dyu(jmt), csu(jmt), h(imt,jmt)
      real forc(imt,jmt), ustuff(imt,jmt), vstuff(imt,jmt)

!-----------------------------------------------------------------------
!     generate arrays of coefficients

!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!     weight "zu" and "zv" by the cell area and take the divergence

      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = h(i,j)*zu(i,j,1)*dyu(j)
          vstuff(i,j) = h(i,j)*zu(i,j,2)*dxu(i)*csu(j)
        enddo
      enddo

      do i=1,imt
        do j=1,jmt
          forc(i,j) = 0.0
        enddo
      enddo

      do i1=-1,0
        do j1=-1,0
          do i=2,imt-1
            do j=2,jmt-1
              forc(i,j) = forc(i,j) + cddxt(i1,j1)*ustuff(i+i1,j+j1)
     &                              + cddyt(i1,j1)*vstuff(i+i1,j+j1)
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine spc9pt (dxu, dyu, csu, h, coef)

!=======================================================================

!     S U R F A C E   P R E S S U R E    C O E F F I C I E N T

!                I N I T I A L I A Z A T I O N

!     inputs:

!     dxu    = width of "u" grid cell (cm)
!     dyu    = height of "u" grid cell (cm)
!     csu    = cosine of "u" grid cell
!     h      = depth at "u,v" cells (cm)

!     outputs:

!     coeff   = 3 x 3 array of coefficients at each (i,j) point
!=======================================================================

      implicit none

      integer i1, j1, i, j, i2, j2

      real c0, p5
      parameter (c0=0.0, p5=0.5)
      real cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real csu(jmt), dxu(imt), dyu(jmt), h(imt,jmt)
      real coef(imt,jmt,-1:1,-1:1)
      real ustuff(imt,jmt), vstuff(imt,jmt)

!-----------------------------------------------------------------------
!     generate arrays of coefficients

!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     compute coefficients for all points
!-----------------------------------------------------------------------

!     initialize all 9 coefficients to zero

      do i1=-1,1
        do j1=-1,1
          do i=1,imt
            do j=1,jmt
              coef(i,j,i1,j1) = c0
            enddo
          enddo
        enddo
      enddo

      do j=1,jmt
        do i=1,imt
          ustuff(i,j) = 0.0
          vstuff(i,j) = 0.0
        enddo
      enddo
      do i=1,imt-1
        do j=1,jmt-1
          ustuff(i,j) = h(i,j)*dyu(j)/(dxu(i)*csu(j))
          vstuff(i,j) = h(i,j)*dxu(i)*csu(j)/dyu(j)
        enddo
      enddo

!     calculate divergence = ddx (ddx (ustuff)) + ddy( ddy (vstuff))

      do i1=0,1
        do j1=0,1
          do i2=-1,0
            do j2=-1,0
              do i=2,imt-1
                do j=2,jmt-1
                  coef(i,j,i1+i2,j1+j2) = coef(i,j,i1+i2,j1+j2)
     &                 + cddxu(i1,j1) * cddxt(i2,j2) * ustuff(i+i2,j+j2)
     &                 + cddyu(i1,j1) * cddyt(i2,j2) * vstuff(i+i2,j+j2)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine filz (fext, cf)
# if defined O_fourfil || defined O_firfil

!=======================================================================
!     subroutine filz sets up input needed for fourier filtering
!     (when the "fourfil" ifdef is defined) -or- symmetric finite
!     impulse response filtering (when the "firfil" ifdef is defined)
!     of "fext" at the specified high latitudes. "fext" is forcing for
!     the external mode.
!=======================================================================

      implicit none

      integer jrow, jj, l, is, ie, ii, i, im, m, n, isv, iev, j, num

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "index.h"
      include "levind.h"
#  if defined O_firfil

      integer jtof(jmt)
#  endif

      real fext(imt,jmt)
      real cf(imt,jmt,3)
#  if defined O_fourfil
      real temp(imt)
#  endif

!=======================================================================

#  if defined O_fourfil
      do jrow=jfrst,jmtm1
        if (jrow.le.jft1 .or. jrow.ge.jft2) then

!---------------------------------------------------------------------
!         fourier filter fext at high latitudes
!---------------------------------------------------------------------

            jj = jrow - jfrst + 1
            if (jrow .ge. jft2) jj = jj - jskpt + 1
            do l=1,lsegf
              is = iszf(jj,l)
              if (is .ne. 0) then
                ie = iezf(jj,l)
                do ii=is,ie
                  i = mod(ii-2,imtm2) + 2
                  temp(ii+1-is) = fext(i,jrow)
                enddo
                im = ie-is+1
#   if defined O_cyclic
                if (im .ne. imtm2) then
                   m = 1
                   n = nint(im*cst(jrow)*cstr(jft0))
                else
                   m = 3
                   n = nint(im*cst(jrow)*cstr(jft0)*p5)
                endif
#   else
                m = 1
                n = nint(im*cst(jrow)*cstr(jft0))
#   endif

                call filtr (temp(1), im, m ,n, 0)

                do ii=is,ie
                  i = mod(ii-2,imtm2)+2
                  fext(i,jrow) = temp(ii+1-is)
                enddo
              endif
            enddo
        endif
      enddo
#  endif
#  if defined O_firfil

!-----------------------------------------------------------------------
!     filter external mode forcing at high latitudes with
!     symmetric finite impulse response filter
!-----------------------------------------------------------------------

!     construct mask and number of jrows to filter "jtof"

      do jrow=1,jmt
        jtof(jrow) = 0
        do i=1,imt
          cf(i,jrow,1) = c0
        enddo
      enddo
      do jrow=jfrst,jmt-1
        if (jrow.le.jft1 .or. jrow.ge.jft2) then
            jj = jrow - jfrst + 1
            if (jrow .ge. jft2) jj = jj - jskpt + 1
            jtof(jrow) = numflt(jj)

!           construct appropriate land/sea mask
!           (1,0) for (filtered,non filtered) points

            do l=1,lseg
              isv = iszf(jj,l)
              if (isv .ne. 0) then
                iev = iezf(jj,l)
                do ii=isv,iev
                  i = mod(ii-2,imt-2) + 2
                  cf(i,jrow,1) = c1
                enddo
              endif
            enddo

        endif
      enddo
      call setbcx (cf(1,1,1), imt, jmt)

!     select points to filter (non filtered points = zero)

      do jrow=1,jmt
        do i=2,imtm1
          cf(i,jrow,2) = fext(i,jrow)*cf(i,jrow,1)
        enddo
      enddo
      call setbcx (cf(1,1,2), imt, jmt)

!     each filtering consists of a double pass

      do j=2,jmtm1
        num = jtof(j)
        do n=1,num
          do i=2,imtm1
            cf(i,j,3) = cf(i,j,1)*(p25*(cf(i-1,j,2) + cf(i+1,j,2)) +
     &               cf(i,j,2)*(c1 - p25*(cf(i-1,j,1) + cf(i+1,j,1))))
          enddo
#  if defined O_cyclic
          cf(1,j,3)   = cf(imtm1,j,3)
          cf(imt,j,3) = cf(2,j,3)
#  else
          cf(1,j,3)   = c0
          cf(imt,j,3) = c0
#  endif
          do i=2,imtm1
            cf(i,j,2) = cf(i,j,1)*(p25*(cf(i-1,j,3) + cf(i+1,j,3)) +
     &               cf(i,j,3)*(c1 - p25*(cf(i-1,j,1) + cf(i+1,j,1))))
          enddo
#  if defined O_cyclic
          cf(1,j,2)   = cf(imtm1,j,2)
          cf(imt,j,2) = cf(2,j,2)
#  else
          cf(1,j,2)   = c0
          cf(imt,j,2) = c0
#  endif
        enddo
      enddo

!     restore filtered "fext" on ocean points

      do jrow=2,jmtm1
        do i=1,imt
          if (cf(i,jrow,1) .ne. c0) fext(i,jrow) = cf(i,jrow,2)
        enddo
      enddo
      call setbcx (fext, imt, jmt)
#  endif

# endif
#endif
      return
      end
