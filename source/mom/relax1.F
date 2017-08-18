      subroutine relax1 (npt, variable, bc_symm
     &,                  guess, dpsi, forc, res
     &,                  cf
     &,                  sor, mxscan, mscan, crit
     &,                  imask, iperm, jperm, iofs, nisle, nippts
     &,                  map
     &,                  converged
     &,                  estimated_error
     &                   )

#if defined O_mom
!=======================================================================
!     MOM 2 Relax using symmetric coefficients as input, but
!     normalizes them as in MOM 1.
!     Normalized coefficients are cfn2, cfs2, etc.
!     Uses parallelization trick to get Gauss/Seidel update
!=======================================================================

!                          O L D   R E L A X

!      solve:

!             A * dpsi = forc

!      for "dpsi" with dirichlet boundary conditions (dpsi=const on
!      each component of the boundary) by a "hypergrid" version of
!      Gauss-Seidel iteration.  In this version, the grid is
!      decomposed into 4 sets, each with the same values of
!      (i mod 2, j mod 2).  All calculations within a set may be
!      done in parallel.

!      inputs:
!              npt   = 5 or 9 (active coefficients)
!              variable = character string identifying solution variable
!              bc_symm = equatorial symmetry type (used only when the
!                        symmetry option is on. otherwise ignore it)
!              guess = initial approximation to solution
!              A     = linear operator (assumed symmetric)
!                      typically A is  grad{(1/h)*grad(dpsi)} -
!                      2dt*acor*{grad(f/h) x grad(dpsi)}
!                      using 5 or 9 pt discretizations
!              cf    = imt x jmt x 3 x 3 array of coefficients of A
!              sor   = over-relaxation multiplier
!              forc  = the sum of all terms evaluated at times tau
!                      or tau-1
!              epsilon = convergence criterion
!              max_iterations = maximum number of iterations
!              imask = shows which land masses have perimeter equations
!              iperm = i coordinate of island perimeter points
!              jperm = j coordinate of island perimeter points
!              iofs  = offset in iperm, jperm for start of perimeter
!                      of land_mass(isle)
!              nisle = actual number of land_masses
!              nippts = number of perimeter ocean points for a land_mass
!      output:
!              dpsi   = answer
!              iterations = actual number of iterations performed
!              converged = logical value
!              estimated_error = estimated maximum error in solution
!                          based on step sizes and convergence rate

!=======================================================================

!      more specifically, the equations to be solved are

!             sum (A(ij,i'j') * dpsi(i'j')) = forc(ij)

!      where the subscripts ij and i'j' range over all "free ocean"
!      T cells ij=(i,j) that are not adjacent to land T cells,
!      and one ij=isle for each boundary component of the ocean.

!      with this choice of variables, in the absence of coriolis terms
!      (acor=0), the operator A is symmetric, i.e.,

!             A(ij,i'j') = A(i'j',ij)

!=======================================================================

      implicit none

      character(16) :: variable
      character(*) :: bc_symm

      integer npt, mscan, j, i, isle, nisle, n

      logical converged

      real c0, c1, sor, mxscan, resmax, absmax, resis, step, step1
      real estimated_error, crit, cfactor, convergence_rate

!     dimensions of local arrays
      include "size.h"

      integer nippts(mnisle), iofs(mnisle), iperm(maxipp), jperm(maxipp)
      integer map(imt,jmt)

      logical imask(-mnisle:mnisle)

      real dpsi(imt,jmt), forc(imt,jmt), res(imt,jmt)
      real cf(imt,jmt,-1:1,-1:1), relmsk(imt,jmt), guess(imt,jmt)
      real rncfdiag(imt,jmt), cfn2(imt,jmt), cfs2(imt,jmt)
      real cfe2(imt,jmt), cfw2(imt,jmt), forc2(imt,jmt)
      real diagsum(mnisle)

!-----------------------------------------------------------------------
!     the parallelization tricks used in relax1 work only for 5 pt
!     operators.  do not use relax1 with 9 point operators.
!-----------------------------------------------------------------------

      if (npt .ne. 5) then
        print '(a)', 'WARNING:  relax1 works only with 5 pt operators'
        mscan = 0
        converged = .false.
        stop '=>relax1'
      endif

!-----------------------------------------------------------------------
!     set locally needed constants
!-----------------------------------------------------------------------

      c0    = 0.0
      c1    = 1.0

!-----------------------------------------------------------------------
!     "normalize" coefficients for "oldrelax" method as in MOM1
!     relmsk is now a locally computed array
!     it is 1 on mid-ocean points, and 0 elsewhere
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          if (map(i,j) .eq. 0) then
            relmsk(i,j) = c1
          else
            relmsk(i,j) = c0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     initialize arrays
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          cfn2(i,j)=c0
          cfs2(i,j)=c0
          cfe2(i,j)=c0
          cfw2(i,j)=c0
          rncfdiag(i,j) = c1
        enddo
      enddo

      do isle=1,nisle
        diagsum(isle) = c0
      enddo

      do j=2,jmt-1
        do i=2,imt-1
          if (map(i,j) .eq. 0) then
            rncfdiag(i,j) =
     &         c1/(cf(i,j,0,1)+cf(i,j,0,-1)+cf(i,j,1,0)+cf(i,j,-1,0))

!           normalize coefficients (mid ocean)

            cfn2(i,j) = cf(i,j, 0, 1)*rncfdiag(i,j)
            cfs2(i,j) = cf(i,j, 0,-1)*rncfdiag(i,j)
            cfe2(i,j) = cf(i,j, 1, 0)*rncfdiag(i,j)
            cfw2(i,j) = cf(i,j,-1, 0)*rncfdiag(i,j)
          endif

!         sum diagonal coefficients on island boundary

          if (map(i,j) .le. -1) then
            isle = -map(i,j)
            if (imask(isle)) then
              diagsum(isle) = diagsum(isle)+cf(i,j,0,0)
            endif
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     normalize coefficients on island boundaries
!-----------------------------------------------------------------------

      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            rncfdiag(i,j) = -c1/diagsum(isle)

!           normalize coefficients (island boundary)

              cfn2(i,j) = cf(i,j, 0, 1)*rncfdiag(i,j)
              cfs2(i,j) = cf(i,j, 0,-1)*rncfdiag(i,j)
              cfe2(i,j) = cf(i,j, 1, 0)*rncfdiag(i,j)
              cfw2(i,j) = cf(i,j,-1, 0)*rncfdiag(i,j)
          enddo
        endif
      enddo

!-----------------------------------------------------------------------
!     pre-multiply all coefficients by sor
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          cfn2(i,j) = cfn2(i,j)*sor
          cfs2(i,j) = cfs2(i,j)*sor
          cfe2(i,j) = cfe2(i,j)*sor
          cfw2(i,j) = cfw2(i,j)*sor
        enddo
      enddo

# if defined O_symmetry

!-----------------------------------------------------------------------
!     check for allowable symmetry conditions at equator
!-----------------------------------------------------------------------

      if (bc_symm .ne. 't even' .and. bc_symm .ne. 't odd') then
        write (*,*)
     &'=>Error: type of bc_symm is ',bc_symm,' in hyper3.F'
     &, '... Only allowed type is "t even" or "t odd"'
        stop
      endif
# endif

!-----------------------------------------------------------------------
!     impose boundary conditions on guess
!     dpsi(0) = guess
!-----------------------------------------------------------------------

      call border(guess, bc_symm)

!-----------------------------------------------------------------------
!     set residuals to zero and normalize forcing
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          res(i,j)  = c0
          forc2(i,j) = forc(i,j)*rncfdiag(i,j)
          dpsi(i,j) = guess(i,j)
        enddo
      enddo

!-----------------------------------------------------------------------
!     begin iteration loop
!-----------------------------------------------------------------------

      do mscan=1,mxscan

!-----------------------------------------------------------------------
!       compute residuals without using updated "dpsi" values to get
!       vector of maximum length
!-----------------------------------------------------------------------

        do j=2,jmt-1
          do i=2,imt-1
          res(i,j) = (cfn2(i,j)*dpsi(i,j+1) +
     &                cfs2(i,j)*dpsi(i,j-1) +
     &                cfe2(i,j)*dpsi(i+1,j) +
     &                cfw2(i,j)*dpsi(i-1,j) -
     &                sor*(dpsi(i,j)+forc2(i,j)))*relmsk(i,j)
          enddo
        enddo

        call border(res, bc_symm)

!-----------------------------------------------------------------------
!       correct southern point using updated "dpsi" to get vectors on "i"
!-----------------------------------------------------------------------

        do j=2,jmt-1
          do i=2,imt-1
            res(i,j) = res(i,j) + cfs2(i,j)*res(i,j-1)*relmsk(i,j)
          enddo

!---------------------------------------------------------------------
!       correct western point using updated "dpsi" to get vectors on "j"
!---------------------------------------------------------------------

          do i=2,imt-1
            res(i,j) = res(i,j) + cfw2(i,j)*res(i-1,j)*relmsk(i,j)
          enddo
        enddo

        call border(res, bc_symm)

!---------------------------------------------------------------------
!       make a correction to dpsi based on the residuals
!---------------------------------------------------------------------

        do j=2,jmt-1
          do i=1,imt
            res(i,j)  = res(i,j)*relmsk(i,j)
            dpsi(i,j) = dpsi(i,j) + res(i,j)
          enddo
        enddo

!---------------------------------------------------------------------
!       find the maximum absolute residual to determine convergence
!---------------------------------------------------------------------

        resmax = absmax(res)

!-----------------------------------------------------------------------
!       do a line integral around each island
!---------------------------------------------------------------------

        do isle=1,nisle
          if (imask(isle)) then
            resis = c0
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              resis = resis +  cfn2(i,j)*dpsi(i  ,j+1)
     &                        +cfs2(i,j)*dpsi(i  ,j-1)
     &                        +cfe2(i,j)*dpsi(i+1,j  )
     &                        +cfw2(i,j)*dpsi(i-1,j  )
     &                        -sor*(          forc2(i,j))
            enddo
            resis = resis - sor*dpsi(i,j)

            resmax = max(abs(resis),resmax)

            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              dpsi(i,j) = dpsi(i,j) + resis
            enddo
          endif
        enddo

        call border(dpsi, bc_symm)

!-----------------------------------------------------------------------
!       test for convergence of the relaxation.
!-----------------------------------------------------------------------

        step = resmax

!-----------------------------------------------------------------------
!       the solver is deemed to have converged when the estimated
!       maximum sum of all future corrections does not exceed
!       crit at any point.
!-----------------------------------------------------------------------

        if (mscan .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. crit) goto 1001
        elseif (step .lt. crit) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(mscan-1))
          estimated_error = step*convergence_rate/(1.0-convergence_rate)
          if (estimated_error  .lt. crit)  goto 1001
        endif
      enddo

!---------------------------------------------------------------------
!     end of iteration loop
!---------------------------------------------------------------------

1001  continue
      if (mscan .lt. mxscan) then
        converged = .true.
      else
         converged = .false.
      endif

!---------------------------------------------------------------------
!     return the last increment to dpsi in the argument res
!-----------------------------------------------------------------------

      do i=1,imt
        do j=1,jmt
          res(i,j) = res(i,j)
        enddo
      enddo
#endif

      return
      end
