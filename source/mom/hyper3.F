      subroutine hyper3 (npt, variable, bc_symm
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
!     MOM 2 "hypergrid" relax using symmetric coefficients
!     It does not normalize symmetric coefficients ala MOM 1.
!     Hyper3 does checkerboard updating as 4 loops of constant stride

!=======================================================================

!                          H Y P E R G R I D

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

      integer j, i, isle, nisle, n, mscan, mxscan, npt, i1, j1

      real c0, c1, sor, resmax, resis, step, step1, estimated_error
      real crit, cfactor, convergence_rate

      logical converged

      include "size.h"

      integer nippts(mnisle), iofs(mnisle), iperm(maxipp), jperm(maxipp)
      integer map(imt,jmt)

      logical imask(-mnisle:mnisle)

      real dpsi(imt,jmt), forc(imt,jmt), res(imt,jmt)
      real cf(imt,jmt,-1:1,-1:1), relmsk(imt,jmt), guess(imt,jmt)
      real rcfdiag(imt,jmt), diagsum(mnisle), resmi(jmt)

!-----------------------------------------------------------------------
!     set locally needed constants
!-----------------------------------------------------------------------

      c0    = 0.0
      c1    = 1.0

!-----------------------------------------------------------------------
!     calculate "normalized" coefficients used in MOM1

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

      do isle=1,nisle
        diagsum(isle) = c0
      enddo

      do j=2,jmt-1
        do i=2,imt-1
          if (cf(i,j,0,0) .eq. 0.0) then
            rcfdiag(i,j) = c0
          elseif (map(i,j) .eq. 0) then
            rcfdiag(i,j) = c1/cf(i,j,0,0)
          else
            rcfdiag(i,j) = c0
          endif

!         sum diagonal coefficients on island boundary

          isle = -map(i,j)
          if (isle .gt. 0 .and. imask(isle)) then
            diagsum(isle) = diagsum(isle)+cf(i,j,0,0)
          endif
        enddo
      enddo

      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            rcfdiag(i,j) = c1/diagsum(isle)
          enddo
        endif
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
!     set residuals to zero and initialize dpsi
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          res(i,j)  = c0
          dpsi(i,j) = guess(i,j)
        enddo
      enddo

!-----------------------------------------------------------------------
!     begin iteration loop
!-----------------------------------------------------------------------

      do mscan=1,mxscan

        do j=2,jmt-1
          resmi(j) = c0
        enddo

!-----------------------------------------------------------------------
!       consider the arrays as being defined on the squares of an
!       "imt by jmt" checkerboard. take four passes: first solve the
!       equation on the black squares in even columns, then on the red
!       squares in even columns, then red squares in odd columns, and
!       finally on black squares in odd columns..
!-----------------------------------------------------------------------

        if (npt .eq. 5) then

!         5 point calculation

          do i1=0,1
            do j1=0,1
              do j=2+j1,jmt-1,2
                do i=2+i1,imt-1,2
                  res(i,j) =  relmsk(i,j) *
     &                      ((forc(i,j)
     &                       -cf(i,j, 0, 1)*dpsi(i,j+1)
     &                       -cf(i,j, 0,-1)*dpsi(i,j-1)
     &                       -cf(i,j, 1, 0)*dpsi(i+1,j)
     &                       -cf(i,j,-1, 0)*dpsi(i-1,j)
     &                       )*rcfdiag(i,j) - dpsi(i,j) )
                enddo

                call border(res, bc_symm)

!               make a correction to dpsi based on the residuals

                do i=2+i1,imt,2
                  dpsi(i,j) = dpsi(i,j) + sor * res(i,j)
                enddo

!               find the maximum absolute residual to determine convergence

                do i=2+i1,imt,2
                  resmi(j) = max(abs(res(i,j)),resmi(j))
                enddo
              enddo
            enddo
            call border(res, bc_symm)
          enddo
        else

!         9 point calculation

          do i1=0,1
            do j1=0,1
              do j=2+j1,jmt-1,2
                do i=2+i1,imt-1,2
                  res(i,j) =  relmsk(i,j) *
     &                      ((forc(i,j)
     &                       -cf(i,j, 0, 1)*dpsi(i,j+1)
     &                       -cf(i,j, 0,-1)*dpsi(i,j-1)
     &                       -cf(i,j, 1, 0)*dpsi(i+1,j)
     &                       -cf(i,j,-1, 0)*dpsi(i-1,j)
     &                       -cf(i,j, 1, 1)*dpsi(i+1,j+1)
     &                       -cf(i,j,-1, 1)*dpsi(i-1,j+1)
     &                       -cf(i,j, 1,-1)*dpsi(i+1,j-1)
     &                       -cf(i,j,-1,-1)*dpsi(i-1,j-1)
     &                       )*rcfdiag(i,j) - dpsi(i,j) )
                enddo

                call border(res, bc_symm)

!               make a correction to dpsi based on the residuals

                do i=2+i1,imt,2
                  dpsi(i,j) = dpsi(i,j) + sor * res(i,j)
                enddo

!               find the maximum absolute residual to determine convergence

                do i=2+i1,imt,2
                  resmi(j) = max(abs(res(i,j)),resmi(j))
                enddo
              enddo
            enddo
            call border(res, bc_symm)
          enddo
        endif

!-----------------------------------------------------------------------
!       find maximum residual
!-----------------------------------------------------------------------

        resmax = c0
        do j=2,jmt-1
          resmax = max(resmi(j),resmax)
        enddo

!-----------------------------------------------------------------------
!       do integration around each island
!-----------------------------------------------------------------------

        do isle=1,nisle
          if (imask(isle)) then
            resis = c0
            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              resis = resis +  forc(i,j)
     &                      -  cf(i,j, 0, 1)*dpsi(i,j+1)
     &                      -  cf(i,j, 0,-1)*dpsi(i,j-1)
     &                      -  cf(i,j, 1, 0)*dpsi(i+1,j)
     &                      -  cf(i,j,-1, 0)*dpsi(i-1,j)
     &                      -  cf(i,j, 1, 1)*dpsi(i+1,j+1)
     &                      -  cf(i,j,-1, 1)*dpsi(i-1,j+1)
     &                      -  cf(i,j, 1,-1)*dpsi(i+1,j-1)
     &                      -  cf(i,j,-1,-1)*dpsi(i-1,j-1)
            enddo
            resis =  resis / diagsum(isle) - dpsi(i,j)

            resmax = max(abs(resis),resmax)

            do n=1,nippts(isle)
              i = iperm(iofs(isle)+n)
              j = jperm(iofs(isle)+n)
              dpsi(i,j) = dpsi(i,j) + sor * resis
            enddo
          endif
        enddo

        call border(dpsi, bc_symm)

!-----------------------------------------------------------------------
!       test for convergence of the relaxation.
!       the solver is deemed to have converged when the estimated
!       maximum sum of all future corrections does not exceed
!       crit at any point.
!-----------------------------------------------------------------------

        step = sor * resmax
        if (mscan .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. crit) goto 1001
        elseif (step .lt. crit) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(mscan-1))
          estimated_error = step*convergence_rate/(1.0-convergence_rate)
          if (estimated_error .lt. crit)  goto 1001
        endif
      enddo

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
          res(i,j) = sor * res(i,j)
        enddo
      enddo
#endif

      return
      end
