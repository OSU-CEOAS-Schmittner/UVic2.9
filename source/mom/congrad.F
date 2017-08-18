      subroutine congr  (npt, variable, bc_symm
     &,                  guess, dpsi, forc, res
     &,                  cf
     &,                  max_iterations, iterations, epsilon
     &,                  imask, iperm, jperm, iofs, nisle, nippts
     &,                  converged
     &,                  estimated_error
     &                  )

# if defined O_mom
!=======================================================================

!                            C O N G R A D

!      solve:

!             A * dpsi = forc

!      for "dpsi" with dirichlet boundary conditions (dpsi=const on
!      each component of the boundary) by a preconditioned conjugate
!      gradient algorithm.

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

!      based on the preconditioned conjugate gradient algorithm given
!      in:

!     A Reformulation and Implementation of the Bryan-Cox-Semtner
!     Ocean Model on the Connection Machine
!     J.K. Dukowicz, R.D. Smith, and R.C. Malone
!     Journal of Atmospheric and Oceanic Technology
!     Vol 10. No. 2 April 1993

!=======================================================================

!      more specifically, the equations to be solved are

!             sum (A(ij,i'j') * dpsi(i'j')) = forc(ij)

!      where the subscripts ij and i'j' range over all "free ocean"
!      T cells ij=(i,j) that are not adjacent to land T cells,
!      and one ij=isle for each boundary component of the ocean.

!      with this choice of variables, in the absence of coriolis terms
!      (acor=0), the operator A is symmetric, i.e.,

!             A(ij,i'j') = A(i'j',ij)

!      the algorithm (essentially executable in Fortran 90) is...

!      subroutine congrad (A, guess, forc, dpsi, iterations)

!      use matrix_module

!      intent (in)     :: A, guess, forc
!      intent (out)    :: dpsi, iterations

!      type(dpsi_type) :: guess, dpsi, Zres, s
!      type(res_type)  :: res, As, forc
!      type(operator)  :: A
!      type(inv_op)    :: Z
!      dimension (0:max_iterations) :: dpsi, res, s, As, beta, alpha

!      dpsi(0) = guess
!      res(0)  = forc - A * dpsi(0)
!      beta(0) = 1
!      s(0)    = zerovector()
!      do k = 1 to max_iterations
!        Zres(k-1) = Z * res(k-1)
!        beta(k)   = res(k-1) * Zres(k-1)
!        s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
!        As(k)     = A * s(k)
!        alpha(k)  = beta(k) / (s(k) * As(k))
!        dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
!        res(k)    = res(k-1) - alpha(k) * As(k)
!        estimated_error = err_est(k, alpha(k), s(k))
!        if (estimated_error) < epsilon) exit
!      enddo
!      if (k > max_iterations) then
!        print *, 'did not converge in ',k,' iterations'
!        stop '=>congrad'
!      endif

!      iterations = k
!      dpsi = dpsi(k)

!      end

!      where...

!      the "vector" and "operator" types used in conjugate gradient
!      are mapped to ordinary 2-dimensional fortran arrays as follows:

!      type(dpsi_type) :: guess, dpsi, Zres, s
!          if ij=(i,j) is a mid-ocean point, map dpsi(ij)-->dpsi(i,j)
!          if ij=isle is an ocean boundary subscript, replicate the
!          value dpsi(isle) in dpsi(i,j) for each (i,j) in the ocean
!          perimeter of land_mass(isle).  the arrays iperm(isle) and
!          jperm(isle), along with iofs(isle) locate these ocean
!          perimeter T cells.
!      type(res_type)  :: res, As, forc
!          if ij=(i,j) is a mid-ocean point, res(ij)-->res(i,j)
!          if ij=isle is an ocean boundary subscript, the value of
!          res(isle) = sum (res(i,j))
!          where the sum is taken over all (i,j) in the ocean perimeter
!          of land_mass(isle).  sometimes, the computed values
!          res(i,j) represent contributions of T cell (i,j) to the
!          component res(isle), and sometimes the values are balanced
!          so that res(i,j)=res(isle)/nippts(isle).  note that, even
!          when balanced, the relation between type(res_type) variables
!          res(isle) and res(i,j) differs from that of type(dpsi_type)
!          variables dpsi(isle) and dpsi(i,j) on T cells in the ocean
!          perimeter.
!      type(operator)  :: A
!          the nearly diagonal quality of the operators used
!          permits a representation as a small collection of
!          2-dimensional arrays.
!          the diagonal, A(ij,ij), is stored in an array cfdiag(i,j)
!          as follows:
!          if ij=(i,j) is a mid-ocean point, A(ij,ij) = cfdiag(i,j)
!          if ij=isle is an ocean boundary subscript,
!          A(isle, isle) = sum (cfdiag(i,j))
!          where the sum is taken over all (i,j) in the ocean perimeter
!          of land_mass(isle).  each cfdiag(i,j) represents the contribution
!          of T cell (i,j) to the island variable diagonal coefficient.
!          the off-diagonal terms A(ij,i`j`) are stored in 4 arrays
!          cfn, cfs, cfe, and cfw if A is a 5-point operator, and in
!          these and 4 additional arrays, cfne, cfnw, cfse, cfsw, if
!          A is a 9-point operator.  For example, if i`=i and j`=j+1,
!          then A(ij,i`j`) is stored in cfn(i,j).
!          if ij=(i,j) is a mid-ocean point and i`j`=isle` is and ocean
!          perimeter subscript, with i`=i and j`=j+1, then
!          cfn(i,j)=A(ij,isle`) is the coefficient of the island
!          variable dpsi(isle`) in the equation for mid-ocean point
!          dpsi(ij)=dpsi(i,j).
!          if ij=isle is an ocean perimeter point and i`j`=(i`,j`) is
!          a mid-ocean point, with i`=i and j`=j-1, then
!          cfs(i,j)=A(isle,i`j`) is the coefficient of the mid-ocean
!          variable dpsi(i`j`)=dpsi(i,j) in the equation for the island
!          variable dpsi(isle).  note that equations for island
!          variables dpsi(isle) are "non-local" in the sense that
!          they usually contain more than 5 or 9 terms, some of which
!          involve values dpsi(i`j`) outside of a compact 5-point
!          or 9-point neighborhood.
!      type(inv_op)    :: Z
!          the approximate inverse operator Z used at present is a
!          diagonal operator Z(ij,ij) = 1/A(ij,ij).
!          if ij=(i,j) is a mid-ocean point,
!            then Z(i,j)=Z(ij)=1/A(ij)=1/cfdiag(i,j)
!          if ij=isle is an ocean perimeter point, then
!          Z(isle) is replicated at each ocean perimeter T cell
!          bordering land_mass(isle).
!            Z(i,j)=Z(isle)=1/A(isle)=1/sum(A(i,j))

!=======================================================================

      implicit none

      character(16) :: variable
      character(*) :: bc_symm

      integer nerror, j, nisle, npt, k, imax, jmax, max_iterations
      integer is, ie, js, je, js1, je1, is1, ie1, jrow, i, kz

      logical converged, diverging

      real zresmax, absvecmax, absmax, epsilon, estimated_error
      real betakm1, betak, dot2, betak_min, smax, step, alpha
      real convergence_rate, betaquot, s_dot_as, step1, cfactor
      real fxa, r2dtuv, f3, atosp, f2, uext, vext, d1, d2
      real diag1, diag0, diag3, diag4, dubdt, dvbdt

      include "size.h"

      integer iperm(maxipp), jperm(maxipp), iofs(mnisle)
      integer iterations, nippts(mnisle)

      logical imask(-mnisle:mnisle)

      real guess(imt,jmt), dpsi(imt,jmt), Zres(imt,jmt)
      real s(imt,jmt), res(imt,jmt), As(imt,jmt), forc(imt,jmt)
      real cf(imt,jmt,-1:1,-1:1), Z(imt,jmt)

# if defined O_symmetry

!-----------------------------------------------------------------------
!     check for allowable symmetry conditions at equator
!-----------------------------------------------------------------------

      if (bc_symm .ne. 't even' .and. bc_symm .ne. 't odd') then
        write (*,*)
     &'=>Error: type of bc_symm is ',bc_symm,' in congrad.F'
     &, '... Only allowed type is "t even" or "t odd"'
        stop
      endif
# endif

!-----------------------------------------------------------------------
!     impose boundary conditions on guess
!     dpsi(0) = guess
!-----------------------------------------------------------------------

      call border(guess, bc_symm)

      do i=1,imt
        do j=1,jmt
          dpsi(i,j) = guess(i,j)
        enddo
      enddo

!-----------------------------------------------------------------------
!     make approximate inverse operator Z (always even symmetry)
!-----------------------------------------------------------------------

      call make_inv (cf, Z,
     &               imask, iperm, jperm, iofs, nisle, nippts)
      call border(Z, 't even')

!-----------------------------------------------------------------------
!     res(0)  = forc - A * dpsi(0)
!     impose cyclic and/or symmetry conditions on res(i,j)
!-----------------------------------------------------------------------

      if (npt .eq. 5) then
        call op5_vec(cf, dpsi, res)
      else
        call op9_vec(cf, dpsi, res)
      endif
      do i=2,imt-1
        do j=2,jmt-1
          res(i,j) = forc(i,j) - res(i,j)
        enddo
      enddo

      call border(res, bc_symm)

!-----------------------------------------------------------------------
!     Zres(k-1) = Z * res(k-1)
!     see if guess is a solution, bail out to avoid division by zero
!-----------------------------------------------------------------------

      k = 0
      call inv_op(Z, res, Zres,
     &            imask, iperm, jperm, iofs, nisle, nippts)

!     set borders of Zres using cyclic/symmetry, if defined.

      call border(Zres, bc_symm)

      Zresmax = absmax(Zres)

      diverging = .false.

!     Assume convergence rate of 0.99 to extrapolate error

      if (100.0 * Zresmax .lt. epsilon) then
        estimated_error = 100.0 * Zresmax
        goto 101
      endif

!-----------------------------------------------------------------------
!     beta(0) = 1
!     s(0)    = zerovector()
!-----------------------------------------------------------------------

      betakm1 = 1.0
      call zero_vec(s)

!-----------------------------------------------------------------------
!     begin iteration loop
!-----------------------------------------------------------------------

      do k = 1,max_iterations

!-----------------------------------------------------------------------
!       Zres(k-1) = Z * res(k-1)
!-----------------------------------------------------------------------

        call inv_op(Z, res, Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)

!       set borders of Zres using cyclic/symmetry, if defined.

        call border(Zres, bc_symm)

!-----------------------------------------------------------------------
!       beta(k)   = res(k-1) * Zres(k-1)
!-----------------------------------------------------------------------

        betak = dot2(Zres, res)
        if (k .eq. 1) then
          betak_min = abs(betak)
        elseif (k .gt. 2) then
          betak_min = min(betak_min, abs(betak))
          if (abs(betak) .gt. 100.0*betak_min) then
            write (*,'(/2(a/))')
     &      'WARNING: conjugate gradient solver terminated because'
     &,     '         correction steps are diverging.'
            write (*,'(/7(a/))')
     &      'PROBABLE CAUSES:'
     &,     '         1. convergence criterion is too tight...'
     &,     '            roundoff error prevents convergence'
     &,     '     or  2. the solution is beginning to blow up...'
     &,     '            if so, it is extremely unlikely that usable'
     &,     '            results can be obtained in subsequent time'
     &,     '            steps.'
            write (*,'(/3(a/))')
     &      'ERROR:   It is assumed that the solution is blowing up.'
     &,     '         It is extremely unlikely that usable results can'
     &,     '         be obtained in subsequent time steps.'
            if (variable .ne. 'surfpres') then
              stop '==>congrad'
            endif
            diverging = .true.
            smax = absmax(s)
            step = abs(alpha) * smax
            estimated_error=step*convergence_rate/(1.0-convergence_rate)
            go to 101
          endif
        endif

!-----------------------------------------------------------------------
!       s(k)      = Zres(k-1) + (beta(k)/beta(k-1)) * s(k-1)
!-----------------------------------------------------------------------

        betaquot = betak/betakm1
        do i=1,imt
          do j=1,jmt
            s(i,j) = Zres(i,j) + betaquot * s(i,j)
          enddo
        enddo

!-----------------------------------------------------------------------
!       As(k)     = A * s(k)
!-----------------------------------------------------------------------

        if (npt .eq. 5) then
          call op5_vec(cf, s, As)
        else
          call op9_vec(cf, s, As)
        endif

        call border(As, bc_symm)

!-----------------------------------------------------------------------
!       If s=0 then the division for alpha(k) gives a float exception.
!       Assume convergence rate of 0.99 to extrapolate error.
!       Also assume alpha(k) ~ 1.
!-----------------------------------------------------------------------

        s_dot_As = dot2(s, As)
        if (abs(s_dot_As) .lt. abs(betak)*1.e-10) then
          smax = absmax(s)
          estimated_error = 100.0 * smax
          goto 101
        endif

!-----------------------------------------------------------------------
!       alpha(k)  = beta(k) / (s(k) * As(k))
!-----------------------------------------------------------------------

        alpha = betak / s_dot_As

!-----------------------------------------------------------------------
!       update values:
!       dpsi(k)   = dpsi(k-1) + alpha(k) * s(k)
!       res(k)    = res(k-1) - alpha(k) * As(k)
!-----------------------------------------------------------------------

        do i=1,imt
          do j=1,jmt
            dpsi (i,j) = dpsi(i,j) + alpha * s(i,j)
            res  (i,j) = res (i,j) - alpha * As(i,j)
          enddo
        enddo
        call avg_dist (res,
     &          imask, iperm, jperm, iofs, nisle, nippts)

        call border(res, bc_symm)

        smax = absmax(s)

!-----------------------------------------------------------------------
!       test for convergence
!       if (estimated_error) < epsilon) exit
!-----------------------------------------------------------------------

        step = abs(alpha) * smax
        if (k .eq. 1) then
          step1 = step
          estimated_error = step
          if (step .lt. epsilon) goto 101
        elseif (step .lt. epsilon) then
          cfactor = log(step/step1)
          convergence_rate = exp(cfactor/(k-1))
          estimated_error = step*convergence_rate/(1.0-convergence_rate)
          if (estimated_error .lt. epsilon) goto 101
        endif

        betakm1 = betak

      enddo

!-----------------------------------------------------------------------
!     end of iteration loop
!-----------------------------------------------------------------------

  101 continue
      if (k .gt. max_iterations) then
        cfactor = log(step/step1)
        convergence_rate = exp(cfactor/(k-1))
        estimated_error = step*convergence_rate/(1.0-convergence_rate)
        converged = .false.
      else
        if (diverging) then
          converged = .false.
        else
          converged = .true.
        endif
      endif

      iterations = k

!-----------------------------------------------------------------------
!     return the last increment of dpsi in the argument res
!-----------------------------------------------------------------------

      if (iterations .eq. 0) then
        do i=1,imt
          do j=1,jmt
            res(i,j) = Zres(i,j)
          enddo
        enddo
      else
        do i=1,imt
          do j=1,jmt
            res(i,j) = alpha * s(i,j)
          enddo
        enddo
      endif

      return
      end

!=======================================================================

!     M A T R I X   M O D U L E   F O R   C O N G R A D

!=======================================================================

      subroutine zero_vec (v)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt)

      do i=1,imt
        do j=1,jmt
          v(i,j) = 0.0
        enddo
      enddo

      return
      end

      subroutine add_vec (v,w,vpw)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt), w(imt,jmt), vpw(imt,jmt)

      do i=1,imt
        do j=1,jmt
          vpw(i,j) = v(i,j) + w(i,j)
        enddo
      enddo

      return
      end

      subroutine sub_vec (v,w,vmw)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt), w(imt,jmt), vmw(imt,jmt)

      do i=1,imt
        do j=1,jmt
          vmw(i,j) = v(i,j) - w(i,j)
        enddo
      enddo

      return
      end

      subroutine mult_vec(v,w,vtw)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt), w(imt,jmt), vtw(imt,jmt)

      do i=1,imt
        do j=1,jmt
          vtw(i,j) = v(i,j) * w(i,j)
        enddo
      enddo

      return
      end

      subroutine div_vec(v,w,vdw)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt), w(imt,jmt), vdw(imt,jmt)

      do i=1,imt
        do j=1,jmt
          if (w(i,j) .ne. 0) then
            vdw(i,j) = v(i,j) / w(i,j)
          else
            vdw(i,j) = 0.0
          endif
        enddo
      enddo

      return
      end

      subroutine scalar_vec (scalar,w,sw)

      implicit none

      integer i, j

      include "size.h"

      real w(imt,jmt), sw(imt,jmt), scalar

      do i=1,imt
        do j=1,jmt
          sw(i,j) = scalar * w(i,j)
        enddo
      enddo

      return
      end

      subroutine neg_vec (v)

      implicit none

      integer i, j

      include "size.h"

      real v(imt,jmt)

      do i=1,imt
        do j=1,jmt
          v(i,j) = -v(i,j)
        enddo
      enddo

      return
      end

      function dot2 (dp_vec, res_vec)

!     this dot product produces the correct answers because for
!     ocean perimeter subscripts, ij=isle, the value on a
!     type(dpsi_type) vector, dp_vec(isle)=dp_vec(i,j), i.e., the true
!     value is replicated, and for a type(res_type) vector,
!     res_vec(isle) = sum (res_vec(i,j)), i.e., the true value is the
!     accumulation of the distributed values.

      implicit none

      integer j, i

      real dot2

      include "size.h"

      real dp_vec(imt,jmt), res_vec(imt,jmt), rowsum (jmt)

      do j=2,jmt-1
        rowsum(j) = 0.0
        do i=2,imt-1
          rowsum(j) = rowsum(j) + dp_vec(i,j) * res_vec(i,j)
        enddo
      enddo

      dot2 = 0.0
      do j=2,jmt-1
        dot2 = dot2 + rowsum(j)
      enddo

      return
      end

      subroutine op5_vec(cf, dpsi, res)

!                       res = A * dpsi

!     this subroutine does not collect the terms of the true value
!     of res(isle) = sum (res(i,j)).  the contributions to the sum
!     remain distributed among the T cells (i,j) that form the
!     ocean perimeter of land_mass(isle).

!     at present, borders are not computed [i=1 or imt] [j=1 or jmt]

      implicit none

      integer i, j

      include "size.h"

      real cf(imt,jmt,-1:1,-1:1), dpsi(imt,jmt), res(imt,jmt)

      do j=2,jmt-1
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i,j)   +
     &               cf(i,j, 0, 1) * dpsi(i,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j)
        enddo
      enddo

      return
      end

      subroutine op9_vec(cf, dpsi, res)

!                       res = A * dpsi

!     this subroutine does not collect the terms of the true value
!     of res(isle) = sum (res(i,j)).  the contributions to the sum
!     remain distributed among the T cells (i,j) that form the
!     ocean perimeter of land_mass(isle).

!     at present, borders are not computed [i=1 or imt] [j=1 or jmt]

      implicit none

      integer i, j

      include "size.h"

      real cf(imt,jmt,-1:1,-1:1), dpsi(imt,jmt), res(imt,jmt)

      do j=2,jmt-1
        do i=2,imt-1
          res(i,j) = cf(i,j, 0, 0) * dpsi(i  ,j  ) +
     &               cf(i,j, 0, 1) * dpsi(i  ,j+1) +
     &               cf(i,j, 0,-1) * dpsi(i  ,j-1) +
     &               cf(i,j, 1, 0) * dpsi(i+1,j  ) +
     &               cf(i,j,-1, 0) * dpsi(i-1,j  ) +
     &               cf(i,j, 1, 1) * dpsi(i+1,j+1) +
     &               cf(i,j,-1, 1) * dpsi(i-1,j+1) +
     &               cf(i,j, 1,-1) * dpsi(i+1,j-1) +
     &               cf(i,j,-1,-1) * dpsi(i-1,j-1)
        enddo
      enddo

      return
      end

      subroutine subset (a, b, nerror)

!     verifies that the set of subscripts for which a(i,j) .ne. 0.0
!     is a subset of the set of subscripts for which b(i,j) .ne. 0.0

      implicit none

      integer i, j, nerror

      include "size.h"

      real a(imt,jmt), b(imt,jmt)

      nerror = 0
      do i=2,imt-1
        do j=2,jmt-1
          if (a(i,j) .ne. 0.0 .and. b(i,j) .eq. 0.0) then
            nerror = nerror + 1
            print '(a,i3,a,i3,a,a)', '(',i,',',j,')'
     &                           ,' forcing is reset to zero'
!           set forcing (i.e., a(i,j)) to zero
            a(i,j) = 0.0
          endif
        enddo
      enddo

      return
      end

      subroutine inv_op(Z, res, Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)

!     apply and approximate inverse Z or the operator A

!     res is type(res_type), i.e., perimeter values res(isle)
!         are the sum of the distributed contributions res(i,j)
!     Zres is type(dpsi_type), i.e., perimeter values Zres(isle)
!         must be replicated at each perimeter point Zres(i,j)

!     borders  of Zres [i=1 or imt] [j=1 or jmt] must be defined
!     and must satisfy cyclic and/or symmetry, if defined.

!     currently, Z is diagonal:  Z(ij) = 1/A(ij)
!     and is stored in type(dpsi_type) format, i.e., Z(isle) is
!     replicated and stored in each Z(i,j) in the perimeter of
!     land_mass(isle).

      implicit none

      integer i, j, nisle

      include "size.h"

      logical imask(-mnisle:mnisle)

      integer iperm(maxipp), jperm(maxipp), iofs(mnisle), nippts(mnisle)

      real Z(imt,jmt), res(imt,jmt), Zres(imt,jmt)

      do i=1,imt
        do j=1,jmt
          Zres(i,j) = Z(i,j) * res(i,j)
        enddo
      enddo

!     sum contributions to Zres(isle)
!     distribute Zres(isle) to all perimeter points

      call sum_dist (Zres,
     &        imask, iperm, jperm, iofs, nisle, nippts)

      return
      end

      function absvecmax(res, imax, jmax)

      implicit none

      integer i, j, imax, jmax

      real absvecmax

      include "size.h"

      real res(imt,jmt)

      absvecmax = 0.0
      do i=2,imt-1
        do j=2,jmt-1
          if (abs(res(i,j)) .gt. absvecmax) then
            absvecmax = abs(res(i,j))
            imax = i
            jmax = j
          endif
        enddo
      enddo
      return
      end

      function absmax (f)

      implicit none

      integer i, j

      real amax, absmax

      include "size.h"

      real f(imt,jmt)

      amax = 0.0
      do i=1,imt
        do j=1,jmt
          amax = max(amax, abs(f(i,j)))
        enddo
      enddo
      absmax = amax
      return
      end

      function absmin (f)

      implicit none

      integer i, j

      real amin, absmin

      include "size.h"

      real f(imt,jmt)

      amin = 1.0e37
      do i=1,imt
        do j=1,jmt
          if (f(i,j) .ne. 0 .and. abs(f(i,j)) .lt. amin) then
            amin = abs(f(i,j))
          endif
        enddo
      enddo
      absmin = amin
      return
      end

      subroutine make_inv (cf, Z,
     &              imask, iperm, jperm, iofs, nisle, nippts)

!     construct an approximate inverse Z to A

!     Z will be diagonal:  Z(ij) = 1/A(ij)
!     and values for ocean perimeter entries Z(isle) will be replicated
!     at all T cells Z(i,j) in the ocean perimeter of land_mass(isle).

!     T cells (i,j) for which there is no diagonal coefficient
!     i.e., A(ij)=A(i,j)=0, are masked off by assigning Z(i,j)=0.
!     there are effectively no equations and no variables dpsi(i,j)
!     at these points.

      implicit none

      integer i, j, nisle, isle, n

      include "size.h"

      integer iperm(maxipp),jperm(maxipp)
      integer iofs(mnisle), nippts(mnisle)

      logical imask(-mnisle:mnisle)

      real cf(imt,jmt,-1:1,-1:1), Z(imt,jmt)

!     copy diagonal coefficients of A to Z

      do i=2,imt-1
        do j=2,jmt-1
          Z(i,j) = cf(i,j,0,0)
        enddo
      enddo

!     for each land_mass(isle),
!     sum the contributions to cfdiag(isle)=A(isle,isle)
!     now stored in Z(i,j) at ocean perimeter T cells and replicate
!     the sum in all Z(i,j) for which (i,j) is in ocean perimeter
!     of land_mass(isle).

      call sum_dist (Z,
     &        imask, iperm, jperm, iofs, nisle, nippts)

!     now invert Z

      do i=2,imt-1
        do j=2,jmt-1
          if (Z(i,j) .ne. 0.0) then
            Z(i,j) = 1/Z(i,j)
          else
            Z(i,j) = 0.0
          endif
        enddo
      enddo

!     make inverse zero on island perimeters that are not integrated

      do isle=1,nisle
        if (.not. imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Z(i,j) = 0.0
          enddo
        endif
      enddo

      return
      end

      subroutine sum_dist (Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)

!     sum contributions to Zres(isle)
!     distribute Zres(isle) to all perimeter points

!     this subroutine converts a type(res_type) vector with
!     distributed contributions to perimeter values
!        Zres(isle) = sum (Zres(i,j))
!     into a type (dpsi_type) vector with replicated values
!     for land_mass perimeters
!        Zres(isle) = Zres(i,j)
!     for all (i,j) in the ocean perimeter of land_mass(isle).

      implicit none

      integer i, j, isle, nisle, n

      include "size.h"

      integer iperm(maxipp),jperm(maxipp)
      integer iofs(mnisle), nippts(mnisle)

      logical imask(-mnisle:mnisle)

      real Zres(imt,jmt), Zresisle(mnisle)

!     sum contributions to Zres(isle)

      do isle=1,nisle
        if (imask(isle)) then
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          enddo
        endif
      enddo

!     distribute Zres(isle) to all perimeter points

      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)
          enddo
        endif
      enddo

      return
      end

      subroutine avg_dist (Zres,
     &              imask, iperm, jperm, iofs, nisle, nippts)

!     avg contributions to Zres(isle)
!     distribute Zres(isle) to all perimeter points

!     this subroutine converts a type(res_type) vector with
!     distributed contributions to perimeter values
!        Zres(isle) = avg (Zres(i,j))
!     into a type (dpsi_type) vector with replicated values
!     for land_mass perimeters
!        Zres(isle) = Zres(i,j)
!     for all (i,j) in the ocean perimeter of land_mass(isle).

      implicit none

      integer i, j, isle, nisle, n

      include "size.h"

      integer iperm(maxipp),jperm(maxipp)
      integer iofs (mnisle), nippts(mnisle)

      logical imask(-mnisle:mnisle)

      real Zres(imt,jmt), Zresisle(mnisle)

!     avg contributions to Zres(isle)

      do isle=1,nisle
        if (imask(isle)) then
!         print *,'isle=',isle,' nisle=',nisle
          Zresisle(isle) = 0.0
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zresisle(isle) = Zresisle(isle) + Zres(i,j)
          enddo
        endif
      enddo

!     distribute Zres(isle) to all perimeter points

      do isle=1,nisle
        if (imask(isle)) then
          do n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            Zres(i,j) = Zresisle(isle)/nippts(isle)
          enddo
        endif
      enddo
#endif

      return
      end
