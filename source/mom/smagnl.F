      subroutine smagnli (hmixset)

#if defined O_mom && defined O_smagnlmix && !defined O_consthmix
!=======================================================================
!     Initialization for the Smagorinsky nonlinear horizontal viscosity
!     as described in Rosati & Miyakoda (jpo,vol 18,#11,1988)
!     see Smagorinsky 1963, Mon Wea Rev, 91, 99-164.
!     Also see Deardorff 1973 J. Fluid Eng. Sep., 429-438.

!     input:
!       hmixset = logical to determine if a vertical mixing scheme was
!                 chosen
!     output:
!       hmixset = logical to determine if a vertical mixing scheme was
!                 chosen
!=======================================================================

      implicit none

      integer io, j, i, k

      logical hmixset

      real c0

      include "size.h"
      include "hmixc.h"
      include "iounit.h"
      include "switch.h"
      include "stdunits.h"

      namelist /smagnl/ diff_c_back

      write (stdout,'(/,20x,a,/)')
     & 'S M A G M I X    I N I T I A L I Z A T I O N'

!-----------------------------------------------------------------------
!     initialize variables (all mixing units are cm**2/sec.)
!-----------------------------------------------------------------------

      c0          = 0.0
      diff_c_back = c0

!-----------------------------------------------------------------------
!     provide for namelist over-ride of above settings + documentation
!-----------------------------------------------------------------------

      call getunit (io, 'control.in'
     &,               'formatted sequential rewind')
      read (io,smagnl,end=100)
100   continue
      write (stdout,smagnl)
      call relunit (io)

!-----------------------------------------------------------------------
!     check for problems
!-----------------------------------------------------------------------

      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "smagnlmix"  cannot be enabled because another   '
     &,'            horizontal mixing scheme has been enabled    '
        stop '=> smagnli'
      else
        hmixset = .true.
      endif

!-----------------------------------------------------------------------
!     initialize arrays
!-----------------------------------------------------------------------

      do j=1,jemw
        do i=1,imt
          do k=1,km
            strain(i,k,j,1) = c0
            strain(i,k,j,2) = c0
            am_lambda(i,k,j)= c0
            am_phi(i,k,j)   = c0
            visc_cnu(i,k,j) = c0
            diff_cnt(i,k,j) = c0
          enddo
        enddo
      enddo
      do j=jsmw,jemw
        do i=1,imt
          do k=1,km
            visc_ceu(i,k,j) = c0
            diff_cet(i,k,j) = c0
          enddo
        enddo
      enddo
      return
      end

      subroutine smagnlc (joff, js, je, is, ie)

      implicit none

      integer istrt, is, iend, ie, jstrt, js, jend, je, j, jrow, joff
      integer k, i, jsc

      real p25, p5, c0, c1, c2, sqrt2r, c14, f1, f2, cphi, deform

      include "size.h"
      include "param.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "switch.h"
      include "stdunits.h"

      real clam(imt)

!=======================================================================
!     compute tension and shearing rates of strain, total deformation
!     and diffusion coefficients for the Smagorinsky nonlinear
!     viscosity as described in Rosati & Miyakoda (jpo,vol 18,#11,1988).

!     input:

!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!      is     = starting longitude index for calculations
!      ie     = ending longitude index for calculations

!     output:

!      strain  = tension (1) and shearing (2) rates of strain
!      defor  = total deformation

!=======================================================================

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      p25    = 0.25
      p5     = 0.5
      c0     = 0.0
      c1     = 1.0
      c2     = 2.0
      sqrt2r = c1/sqrt(c2)
      c14    = 0.14

      istrt  = max(2,is)
      iend   = min(imt-1,ie)
      jstrt  = max(1,js-1)
      jend   = je-1

!-----------------------------------------------------------------------
!     compute rates of strain: tension and shear on north face of
!     "u" cell. Rosati & Miyakoda (jpo,vol 18,#11,1988) Eqn 2.24 & 2.25
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = j + joff
        f1   = p5*cstr(jrow+1)
        f2   = dytr(jrow+1)*cst(jrow+1)
        do k=1,km
          do i=istrt,iend
            strain(i,k,j,1) = ((u(i+1,k,j,1,taum1)+u(i+1,k,j+1,1,taum1))
     &                      - (u(i-1,k,j,1,taum1)+u(i-1,k,j+1,1,taum1)))
     &                          *f1*dxu2r(i)
     &                      - (csur(jrow+1)*u(i,k,j+1,2,taum1) -
     &                         csur(jrow  )*u(i,k,j,2,taum1))*f2
            strain(i,k,j,2) = ((u(i+1,k,j,2,taum1)+u(i+1,k,j+1,2,taum1))
     &                      - (u(i-1,k,j,2,taum1)+u(i-1,k,j+1,2,taum1)))
     &                          *f1*dxu2r(i)
     &                      + (csur(jrow+1)*u(i,k,j+1,1,taum1) -
     &                         csur(jrow  )*u(i,k,j,1,taum1))*f2
          enddo
        enddo
        call setbcx (strain(1,1,j,1), imt, km)
        call setbcx (strain(1,1,j,2), imt, km)
      enddo

!-----------------------------------------------------------------------
!     compute effective anisentropic wavenumber of diffusing turbulence
!     for effective 2d isentropic wavenumber set coeffs
!     to: (c14*sqrt(csu(jrow)*dxu(i)*dyu(jrow)))**2*sqrt2r
!     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.28 & 2.29
!     (note: Eqn 2.28 should have m**-1)
!     compute total deformation + viscosity coefficients on north face
!     of "u" cell
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow     = joff + j
        cphi  = (c14*dyu(jrow))**2*sqrt2r
        do i=istrt,iend
          clam(i) = (c14*csu(jrow)*dxu(i))**2*sqrt2r
        enddo
        do k=1,km
          do i=istrt,iend
            deform     = sqrt(c2*(strain(i,k,j,1)**2
     &                                + strain(i,k,j,2)**2))
            am_lambda(i,k,j) = clam(i)*deform
            am_phi(i,k,j)    = cphi*deform
          enddo
        enddo
        call setbcx (am_lambda(1,1,j), imt, km)
        call setbcx (am_phi(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     set j index for calculating diffusive coefficient on north face
!     of cells.
!-----------------------------------------------------------------------

      jsc = max(jsmw,js-1)

# if defined O_matrix_sections
      if (prxzts .and. eots) then

!-----------------------------------------------------------------------
!       calculate the mixing coefficients for momentum and diffusion
!       Rosati & Miyakoda(jpo,vol 18,#11,1988)
!       Eqn: 2.26 & 2.27  mixing coeff for momentum
!       Both "visc_cnu" and "visc_ceu" are purely diagnostic here since
!       momentum flux "diff_fe" and "diff_fn" is calculated in terms of
!       shearing and tension rates of strain in subroutine "smagnlm"
!-----------------------------------------------------------------------

!       viscosity coeff on north face of "u" cells

        do j=jstrt,jend
          do k=1,km
            do i=istrt,iend
              visc_cnu(i,k,j) = am_phi(i,k,j)
            enddo
          enddo
          call setbcx (visc_cnu(1,1,j), imt, km)
        enddo

!       viscosity coeff on east face of "u" cells

        do j=jsc,jend
          do k=1,km
            do i=istrt-1,iend
              visc_ceu(i,k,j) = p25*(am_lambda(i,k,j)
     &                             + am_lambda(i+1,k,j)
     &                             + am_lambda(i,k,j-1)
     &                             + am_lambda(i+1,k,j-1))
            enddo
          enddo
          call setbcx (visc_ceu(1,1,j), imt, km)
        enddo
      endif
# endif

!-----------------------------------------------------------------------
!     calculate the diffusion coefficients for "t" cells
!     Rosati & Miyakoda(jpo,vol 18,#11,1988)
!     Eqn: 2.35 & 2.36  mixing coeff for tracers
!-----------------------------------------------------------------------

!     diffusion coeff on east face of "t" cells

      do j=jsc,jend
        do k=1,km
          do i=istrt-1,iend
            diff_cet(i,k,j) = am_lambda(i,k,j-1) + diff_c_back
          enddo
        enddo
        call setbcx (diff_cet(1,1,j), imt, km)
      enddo

!     diffusion coeff on north face of "t" cells. Index jsc is used
!     because "j-1" is needed. The southern wall is taken care of by
!     masking in the diffusion operator.

      do j=jsc,jend
        jrow = j + joff
        do k=1,km
          do i=istrt,iend
            diff_cnt(i,k,j) = p25*(am_phi(i,k,j) + am_phi(i-1,k,j)
     &                           + am_phi(i,k,j-1) + am_phi(i-1,k,j-1))
     &                           + diff_c_back
          enddo
        enddo
        call setbcx (diff_cnt(1,1,j), imt, km)
      enddo

# if defined O_matrix_sections
      if (prxzts .and. eots) then
        call diagnl (joff, jsc, jend)
      endif
# endif

      return
      end

      subroutine smagnlm (joff, js, je, is, ie, n)

      implicit none

      integer istrt, is, iend, ie, jstrt, js, jend, je, j, k, i, n
      integer jrow, joff

      real p25, p5, c0, rmsq, cstsq, cstsqm, f1

      include "size.h"
      include "param.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "stdunits.h"

!=======================================================================
!     compute diffusive flux across north and east face of "u" cells
!     for velocity component "n"

!     input:

!      joff   = offset between row j in the MW and latitude jrow
!      js     = starting row for calculations
!      je     = ending row for calculations
!      is     = starting longitude index for calculations
!      ie     = ending longitude index for calculations
!      n      = velocity component: 1 = "u", 2 = "v"
!        + quantities from subroutine "dform"

!     output:

!      diff_fe  = diffusive flux across eastern face of "u" cells
!      diff_fn  = diffusive flux across northern face of "u" cells

!=======================================================================

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      p25    = 0.25
      p5     = 0.5
      c0     = 0.0

      istrt  = max(2,is)
      iend   = min(imt-1,ie)
      jstrt  = max(2,js)
      jend   = min(jemw,je)

!-----------------------------------------------------------------------
!     compute zonal flux components of the stress tensor on the eastern
!     face of the "u" cell
!     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.18
!-----------------------------------------------------------------------

      do j=jstrt,jend
        do k=1,km
          do i=istrt-1,ie
            diff_fe(i,k,j) = p25*(
     &                am_lambda(i,k,j)*strain(i,k,j,n)
     &              + am_lambda(i,k,j-1)*strain(i,k,j-1,n)
     &              + am_lambda(i+1,k,j)*strain(i+1,k,j,n)
     &              + am_lambda(i+1,k,j-1)*strain(i+1,k,j-1,n))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute meridional flux components of the stress tensor on the
!     north face of the "u" cell
!     Rosati & Miyakoda(jpo,vol 18,#11,1988). Eqn: 2.19
!-----------------------------------------------------------------------

      if (n .eq. 1) then

!       northward flux term for zonal memoentum eqn is zero

        do j=jstrt-1,jend
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = c0
            enddo
          enddo
        enddo

!       compute second term (which is not a flux term) when working on
!       the zonal momentum eqn.

        do j=jstrt,jend
          jrow   = joff + j
          rmsq   = csur(jrow)**2*dyur(jrow)
          cstsq  = cst(jrow+1)**2
          cstsqm = cst(jrow)**2
          do k=1,km
            do i=istrt,iend
              smag_metric(i,k,j) = rmsq*
     &                     (am_phi(i,k,j)*strain(i,k,j,2)*cstsq-
     &                      am_phi(i,k,j-1)*strain(i,k,j-1,2)*cstsqm)
            enddo
          enddo
        enddo
      elseif (n .eq. 2) then

!       northward flux for meridional momentum equation.

        do j=jstrt-1,jend
          jrow = joff + j
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =-cst(jrow+1)*am_phi(i,k,j)*strain(i,k,j,1)
            enddo
          enddo
        enddo

!       compute third term (which is not a flux term) when working on
!       the meridional component of the momentum eqn.

        do j=jstrt,jend
          jrow = joff + j
          f1   = csur(jrow)*sine(jrow)*p5/radius
          do k=1,km
            do i=istrt,iend
              smag_metric(i,k,j) = f1*(am_lambda(i,k,j)*strain(i,k,j,1)
     &                          +  am_lambda(i,k,j-1)*strain(i,k,j-1,1))
            enddo
          enddo
        enddo
      endif

      return
      end

# if defined O_matrix_sections
      subroutine diagnl (joff, js, je)

      implicit none

      integer j, js, je, jrow, joff, jlat, jj, indp, is, ie, ks, ke

      real reltim, fx, scl

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cprnts.h"
      include "hmixc.h"
      include "iounit.h"
      include "switch.h"
      include "tmngr.h"

!=======================================================================
!     Show some diagnostics

!     input:

!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!=======================================================================

      do j=js,je
        jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2

!           write out the diffusion coeffs for tracers

            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cet', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cet(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cet ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cet(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cet(1,1,j), imt*km)
              call relunit (io)
            endif

            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cnt', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cnt(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cnt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cnt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cnt(1,1,j), imt*km)
              call relunit (io)
            endif

!           write out the viscosity coeffs for momentum

            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_ceu', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_ceu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => visc_ceu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_ceu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_ceu(1,1,j), imt*km)
              call relunit (io)
            endif

            scl = 1.e7
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cnu', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cnu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                    'unformatted sequential append ieee')
              write (stdout,*) ' => visc_cnu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cnu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cnu(1,1,j), imt*km)
              call relunit (io)
            endif
          endif
        enddo
      enddo

9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
# endif
#endif

      return
      end
