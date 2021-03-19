      subroutine clinic (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     compute internal mode velocity components for rows js through je
!     in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer istrt, iend
      integer i, k, j, jrow, n, js, je, limit, joff, kb, is, ie

      real adv_ux, adv_uy, adv_uz, adv_metric, diff_ux, diff_uz, fx
      real diff_uy, diff_metric, coriolis, aprime, grav_rho0r, fxa
      real fxb, t1, t2, ambi_csur, del2, ambi_cst_dytr, detmr, unep

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
# if defined O_neptune
      include "cnep.h"
# endif
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "hmixc.h"
      include "emode.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"
      include "fdifm.h"

      parameter (istrt=2, iend=imt-1)

      real tempik(imt,km,jsmw:jmw)
      real baru(imt,jsmw:jemw,2)

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: these are currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

# if defined O_biharmonic
      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
# else
      do j=js,je
# endif
        jrow = j + joff
# if defined O_anisotropic_viscosity
        do k=1,km
           do i=istrt-1,iend
              csudxur(i,j)  = csur(jrow)*dxur(i)
              csudxu2r(i,j) = csur(jrow)*dxur(i)*p5
              am_csudxtr(i,k,j) = visc_ceu(i,k,j)*csur(jrow)*dxtr(i+1)
           enddo
        enddo
# else
        do i=istrt-1,iend
          csudxur(i,j)  = csur(jrow)*dxur(i)
          csudxu2r(i,j) = csur(jrow)*dxur(i)*p5
#  if defined O_consthmix
          am_csudxtr(i,j)  = am*csur(jrow)*dxtr(i+1)
#  endif
        enddo
# endif
      enddo

!-----------------------------------------------------------------------
!     construct the hydrostatic pressure gradients: 1 = dp/dx; 2 = dp/dy
!-----------------------------------------------------------------------

!     compute horizontal pressure gradient at the first level

# if defined O_pressure_gradient_average

!     construct density as in Brown and Campana (1978) (see manual)

      call state (t(1,1,1,1,taum1), t(1,1,1,2,taum1), rhotaum1(1,1,jsmw)
     &,           js, je+1, istrt-1, iend+1)
      call state (t(1,1,1,1,taup1), t(1,1,1,2,taup1), rhotaup1(1,1,jsmw)
     &,           js, je+1, istrt-1, iend+1)
      aprime = 0.25
      do j=js,je+1
        do k=1,km
          do i=1,imt
            rhotilde(i,k,j) = aprime*(rhotaum1(i,k,j) + rhotaup1(i,k,j))
     &       + (1.0-2.0*aprime)*rho(i,k,j)
          enddo
        enddo
      enddo
# endif

      grav_rho0r = grav*rho0r
      do j=js,je
        jrow = j + joff
        fxa  = grav_rho0r*dzw(0)*csur(jrow)
        fxb  = grav_rho0r*dzw(0)*dyu2r(jrow)
        do i=istrt-1,iend
# if defined O_pressure_gradient_average
          t1              = rhotilde(i+1,1,j+1) - rhotilde(i  ,1,j)
          t2              = rhotilde(i  ,1,j+1) - rhotilde(i+1,1,j)
# else
          t1              = rho(i+1,1,j+1) - rho(i  ,1,j)
          t2              = rho(i  ,1,j+1) - rho(i+1,1,j)
# endif
          grad_p(i,1,j,1) = (t1-t2)*fxa*dxu2r(i)
          grad_p(i,1,j,2) = (t1+t2)*fxb
        enddo
      enddo

!     compute the change in pressure gradient between levels

      do j=js,je+1
        do k=2,km
          do i=istrt-1,iend+1
# if defined O_pressure_gradient_average
            tempik(i,k,j) = rhotilde(i,k-1,j) + rhotilde(i,k,j)
# else
            tempik(i,k,j) = rho(i,k-1,j) + rho(i,k,j)
# endif
          enddo
        enddo
      enddo

      do j=js,je
        jrow = j + joff
        fxa = grav_rho0r*csur(jrow)*p5
        fxb = grav_rho0r*dyu4r(jrow)
        do k=2,km
          do i=istrt-1,iend
            t1              = tempik(i+1,k,j+1) - tempik(i  ,k,j)
            t2              = tempik(i  ,k,j+1) - tempik(i+1,k,j)
            grad_p(i,k,j,1) = fxa*(t1-t2)*dzw(k-1)*dxu2r(i)
            grad_p(i,k,j,2) = fxb*(t1+t2)*dzw(k-1)
          enddo
        enddo
      enddo

!     integrate downward from the first level

      do j=js,je
        do k=1,kmm1
          do i=istrt-1,iend
            grad_p(i,k+1,j,1) = grad_p(i,k,j,1) + grad_p(i,k+1,j,1)
            grad_p(i,k+1,j,2) = grad_p(i,k,j,2) + grad_p(i,k+1,j,2)
          enddo
        enddo
      enddo

      do j=js,je
        call setbcx (grad_p(1,1,j,1), imt, km)
        call setbcx (grad_p(1,1,j,2), imt, km)
      enddo

!-----------------------------------------------------------------------
!     solve for one component of velocity at a time
!     n = 1 => zonal component
!     n = 2 => meridional component
!-----------------------------------------------------------------------

      do n=1,2

# if !defined O_linearized_advection

!-----------------------------------------------------------------------
!       calculate 2*advective flux (for speed) across east face of
!       "u" cells.
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              adv_fe(i,k,j) = adv_veu(i,k,j)*(u(i,  k,j,n,tau) +
     &                                        u(i+1,k,j,n,tau))
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       2*advective flux across northern face of "u" cells is built
!       into ADV_Uy. (It's done this way for performance issues)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       diffusive flux across east face of "u" cell
!       diffusive flux across north face of "u" cell
!-----------------------------------------------------------------------

#  if defined O_consthmix && !defined O_biharmonic

!       build diffusive flux on eastern face of "u" cells

        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend
#   if defined O_anisotropic_viscosity
              diff_fe(i,k,j) = am_csudxtr(i,k,j)*
#   else
              diff_fe(i,k,j) = am_csudxtr(i,j)*
#   endif
     &                        (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1)
#   if defined O_neptune
     &                        - unep(i+1,jrow,n)*umask(i+1,k,j) +
     &                          unep(i,jrow,n)*umask(i,k,j)
#   endif
     &                        )
            enddo
          enddo
        enddo

!       diffusive flux on northern face of "u" cells is built
!       into DIFF_Uy
#  endif
# endif
# if defined O_consthmix && defined O_biharmonic

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "u" cell
!       diffusive flux across northern face of "u" cell
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          ambi_csur = visc_ceu*csur(jrow)
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = ambi_csur*dxtr(i+1)*
     &                        (del2(i+1,k,j,n) - del2(i,k,j,n))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          ambi_cst_dytr = visc_cnu*cst(jrow+1)*dytr(jrow+1)
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = ambi_cst_dytr*
     &                        (del2(i,k,j+1,n) - del2(i,k,j,n))
            enddo
          enddo
        enddo
# endif
# if defined O_smagnlmix && !defined O_consthmix

!       calculate diffusive flux on eastern and northern faces of
!      "u" cells

        call smagnlm (joff, js, je, istrt, iend, n)
# endif

# if !defined O_linearized_advection

!-----------------------------------------------------------------------
!       calculate 2*advective flux (for speed) on bottom face of
!       "u" cell. also diffusive flux on bottom face of "u" cell
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,kmm1
            do i=istrt,iend
              adv_fb(i,k,j)  = adv_vbu(i,k,j)*(u(i,k,  j,n,tau) +
     &                                        u(i,k+1,j,n,tau))
              diff_fb(i,k,j) = visc_cbu(i,k,j)*dzwr(k)*
     &                         (u(i,k,j,n,taum1) - u(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo
# endif

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "u" cells for mixing
!       and advection.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            kb              = kmu(i,jrow)
            diff_fb(i,0,j)  = smf(i,j,n)
            diff_fb(i,kb,j) = bmf(i,j,n)
            adv_fb(i,0,j)   = adv_vbu(i,0,j)*(u(i,1,j,n,tau) +
     &                                        u(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbu(i,km,j)*u(i,km,j,n,tau)
          enddo
        enddo

# if defined O_source_term || defined O_mobi || defined O_carbon_14

!-----------------------------------------------------------------------
!       set source term for "u" cell
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km
            do i=istrt,iend
              source(i,k,j) = c0
            enddo
          enddo
        enddo
# endif

!-----------------------------------------------------------------------
!       solve for the internal mode part of du/dt at center of
!       "u" cells by neglecting the surface pressure gradients. use
!       statement functions to represent each component of the
!       calculation.
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) =
     &          (DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j,jrow,n) + DIFF_Uz(i,k,j)
     &          + DIFF_metric(i,k,j,jrow,n)
# if !defined O_linearized_advection
     &          - ADV_Ux(i,k,j) - ADV_Uy(i,k,j,jrow,n) - ADV_Uz(i,k,j)
     &          + ADV_metric(i,k,j,jrow,n)
# endif
     &          - grad_p(i,k,j,n) + CORIOLIS(i,k,j,jrow,n)
# if defined O_source_term || defined O_mobi || defined O_carbon_14
     &          + source(i,k,j)
# endif
     &          )*umask(i,k,j)
            enddo
          enddo
        enddo
# if defined O_implicitvmix

!-----------------------------------------------------------------------
!       add in du/dt component due to implicit vertical diffusion
!-----------------------------------------------------------------------

        call ivdifu (joff, js, je, istrt, iend, n)
# endif

!-----------------------------------------------------------------------
!       construct diagnostics associated with velocity component "n"
!-----------------------------------------------------------------------

        call diagc1 (joff, js, je, istrt, iend, n)

!-----------------------------------------------------------------------
!       construct the vertical average of du/dt for forcing
!       the barotropic equation
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = c0
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km
            fx = dzt(k)
            do i=istrt,iend
              zu(i,jrow,n) = zu(i,jrow,n) + u(i,k,j,n,taup1)*fx
            enddo
          enddo
        enddo

        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            zu(i,jrow,n) = zu(i,jrow,n)*hr(i,jrow)
          enddo
        enddo
# if defined O_symmetry
        do j=js,je
          jrow = j + joff
          if (jrow .eq. jmtm1 .and. n .eq. 2) then
            do i=istrt,iend
              zu(i,jrow,2)   =  c0
              zu(i,jrow+1,2) = -zu(i,jrow-1,2)
              zu(i,jrow+1,1) =  zu(i,jrow-1,1)
            enddo
          endif
        enddo
# endif

!-----------------------------------------------------------------------
!       end of velocity component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     compute "tau+1" velocities accounting for implicit part of the
!     coriolis term if treated implicitly. velocities are in error by an
!     arbitrary constant related to neglecting the unknown surface
!     pressure gradients
!-----------------------------------------------------------------------

# if defined O_damp_inertial_oscillation
      do j=js,je
        jrow  = j + joff
        do k=1,km
          do i=istrt,iend
            fx    = c2dtuv*acor*cori(i,jrow,1)
            detmr = c1/(c1 + fx*fx)
            t1 = (u(i,k,j,1,taup1) + fx*u(i,k,j,2,taup1))*detmr
            t2 = (u(i,k,j,2,taup1) - fx*u(i,k,j,1,taup1))*detmr
            u(i,k,j,1,taup1) = u(i,k,j,1,taum1) + c2dtuv*t1
            u(i,k,j,2,taup1) = u(i,k,j,2,taum1) + c2dtuv*t2
          enddo
        enddo
      enddo
# else
      do n=1,2
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taum1)
     &                            + c2dtuv*u(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      enddo
# endif

!-----------------------------------------------------------------------
!     subtract incorrect vertical means (related to ignoring horizontal
!     gradients of the surface pressure) to get pure internal modes.
!-----------------------------------------------------------------------

      do n=1,2
        do j=js,je
          do i=istrt,iend
            baru(i,j,n) = c0
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              baru(i,j,n) = baru(i,j,n) + u(i,k,j,n,taup1)*dzt(k)
            enddo
          enddo
        enddo
        do j=js,je
          jrow  = j + joff
          do i=istrt,iend
            baru(i,j,n) = baru(i,j,n)*hr(i,jrow)
          enddo
        enddo
        do j=js,je
          do k=1,km
            do i=istrt,iend
              u(i,k,j,n,taup1) = u(i,k,j,n,taup1)
     &                          - umask(i,k,j)*baru(i,j,n)
            enddo
          enddo
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics involving internal mode velocity at "tau+1"
!-----------------------------------------------------------------------

      call diagc2 (joff, js, je, is, ie)

# if defined O_fourfil || defined O_firfil

!-----------------------------------------------------------------------
!     filter velocity components at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filuv (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in clinic'
        stop '=>clinic'
      endif
# endif
      do j=js,je
        call setbcx (u(1,1,j,1,taup1), imt, km)
        call setbcx (u(1,1,j,2,taup1), imt, km)
      enddo

# if defined O_ice_evp

!-----------------------------------------------------------------------
!     if needed, construct the Ice S.B.C.(surface boundary conditions)
!     averaged over this segment
!-----------------------------------------------------------------------

      call isbcu (joff, js, je, istrt, iend, igu, igv)
      call asbcu (joff, js, je, istrt, iend, isu, isv)
# endif

      return
      end

      subroutine diagc1 (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     construct diagnostics which don`t require internal mode velocity
!     at "tau+1" for each velocity component "n"

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = (1,2) = (u,v) velocity component
!-----------------------------------------------------------------------

      implicit none

      integer n, j, js, je, jrow, joff, k, i, is, ie

      real dudx, ce, dudy, cn, dudz, cb, fx, weight

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "diag.h"
      include "diaga.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

      real temp(imt,km)

# if defined O_save_mixing_coeff

!-----------------------------------------------------------------------
!     diagnostic: estimate mixing coefficients on east, north, and
!                 bottom face of U cells from the flux
!-----------------------------------------------------------------------

      if (cmixts .and. n .eq. 1 .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dudx = (u(i+1,k,j,1,taum1)-u(i,k,j,1,taum1))
     &                *csur(jrow)*dxtr(i+1) + epsln
              ce(i,k,j,1) = diff_fe(i,k,j)/dudx
#  if !defined O_consthmix || defined O_biharmonic
              dudy = (u(i,k,j+1,1,taum1)-u(i,k,j,1,taum1))
     &                *dytr(jrow+1) + epsln
              cn(i,k,j,1) = cstr(jrow+1)*diff_fn(i,k,j)/dudy
#  else
              cn(i,k,j,1) = cstr(jrow+1)*am
#  endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dudz = (u(i,k,j,1,taum1)-u(i,k+1,j,1,taum1))
     &                *dzwr(k) + epsln
              cb(i,k,j,1) = diff_fb(i,k,j)/dudz
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,1) = 0.0
          enddo
        enddo

        do j=js,je
          call setbcx (ce(1,1,j,1), imt, km)
          call setbcx (cn(1,1,j,1), imt, km)
          call setbcx (cb(1,1,j,1), imt, km)
        enddo
      endif
# endif

# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     diagnostic: accumulate global kinetic energy on "tau" velocity
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          fx = rho0*p5*csu(jrow)*dyu(jrow)
#  if defined O_symmetry
          if (jrow .eq. jmtm1) fx = fx*p5
#  endif
          do k=1,km
            do i=is,ie
              weight    = fx*dzt(k)*dxu(i)
              temp(i,k) = u(i,k,j,n,tau)**2*weight
            enddo
            do i=is,ie
              ektot(k,jrow) = ektot(k,jrow) + temp(i,k)
            enddo
          enddo
        enddo
      endif
# endif

# if defined O_energy_analysis

!-----------------------------------------------------------------------
!     diagnostic: integrate work done by the r.h.s. terms in the
!                  momentum equations.
!-----------------------------------------------------------------------

      if (glents .and. eots) call ge1 (joff, js, je, is, ie, n)
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the momentum equations
!                 over specified regional volumes
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call utb1 (joff, js, je, is, ie, n)
# endif
# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s terms in the momentum equation
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call uxbt1 (joff, js, je, n)
# endif
      return
      end

      subroutine diagc2 (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!     construct diagnostics requiring internal mode velocity at "tau+1"
!     and those not dependent on velocity component fluxes.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!-----------------------------------------------------------------------

      implicit none

      integer joff, js, je, is, ie

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "scalar.h"
      include "switch.h"

# if defined O_energy_analysis

!-----------------------------------------------------------------------
!     diagnostic: integrate work done by du/dt in the momentum equations
!                 the external mode part of "u" at "tau+1" will be
!                 accounted for after the external mode is solved.
!                 also, integrate the work done by buoyancy.
!-----------------------------------------------------------------------

      if (glents .and. eots) then
        call ge2 (joff, js, je, is, ie, kmt, kmu, c2dtuv, grav, rho0r)
      endif
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: add du/dt and implicit coriolis terms to the integrals
!                 over specified volumes. the external mode parts will
!                 be accounted for after the external mode is solved.
!-----------------------------------------------------------------------

      if (trmbts .and. eots) then
        call utb2 (joff, js, je, is, ie, c2dtuv, acor)
      endif
# endif

# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate du/dt and implicit coriolis terms from the
!                 momentum equations
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call uxbt2 (joff, js, je, c2dtuv, acor)
# endif
      return
      end

      subroutine asbcu (joff, js, je, is, ie, iu, iv)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C.(surface boundary conditions)

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       iu   = index for u component
!       iv   = index for v component

!     reference: Pacanowski, R.C., Effect of Equatorial Currents
!                on Surface Stress (JPO, Vol 17, No. 6, June 1987)
!-----------------------------------------------------------------------

      implicit none

      integer iu, iv, j, js, je, jrow, joff, i, is, ie

      real rts

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!     initialize S.B.C. at the beginning of each ocean segment
!     (do not alter values in land)

      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = c0
              sbc(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif

!     accumulate surface currents for the Atmos S.B.C. every time step

      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbc(i,jrow,iu) = sbc(i,jrow,iu) + p25*(
     &                       u(i,1,j,1,tau) + u(i-1,1,j,1,tau)
     &                     + u(i,1,j-1,1,tau) + u(i-1,1,j-1,1,tau))
            sbc(i,jrow,iv) = sbc(i,jrow,iv) + p25*(
     &                       u(i,1,j,2,tau) + u(i-1,1,j,2,tau)
     &                     + u(i,1,j-1,2,tau) + u(i-1,1,j-1,2,tau))
          enddo
        enddo
      endif

!     average the surface currents for the Atmos S.B.C. at the end
!     of each ocean segment. (do not alter values in land)

      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = rts*sbc(i,jrow,iu)
              sbc(i,jrow,iv) = rts*sbc(i,jrow,iv)
            endif
          enddo
        enddo
      endif
#endif

      return
      end

#if defined O_mom && defined O_ice_evp
      subroutine isbcu (joff, js, je, is, ie, iu, iv)

!-----------------------------------------------------------------------
!     construct the Ice S.B.C.(surface boundary conditions)

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       iu   = index for u component
!       iv   = index for v component

!-----------------------------------------------------------------------

      implicit none

      integer iu, iv, j, js, je, jrow, joff, i, is, ie

      real rts

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!     initialize S.B.C. at the beginning of each ocean segment
!     (do not alter values in land)

      if (eots .and. osegs .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = c0
              sbc(i,jrow,iv) = c0
            endif
          enddo
        enddo
      endif

!     accumulate geostrophic currents (level 2) for the Ice S.B.C.
!     every time step

      if (eots .and. iu .ne. 0 .and. iv .ne. 0) then
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            sbc(i,jrow,iu) = sbc(i,jrow,iu) + u(i,2,j,1,tau)
            sbc(i,jrow,iv) = sbc(i,jrow,iv) + u(i,2,j,2,tau)
          enddo
        enddo
      endif

!     average the currents for the Ice S.B.C. at the end of each ocean
!     segment. (do not alter values in land)

      if (eots .and. osege .and. iu .ne. 0 .and. iv .ne. 0) then
        rts = c1/ntspos
        do j=js,je
          jrow  = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) then
              sbc(i,jrow,iu) = rts*sbc(i,jrow,iu)
              sbc(i,jrow,iv) = rts*sbc(i,jrow,iv)
            endif
          enddo
        enddo
      endif

      return
      end
#endif

#if defined O_implicitvmix
      subroutine ivdifu (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     solve vertical diffusion of velocity implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = velocity component
!       twodt = (2*dtuv, dtuv) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real c1, up1, r2dtuv, vvca

      include "size.h"
      include "param.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

      real twodt(km)

!     set some constants

      c1 = 1.0

!     construct the "tau+1" velocity without implicit vertical diffusion

      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) = u(i,k,j,n,taum1)+c2dtuv*u(i,k,j,n,taup1)
          enddo
        enddo
      enddo
# if defined O_xbts || defined O_energy_analysis || defined O_term_balances

!     store terms to compute implicit vertical diffusion on
!     diagnostic time steps

      if ((xbtperts .or. glents .or. trmbts) .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = u(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
# endif

!     add in the implicit vertical diffusion

      do k=1,km
        twodt(k) = c2dtuv
      enddo
      call invtri (u(1,1,1,n,taup1), smf(1,1,n), bmf(1,1,n)
     &, visc_cbu(1,1,jsmw), twodt, kmu, umask(1,1,1), is, ie
     &, joff, js, je)

      r2dtuv = c1/c2dtuv
# if defined O_xbts || defined O_energy_analysis || defined O_term_balances

!     compute residual implicit vertical diffusion for diagnostics

      if ((xbtperts .or. glents .or. trmbts) .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = r2dtuv*(u(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
# endif

!     convert back to time change of velocity

      do j=js,je
        do k=1,km
          do i=is,ie
            u(i,k,j,n,taup1) =r2dtuv*(u(i,k,j,n,taup1)-u(i,k,j,n,taum1))
          enddo
        enddo
      enddo

      return
      end
#endif
