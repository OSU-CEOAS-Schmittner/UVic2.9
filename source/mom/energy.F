      subroutine ge1 (joff, js, je, is, ie, n)

#if defined O_mom && defined O_energy_analysis
!-----------------------------------------------------------------------
!     compute global energetics by taking u dot the momentum equations
!     and integrating over the ocean volume

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = velocity component
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, jrow, n, js, je, joff, is, ie

      real adv_ux, adv_uy, adv_uz, adv_metric, diff_ux, diff_uz
      real diff_uy, diff_metric, coriolis, acor
      real diag1, diag0, fx, uext, uint, boxvol, term, dhdt, gcor

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "emode.h"
      include "diag.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "vmixc.h"
      include "fdifm.h"

      real ext(imt,2)

      do j=js,je
        jrow = j + joff
# if defined O_stream_function

!-----------------------------------------------------------------------
!       compute the external mode velocities
!-----------------------------------------------------------------------

        do i=is,ie
          diag1    = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
          diag0    = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
          ext(i,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
          ext(i,2) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
        enddo
# endif

!-----------------------------------------------------------------------
!       compute work done by each term on the internal and external
!       modes in the momentum equations
!-----------------------------------------------------------------------

        fx = csu(jrow)*dyu(jrow)
# if defined O_symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif

        do k=1,km
          do i=is,ie
# if defined O_stream_function
            uext   = ext(i,n)
# else
            uext   = ubar(i,jrow,n)
# endif
            uint   = u(i,k,j,n,tau) - uext
            boxvol = fx*dxu(i)*dzt(k)

!-----------------------------------------------------------------------
!           pressure term
!-----------------------------------------------------------------------

            term = -umask(i,k,j)*grad_p(i,k,j,n)
            call addto (engint(k,6,jrow), uint*term*boxvol)
            call addto (engext(6,jrow),   uext*term*boxvol)

!-----------------------------------------------------------------------
!           coriolis term does no work: u*fv - v*fu = 0
!           implicit coriolis work will be reflected in the imbalance
!           between horizontal pressure forces and buoyancy
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!           zonal and meridional advection of momentum + metric term
!-----------------------------------------------------------------------

            term =-umask(i,k,j)*(ADV_Ux(i,k,j) + ADV_Uy(i,k,j,jrow,n)
     &           + ADV_metric(i,k,j,jrow,n))
            call addto (engint(k,2,jrow), uint*term*boxvol)
            call addto (engext(2,jrow),   uext*term*boxvol)

!-----------------------------------------------------------------------
!           vertical advection of momentum
!-----------------------------------------------------------------------

            term = -umask(i,k,j)*ADV_Uz(i,k,j)
            call addto (engint(k,3,jrow), uint*term*boxvol)
            call addto (engext(3,jrow),   uext*term*boxvol)
# if defined O_implicit_free_surface

!           add effect due to change in volume of top layer

            if (k .eq. 1) then
              dhdt = p5*umask(i,1,j)*fx*dxu(i)*adv_vbu(i,0,j)
              call addto (engext(3,jrow),  dhdt*u(i,1,j,n,tau)**2)
            endif
# endif

!-----------------------------------------------------------------------
!           zonal and meridional diffusion of momentum + metric term
!-----------------------------------------------------------------------

            term = umask(i,k,j)*(DIFF_Ux(i,k,j) + DIFF_Uy(i,k,j,jrow,n)
     &           + DIFF_metric(i,k,j,jrow,n))
            call addto (engint(k,4,jrow), uint*term*boxvol)
            call addto (engext(4,jrow),   uext*term*boxvol)

!-----------------------------------------------------------------------
!           vertical diffusion of momentum
!-----------------------------------------------------------------------

            term = umask(i,k,j)*DIFF_Uz(i,k,j)
# if defined O_implicitvmix
     &                                    + umask(i,k,j)*zzi(i,k,j)
# endif
            call addto (engint(k,5,jrow), uint*term*boxvol)
            call addto (engext(5,jrow),   uext*term*boxvol)
          enddo
        enddo

!-----------------------------------------------------------------------
!       work done by wind stress
!-----------------------------------------------------------------------

        k = 1
        do i=is,ie
# if defined O_stream_function
          uext = ext(i,n)
# else
          uext = ubar(i,jrow,n)
# endif
          uint = u(i,k,j,n,tau) - uext
          term = umask(i,k,j)*smf(i,j,n)
          call addto (engint(k,7,jrow), uint*term*fx*dxu(i))
          call addto (engext(7,jrow),   uext*term*fx*dxu(i))
        enddo

!-----------------------------------------------------------------------
!       work done by bottom drag
!-----------------------------------------------------------------------

        do i=is,ie
# if defined O_stream_function
          uext = ext(i,n)
# else
          uext = ubar(i,jrow,n)
# endif
          k = kmu(i,jrow)
          if (k .ne. 0) then
            uint = u(i,k,j,n,tau) - uext
            term = -umask(i,k,j)*bmf(i,j,n)
            call addto (engint(k,8,jrow), uint*term*fx*dxu(i))
            call addto (engext(8,jrow),   uext*term*fx*dxu(i))
          endif
        enddo
      enddo

      return
      end

      subroutine ge2 (joff, js, je, is, ie, kmt, kmu, c2dtuv, grav
     &,               rho0r)

!-----------------------------------------------------------------------
!     compute global energetics by taking u dot the momentum equations
!     and integrating over the entire ocean volume

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, jrow, joff, i, is, ie, n, k, kz

      real fx, diag1, diag0, r2dt, c2dtuv, uext, uint, boxvol, term
      real f1, area, grav, rho0r, udxdy, tdxdy

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "emode.h"
      include "diag.h"
      include "grdvar.h"
      include "mw.h"

      integer kmt(imt,jmt), kmu(imt,jmt)

      real ext(imt,2)

      do j=js,je
        jrow  = j + joff

!-----------------------------------------------------------------------
!       set local constants
!-----------------------------------------------------------------------

        fx = csu(jrow)*dyu(jrow)
# if defined O_symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif
# if defined O_stream_function

!-----------------------------------------------------------------------
!       compute the external mode velocities
!-----------------------------------------------------------------------

        do i=is,ie
          diag1    = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
          diag0    = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
          ext(i,1) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
          ext(i,2) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
        enddo
# endif

!-----------------------------------------------------------------------
!       compute work done by each term on the internal and external
!       modes in the momentum equations
!-----------------------------------------------------------------------

        r2dt = c1/c2dtuv
        do n=1,2
          do k=1,km
            do i=is,ie
# if defined O_stream_function
              uext = ext(i,n)
# else
              uext = ubar(i,jrow,n)
# endif
              uint   = u(i,k,j,n,tau) - uext
              boxvol = fx*dxu(i)*dzt(k)

!-----------------------------------------------------------------------
!             total change in kinetic energy.
!-----------------------------------------------------------------------

              term = umask(i,k,j)*(u(i,k,j,n,taup1)-
     &                             u(i,k,j,n,taum1))*r2dt
              call addto (engint(k,1,jrow), uint*term*boxvol)
               call addto (engext(1,jrow),   uext*term*boxvol)
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute the work done by buoyancy integrated over the entire
!     ocean volume.
!-----------------------------------------------------------------------

      do j=js,je
        jrow  = j + joff
        f1 = cst(jrow)*dyt(jrow)
        do i=is,ie
          kz = kmt(i,jrow)
          if (kz .ne. 0) then
            area = f1*dxt(i)
            fx   = area*grav*rho0r*p5
# if defined O_pressure_gradient_average
            term =  - c2*fx*dzw(0)*adv_vbt(i,0,j)*rhotilde(i,1,j)
# else
            term =  - c2*fx*dzw(0)*adv_vbt(i,0,j)*rho(i,1,j)
# endif
            call addto (buoy(1,jrow), term)
            do k=2,kz
              term =-fx*dzw(k-1)*adv_vbt(i,k-1,j)*
# if defined O_pressure_gradient_average
     &               (rhotilde(i,k-1,j) + rhotilde(i,k,j))
# else
     &               (rho(i,k-1,j) + rho(i,k,j))
# endif
              call addto (buoy(k,jrow), term)
            enddo
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     find maximum error in continuity for "t" cells and "u" cells
!-----------------------------------------------------------------------

      do j=js,je
        jrow  = j + joff
        do k=1,km
          do i=is,ie
            term =
     &        ((adv_vet(i,k,j) - adv_vet(i-1,k,j))*cstr(jrow)*dxtr(i)
     &        +(adv_vnt(i,k,j) - adv_vnt(i,k,j-1))*cstr(jrow)*dytr(jrow)
     &        +(adv_vbt(i,k-1,j) - adv_vbt(i,k,j))*dztr(k))*tmask(i,k,j)
            if (abs(term) .gt. abs(tcerr(jrow))) then
              tcerr(jrow) = term
              itcerr(jrow) = i
              jtcerr(jrow) = jrow
              ktcerr(jrow) = k
            endif

            term =
     &        ((adv_veu(i,k,j) - adv_veu(i-1,k,j))*csur(jrow)*dxur(i)
     &        +(adv_vnu(i,k,j) - adv_vnu(i,k,j-1))*csur(jrow)*dyur(jrow)
     &        +(adv_vbu(i,k-1,j) - adv_vbu(i,k,j))*dztr(k))*umask(i,k,j)
            if (abs(term) .gt. abs(ucerr(jrow))) then
              ucerr(jrow) = term
              iucerr(jrow) = i
              jucerr(jrow) = jrow
              kucerr(jrow) = k
            endif
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     find max error in "adv_vbt" at bottom and max "adv_vbu" at bottom
!-----------------------------------------------------------------------

      do j=js,je
        jrow  = j + joff
        do i=is,ie
          k = kmt(i,jrow)
          if (k.ne.0 .and.
     &        (abs(adv_vbt(i,k,j)) .gt. abs(wtbot(jrow)))) then
            wtbot(jrow)  = adv_vbt(i,k,j)
            iwtbot(jrow) = i
            jwtbot(jrow) = jrow
            kwtbot(jrow) = k
          endif
        enddo
      enddo

      do j=js,je
        jrow  = j + joff
        do i=is,ie
          k = kmu(i,jrow)
          if (k.ne.0 .and.
     &        (abs(adv_vbu(i,k,j)) .gt. abs(wubot(jrow))))then
            wubot(jrow)  = adv_vbu(i,k,j)
            iwubot(jrow) = i
            jwubot(jrow) = jrow
            kwubot(jrow) = k
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     integrate "adv_vbt" and "adv_vbu" for all lat and lon
!-----------------------------------------------------------------------

      do j=js,je
        jrow  = j + joff
        do k=1,km
          do i=is,ie
            udxdy = umask(i,k,j)*dxu(i)*csu(jrow)*dyu(jrow)
            tdxdy = tmask(i,k,j)*dxt(i)*cst(jrow)*dyt(jrow)
            wtlev(k,jrow)  = wtlev(k,jrow) + adv_vbt(i,k,j)*tdxdy
            wulev(k,jrow)  = wulev(k,jrow) + adv_vbu(i,k,j)*udxdy
          enddo
        enddo
      enddo

      return
      end

      subroutine ge3 (c2dtuv)

!-----------------------------------------------------------------------
!     compute global energetics by taking u dot the momentum equations
!     and integrating over the entire ocean volume
!-----------------------------------------------------------------------

      implicit none

      integer i, jrow, m

      real  fx, diag1, diag0, r2dt, c2dtuv, uext, vext, uextn
      real vextn, boxvol, term

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "emode.h"
      include "diag.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"

      real  ext(imt,2,2)

      do jrow=2,jmt-1
        fx = csu(jrow)*dyu(jrow)
# if defined O_symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif
# if defined O_stream_function

!-----------------------------------------------------------------------
!       compute the external mode velocities after external mode has
!       been updated. since the external mode has been updated,
!       m = 1 is "tau+1"
!       m = 2 is "tau"
!-----------------------------------------------------------------------

        do m=1,2
          do i=2,imtm1
            diag1      = psi(i+1,jrow+1,m) - psi(i  ,jrow,m)
            diag0      = psi(i  ,jrow+1,m) - psi(i+1,jrow,m)
            ext(i,1,m) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
            ext(i,2,m) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
          enddo
        enddo
# endif

!-----------------------------------------------------------------------
!       compute external mode part of work done by "tau+1"
!       component of d/dt. internal mode part is zero. since the
!       external mode has been updated:
!       ubarm1 is "tau"
!       ubar   is "tau+1"
!       note: if using 5 point numerics when solving for the stream
!             function, the total integral is not conserved  for the
!             external mode which shows itself as the "ficticious" term
!             in the global energy integrals.
!-----------------------------------------------------------------------

        r2dt = c1/c2dtuv
        do i=2,imtm1
# if defined O_stream_function
          uext   = ext(i,1,2)
          vext   = ext(i,2,2)
          uextn  = ext(i,1,1)
          vextn  = ext(i,2,1)
# else
          uext   = ubarm1(i,jrow,1)
          vext   = ubarm1(i,jrow,2)
          uextn  = ubar(i,jrow,1)
          vextn  = ubar(i,jrow,2)
# endif
          boxvol = fx*dxu(i)*h(i,jrow)
          term = (uext*uextn + vext*vextn)*r2dt
           call addto (engext(1,jrow),   term*boxvol)
        enddo

      enddo

      return
      end

      subroutine addto (sum, term)

      implicit none

      real sum, term

      sum = sum + term
#endif
      return
      end
