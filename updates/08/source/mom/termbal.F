      subroutine utb1 (joff, js, je, is, ie, n)

#if defined O_mom && defined O_term_balances
!=======================================================================
!     accumulate terms in the momentum equations over the
!     volume of the specified regions

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       n    = velocity component
!=======================================================================

      implicit none

      integer i, k, j, n, jrow, js, je, joff, is, ie, nreg

      real adv_ux, adv_uy, adv_uz, adv_metric, diff_ux, diff_uz
      real diff_uy, diff_metric, coriolis, fx, boxvol, term, dudx
      real dvdy, dwdz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "vmixc.h"
      include "fdifm.h"

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        fx = csu(jrow)*dyu(jrow)
# if defined O_symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif

!-----------------------------------------------------------------------
!       accumulate terms for all regions within the current jrow
!-----------------------------------------------------------------------

        do k=1,km
          do i=is,ie
            nreg = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
            if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
              boxvol = fx*dxu(i)*dzt(k)

!-----------------------------------------------------------------------
!             pressure term
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*grad_p(i,k,j,n)
              call addto (termbm(k,2,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             zonal advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Ux(i,k,j)
              call addto (termbm(k,3,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             pure zonal advection of momentum
!-----------------------------------------------------------------------

!             - U(U)x = U(U)x - (UU)x (when n=1)
!             - U(V)x = V(U)x - (UV)x (when n=2)

              dudx = (adv_veu(i,k,j)-adv_veu(i-1,k,j))*dxur(i)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dudx - ADV_Ux(i,k,j))
              call addto (termbm(k,14,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             advective metric term
!-----------------------------------------------------------------------

              term = ADV_metric(i,k,j,jrow,n)
              call addto (termbm(k,13,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             meridional advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Uy(i,k,j,jrow,n)
              call addto (termbm(k,4,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             pure meridional advection of momentum
!-----------------------------------------------------------------------

!             - V(U)y = U(V)y - (VU)y (when n=1)
!             - V(V)y = V(V)y - (VV)y (when n=2)

              dvdy = (adv_vnu(i,k,j)-adv_vnu(i,k,j-1))*dyur(jrow)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dvdy
     &             - ADV_Uy(i,k,j,jrow,n))
              call addto (termbm(k,15,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             vertical advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Uz(i,k,j)
              call addto (termbm(k,5,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             pure vertical advection of momentum
!-----------------------------------------------------------------------

!             - W(U)z = U(W)z - (WU)z (when n=1)
!             - W(V)z = V(W)z - (WV)z (when n=2)

              dwdz = (adv_vbu(i,k-1,j)-adv_vbu(i,k,j))*dztr(k)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dwdz - ADV_Uz(i,k,j))
              call addto (termbm(k,16,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             zonal diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Ux(i,k,j)
              call addto (termbm(k,6,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             meridional diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Uy(i,k,j,jrow,n)
              call addto (termbm(k,7,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             diffusive metric term
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_metric(i,k,j,jrow,n)
              call addto (termbm(k,9,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             vertical diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Uz(i,k,j)
# if defined O_implicitvmix
     &              +umask(i,k,j)*zzi(i,k,j)
# endif
              call addto (termbm(k,8,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             coriolis term
!-----------------------------------------------------------------------

              term = umask(i,k,j)*CORIOLIS(i,k,j,jrow,n)
              call addto (termbm(k,10,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!             accumulate the source term
!-----------------------------------------------------------------------

# if defined O_source_term || defined O_mobi || defined O_carbon_14
              term = umask(i,k,j)*source(i,k,j)
              call addto (termbm(k,11,n,nreg), term*boxvol)
# else
              termbm(k,11,n,nreg) = 0.0
# endif

!-----------------------------------------------------------------------
!             accumulate u, v, and w
!-----------------------------------------------------------------------

              term = umask(i,k,j)*u(i,k,j,n,tau)
              call addto (termbm(k,17,n,nreg), term*boxvol)

              if (n .eq. 2) then
                term = p5*(adv_vbu(i,k,j)+adv_vbu(i,k-1,j))*umask(i,k,j)
                call addto (avgw(nreg), term*boxvol)
              endif

!-----------------------------------------------------------------------
!             accumulate the surface momentum flux
!-----------------------------------------------------------------------

              if (k .eq. 1) then
                term = umask(i,k,j)*smf(i,j,n)
                call addto (smflx(n,nreg), term*fx*dxu(i))
              endif
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine utb2 (joff, js, je, is, ie, c2dtuv, acor)

!=======================================================================
!     accumulate external mode parts of d/dt and the implicit coriolis
!     term in the momentum equations over the volume in the specified
!     regions

!     input:
!       joff   = offset relating "j" in the MW to latitude "jrow"
!       js     = starting row in the MW
!       je     = ending row in the MW
!       is     = starting longitude index in the MW
!       ie     = ending longitude index in the MW
!       c2dtuv = (2*dtuv,dtuv) on (lpfrod,mixing) time steps
!       acor   = implicit factor
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, n, k, i, is, ie, nreg

      real r2dt, c2dtuv, fx, boxvol, term, acor

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "grdvar.h"
      include "mw.h"

!-----------------------------------------------------------------------
!       local constants
!-----------------------------------------------------------------------

      r2dt = c1/c2dtuv

      do j=js,je
        jrow = j + joff
        fx   = csu(jrow)*dyu(jrow)
# if defined O_symmetry
        if (jrow .eq. jmtm1) fx = fx*p5
# endif
        do n=1,2
          do k=1,km
            do i=is,ie
              nreg   = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
                boxvol = fx*dxu(i)*dzt(k)

!-----------------------------------------------------------------------
!               d/dt of velocity (external mode part of tau+1 will be
!               added later when the external mode is solved)
!-----------------------------------------------------------------------

                term = umask(i,k,j)*(u(i,k,j,n,taup1) -
     &                               u(i,k,j,n,taum1))*r2dt
                call addto (termbm(k,1,n,nreg), term*boxvol)

!-----------------------------------------------------------------------
!               implicit coriolis term (external mode part will be added
!               later when external mode is solved)
!-----------------------------------------------------------------------

                if (acor .ne. c0) then
                  term = umask(i,k,j)*acor*cori(i,jrow,n)*
     &                       (u(i,k,j,3-n,taup1) - u(i,k,j,3-n,taum1))
                  call addto (termbm(k,10,n,nreg), term*boxvol)
                endif
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine utb3

!=======================================================================
!     accumulate external mode parts of d/dt, the implicit coriolis
!     term and the surface pressure gradientsover the volume in the
!     specified regions.
!=======================================================================

      implicit none

      integer is, ie, js, je, jrow, i, kz, k, n

      real fddt, fspr, atosp, f1, uext, vext, boxfac, boxspr, boxacr

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "diag.h"

      parameter (is=1, ie=1, js=1, je=1)
      real psgrad(is:ie,js:je,2)

      do jrow=1,jmt-1
        fddt  = csu(jrow)*dyu(jrow)/c2dtuv
        fspr  = csu(jrow)*dyu(jrow)
        do i=2,imtm1
          atosp = acor*cori(i,jrow,1)
          f1    = atosp*csu(jrow)*dyu(jrow)
          kz = kmu(i,jrow)
          if (kz .ne. 0) then
            do k=1,kz
              n = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (n .gt. 0 .and. mskhr(i,jrow) .gt. 0) then

!               construct the surface pressure gradients for pt (i,jrow)

                if (k .eq. 1) then
                  call calc_psgrad(psgrad, uext, vext, jrow, jrow, i, i)
                endif
                boxfac = fddt*dxu(i)*dzt(k)
                boxspr = fspr*dxu(i)*dzt(k)
                termbm(k,1,1,n)  = termbm(k,1,1,n)  + uext*boxfac
                termbm(k,1,2,n)  = termbm(k,1,2,n)  + vext*boxfac
                termbm(k,12,1,n) = termbm(k,12,1,n) -
     &                              psgrad(is,js,1)*boxspr
                termbm(k,12,2,n) = termbm(k,12,2,n) -
     &                              psgrad(is,js,2)*boxspr
                boxacr = f1*dxu(i)*dzt(k)
                termbm(k,10,1,n) = termbm(k,10,1,n) + vext*boxacr
                termbm(k,10,2,n) = termbm(k,10,2,n) - uext*boxacr
              endif
            enddo
          endif
        enddo
      enddo

      return
      end

      subroutine ttb1 (joff, js, je, is, ie, n)

!=======================================================================
!     accumulate terms in the tracer equations over the volume in the
!     specified regions

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, istrt, is, iend, ie, js
      integer je, joff, nreg

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, fx, area, boxvol
      real term, r2dt, dudx, dvdy, dwdz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "vmixc.h"

# if defined O_isopycmix || defined O_isneutralmix
      include "isopyc.h"
# endif
      include "fdift.h"

!-----------------------------------------------------------------------
!     limit the longitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

      do j=js,je
        jrow = j + joff

!-----------------------------------------------------------------------
!       set local constants
!-----------------------------------------------------------------------

        fx   = cst(jrow)*dyt(jrow)
        do k=1,km
          do i=istrt,iend
            nreg = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
            if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
              area   = fx*dxt(i)
              boxvol = area*dzt(k)

!-----------------------------------------------------------------------
!             tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*t(i,k,j,n,tau)
              call addto (termbt(k,15,n,nreg), term*boxvol)
              call addto (termbt(k,15,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             d(tracer)/dt
!-----------------------------------------------------------------------

              r2dt = c1/(c2dtts*dtxcel(k))
              term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                             t(i,k,j,n,taum1))*r2dt
              call addto (termbt(k,9,n,nreg), term*boxvol)
              call addto (termbt(k,9,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             zonal advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Tx(i,k,j)
# if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Txiso(i,k,j,n)
# endif
              call addto (termbt(k,2,n,nreg), term*boxvol)
              call addto (termbt(k,2,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             pure zonal advection of tracer
!-----------------------------------------------------------------------

!             - U(T)x = T(U)x - (UT)x

              dudx = (adv_vet(i,k,j)-adv_vet(i-1,k,j))*dxtr(i)
     &               *cstr(jrow)
# if defined O_gent_mcwilliams && !defined O_fct
     &              +(adv_vetiso(i,k,j)-adv_vetiso(i-1,k,j))*dxtr(i)
     &               *cstr(jrow)
# endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dudx - ADV_Tx(i,k,j))
# if defined O_gent_mcwilliams && !defined O_fct
     &              -tmask(i,k,j)*ADV_Txiso(i,k,j,n)
# endif
              call addto (termbt(k,11,n,nreg), term*boxvol)
              call addto (termbt(k,11,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             meridional advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Ty(i,k,j,jrow,n)
# if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Tyiso(i,k,j,jrow,n)
# endif
              call addto (termbt(k,3,n,nreg), term*boxvol)
              call addto (termbt(k,3,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             pure meridional advection of tracer
!-----------------------------------------------------------------------

!             - V(T)y = T(V)y - (VT)y

              dvdy = (adv_vnt(i,k,j)-adv_vnt(i,k,j-1))*dytr(jrow)
     &               *cstr(jrow)
# if defined O_gent_mcwilliams && !defined O_fct
     &             + (adv_vntiso(i,k,j)-adv_vntiso(i,k,j-1))*dytr(jrow)
     &               *cstr(jrow)
# endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dvdy
     &             - ADV_Ty(i,k,j,jrow,n))
# if defined O_gent_mcwilliams && !defined O_fct
     &              -tmask(i,k,j)*ADV_Tyiso(i,k,j,jrow,n)
# endif
              call addto (termbt(k,12,n,nreg), term*boxvol)
              call addto (termbt(k,12,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             vertical advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Tz(i,k,j)
#  if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Tziso(i,k,j)
#  endif
              call addto (termbt(k,4,n,nreg), term*boxvol)
              call addto (termbt(k,4,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             pure vertical advection of tracer
!-----------------------------------------------------------------------

!             - W(T)z = T(W)z - (WT)z

              dwdz = (adv_vbt(i,k-1,j)-adv_vbt(i,k,j))*dztr(k)
#  if defined O_gent_mcwilliams && !defined O_fct
     &              + (adv_vbtiso(i,k-1,j)-adv_vbtiso(i,k,j))*dztr(k)
#  endif
              term = tmask(i,k,j)*(t(i,k,j,n,tau)*dwdz - ADV_Tz(i,k,j))
#  if defined O_gent_mcwilliams && !defined O_fct
     &              -tmask(i,k,j)*ADV_Tziso(i,k,j)
#  endif
              call addto (termbt(k,13,n,nreg), term*boxvol)
              call addto (termbt(k,13,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             zonal diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Tx(i,k,j)
              call addto (termbt(k,5,n,nreg), term*boxvol)
              call addto (termbt(k,5,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             meridional diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Ty(i,k,j,jrow,n)
              call addto (termbt(k,6,n,nreg), term*boxvol)
              call addto (termbt(k,6,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             vertical diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Tz(i,k,j)
# if defined O_implicitvmix || defined O_isopycmix
     &               + tmask(i,k,j)*zzi(i,k,j)
# endif
              call addto (termbt(k,7,n,nreg), term*boxvol)
              call addto (termbt(k,7,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!             tracer source term
!-----------------------------------------------------------------------

# if defined O_source_term || defined O_mobi || defined O_carbon_14
              term = tmask(i,k,j)*source(i,k,j)
              call addto (termbt(k,8,n,nreg), term*boxvol)
              call addto (termbt(k,8,n,0),    term*boxvol)
# else
              termbt(k,8,n,nreg) = 0.0
              termbt(k,8,n,0)    = 0.0
# endif

              if (k .eq. 1) then

!-----------------------------------------------------------------------
!               surface tracer
!-----------------------------------------------------------------------

                term = tmask(i,k,j)*t(i,k,j,n,tau)
                call addto (asst(n,nreg), term*area)
                call addto (asst(n,0),    term*area)

!-----------------------------------------------------------------------
!               surface tracer flux
!-----------------------------------------------------------------------

                term = tmask(i,k,j)*stf(i,j,n)
                call addto (stflx(n,nreg), term*area)
                call addto (stflx(n,0), term*area)
              endif
            endif
          enddo
        enddo
      enddo

      return
      end

      subroutine ttb2 (joff, js, je, is, ie, iterm)

!=======================================================================
!     accumulate d/dt and change in tracer in the tracer equations over
!     the volume in the specified regions

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW

!     iterm = 1  => total change
!     iterm = 10 => change due to filtering
!=======================================================================

      implicit none

      integer iterm, j, js, je, jrow, joff, n, k, i, is, ie, nreg

      real fx, r2dt, area, boxvol, term

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "grdvar.h"
      include "mw.h"
      include "scalar.h"

       if (iterm .ne. 1 .and. iterm .ne. 10) then
        write (stdout,*) '=>Error: iterm=',iterm,' in ttb2'
        stop '=>ttb2'
      endif

      do j=js,je
        jrow = j + joff
        fx   = cst(jrow)*dyt(jrow)
        do n=1,nt
          do k=1,km
            r2dt = c1/(c2dtts*dtxcel(k))
            do i=is,ie
              nreg   = nhreg*(mskvr(k)-1) + mskhr(i,jrow)
              if (nreg .gt. 0 .and. mskhr(i,jrow) .gt. 0) then
                area   = fx*dxt(i)
                boxvol = area*dzt(k)

!-----------------------------------------------------------------------
!               d/dt(tracer)
!-----------------------------------------------------------------------

                term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                               t(i,k,j,n,taum1))*r2dt
                call addto (termbt(k,iterm,n,nreg), term*boxvol)
                call addto (termbt(k,iterm,n,0),    term*boxvol)

!-----------------------------------------------------------------------
!               change in variance of tracer
!-----------------------------------------------------------------------

                if (iterm .eq. 1) then
                  term = tmask(i,k,j)*(t(i,k,j,n,taup1)**2-
     &                                 t(i,k,j,n,taum1)**2)
                  call addto (termbt(k,14,n,nreg), term*boxvol)
                  call addto (termbt(k,14,n,0),    term*boxvol)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
#endif

      return
      end
