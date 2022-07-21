      subroutine adv_flux (joff, js, je, is, ie, n)

#if defined O_mom
# if defined O_linearized_advection
!=======================================================================
!     Linearized advective tracer flux

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer istrt, iend, j, js, je, k, i, n, joff, ie, is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      parameter (istrt=2, iend=imt-1)


!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!     istrt = max(2,is)
!     iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     advective flux across eastern and northern face of "T" cells.
!     is zero due to linearization about state of no motion
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     calculate 2*advective flux across bottom face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km-1
          do i=istrt,iend
            adv_fb(i,k,j)  = adv_vbt(i,k,j)*(tbarz(k,n) + tbarz(k+1,n))
          enddo
        enddo
      enddo

# elif defined O_quicker
!=======================================================================
!     3rd order advective tracer flux

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!        n   = tracer
!=======================================================================

      implicit none

      integer istrt, iend, i, k, j, ip, kr, jq, lag, js, je, ip2, n
      integer joff, jrow, jp2, ie, is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "mw.h"

      parameter (istrt=2, iend=imt-1)

      real totvel, upos, uneg, vpos, vneg, wpos, wneg


#  if defined O_gent_mcwilliams
      include "isopyc.h"
#  endif

#  if defined O_ncar_upwind3
      lag = tau
#  else
      lag = taum1
#  endif

!-----------------------------------------------------------------------
!     calculate 2*advective flux across eastern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km
          do i=istrt,iend-1
            ip2 = i + 2
            totvel = adv_vet(i,k,j)
#  if defined O_gent_mcwilliams
     &              + adv_vetiso(i,k,j)
#  endif
            upos = p5*(totvel + abs(totvel))
     &             *tmask(i-1,k,j)*tmask(i,k,j)*tmask(i+1,k,j)
            uneg = p5*(totvel - abs(totvel))
     &             *tmask(ip2,k,j)*tmask(i+1,k,j)*tmask(i,k,j)

            adv_fe(i,k,j) = totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau)
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)*t(i-1,k,j,n,lag))
     &                  - uneg*(curv_xn(i,1)*t(ip2,k,j,n,lag)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
          enddo
        enddo

        do k=1,km
          i=iend
          ip2 = 3
          totvel = adv_vet(i,k,j)
#  if defined O_gent_mcwilliams
     &              + adv_vetiso(i,k,j)
#  endif
          upos = p5*(totvel + abs(totvel))
     &             *tmask(i-1,k,j)*tmask(i,k,j)*tmask(i+1,k,j)
          uneg = p5*(totvel - abs(totvel))
     &             *tmask(ip2,k,j)*tmask(i+1,k,j)*tmask(i,k,j)

          adv_fe(i,k,j) = totvel*(
     &                          quick_x(i,1)*t(i,  k,j,n,tau)
     &                        + quick_x(i,2)*t(i+1,k,j,n,tau))
     &                  - upos*(curv_xp(i,1)*t(i+1,k,j,n,lag)
     &                         +curv_xp(i,2)*t(i  ,k,j,n,lag)
     &                         +curv_xp(i,3)*t(i-1,k,j,n,lag))
     &                  - uneg*(curv_xn(i,1)*t(ip2,k,j,n,lag)
     &                         +curv_xn(i,2)*t(i+1,k,j,n,lag)
     &                         +curv_xn(i,3)*t(i  ,k,j,n,lag))
        enddo
        call setbcx (adv_fe(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across northern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      if (joff +js .eq. 2) then
        do j=1,1
          do k=1,km
            do i=2,imt-1
              adv_f4n(i,k,j,n) = c0
            enddo
          enddo
        enddo
      endif
      do j=js,je
        jrow = j + joff
        jp2 = min(j+2+joff,jmt) - joff
        do k=1,km
          do i=istrt,iend
            totvel = adv_vnt(i,k,j)
#  if defined O_gent_mcwilliams
     &              + adv_vntiso(i,k,j)
#  endif
            vpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k,j-1)*tmask(i,k,j)*tmask(i,k,j+1)
            vneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k,jp2)*tmask(i,k,j+1)*tmask(i,k,j)

            adv_f4n(i,k,j,n) = totvel*(
     &                          quick_y(jrow,1)*t(i,k,j  ,n,tau)
     &                        + quick_y(jrow,2)*t(i,k,j+1,n,tau))
     &                  - vpos*(curv_yp(jrow,1)*t(i,k,j+1,n,lag)
     &                         +curv_yp(jrow,2)*t(i,k,j  ,n,lag)
     &                         +curv_yp(jrow,3)*t(i,k,j-1,n,lag))
     &                  - vneg*(curv_yn(jrow,1)*t(i,k,jp2,n,lag)
     &                         +curv_yn(jrow,2)*t(i,k,j+1,n,lag)
     &                         +curv_yn(jrow,3)*t(i,k,j  ,n,lag))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across bottom face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=2,km-2
          do i=istrt,iend
            totvel = adv_vbt(i,k,j)
#  if defined O_gent_mcwilliams
     &              + adv_vbtiso(i,k,j)
#  endif
            wpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k+2,j)*tmask(i,k+1,j)*tmask(i,k,j)
            wneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k-1,j)*tmask(i,k,j)*tmask(i,k+1,j)

            adv_fb(i,k,j)  = totvel*(
     &                          quick_z(k,1)*t(i,k  ,j,n,tau)
     &                        + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                  - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                         +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                         +curv_zp(k,3)*t(i,k-1,j,n,lag))
     &                  - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                         +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                         +curv_zn(k,3)*t(i,k  ,j,n,lag))
          enddo
        enddo
        k=1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
#  if defined O_gent_mcwilliams
     &            + adv_vbtiso(i,k,j)
#  endif
          wpos = p5*(totvel + abs(totvel))
     &             *tmask(i,k+2,j)*tmask(i,k+1,j)*tmask(i,k,j)

          adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau)
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wpos*(curv_zn(k,1)*t(i,k+2,j,n,lag)
     &                       +curv_zn(k,2)*t(i,k+1,j,n,lag)
     &                       +curv_zn(k,3)*t(i,k  ,j,n,lag))
        enddo
        k=km-1
        do i=istrt,iend
          totvel = adv_vbt(i,k,j)
#  if defined O_gent_mcwilliams
     &            + adv_vbtiso(i,k,j)
#  endif
          wneg = p5*(totvel - abs(totvel))
     &             *tmask(i,k-1,j)*tmask(i,k,j)*tmask(i,k+1,j)

          adv_fb(i,k,j)  = totvel*(
     &                        quick_z(k,1)*t(i,k  ,j,n,tau)
     &                      + quick_z(k,2)*t(i,k+1,j,n,tau))
     &                - wneg*(curv_zp(k,1)*t(i,k+1,j,n,lag)
     &                       +curv_zp(k,2)*t(i,k  ,j,n,lag)
     &                       +curv_zp(k,3)*t(i,k-1,j,n,lag))
        enddo
      enddo

# elif defined O_fourth_order_tracer_advection
!=======================================================================
!     4th order advective tracer flux

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!        n   = tracer
!=======================================================================

      implicit none

      integer istrt, iend, j, js, je, k, i, mask, n, ip2, joff, jp2, m
      integer ie, is
      real a2nd, b2nd, a4th, b4th

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      parameter (istrt=2, iend=imt-1)


!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!     istrt = max(2,is)
!     iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     calculate 2*advective flux across eastern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      a2nd = 1.0
      b2nd = 0.0
      a4th = 7.0/6.0
      b4th = -1.0/6.0

      do j=js,je
        do k=1,km
          do i=istrt,iend-1
            mask = tmask(i-1,k,j)*tmask(i+2,k,j)
            adv_fe(i,k,j) = adv_vet(i,k,j)*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) +
     &                                      t(i+1,k,j,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i-1,  k,j,n,tau) +
     &                                      t(i+2,k,j,n,tau)))
          enddo
          i = iend
#  if defined O_cyclic
          ip2 = 3
#  else
          ip2 = imt
#  endif
          mask = tmask(i-1,k,j)*tmask(ip2,k,j)
          adv_fe(i,k,j) = adv_vet(i,k,j)*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) +
     &                                      t(i+1,k,j,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i-1,  k,j,n,tau) +
     &                                      t(ip2,k,j,n,tau)))
          adv_fe(1,k,j) = adv_fe(imt-1,k,j)
        enddo
      enddo

!-----------------------------------------------------------------------
!     2*advective flux across northern face of "T" cells is built
!     into ADV_Ty. (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      if (joff + js .eq. 2) then
        do j=1,1
          do k=1,km
            do i=2,imt-1
              adv_f4n(i,k,j,n) = 0.0
            enddo
          enddo
        enddo
      endif
      do j=js,je
        jp2 = min(j+2+joff,jmt) - joff
        do k=1,km
          do i=istrt,iend
            mask = tmask(i,k,j-1)*tmask(i,k,jp2)
            adv_f4n(i,k,j,n) = adv_vnt(i,k,j)*(
     &       (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) +
     &                                      t(i,k,j+1,n,tau))+
     &       (b2nd*(1.0-mask) + b4th*mask)*(t(i,  k,j-1,n,tau) +
     &                                      t(i,k,jp2,n,tau)))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across bottom face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

        do j=js,je
          do k=2,km-2
            do i=istrt,iend
              mask = tmask(i,k-1,j)*tmask(i,k+2,j)
              adv_fb(i,k,j) = adv_vbt(i,k,j)*(
     &         (a2nd*(1.0-mask) + a4th*mask)*(t(i,  k,j,n,tau) +
     &                                        t(i,k+1,j,n,tau))+
     &         (b2nd*(1.0-mask) + b4th*mask)*(t(i,k-1,j,n,tau) +
     &                                        t(i,k+2,j,n,tau)))
            enddo
          enddo
          k = 1
          m = km-1
          do i=istrt,iend
            adv_fb(i,k,j) = adv_vbt(i,k,j)*(t(i,k  ,j,n,tau) +
     &                                      t(i,k+1,j,n,tau))
            adv_fb(i,m,j) = adv_vbt(i,m,j)*(t(i,m  ,j,n,tau) +
     &                                      t(i,m+1,j,n,tau))
          enddo
        enddo

# elif defined O_fct
!=======================================================================
!     computes advective fluxes using a flux-corrected transport scheme

!        for reference see
!        Gerdes, R., C. Koeberle and J. Willebrandt, 1991
!        the influence of numerical advection schemes on the results of
!        ocean general circulation models. Clim Dynamics 5, 211-226
!        and
!        Zalesak, S. T., 1979: Fully multidimensional flux-corrected
!        transport algorithms for fluids. J. Comp. Phys. 31, 335-362.

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       jstrt = starting row in the MW for fct
!       jend  = ending row in the MW for fct
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       istrt = max(2,starting longitude index in the MW)
!       iend  = min(imt-1,ending longitude index in the MW)
!       n     = tracer index

!     output: ( via common mw in "mw.h" )
!       adv_fn = 2*advective flux across northern face of T-cell
!       adv_fe = 2*advective flux across eastern face of T-cell
!       adv_fb = 2*advective flux across bottom face of T-cell
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, istrt, is, iend, ie
      integer istrtm1, iendp1, jstrt, js, jend, je, joff, jlast, jp2
      integer jp1

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, aidif, totadv
      real fxa, fxb, pplus, pminus, qplus, qminus, reltim

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
#  if defined O_gent_mcwilliams
      include "isopyc.h"
#  endif

!---------------------------------------------------------------------
!     dimension local data
!---------------------------------------------------------------------

      real twodt(km), dcf(imt), Trmin(imt), Trmax(imt)
      real Cpos(imt), Cneg(imt), flxlft(imt), flxrgt(imt)
      real Rpl(imt,km), Rmn(imt,km), tmaski(imt,km,jmw)
#  if defined O_fct_3d
      real tmpext(imt,km,2)
#  endif
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
     &,         t_lo(imt,km)
#  endif
      include "fdift.h"

!-----------------------------------------------------------------------
!     limit the indices based on those from the argument list
!-----------------------------------------------------------------------

      istrt   = max(2,is)
      iend    = min(imt-1,ie)
      istrtm1 = istrt - 1
      iendp1  = iend + 1
      jstrt   = js
      jend    = min(je,jmt-1-joff)

!-----------------------------------------------------------------------
!     initialization when calculating jrow 2
!-----------------------------------------------------------------------

      if (joff + js .eq. 2) then
        jstrt = js - 1
        do k=1,km
          do i=istrt-1,iend
            adv_fn(i,k,1) = c0
            anti_fn(i,k,1,n) = c0
            R_plusY(i,k,1,n) = c0
            R_minusY(i,k,1,n) = c0
#  if defined O_fct_3d
            R_plus3(i,k,1,n) = c0
            R_minus3(i,k,1,n) = c0
#  endif
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     create an inverse land mask
!-----------------------------------------------------------------------

      do j=1,jmw
        do k=1,km
          do i=1,imt
            tmaski(i,k,j) = c1 - tmask(i,k,j)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective (low order scheme) flux across northern,
!     eastern and bottom faces of "T" cells
!-----------------------------------------------------------------------

      jlast = min(jend+1+joff,jmt-1) - joff
      do j=js-1,jlast
        do k=1,km
          do i=istrt,iend
            totadv = adv_vnt(i,k,j)
#  if defined O_gent_mcwilliams
     &             + adv_vntiso(i,k,j)
#  endif
            adv_fn(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i,k,j+1,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i,k,j+1,n,taum1))
          enddo
        enddo
      enddo

      jlast = min(jend+1+joff,jmt-1) - joff
      do j=js,jlast
        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,j)
#  if defined O_gent_mcwilliams
     &             + adv_vetiso(i,k,j)
#  endif
            adv_fe(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i+1,k,j,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i+1,k,j,n,taum1))
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,j)
#  if defined O_gent_mcwilliams
     &             + adv_vbtiso(i,k,j)
#  endif
            adv_fb(i,k,j) = totadv*
     &                      (t(i,k+1,j,n,taum1) + t(i,k,j,n,taum1))
     &                    + abs(totadv)*
     &                      (t(i,k+1,j,n,taum1) - t(i,k,j,n,taum1))
          enddo
        enddo
        do i=istrt,iend
          adv_fb(i,0,j) = adv_vbt(i,0,j)*c2*t(i,1,j,n,taum1)
          adv_fb(i,km,j) = c0
        enddo
      enddo

!-----------------------------------------------------------------------
!     main j loop
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = (j+1) + joff
        jp2  = min(j+2+joff,jmt) - joff
        jp1  = min(j+1+joff,jmt-1) - joff

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer at center of "T" cells in row j+1
!       - low order solution -
!-----------------------------------------------------------------------

        do k=1,km
          twodt(k) = c2dtts*dtxcel(k)
          do i=istrt,iend
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            t_lo(i,k) = (t(i,k,j+1,n,taum1) - twodt(k)
#  else
            t_lo(i,k,j+1,n) = (t(i,k,j+1,n,taum1) - twodt(k)
#  endif
     &        *(ADV_Tx(i,k,jp1) + ADV_Ty(i,k,jp1,jrow,n) +
     &        ADV_Tz(i,k,jp1))*tmask(i,k,j+1))
          enddo
        enddo
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
        call setbcx (t_lo(1,1), imt, km)
#  else
        call setbcx (t_lo(1,1,j+1,n), imt, km)
#  endif

!-----------------------------------------------------------------------
!       next calculate raw antidiffusive fluxes, that is high order
!       scheme flux (leap frog) minus the low order (upstream)
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,jp1)
#  if defined O_gent_mcwilliams
     &             + adv_vetiso(i,k,jp1)
#  endif
            anti_fe(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i+1,k,j+1,n,tau)) - adv_fe(i,k,jp1)
          enddo
          do i=istrt,iend
            totadv = adv_vnt(i,k,jp1)
#  if defined O_gent_mcwilliams
     &             + adv_vntiso(i,k,jp1)
#  endif
            anti_fn(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i,k,jp2,n,tau)) - adv_fn(i,k,jp1)
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,jp1)
#  if defined O_gent_mcwilliams
     &             + adv_vbtiso(i,k,jp1)
#  endif
            anti_fb(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i,k+1,j+1,n,tau)) - adv_fb(i,k,jp1)
     &                           *tmask(i,k,j+1)
          enddo
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n) = adv_vbt(i,0,j+1)*c2*t(i,1,j+1,n,taum1)
          anti_fb(i,km,j+1,n) = c0
        enddo

!-----------------------------------------------------------------------
!       now calculate and apply one-dimensional delimiters to these
!       raw antidiffusive fluxes

!       1) calculate T*, that are all halfway neighbors of T
!       2) calculate ratio R+- of Q+- to P+-, that is maximal/minimal
!          possible change of T if no limit would be active,
!          must be at least 1
!       3) choose correct ratio depending on direction of flow as a
!          delimiter
!       4) apply this delimiter to raw antidiffusive flux
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       delimit x-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function
#  if defined O_fct_dlm1 || !defined O_fct_dlm2

!         running mean of two adjacent points

          do i=istrt,iendp1
            Trmax(i) = p5*(t(i-1,k,j+1,n,tau) + t(i,k,j+1,n,tau))
          enddo
#  endif

!         extremum of low order solution central point and adjacent
!         halfway neighbours; check for land

          do i=istrt,iend
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            fxa = tmask(i-1,k,j+1)*Trmax(i) +
     &            tmaski(i-1,k,j+1)*t_lo(i,k)
            fxb = tmask(i+1,k,j+1)*Trmax(i+1) +
     &            tmaski(i+1,k,j+1)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
#  else
            fxa = tmask(i-1,k,j+1)*t_lo(i-1,k,j+1,n) +
     &            tmaski(i-1,k,j+1)*t_lo(i,k,j+1,n)
            fxb = tmask(i+1,k,j+1)*t_lo(i+1,k,j+1,n) +
     &            tmaski(i+1,k,j+1)*t_lo(i,k,j+1,n)
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
#  endif
#  if defined O_fct_3d
            tmpext(i,k,1) = Trmax(i)
            tmpext(i,k,2) = Trmin(i)
#  endif
            dcf(i) = cstdxt2r(i,j+1)
            flxlft(i) = anti_fe(i-1,k,j+1,n)
            flxrgt(i) = anti_fe(i,k,j+1,n)
          enddo

!         calculate ratio R

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            Qplus  = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
#  else
            Qplus  = Trmax(i) - t_lo(i,k,j+1,n)
            Qminus = t_lo(i,k,j+1,n) - Trmin(i)
#  endif
            Rpl(i,k) = min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            Rmn(i,k) = min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo
          call setbcx (Rpl, imt, km)
          call setbcx (Rmn, imt, km)

!         calculate delimiter using ratio at adjacent points

          do i=istrt,iendp1
            Cpos(i-1) = min(Rpl(i,k),Rmn(i-1,k))
            Cneg(i-1) = min(Rpl(i-1,k),Rmn(i,k))
          enddo

!         finally apply appropriate delimiter to flux

          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fe(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fe(i,k,j+1,n)))
          enddo
        enddo

!-----------------------------------------------------------------------
!       delimit y-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

          do i=istrt,iend
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            fxa = p5*tmask(i,k,j)*(t(i,k,j,n,tau) +
     &            t(i,k,j+1,n,tau)) +
     &            tmaski(i,k,j)*t_lo(i,k)
            fxb = p5*tmask(i,k,jp2)*(t(i,k,j+1,n,tau) +
     &            t(i,k,jp2,n,tau)) +
     &            tmaski(i,k,jp2)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
#  else
            fxa = tmask(i,k,j)*t_lo(i,k,j,n) +
     &            tmaski(i,k,j)*t_lo(i,k,j+1,n)
            fxb = p5*tmask(i,k,jp2)*(t(i,k,j+1,n,tau) +
     &            t(i,k,jp2,n,tau)) +
     &            tmaski(i,k,jp2)*t_lo(i,k,j+1,n)
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
#  endif
#  if defined O_fct_3d
            tmpext(i,k,1) = max(Trmax(i),tmpext(i,k,1))
            tmpext(i,k,2) = min(Trmin(i),tmpext(i,k,2))
#  endif
            dcf(i) = cstdyt2r(jrow)
            flxlft(i) = anti_fn(i,k,j,n)
            flxrgt(i) = anti_fn(i,k,j+1,n)
          enddo

!         calculate ratio R, related to a point

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            Qplus  = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
#  else
            Qplus  = Trmax(i) - t_lo(i,k,j+1,n)
            Qminus = t_lo(i,k,j+1,n) - Trmin(i)
#  endif
            R_plusY(i,k,j+1,n) =
     &        min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            R_minusY(i,k,j+1,n) =
     &        min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo

!         calculate delimiter using ratio at adjacent points

          do i=istrt,iend
            Cpos(i) = min(R_plusY(i,k,j+1,n),R_minusY(i,k,j,n))
            Cneg(i) = min(R_plusY(i,k,j,n),R_minusY(i,k,j+1,n))
          enddo

!         finally get delimiter c dependent on direction of flux and
!         apply it to raw antidiffusive flux

          do i=istrt,iend
            anti_fn(i,k,j,n) = p5*((Cpos(i) + Cneg(i))
     &                             *anti_fn(i,k,j,n) +
     &                             (Cpos(i) - Cneg(i))
     &                             *abs(anti_fn(i,k,j,n)))
          enddo
        enddo

!-----------------------------------------------------------------------
!       delimit z-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

          do i=istrt,iend
            dcf(i) = dzt2r(k)
            flxlft(i) = anti_fb(i,k,j+1,n)
            flxrgt(i) = anti_fb(i,k-1,j+1,n)
            if (k .gt. 1)then
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
              fxa = p5*tmask(i,k-1,j+1)*
     &              (t(i,k-1,j+1,n,tau) + t(i,k,j+1,n,tau)) +
     &              tmaski(i,k-1,j+1)*t_lo(i,k)
#  else
              fxa = tmask(i,k-1,j+1)*t_lo(i,k-1,j+1,n) +
     &              tmaski(i,k-1,j+1)*t_lo(i,k,j+1,n)
#  endif
            else
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
              fxa = t_lo(i,k)
#  else
              fxa = t_lo(i,k,j+1,n)
#  endif
            endif
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            if (k .lt. km) then
              fxb = p5*tmask(i,k+1,j+1)*
     &              (t(i,k,j+1,n,tau)+t(i,k+1,j+1,n,tau)) +
     &              tmaski(i,k+1,j+1)*t_lo(i,k)
#  else
            if (k .lt. km) then
              fxb = tmask(i,k+1,j+1)*t_lo(i,k+1,j+1,n) +
     &              tmaski(i,k+1,j+1)*t_lo(i,k,j+1,n)
#  endif
            else
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
              fxb = t_lo(i,k)
#  else
              fxb = t_lo(i,k,j+1,n)
#  endif
            endif
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
#  else
            Trmax(i) = max(fxa,fxb,t_lo(i,k,j+1,n))
            Trmin(i) = min(fxa,fxb,t_lo(i,k,j+1,n))
#  endif
#  if defined O_fct_3d
            tmpext(i,k,1) = max(Trmax(i),tmpext(i,k,1))
            tmpext(i,k,2) = min(Trmin(i),tmpext(i,k,2))
#  endif
          enddo

!         calculate delimiter using ratio at adjacent points
!         this variable is related to an arc (between two points,
!         the same way as fluxes are defined.)

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            Qplus = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
#  else
            Qplus = Trmax(i) - t_lo(i,k,j+1,n)
            Qminus = t_lo(i,k,j+1,n) - Trmin(i)
#  endif
            Rpl(i,k) = min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            Rmn(i,k) = min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo

        enddo

        do k=1,kmm1

!         calculate delimiter using ratio at adjacent points
!         this variable is related to an arc (between two points,
!         the same way as fluxes are defined.)

          do i=istrt,iend
            Cneg(i) = min(Rpl(i,k+1),Rmn(i,k))
            Cpos(i) = min(Rpl(i,k),Rmn(i,k+1))
          enddo

!         finally get delimiter c dependent on direction of flux and
!         apply it to raw antidiffusive flux

          do i=istrt,iend
            anti_fb(i,k,j+1,n) = p5*((Cpos(i)+Cneg(i))
     &                               *anti_fb(i,k,j+1,n) +
     &                               (Cpos(i)-Cneg(i))
     &                               *abs(anti_fb(i,k,j+1,n)))
          enddo
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n)  = c0
          anti_fb(i,km,j+1,n) = c0
        enddo
#  if defined O_fct_3d

!-----------------------------------------------------------------------
!     then calculate and apply 3-d delimiter to just delimited
!     antidiffusive fluxes
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

          do i=istrt,iend
            Trmax(i) = tmpext(i,k,1)
            Trmin(i) = tmpext(i,k,2)
          enddo

          do i=istrt,iend
#  if defined O_fct_dlm1 || !defined O_fct_dlm2
            Qplus  = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
#  else
            Qplus  = Trmax(i) - t_lo(i,k,j+1,n)
            Qminus = t_lo(i,k,j+1,n) - Trmin(i)
#  endif
            R_plus3(i,k,j+1,n) = min(1.,tmask(i,k,j+1)*Qplus/
     &                           (epsln+c2dtts*(
     &                           cstdxt2r(i,j+1)
     &                             *(max(c0,anti_fe(i-1,k,j+1,n)) -
     &                               min(c0,anti_fe(i,k,j+1,n))) +
     &                           cstdyt2r(jrow)
     &                             *(max(c0,anti_fn(i,k,j,n)) -
     &                               min(c0,anti_fn(i,k,j+1,n))) +
     &                           dzt2r(k)
     &                             *(max(c0,anti_fb(i,k,j+1,n)) -
     &                               min(c0,anti_fb(i,k-1,j+1,n)))
     &                           )))

            R_minus3(i,k,j+1,n) = min(1.,tmask(i,k,j+1)*Qminus/
     &                            (epsln+c2dtts*(
     &                            cstdxt2r(i,j+1)
     &                              *(max(c0,anti_fe(i,k,j+1,n)) -
     &                                min(c0,anti_fe(i-1,k,j+1,n))) +
     &                            cstdyt2r(jrow)
     &                              *(max(c0,anti_fn(i,k,j+1,n)) -
     &                                min(c0,anti_fn(i,k,j,n))) +
     &                            dzt2r(k)
     &                              *(max(c0,anti_fb(i,k-1,j+1,n)) -
     &                                min(c0,anti_fb(i,k,j+1,n)))
     &                            )))

          enddo
        enddo
        call setbcx (R_plus3(1,1,j+1,n), imt, km)
        call setbcx (R_minus3(1,1,j+1,n), imt, km)

!       finally apply 3-d delimiters to precorrected fluxes

        do k=1,km
          do i=istrt,iendp1
            Cpos(i-1) = min(R_plus3(i,k,j+1,n),R_minus3(i-1,k,j+1,n))
            Cneg(i-1) = min(R_plus3(i-1,k,j+1,n),R_minus3(i,k,j+1,n))
          enddo
          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fe(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fe(i,k,j+1,n)))
          enddo

          do i=istrt,iend
            Cpos(i) = min(R_plus3(i,k,j+1,n),R_minus3(i,k,j,n))
            Cneg(i) = min(R_plus3(i,k,j,n),R_minus3(i,k,j+1,n))
          enddo
          do i=istrt,iend
            anti_fn(i,k,j,n) = p5*((Cpos(i) + Cneg(i))
     &                             *anti_fn(i,k,j,n) +
     &                             (Cpos(i) - Cneg(i))
     &                             *abs(anti_fn(i,k,j,n)))
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            Cneg(i) = min(R_plus3(i,k+1,j+1,n),R_minus3(i,k,j+1,n))
            Cpos(i) = min(R_plus3(i,k,j+1,n),R_minus3(i,k+1,j+1,n))
          enddo
          do i=istrt,iend
            anti_fb(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fb(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fb(i,k,j+1,n)))
          enddo
        enddo

#  endif

!-----------------------------------------------------------------------
!       complete advective fluxes by adding low order fluxes to
!       delimited antidiffusive fluxes
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = anti_fe(i,k,j+1,n) + adv_fe(i,k,jp1)
          enddo
          do i=istrt,iend
            anti_fn(i,k,j,n) = (anti_fn(i,k,j,n) + adv_fn(i,k,j))
     &                         *tmask(i,k,j)
            anti_fb(i,k,j+1,n) = (anti_fb(i,k,j+1,n) + adv_fb(i,k,jp1))
     &                           *tmask(i,k,j+1)
          enddo
        enddo

      enddo

!-----------------------------------------------------------------------
!     set 2*corrected advective fluxes across northern, eastern and
!     bottom faces of "T" cells
!-----------------------------------------------------------------------

      do j=js-1,jend
        do k=1,km
          do i=istrt,iend
            adv_fn(i,k,j) = anti_fn(i,k,j,n)
          enddo
        enddo
      enddo

      do j=js,jend
        do k=1,km
          do i=istrtm1,iend
            adv_fe(i,k,j) = anti_fe(i,k,j,n)
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            adv_fb(i,k,j) = anti_fb(i,k,j,n)
          enddo
        enddo
      enddo

# else
!=======================================================================
!     2nd order advective tracer flux

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!        n   = tracer
!=======================================================================

      implicit none

      integer istrt, iend, j, js, je, k, i, n, joff, ie, is, ip, kr, jq

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      parameter (istrt=2, iend=imt-1)

!-----------------------------------------------------------------------
!     calculate 2*advective flux across eastern face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km
          do i=istrt-1,iend
            adv_fe(i,k,j) = adv_vet(i,k,j)*(t(i,  k,j,n,tau) +
     &                                      t(i+1,k,j,n,tau))
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective flux across bottom face of "T" cells.
!     (It`s done this way for performance issues)
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km-1
          do i=istrt,iend
            adv_fb(i,k,j)  = adv_vbt(i,k,j)*(t(i,k,  j,n,tau) +
     &                                       t(i,k+1,j,n,tau))
          enddo
        enddo
      enddo

# endif
#endif
      return
      end
