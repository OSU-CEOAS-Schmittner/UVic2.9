      subroutine xbti

#if defined O_mom && defined O_xbts
!=======================================================================
!     initialize all  XBT positions and arrays

!     inputs:

!     xbtlon = longiutde of XBT
!     xbtlat = latitude of XBT
!     xbtdpt = depth of XBT

!     outputs:

!     numxbt = number of XBTs
!     xname  = names of quantities measured by the XBTs
!     nxbtts = time step counter for averaging XBTs
!     ixbt   = nearest model grid point to xbtlon
!     jxbt   = nearest model grid point to xbtlat
!     kxbt   = nearest model grid point to xbtdpt
!     nsxbt  = starting XBT number for jrow
!     nexbt  = ending XBT number for jrow
!     axbt   = space for accumulating XBT data for averaging
!=======================================================================

      implicit none

      integer n, i, maxkm, indp, nsort, isort, jrow, num, k, m

      logical errors
      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "iounit.h"
      include "xbt.h"

      integer numb(maxxbt), numrow(jmt)

!-----------------------------------------------------------------------
!     initialize positions of XBT stations. there may be up to
!     "maxxbt" stations specified.

!      example: deploy 15 XBT`s at latitude 10S every 10 degs in the
!               pacific starting at 135E and drop them down to 4000m.
!               then deploy 15 more at the equator.

!     10S pacific deployment

!      n = 0
!      do  i=1,15
!        n = n + 1
!        xbtlon(n) = 135.0 + (i-1)*10.0
!         xbtlat(n) = -10.0
!         xbtdpt(n) = 4000.0e2
!      enddo

!     equatorial pacific deployment

!      do  i=1,15
!        n = n + 1
!        xbtlon(n) = 135.0 + (i-1)*10.0
!         xbtlat(n) = 0.0
!         xbtdpt(n) = 4000.0e2
!      enddo
!-----------------------------------------------------------------------

!     initialize all XBTs to "not set" condition

      do n=1,maxxbt
        xbtlon(n) = epsln
        xbtlat(n) = epsln
        xbtdpt(n) = epsln
      enddo

!-----------------------------------------------------------------------
!     begin user deployment of XBTs
!-----------------------------------------------------------------------

!     test deployment of 2 XBTs near equator to depth of near 40m

      n = 0
      do i=1,2
        n = n + 1
        xbtlon(n) = 180.0 + (i-1)*10.0
        xbtlat(n) = 0.0 - 3*i
        xbtdpt(n) = 40.0e2
      enddo

!-----------------------------------------------------------------------
!     end user deployment of XBTs
!-----------------------------------------------------------------------

!     verify that the number of XBTs doesn`t exceed the "maxxbt" limit

      if (n .gt. maxxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Error: ==> number of XBT stations (',n
     &,                  ') exceeds parameter "maxxbt"'
        stop "=>xbti"
      elseif (n .lt. maxxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Warning: ==> save space by reducing'
     &,                  ' parameter "maxxbt" to ',n,' in "xbt.h"'
      endif
      if (n .eq. 0) then
        write (stdout,*) ' '
        write (stdout,*) ' Error: ==> no XBT locations were specified.'
     &,                  ' set them in routine "xbti"'
        stop "=>xbti"
      endif

!-----------------------------------------------------------------------
!     specify names of quantities measured by the XBT
!-----------------------------------------------------------------------

      xnameu(1,1) = ' dU/dt      '
      xnameu(1,2) = ' dV/dt      '
      xnameu(2,1) = ' -(P)x      '
      xnameu(2,2) = ' -(P)y      '
      xnameu(3,1) = ' -(UU)x     '
      xnameu(3,2) = ' -(UV)x     '
      xnameu(4,1) = ' -(VU)y     '
      xnameu(4,2) = ' -(VV)y     '
      xnameu(5,1) = ' -(WU)z     '
      xnameu(5,2) = ' -(WV)z     '
      xnameu(6,1) = ' DIFF_Ux    '
      xnameu(6,2) = ' DIFF_Vx    '
      xnameu(7,1) = ' DIFF_Uy    '
      xnameu(7,2) = ' DIFF_Vy    '
      xnameu(8,1) = ' DIFF_Uz    '
      xnameu(8,2) = ' DIFF_Vz    '
      xnameu(9,1) = ' DIFF_Umet  '
      xnameu(9,2) = ' DIFF_Vmet'
      xnameu(10,1) = ' fV         '
      xnameu(10,2) = ' -fU        '
      xnameu(11,1) = '   source   '
      xnameu(11,2) = '   source   '
      xnameu(12,1) = ' -(surf P)x '
      xnameu(12,2) = ' -(surf P)y '
      xnameu(13,1) = 'ADV_Umetric '
      xnameu(13,2) = '-ADV_Vmetric'
      xnameu(14,1) = ' -U(U)x     '
      xnameu(14,2) = ' -U(V)x     '
      xnameu(15,1) = ' -V(U)y     '
      xnameu(15,2) = ' -V(V)y     '
      xnameu(16,1) = ' -W(U)z     '
      xnameu(16,2) = ' -W(V)z     '
      xnameu(17,1) = ' average U  '
      xnameu(17,2) = ' average V  '

      xnamex(1) = ' Surf T flux'
      xnamex(2) = ' Taux       '
      xnamex(3) = ' Tauy       '
      xnamex(4) = ' average W  '

      xnamet(1) = ' dT/dt      '
      xnamet(2) = ' -(UT)x     '
      xnamet(3) = ' -(VT)y     '
      xnamet(4) = ' -(WT)z     '
      xnamet(5) = ' DIFF_Tx    '
      xnamet(6) = ' DIFF_Ty    '
      xnamet(7) = ' DIFF_Tz    '
      xnamet(8) = ' T source   '
      xnamet(9) = ' T convect  '
      xnamet(10) = ' T filter   '
      xnamet(11) = ' -U(T)x     '
      xnamet(12) = ' -V(T)y     '
      xnamet(13) = ' -W(T)z     '
      xnamet(14) = ' chg var  T '
      xnamet(15) = ' average T  '

!-----------------------------------------------------------------------
!     initialize counter for the number of time steps per average
!-----------------------------------------------------------------------

      nxbtts = 0

!-----------------------------------------------------------------------
!     convert XBT positions to nearest model grid points and set
!     unused XBTs to arbitrary value. don`t allow more than "kmxbt"
!     levels to be sampled at any position.
!-----------------------------------------------------------------------

      errors = .false.
      maxkm = 1
      do n=1,maxxbt
        if (xbtlon(n) .ne. epsln) then
          numxbt  = n
          ixbt(n) = indp(xbtlon(n), xt(2), imt-2) + 1
          jxbt(n) = indp(xbtlat(n), yt(2), jmt-2) + 1
          kxbt(n) = indp(xbtdpt(n), zt, km)
          maxkm   = max(kxbt(n),maxkm)
          numb(n) = (jxbt(n)-1)*imt + ixbt(n)
          if (kxbt(n) .gt. kmxbt) then
            errors = .true.
            write (stdout,'(/a,i4,a,f10.2,a/a,i3,a,i3,a,i3/)')
     &      '=> Error: XBT station #',n, ' is set for a depth of '
     &,                xbtdpt(n),' cm in xbti.F'
     &,     '          The nearest model level is k=',kxbt(n)
     &,                '. Increase "kmxbt" in xbt.h from ', kmxbt
     &,                ' to ',kxbt(n)
          endif
        else
          ixbt(n) = 0
          jxbt(n) = 0
          kxbt(n) = 0
          numb(n) = 0
        endif
      enddo
      if (maxkm .lt. kmxbt) then
        write (stdout,*) ' '
        write (stdout,*) ' Warning: ==> save space by reducing'
     &,                  ' parameter "kmxbt" to ',maxkm,' in "xbt.h"'
      endif
      if (errors) then
        stop '=>xbt.F'
      endif

!-----------------------------------------------------------------------
!     sort the XBTs from south to north and west to east
!-----------------------------------------------------------------------

      do nsort=1,1000
        isort = 0
        do n=2,numxbt
          if (numb(n) .lt. numb(n-1)) then
            isort = 1
            call iswapx (ixbt(n), ixbt(n-1))
            call iswapx (jxbt(n), jxbt(n-1))
            call iswapx (kxbt(n), kxbt(n-1))
            call iswapx (numb(n), numb(n-1))
          endif
        enddo
        if (isort .eq. 0) go to 10
      enddo
10    continue
      write (stdout,8900)
      do n=1,numxbt
        write (stdout,9000) n, yt(jxbt(n)), xt(ixbt(n))
     &,                       kxbt(n), zt(kxbt(n))*0.01
      enddo

!-----------------------------------------------------------------------
!     count the number of XBTs on each model latitude.
!     nsxbt is the starting XBT number
!     nexbt is the ending XBT number
!-----------------------------------------------------------------------

      do jrow=1,jmt
        numrow(jrow) = 0
      enddo
      do n=1,numxbt
        numrow(jxbt(n)) = numrow(jxbt(n)) + 1
      enddo
      n = 0
      do jrow=1,jmt
        if (numrow(jrow) .ne. 0) then
          nsxbt(jrow) = n + 1
          nexbt(jrow) = nsxbt(jrow) + numrow(jrow) - 1
          n           = n + numrow(jrow)
        else
          nsxbt(jrow) = 0
          nexbt(jrow) = 0
        endif
      enddo

!-----------------------------------------------------------------------
!     initialize all time accumulators to zero
!-----------------------------------------------------------------------

      do num=1,numxbt
        do k=1,kmxbt
          do n=1,nt
            do m=1,ntxbt
              txbt(k,m,n,num) = c0
            enddo
          enddo
          do n=1,2
            do m=1,nuxbt
              uxbt(k,m,n,num) = c0
            enddo
          enddo
          xbtw(k,num) = c0
        enddo
        uxbtsf(1,num) = c0
        uxbtsf(2,num) = c0
        do n=1,nt
          txbtsf(n,num) = c0
        enddo
      enddo

      return
8900  format (//,20x,' X B T    S T A T I O N    L O C A T I O N S'/)
9000  format (1x, ' XBT station #',i4, ' is at lat =',f6.2,', lon ='
     &,      f6.2, ', for ',i3,' levels down to a depth of ',f6.1,' m')
      end
      subroutine iswapx (i, j)

      implicit none

      integer i, j, itemp

      itemp = i
      i     = j
      j     = itemp
      return
      end

      subroutine txbt1 (joff, js, je, n)

!=======================================================================
!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       n     = tracer component
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, nth, lev

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, term, r2dt, dudx
      real dvdy, dwdz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "xbt.h"
      include "scalar.h"
      include "vmixc.h"

# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "fdift.h"

      if (n .gt. nt) then
        write (stdout,*) '=> Error: n=',n,' in txbt1.F'
        stop '=>txbt1'
      endif

      do j=js,je
        jrow = j + joff
        if (nsxbt(jrow) .ne. 0) then

!-----------------------------------------------------------------------
!         accumulate data for the nth XBT for one time step
!-----------------------------------------------------------------------

          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
            lev = kxbt(nth)
            do k=1,lev

!-----------------------------------------------------------------------
!             tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*t(i,k,j,n,tau)
              txbt(k,15,n,nth) = txbt(k,15,n,nth) + term

!-----------------------------------------------------------------------
!             d(tracer)/dt
!-----------------------------------------------------------------------

              r2dt = c1/(c2dtts*dtxcel(k))
              term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                             t(i,k,j,n,taum1))*r2dt
              txbt(k,9,n,nth) = txbt(k,9,n,nth) + term

!-----------------------------------------------------------------------
!             zonal advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Tx(i,k,j)
# if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Txiso(i,k,j,n)
# endif
              txbt(k,2,n,nth) = txbt(k,2,n,nth) + term

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
              txbt(k,11,n,nth) = txbt(k,11,n,nth) + term

!-----------------------------------------------------------------------
!             meridional advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Ty(i,k,j,jrow,n)
# if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Tyiso(i,k,j,jrow,n)
# endif
              txbt(k,3,n,nth) = txbt(k,3,n,nth) + term

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
              txbt(k,12,n,nth) = txbt(k,12,n,nth) + term

!-----------------------------------------------------------------------
!             vertical advection (flux form) of tracer
!-----------------------------------------------------------------------

              term = -tmask(i,k,j)*ADV_Tz(i,k,j)
#  if defined O_gent_mcwilliams && !defined O_fct
     &               -tmask(i,k,j)*ADV_Tziso(i,k,j)
#  endif
              txbt(k,4,n,nth) = txbt(k,4,n,nth) + term

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
              txbt(k,13,n,nth) = txbt(k,13,n,nth) + term

!-----------------------------------------------------------------------
!             zonal diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Tx(i,k,j)
              txbt(k,5,n,nth) = txbt(k,5,n,nth) + term

!-----------------------------------------------------------------------
!             meridional diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Ty(i,k,j,jrow,n)
              txbt(k,6,n,nth) = txbt(k,6,n,nth) + term

!-----------------------------------------------------------------------
!             vertical diffusion of tracer
!-----------------------------------------------------------------------

              term = tmask(i,k,j)*DIFF_Tz(i,k,j)
# if defined O_implicitvmix || defined O_isopycmix
     &               + tmask(i,k,j)*zzi(i,k,j)
# endif
              txbt(k,7,n,nth) = txbt(k,7,n,nth) + term

!-----------------------------------------------------------------------
!             tracer source term
!-----------------------------------------------------------------------

# if defined O_source_term || defined O_npzd || defined O_carbon_14
              term = tmask(i,k,j)*source(i,k,j)
              txbt(k,8,n,nth) = txbt(k,8,n,nth) + term
# else
              txbt(k,8,n,nth) = 0.0
# endif

            enddo

!-----------------------------------------------------------------------
!           surface tracer flux
!-----------------------------------------------------------------------

            k    = 1
            term = tmask(i,k,j)*stf(i,j,n)
            txbtsf(n,nth) = txbtsf(n,nth) + term
          enddo
        endif
      enddo

      return
      end

      subroutine txbt2 (joff, js, je, iterm)

!=======================================================================
!     Accumulate parts of d(tracer)/dt and change in tracer variance

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       iterm = 1  => total change
!       iterm = 10 => change due to filtering
!=======================================================================

      implicit none

      integer iterm, j, js, je, jrow, joff, nth, i, lev, n, k

      real r2dt, term

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "mw.h"
      include "scalar.h"
      include "xbt.h"

      if (iterm .ne. 1 .and. iterm .ne. 10) then
        write (stdout,*) '=>Error: iterm=',iterm,' in txbt2'
        stop '=>txbt2'
      endif

      do j=js,je
        jrow  = j + joff
        if (nsxbt(jrow) .ne. 0) then
          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
            lev = kxbt(nth)
            do n=1,nt
              do k=1,lev

!-----------------------------------------------------------------------
!               d/dt(tracer)
!-----------------------------------------------------------------------

                r2dt = c1/(c2dtts*dtxcel(k))
                term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                               t(i,k,j,n,taum1))*r2dt
                txbt(k,iterm,n,nth) = txbt(k,iterm,n,nth) + term

!-----------------------------------------------------------------------
!               change in variance of tracer
!-----------------------------------------------------------------------

                if (iterm .eq. 1) then
                  term = tmask(i,k,j)*(t(i,k,j,n,taup1)**2-
     &                                 t(i,k,j,n,taum1)**2)
                  txbt(k,14,n,nth) = txbt(k,14,n,nth) + term
                endif
              enddo
            enddo
          enddo
        endif
      enddo

      return
      end

      subroutine uxbt1 (joff, js, je, n)

!=======================================================================

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       n    = (1,2) for (u,v) component of velocity
!=======================================================================

      implicit none

      integer i, k, jrow, j, n, js, je, joff, nth, lev

      real adv_ux, adv_uy, adv_uz, adv_metric, diff_ux, diff_uz
      real diff_uy, diff_metric, coriolis, term, dudx, dvdy, dwdz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "xbt.h"
      include "vmixc.h"
      include "fdifm.h"

      do j=js,je
        jrow = j + joff
        if (nsxbt(jrow) .ne. 0) then

!         only allow n = 1 (u component) or 2 (v component)

          if (n .gt. 2) then
            write (stdout,*) '=> Error: n=',n,' in uxbt1.F'
            stop '=>uxbt1'
          endif

!-----------------------------------------------------------------------
!         accumulate data for the nth XBT for one time step
!-----------------------------------------------------------------------

          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
            lev = kxbt(nth)
            do k=1,lev

!-----------------------------------------------------------------------
!             pressure term
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*grad_p(i,k,j,n)
              uxbt(k,2,n,nth) = uxbt(k,2,n,nth) + term

!-----------------------------------------------------------------------
!             zonal advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Ux(i,k,j)
              uxbt(k,3,n,nth) = uxbt(k,3,n,nth) + term

!-----------------------------------------------------------------------
!             pure zonal advection of momentum
!-----------------------------------------------------------------------

!             - U(U)x = U(U)x - (UU)x (when n=1)
!             - U(V)x = V(U)x - (UV)x (when n=2)

              dudx = (adv_veu(i,k,j)-adv_veu(i-1,k,j))*dxur(i)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dudx - ADV_Ux(i,k,j))
              uxbt(k,14,n,nth) = uxbt(k,14,n,nth) + term

!-----------------------------------------------------------------------
!             advective metric term
!-----------------------------------------------------------------------

              term = ADV_metric(i,k,j,jrow,n)
              uxbt(k,13,n,nth) = uxbt(k,13,n,nth) + term

!-----------------------------------------------------------------------
!             meridional advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Uy(i,k,j,jrow,n)
              uxbt(k,4,n,nth) = uxbt(k,4,n,nth) + term

!-----------------------------------------------------------------------
!             pure meridional advection of momentum
!-----------------------------------------------------------------------

!             - V(U)y = U(V)y - (VU)y (when n=1)
!             - V(V)y = V(V)y - (VV)y (when n=2)

              dvdy = (adv_vnu(i,k,j)-adv_vnu(i,k,j-1))*dyur(jrow)
     &               *csur(jrow)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dvdy
     &             - ADV_Uy(i,k,j,jrow,n))
              uxbt(k,15,n,nth) = uxbt(k,15,n,nth) + term

!-----------------------------------------------------------------------
!             vertical advection (flux form) of momentum
!-----------------------------------------------------------------------

              term = -umask(i,k,j)*ADV_Uz(i,k,j)
              uxbt(k,5,n,nth) = uxbt(k,5,n,nth) + term

!-----------------------------------------------------------------------
!             pure vertical advection of momentum
!-----------------------------------------------------------------------

!             - W(U)z = U(W)z - (WU)z (when n=1)
!             - W(V)z = V(W)z - (WV)z (when n=2)

              dwdz = (adv_vbu(i,k-1,j)-adv_vbu(i,k,j))*dztr(k)
              term = umask(i,k,j)*(u(i,k,j,n,tau)*dwdz - ADV_Uz(i,k,j))
              uxbt(k,16,n,nth) = uxbt(k,16,n,nth) + term

!-----------------------------------------------------------------------
!             zonal diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Ux(i,k,j)
              uxbt(k,6,n,nth) = uxbt(k,6,n,nth) + term

!-----------------------------------------------------------------------
!             meridional diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Uy(i,k,j,jrow,n)
              uxbt(k,7,n,nth) = uxbt(k,7,n,nth) + term

!-----------------------------------------------------------------------
!             diffusive metric term
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_metric(i,k,j,jrow,n)
              uxbt(k,9,n,nth) = uxbt(k,9,n,nth) + term

!-----------------------------------------------------------------------
!             vertical diffusion of momentum
!-----------------------------------------------------------------------

              term = umask(i,k,j)*DIFF_Uz(i,k,j)
# if defined O_implicitvmix
     &              +umask(i,k,j)*zzi(i,k,j)
# endif
              uxbt(k,8,n,nth) = uxbt(k,8,n,nth) + term

!-----------------------------------------------------------------------
!             coriolis term
!-----------------------------------------------------------------------

              term = umask(i,k,j)*CORIOLIS(i,k,j,jrow,n)
              uxbt(k,10,n,nth) = uxbt(k,10,n,nth) + term

!-----------------------------------------------------------------------
!             accumulate the source term
!-----------------------------------------------------------------------

# if defined O_source_term || defined O_npzd || defined O_carbon_14
              term = umask(i,k,j)*source(i,k,j)
              uxbt(k,11,n,nth) = uxbt(k,11,n,nth) + term
# else
              uxbt(k,11,n,nth) = 0.0
# endif

!-----------------------------------------------------------------------
!             accumulate u, v, and w
!-----------------------------------------------------------------------

              term = umask(i,k,j)*u(i,k,j,n,tau)
              uxbt(k,17,n,nth) = uxbt(k,17,n,nth) + term

              if (n .eq. 2) then
                term = p5*(adv_vbu(i,k,j)+adv_vbu(i,k-1,j))*umask(i,k,j)
                xbtw(k,nth) = xbtw(k,nth) + term
              endif
            enddo

!-----------------------------------------------------------------------
!           accumulate the surface momentum flux
!-----------------------------------------------------------------------

            k    = 1
            term = umask(i,k,j)*smf(i,j,n)
            uxbtsf(n,nth) = uxbtsf(n,nth) + term

!-----------------------------------------------------------------------
!           update accumulation counter once per time step
!-----------------------------------------------------------------------

            if (nth .eq. numxbt .and. (n .eq. 2))  nxbtts = nxbtts + 1
          enddo
        endif
      enddo

      return
      end

      subroutine uxbt2 (joff, js, je, c2dtuv, acor)

!=======================================================================
!     Accumulate  d/dt and the implicit coriolis terms in the
!     momentum equations for XBTs

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       c2dtuv = (2*dtuv,dtuv) on (lpfrod,mixing) time steps
!       acor   = implicit factor
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, nth, i, lev, n, k

      real r2dt, c2dtuv, term, acor

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"
      include "xbt.h"

!-----------------------------------------------------------------------
!     local constants
!-----------------------------------------------------------------------

      r2dt = c1/c2dtuv

!-----------------------------------------------------------------------
!     d/dt of velocity (external mode part of tau+1 will be added
!     later when the external mode is solved)
!-----------------------------------------------------------------------

      do j=js,je
        jrow  = j + joff
        if (nsxbt(jrow) .ne. 0) then
          do nth=nsxbt(jrow),nexbt(jrow)
            i   = ixbt(nth)
            lev = kxbt(nth)
            do n=1,2
              do k=1,lev

!-----------------------------------------------------------------------
!               d/dt of velocity (external mode part of tau+1 will be
!               added later when the external mode is solved)
!-----------------------------------------------------------------------

                term = umask(i,k,j)*(u(i,k,j,n,taup1) -
     &                               u(i,k,j,n,taum1))*r2dt
                uxbt(k,1,n,nth) = uxbt(k,1,n,nth) + term

!-----------------------------------------------------------------------
!               implicit coriolis term (external mode part will be
!               added later when external mode is solved)
!-----------------------------------------------------------------------

                if (acor .ne. c0) then
                  term = umask(i,k,j)*acor*cori(i,jrow,n)*
     &                     (u(i,k,j,3-n,taup1) - u(i,k,j,3-n,taum1))
                  uxbt(k,10,n,nth) = uxbt(k,10,n,nth) + term
                endif
              enddo
            enddo
          enddo
        endif
      enddo

      return
      end

      subroutine uxbt3

!=======================================================================
!     Accumulate surface pressure gradients, external mode
!     part of d/dt and external moe part of implicit coriolis term
!     for XBTs
!=======================================================================

      implicit none

      integer is, ie, js, je, n, i, jrow, lev, kz, kl, k
      parameter (is=1, ie=1, js=1, je=1)

      real r2dtuv, uext, vext, atosp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "xbt.h"

      real psgrad(is:ie,js:je,2)

      r2dtuv = c1/c2dtuv
      do n=1,numxbt
        i      = ixbt(n)
        jrow   = jxbt(n)
        lev    = kxbt(n)
        kz     = kmu(i,jrow)
        if (kz .ne. 0) then

!         construct the surface pressure gradients for pt (i,jrow)
!         note: it will be stored in psgrad(is,js,:)

          call calc_psgrad (psgrad, uext, vext, jrow, jrow, i, i)
          atosp = acor*cori(i,jrow,1)
          kl = min(kz,lev)
          do k=1,kl
            uxbt(k,1,1,n)  = uxbt(k,1,1,n)  + uext*r2dtuv
            uxbt(k,1,2,n)  = uxbt(k,1,2,n)  + vext*r2dtuv
            uxbt(k,12,1,n) = uxbt(k,12,1,n) - psgrad(is,js,1)
            uxbt(k,12,2,n) = uxbt(k,12,2,n) - psgrad(is,js,2)
            uxbt(k,10,1,n) = uxbt(k,10,1,n) + vext*atosp
            uxbt(k,10,2,n) = uxbt(k,10,2,n) - uext*atosp
          enddo
        endif
      enddo

      return
      end

      subroutine xbto

!=======================================================================
!     XBT output: average and save all XBT data
!=======================================================================

      implicit none

      integer io, n, i, jrow, lev, k, m, lll

      real reltim, period, rnavg, tconv, tfilt, erru, errv, errt(nt)

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "iounit.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "xbt.h"

!     reltim defines the end of the averaged data in years

      reltim = prelyr
      period = dtts*nxbtts

      if (ioxbt .ne. stdout .or. ioxbt .lt. 0) then

        call getunit (io, 'xbt.dta'
     &,               'unformatted sequential append ieee')

        write (stdout,9999) numxbt, itt, period*secday, stamp
        iotext ='read (ioxbt) reltim, period, numxbt, ntxbt, nuxbt'
        write (io) pstamp, iotext, expnam
        write (io) reltim, period, numxbt, ntxbt, nuxbt

        iotext ='read (ioxbt) (xnamet(1:12)(n),n=1,ntxbt)'
        write (io) pstamp, iotext, expnam
        write (io) xnamet

        iotext ='read (ioxbt) ((xnameu(1:12)(n,m),n=1,nuxbt),m=1,2)'
        write (io) pstamp, iotext, expnam
        write (io) xnameu

        iotext ='read (ioxbt) (xnamex(1:12)(n),n=1,4)'
        write (io) pstamp, iotext, expnam
        write (io) xnamex
      endif

      if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
        write (stdout,'(//,20x,a,/)')
     &  'A V E R A G E D    X B T   S T A T I O N S'
        write (stdout,8800) numxbt, itt, period*secday, stamp
      endif

      do n=1,numxbt
        i     = ixbt(n)
        jrow  = jxbt(n)
        lev   = kxbt(n)

        if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
          write (stdout,8900)  n, yt(jrow), xt(i)
        endif

        if (ioxbt .ne. stdout .or. ioxbt .lt. 0) then
          iotext ='read (ioxbt) lev, rlon, rlat, tdx, tdy, udx, udy'
          write (io) pstamp, iotext, expnam
          write (io) lev, xt(i), yt(jrow), dxt(i)*cst(jrow), dyt(jrow)
     &,               dxu(i)*csu(jrow), dyu(jrow)

          iotext ='read (ioxbt) (dzt(k)=1,lev)'
          write (io) pstamp, iotext, expnam
          write (io) (dzt(k),k=1,lev)
        endif

!-----------------------------------------------------------------------
!       construct time mean quantities
!-----------------------------------------------------------------------

        rnavg = c1/nxbtts

!-----------------------------------------------------------------------
!       average the data for xbt number n, write it out, then zero
!       the accumulators
!-----------------------------------------------------------------------

        do i=1,nt
          do m=1,ntxbt
            do k=1,lev
              txbt(k,m,i,n) = rnavg*txbt(k,m,i,n)
            enddo
          enddo
          txbtsf(i,n) = rnavg*txbtsf(i,n)
        enddo
        do i=1,2
          do m=1,nuxbt
            do k=1,lev
              uxbt(k,m,i,n) = rnavg*uxbt(k,m,i,n)
            enddo
          enddo
          uxbtsf(i,n) = rnavg*uxbtsf(i,n)
        enddo
        do k=1,lev
          xbtw(k,n) = rnavg*xbtw(k,n)
        enddo

!       construct change due to convection and filtering

        do i=1,nt
          do k=1,lev
            tconv = txbt(k,10,i,n) - txbt(k,9,i,n)
            tfilt = txbt(k,1,i,n) - txbt(k,10,i,n)
            txbt(k,9,i,n)  = tconv
            txbt(k,10,i,n) = tfilt
          enddo
        enddo

!-----------------------------------------------------------------------
!       write out the results
!-----------------------------------------------------------------------

        if (ioxbt .eq. stdout .or. ioxbt .lt. 0) then
          do k=1,lev
            write (stdout,8700) k, zt(k)*0.01
            do m=1,nuxbt
              write (stdout,9000) xnameu(m,1), uxbt(k,m,1,n)
     &,                           xnameu(m,2), uxbt(k,m,2,n)
            enddo
            erru = c0
            errv = c0
            do lll=2,13
              erru = erru + uxbt(k,lll,1,n)
              errv = errv + uxbt(k,lll,2,n)
            enddo
            erru = uxbt(k,1,1,n) - erru
            errv = uxbt(k,1,2,n) - errv
            write (stdout,9000) ' error      ', erru
     &,                         ' error      ', errv
            if (k .eq. 1) then
              write (stdout,9000) xnamex(2),uxbtsf(1,n)
     &,                           xnamex(3),uxbtsf(2,n)
            endif
            write (stdout,*) ' '
            write (stdout,*) '           ',xnamex(4),' = ',xbtw(k,n)

            write (stdout,9001) (i,i=1,nt)
            do m=1,ntxbt
              write (stdout,9002) xnamet(m), (txbt(k,m,i,n),i=1,nt)
            enddo
            do i=1,nt
              errt(i) = c0
              do lll=2,10
                errt(i) = errt(i) + txbt(k,lll,i,n)
              enddo
              errt(i) = errt(i) - txbt(k,1,i,n)
            enddo
            write (stdout,9002) '  error     ', (errt(i),i=1,nt)
            if (k .eq. 1) then
              write (stdout,9002) xnamex(1), (uxbtsf(i,n),i=1,nt)
            endif
          enddo
        endif

        if (ioxbt .ne. stdout .or. ioxbt .lt. 0) then
          iotext =
     &   'read (ioxbt) (((uxbt(k,m,i),k=1,lev),m=1,nuxbt),i=1,2)'
          write (io) pstamp, iotext, expnam
          write (io) (((uxbt(k,m,i,n),k=1,lev),m=1,nuxbt),i=1,2)

          iotext = 'read (ioxbt) (xbtw(k),k=1,lev)'
          write (io) pstamp, iotext, expnam
          write (io) (xbtw(k,n),k=1,lev)

          iotext = 'read (ioxbt) (uxbtsf(i),i=1,2)'
          write (io) pstamp, iotext, expnam
          write (io) (uxbtsf(i,n),i=1,2)

          iotext =
     &    'read (ioxbt) (((txbt(k,m,i),k=1,lev),m=1,ntxbt),i=1,nt)'
          write (io) pstamp, iotext, expnam
          write (io) (((txbt(k,m,i,n),k=1,lev),m=1,ntxbt),i=1,nt)

          iotext = 'read (ioxbt) (txbtsf(i),i=1,nt)'
          write (io) pstamp, iotext, expnam
          write (io) (txbtsf(i,n),i=1,nt)
        endif

      enddo

!-----------------------------------------------------------------------
!       zero the accumulators
!-----------------------------------------------------------------------

      do n=1,numxbt
        do i=1,nt
          do m=1,ntxbt
            do k=1,lev
              txbt(k,m,i,n) = c0
            enddo
          enddo
          txbtsf(i,n) = c0
        enddo
        do i=1,2
          do m=1,nuxbt
            do k=1,lev
              uxbt(k,m,i,n) = c0
            enddo
          enddo
          uxbtsf(i,n) = c0
        enddo
        do k=1,lev
          xbtw(k,n) = c0
        enddo
      enddo

!-----------------------------------------------------------------------
!     zero the "averaging" counter for the next averaging period
!-----------------------------------------------------------------------

      nxbtts = 0

      if (ioxbt .ne. stdout .or. ioxbt .lt. 0) then
        call relunit (io)
      endif

8700  format (/,1x,' at model level ',i3,': depth = '
     &,f8.2,'m'/)
8800  format (//1x,'===> XBT diagnostic: ',i4, ' stations at ts = ',i7
     &,', averaged over ', f8.2,' days ending on ',a32)
8900  format (//1x,' XBT station #',i4,' location:  lat =',f6.2
     &,       ', lon =',f6.2)
9000  format (1x, a12,' = ',1pe14.7,2x, a12,' = ',1pe14.7,2x
     &,           a12,' = ',1pe14.7,2x, a12,' = ',1pe14.7,2x)
9001  format (/1x,14x, 3x,'tracer # ',i2,2x, 3x,'tracer # ',i2,2x
     &,                3x,'tracer # ',i2,2x, 3x,'tracer # ',i2,2x)
9002  format (1x,a12, 2x,1pe14.7, 2x,1e14.7, 2x,1e14.7, 2x,1e14.7)
9999  format (/1x
     & ,'===> ',i4,' XBT"s written to file xbt.dta on ts = ',i10
     &,        ', averaged over = ',f7.2, 'days ', a,/)
#endif

      return
      end
