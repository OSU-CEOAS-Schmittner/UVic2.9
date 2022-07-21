      subroutine stabi

#if defined O_mom && defined O_stability_tests
!-----------------------------------------------------------------------
!     initialization for stability testing
!-----------------------------------------------------------------------

      implicit none

      include "stab.h"
      tdig     = 1.e-4
      cflcrt   = 1.5
      maxcfl   = 3
      cflons   = 0.0
      cflone   = 360.0
      cflats   = -90.0
      cflate   = 90.0
      cfldps   = 0.0
      cfldpe   = 6000.0e2
      return
      end

      subroutine stab (j, jrow)

!-----------------------------------------------------------------------
!     test for various measures of stability
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, jrow, is, ie, ks, ke, m, kp1, kk
      integer n, jj, ii, iobadt, iobads

      real cl, cosur, dtmax, f1, f2, f3, cflu, cflv, cflwu, cflwt
      real umax, pcflu, vmax, pcflv, wmax, pcflwu, pcflwt, scl, fx
      real dtdz, ramb, rahb, rame, ramn, rahe, rahn, reyx, reyy
      real reyz, pecx, pecy, pecz, tbig, tsml, tcrit, tlocal, slocal

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "docnam.h"
      include "grdvar.h"
      include "hmixc.h"
      include "iounit.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "stab.h"
      include "state.h"
      include "tmngr.h"
      include "vmixc.h"

!-----------------------------------------------------------------------
!     perform stability tests only within specified lat,lon and depth
!     region
!-----------------------------------------------------------------------

      if (jrow .lt. jscfl .or. jrow .gt. jecfl) return

!-----------------------------------------------------------------------
!     scan for CFL violations. save locations of closest approach to
!     local CFL limit.
!-----------------------------------------------------------------------

      cl = cflcrt * p5
      cosur = max(csur(jrow),epsln)
      do k=kscfl,kecfl
        dtmax = max(dtuv, dtts*dtxcel(k))
        f1    = dtmax*dyur(jrow)
        f2    = dtmax*dzwr(k)
        f3    = dtmax*cosur
        do i=iscfl,iecfl
          cflu  = abs(f3*dxur(i)*u(i,k,j,1,tau))
          cflv  = abs(f1*u(i,k,j,2,tau))
          cflwu = abs(f2*adv_vbu(i,k,j))*umask(i,k,j)
          cflwt = abs(f2*adv_vbt(i,k,j))*tmask(i,k,j)
          if (cflu .ge. cl .or. cflv .ge. cl .or. cflwu .ge. cl .or.
     &        cflwt .ge. cl) then
            write (stdout,'(/,a33,i4,a1,i3,a1,i3,a13,f6.3,/)')
     &       ' ==> CFL exceeded at coordinate (i,j,k) = (',i
     &,      ',',jrow,',',k,') by factor =',cflcrt
            umax  = p5*csu(jrow)*dxu(i)/dtmax
            pcflu = abs(100.0*u(i,k,j,1,tau)/umax)
            vmax  = p5*dyu(jrow)/dtmax
            pcflv = abs(100.0*u(i,k,j,2,tau)/vmax)
            wmax  = p5*dzw(k)/dtmax
            pcflwu= abs(100.0*adv_vbu(i,k,j)/wmax)
            pcflwt= abs(100.0*adv_vbt(i,k,j)/wmax)
            write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' u reached   ', pcflu,' % of the CFL limit (',umax,')'
            write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' v reached   ', pcflv,' % of the CFL limit (',vmax,')'
            write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' adv_vbu reached', pcflwu,' % of the CFL limit (',wmax,')'
            write (stdout,'(a,f8.2,a,g15.8,a)')
     &       ' adv_vbt reached', pcflwt,' % of the CFL limit (',wmax,')'
            is = max(1,i-3)
            ie = min(imt,i+3)
            ks = max(1,k-3)
            ke = min(km,k+3)
            scl = c0
            fx  = 1.0e-2
            write (stdout,9100)'u velocity', itt
     &,      jrow, yu(jrow), xu(is), xu(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (u(1,1,j,1,tau), imt, is, ie, ks, ke, scl)

            write (stdout,9100) 'v velocity', itt
     &,      jrow, yu(jrow), xu(is), xu(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (u(1,1,j,2,tau), imt, is, ie, ks, ke, scl)

            write (stdout,9100)  'adv_vbu ', itt
     &,      jrow, yu(jrow), xu(is), xu(ie), fx*zw(ks), fx*zw(ke), scl
            call matrix (adv_vbu(1,1,j), imt, is, ie, ks, ke, scl)

            write (stdout,9100)  'adv_vbt ', itt
     &,      jrow, yt(jrow), xt(is), xt(ie), fx*zw(ks), fx*zw(ke), scl
            call matrix (adv_vbt(1,1,j), imt, is, ie, ks, ke, scl)

            do m=1,nt
              write (stdout,9100) trname(m), itt
     &,        jrow, yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (t(1,1,j,m,tau), imt, is, ie, ks, ke, scl)
            enddo
            numcfl = numcfl + 1

!           turn off checking on this time step when max is exceeded

            if (numcfl .gt. maxcfl) stabts = .false.
          endif
        enddo
      enddo
      do k=kscfl,kecfl
        dtmax = max(dtuv, dtts*dtxcel(k))
        vmax  = p5*dyu(jrow)/dtmax
        wmax  = p5*dzw(k)/dtmax
        do i=iscfl,iecfl
          umax  = p5*csu(jrow)*dxu(i)/dtmax
          if (abs(100.0*u(i,k,j,1,tau)/umax) .gt. cflup) then
            cflup = abs(100.0*u(i,k,j,1,tau)/umax)
            cflum = u(i,k,j,1,tau)
            icflu = i
            jcflu = jrow
            kcflu = k
          endif
          if (abs(100.0*u(i,k,j,2,tau)/vmax) .gt. cflvp) then
            cflvp = abs(100.0*u(i,k,j,2,tau)/vmax)
            cflvm = u(i,k,j,2,tau)
            icflv = i
            jcflv = jrow
            kcflv = k
          endif
          if (abs(100.0*umask(i,k,j)*adv_vbu(i,k,j)/wmax) .gt. cflwup)
     &      then
            cflwup = abs(100.0*adv_vbu(i,k,j)/wmax)
            cflwum = adv_vbu(i,k,j)
            icflwu = i
            jcflwu = jrow
            kcflwu = k
          endif
          if (abs(100.0*tmask(i,k,j)*adv_vbt(i,k,j)/wmax) .gt. cflwtp)
     &      then
            cflwtp = abs(100.0*adv_vbt(i,k,j)/wmax)
            cflwtm = adv_vbt(i,k,j)
            icflwt = i
            jcflwt = jrow
            kcflwt = k
          endif
       enddo
      enddo

!-----------------------------------------------------------------------
!     look for max peclet numbers using velocities at "tau".
!     look for max reynolds numbers using velocities at "tau".
!-----------------------------------------------------------------------

      do k=kscfl,kecfl
        do i=iscfl,iecfl
          kp1 = min(k+1,km)
          dtdz = (t(i,k,j,1,taum1)-t(i,kp1,j,1,taum1))*dzwr(k)
          ramb = c1/(visc_cbu(i,k,j) + epsln)
# if defined O_isopycmix
          rahb = dtdz/(diff_fbiso(i,k,j) + diff_fb(i,k,j) + epsln)
# else
          rahb = dtdz/(diff_cbt(i,k,j)*dtdz + epsln)
# endif

# if defined O_consthmix
#  if defined O_anisotropic_viscosity
          ramn = c1/(visc_cnu(i,k,jrow) + epsln)
          rame = c1/(visc_ceu(i,k,j) + epsln)
#  else
          rame = c1/(visc_ceu + epsln)
          ramn = c1/(visc_cnu + epsln)
#  endif
#  if defined O_isopycmix
          rahe = c1/(ah+ahisop*fisop(i,jrow,k))
          rahn = c1/(ah+ahisop*fisop(i,jrow,k))
#  else
#   if defined O_bryan_lewis_horizontal
          rahe = c1/(diff_cet(k) + epsln)
          rahn = c1/(diff_cnt(k) + epsln)
#   else
          rahe = c1/(diff_cet + epsln)
          rahn = c1/(diff_cnt + epsln)
#   endif
#  endif
# else
          rame = c1/(visc_ceu(i,k,j) + epsln)
          ramn = c1/(visc_cnu(i,k,j) + epsln)
          rahe = c1/(diff_cet(i,k,j) + epsln)
          rahn = c1/(diff_cnt(i,k,j) + epsln)
# endif
          reyx = abs(u(i,k,j,1,tau)*dxu(i))*rame
          if (reyx .gt. reynx) then
            ireynx = i
            jreynx = jrow
            kreynx = k
            reynx  = reyx
            reynu  = u(i,k,j,1,tau)
            reynmu = c1/rame
          endif
          reyy = abs(u(i,k,j,2,tau)*dyu(jrow))*ramn
          if (reyy .gt. reyny) then
            ireyny = i
            jreyny = jrow
            kreyny = k
            reyny  = reyy
            reynv  = u(i,k,j,2,tau)
            reynmv = c1/ramn
          endif
          kk = min(k+1,km)
          if (k .ge. kmu(i,jrow)) then
            reyz = c0
          else
            reyz =umask(i,kk,j)*abs(adv_vbu(i,k,j)*dzw(k))*ramb
          endif
          if (reyz .gt. reynz) then
            ireynz = i
            jreynz = jrow
            kreynz = k
            reynz  = reyz
            reynw  = adv_vbu(i,k,j)
            reynmw = c1/ramb
          endif
          pecx = abs(u(i,k,j,1,tau)*dxu(i))*rahe
          if (pecx .gt. peclx) then
            ipeclx = i
            jpeclx = jrow
            kpeclx = k
            peclx  = pecx
            peclu  = u(i,k,j,1,tau)
            peclmu = c1/rahe
          endif
          pecy = abs(u(i,k,j,2,tau)*dyu(jrow))*rahn
          if (pecy .gt. pecly) then
            ipecly = i
            jpecly = jrow
            kpecly = k
            pecly  = pecy
            peclv  = u(i,k,j,2,tau)
            peclmv = c1/rahn
          endif
          kk = min(k+1,km)
          if (k .ge. kmt(i,jrow)) then
            pecz = 0.0
          else
            pecz =tmask(i,kk,j)*abs(adv_vbt(i,k,j)*dzw(k))*rahb
          endif
          if (pecz .gt. peclz) then
            ipeclz = i
            jpeclz = jrow
            kpeclz = k
            peclz  = pecz
            peclw  = adv_vbt(i,k,j)
            peclmw = c1/rahb
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     look for ficticious creation of local extremum for tracers
!     by finding local min and max tracer at "tau" and comparing to
!     tracer at "tau+1"
!-----------------------------------------------------------------------

      call getunit (iostab, 'iostab'
     &,             'formatted sequential append')
      ks = max(2,kscfl)
      ke = min(km-1,kecfl)
      is = max(2,iscfl)
      ie = min(imt-1,iecfl)
      do n=1,nt
        do k=ks,ke
          do i=is,ie
            if (tmask(i,k,j) .ne. c0) then
              tbig = max(t(i,k,j,n,tau),t(i,k,j,n,taum1))
              tsml = min(t(i,k,j,n,tau),t(i,k,j,n,taum1))
              do jj=j-1,j+1,2
                if (tmask(i,k,jj) .ne. c0) then
                  tbig = max(tbig,t(i,k,jj,n,tau),t(i,k,jj,n,taum1))
                  tsml = min(tsml,t(i,k,jj,n,tau),t(i,k,jj,n,taum1))
                endif
              enddo

              do ii=i-1,i+1,2
                if (tmask(ii,k,j) .ne. c0) then
                  tbig = max(tbig,t(ii,k,j,n,tau),t(ii,k,j,n,taum1))
                  tsml = min(tsml,t(ii,k,j,n,tau),t(ii,k,j,n,taum1))
                endif
              enddo

              do kk=k-1,k+1,2
                if (tmask(i,kk,j) .ne. c0) then
                  tbig = max(tbig,t(i,kk,j,n,tau),t(i,kk,j,n,taum1))
                  tsml = min(tsml,t(i,kk,j,n,tau),t(i,kk,j,n,taum1))
                endif
              enddo

              tcrit = tdig*abs(t(i,k,j,n,taup1))
              if (tmask(i,k,j) .ne. c0 .and.
     &          ((t(i,k,j,n,taup1) .gt. tbig + tcrit)
     &                .or. (t(i,k,j,n,taup1) .lt. tsml - tcrit))) then
                write (iostab,'(i4, i4, i4, i2, 3g14.7)')
     &           i, k, jrow, n, t(i,k,j,n,taup1), tsml, tbig
              endif
            endif
          enddo
        enddo
      enddo
      call relunit (iostab)

!-----------------------------------------------------------------------
!     look for temperatures and salinities outside ranges used for
!     specifying density coefficients
!-----------------------------------------------------------------------

      call getunit (iobadt, 'iobadt'
     &,             'formatted sequential append')
      call getunit (iobads, 'iobads'
     &,             'formatted sequential append')
      ks = max(1,kscfl)
      ke = min(km,kecfl)
      is = max(2,iscfl)
      ie = min(imt-1,iecfl)
      do k=1,km
        do i=is,ie
          if (tmask(i,k,j) .ne. c0) then
            tlocal = t(i,k,j,1,taup1)
            slocal = t(i,k,j,2,taup1)*1000.0 + 35.0
            n = 1
            if (tlocal .lt. tmink(k) .or. tlocal .gt. tmaxk(k)) then
                write (iobadt,'(i4, i4, i4, i2, 3g14.7)')
     &           i, k, jrow, n, tlocal, tmink(k), tmaxk(k)
            endif
            n = 2
            if (slocal .lt. smink(k) .or. slocal .gt. smaxk(k)) then
                write (iobads,'(i4, i4, i4, i2, 3g14.7)')
     &           i, k, jrow, n, slocal, smink(k), smaxk(k)
            endif
          endif
        enddo
      enddo
      call relunit (iobadt)
      call relunit (iobads)

9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
#endif
      return
      end
