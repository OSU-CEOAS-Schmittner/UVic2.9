      subroutine therm  (is, ie, js, je)

#if defined O_ice && defined O_embm
!=======================================================================
!     thermodynamic ice model

!     Note: if run with embm this routine must be called after "fluxes"
!           and before "solve"

!     calculates ice and open water growth rates based on the surface
!     energy budget. see Parkinson and Washington, JGR, Vol.84, C1,
!     311-337, 1979 and Hibler, JPO, Vol.9, 815-846, 1979

!     heat and fresh water fluxes between the ocean and atmosphere
!     are adjusted depending on ice growth or melt. ice thickness is
!     changed by the amount of growth or melt and sublimation
!=======================================================================

      implicit none

      integer i, ie, iem1, imax, index, iter, is, isp1, j, je, jem1
      integer jmax, js, jsp1, maxit

      real ai, aice3, al, amin, ao, as, as_agric, ca,  delta, df, dh
      real dha, dhflxi, dhflxs, dhi, hice3, dho, dhs, hsno3, dhss
      real dhstot, dhtot, dswr, dt, dultnt, dulwr, dusens, emax, f, fa
      real fas, fb, fbot, fcond, fd, fds, fe, ff, ffs, fh, fl, fls, fm
      real fn, fptf, fpts, ftopi, ftopo, ho, hsextra, qair,qice, sla
      real sub, tair, tcdh, ti, tiold, tol, tol2, ultnt, ulwr, usens
      real wspd, zintfc, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "coord.h"
      include "grdvar.h"
      include "veg.h"
# if defined O_mtlm
      include "mtlm.h"
# endif

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      fa = dts/(rhoice*flice)
      fb = 0.94*rhoatm*cpatm
      fd = rhoatm/rhoice
      fe = rhoatm*slice
      ff = rhoice*flice
      fh = 21.8746*265.5
      fas = dts/(rhosno*flice)
      fds = rhoatm/rhosno
      ffs = rhosno*flice
      ho = 1.0
      amin = 0.15
      maxit = 10
      tol = 0.01
      emax = 0.0
      imax = 0
      jmax = 0
      fptf = 0.0
# if defined O_ice_evp
      index = 1
# else
      index = lf
# endif
      sla = zw(1)*secday/dampice/2.389e-8
      tol2 = tol*2.0
      C2K = 273.15

      do j=jsp1,jem1
        do i=isp1,iem1

          hsno3 = hsno(i,j,index)
          hice3 = hice(i,j,index)
          aice3 = 0

!-----------------------------------------------------------------------
!         set the incoming shortwave over snow and ice
!-----------------------------------------------------------------------
!         if snow is less than 25 cm linearly reduce to ice albedo
          as = min(hsno(i,j,2)*0.04/(aice(i,j,2) + epsln), c1)
          ca = ice_calb*(c1 - as) + sno_calb*as
# if defined O_sulphate_data || defined O_sulphate_data_transient
          dswr = solins(i,j)*sbc(i,j,iaca)*pass
     &           *max(c0, ca - sulph(i,j,2))
# else
          dswr = solins(i,j)*sbc(i,j,iaca)*pass*ca
# endif
          wspd = sbc(i,j,iws)

# if defined O_ice_evp
          aice(i,j,index) = min(c1, aice(i,j,index))
# endif
          ai = aice(i,j,2)

# if  defined O_mtlm
          if (tmsk(i,j) .lt. 0.5 .and. land_map(i,j) .eq. 0) then
# else
          if (tmsk(i,j) .lt. 0.5) then
# endif

!-----------------------------------------------------------------------
!           land points
!-----------------------------------------------------------------------
            tair = at(i,j,2,isat) - elev(i,j)*rlapse
# if defined O_landice_data
     &           - hicel(i,j,2)*rlapse
# endif
# if defined O_sealev || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
# endif
            aice3 = 0.0
!           snow masking distance for different vegetation types
            as = hsno(i,j,2)/(100.0*veg_smd(iveg(i,j)))
# if defined O_crop_data || defined O_pasture_data || defined O_agric_data
            as_agric = hsno(i,j,2)/(100.0*veg_smd(iagric))
            as_agric = min(c1, max(c0, as_agric))
            as = (1.-agric(i,j,2))*as + agric(i,j,2)*as_agric
# endif
!           limit snow coverage between 0 and 1
            aice3 = min(c1, max(c0, as))
# if defined O_landice_data
            aice3 = max(aice3, aicel(i,j,2))
# endif

!-----------------------------------------------------------------------
!           start loop for grid points with snow or ice
!-----------------------------------------------------------------------
            if (ai .gt. c0) then

!-----------------------------------------------------------------------
!             find snow temperature by balancing the surface heat budget
!               dwsr = ultnt + usens + ulwr
!             using Newton's method:
!               t(i+1) = t(i) - f(t(i))/df(t(i))
!             where:
!               f(t(i)) = dwsr - ultnt - usens - ulwr
!               -df(t(i)) = dultnt - dusens - dulwr
!-----------------------------------------------------------------------
              al = 1.0 - ai
              ti = tice(i,j)
              tiold = tice(i,j)
              fm = esatm*(tair + C2K)**4
              fls = fe*dalt_i*wspd
              dusens = fb*dalt_i*wspd
              qair = rh(i,j)*cssh*exp(17.67*tair/(tair + 243.5))
              iter = 0
              delta = tol2
              do while (abs(delta) .gt. tol .and. iter .le. maxit)
                iter = iter + 1
                dt = ti - tair
                qice = cssh*exp(21.8746*ti/(ti + 265.5))
                if (qice .gt. qair) then
                  ultnt = fls*(qice - qair)
                  dultnt = fls*qice*21.8746*265.5/(ti + 265.5)**2
                else
                  ultnt = 0.0
                  dultnt = 0.0
                endif
                usens = dusens*dt
                ulwr = esice*(ti + C2K)**4 - fm
                dulwr = 4.0*esice*(ti + C2K)**3
                f = dswr - ultnt - usens - ulwr
                df = dultnt + dusens + dulwr
                delta = f/df
                ti = ti + delta
              enddo
              if (iter .gt. maxit) then
!               if not converged, set to last converged temperature
                if (abs(delta) .gt. emax .and. ti .lt. fptf) then
                  emax = abs(delta)
                  imax = i
                  jmax = j
                  ti = tiold
                endif
              endif

!-----------------------------------------------------------------------
!             set maximum tice to freezing and calculate fluxes
!-----------------------------------------------------------------------
              ti = min(ti, fptf)
              dt = ti - tair
              qice = cssh*exp(21.8746*ti/(ti + 265.5))
              sub = max(c0, fds*dalt_i*wspd*(qice - qair))
              usens = dusens*dt
              ulwr = esice*(ti + C2K)**4 - fm

!             ensure that snow sublimated does not exceed hsno
              dha = -dts*sub*ai
              if (-dha .gt. hsno3) then
                dha = -hsno3
                sub = -dha/(ai*dts)
              endif
              ultnt = rhosno*slice*sub
              sub = sub*rhosno

              ftopi = dswr - ulwr - usens - ultnt

!-----------------------------------------------------------------------
!             add ice/snow covered area fluxes to land fluxes
!-----------------------------------------------------------------------
              tice(i,j) = ti
              evap(i,j) = evap(i,j)*al + sub*ai
              dnswr(i,j) = dnswr(i,j)*al + dswr*ai
              upltnt(i,j) = upltnt(i,j)*al + ultnt*ai
              upsens(i,j) = upsens(i,j)*al + usens*ai
              uplwr(i,j) = uplwr(i,j)*al + ulwr*ai

!-----------------------------------------------------------------------
!             calculate total change in snow volume on land
!             allocate snowmelt to flux for input to bucket
!-----------------------------------------------------------------------
              dhs = 0.0
              dhtot = 0.0
              if (tice(i,j) .ge. fptf .and. ftopi .gt. 0.0)
     &          dhs = -ai*fas*ftopi
              dhs = min(0.0, max(-(hsno3 + dha), dhs))
              dhtot = dhs + dha
              hsno3 = hsno3 + dhtot

              flux(i,j,ishum) = flux(i,j,ishum) + dhs*rhosno/dts

!             ensure fluxes are balanced since land can't absorb error
              upsens(i,j) = dnswr(i,j) - upltnt(i,j) - uplwr(i,j)
     &                   + dhs*ffs/dts

!-----------------------------------------------------------------------
!     if no snow on land
!-----------------------------------------------------------------------
            else
              tice(i,j) = 0.0
            endif

# if !defined O_ice_cpts
#  if  defined O_mtlm
          elseif (tmsk(i,j) .ge. 0.5) then
#  else
          else
#  endif

!-----------------------------------------------------------------------
!           ocean points
!-----------------------------------------------------------------------
            ao = 1.0 - ai
            tair = at(i,j,2,isat)
#  if defined O_sealev_data || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
#  endif

!-----------------------------------------------------------------------
!           calculate fluxes to and from the ocean (without ice)
!-----------------------------------------------------------------------
            ftopo = dnswr(i,j) - uplwr(i,j) - upsens(i,j) - upltnt(i,j)
            fbot = sla*(frzpt(i,j) - sbc(i,j,isst))
!           calculate growth of ice in open water areas
            dho = fa*(fbot - ftopo)

            if (ai .ne. 0.0) then

!-----------------------------------------------------------------------
!             find ice temperature by balancing the surface heat budget:
!               tcdh*(ti - fpts) = dswr - ultnt - usens - ulwr
!             using Newton's method:
!               t(i+1) = t(i) - f(t(i))/df(t(i))
!             where:
!               f(t(i)) = dswr - ultnt - usens - ulwr - tcdh*(ti - fpts)
!               -df(t(i)) = dultnt + dusens + dulwr + tcdh
!-----------------------------------------------------------------------
              tcdh = condice/(hice(i,j,2) + 6.5*hsno(i,j,2))
              ti = tice(i,j)
              tiold = tice(i,j)
              fpts = frzpt(i,j)
              fm = esatm*(tair + C2K)**4
              fn = 4.0*esice

              fl = fe*dalt_i*wspd
              dusens = fb*dalt_i*wspd
              qair = at(i,j,2,ishum)
              iter = 0
              delta = tol2
              do while (abs(delta) .gt. tol .and. iter .le. maxit)
                iter = iter + 1
                dt = ti - tair
                qice = cssh*exp(21.8746*ti/(ti + 265.5))
                if (qice .gt. qair) then
                  ultnt = fl*(qice - qair)
                  dultnt = fl*qice*fh/(ti + 265.5)**2
                else
                  ultnt = 0.0
                  dultnt = 0.0
                endif
                usens = dusens*dt
                ulwr = esice*(ti + C2K)**4 - fm
                dulwr = fn*(ti + C2K)**3
                f = dswr - ultnt - usens - ulwr - tcdh*(ti - fpts)
                df = dultnt + dusens + dulwr + tcdh
                delta = f/df
                ti = ti + delta
              enddo
              if (iter .gt. maxit) then
!               if not converged, set to last converged temperature
                if (abs(delta) .gt. emax .and. ti .lt. fptf) then
                  emax = abs(delta)
                  imax = i
                  jmax = j
                  ti = tiold
                endif
              endif

!-----------------------------------------------------------------------
!             set maximum tice to freezing and calculate fluxes
!-----------------------------------------------------------------------
              ti = min(ti, fptf)
              dt = ti - tair
              qice = cssh*exp(21.8746*ti/(ti + 265.5))
              sub = max(c0, dalt_i*wspd*(qice - qair))
              ultnt = fe*sub
              fcond = tcdh*(ti - fpts)
              if (hsno(i,j,index) .gt. 0.0) then
                sub = fds*sub
                dha = -dts*sub
                sub = sub*ai*rhosno
              else
                sub = fd*sub
                dha = -dts*sub
                sub = sub*ai*rhoice
              endif

              usens = dusens*dt
              ulwr = esice*(ti + C2K)**4 - fm

!-----------------------------------------------------------------------
!             add ice covered area fluxes to ocean area fluxes
!-----------------------------------------------------------------------
              tice(i,j) = ti
              dnswr(i,j) = dnswr(i,j)*ao + dswr*ai
              upltnt(i,j) = upltnt(i,j)*ao + ultnt*ai
              upsens(i,j) = upsens(i,j)*ao + usens*ai
              uplwr(i,j) = uplwr(i,j)*ao + ulwr*ai
              ftopi = dswr - ulwr - usens - ultnt

!-----------------------------------------------------------------------
!             calculate change in ice volume due to sublimation of
!             ice (dha). adjust evaporation to the atmosphere to
!             account for sublimation from ice. subtract this
!             adjustment from the ocean freshwater flux
!-----------------------------------------------------------------------
              if (addflxa) flux(i,j,ishum) = flux(i,j,ishum) + dts*sub
              evap(i,j) = evap(i,j)*ao + sub
            else
              tice(i,j) = sbc(i,j,isst)
              ftopi = 0.0
              fcond = 0.0
              dha = 0.0
            endif

!-----------------------------------------------------------------------
!           calculate total change in ice volume (dh)
!-----------------------------------------------------------------------
            dha = dha*ai
            dhflxs = 0.0
            dhs = 0.0
            if (hsno(i,j,index) .le. 0.0) then
!             total growth of ice from the ocean
              dhi = ai*fa*(fbot - ftopi) + ao*dho
!             total growth (loss + sublimation limited to total amount)
              dh = max(-hice(i,j,index), dhi + dha)
!             adjust ocean fluxes for ice growth or melt + sublimation
              dhflxi = dh - dha
            else
!             total growth of ice from the ocean
              dhi = ai*fa*(fbot - fcond) + ao*dho
!             loss of snow due to melt
              if (tice(i,j) .ge. fptf) dhs = ai*fas*(fcond - ftopi)
!             total loss of snow including sublimation
              dhs = dhs + dha
!             check if snow loss greater than total
              if (-dhs .gt. hsno(i,j,index)) then
!               take extra melt from ice
                dhi = dhi + rhosno/rhoice*(dhs + hsno(i,j,index))
!               remove all snow
                dhs = -hsno(i,j,index)
              endif
!             adjust ocean fluxes for snow melt and sublimation
              dhflxs = dhs - dha
!             ice loss limited to total
              dh = max(-hice(i,j,index), dhi)
!             adjust ocean fluxes for ice growth or melt
              dhflxi = dh
            endif

!-----------------------------------------------------------------------
!           calculate new area and thickness from thermodynamics
!-----------------------------------------------------------------------
!           use minimum area (amin) of open water
            ai = max(amin, aice(i,j,index))
            aice3 = aice(i,j,index) + ((1.0 - ai)*max(c0, dho)/ho
     &         + 0.5*min(c0, dhi)*ai/(hice(i,j,index) + epsln))
!           update ice and snow thickness
            hice3 = hice(i,j,index) + dh
            hsno3 = hsno(i,j,index) + dhs

!           lower ice area where thickness is < 1 cm
            aice3 = min(aice3, hice3)
!           max ice area where thickness is > 10 m
            aice3 = max(aice3, hice3*0.001)
            aice3 = max(c0, min(c1, aice3))
            if (aice3 .eq. 0.0) then
              dhflxs = dhflxs - hsno3
              hsno3 = 0.0
            endif

!           check if the weight of the snow pushes the ice/snow
!           interface below the waterline (if so, change snow to ice)
            zintfc = hice3 - (rhosno*hsno3 + rhoice*hice3)/rhoocn
            if (zintfc .lt. 0.0) then
              dhss = rhoice/rhosno*zintfc
              if (-dhss .gt. hsno3) then
                write(*,*) '==> Warning: dhss is too large: ',dhss
                dhss = -hsno3
              endif
              hice3 = hice3 - rhosno/rhoice*dhss
              hsno3 = hsno3 + dhss
            endif
            hsno3 = max(hsno3, c0)

!-----------------------------------------------------------------------
!           adjust fluxes to the ocean due to ice melt or growth
!-----------------------------------------------------------------------
            if (addflxa) then
              flux(i,j,isat) = flux(i,j,isat) + ff*dhflxi + ffs*dhflxs
              flux(i,j,ishum) = flux(i,j,ishum) - rhoice*dhflxi
     &                        - rhosno*dhflxs
#  if defined O_convect_brine
              ao = 1.0 - aice(i,j,2)
              dho = ao*fa*(fbot - ftopo)
              if (dho .gt. 0.) then
                flux(i,j,ishum) = flux(i,j,ishum) + rhoice*dho
                cbf(i,j,0) = cbf(i,j,0) - rhoice*dho
                cba(i,j,0) = cba(i,j,0) + ao*dts
              endif
              ai = aice(i,j,2)
              dhi = ai*fa*(fbot - ftopi)
              if (dhi .gt. 0.) then
                flux(i,j,ishum) = flux(i,j,ishum) + rhoice*dhi
                cbf(i,j,1) = cbf(i,j,1) - rhoice*dhi
                cba(i,j,1) = cba(i,j,1) + ai*dts
              endif
#   endif
            endif
#  if defined O_plume_brine
            subflux(i,j,ishum) = subflux(i,j,ishum) - rhoice*dhflxi
#  endif
# endif

          endif

!-----------------------------------------------------------------------
!         shuffle time levels
!-----------------------------------------------------------------------
          hice(i,j,1) = hice(i,j,2)
          hice(i,j,2) = hice3
          aice(i,j,1) = aice(i,j,2)
          aice(i,j,2) = aice3
          hsno(i,j,1) = hsno(i,j,2)
          hsno(i,j,2) = hsno3

        enddo
      enddo

      if (emax .gt. 0.0) write (stdout,*)
     &  '==> Warning: ice temperature not converging: emax, i, j:'
     &, emax, imax, jmax

!-----------------------------------------------------------------------
!     set boundary conditions
!-----------------------------------------------------------------------
      call embmbc (hice(1,1,2))
      call embmbc (aice(1,1,2))
      call embmbc (hsno(1,1,2))

#endif
      return
      end
