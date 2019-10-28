!=======================================================================
      subroutine sed (is, ie, js, je)

#if defined O_sed
!-----------------------------------------------------------------------
!     sediment model
!-----------------------------------------------------------------------

      implicit none

      integer ie, is, je, js, i, ip, i_year, j, k, n_control, n_debug
      integer n_year, nsed

      real rtsed, r13or, r14or, r13ca, r14ca, rainr, tmp

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "tmngr.h"
#include "switch.h"
#include "sed.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"

      integer n_buried_depths(ipmax)

!     results from sediment
      real ttrtc(ipmax), ttral(ipmax), difal(ipmax), diftc(ipmax)

!     internal variables to sediment
      real temp_p(ipmax), sal_p(ipmax), alk_p(ipmax), tco2_p(ipmax)
      real o2_bw_p(ipmax), rain_clay_p(ipmax), co2_p(ipmax)
      real hco3_p(ipmax), resp_c(ipmax,nzmax,3), cal_c(ipmax,nzmax)

      real dmsk(imt,jmt), tmpij(imt,jmt), tmpip(ipmax)

!-----------------------------------------------------------------------
!     n_control = 0  initialize without setting calcite concentration
!     n_control = 1  find steady state calcite concentration
!     n_control = 2  step ahead time-dependent calcite concentration
!     n_control = 3  maintain constant calcite concentration
!     n_control = 4  return all calcite rain to the overlying water
!     n_control = 5  use old reaction rates to save time but update bury
!-----------------------------------------------------------------------
      n_control = 2

!-----------------------------------------------------------------------
!     you must uncomment lines in sed.f and sediment.f when debugging.
!     Look for: "!!     uncomment for debugging"
!     n_debug = 0  if things are working well
!     n_debug = 1  watch the convergence in stdout
!     n_debug = 2  when convergence is slow, see one a site converging
!     n_debug = 3  when code is blowing up, see everything for all sites
!     n_debug > 3  specific site index
!-----------------------------------------------------------------------
      n_debug = 0

!-----------------------------------------------------------------------
!     transfer 2d ocean bottom variables to 1d arrays
!     convert concentrations from mol m-3 (or umol cm-3) to mol l-1
!-----------------------------------------------------------------------
      rtsed = 1.
!     set the ratio between calcite and organic rain rates
      if (atsed .gt. 0) rtsed = 1./atsed
      do ip=1,ipsed
        i = imap(ip)
        j = jmap(ip)
        temp_p(ip) = sbc(i,j,ibtemp)*rtsed
!       convert from mom units of salinity to psu
        sal_p(ip) = sbc(i,j,ibsalt)*rtsed*1000. + 35.
        alk_p(ip) = sbc(i,j,ibalk)*rtsed*0.001
        tco2_p(ip) = sbc(i,j,ibdic)*rtsed*0.001
        o2_bw_p(ip) = max(sbc(i,j,ibo2)*rtsed*0.001, 1.e-6)
!       convert from umol/cm2/s to mol/cm2/dtsed
        rain_cal_p(ip) = sbc(i,j,ircal)*rtsed*1.e-6*dtsed
        rain_org_p(ip) = sbc(i,j,irorg)*rtsed*1.e-6*dtsed
# if defined O_sed_constrain_rainr
!       set limits on rain ratio by limiting effective organic rain
        rain_org_p(ip) = min(rain_org_p(ip),rain_cal_p(ip)*10.0)
        rain_org_p(ip) = max(rain_org_p(ip),rain_cal_p(ip)*0.1)
!       bring the rain ratio closer to 1 by scaling effective organic rain
        if (rain_org_p(ip) .lt. rain_cal_p(ip)) then
          rain_org_p(ip) = (rain_org_p(ip) - rain_cal_p(ip))*.5
     &                   + rain_cal_p(ip)
        else
          rain_org_p(ip) = 1./((1./rain_org_p(ip)-1./rain_cal_p(ip))*.5
     &                   + 1./rain_cal_p(ip))
        endif
# endif
!       parametrize clay rain from calcite
        rain_clay_p(ip) = 9.1e-6*(rain_org_p(ip)*1.e6)**1.41
# if defined O_npzd_alk
        sbc(i,j,ibalkfx) = 0.
# endif
# if defined O_carbon
        sbc(i,j,ibdicfx) = 0.
# endif
      enddo
# if defined O_time_averages

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_sed_tavg (is, ie, js, je, 0)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_sed_tsi (is, ie, js, je, 0)
# endif

!!     uncomment for debugging
!      if (n_debug .gt. 3) then
!        print*, "input at: ", n_debug
!        print*, "temp_p: ", temp_p(n_debug)
!        print*, "sal_p: ", sal_p(n_debug)
!        print*, "alk_p: ", alk_p(n_debug)
!        print*, "tco2_p: ", tco2_p(n_debug)
!        print*, "o2_bw_p: ", o2_bw_p(n_debug)
!        print*, "rain_cal_p: ", rain_cal_p(n_debug)
!        print*, "rain_org_p: ", rain_org_p(n_debug)
!        print*, "rain_clay_p: ", rain_clay_p(n_debug)
!      endif

!-----------------------------------------------------------------------
!     update sediment year for sediment profiles
!-----------------------------------------------------------------------
# if defined O_sed_profile
      sed_year = sed_year + dtsed
# else
      sed_year = 0.
# endif
      n_year = nint(sed_year)

!-----------------------------------------------------------------------
!     loop through sediment model nsedacc times (>1 during spinup)
!-----------------------------------------------------------------------
      do nsed=1,nsedacc

        call calc_k (temp_p, sal_p, water_z_p, k1, k2, k3, csat, ipsed)

        call calc_buff (alk_p, tco2_p, sal_p, k1, k2, k3, co2_p, hco3_p
     &,                 co3_p, ipsed)

        call setup_pw (co2_p, hco3_p, co3_p, o2_bw_p, carb, o2, csat
     &,                kmax, ipsed, ipmax, nzmax)

        call estimate_rc (rain_org_p, rc, ipsed)

        if (n_control .le. 0) then

!-----------------------------------------------------------------------
!         then initialize without setting calcite concentration
!-----------------------------------------------------------------------
          do ip=1,ipsed
            buried_mass(1,ip) = 500.
            buried_calfrac(1,ip) = calgg(kmax,ip)
            do k=2,ibmax
               buried_mass(k,ip) = 0.
               buried_calfrac(k,ip) = calgg(kmax,ip)
            enddo
          enddo
          call set_pore (calgg, zsed, pore, kmax, ipsed, ipmax, nzmax)
          call get_sed_ml_mass (delz, pore, sed_ml_mass, nzmax, ipsed
     &,                         ipmax, kmax)

        elseif (n_control .eq. 1) then

!-----------------------------------------------------------------------
!         find steady state calcite concentration
!-----------------------------------------------------------------------
!!         uncomment for debugging
!          if (n_debug .ge. 1) write(6,*) 'Starting calss'

          call sed_ss (zsed, delz, form, pore, kmax, o2, zrct, carb
     &,                orgml, orggg, calml, calgg, rain_cal_p
     &,                rain_org_p, rain_clay_p, r13or, r14or, r13ca
     &,                r14ca, rc, dissc, dissn, csat, k1, k2, dopls
     &,                domin, dcpls, dcmin, dbpls, dbmin, resp_c, cal_c
     &,                ttrorg, ttrcal, ttrtc, ttral, diftc, difal
     &,                c_advect, buried_mass, buried_calfrac
     &,                n_buried_depths, ipsed, ipmax, nzmax, ibmax
     &,                n_control, n_debug)

!         reinitialize without setting calcite concentration
          do ip=1,ipsed
            buried_mass(1,ip) = 500.
            buried_calfrac(1,ip) = calgg(kmax,ip)
            do k=2,ibmax
               buried_mass(k,ip) = 0.
               buried_calfrac(k,ip) = calgg(kmax,ip)
            enddo
          enddo
          call set_pore (calgg, zsed, pore, kmax, ipsed, ipmax, nzmax)
          call get_sed_ml_mass (delz, pore, sed_ml_mass, nzmax, ipsed
     &,                         ipmax, kmax)

        elseif (n_control .eq. 2) then

!-----------------------------------------------------------------------
!         step ahead time-dependent calcite concentration
!-----------------------------------------------------------------------
          call sed_const_cal (zsed, delz, form, pore, kmax, o2, zrct
     &,                       carb, orgml, orggg, calml, calgg
     &,                       rain_cal_p, rain_org_p, rain_clay_p, r13or
     &,                       r14or, r13ca, r14ca, rc, dissc, dissn
     &,                       csat, k1, k2, dopls, domin, dcpls, dcmin
     &,                       dbpls, dbmin, resp_c, cal_c, ttrorg
     &,                       ttrcal, ttrtc, ttral, diftc, difal
     &,                       c_advect, ipsed, ipmax, nzmax, n_control
     &,                       n_debug)

          call bury (n_control, zsed, delz, pore, kmax, calgg, orggg
     &,              sed_ml_mass, rain_cal_p, ttrcal, c_advect
     &,              rain_clay_p, n_year, buried_mass, buried_calfrac
     &,              depth_age, ipsed, ipmax, nzmax, ibmax)

        elseif (n_control .eq. 3) then

!-----------------------------------------------------------------------
!        maintain constant calcite concentration (update diss but don't bury)
!-----------------------------------------------------------------------
         call sed_const_cal (zsed, delz, form, pore, kmax, o2, zrct
     &,                      carb, orgml, orggg, calml, calgg
     &,                      rain_cal_p, rain_org_p, rain_clay_p, r13or
     &,                      r14or, r13ca, r14ca, rc, dissc, dissn
     &,                      csat, k1, k2, dopls, domin, dcpls, dcmin
     &,                      dbpls, dbmin, resp_c, cal_c, ttrorg
     &,                      ttrcal, ttrtc, ttral, diftc, difal
     &,                      c_advect, ipsed, ipmax, nzmax, n_control
     &,                      n_debug)

        elseif (n_control .eq. 4) then

!-----------------------------------------------------------------------
!         return all calcite rain to the overlying water
!-----------------------------------------------------------------------
          do ip=1,ipsed
            ttrcal(ip) = rain_cal_p(ip)
          enddo

        elseif (n_control .eq. 5) then

!-----------------------------------------------------------------------
!         use old reaction rates but update bury
!-----------------------------------------------------------------------
          call bury (n_control, zsed, delz, pore, kmax, calgg, orggg
     &,              sed_ml_mass, rain_cal_p, ttrcal, c_advect
     &,              rain_clay_p, n_year, buried_mass, buried_calfrac
     &,              depth_age, ipsed, ipmax, nzmax, ibmax)

        endif

      enddo

!-----------------------------------------------------------------------
!     set variables for coupling back to the ocean model
!-----------------------------------------------------------------------
# if defined O_sed_weath_diag
      weathflx = 0.
# endif
# if defined O_npzd_alk
      sbc(:,:,ibalkfx) = 0.
# endif
# if defined O_carbon
      sbc(:,:,ibdicfx) = 0.
# endif
      if (atsed .gt. 0) then
        rtsed = 1.e6/dtsed
        do ip=1,ipsed
          i = imap(ip)
          j = jmap(ip)
          tmp =  (ttrcal(ip) - rain_cal_p(ip))*rtsed
# if defined O_sed_weath_diag
!         diagnose weathering flux as he total burial flux (umol/s)
          weathflx = weathflx - tmp*dxt(i)*dyt(j)*cst(j)
# endif
!         set source terms, convert to umol cm-3 s-1
!         the source terms are really a correction to assumed zero burial
!         during the previous coupling time and is applied over the next
          tmp = tmp*dztr(kmt(i,j))
# if defined O_npzd_alk
          sbc(i,j,ibalkfx) = 2.*tmp
# endif
# if defined O_carbon
          sbc(i,j,ibdicfx) = tmp
# endif
!         zero accumulators
          sbc(i,j,ibtemp) = 0.
          sbc(i,j,ibsalt) = 0.
# if defined O_npzd_alk
          sbc(i,j,ibalk) = 0.
# endif
# if defined O_carbon
          sbc(i,j,ibdic) = 0.
# endif
# if defined O_npzd_o2
          sbc(i,j,ibo2) = 0.
# endif
          sbc(i,j,ircal) = 0.
          sbc(i,j,irorg) = 0.
        enddo
        atsed = 0.
      endif
#endif

      return
      end
