      subroutine glsbc (is, ie, js, je)

#if defined O_mtlm && defined O_embm
!-----------------------------------------------------------------------
!     Get the boundary conditions for the land model.
!-----------------------------------------------------------------------

       implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cembm.h"
      include "csbc.h"
      include "calendar.h"
      include "tmngr.h"
      include "switch.h"
      include "levind.h"
      include "mtlm.h"
# if defined O_mtlm_carbon_13
      include "mtlmc13.h"
# endif
# if defined O_mtlm_carbon_14
      include "mtlmc14.h"
# endif
      include "insolation.h"
      include "atm.h"
# if defined O_crop_data_transient || defined O_pasture_data_transient || defined O_agric_data_transient
      include "veg.h"
# endif

      integer ie, is, je, js
      integer i, j, L, n

      real calday, fc, fd, fe, ff, cosz(POINTS), dt, t, pi, degrad, C2K

      fc = 1.0/atatm
      fd = 1.e-2/atatm
      fe = 1.e-3/atatm
      ff = 10./atatm
      atlnd = 0.
      pi = 4.*atan(1.)
      degrad = pi/180.
      C2K = 273.15

!----------------------------------------------------------------------
!     Calculate the diurnal cycle in the SW radiation
!----------------------------------------------------------------------
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      dt = SEC_DAY/STEP_DAY
      do n=1,STEP_DAY
        t = real(n-1)*dt
        call zenith (POINTS, t, dt, SEC_DAY, LAT, LONG, sindec, cosz)
        SUN(:,n) = solarconst*eccf*cosz(:)
      enddo
      do j=2,jmt-1
        do i=2,imt-1
          L = land_map(i,j)
          if (L .ne. 0) then
            t = 0.
            do n=1,STEP_DAY
              t = t + SUN(L,n)
            enddo
!           make sure the daily insolations agree
            SUN(L,:) = SUN(L,:)*solins(i,j)*STEP_DAY/(t + epsln)
!           calculate downward shortwave at the surface
# if defined O_sulphate_data || defined O_sulphate_data_transient
            SUN(L,:) = SUN(L,:)*1.e-3*sbc(i,j,iaca)*pass
     &                 *(sbc(i,j,isca) - sulph(i,j,2))
# else
            SUN(L,:) = SUN(L,:)*1.e-3*sbc(i,j,iaca)*pass*sbc(i,j,isca)
# endif
          endif
        enddo
      enddo

!----------------------------------------------------------------------
!     Calculate the time of maximum temperature. Assume at local noon.
!----------------------------------------------------------------------
      do L=1,POINTS
        TIME_MAX(L) = SEC_DAY*(0.5 + LONG(L)/360.)
        if (TIME_MAX(L) .lt. 0) then
          TIME_MAX(L) = TIME_MAX(L)+SEC_DAY*int(1.+TIME_MAX(L)/SEC_DAY)
        elseif (TIME_MAX(L) .gt. SEC_DAY) then
          TIME_MAX(L) = TIME_MAX(L)-SEC_DAY*int(TIME_MAX(L)/SEC_DAY)
        endif
      enddo

!----------------------------------------------------------------------
!     Set other boundary conditions and zero accumulators
!----------------------------------------------------------------------
      do j=2,jmtm1
        do i=2,imtm1
          L = land_map(i,j)
          if (L .ne. 0) then
!           set boundary conditions for land
            RAIN(L) = ff*sbc(i,j,ipr)
            SNOW(L) = ff*sbc(i,j,ips)
            SW_C(L) = fe*sbc(i,j,iswr)
            T_C(L) = fc*sbc(i,j,iat) + C2K
            WIND(L) = fd*sbc(i,j,iaws)
            RH_C(L) = fc*sbc(i,j,irh)
            DTEMP_DAY(L) = sbc(i,j,idtr)
            TSTAR(L,:) = T_C(L)
            TSOIL(L) = T_C(L)
# if defined O_carbon_co2_2d
            CO2(L) = 1.0E-6*at(i,j,2,ico2)*EPCO2
# else
            CO2(L) = 1.0E-6*co2ccn*EPCO2
# endif
# if defined O_crop_data_transient || defined O_pasture_data_transient || defined O_agric_data_transient
            FRACA(L) = agric(i,j,2)
# endif
!           zero atmosphere accumulation variables
            sbc(i,j,iro) = 0.
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
            sbc(i,j,isca) = 0.
# if defined O_carbon
            sbc(i,j,inpp) = 0.
            sbc(i,j,isr) = 0.
# endif
# if defined O_mtlm_carbon_13
            sbc(i,j,inpp13) = 0.
            sbc(i,j,isr13) = 0.
            sbc(i,j,iburn13) = 0.
# endif
# if defined O_mtlm_carbon_14
            sbc(i,j,inpp14) = 0.
            sbc(i,j,isr14) = 0.
            sbc(i,j,iburn14) = 0.
# endif
          endif
        enddo
      enddo
# if defined O_mtlm_carbon_13
      rc13a = (dc13ccn/1000. + 1.)*rc13std
c      print*,'dc13ccn,rc13a',dc13ccn,rc13a
# endif
# if defined O_mtlm_carbon_14
c      rc14a = c14ccn/(co2ccn-c14ccn)
      rc14a = (dc14ccn/1000. + 1.)*rc14std
c      print*,'dc14ccn,rc14a',dc14ccn,rc14a
# endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_mtlm_tavg (is, ie, js, je, 0)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_mtlm_tsi (is, ie, js, je, 0)
# endif
#endif

      return
      end
