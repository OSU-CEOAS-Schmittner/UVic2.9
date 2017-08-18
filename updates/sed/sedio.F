       subroutine sedout (is, ie, js, je)

#if defined O_sed
!=======================================================================
!     output routine for sediment model

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name
      character(32) :: nstamp

      integer ie, is, je, js, ntrec, nyear, nmonth, nday, nhour, nmin
      integer nsec

      real avgper, time, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "grdvar.h"
      include "levind.h"
      include "sed.h"
      include "iounit.h"
      include "switch.h"
      include "tmngr.h"

!-----------------------------------------------------------------------
!     write sediment model diagnostics
!-----------------------------------------------------------------------

# if defined O_time_step_monitor
!     accumulate time step integrals
      call ta_sed_tsi (is, ie, js, je, 1)

      if (tsits .and. ntatis .ne. 0) then

        call ta_sed_tsi (is, ie, js, je, 2)

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

        call def_tsi
        call def_tsi_sed (fname)
#  if defined O_save_carbon_totals
!       convert from mol cm-2 s-1 to Pg year-1
        tai_cfo2s = (tai_rain_cal - tai_ttrcal)*sedsa*12.e-15
     &               *yrlen*daylen
!       convert from kg s-1 to Pg year-1
        tai_cfl2o = tai_weathflx*1.e-12*yrlen*daylen
#  else
        tai_cfo2s = 0.
        tai_cfl2o = 0.
#  endif

#  if !defined O_save_time_relyear0
!       make output time relative to year 1
        time = time - 1.
#  endif
        avgper = tsiper*accel
        if (avgper .le. 1e-6) avgper = 0.
#  if defined O_save_time_endper
        tmp = 0.
#  elif defined O_save_time_startper
        tmp = 1.
#  else
        tmp = 0.5
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
        time = time - tmp*avgper/360.
#   elif defined O_calendar_gregorian
        time = time - tmp*avgper/365.25
#   else
        time = time - tmp*avgper/365.
#   endif
#  else
#   if defined O_calendar_360_day
        time = time*360. - tmp*avgper
#   elif defined O_calendar_gregorian
        time = time*365.25 - tmp*avgper
#   else
        time = time*365. - tmp*avgper
#   endif
#  endif

        call sed_tsi_out (fname, avgper, time, stamp, tai_ttrcal
     &,   tai_rain_cal, tai_cal, tai_calmass, tai_calmass_bur
     &,   tai_co3, tai_co3sat, tai_weathflx, tai_rainr, tai_carblith
     &,   tai_cfo2s, tai_cfl2o, ntrec)

        call ta_sed_tsi (is, ie, js, je, 0)

      endif

# endif
# if defined O_time_averages
!     accumulate time averages
      call ta_sed_tavg (is, ie, js, je, 1)

      if (timavgts .and. ntatss .ne. 0) then

!-----------------------------------------------------------------------
!       write sediment model time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_sed_tavg (is, ie, js, je, 2)

!       write time averaged data

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        call def_tavg
        call def_tavg_sed (fname)

#  if !defined O_save_time_relyear0
!       make output time relative to year 1
        time = time - 1.
#  endif
        avgper = timavgper*accel
        if (avgper .le. 1e-6) avgper = 0.
#  if defined O_save_time_endper
        tmp = 0.
#  elif defined O_save_time_startper
        tmp = 1.
#  else
        tmp = 0.5
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
        time = time - tmp*avgper/360.
#   elif defined O_calendar_gregorian
        time = time - tmp*avgper/365.25
#   else
        time = time - tmp*avgper/365.
#   endif
#  else
#   if defined O_calendar_360_day
        time = time*360. - tmp*avgper
#   elif defined O_calendar_gregorian
        time = time*365.25 - tmp*avgper
#   else
        time = time*365. - tmp*avgper
#   endif
#  endif

        call sed_tavg_out (fname, is, ie, js, je, imt, jmt, xt, yt, xu
     &,   yu, dxt, dyt, dxu, dyu, avgper, time, stamp, ta_ttrcal
     &,   ta_rain_cal, ta_cal, ta_calmass, ta_calmass_bur, ta_co3
     &,   ta_co3sat, ta_rainr, map_sed, kmt, tlat, tlon, tgarea, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Sed time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

!       zero time average accumulators

        call ta_sed_tavg (is, ie, js, je, 0)

      endif

# endif
!-----------------------------------------------------------------------
!       write sediment model restart
!-----------------------------------------------------------------------

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_sed (0, fname)
          call sed_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_sed (1, fname)
          call sed_rest_out (fname, is, ie, js, je)
        endif
      endif

      return
      end

# if defined O_sed && defined O_time_averages
      subroutine ta_sed_tavg (is, ie, js, je, iflag)

!=======================================================================
!     sediment data time averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer i, ie, iflag, ip, is, j, je, js, k

      real rntatss, rtsed

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "sed.h"

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatss = 0
        ta_ttrcal(is:ie,js:je) = 0.
        ta_rain_cal(is:ie,js:je) = 0.
        ta_cal(is:ie,js:je) = 0.
        ta_calmass(is:ie,js:je) = 0.
        ta_calmass_bur(is:ie,js:je) = 0.
        ta_co3(is:ie,js:je) = 0.
        ta_co3sat(is:ie,js:je) = 0.
        ta_rainr(is:ie,js:je) = 0.

      elseif (iflag .eq. 1) then

!       accumulate
        ntatss = ntatss + 1
        rtsed = 1./dtsed
        do j=js,je
          do i=is,ie
            ip = map_sed(i,j)
            if (ip .gt. 0) then
!             convert fluxes from per dtsed to per second with rtsed
              ta_ttrcal(i,j) = ta_ttrcal(i,j) + ttrcal(ip)*rtsed
              ta_rain_cal(i,j) = ta_rain_cal(i,j)
     &                         + rain_cal_p(ip)*rtsed
              ta_cal(i,j) = ta_cal(i,j) + calgg(kmax,ip)
              ta_calmass(i,j) = ta_calmass(i,j)
     &                        + calgg(kmax,ip)*sed_ml_mass(ip)
              do k=1,ibmax
                ta_calmass_bur(i,j) = ta_calmass_bur(i,j)
     &                              + buried_mass(k,ip)
     &                               *buried_calfrac(k,ip)
              enddo
              ta_co3(i,j) = ta_co3(i,j) + co3_p(ip)
              ta_co3sat(i,j) = ta_co3sat(i,j) + csat(ip)
              ta_rainr(i,j) = ta_rainr(i,j) + rain_cal_p(ip)
     &                        /(rain_org_p(ip) + 1.e-20)
            endif
          enddo
        enddo

      elseif (iflag .eq. 2 .and. ntatss .ne. 0) then

!       average
        rntatss = 1./float(ntatss)
        ta_ttrcal(is:ie,js:je) = ta_ttrcal(is:ie,js:je)*rntatss
        ta_rain_cal(is:ie,js:je) = ta_rain_cal(is:ie,js:je)*rntatss
        ta_cal(is:ie,js:je) = ta_cal(is:ie,js:je)*rntatss
!       convert from g CaCO3 cm-2 to kg C m-2 (1e4*1e-3*12/100)
        ta_calmass(is:ie,js:je) = ta_calmass(is:ie,js:je)*rntatss*1.2
!       convert from g CaCO3 cm-2 to kg C m-2 (1e4*1e-3*12/100)
        ta_calmass_bur(is:ie,js:je) = ta_calmass_bur(is:ie,js:je)
     &                                *rntatss*1.2
        ta_co3(is:ie,js:je) = ta_co3(is:ie,js:je)*rntatss
        ta_co3sat(is:ie,js:je) = ta_co3sat(is:ie,js:je)*rntatss
        ta_rainr(is:ie,js:je) = ta_rainr(is:ie,js:je)*rntatss
      endif

      return
      end
# endif

# if defined O_time_step_monitor
      subroutine ta_sed_tsi (is, ie, js, je, iflag)

!=======================================================================
!     sediment data time integral averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer ie, iflag, ip, is, je, js, k

      real rntatis, tmp, rtsed

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "grdvar.h"
      include "sed.h"

      real dmsk(imt,jmt), tmpij(imt,jmt), tmpip(ipmax)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatis = 0
        tai_ttrcal = 0.
        tai_rain_cal = 0.
        tai_cal = 0.
        tai_calmass = 0.
        tai_calmass_bur = 0.
        tai_co3 = 0.
        tai_co3sat = 0.
        tai_weathflx = 0.
        tai_rainr = 0.
        tai_carblith = 0.

      elseif (iflag .eq. 1) then

!       accumulate
        ntatis = ntatis + 1
        dmsk(:,:) = 0.
        tmpip(:) = 0.
        where (map_sed(is:ie,js:je) .gt. 0.) dmsk(:,:) = 1.
        rtsed = 1./dtsed
!       convert fluxes from per dtsed to per second with rtsed
        tmpip(1:ipmax) = ttrcal(1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_ttrcal = tai_ttrcal + tmp*rtsed
        tmpip(1:ipmax) = rain_cal_p(1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_rain_cal = tai_rain_cal + tmp*rtsed
!       proportion to percent
        tmpip(1:ipmax) = calgg(kmax,1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_cal = tai_cal + tmp
        do ip=1,ipsed
          tmpip(ip) = calgg(kmax,ip)*sed_ml_mass(ip)
        enddo
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areatot (tmpij, dmsk, tmp)
!       convert from g CaCO3 to kg C (1e-3*12/100)
        tai_calmass = tai_calmass + tmp*12.e-5
        tmpip(:) = 0.
        do k=1,ibmax
          do ip=1,ipsed
            tmpip(ip) = tmpip(ip)+buried_mass(k,ip)*buried_calfrac(k,ip)
          enddo
        enddo
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areatot (tmpij, dmsk, tmp)
!       convert from g CaCO3 to kg C (1e-3*12/100)
        tai_calmass_bur = tai_calmass_bur + tmp*12.e-5
        tmpip(1:ipmax) = co3_p(1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_co3 = tai_co3 + tmp
        tmpip(1:ipmax) = csat(1:ipmax)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_co3sat = tai_co3sat + tmp
!       convert umol C s-1 to kg C s-1
        tai_weathflx = tai_weathflx + weathflx*12.e-9
!       convert umol to Pg
        tai_carblith = tai_carblith + carblith*12.e-21
        tmpip(1:ipmax) = rain_cal_p(1:ipmax)
     &                   /(rain_org_p(1:ipmax) + 1.e-20)
        call unloadsed (ipmax, tmpip, imt, jmt, map_sed, tmpij)
        call areaavg (tmpij, dmsk, tmp)
        tai_rainr = tai_rainr + tmp

      elseif (iflag .eq. 2 .and. ntatis .ne. 0) then

!       average
        rntatis = 0.
        if (ntatis .gt. 0.) rntatis = 1./float(ntatis)
        tai_ttrcal = tai_ttrcal*rntatis
        tai_rain_cal = tai_rain_cal*rntatis
        tai_cal = tai_cal*rntatis
        tai_calmass = tai_calmass*rntatis
        tai_calmass_bur = tai_calmass_bur*rntatis
        tai_co3 = tai_co3*rntatis
        tai_co3sat = tai_co3sat*rntatis
        tai_weathflx = tai_weathflx*rntatis
        tai_carblith = tai_carblith*rntatis
        tai_rainr = tai_rainr*rntatis
      endif

      return
      end
# endif

!=======================================================================
      subroutine unloadsed (kd, dk, id, jd, map, dij)

!-----------------------------------------------------------------------
!     load 1d sediment array into 2d array
!-----------------------------------------------------------------------

      implicit none

      integer i, id, j, jd, k, kd, map(id,jd)

      real dk(kd), dij(id,jd)

      dij(:,:) = 0.
      dij(:,:) = 0.
      do j=1,jd
        do i=1,id
          k = map(i,j)
          if (k .ge. 1 .and. k .le. kd) dij(i,j) = dk(k)
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine loadsed (kd, dk, id, jd, map, dij)

!-----------------------------------------------------------------------
!     load 2d array into 1d sediment array
!-----------------------------------------------------------------------

      implicit none

      integer i, id, j, jd, k, kd, map(id,jd)

      real dk(kd), dij(id,jd)

      dk(:) = 0.
      do j=1,jd
        do i=1,id
          k = map(i,j)
          if (k .ge. 1 .and. k .le. kd) dk(k) = dij(i,j)
        enddo
      enddo
#endif

      return
      end
