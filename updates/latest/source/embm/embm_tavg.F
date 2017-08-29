      subroutine embm_tavg_def (fname, imt, jmt, nat, ncat, xt, yt
     &,                         calendar, expnam, runstamp, mapat)

#if defined O_embm
!=======================================================================
!     definition routine for atmospheric time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapat        = tracer map
!=======================================================================

      implicit none

      integer iou, j, n, imt, jmt, nat, ncat, igs, ige, ig, jgs, jge
      integer jg, it(10), iu(10), id_time, id_xt, id_xu, id_yt, id_yu
      integer id_cat, id_cat_e, id_xt_e, id_xu_e, id_yt_e, id_yu_e

      character(*) :: fname, calendar, expnam, runstamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), yt(jmt)
      real c0, c1, c100, c500, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('longitude', iou, ig, id_xt)
      call defdim ('latitude', iou, jg, id_yt)
      call defdim ('longitude_V', iou, ig, id_xu)
      call defdim ('latitude_V', iou, jg, id_yu)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('longitude_V_edges', iou, ig+1, id_xu_e)
      call defdim ('latitude_V_edges', iou, jg+1, id_yu_e)
# if defined O_ice_cpts && defined O_ice
      call defdim ('cat', iou, ncat, id_cat)
      call defdim ('cat_edges', iou, ncat+1, id_cat_e)
# endif

!-----------------------------------------------------------------------
!       define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
# if defined O_units_time_years
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'years since 1-1-1')
#  else
     &, 'time', 'time', 'years since 0-1-1')
#  endif
# else
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'days since 1-1-1')
#  else
     &, 'time', 'time', 'days since 0-1-1')
#  endif
# endif
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_avgper', iou, 1, it, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','day')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      call defvar ('G_dxt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'width t grid', ' ', 'm')
      it(1) = id_yt
      call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      call defvar ('G_dyt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'height t grid', ' ', 'm')
      it(1) = id_xu
      call defvar ('longitude_V', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude', 'longitude', 'degrees_east')
      call defvar ('G_dxu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'width u grid', ' ', 'm')
      it(1) = id_yu
      call defvar ('latitude_V', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude', 'latitude', 'degrees_north')
      call defvar ('G_dyu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'height u grid', ' ', 'm')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'longitude', 'degrees_east')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'longitude', 'degrees_east')
# if defined O_ice_cpts && defined O_ice
      it(1) = id_cat
      call defvar ('cat', iou, 1, it, c0, c0, 'Z', 'F'
     &, 'ice category', ' ',' ')
      it(1) = id_cat_e
      call defvar ('cat_edges', iou, 1, it, c0, c0, 'Z', 'F'
     &, 'ice category edges', ' ',' ')
# endif

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('L_elev', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'land elevation and ocean depth', 'surface_altitude', 'm')
      call defvar ('L_rivers', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'river basin number', ' ' ,'1')
      call defvar ('G_mskt', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'ocean mask', ' ' ,'1')
      call defvar ('G_mskhr', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'horizontal region mask', ' ' ,'1')
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degrees_east')
      call defvar ('G_latU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid longitude', 'longitude', 'degrees_east')
      call defvar ('G_areaT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid area', ' ', 'm2')
      call defvar ('G_areaU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid area', ' ', 'm2')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('A_slat', iou, 3, it, -c100, c500, ' ', 'F'
     &,     'sea level atmospheric temperature'
# if defined O_units_temperature_Celsius
     &,     'air_temperature', 'C')
# else
     &,     'air_temperature', 'K')
# endif
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('A_shum', iou, 3, it, -c100, c100, ' ', 'F'
     &,     'atmospheric surface specific humidity'
     &,     'specific_humidity', '1')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('A_co2', iou, 3, it, c0, c1e6, ' ', 'F'
     &,     'atmospheric co2', 'atmospheric_co2', 'ppm')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('A_tracer'//trim(a3), iou ,3, it, -c1e6, c1e6
     &,     ' ', 'F', 'tracer '//trim(a3)
     &,     'tracer_'//trim(a3), 'unknown')
        endif
      enddo
      call defvar ('A_sat', iou, 3, it, -c100, c500, ' ', 'F'
# if defined O_units_temperature_Celsius
     &, 'atmospheric surface temperature', 'air_temperature', 'C')
# else
     &, 'atmospheric surface temperature', 'air_temperature', 'K')
# endif
      call defvar ('F_precip', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation (includes snow in water equivalent)'
     &, 'precipitation_flux', 'kg m-2 s-1')
      call defvar ('F_snow', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation as snow'
     &, 'snowfall_flux', 'kg m-2 s-1')
      call defvar ('F_evap', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'upward evaporation plus sublimation'
     &, 'water_evaporation_flux', 'kg m-2 s-1')
      call defvar ('F_rivdis', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'river discharge', 'river_discharge_flux', 'kg m-2 s-1')
      call defvar ('F_virtual', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'normalized virtual flux', 'virtual_flux', 'm s-1')
      call defvar ('F_outlwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'outgoing longwave at top of atmosphere'
     &, 'toa_outgoing_longwave_flux', 'W m-2')
      call defvar ('F_uplwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface net upward longwave'
     &, 'surface_net_upward_longwave_flux', 'W m-2')
      call defvar ('F_upsens', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward sensible heat'
     &, 'surface_upward_sensible_heat_flux', 'W m-2')
      call defvar ('F_dnswr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface net downward shortwave (absorbed)'
     &, 'surface_net_downward_shortwave_flux', 'W/m^2')
      call defvar ('F_solins', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'incoming shortwave at top of atmosphere'
     &, 'toa_incoming_shortwave_flux', 'W m-2')
      call defvar ('F_outswr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'outgoing shortwave at top of atmosphere'
     &, 'toa_outgoing_shortwave_flux', 'W m-2')
      call defvar ('F_netrad', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'net radiation at top of atmosphere'
     &, 'toa_net_radiation_flux', 'W m-2')
      call defvar ('A_albplt', iou, 3, it, c0, c1, ' ', 'F'
     &, 'planetary albedo', 'planetary_albedo', '1')
      call defvar ('A_albatm', iou, 3, it, c0, c1, ' ', 'F'
     &, 'atmospheric albedo', ' ', '1')
      call defvar ('A_albsur', iou, 3, it, c0, c1, ' ', 'F'
     &, 'surface albedo', 'surface_albedo', '1')
      call defvar ('F_runoff', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'soil runoff', 'runoff_flux', 'kg m-2 s-1')
      call defvar ('A_windspd', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'surface wind speed', 'wind_speed', 'm s-1')
# if defined O_save_embm_wind
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('A_windtX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward wind for advection of temperature'
     &,     'eastward_wind', 'm s-1')
          call defvar ('A_windtY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward wind for advection of temperature'
     &,     'northward_wind', 'm s-1')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('A_windqX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward wind for advection of humidity'
     &,     'eastward_wind', 'm s-1')
          call defvar ('A_windqY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward wind for advection of humidity'
     &,     'northward_wind', 'm s-1')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('A_windcX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward wind for advection of carbon'
     &,     'eastward_wind', 'm s-1')
          call defvar ('A_windcY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward wind for advection of carbon'
     &,     'northward_wind', 'm s-1')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('A_wind'//trim(a3)//'X', iou ,3, it, -c1e6, c1e6
     &,     ' ', 'F', 'eastward wind for tracer '//trim(a3)
     &,     'eastward_wind_for_tracer_'//trim(a3), 'unknown')
          call defvar ('A_wind'//trim(a3)//'Y', iou ,3, it, -c1e6, c1e6
     &,     ' ', 'F', 'northward wind for tracer '//trim(a3)
     &,     'northward_wind_for_tracer_'//trim(a3), 'unknown')
        endif
      enddo
# endif
# if defined O_embm_awind
      call defvar ('A_awindX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward wind  anomaly', 'eastward_wind_anomaly', 'm s-1')
      call defvar ('A_awindY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward wind anomaly', 'northward_wind_anomaly', 'm s-1')
      call defvar ('A_avgslat', iou, 3, iu, -c100, c500, ' ', 'F'
     &, 'running average sea level atmospheric temperature'
#  if defined O_units_temperature_Celsius
     &, 'air_temperature', 'C')
#  else
     &, 'air_temperature', 'K')
#  endif
      call defvar ('A_apress', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'anomalous surface pressure', 'pressure', 'Pa')
# endif
      call defvar ('L_soilmois', iou, 3, it, c0, c1e3, ' ', 'F'
     &, 'soil moisture', 'soil_moisture_content', 'kg m-2')
      call defvar ('L_tempsur', iou, 3, it, -c100, c500, ' ', 'F'
# if defined O_units_temperature_Celsius
     &, 'land surface temperature', 'surface_temperature', 'C')
# else
     &, 'land surface temperature', 'surface_temperature', 'K')
# endif
# if defined O_save_flxadj
      call defvar ('F_adjtemp', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'surface downward heat flux adjustment'
     &, ' ', 'W m-2')
      call defvar ('F_adjsal', iou,3, it, -c100, c100, ' ', 'F'
     &, 'surface downward salt flux adjustment'
     &, ' ', 'kg m-2 s-1')
# endif
# if defined O_save_embm_diff
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('A_difftX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward diffusion for temperature'
     &,     'eastward_diffusion', 'm2 s-1')
          call defvar ('A_difftY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward diffusion for temperature'
     &,     'northward_diffusion', 'm2 s-1')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('A_diffqX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward diffusion for humidity'
     &,     'eastward_diffusion', 'm2 s-1')
          call defvar ('A_diffqY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward diffusion for humidity'
     &,     'northward_diffusion', 'm2 s-1')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('A_diffcX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'eastward diffusion for carbon'
     &,     'eastward_diffusion', 'm2 s-1')
          call defvar ('A_diffcY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &,     'northward diffusion for carbon'
     &,     'northward_diffusion', 'm2 s-1')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('A_diff'//trim(a3)//'X', iou ,3, it, -c1e6, c1e6
     &,     ' ', 'F', 'eastward diffusion for tracer '//trim(a3)
     &,     'eastward_diffusion_for_tracer_'//trim(a3), 'unknown')
          call defvar ('A_diff'//trim(a3)//'Y', iou ,3, it, -c1e6, c1e6
     &,     ' ', 'F', 'northward diffusion for tracer '//trim(a3)
     &,     'northward_diffusion_for_tracer_'//trim(a3), 'unknown')
        endif
      enddo
# endif
# if defined O_landice_data
      call defvar ('L_icefra', iou, 3, it, c0, c100, ' ', 'F'
     &, 'ice sheet area fraction', 'land_ice_area_fraction', '1')
      call defvar ('L_icethk', iou, 3, it, c0, c1e6, ' ', 'F'
     &, 'ice sheet thickness anomaly', ' ', 'm')
# endif
# if defined O_carbon && defined O_carbon_co2_2d
      call defvar ('F_co2', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'surface downward flux of co2', ' ', 'kg m-2 s-1')
#  if defined O_co2emit_data || defined O_co2emit_data_transient
      call defvar ('F_co2emit', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'emissions of co2', ' ', 'kg m-2 s-1')
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      call defvar ('A_sulphfor', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward shortwave flux due to sulphates'
     &, 'surface_upward_shortwave_flux_due_to_sulphates', 'W m-2')
# endif
# if defined O_ice
      call defvar ('O_icetemp', iou, 3, it, -c100, c500, ' ', 'F'
     &, 'surface ice temperature', 'surface_temperature'
#  if defined O_units_temperature_Celsius
     &, 'C')
#  else
     &, 'K')
#  endif
      call defvar ('O_icethk', iou, 3, it, c0, c1e6, ' ', 'F'
     &, 'ice thickness', 'sea_ice_thickness', 'm')
      call defvar ('O_icefra', iou, 3, it, c0, c100, ' ', 'F'
     &, 'ice area fraction (includes land ice area fraction)'
     &, 'sea_ice_area_fraction', '1')
      call defvar ('O_snothk', iou, 3, it, -0.01, c1e6, ' ', 'F'
     &, 'surface snow thickness', 'surface_snow_thickness', 'm')
# endif
# if defined O_ice_evp && defined O_ice
      call defvar ('O_icevelX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward ice velocity', 'eastward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('O_icevelY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward ice velocity', 'northward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('O_iceintX', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'eastward ice interaction'
     &, 'downward_eastward_stress_at_sea_ice_base', 'Pa')
      call defvar ('O_iceintY', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'northward ice interaction'
     &, 'downward_northward_stress_at_sea_ice_base', 'Pa')
# endif
# if defined O_ice_cpts && defined O_ice

!-----------------------------------------------------------------------
!       define time dependent 4d data (x,y,c,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_cat
      iu(3) = id_cat
      it(4) = id_time
      iu(4) = id_time
      call defvar ('O_icetempc', iou, 4, it, -c100, c500, ' ', 'F'
     &, 'multi-category ice temperature'
#  if defined O_units_temperature_Celsius
     &, 'surface_sea_ice_temperature', 'C')
#  else
     &, 'surface_sea_ice_temperature', 'K')
#  endif
      call defvar ('O_icethkc', iou, 4, it, c0, c1e6, ' ', 'F'
     &, 'multi-category ice thickness', 'sea_ice_thickness', 'm')
      call defvar ('O_icefrac', iou, 4, it, c0, c100, ' ', 'F'
     &, 'multi-category ice area', 'sea_ice_area_fraction', '1')
      call defvar ('O_snothkc', iou, 4, it, c0, c1e6, ' ', 'F'
     &, 'multi-category snow thickness', 'surface_snow_thickness'
     &, 'm')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine embm_tavg_out (fname, ids, ide, jds, jde, imt, jmt, nat
     &,                         ncat, xt, yt, xu, yu, dxt, dyt, dxu, dyu
     &,                         avgper, time, stamp, mapat, at, sat
     &,                         precip, evap, disch, vflux, outlwr
     &,                         uplwr, upsens, dnswr, solins, outswr
     &,                         netrad, p_alb, a_alb, s_alb, elev, psno
     &,                         ws, runoff
# if defined O_save_embm_wind
     &,                         wx, wy
# endif
# if defined O_embm_awind
     &,                         awx, awy, rtbar, apress
# endif
     &,                         soilm, surf
# if defined O_ice
     &,                         tice, hice, aice, hsno
# endif
# if defined O_ice_cpts && defined O_ice
     &,                         Ts, heff, A, hseff
# endif
# if defined O_ice_evp && defined O_ice
     &,                         uice, vice, xint, yint
# endif
     &,                         tlat, tlon, ulat, ulon, tgarea, ugarea
# if defined O_save_flxadj
     &,                         flxadj_t, flxadj_s
# endif
# if defined O_save_embm_diff
     &,                         dn, de
# endif
# if defined O_landice_data
     &,                         aicel, hicel
# endif
# if defined O_carbon_co2_2d
     &,                         flux_co2
#  if defined O_co2emit_data || defined O_co2emit_data_transient
     &,                         co2emit
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
     &,                         sulph
# endif
     &,                         tmsk, mskhr, nriv, ntrec)
!=======================================================================
!     output routine for atmospheric time averages

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths

!     time         = time in years
!     stamp        = time stamp
!     at, ...      = data to be written

!   outputs:
!     ntrec        = number of time record in file
!=======================================================================

      implicit none

      integer iou, j, k, ln, n, ntrec, imt, jmt, nat, ncat
      integer ids, ide, jds, jde, igs, ige, ig, jgs, jge, jg
      integer ils, ile, jls, jle, ib(10), ic(10)
      integer mskhr(ids:ide,jds:jde), nriv(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      character(*) :: fname, stamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), xu(imt), yt(jmt), yu(jmt)
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt)
      real avgper, at(ids:ide,jds:jde,nat)
      real sat(ids:ide,jds:jde), precip(ids:ide,jds:jde)
      real evap(ids:ide,jds:jde), disch(ids:ide,jds:jde)
      real vflux(ids:ide,jds:jde), outlwr(ids:ide,jds:jde)
      real uplwr(ids:ide,jds:jde), upsens(ids:ide,jds:jde)
      real dnswr(ids:ide,jds:jde), solins(ids:ide,jds:jde)
      real outswr(ids:ide,jds:jde), netrad(ids:ide,jds:jde)
      real p_alb(ids:ide,jds:jde), a_alb(ids:ide,jds:jde)
      real s_alb(ids:ide,jds:jde), elev(ids:ide,jds:jde)
      real psno(ids:ide,jds:jde), ws(ids:ide,jds:jde)
      real runoff(ids:ide,jds:jde)
# if defined O_save_embm_wind
      real  wx(ids:ide,jds:jde,nat), wy(ids:ide,jds:jde,nat)
# endif
# if defined O_embm_awind
      real awx(ids:ide,jds:jde), awy(ids:ide,jds:jde)
      real rtbar(ids:ide,jds:jde), apress(ids:ide,jds:jde)
# endif
      real soilm(ids:ide,jds:jde), surf(ids:ide,jds:jde)
# if defined O_ice
      real tice(ids:ide,jds:jde), hice(ids:ide,jds:jde)
      real aice(ids:ide,jds:jde), hsno(ids:ide,jds:jde)
# endif
# if defined O_ice_cpts && defined O_ice
      real cat(ncat), cat_e(ncat+1)
      real Ts(ids:ide,jds:jde,ncat), heff(ids:ide,jds:jde,ncat)
      real A(ids:ide,jds:jde,ncat), hseff(ids:ide,jds:jde,ncat)
# endif
# if defined O_ice_evp && defined O_ice
      real uice(ids:ide,jds:jde), vice(ids:ide,jds:jde)
      real xint(ids:ide,jds:jde), yint(ids:ide,jds:jde)
# endif
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
      real tgarea(ids:ide,jds:jde), ugarea(ids:ide,jds:jde)
# if defined O_save_flxadj
      real flxadj_t(ids:ide,jds:jde), flxadj_s(ids:ide,jds:jde)
# endif
# if defined O_save_embm_diff
      real dn(ids:ide,jds:jde,nat), de(ids:ide,jds:jde,nat)
# endif
# if defined O_landice_data
      real aicel(ids:ide,jds:jde), hicel(ids:ide,jds:jde)
# endif
# if defined O_carbon_co2_2d
      real flux_co2(ids:ide,jds:jde)
#  if defined O_co2emit_data || defined O_co2emit_data_transient
      real co2emit(ids:ide,jds:jde)
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      real sulph(ids:ide,jds:jde)
# endif
      real tmsk(ids:ide,jds:jde)
      real time, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), tmpmask(ids:ide,jds:jde)
      real c0, c1, c10, c100, c1e3, c1e4, p1, C2K, cal2J
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      p1 = 0.1
      C2K = 273.15
      cal2J = 2.389e-05

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, time, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)

      if (ntrec .eq. 1) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )

        ib(1) = 1
        ic(1) = ig
        tmpi(igs:ige) = xt(igs:ige)
        call putvara ('longitude', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxt(igs:ige)
        call putvara ('G_dxt', iou, ig, ib, ic, tmpi, c100, c0)
        tmpi(igs:ige) = xu(igs:ige)
        call putvara ('longitude_V', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxu(igs:ige)
        call putvara ('G_dxu', iou, ig, ib, ic, tmpi, c100, c0)

        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('latitude', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyt(jgs:jge)
        call putvara ('G_dyt', iou, jg, ib, ic, tmpj, c100, c0)
        tmpj(jgs:jge) = yu(jgs:jge)
        call putvara ('latitude_V', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyu(jgs:jge)
        call putvara ('G_dyu', iou, jg, ib, ic, tmpj, c100, c0)

        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('longitude_edges', iou, ig+1, ib, ic, tmpie
     &,   c1, c0)
        call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xu_e(igs:ige+1)
        call putvara ('longitude_V_edges', iou, ig+1, ib, ic, tmpie
     &,   c1, c0)

        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('latitude_edges', iou, jg+1, ib, ic, tmpje
     &,   c1, c0)
        call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
        call putvara ('latitude_V_edges', iou, jg+1, ib, ic, tmpje
     &,   c1, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpie )
        deallocate ( tmpje )
# if defined O_ice_cpts && defined O_ice

        do n=1,ncat
          cat(n) = float(n)
          cat_e(n) = float(n) - 0.5
        enddo
        cat_e(ncat+1) = cat(ncat) + 0.5
        ic(1) = ncat
        call putvara ('cat', iou, ncat, ib, ic, cat, 1., 0.)
        ic(1) = ncat+1
        call putvara ('cat_edges', iou, ncat+1, ib, ic, cat_e, 1., 0.)
# endif

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = elev(ils:ile,jls:jle)
        call putvara ('L_elev', iou, ln, ib, ic, tmpij, c100, c0)
        tmpij(ils:ile,jls:jle) = nriv(ils:ile,jls:jle)
        call putvara ('L_rivers', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tmsk(ils:ile,jls:jle)
        call putvara ('G_mskt', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = mskhr(ils:ile,jls:jle)
        call putvara ('G_mskhr', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('G_latT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('G_lonT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulat(ils:ile,jls:jle)
        call putvara ('G_latU', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulon(ils:ile,jls:jle)
        call putvara ('G_lonU', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tgarea(ils:ile,jls:jle)
        call putvara ('G_areaT', iou, ln, ib, ic, tmpij, c1e4, c0)
        tmpij(ils:ile,jls:jle) = ugarea(ils:ile,jls:jle)
        call putvara ('G_areaU', iou, ln, ib, ic, tmpij, c1e4, c0)

      endif

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
# if defined O_units_temperature_Celsius
          call putvara ('A_slat', iou, ln, ib, ic, tmpij, c1, c0)
# else
          call putvara ('A_slat', iou, ln, ib, ic, tmpij, c1, -C2K)
# endif
        elseif (trim(mapat(n)) .eq. 'shum') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_shum', iou, ln, ib, ic, tmpij, c1, c0)
        elseif (trim(mapat(n)) .eq. 'co2') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_co2', iou, ln, ib, ic, tmpij, c1, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_tracer'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
        endif
      enddo
      tmpij(ils:ile,jls:jle) = sat(ils:ile,jls:jle)
# if defined O_units_temperature_Celsius
      call putvara('A_sat', iou, ln, ib, ic, tmpij, c1, c0)
# else
      call putvara('A_sat', iou, ln, ib, ic, tmpij, c1, -C2K)
# endif
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call putvara ('F_precip', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = psno(ils:ile,jls:jle)
      call putvara ('F_snow', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = evap(ils:ile,jls:jle)
      call putvara ('F_evap', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = disch(ils:ile,jls:jle)
      call putvara ('F_rivdis', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = vflux(ils:ile,jls:jle)
      call putvara ('F_virtual', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = outlwr(ils:ile,jls:jle)
      call putvara ('F_outlwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = uplwr(ils:ile,jls:jle)
      call putvara ('F_uplwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = upsens(ils:ile,jls:jle)
      call putvara ('F_upsens', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = dnswr(ils:ile,jls:jle)
      call putvara ('F_dnswr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = solins(ils:ile,jls:jle)
      call putvara ('F_solins', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = outswr(ils:ile,jls:jle)
      call putvara ('F_outswr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = netrad(ils:ile,jls:jle)
      call putvara ('F_netrad', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpmask(:,:) = 1.
      where (p_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = p_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albplt', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpmask(:,:) = 1.
      where (a_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = a_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albatm', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpmask(:,:) = 1.
      where (s_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = s_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = runoff(ils:ile,jls:jle)
      call putvara ('F_runoff', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = ws(ils:ile,jls:jle)
      call putvara ('A_windspd', iou, ln, ib, ic, tmpij, c100, c0)
# if defined O_save_embm_wind
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          tmpij(ils:ile,jls:jle) = wx(ils:ile,jls:jle,n)
          call putvara ('A_windtX', iou, ln, ib, ic, tmpij, c100, c0)
          tmpij(ils:ile,jls:jle) = wy(ils:ile,jls:jle,n)
          call putvara ('A_windtY', iou, ln, ib, ic, tmpij, c100, c0)
        elseif (trim(mapat(n)) .eq. 'shum') then
          tmpij(ils:ile,jls:jle) = wx(ils:ile,jls:jle,n)
          call putvara ('A_windqX', iou, ln, ib, ic, tmpij, c100, c0)
          tmpij(ils:ile,jls:jle) = wy(ils:ile,jls:jle,n)
          call putvara ('A_windqY', iou, ln, ib, ic, tmpij, c100, c0)
        elseif (trim(mapat(n)) .eq. 'co2') then
          tmpij(ils:ile,jls:jle) = wx(ils:ile,jls:jle,n)
          call putvara ('A_windcX', iou, ln, ib, ic, tmpij, c100, c0)
          tmpij(ils:ile,jls:jle) = wy(ils:ile,jls:jle,n)
          call putvara ('A_windcY', iou, ln, ib, ic, tmpij, c100, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          tmpij(ils:ile,jls:jle) = wx(ils:ile,jls:jle,n)
          call putvara ('A_wind'//trim(a3)//'X', iou, ln, ib, ic, tmpij
     &,     c1, c0)
          tmpij(ils:ile,jls:jle) = wy(ils:ile,jls:jle,n)
          call putvara ('A_wind'//trim(a3)//'Y', iou, ln, ib, ic, tmpij
     &,     c1, c0)
        endif
      enddo
# endif
# if defined O_embm_awind
      tmpij(ils:ile,jls:jle) = awx(ils:ile,jls:jle)
      call putvara ('A_awindX', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) =  awy(ils:ile,jls:jle)
      call putvara ('A_awindY', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call putvara ('A_avgslat', iou, ln, ib, ic, tmpij, c1, c0)
#  if defined O_units_temperature_Celsius
      call putvara ('A_avgslat', iou, ln, ib, ic, tmpij, c1, c0)
#  else
      call putvara ('A_avgslat', iou, ln, ib, ic, tmpij, c1, -C2K)
#  endif
      tmpij(ils:ile,jls:jle) =  apress(ils:ile,jls:jle)
      call putvara ('A_apress', iou, ln, ib, ic, tmpij, c10, c0)
# endif
      tmpmask(:,:) = 0.
      where (tmsk(:,:) .lt. 0.5) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle)
      call putvaramsk ('L_soilmois', iou, ln, ib, ic, tmpij, tmpijm
     &, p1, c0)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call putvaramsk ('L_tempsur', iou, ln, ib, ic, tmpij, tmpijm
# if defined O_units_temperature_Celsius
     &, c1, c0)
# else
     &, c1, -C2K)
# endif
# if defined O_ice
      tmpmask(:,:) = 0.
      where (hice(:,:) .gt. 0) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = tice(ils:ile,jls:jle)
      call putvaramsk ('O_icetemp', iou, ln, ib, ic, tmpij, tmpijm
#  if defined O_units_temperature_Celsius
     &, c1, c0)
#  else
     &, c1, -C2K)
#  endif
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle)
      call putvara ('O_icethk', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle)
      call putvara ('O_icefra', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle)
      call putvara ('O_snothk', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if defined O_ice_evp && defined O_ice
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call putvara ('O_icevelX', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call putvara ('O_icevelY', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = xint(ils:ile,jls:jle)
      call putvara ('O_iceintX', iou, ln, ib, ic, tmpij, c10, c0)
      tmpij(ils:ile,jls:jle) = yint(ils:ile,jls:jle)
      call putvara ('O_iceintY', iou, ln, ib, ic, tmpij, c10, c0)
# endif
# if defined O_save_flxadj
      tmpij(ils:ile,jls:jle) = flxadj_t(ils:ile,jls:jle)
      call putvara('F_adjtemp', iou, ln, ib, ic, tmpij, cal2J, c0)
      tmpij(ils:ile,jls:jle) = flxadj_s(ils:ile,jls:jle)
      call putvara('F_adjsal', iou, ln, ib, ic, tmpij, p1, c0)
# endif
# if defined O_save_embm_diff
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          tmpij(ils:ile,jls:jle) = de(ils:ile,jls:jle,n)
          call putvara ('A_difftX', iou, ln, ib, ic, tmpij, c1e4, c0)
          tmpij(ils:ile,jls:jle) = dn(ils:ile,jls:jle,n)
          call putvara ('A_difftY', iou, ln, ib, ic, tmpij, c1e4, c0)
        elseif (trim(mapat(n)) .eq. 'shum') then
          tmpij(ils:ile,jls:jle) = de(ils:ile,jls:jle,n)
          call putvara ('A_diffqX', iou, ln, ib, ic, tmpij, c1e4, c0)
          tmpij(ils:ile,jls:jle) = dn(ils:ile,jls:jle,n)
          call putvara ('A_diffqY', iou, ln, ib, ic, tmpij, c1e4, c0)
        elseif (trim(mapat(n)) .eq. 'co2') then
          tmpij(ils:ile,jls:jle) = de(ils:ile,jls:jle,n)
          call putvara ('A_diffcX', iou, ln, ib, ic, tmpij, c1e4, c0)
          tmpij(ils:ile,jls:jle) = dn(ils:ile,jls:jle,n)
          call putvara ('A_diffcY', iou, ln, ib, ic, tmpij, c1e4, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          tmpij(ils:ile,jls:jle) = de(ils:ile,jls:jle,n)
          call putvara('A_diff'//trim(a3)//'X', iou, ln, ib, ic, tmpij
     &,     c1e4, c0)
          tmpij(ils:ile,jls:jle) = dn(ils:ile,jls:jle,n)
          call putvara('A_diff'//trim(a3)//'Y', iou, ln, ib, ic, tmpij
     &,     c1e4, c0)
        endif
      enddo
# endif
# if defined O_landice_data
      tmpij(ils:ile,jls:jle) = aicel(ils:ile,jls:jle)
      call putvara('L_icefra', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hicel(ils:ile,jls:jle)
      call putvara('L_icethk', iou, ln, ib, ic, tmpij, c100, c0)
# endif
# if O_carbon_co2_2d
      tmpij(ils:ile,jls:jle) = flux_co2(ils:ile,jls:jle)
      call putvara('F_co2', iou, ln, ib, ic, tmpij, c10, c0)
#  if defined O_co2emit_data || defined O_co2emit_data_transient
      tmpij(ils:ile,jls:jle) = co2emit(ils:ile,jls:jle)
      call putvara('F_co2emit', iou, ln, ib, ic, tmpij, c10, c0)
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      tmpij(ils:ile,jls:jle) = sulph(ils:ile,jls:jle)
      call putvara('A_sulphfor', iou, ln, ib, ic, tmpij, c1e3, c0)
# endif
# if defined O_ice_cpts && defined O_ice

!-----------------------------------------------------------------------
!     write 4d data (x,y,c,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpmask(:,:) = 0.
      where (tmsk(:,:) .ge. 0.5) tmpmask(:,:) = 1.
      do k=1,ncat
        ib(3) = k
        tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
        tmpij(ils:ile,jls:jle) = Ts(ils:ile,jls:jle,k)
        call putvaramsk ('O_icetempc', iou, ln, ib, ic, tmpij, tmpijm
#  if defined O_units_temperature_Celsius
     &,   c1, c0)
#  else
     &,   c1, -C2K)
#  endif
        tmpij(ils:ile,jls:jle) = heff(ils:ile,jls:jle,k)
        call putvaramsk ('O_icethkc', iou, ln, ib, ic, tmpij, tmpijm
     &,   c100, c0)
        tmpij(ils:ile,jls:jle) = A(ils:ile,jls:jle,k)
        call putvaramsk ('O_icefrac', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        tmpij(ils:ile,jls:jle) = hseff(ils:ile,jls:jle,k)
        call putvaramsk ('O_snothkc', iou, ln, ib, ic, tmpij, tmpijm
     &,   c100, c0)
      enddo
# endif

      deallocate ( tmpij )
      deallocate ( tmpijm )
#endif

      return
      end
