      subroutine embm_tsi_def (fname, calendar, expnam, runstamp)

#if defined O_embm
!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou

      real c0, c1, c100, c500, c1e3, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

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
      id(1) = id_time

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'
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
      call defvar ('T_avgper', iou, 1, id, c0, c0, ' ', 'F'
     &,   'averaging period', ' ','day')
      call defvar ('A_sat', iou, 1, id, -c100, c500, ' ', 'F'
# if defined O_units_temperature_Celsius
     &,   'global average surface air temperature', ' ', 'C')
# else
     &,   'global average surface air temperature', ' ', 'K')
# endif
      call defvar ('A_shum', iou, 1, id, -c100, c100, ' ', 'F'
     &,   'global average surface specific humidity', ' ', '1')
      call defvar ('F_precip', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_evap', iou, 1, id, -c100, c100, ' '
     &,   'F', 'global average evaporation', ' ','kg m-2 s-1')
      call defvar ('A_co2', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global average CO2 concentration', ' ','ppm')
# if defined O_save_carbon_totals
      call defvar ('A_totcarb', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global total atmosphere carbon', ' ', 'Pg')
      call defvar ('A_carbemit', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global total accumulated carbon emissions', ' ', 'Pg')
      call defvar ('F_co2emit', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global total CO2 emissions', ' ','Pg yr-1')
# else
      call defvar ('F_co2emit', iou, 1, id, c0, c1e3, ' '
     &,   'F', 'global total CO2 emissions', ' ','kg s-1')
# endif
      call defvar ('A_maxiter', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'maximum atmosphere solver iterations', ' ','1')
# if defined O_carbon_13
      call defvar ('A_dc13', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &,   'global average delta C13 concentration', ' ','permil')
# endif
# if defined O_carbon_14
      call defvar ('A_dc14', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &,   'global average delta C14 concentration', ' ','permil')
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call defvar ('A_cfc11', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'global average CFC11 concentration', ' ','ppt')
      call defvar ('A_cfc12', iou, 1, id, -c1e3, c1e3, ' '
     &,   'F', 'global average CFC12 concentration', ' ','ppt')
# endif
# if defined O_ice
      call defvar ('O_snovol', iou, 1, id, c0, c1e20, ' '
     &,   'F', 'global snow volume', ' ', 'm3')
      call defvar ('O_icevol', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice volume', ' ', 'm3')
      call defvar ('O_icearea', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global sea ice area', ' ', 'm2')
#   if defined O_landice_data_transient || defined O_ism
      call defvar ('L_icevol', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'global land ice volume', ' ', 'm3')
      call defvar ('L_icearea', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'global land ice area', ' ', 'm2')
#  endif
# endif
# if defined O_tai_ns
      call defvar ('A_satN', iou, 1, id, -c100, c500, ' '
#  if defined O_units_temperature_Celsius
     &,   'F', 'NH average surface air temperature', ' ', 'C')
      call defvar ('A_satS', iou, 1, id, -c100, c500, ' '
     &,   'F', 'SH average surface air temperature', ' ', 'C')
#  else
     &,   'F', 'NH average surface air temperature', ' ', 'K')
      call defvar ('A_satS', iou, 1, id, -c100, c500, ' '
     &,   'F', 'SH average surface air temperature', ' ', 'K')
#  endif
      call defvar ('A_shumN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average surface specific humidity', ' ', '1')
      call defvar ('A_shumS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average surface specific humidity', ' ', '1')
      call defvar ('F_precipN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_precipS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average precipitation', ' ','kg m-2 s-1')
      call defvar ('F_evapN', iou, 1, id, -c100, c100, ' '
     &,   'F', 'NH average evaporation', ' ','kg m-2 s-1')
      call defvar ('F_evapS', iou, 1, id, -c100, c100, ' '
     &,   'F', 'SH average evaporation', ' ','kg m-2 s-1')
#  if defined O_ice
      call defvar ('O_snovolN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH snow volume', ' ', 'm3')
      call defvar ('O_snovolS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH snow volume', ' ', 'm3')
      call defvar ('O_icevolN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice volume', ' ', 'm3')
      call defvar ('O_icevolS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice volume', ' ', 'm3')
      call defvar ('O_iceareaN', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'NH sea ice area', ' ', 'm2')
      call defvar ('O_iceareaS', iou, 1, id, c0, c1e20
     &,   ' ', 'F', 'SH sea ice area', ' ', 'm2')
#   if defined O_landice_data_transient || defined O_ism
      call defvar ('L_icevolN', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'NH land ice volume', ' ', 'm3')
      call defvar ('L_icevolS', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'SH land ice volume', ' ', 'm3')
      call defvar ('L_iceareaN', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'NH land ice area', ' ', 'm2')
      call defvar ('L_iceareaS', iou, 1, id, c0
     &,  c1e20, ' ', 'F', 'SH land ice area', ' ', 'm2')
#   endif
#  endif
# endif
# if defined O_tai_lo
#  if defined O_units_temperature_Celsius
      call defvar ('A_satL', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over land', ' ', 'C')
      call defvar ('A_satO', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over ocean', ' ', 'C')
#  else
      call defvar ('A_satL', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over land', ' ', 'K')
      call defvar ('A_satO', iou, 1, id, -c100, c500, ' '
     &,   'F', 'average surface air temperature over ocean', ' ', 'K')
#  endif
      call defvar ('F_precipL', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average precipitation over land', ' ','kg m-2 s-1')
      call defvar ('F_precipO', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average precipitation over ocean', ' ','kg m-2 s-1')
      call defvar ('F_evapL', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average evaporation over land', ' ','kg m-2 s-1')
      call defvar ('F_evapO', iou, 1, id, -c100, c100, ' '
     &,   'F', 'average evaporation over ocean', ' ','kg m-2 s-1')

# endif
# if defined O_tai_rad
      call defvar ('F_solins', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'incoming solar insolation', ' ', 'W m-2')
      call defvar ('F_upsens', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'surface upward sensible heat', ' ', 'W m-2')
      call defvar ('F_uplwr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'surface net upward longwave', ' ', 'W m-2')
      call defvar ('F_outlwr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'TOA outgoing longwave', ' ', 'W m-2')
      call defvar ('F_dnswr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net surface downward shortwave (abs.)', ' ', 'W m-2')
      call defvar ('F_absswr', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net absorbed shortwave radiation', ' ', 'W m-2')
      call defvar ('F_netrad', iou, 1, id, -c1e3, c1e3, ' '
     &,    'F', 'net top of the atmosphere radiation', ' ', 'W m-2')
      call defvar ('A_albplt', iou, 1, id, c0, c1, ' '
     &,    'F', 'plantary albedo', ' ', '1')
      call defvar ('A_albatm', iou, 1, id, c0, c1, ' '
     &,    'F', 'atmospheric albedo', ' ', '1')
      call defvar ('A_albsur', iou, 1, id, c0, c1, ' '
     &,    'F', 'surface albedo', ' ', '1')
#  if defined O_tai_lo
      call defvar ('A_albsurL', iou, 1, id, c0, c1, ' '
     &,    'F', 'land surface albedo', ' ', '1')
      call defvar ('A_albsurO', iou, 1, id, c0, c1, ' '
     &,    'F', 'ocean surface albedo', ' ', '1')
#  endif
# endif
      call defvar ('O_tempsur', iou, 1, id, -c100, c500, ' '
# if defined O_units_temperature_Celsius
     &,   'F', 'global average sea surface temperature', ' ','C')
# else
     &,   'F', 'global average sea surface temperature', ' ','K')
# endif
      call defvar ('O_salsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface salinity', ' ','psu')
# if defined O_carbon
      call defvar ('O_dicsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface dic', ' ','mol m-3')
#  if defined O_carbon_13
      call defvar ('O_dic13sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface dic13', ' ','mol m-3')
#  endif
#  if defined O_carbon_14
      call defvar ('O_c14sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface c14', ' ','mol m-3')
#  endif
# endif
# if defined O_mobi_silicon
      call defvar ('O_silsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface sil', ' ','mol m-3')
# endif
# if defined O_mobi_alk
      call defvar ('O_alksur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface alkalinity', ' ','mol m-3')
# endif
# if defined O_mobi_o2
      call defvar ('O_o2sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface oxygen', ' ','mol m-3')
# endif
# if defined O_mobi
      call defvar ('O_po4sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface phosphate', ' ','mol m-3')
#  if defined O_mobi_nitrogen
      call defvar ('O_dopsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface DOP', ' ','mol m-3')
      call defvar ('O_no3sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface nitrate', ' ','mol m-3')
      call defvar ('O_donsur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface DON', ' ','mol m-3')
#   if defined O_mobi_nitrogen_15
      call defvar ('O_din15sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface nitrate 15', ' ','mol m-3')
      call defvar ('O_don15sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface DON15', ' ','mol m-3')
#   endif
#  endif
!juan
#  if defined O_mobi_iron
      call defvar ('O_dfesur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface iron', ' ','mol m-3')
#  endif
!
#  if defined O_carbon_13
      call defvar ('O_doc13sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface DOC13', ' ','mol m-3')
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call defvar ('O_cfc11sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface cfc11', ' ','mol m-3')
      call defvar ('O_cfc12sur', iou, 1, id, c0, c100, ' '
     &,   'F', 'global average sea surface cfc12', ' ','mol m-3')
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      call defvar ('A_sulphfor', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average upward shortwave flux due to sulphates'
     &, ' ', 'W m-2')
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
      call defvar ('A_volcfor', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average outward shortwave flux due to volcanos'
     &, ' ', 'W m-2')
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
      call defvar ('A_aggfor', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average longwave reduction due to additional gases'
     &, ' ', 'W m-2')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine embm_tsi_out (fname, avgper, time, stamp, sat, shum
     &,                        precip, evap, v_oice, a_oice, v_snow
     &,                        v_lice, a_lice, co2ccn, co2emit, dc14ccn
     &,                        dc13ccn
     &,                        cfc11ccn, cfc12ccn, scan, nsat, ssat
     &,                        nshum, sshum, nprecip, sprecip, nevap
     &,                        sevap, v_noice, v_soice, a_noice
     &,                        a_soice, v_nsnow, v_ssnow, v_nlice
     &,                        v_slice, a_nlice, a_slice, lsat, osat
     &,                        lprecip, oprecip, levap, oevap, solins
     &,                        upsens, uplwr, outlwr, dnswr, absswr
     &,                        netrad, palb, aalb, salb, lsalb, osalb
     &,                        sst, sss, ssdic, ssdic13, ssc14, ssalk
     &,                        sso2, sspo4, ssdop, ssno3, ssdon, ssdin15
     &,                        ssdon15, ssdoc13, sscfc11, sscfc12, sulph
     &,                        volc, agg, catm, carbemit, ntrec, ssdfe
     &,                        sssil
     &                         )     
!=======================================================================
!     output routine for atmospheric time step integrals

!   inputs:
!     fname    = file name
!     avgper   = length of averaging period
!     time     = time in years
!     stamp    = time stamp
!     sat, ... = data to be written

!   outputs:
!     ntrec    = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real sat, shum, precip, evap, v_oice, a_oice, v_snow, v_lice
      real a_lice, co2ccn, co2emit, scan, dc14ccn, cfc11ccn, cfc12ccn
      real dc13ccn
      real nsat, ssat, nshum, sshum, nprecip, sprecip, nevap, sevap
      real v_noice, v_soice, a_noice, a_soice, v_nsnow, v_ssnow
      real v_nlice, v_slice, a_nlice, a_slice, lsat, osat, lprecip
      real oprecip, levap, oevap, avgper, solins, upsens, uplwr
      real outlwr, dnswr, absswr, netrad, palb, aalb, salb, lsalb
      real osalb, sst, sss, ssdic, ssc14, ssalk, sso2, sspo4, ssno3
      real sscfc11, sscfc12, sulph, volc, agg, catm, carbemit, time, tmp
      real c0, c1, c100, c1e3, c1e4, c1e6, p1, p001, p035, cal2J, C2K
      real kgsPgyr, ssdop, ssdon, ssdin15, ssdon15, ssdic13, ssdoc13
      real ssdfe, sssil

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      C2K = 273.15
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      cal2J = 2.389e-05
      kgsPgyr = 1.e12/(86400.*365.)

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, time, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
# if defined O_units_temperature_Celsius
      call putvars ('A_sat', iou, ntrec, sat, c1, c0)
# else
      call putvars ('A_sat', iou, ntrec, sat, c1, -C2K)
# endif
      call putvars ('A_shum', iou, ntrec, shum, c1, c0)
      call putvars ('F_precip', iou, ntrec, precip, p1, c0)
      call putvars ('F_evap', iou, ntrec, evap, p1, c0)
      call putvars ('A_co2', iou, ntrec, co2ccn, c1, c0)
# if defined O_save_carbon_totals
      call putvars ('A_totcarb', iou, ntrec, catm, c1, c0)
      call putvars ('A_carbemit', iou, ntrec, carbemit, c1, c0)
      call putvars ('F_co2emit', iou, ntrec, co2emit, kgsPgyr, c0)
# else
      call putvars ('F_co2emit', iou, ntrec, co2emit, c1, c0)
# endif
      call putvars ('A_maxiter', iou, ntrec, scan, c1, c0)
# if defined O_carbon_13
      call putvars ('A_dc13', iou, ntrec, dc13ccn, c1, c0)
# endif
# if defined O_carbon_14
      call putvars ('A_dc14', iou, ntrec, dc14ccn, c1, c0)
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call putvars ('A_cfc11', iou, ntrec, cfc11ccn, c1, c0)
      call putvars ('A_cfc12', iou, ntrec, cfc12ccn, c1, c0)
# endif
# if defined O_ice
      call putvars ('O_icevol', iou, ntrec, v_oice, c1e6, c0)
      call putvars ('O_icearea', iou, ntrec, a_oice, c1e4, c0)
      call putvars ('O_snovol', iou, ntrec, v_snow, c1e6, c0)
#  if defined O_landice_data_transient || defined O_ism
      call putvars ('L_icevol', iou, ntrec, v_lice, c1e6, c0)
      call putvars ('L_icearea', iou, ntrec, a_lice, c1e4, c0)
#  endif
# endif
# if defined O_tai_ns
#  if defined O_units_temperature_Celsius
      call putvars ('A_satN', iou, ntrec, nsat, c1, c0)
      call putvars ('A_satS', iou, ntrec, ssat, c1, c0)
#  else
      call putvars ('A_satN', iou, ntrec, nsat, c1, -C2K)
      call putvars ('A_satS', iou, ntrec, ssat, c1, -C2K)
#  endif
      call putvars ('A_shumN', iou, ntrec, nshum, c1, c0)
      call putvars ('A_shumS', iou, ntrec, sshum, c1, c0)
      call putvars ('F_precipN', iou, ntrec, nprecip, p1, c0)
      call putvars ('F_precipS', iou, ntrec, sprecip, p1, c0)
      call putvars ('F_evapN', iou, ntrec, nevap, p1, c0)
      call putvars ('F_evapS', iou, ntrec, sevap, p1, c0)
#  if defined O_ice
      call putvars ('O_icevolN', iou, ntrec, v_noice, c1e6, c0)
      call putvars ('O_icevolS', iou, ntrec, v_soice, c1e6, c0)
      call putvars ('O_iceareaN', iou, ntrec, a_noice, c1e4, c0)
      call putvars ('O_iceareaS', iou, ntrec, a_soice, c1e4, c0)
      call putvars ('O_snovolN', iou, ntrec, v_nsnow, c1e6, c0)
      call putvars ('O_snovolS', iou, ntrec, v_ssnow, c1e6, c0)
#   if defined O_landice_data_transient || defined O_ism
      call putvars ('L_icevolN', iou, ntrec, v_nlice, c1e6, c0)
      call putvars ('L_icevolS', iou, ntrec, v_slice, c1e6, c0)
      call putvars ('L_iceareaN', iou, ntrec, a_nlice, c1e4, c0)
      call putvars ('L_iceareaS', iou, ntrec, a_slice, c1e4, c0)
#   endif
#  endif
# endif
# if defined O_tai_lo
#  if defined O_units_temperature_Celsius
      call putvars ('A_satL', iou, ntrec, lsat, c1, c0)
      call putvars ('A_satO', iou, ntrec, osat, c1, c0)
#  else
      call putvars ('A_satL', iou, ntrec, lsat, c1, -C2K)
      call putvars ('A_satO', iou, ntrec, osat, c1, -C2K)
#  endif
      call putvars ('F_precipL', iou, ntrec, lprecip, p1, c0)
      call putvars ('F_precipO', iou, ntrec, oprecip, p1, c0)
      call putvars ('F_evapL', iou, ntrec, levap, p1, c0)
      call putvars ('F_evapO', iou, ntrec, oevap, p1, c0)
# endif
# if defined O_tai_rad
      call putvars ('F_solins', iou, ntrec, solins, c1e3, c0)
      call putvars ('F_upsens', iou, ntrec, upsens, c1e3, c0)
      call putvars ('F_uplwr', iou, ntrec, uplwr, c1e3, c0)
      call putvars ('F_outlwr', iou, ntrec, outlwr, c1e3, c0)
      call putvars ('F_dnswr', iou, ntrec, dnswr, c1e3, c0)
      call putvars ('F_absswr', iou, ntrec, absswr, c1e3, c0)
      call putvars ('F_netrad', iou, ntrec, netrad, c1e3, c0)
      call putvars ('A_albplt', iou, ntrec, palb, c1, c0)
      call putvars ('A_albatm', iou, ntrec, aalb, c1, c0)
      call putvars ('A_albsur', iou, ntrec, salb, c1, c0)
#  if defined O_tai_lo
      call putvars ('A_albsurL', iou, ntrec, lsalb, c1, c0)
      call putvars ('A_albsurO', iou, ntrec, osalb, c1, c0)
#  endif
# endif
# if defined O_units_temperature_Celsius
      call putvars ('O_tempsur', iou, ntrec, sst, c1, c0)
# else
      call putvars ('O_tempsur', iou, ntrec, sst, c1, -C2K)
# endif
      call putvars ('O_salsur', iou, ntrec, sss, p001, -p035)
# if defined O_carbon
      call putvars ('O_dicsur', iou, ntrec, ssdic, c1, c0)
#  if defined O_carbon_13
      call putvars ('O_dic13sur', iou, ntrec, ssdic13, c1, c0)
#  endif
#  if defined O_carbon_14
      call putvars ('O_c14sur', iou, ntrec, ssc14, c1, c0)
#  endif
# endif
# if defined O_mobi_silicon
      call putvars ('O_silsur', iou, ntrec, sssil, c1, c0)
# endif
# if defined O_mobi_alk
      call putvars ('O_alksur', iou, ntrec, ssalk, c1, c0)
# endif
# if defined O_mobi_o2
      call putvars ('O_o2sur', iou, ntrec, sso2, c1, c0)
# endif
# if defined O_mobi
      call putvars ('O_po4sur', iou, ntrec, sspo4, c1e3, c0)
#  if defined O_mobi_nitrogen
      call putvars ('O_dopsur', iou, ntrec, ssdop, c1e3, c0)
      call putvars ('O_no3sur', iou, ntrec, ssno3, c1e3, c0)
      call putvars ('O_donsur', iou, ntrec, ssdon, c1e3, c0)
#   if defined O_mobi_nitrogen_15
      call putvars ('O_din15sur', iou, ntrec, ssdin15, c1e3, c0)
      call putvars ('O_don15sur', iou, ntrec, ssdon15, c1e3, c0)
#   endif
#  endif
!juan
#  if defined O_mobi_iron
      call putvars ('O_dfesur', iou, ntrec, ssdfe, c1e3, c0)
#  endif
!
#  if defined O_carbon_13
      call putvars ('O_doc13sur', iou, ntrec, ssdoc13, c1, c0)
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call putvars ('O_cfc11sur', iou, ntrec, sscfc11, c1, c0)
      call putvars ('O_cfc12sur', iou, ntrec, sscfc12, c1, c0)
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      call putvars ('A_sulphfor', iou, ntrec, sulph, c1e3, c0)
# endif
# if defined O_volcano_data || defined O_volcano_data_transient
      call putvars ('A_volcfor', iou, ntrec, volc, c1e3, c0)
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient
      call putvars ('A_aggfor', iou, ntrec, agg, c1e3, c0)
# endif
#endif

      return
      end
