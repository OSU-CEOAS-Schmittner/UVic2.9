      subroutine mom_tsi_def (fname, calendar, expnam, runstamp)

#if defined O_mom
!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou, ntrec

      real c0, c1, c10, c100, c500, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c10 = 10.
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
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global atributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      id = id_time

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
     &, 'averaging period', ' ','day')
      call defvar ('O_ke', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'kinetic energy per unit volume', ' ', 'J m-3')
      call defvar ('O_temp', iou, 1, id, -c100, c500, ' ', 'F'
#  if defined O_units_temperature_Celsius
     &, 'global average ocean temperature', ' ', 'C')
      call defvar ('O_sal', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean salinity', ' ', '1e-3')
      call defvar ('O_tempvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean temperature', ' ', 'C2')
      call defvar ('O_salvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean salinity', ' ', '1e-6')
      call defvar ('O_absdtemp', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'rate of change of ocean temperature', ' ', 'C s-1')
#  else
     &, 'global average ocean temperature', ' ', 'K')
      call defvar ('O_sal', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean salinity', ' ', '1e-3')
      call defvar ('O_tempvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean temperature', ' ', 'K2')
      call defvar ('O_salvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean salinity', ' ', '1e-6')
      call defvar ('O_absdtemp', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'rate of change of ocean temperature', ' ', 'K s-1')
#  endif
      call defvar ('O_absdsal', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'rate of change of ocean salinity', ' ', '1e-3 s-1')
      call defvar ('O_avgiter', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'ocean solver iterations', ' ', '1')
      call defvar ('F_heat', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean heat flux', ' ','W m-2')
      call defvar ('F_salt', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean salt flux', ' ','kg m-2 s-1')
# if defined O_carbon
      call defvar ('O_dic', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean dic', ' ', 'mol m-3')
      call defvar ('F_dic', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean dic flux', ' ', 'mol m-2 s-1')
#  if defined O_carbon_13
      call defvar ('O_dic13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average dic13', ' ', 'mol m-3')
      call defvar ('F_dic13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average c13 flux', ' ', 'mol m-2 s-1')
#  endif
#  if defined O_carbon_14
      call defvar ('O_c14', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average c14', ' ', 'mol m-3')
      call defvar ('O_dc14', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean delta c14', ' ', 'permil')
      call defvar ('F_c14', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean c14 flux', ' ', 'mol m-2 s-1')
#  endif
# endif
# if defined O_mobi_alk
      call defvar ('O_alk', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean alkalinity', ' ', 'mol m-3')
# endif
# if defined O_mobi_silicon
      call defvar ('O_sil', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean sil', ' ', 'mol m-3')
      call defvar ('F_sil', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean sil flux', ' ', 'mol m-2 s-1')      
      call defvar ('O_opl', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean opal', ' ', 'mol m-3')
# endif
# if defined O_mobi_o2
      call defvar ('O_o2', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean oxygen', ' ', 'mol m-3')
      call defvar ('F_o2', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean oxygen flux', ' ', 'mol m-2 s-1')
# endif
# if defined O_mobi
      call defvar ('O_po4', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean phosphate', ' ', 'mol m-3')
      call defvar ('O_phyt', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean phytoplankton', ' ', 'mol N m-3')
      call defvar ('O_zoop', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean zooplankton', ' ', 'mol N m-3')
      call defvar ('O_detr', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean detritus', ' ', 'mol N m-3')
#  if defined O_mobi_caco3
      call defvar ('O_caco3', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean dead calcite', ' ', 'mol C m-3')
#  endif
#  if defined O_kk_ballast
      call defvar ('O_detr_B', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean ballasted detritus', ' ', 'mol N m-3')
#  endif
#  if defined O_mobi_silicon
      call defvar ('O_diat', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean diatoms', ' ', 'mol N m-3')
#  endif
# endif
# if defined O_mobi_nitrogen
      call defvar ('O_dop', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean DOP', ' ', 'mol m-3')
      call defvar ('O_no3', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean nitrate', ' ', 'mol m-3')
      call defvar ('O_don', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean DON', ' ', 'mol m-3')
      call defvar ('O_diaz', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean diazotrophs', ' ', 'mol N m-3')
#  if defined O_mobi_nitrogen_15
      call defvar ('O_din15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean nitrate 15', ' ', 'mol m-3')
      call defvar ('O_don15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean DON 15', ' ', 'mol m-3')
      call defvar ('O_phytn15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean phytoplankton n15', ' ', 'mol N m-3')
#   if defined O_mobi_silicon
      call defvar ('O_diatn15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean diat n15', ' ', 'mol N m-3')
#   endif      
      call defvar ('O_zoopn15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean zooplankton n15', ' ', 'mol N m-3')
      call defvar ('O_detrn15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean detritus n15', ' ', 'mol N m-3')
      call defvar ('O_diazn15', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average ocean diazotroph n15', ' ', 'mol N m-3')
#  endif
# endif
!juan
# if defined O_mobi_iron
      call defvar ('O_dfe', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean iron', ' ', 'mol Fe m-3')
      call defvar ('O_detrfe', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean particulate iron', ' ', 'mol Fe m-3')
# endif  
!
# if defined O_carbon_13
      call defvar ('O_doc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average DOC13', ' ', 'mol m-3')
      call defvar ('O_phytc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average phytoplankton c13', ' ', 'mol m-3')
#  if defined O_mobi_silicon
      call defvar ('O_diatc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average diatom c13', ' ', 'mol m-3')
#  endif      
#  if defined O_mobi_caco3
      call defvar ('O_caco3c13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average CaCO3 c13', ' ', 'mol m-3')
#  endif      
      call defvar ('O_zoopc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average zooplankton c13', ' ', 'mol m-3')
      call defvar ('O_detrc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average detrius c13', ' ', 'mol m-3')
#  if defined O_mobi_nitrogen
      call defvar ('O_diazc13', iou, 1, id, c0, c100, ' ', 'F'
     &, 'Global average diazotroph c13', ' ', 'mol m-3')
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call defvar ('O_cfc11', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean CFC11', ' ', 'mol m-3')
      call defvar ('F_cfc11', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean CFC11 flux', ' ', 'mol m-2 s-1')
      call defvar ('O_cfc12', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean CFC12', ' ', 'mol m-3')
      call defvar ('F_cfc12', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average ocean CFC12 flux', ' ', 'mol m-2 s-1')
# endif
# if defined O_tai_otsf
      call defvar ('O_motmax', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'maximum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
      call defvar ('O_motmin', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'minimum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
# endif
# if defined O_tai_slh
      call defvar ('O_dsealev', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'relative sea level height', ' ', 'm')
# endif
# if defined O_save_carbon_carbonate_chem
      call defvar ('O_phsur', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average sea surface pH', ' ', '1')
      call defvar ('O_co3sur', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average sea surface CO3', ' ', 'mol m-3')
      call defvar ('O_ocalcsur', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average sea surface omega calcite', ' ', '1')
      call defvar ('O_oaragsur', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average sea surface omega aragonite', ' ', '1')
      call defvar ('O_pco2sur', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average sea surface pCO2', ' ', 'ppmv')
# endif
# if defined O_save_carbon_totals
      call defvar ('O_totcarb', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average total ocean carbon', ' ', 'Pg')
      call defvar ('F_carba2o', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average total atmosphere to ocean carbon flux'
     &, ' ', 'Pg year-1')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mom_tsi_out (fname, avgper, time, stamp, ektot, tbar1
     &,                       tbar2, travar1, travar2, dtabs1, dtabs2
     &,                       scan, hflx, sflx, tbar_dic, dicflx
     &,                       tbar_dic13, dic13flx
     &,                       tbar_alk, tbar_o2, o2flx, tbar_po4
     &,                       tbar_dop, tbar_phyt, tbar_zoop, tbar_detr
     &,                       tbar_no3, tbar_don, tbar_diaz, tbar_din15
     &,                       tbar_don15, tbar_phytn15, tbar_zoopn15
     &,                       tbar_detrn15, tbar_diazn15, tbar_doc13
     &,                       tbar_phytc13, tbar_zoopc13, tbar_detrc13
     &,                       tbar_diazc13, tbar_c14
     &,                       tbar_dc14, c14flx, tbar_cfc11, cfc11flx
     &,                       tbar_cfc12, cfc12flx, otmax, otmin, slh
     &,                       sspH, ssCO3, ssOc, ssOa, sspCO2, cocn
     &,                       cfa2o, ntrec, tbar_dfe, tbar_ddfe
# if defined O_kk_ballast
     &,                       tbar_d_B
# endif
# if defined O_mobi_caco3
     &,                       tbar_caco3, tbar_caco3c13
# endif
# if defined O_mobi_silicon
     &,                       tbar_diat
     &,                       tbar_diatn15, tbar_diatc13
     &,                       tbar_sil, silflx, tbar_opl
# endif
     &     )
!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     avgper     = length of averaging period
!     time       = time in years
!     stamp      = time stamp
!     ektot, ... = data to be written

!   outputs:
!     ntrec      = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real ektot, tbar1, tbar2, travar1, travar2, dtabs1, dtabs2, scan
      real hflx, sflx, tbar_dic, dicflx, tbar_alk, tbar_o2, o2flx
      real tbar_po4, tbar_phyt, tbar_zoop, tbar_detr, tbar_no3
      real tbar_diaz, tbar_c14, tbar_dc14, c14flx, tbar_cfc11
      real cfc11flx, tbar_cfc12, cfc12flx
      real otmax, otmin, slh, sspH, ssCO3, ssOc, ssOa, sspCO2, cocn
      real cfa2o, avgper, time, tmp, c0, c1, c10, c100, c1e3, c1e6
      real p001, p035, p1, C2K, cal2J, tbar_dop, tbar_don, tbar_din15
      real tbar_don15, tbar_phytn15, tbar_zoopn15, tbar_detrn15
      real tbar_diazn15, tbar_diatn15
      real tbar_dic13, dic13flx, tbar_doc13, tbar_phytc13, tbar_zoopc13
      real tbar_detrc13, tbar_diazc13, tbar_dfe, tbar_ddfe
      real tbar_d_B, tbar_c, tbar_caco3, tbar_diatc13, tbar_caco3c13
      real tbar_diat, tbar_sil, silflx, tbar_opl
      
      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e6 = 1.e6
      p001 = 0.001
      p035 = 0.035
      p1 = 0.1
      C2K = 273.15
      cal2J = 2.389e-05

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
      tmp = nyear
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
      call putvars ('O_ke', iou, ntrec, ektot, c10, c0)
#  if defined O_units_temperature_Celsius
      call putvars ('O_temp', iou, ntrec, tbar1, c1, c0)
      call putvars ('O_sal', iou, ntrec, tbar2, p001, -p035)
      call putvars ('O_tempvar', iou, ntrec, travar1, c1, c0)
#  else
      call putvars ('O_temp', iou, ntrec, tbar1, c1, -C2K)
      call putvars ('O_sal', iou, ntrec, tbar2, p001, -p035)
      tmp = 0.0
      if (travar1 .gt. 0.0) tmp = (sqrt(travar1) + C2K)**2
      call putvars ('O_tempvar', iou, ntrec, tmp, c1, c0)
#  endif
      call putvars ('O_salvar', iou, ntrec, tmp, c1, c0)
      call putvars ('O_absdtemp', iou, ntrec, dtabs1, c1, c0)
      call putvars ('O_absdsal', iou, ntrec, dtabs2, p001, c0)
      call putvars ('O_avgiter', iou, ntrec, scan, c1, c0)
      call putvars ('F_heat', iou, ntrec, hflx, cal2J, c0)
      call putvars ('F_salt', iou, ntrec, sflx, p1, c0)
# if defined O_carbon
      call putvars ('O_dic', iou, ntrec, tbar_dic, c1, c0)
      call putvars ('F_dic', iou, ntrec, dicflx, c100, c0)
#  if defined O_carbon_13
      call putvars ('O_dic13', iou, ntrec, tbar_dic13, c1, c0)
      call putvars ('F_dic13', iou, ntrec, dic13flx, c100, c0)
#  endif
#  if defined O_carbon_14
      call putvars ('O_c14', iou, ntrec, tbar_c14, c1, c0)
      call putvars ('O_dc14', iou, ntrec, tbar_dc14, c1, c0)
      call putvars ('F_c14', iou, ntrec, c14flx, c100, c0)
#  endif
# endif
# if defined O_mobi_alk
      call putvars ('O_alk', iou, ntrec, tbar_alk, c1, c0)
# endif
# if defined O_mobi_silicon
      call putvars ('O_sil', iou, ntrec, tbar_sil, c1, c0)
      call putvars ('F_sil', iou, ntrec, silflx, c100, c0)     
      call putvars ('O_opl', iou, ntrec, tbar_opl, c1, c0)
# endif
# if defined O_mobi_o2
      call putvars ('O_o2', iou, ntrec, tbar_o2, c1, c0)
      call putvars ('F_o2', iou, ntrec, o2flx, c100, c0)
# endif
# if defined O_mobi
      call putvars ('O_po4', iou, ntrec, tbar_po4, c1e3, c0)
      call putvars ('O_phyt', iou, ntrec, tbar_phyt, c1e3, c0)
      call putvars ('O_zoop', iou, ntrec, tbar_zoop, c1e3, c0)
      call putvars ('O_detr', iou, ntrec, tbar_detr, c1e3, c0)
#  if defined O_mobi_caco3
      call putvars ('O_caco3', iou, ntrec, tbar_caco3, c1e3, c0)
#  endif
#  if defined O_kk_ballast
      call putvars ('O_detr_B', iou, ntrec, tbar_d_B, c1e3, c0)
#  endif
#  if defined O_mobi_silicon
      call putvars ('O_diat', iou, ntrec, tbar_diat, c1e3, c0)
#  endif
# endif
# if defined O_mobi_nitrogen
      call putvars ('O_dop', iou, ntrec, tbar_dop, c1e3, c0)
      call putvars ('O_no3', iou, ntrec, tbar_no3, c1e3, c0)
      call putvars ('O_don', iou, ntrec, tbar_don, c1e3, c0)
      call putvars ('O_diaz', iou, ntrec, tbar_diaz, c1e3, c0)
#  if defined O_mobi_nitrogen_15
      call putvars ('O_din15', iou, ntrec, tbar_din15, c1e3, c0)
      call putvars ('O_don15', iou, ntrec, tbar_don15, c1e3, c0)
      call putvars ('O_phytn15', iou, ntrec, tbar_phytn15, c1e3, c0)
#   if defined O_mobi_silicon
      call putvars ('O_diatn15', iou, ntrec, tbar_diatn15, c1e3, c0)
#   endif      
      call putvars ('O_zoopn15', iou, ntrec, tbar_zoopn15, c1e3, c0)
      call putvars ('O_detrn15', iou, ntrec, tbar_detrn15, c1e3, c0)
      call putvars ('O_diazn15', iou, ntrec, tbar_diazn15, c1e3, c0)
#  endif
# endif
!juan
# if defined O_mobi_iron
      call putvars ('O_dfe', iou, ntrec, tbar_dfe, c1e3, c0)
      call putvars ('O_detrfe', iou, ntrec, tbar_ddfe, c1e3, c0)
# endif
!
# if defined O_carbon_13
      call putvars ('O_doc13', iou, ntrec, tbar_doc13, c1, c0)
      call putvars ('O_phytc13', iou, ntrec, tbar_phytc13, c1, c0)
#  if defined O_mobi_silicon
      call putvars ('O_diatc13', iou, ntrec, tbar_diatc13, c1, c0)
#  endif      
#  if defined O_mobi_caco3
      call putvars ('O_caco3c13', iou, ntrec, tbar_caco3c13, c1, c0)
#  endif      
      call putvars ('O_zoopc13', iou, ntrec, tbar_zoopc13, c1, c0)
      call putvars ('O_detrc13', iou, ntrec, tbar_detrc13, c1, c0)
#  if defined O_mobi_nitrogen
      call putvars ('O_diazc13', iou, ntrec, tbar_diazc13, c1, c0)
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call putvars ('O_cfc11', iou, ntrec, tbar_cfc11, c1, c0)
      call putvars ('F_cfc11', iou, ntrec, cfc11flx, c100, c0)
      call putvars ('O_cfc12', iou, ntrec, tbar_cfc12, c1, c0)
      call putvars ('F_cfc12', iou, ntrec, cfc12flx, c100, c0)
# endif
# if defined O_tai_otsf
      call putvars ('O_motmax', iou, ntrec, otmax, c1e6, c0)
      call putvars ('O_motmin', iou, ntrec, otmin, c1e6, c0)
# endif
# if defined O_tai_slh
      call putvars ('O_dsealev', iou, ntrec, slh, c100, c0)
# endif
# if defined O_save_carbon_carbonate_chem
      call putvars ('O_phsur', iou, ntrec, sspH, c1, c0)
      call putvars ('O_co3sur', iou, ntrec, ssCO3, c1, c0)
      call putvars ('O_ocalcsur', iou, ntrec, ssOc, c1, c0)
      call putvars ('O_oaragsur', iou, ntrec, ssOa, c1, c0)
      call putvars ('O_pco2sur', iou, ntrec, sspCO2, c1, c0)
# endif
# if defined O_save_carbon_totals
      call putvars ('O_totcarb', iou, ntrec, cocn, c1, c0)
      call putvars ('F_carba2o', iou, ntrec, cfa2o, c1, c0)
# endif
#endif

      return
      end
