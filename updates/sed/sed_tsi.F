      subroutine sed_tsi_def (fname, calendar, expnam, runstamp)

#if defined O_sed
!=======================================================================
!     output routine for sediment time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou

      real c0, c1, c100, c400, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c400 = 400.
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
      call defvar ('F_caco3dis', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average sediment upward flux of calcite', ' '
     &, 'mol m-2 s-1')
      call defvar ('F_caco3rai', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average sediment downward flux of calcite', ' '
     &, 'mol m-2 s-1')
      call defvar ('S_caco3per', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average CaCO3 pore layer portion', ' ', 'percent')
# if defined O_save_carbon_totals
      call defvar ('S_carblith', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average change in lithosphere carbon', ' ', 'Pg')
      call defvar ('F_carbo2s', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average total ocean to sediment carbon flux'
     &, ' ', 'Pg year-1')
      call defvar ('F_carbl2o', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average total land to ocean carbon flux (weathering)'
     &, ' ', 'Pg year-1')
      call defvar ('S_caco3mas', iou, 1, id, -c1e3, c1e20, ' ', 'F'
     &, 'total CaCO3 pore layer mass', ' ', 'Pg C')
      call defvar ('S_caco3bur', iou, 1, id, -c1e3, c1e20, ' ', 'F'
     &, 'total CaCO3 buried mass', ' ', 'Pg C')
# else
      call defvar ('S_caco3mas', iou, 1, id, -c1e3, c1e20, ' ', 'F'
     &, 'total CaCO3 pore layer mass', ' ', 'kg C')
      call defvar ('S_caco3bur', iou, 1, id, -c1e3, c1e20, ' ', 'F'
     &, 'total CaCO3 buried mass', ' ', 'kg C')
# endif
      call defvar ('S_co3', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average sediment surface CO3 concentration', ' '
     &, 'mol m-3')
      call defvar ('S_co3sat', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global average sediment surface CO3 saturation ', ' '
     &, 'mol m-3')
      call defvar ('F_weath', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'global total weathering flux', ' ', 'kg C s-1')
      call defvar ('S_rainrat', iou, 1, id, -c1e3, c1e3, ' ', 'F'
     &, 'rain ratio', ' ', '1')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine sed_tsi_out (fname, avgper, time, stamp, ttrcal
     &,                       rain_cal, cal, calmass, calmass_bur, co3
     &,                       co3sat, weathflx, rainr, carblith, cfo2s
     &,                       cfl2o, ntrec)
!=======================================================================
!     output routine for sediment time step integrals

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

      real ttrcal, rain_cal, cal, calmass, calmass_bur, co3, co3sat
      real weathflx, rainr, carblith, cfo2s, cfl2o, avgper, time, tmp
      real c0, c1, c1e3, c1e4, c12e6, c1e12, p01, p001, p0001

      c0 = 0.
      c1 = 1.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c12e6 = 12.e6
      c1e12 = 1.e12
      p01 = 0.01
      p001 = 0.001
      p0001 = 0.0001

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
      call putvars ('F_caco3dis', iou, ntrec, ttrcal, p0001, c0)
      call putvars ('F_caco3rai', iou, ntrec, rain_cal, p0001, c0)
      call putvars ('S_caco3per', iou, ntrec, cal, p01, c0)
# if defined O_save_carbon_totals
      call putvars ('S_carblith', iou, ntrec, carblith, c1, c0)
      call putvars ('F_carbo2s', iou, ntrec, cfo2s, c1, c0)
      call putvars ('F_carbl2o', iou, ntrec, cfl2o, c1, c0)
      call putvars ('S_caco3mas', iou, ntrec, calmass, c1e12, c0)
      call putvars ('S_caco3bur', iou, ntrec, calmass_bur, c1e12, c0)
# else
      call putvars ('S_caco3mas', iou, ntrec, calmass, c1, c0)
      call putvars ('S_caco3bur', iou, ntrec, calmass_bur, c1, c0)
# endif
      call putvars ('S_co3', iou, ntrec, co3, p001, c0)
      call putvars ('S_co3sat', iou, ntrec, co3sat, p001, c0)
      call putvars ('F_weath', iou, ntrec, weathflx, c1, c0)
      call putvars ('S_rainrat', iou, ntrec, rainr, c1, c0)
#endif

      return
      end
