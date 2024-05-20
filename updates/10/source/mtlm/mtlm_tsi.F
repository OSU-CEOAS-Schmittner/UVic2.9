      subroutine mtlm_tsi_def (fname, calendar, expnam, runstamp)

#if defined O_mtlm
!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou, ntrec

      real c0, c1, c1e20

      c0 = 0.
      c1 = 1.
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
     &, 'averaging period', ' ','day')
      call defvar ('L_soiltemp', iou, 1, id, -c1e20, c1e20, ' '
# if defined O_units_temperature_Celsius
     &, 'F', 'global average soil temperature', ' ', 'C')
# else
     &, 'F', 'global average soil temperature', ' ', 'K')
# endif
# if defined O_save_carbon_totals
      call defvar ('L_totcarb', iou, 1, id, -c1e20, c1e20, ' ', 'F'
     &, 'global total land carbon', ' ', 'Pg')
      call defvar ('F_carba2l', iou, 1, id, -c1e20, c1e20, ' ', 'F'
     &, 'global total atmopsphere to land carbon flux'
     &, ' ', 'Pg year-1')
      call defvar ('L_soilcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon', ' ', 'Pg')
      call defvar ('L_soilresp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration flux', ' ', 'Pg yr-1')
      call defvar ('L_veglit', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total leaf litter flux', ' ', 'Pg yr-1')
      call defvar ('L_vegburn', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning flux', ' ', 'Pg yr-1')
      call defvar ('L_vegcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon', ' ', 'Pg')
      call defvar ('L_vegnpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity', ' ', 'Pg yr-1')
      call defvar ('L_veggpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total gross primary productivity', ' ', 'Pg yr-1')
#  if defined O_mtlm_carbon_13
      call defvar ('L_soilcarb13', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon 13', ' ', 'Pg')
      call defvar ('L_soilresp13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration 13 flux', ' ', 'Pg yr-1')
      call defvar ('L_vegburn13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning 13 flux', ' ', 'Pg yr-1')
      call defvar ('L_vegcarb13', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon 13', ' ', 'Pg')
      call defvar ('L_vegnpp13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity 13', ' ', 'Pg yr-1')
#  endif
#  if defined O_mtlm_carbon_14
      call defvar ('L_soilcarb14', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon 14', ' ', 'Pg')
      call defvar ('L_soilresp14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration 14 flux', ' ', 'Pg yr-1')
      call defvar ('L_vegburn14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning 14 flux', ' ', 'Pg yr-1')
      call defvar ('L_vegcarb14', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon 14', ' ', 'Pg')
      call defvar ('L_vegnpp14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity 14', ' ', 'Pg yr-1')
#  endif
# else
      call defvar ('L_soilcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon', ' ', 'kg')
      call defvar ('L_soilresp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration flux', ' ', 'kg s-1')
      call defvar ('L_veglit', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total leaf litter flux', ' ', 'kg s-1')
      call defvar ('L_vegburn', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning flux', ' ', 'kg s-1')
      call defvar ('L_vegcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon', ' ', 'kg')
      call defvar ('L_vegnpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity', ' ', 'kg s-1')
      call defvar ('L_veggpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total gross primary productivity', ' ', 'kg s-1')
#  if defined O_mtlm_carbon_13
      call defvar ('L_soilcarb13', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon 13', ' ', 'kg')
      call defvar ('L_soilresp13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration 13 flux', ' ', 'kg s-1')
      call defvar ('L_vegburn13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning 13 flux', ' ', 'kg s-1')
      call defvar ('L_vegcarb13', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon 13', ' ', 'kg')
      call defvar ('L_vegnpp13', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity 13', ' ', 'kg s-1')
#  endif
#  if defined O_mtlm_carbon_14
      call defvar ('L_soilcarb14', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon 14', ' ', 'kg')
      call defvar ('L_soilresp14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration 14 flux', ' ', 'kg s-1')
      call defvar ('L_vegburn14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning 14 flux', ' ', 'kg s-1')
      call defvar ('L_vegcarb14', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon 14', ' ', 'kg')
      call defvar ('L_vegnpp14', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity 14', ' ', 'kg s-1')
#  endif
# endif
      call defvar ('L_veghgt', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average vegetation height', ' ', 'm')
      call defvar ('L_veglai', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average leaf area index', ' ', '1')
# if !defined O_embm
      call defvar ('L_snothk', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total snow mass', ' ', 'kg')
      call defvar ('L_skintemp', iou, 1, id, -c1e20, c1e20, ' '
#  if defined O_units_temperature_Celsius
     &, 'F', 'global average skin temperature', ' ', 'C')
#  else
     &, 'F', 'global average skin temperature', ' ', 'K')
#  endif
      call defvar ('L_soilmois', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average soil moisture', ' ', 'Kg m-2')
      call defvar ('F_evapLtrans', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average evapotranspiration', ' ', 'kg m-2 s-1')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mtlm_tsi_out (fname, avgper, time, stamp, CS, RESP_S
     &,                        LIT_C_T, BURN, CV, NPP, GPP, HT, LAI
     &,                        LYING_SNOW, TSOIL, TSTAR, M, ET, clnd
     &,                        cfa2l, ntrec
# if defined O_mtlm_carbon_13
     &,                        RESP_S13, NPP13, BURN13
     &,                        CS13, CV13
# endif
# if defined O_mtlm_carbon_14
     &,                        RESP_S14, NPP14, BURN14
     &,                        CS14, CV14
# endif
     &                         )
!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname   = file name
!     avgper  = length of averaging period
!     time    = time in years
!     stamp   = time stamp
!     CS, ... = data to be written

!   outputs:
!     ntrec   = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real CS, RESP_S, LIT_C_T, BURN, CV, NPP, GPP, HT, LAI, LYING_SNOW
      real TSOIL, TSTAR, M, ET, clnd, cfa2l, avgper, time, tmp
      real c0, c1, C2K, c1e12, kgsPgyr
# if defined O_mtlm_carbon_13
      real CS13, CV13, RESP_S13, NPP13, BURN13
# endif
# if defined O_mtlm_carbon_14
      real CS14, CV14, RESP_S14, NPP14, BURN14
# endif
      c0 = 0.
      c1 = 1.
      C2K = 273.15
      c1e12 = 1.e12
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
      call putvars ('L_soiltemp', iou, ntrec, TSOIL, c1, C2K)
# else
      call putvars ('L_soiltemp', iou, ntrec, TSOIL, c1, c0)
# endif
# if defined O_save_carbon_totals
      call putvars ('L_totcarb', iou, ntrec, clnd, c1, c0)
      call putvars ('F_carba2l', iou, ntrec, cfa2l, c1, c0)
      call putvars ('L_soilcarb', iou, ntrec, CS, c1e12, c0)
      call putvars ('L_soilresp', iou, ntrec, RESP_S, kgsPgyr, c0)
      call putvars ('L_veglit', iou, ntrec, LIT_C_T, kgsPgyr, c0)
      call putvars ('L_vegburn', iou, ntrec, BURN, kgsPgyr, c0)
      call putvars ('L_vegcarb', iou, ntrec, CV, c1e12, c0)
      call putvars ('L_vegnpp', iou, ntrec, NPP, kgsPgyr, c0)
      call putvars ('L_veggpp', iou, ntrec, GPP, kgsPgyr, c0)
#  if defined O_mtlm_carbon_13
      call putvars ('L_soilcarb13', iou, ntrec, CS13, c1e12, c0)
      call putvars ('L_soilresp13', iou, ntrec, RESP_S13, kgsPgyr, c0)
      call putvars ('L_vegburn13', iou, ntrec, BURN13, kgsPgyr, c0)
      call putvars ('L_vegcarb13', iou, ntrec, CV13, c1e12, c0)
      call putvars ('L_vegnpp13', iou, ntrec, NPP13, kgsPgyr, c0)
#  endif
#  if defined O_mtlm_carbon_14
      call putvars ('L_soilcarb14', iou, ntrec, CS14, c1e12, c0)
      call putvars ('L_soilresp14', iou, ntrec, RESP_S14, kgsPgyr, c0)
      call putvars ('L_vegburn14', iou, ntrec, BURN14, kgsPgyr, c0)
      call putvars ('L_vegcarb14', iou, ntrec, CV14, c1e12, c0)
      call putvars ('L_vegnpp14', iou, ntrec, NPP14, kgsPgyr, c0)
#  endif
# else
      call putvars ('L_soilcarb', iou, ntrec, CS, c1, c0)
      call putvars ('L_soilresp', iou, ntrec, RESP_S, c1, c0)
      call putvars ('L_veglit', iou, ntrec, LIT_C_T, c1, c0)
      call putvars ('L_vegburn', iou, ntrec, BURN, c1, c0)
      call putvars ('L_vegcarb', iou, ntrec, CV, c1, c0)
      call putvars ('L_vegnpp', iou, ntrec, NPP, c1, c0)
      call putvars ('L_veggpp', iou, ntrec, GPP, c1, c0)
#  if defined O_mtlm_carbon_13
      call putvars ('L_soilcarb13', iou, ntrec, CS13, c1, c0)
      call putvars ('L_soilresp13', iou, ntrec, RESP_S13, c1, c0)
      call putvars ('L_vegburn13', iou, ntrec, BURN13, c1, c0)
      call putvars ('L_vegcarb13', iou, ntrec, CV13, c1, c0)
      call putvars ('L_vegnpp13', iou, ntrec, NPP13, c1, c0)
#  endif
#  if defined O_mtlm_carbon_14
      call putvars ('L_soilcarb14', iou, ntrec, CS14, c1, c0)
      call putvars ('L_soilresp14', iou, ntrec, RESP_S14, c1, c0)
      call putvars ('L_vegburn14', iou, ntrec, BURN14, c1, c0)
      call putvars ('L_vegcarb14', iou, ntrec, CV14, c1, c0)
      call putvars ('L_vegnpp14', iou, ntrec, NPP14, c1, c0)
#  endif
# endif
      call putvars ('L_veghgt', iou, ntrec, HT, c1, c0)
      call putvars ('L_veglai', iou, ntrec, LAI, c1, c0)
# if !defined O_embm
      call putvars ('L_snothk', iou, ntrec, LYING_SNOW, c1, c0)
#  if defined O_units_temperature_Celsius
      call putvars ('L_skintemp', iou, ntrec, TSTAR, c1, C2K)
#  else
      call putvars ('L_skintemp', iou, ntrec, TSTAR, c1, c0)
#  endif
      call putvars ('L_soilmois', iou, ntrec, M, c1, c0)
      call putvars ('F_evapLtrans', iou, ntrec, ET, c1, c0)
# endif
#endif

      return
      end
