      subroutine mom_tavg_def (fname, imt, jmt, km, nt, kpzd, xt, yt
     &,                        calendar, expnam, runstamp, mapt)

#if defined O_mom
!=======================================================================
!     definition routine for ocean time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapt         = tracer map
!=======================================================================

      implicit none

      integer iou, j, n, imt, jmt, km, nt, kpzd, igs, ige, ig, jgs, jge
      integer jg, kgs, kge, kg, lgs, lge, lg, it(10), iu(10), id_time
      integer id_xt, id_xu, id_yt, id_yu, id_zt, id_zw, id_zl, id_xt_e
      integer id_xu_e, id_yt_e, id_yu_e, id_zt_e, id_zw_e, id_zl_e

      save iou

      character(*) :: fname, calendar, expnam, runstamp
      character(3) :: a3
      character(10) :: mapt(nt)

      real xt(imt), yt(jmt)
      real c0, c1, c100, c500, c1e3, c1e4, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
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
      kgs = 1
      kge = km
      kg  = kge-kgs+1
# if defined O_save_npzd
      lgs = max(1,kgs)
      lge = min(kpzd,kg)
      lg  = lge-lgs+1
# endif

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
      call defdim ('longitude', iou, ig, id_xt)
      call defdim ('latitude', iou, jg, id_yt)
      call defdim ('longitude_V', iou, ig, id_xu)
      call defdim ('latitude_V', iou, jg, id_yu)
      call defdim ('depth', iou, kg, id_zt)
      call defdim ('depth_W', iou, kg, id_zw)
      call defdim ('longitude_edges', iou, ig+1, id_xt_e)
      call defdim ('latitude_edges', iou, jg+1, id_yt_e)
      call defdim ('longitude_V_edges', iou, ig+1, id_xu_e)
      call defdim ('latitude_V_edges', iou, jg+1, id_yu_e)
      call defdim ('depth_edges', iou, kg+1, id_zt_e)
      call defdim ('depth_W_edges', iou, kg+1, id_zw_e)
# if defined O_save_npzd
      call defdim ('npzdlev', iou, lg, id_zl)
      call defdim ('npzdlev_edges', iou, lg+1, id_zl_e)
# endif

!-----------------------------------------------------------------------
!     define 1d data (t)
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
      it(1) = id_zt
      call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      call defvar ('G_dzt', iou, 1, it, c0, c0, ' ', 'D'
     &, 'thickness t grid', ' ', 'm')
      it(1) = id_zw
      call defvar ('depth_W', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      it(1) = id_xt_e
      call defvar ('longitude_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yt_e
      call defvar ('latitude_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'latitude', 'degrees_north')
      it(1) = id_xu_e
      call defvar ('longitude_V_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude edges', 'longitude', 'degrees_east')
      it(1) = id_yu_e
      call defvar ('latitude_V_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude edges', 'latitude', 'degrees_north')
        it(1) = id_zt_e
      call defvar ('depth_edges', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth edges', 'depth', 'm')
      it(1) = id_zw_e
      call defvar ('depth_W_edges', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth of w grid edges', 'depth', 'm')
# if defined O_save_npzd
      it(1) = id_zl
      call defvar ('npzdlev', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth', 'depth', 'm')
      it(1) = id_zl_e
      call defvar ('npzdlev_edges', iou, 1, it, c0, c0, 'Z', 'D'
     &, 'depth edges', ' ', 'm')
# endif

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('G_kmt', iou, 2, it, c0, c1e6, ' ', 'I'
     &, 'ocean grid depth level', 'model_level_number' ,'1')
      call defvar ('G_mskhr', iou, 2, it, c0, c1e6, ' ', 'I'
     &, 'horizontal region mask', ' ' ,'1')
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degrees_east')
      call defvar ('G_latU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid latitude', 'latitude', 'degrees_north')
      call defvar ('G_lonU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid longitude', 'longitude', 'degrees_east')
# if defined O_save_anisotropic_viscosity
      call defvar ('O_mvisc', iou, 2, iu, c0, c1e20,' ', 'F'
     &, 'meridional viscosity', 'viscosity', 'cm2 s-1')
      call defvar ('O_zvisc', iou, 2, iu, c0, c1e20, ' ', 'F'
     &, 'zonal viscosity', 'viscosity', 'cm2 s-1')
# endif
      call defvar ('G_areaT', iou, 2, it, -c1e6, c1e6, ' ', 'D'
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
      do n=1,nt
        if (trim(mapt(n)) .eq. 'temp') then
          call defvar ('F_heat', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &,     'surface downward heat flux', ' ', 'W m-2')
        elseif (trim(mapt(n)) .eq. 'salt') then
          call defvar ('F_salt', iou,3, it, -c100, c100, ' ', 'F'
     &,     'surface downward salt flux', ' ', 'kg m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'dic') then
          call defvar ('F_dic', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward carbon flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'alk') then
          call defvar ('F_alk', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward alkalinity flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'o2') then
          call defvar ('F_o2', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward oxygen flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'po4') then
          call defvar ('F_po4', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward phosphate flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'no3') then
          call defvar ('F_no3', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward nitrate flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'c14') then
          call defvar ('F_c14', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward carbon 14 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call defvar ('F_cfc11', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC11 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call defvar ('F_cfc12', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC12 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'phyt') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        else
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call defvar ('F_'//trim(a3), iou ,3, it, -c1e6, c1e6, ' '
     &,     'F', 'tracer flux '//trim(a3)
     &,     'tracer_flux_'//trim(a3), 'unknown')
        endif
      enddo
      call defvar ('O_tauX', iou, 3, iu, -c1e6, c1e6, ' ', 'F'
     &, 'surface eastward momentum flux'
     &, 'surface_downward_eastward_stress', 'Pa')
      call defvar ('O_tauY', iou, 3, iu, -c1e6, c1e6, ' ', 'F'
     &, 'surface northward momentum flux'
     &, 'surface_downward_northward_stress', 'Pa')
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      call defvar ('O_presssur', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'surface pressure', ' ', 'Pa')
# endif
# if defined O_stream_function
      call defvar ('O_psi', iou, 3, it, -c1e20, c1e20, ' ', 'F'
     &, 'transport streamfunction'
     &, 'ocean_barotropic_streamfunction', 'm3 s-1')
# endif
# if defined O_save_convection
      call defvar ('O_convlev', iou, 3, it, -c1, c1e6, ' ', 'F'
     &, 'number of convected levels', ' ', '1')
      call defvar ('O_ventdep', iou, 3, it, -c1, c1e6, ' ', 'F'
     &, 'ventelation depth', ' ', 'm')
      call defvar ('O_convpe', iou, 3, it, -c1, c1e6, ' ', 'F'
     &, 'potential energy lost due to convection', ' ', 'W m-2')
# endif
# if defined O_save_carbon_carbonate_chem
      call defvar ('O_phsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface pH', ' ', '1')
      call defvar ('O_co3sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface CO3', ' ', 'mol m-3')
      call defvar ('O_ocalcsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega calcite', ' ', '1')
      call defvar ('O_oaragsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega aragonite', ' ', '1')
      call defvar ('O_pco2sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface pCO2', ' ', '1')
# endif
# if defined O_save_npzd
      call defvar ('O_caco3pro', iou, 3, it, c0, c100, ' ', 'F'
     &, 'calcite production', ' ', 'mol C m-2 s-1')
# endif

!-----------------------------------------------------------------------
!     define time dependent 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_zt
      iu(3) = id_zt
      it(4) = id_time
      iu(4) = id_time
      do n=1,nt
        if (trim(mapt(n)) .eq. 'temp') then
          call defvar ('O_temp', iou, 4, it, -c100, c500, ' ', 'F'
     &,     'ocean potential temperature'
#  if defined O_units_temperature_Celsius
     &,     'sea_water_potential_temperature', 'C')
#  else
     &,     'sea_water_potential_temperature', 'K')
#  endif
        elseif (trim(mapt(n)) .eq. 'salt') then
          call defvar ('O_sal', iou, 4, it, c0, c100, ' ', 'F'
     &,     'ocean salinity', 'sea_water_salinity', '1e-3')
        elseif (trim(mapt(n)) .eq. 'dic') then
          call defvar ('O_dic', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'dissolved inorganic carbon', ' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'alk') then
          call defvar ('O_alk', iou, 4, it, -c1, c500, ' ', 'F'
     &,     'total alkalinity',' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'o2') then
          call defvar ('O_o2', iou, 4, it, -c1, c500, ' ', 'F'
     &,     'O_o2',' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'po4') then
          call defvar ('O_po4', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'O_po4', ' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'phyt') then
          call defvar ('O_phyt', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'phytoplankton', ' ', 'mol N m-3')
        elseif (trim(mapt(n)) .eq. 'zoop') then
          call defvar ('O_zoop', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'zooplankton', ' ', 'mol N m-3')
        elseif (trim(mapt(n)) .eq. 'detr') then
          call defvar ('O_detr', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'detritus', ' ', 'mol N m-3')
        elseif (trim(mapt(n)) .eq. 'no3') then
          call defvar ('O_no3', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'O_no3', ' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'diaz') then
          call defvar ('O_diaz', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'diazotrophs', ' ', 'mol N m-3')
        elseif (trim(mapt(n)) .eq. 'c14') then
          call defvar ('O_c14', iou, 4, it, -c1e3, c1e3, ' ', 'F'
     &,     'carbon 14', ' ', 'mol m-3')
          call defvar ('O_dc14', iou, 4, it, -c1e3, c1e3, ' ', 'F'
     &,     'delta carbon 14', ' ', 'permil')
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call defvar ('O_cfc11', iou, 4, it, -c1, c500, ' ', 'F'
     &,     'O_cfc11', ' ', 'mol m-3')
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call defvar ('O_cfc12', iou, 4, it, -c1, c500, ' ', 'F'
     &,     'O_cfc12', ' ', 'mol m-3')
        else
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call defvar ('tracer_'//trim(a3), iou ,4, it, -c1e6, c1e6
     &,     ' ', 'F', 'tracer '//trim(a3)
     &,     'tracer_'//trim(a3), 'unknown')
        endif
      enddo
# if defined O_save_carbon_totals
       call defvar ('O_totcarb', iou, 4, it, -c1, c100, ' ', 'F'
     &,     'total carbon', ' ', 'mol m-3')
# endif
      iu(3) = id_zt
      call defvar ('O_velX', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'eastward ocean velocity', 'eastward_sea_water_velocity'
     &, 'm s-1')
      call defvar ('O_velY', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'northward ocean velocity', 'northward_sea_water_velocity'
     &, 'm s-1')
        iu(3) = id_zw
      call defvar ('O_velZ', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'upward ocean velocity', 'upward_sea_water_velocity'
     &, 'm s-1')
# if defined O_gent_mcwilliams
        iu(1) = id_xu
        iu(2) = id_yt
        iu(3) = id_zt
      call defvar ('O_gmvelX', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'eastward GM velocity', ' ', 'm s-1')
        iu(1) = id_xt
        iu(2) = id_yu
        iu(3) = id_zt
      call defvar ('O_gmvelY', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'northward GM velocity', ' ', 'm s-1')
        iu(1) = id_xt
        iu(2) = id_yt
        iu(3) = id_zw
      call defvar ('O_gmvelZ', iou, 4, iu, -c100, c100, ' ', 'F'
     &, 'upward GM velocity', ' ', 'm s-1')
# endif
# if defined O_save_kv
      iu(1) = id_xu
      iu(2) = id_yu
      iu(3) = id_zw
      call defvar ('O_diffZ', iou, 4, iu, c0, c1, ' ', 'F'
     &, 'diapycnal diffusivity', ' ', 'm2 s-1')
# endif
# if defined O_save_npzd
      it(3) = id_zl
# else
      it(3) = id_zt
# endif
# if defined O_save_npzd
      call defvar ('O_phytnpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_phytgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on P', ' ', 'mol N m-3 s-1')
      call defvar ('O_zoograz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on Z', ' ', 'mol N m-3 s-1')
      call defvar ('O_detgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on detritus', ' ', 'mol N m-3 s-1')
      call defvar ('O_phytmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'phytoplankton mortality rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_phytrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'microbial fast recycling rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_zoopmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'zooplankton mortality rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_excret', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'excretion rate', ' ', 'mol N m-3 s-1')
# if defined O_npzd_extra_diagnostics
      call defvar ('O_avej', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant phyt. growth rate', ' ', 's-1')
      call defvar ('O_avej_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant Diaz growth rate', ' ', 's-1')
      call defvar ('O_gmax', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'temp-dependant zoo growth rate', ' ', 's-1')
      call defvar ('O_no3P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'no3-dependant phyt. growth rate', ' ', 's-1')
      call defvar ('O_po4P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant phyt. growth rate', ' ', 's-1')
      call defvar ('O_po4_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant Diaz growth rate', ' ', 's-1')
# endif
#  if defined O_npzd_nitrogen
      call defvar ('O_diaznpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'net primary production diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_diazgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_diazmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'mortality diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_no3fix', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'nitrogen fixation', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_denitr', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'denitrification', ' ', 'mol N m-3 s-1')
#  endif
      it(3) = id_zt
      call defvar ('O_detrremi', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'remineralisation rate', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_detrexp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'detrital export rate', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_caco3exp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'calcite export rate', ' ', 'mol C m-3 s-1')
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mom_tavg_out (fname, ids, ide, jds, jde, imt, jmt, km
     &,                        nt, xt, yt, zt, xu, yu, zw, dxt, dyt, dzt
     &,                        dxu, dyu, dzw, avgper, time, stamp, mapt
     &,                        t, u, v, adv_vbt, stf, taux, tauy, redctn
# if defined O_gent_mcwilliams
     &,                        adv_vetiso, adv_vntiso, adv_vbtiso
# endif
# if defined O_carbon && defined O_carbon_14
     &,                        dc14
# endif
# if defined O_save_kv
     &,                        diff_cbt
# endif
# if defined O_save_npzd
     &,                        kpzd, npp, graz, morp, morpt, morz, remi
     &,                        excr, expo, expocal, procal, graz_Z
     &,                        graz_Det
#  if defined O_npzd_nitrogen
     &,                        npp_D, graz_D, morp_D, nfix, deni
#  endif
# if defined O_npzd_extra_diagnostics
     &,                        avej, avej_D, gmax, no3P, po4P, po4_D
# endif
# endif
# if defined O_save_convection
     &,                        totalk, vdepth, pe
# endif
#if defined O_save_carbon_carbonate_chem
     &,                        sspH, ssCO3, ssOc, ssOa, sspCO2
#endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
     &,                        ps
# endif
# if defined O_stream_function
     &,                        psi
# endif
     &,                        kmt, mskhr, tm, um
     &,                        tlat, tlon, ulat, ulon
# if defined O_save_anisotropic_viscosity
     &,                        mvisc, zvisc
# endif  
     &,                        tgarea, ugarea
     &,                        ntrec)
!=======================================================================
!     output routine for ocean time averages

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if  no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths
!     avgper       = length of averaging period
!     time         = time in years
!     t, ...       = data to be written

!   outputs:
!     ntrec        = number of time record in file
!=======================================================================

      implicit none

      integer iou, j, ln, n, ntrec, imt, jmt, km, nt, kpzd, ids, ide
      integer jds, jde, igs, ige, ig, jgs, jge, jg, kgs, kge, kg, lgs
      integer lge, lg, ils, ile, jls, jle, kls, kle, lls, lle, ib(10)
      integer ic(10), kmt(ids:ide,jds:jde), mskhr(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      save iou

      character(*) :: fname, stamp
      character(3) :: a3
      character(10) :: mapt(nt)

      real xt(imt), xu(imt), yt(jmt), yu(jmt), zt(km), zw(km), avgper
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt), dzt(km), dzw(km)
      real t(ids:ide,jds:jde,km,nt), u(ids:ide,jds:jde,km)
      real v(ids:ide,jds:jde,km), adv_vbt(ids:ide,jds:jde,km)
      real stf(ids:ide,jds:jde,nt), taux(ids:ide,jds:jde)
      real tauy(ids:ide,jds:jde), redctn
# if defined O_gent_mcwilliams
      real adv_vetiso(ids:ide,jds:jde,km)
      real adv_vntiso(ids:ide,jds:jde,km)
      real adv_vbtiso(ids:ide,jds:jde,km)
# endif
# if defined O_save_carbon_totals
      real toc(ids:ide,jds:jde,km)
# endif
# if defined O_carbon && defined O_carbon_14
      real dc14(ids:ide,jds:jde,km)
# endif
# if defined O_save_kv
      real diff_cbt(ids:ide,jds:jde,km)
# endif
# if defined O_save_npzd
      real npp(ids:ide,jds:jde,kpzd), graz(ids:ide,jds:jde,kpzd)
      real graz_Z(ids:ide,jds:jde,kpzd) 
      real graz_Det(ids:ide,jds:jde,kpzd)
      real morp(ids:ide,jds:jde,kpzd), morz(ids:ide,jds:jde,kpzd)
      real excr(ids:ide,jds:jde,kpzd), morpt(ids:ide,jds:jde,kpzd)
      real remi(ids:ide,jds:jde,km), expo(ids:ide,jds:jde,km)
      real expocal(ids:ide,jds:jde,km), procal(ids:ide,jds:jde)
#  if defined O_npzd_nitrogen
      real npp_D(ids:ide,jds:jde,kpzd), graz_D(ids:ide,jds:jde,kpzd)
      real morp_D(ids:ide,jds:jde,kpzd), nfix(ids:ide,jds:jde,kpzd)
      real deni(ids:ide,jds:jde,km)
#  endif
# if defined O_npzd_extra_diagnostics
      real avej(ids:ide,jds:jde,kpzd), avej_D(ids:ide,jds:jde,kpzd)
      real gmax(ids:ide,jds:jde,kpzd), no3P(ids:ide,jds:jde,kpzd)
      real po4P(ids:ide,jds:jde,kpzd), po4_D(ids:ide,jds:jde,kpzd)
# endif
# endif
# if defined O_save_convection
      real totalk(ids:ide,jds:jde), vdepth(ids:ide,jds:jde)
      real pe(ids:ide,jds:jde)
# endif
# if defined O_save_carbon_carbonate_chem
      real sspH(ids:ide,jds:jde), ssCO3(ids:ide,jds:jde)
      real ssOc(ids:ide,jds:jde), ssOa(ids:ide,jds:jde)
      real sspCO2(ids:ide,jds:jde)
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      real ps(ids:ide,jds:jde)
# endif
# if defined O_stream_function
      real psi(ids:ide,jds:jde)
# endif
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
# if defined O_save_anisotropic_viscosity
      real mvisc(ids:ide,jds:jde), zvisc(ids:ide,jds:jde)
# endif
      real tgarea(ids:ide,jds:jde), ugarea(ids:ide,jds:jde)
      real tm(ids:ide,jds:jde,km), um(ids:ide,jds:jde,km)
      real time, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), zt_e(km+1), zw_e(km+1), c0, c1, c10, c100
      real c1e3, c1e4, c1e5, c1e6, p1, p001, p035, C2K, cal2J
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpijk(:,:,:), tmpijkm(:,:,:)
      real, allocatable :: tmpijl(:,:,:), tmpijlm(:,:,:)
      real, allocatable :: tmpi(:), tmpj(:), tmpk(:), tmpl(:)
      real, allocatable :: tmpie(:), tmpje(:), tmpke(:), tmple(:)

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e5 = 1.e5
      c1e6 = 1.e6
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      C2K = 273.15
      cal2J = 2.389e-05

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      if (ntrec .le. 0) call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
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
      kgs = 1
      kge = km
      kg  = kge-kgs+1
# if defined O_save_npzd
      lgs = max(1,kgs)
      lge = min(kpzd,kg)
      lg  = lge-lgs+1
# endif

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)
      kls = max(1,kgs)
      kle = min(km,kge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )
      allocate ( tmpijk(ils:ile,jls:jle,kls:kle) )
      allocate ( tmpijkm(ils:ile,jls:jle,kls:kle) )

# if defined O_save_npzd
      lls = max(1,lgs)
      lle = min(km,lge)

      allocate ( tmpijl(ils:ile,jls:jle,lls:lle) )
      allocate ( tmpijlm(ils:ile,jls:jle,lls:lle) )
# endif
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
        allocate ( tmpk(kgs:kge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )
        allocate ( tmpke(kgs:kge+1) )

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

        ic(1) = kg
        tmpk(kgs:kge) = zt(kgs:kge)
        call putvara ('depth', iou, kg, ib, ic, tmpk, c100, c0)
        tmpk(kgs:kge) = dzt(kgs:kge)
        call putvara ('G_dzt', iou, kg, ib, ic, tmpk, c100, c0)
        tmpk(kgs:kge) = zw(kgs:kge)
        call putvara ('depth_W', iou, kg, ib, ic, tmpk, c100, c0)

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

        ic(1) = kg + 1
        call edge_maker (1, zt_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zt_e(kgs:kge+1)
        call putvara ('depth_edges', iou, kg+1, ib, ic, tmpke
     &,   c100, c0)
        call edge_maker (2, zw_e, zt, dzt, zw, dzw, km)
        tmpke(kgs:kge+1) = zw_e(kgs:kge+1)
        call putvara ('depth_W_edges', iou, kg+1, ib, ic, tmpke
     &,   c100, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpk )
        deallocate ( tmpie )
        deallocate ( tmpje )
        deallocate ( tmpke )

# if defined O_save_npzd
        allocate ( tmpl(lgs:lge) )
        allocate ( tmple(lgs:lge+1) )
        ic(1) = lg
        tmpl(lgs:lge) = zt(lgs:lge)
        call putvara ('npzdlev', iou, lg, ib, ic, tmpl, c100, c0)
        ic(1) = lg + 1
        tmple(lgs:lge+1) = zt_e(lgs:lge+1)
        call putvara ('npzdlev_edges', iou, lg+1, ib, ic, tmple
     &,   c100, c0)
        deallocate ( tmpl )
        deallocate ( tmple )

# endif
!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = kmt(ils:ile,jls:jle)
        call putvara ('G_kmt', iou, ln, ib, ic, tmpij, c1, c0)
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
# if defined O_save_anisotropic_viscosity
        tmpij(ils:ile,jls:jle) = mvisc(ils:ile,jls:jle)
        call putvara('O_mvisc', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = zvisc(ils:ile,jls:jle)
        call putvara('O_zvisc', iou, ln, ib, ic, tmpij, c1, c0)
# endif
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
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      do n=1,nt
        tmpij(ils:ile,jls:jle) = stf(ils:ile,jls:jle,n)
        if (trim(mapt(n)) .eq. 'temp') then
          call putvaramsk('F_heat', iou, ln, ib, ic, tmpij, tmpijm
     &,     cal2J, c0)
        elseif (trim(mapt(n)) .eq. 'salt') then
          call putvaramsk('F_salt', iou, ln, ib, ic, tmpij, tmpijm
     &,     p1, c0)
        elseif (trim(mapt(n)) .eq. 'dic') then
          call putvaramsk('F_dic', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'alk') then
          call putvaramsk('F_alk', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'o2') then
          call putvaramsk('F_o2', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'po4') then
          call putvaramsk('F_po4', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'no3') then
          call putvaramsk('F_no3', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'c14') then
          call putvaramsk('F_c14', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call putvaramsk('F_cfc11', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call putvaramsk('F_cfc12', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'phyt') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call putvaramsk('F_'//trim(a3), iou, ln, ib, ic, tmpij
     &,     tmpijm, c1, c0)
        endif
      enddo
      tmpijm(ils:ile,jls:jle) = um(ils:ile,jls:jle,1)
      tmpij(ils:ile,jls:jle) = taux(ils:ile,jls:jle)
      call putvaramsk ('O_tauX', iou, ln, ib, ic, tmpij, tmpijm
     &, c10, c0)
      tmpij(ils:ile,jls:jle) = tauy(ils:ile,jls:jle)
      call putvaramsk ('O_tauY', iou, ln, ib, ic, tmpij, tmpijm
     &, c10, c0)
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      tmpij(ils:ile,jls:jle) = ps(ils:ile,jls:jle)
      call putvaramsk ('O_presssur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
# endif
# if defined O_stream_function
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      tmpij(ils:ile,jls:jle) = psi(ils:ile,jls:jle)
      call putvaramsk ('O_psi', iou, ln, ib, ic, tmpij, tmpijm
     &, c1e6, c0)
# endif
# if defined O_save_convection
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      tmpij(ils:ile,jls:jle) = totalk(ils:ile,jls:jle)
      call putvaramsk ('O_convlev', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = vdepth(ils:ile,jls:jle)
      call putvaramsk ('O_ventdep', iou, ln, ib, ic, tmpij, tmpijm
     &, c100, c0)
      tmpij(ils:ile,jls:jle) = pe(ils:ile,jls:jle)
      call putvaramsk ('O_convpe', iou, ln, ib, ic, tmpij, tmpijm
     &, c1e3, c0)
# endif
# if defined O_save_carbon_carbonate_chem
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      where (sspH(ils:ile,jls:jle) .le. 0.) tmpijm(ils:ile,jls:jle) = 0.
      tmpij(ils:ile,jls:jle) = sspH(ils:ile,jls:jle)
      call putvaramsk ('O_phsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = ssCO3(ils:ile,jls:jle)
      call putvaramsk ('O_co3sur', iou, ln, ib, ic, tmpij, tmpijm, c1
     &, c0)
      tmpij(ils:ile,jls:jle) = ssOc(ils:ile,jls:jle)
      call putvaramsk ('O_ocalcsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = ssOa(ils:ile,jls:jle)
      call putvaramsk ('O_oaragsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = sspCO2(ils:ile,jls:jle)
      call putvaramsk ('O_pco2sur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
# endif
# if defined O_save_npzd
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      tmpij(ils:ile,jls:jle) = procal(ils:ile,jls:jle)
      call putvaramsk ('O_caco3pro', iou, ln, ib, ic, tmpij, tmpijm
     &, c100, c0)
# endif

!-----------------------------------------------------------------------
!     write 4d data (x,y,z,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
# if defined O_save_carbon_totals
      toc(:,:,:) = 0.
# endif
      do n=1,nt
        tmpijk(ils:ile,jls:jle,kls:kle) = t(ils:ile,jls:jle,kls:kle,n)
        if (trim(mapt(n)) .eq. 'temp') then
          call putvaramsk('O_temp', iou, ln, ib, ic, tmpijk
# if defined O_units_temperature_Celsius
     &,     tmpijkm, c1, c0)
# else
     &,     tmpijkm, c1, -C2K)
# endif
        elseif (trim(mapt(n)) .eq. 'salt') then
          call putvaramsk('O_sal', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     p001, -p035)
        elseif (trim(mapt(n)) .eq. 'dic') then
          call putvaramsk('O_dic', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
# if defined O_save_carbon_totals
          toc(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
     &      + tmpijk(ils:ile,jls:jle,kls:kle)
# endif
        elseif (trim(mapt(n)) .eq. 'alk') then
          call putvaramsk('O_alk', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
        elseif (trim(mapt(n)) .eq. 'o2') then
          call putvaramsk('O_o2', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
        elseif (trim(mapt(n)) .eq. 'po4') then
          call putvaramsk('O_po4', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
        elseif (trim(mapt(n)) .eq. 'phyt') then
          call putvaramsk('O_phyt', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
# if defined O_save_carbon_totals
          toc(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
     &      + tmpijk(ils:ile,jls:jle,kls:kle)*redctn
# endif
        elseif (trim(mapt(n)) .eq. 'zoop') then
          call putvaramsk('O_zoop', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
# if defined O_save_carbon_totals
          toc(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
     &      + tmpijk(ils:ile,jls:jle,kls:kle)*redctn
# endif
        elseif (trim(mapt(n)) .eq. 'detr') then
          call putvaramsk('O_detr', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
# if defined O_save_carbon_totals
          toc(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
     &      + tmpijk(ils:ile,jls:jle,kls:kle)*redctn
# endif
        elseif (trim(mapt(n)) .eq. 'no3') then
          call putvaramsk('O_no3', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
        elseif (trim(mapt(n)) .eq. 'diaz') then
          call putvaramsk('O_diaz', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1e3, c0)
# if defined O_save_carbon_totals
          toc(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
     &      + tmpijk(ils:ile,jls:jle,kls:kle)*redctn
# endif
        elseif (trim(mapt(n)) .eq. 'c14') then
          call putvaramsk('O_c14', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call putvaramsk('O_cfc11', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call putvaramsk('O_cfc12', iou, ln, ib, ic, tmpijk, tmpijkm
     &,     c1, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call putvaramsk('tracer_'//trim(a3), iou, ln, ib, ic
     &,     tmpijk, tmpijkm, c1, c0)
        endif
      enddo
# if defined O_save_carbon_totals
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = toc(ils:ile,jls:jle,kls:kle)
      call putvaramsk('O_totcarb', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1, c0)
# endif
      tmpijkm(ils:ile,jls:jle,kls:kle) = um(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = u(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_velX', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c100, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = v(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_velY', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c100, c0)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = adv_vbt(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_velZ', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c100, c0)
# if defined O_gent_mcwilliams
      tmpijkm(ils:ile,jls:jle,kls:kle) = um(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) =
     &  adv_vetiso(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_gmvelX', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c100, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) =
     &  adv_vntiso(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_gmvelY', iou, ln, ib, ic, tmpijk, tmpijkm
     &,  c100, c0)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) =
     &  adv_vbtiso(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_gmvelZ', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c100, c0)
# endif
# if defined O_carbon && defined O_carbon_14
      tmpijk(ils:ile,jls:jle,kls:kle) = dc14(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_dc14', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1, c0)
# endif
# if defined O_save_kv
      tmpijk(ils:ile,jls:jle,kls:kle) =
     &  diff_cbt(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_diffZ', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e4, c0)
# endif
# if defined O_save_npzd
      tmpijlm(ils:ile,jls:jle,lls:lle) = tm(ils:ile,jls:jle,lls:lle)
      ib(3) = lls-lgs+1
      ic(3) = lle-lls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijl(ils:ile,jls:jle,lls:lle) = npp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytnpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz_Z(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoograz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = 
     & graz_Det(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morpt(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytrecy', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoopmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = excr(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_excret', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
# if defined O_npzd_extra_diagnostics
      tmpijl(ils:ile,jls:jle,lls:lle) = avej(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = avej_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = gmax(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_gmax', iou, ln, ib, ic, tmpijl, tmpijlm 
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = no3P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_no3P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = po4P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = po4_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
# endif
#  if defined O_npzd_nitrogen
      tmpijl(ils:ile,jls:jle,lls:lle) = npp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diaznpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = nfix(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_no3fix', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = deni(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_denitr', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
#  endif
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = remi(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrremi', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = expo(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrexp', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = expocal(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_caco3exp', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1, c0)
# endif

      deallocate ( tmpij )
      deallocate ( tmpijm )
      deallocate ( tmpijk )
      deallocate ( tmpijkm )
# if defined O_save_npzd
      deallocate ( tmpijl )
      deallocate ( tmpijlm )
# endif
#endif

      return
      end
