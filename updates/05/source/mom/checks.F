      subroutine checks (errorc, vmixset, hmixset)

# if defined O_mom
      implicit none

      integer i, k, j, ip, kr, jq, n, num_mw, jb, jjs, jj, ncrow, jw
      integer je, in, is, numk

      real critv, t1, dymin, dxmin, jrow, dzmin, xlmax, dtxl, num
      real dxdymn, fimax, fmax, dysq, dxsq, clix, h1, h2, hx

      logical errorc, vmixset, hmixset
      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "hmixc.h"
      include "iounit.h"
      include "levind.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"
# if defined O_xbts
      include "xbt.h"
# endif

!-----------------------------------------------------------------------
!     do consistency checks
!-----------------------------------------------------------------------

      write (stdout,'(/,20x,a,/)')
     &         'G E N E R A L    C O N S I S T E N C Y    C H E C K S'

      if (imt .lt. 3) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "imt" less than 3 is not allowed         '
        errorc = .true.
      endif

      if (jmt .lt. 4) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: parameter "jmt" less than 4 is not allowed          '
        errorc = .true.
      endif

# if defined O_pressure_gradient_average
      if (.not. eb) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: option "pressure_gradient_average" may need       '
     &,'             euler backward mixing timesteps to suppress the   '
     &,'             computational modes.                              '
      endif
# endif

# if defined O_linearized_advection
#  if defined O_fct
      write (stdout,'(/,(1x,a))')
     & '==> Error: fct" is incompatable with "linearized_advection".   '
      errorc = .true.
#  endif
#  if defined O_fourth_order_tracer_advection
      write (stdout,'(/,(1x,a))')
     & '==> Error: "fourth_order_tracer_advection" is incompatable     '
     &,'            with "linearized_advection".                       '
      errorc = .true.
#  endif
#  if defined O_quicker
      write (stdout,'(/,(1x,a))')
     & '==> Error: "quicker" is incompatable with                      '
     &,'           "linearized_advection".                             '
      errorc = .true.
#  endif
#  if !defined O_linearized_density
      write (stdout,'(/,(1x,a))')
     & '==> Error: "linearized_density" must be enabled when using     '
     &,'           "linearized_advection".                             '
      errorc = .true.
#  endif
# endif

# if defined O_consthmix
      if (hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "consthmix"  cannot be enabled because another      '
     &,'            horizontal mixing scheme has been enabled          '
        errorc = .true.
      else
        hmixset = .true.
      endif
# endif
      if (.not.hmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No horizontal mixing scheme has been enabled        '
        errorc = .true.
      endif

# if defined O_constvmix
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "constvmix"  cannot be enabled because another      '
     &,'            vertical mixing scheme has been enabled            '
        errorc = .true.
      else

!       set vmixset = true for enabeling "constvmix"

        vmixset = .true.
      endif
# endif
# if defined O_bryan_lewis_vertical
      if (.not. vmixset) then
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a/1x,a/1x,a))')
     & '==> Error: "bryan_lewis_vertical" tracer diffusion requires    '
     &,'           "constvmix" or "ppvmix" to be enabled.              '
     &,'           If "constvmix" is enabled, "bryan_lewis_vertical"   '
     &,'           coefficients replace kappa_h. If "ppvmix" is        '
     &,'           is enabled, "bryan_lewis_vertical" coefficients     '
     &,'           are added to the predicted coeffs.                  '
        errorc = .true.
      endif
# endif
# if defined O_bryan_lewis_horizontal
#  if !defined O_consthmix
      write (stdout,'(/,(1x,a/1x,a))')
     & '==> Error: "bryan_lewis_horizontal" tracer diffusion requires  '
     &,'           "consthmix" to be enabled.                          '
      errorc = .true.
#  endif
# endif

      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: No vertical mixing scheme has been enabled          '
        errorc = .true.
      endif
      if (.not.vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: there is no vertical mixing scheme enabled          '
        errorc = .true.
      endif
# if defined O_rigid_lid_surface_pressure
      write (stdout,9000)
     & '==> Warning: when using the "rigid_lid_surface_pressure" method'
     &,'             a checkerboard may develop in "ps". Using option  '
     &,'             "remove_ps_checkerboard" may help to remove it.   '
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  if defined O_hypergrid || defined O_oldrelax
      write (stdout,9000)
     & '==> Error: "hypergrid" or "oldrelax" is not allowed when       '
     &,'           using "implicit_free_surface" or                    '
     &,'           "rigid_lid_surface_pressure" options.               '
      errorc = .true.
#  endif
#  if defined O_sf_5_point
      write (stdout,9000)
     & '==> Warning: "sf_5_point" is ignored when using the            '
     &,'             "implicit_free_surface" or                        '
     &,'             "rigid_lid_surface_pressure" option.              '
     &,'             conjugate_gradient with 9 point numerics is used. '
#  endif
# endif

# if !defined O_stream_function && !defined O_implicit_free_surface
#  if !defined O_rigid_lid_surface_pressure
      write (stdout,9000)
     & '==> Error: choose "stream_function", "implicit_free_surface"   '
     &,'           or the "rigid_lid_surface_pressure"  method         '
      errorc = .true.
#  endif
# endif

# if defined O_stream_function
#  if !defined O_sf_5_point && !defined O_sf_9_point
      write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using      '
     &,'           the stream_function method                          '
      errorc = .true.
#  endif
#  if defined O_sf_5_point && defined O_sf_9_point
      write (stdout,9000)
     & '==> Error: choose "sf_5_point" or "sf_9_point" when using      '
     &,'           the stream_function method (not both)               '
        errorc = .true.
#  endif
# endif
# if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
      if (jmw .lt. 4) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 4 rows when using any'
     &,'           fourth order options                                '
        write (stdout,*)'           you have set jmw=',jmw
        errorc = .true.
      endif
      if (jmw .gt. 4) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 4 ("jmw"=4 will use the minimum memory)   '
        write (stdout,*)'             you have set jmw=',jmw
      endif
# else
      if (jmw .lt. 3) then
        write (stdout,9000)
     & '==> Error: the MW can not have fewer than 3 rows (jmw >2)      '
        write (stdout,*)'           you have set jmw=',jmw
        errorc = .true.
      endif
      if (jmw .gt. 3) then
        write (stdout,9000)
     & '==> Warning: "jmw" > 3 ("jmw"=3 will use the minimum memory)   '
      endif
# endif
      if (jmw .gt. jmt) then
        write (stdout,9000)
     & '==> Error: the MW can not have more rows than "jmt"            '
        write (stdout,*)'           you have set jmw=',jmw, ', jmt=',jmt
        errorc = .true.
      endif
      if (jmw .eq. jmt) then
        write (stdout,9000)
     & '==> Warning: The MW is open all the way ("jmw" = "jmt") which  '
     &,'             is the maximum memory configuration. Note that    '
     &,'             latitude rows are kept in the MW and not on disk! '
      endif

      if (nkflds .lt. 2) then
        write (stdout,9000)
     & '==> Error: "nkflds" must be at least 2                         '
        write (stdout,*)'           nkflds is set = ',nkflds
        errorc = .true.
      endif

# if defined O_restorst
      if (dampts(1) .eq. c0 .or. dampts(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping time scale "dampts" must be set > 0     '
     &,'           when enabeling option "restorst"                    '
        errorc = .true.
      endif
      if (dampts(1) .ne. dampts(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping time scale "dampts"'
     &,'             is set differently for T and S in option          '
     &,'             "restorst"? (they can be set the same in the      '
     &,'             namelist)                                         '
      endif
      if (dampdz(1) .eq. c0 .or. dampdz(2) .eq. c0) then
        write (stdout,9000)
     & '==> Error: the damping thickness "dampdz" must be set nonzero  '
     &,'           when enabeling the "restorst" otpion                '
        errorc = .true.
      endif
      if (dampdz(1) .ne. dampdz(2)) then
        write (stdout,9000)
     & '==> Warning: are you aware that the damping thickness "dampdz" '
     &,'             is set differently for T and S in option          '
     &,'             "restorst"?                                       '
      endif
# else
      if (dampts(1) .ne. c0 .or. dampts(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping time scale "dampts" is > zero but     '
     &,'             the "restorst" otpion is not enabled              '
      endif
      if (dampdz(1) .ne. c0 .or. dampdz(2) .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the damping thickness "dampdz" is > zero but      '
     &,'             the "restorst" otpion is not enabled              '
      endif
# endif

# if defined O_fourth_order_tracer_advection
      write (stdout,9000)
     & '==> Warning: "fourth_order_tracer_advection" has a more        '
     &,'             stringent stability requirement than the standard '
     &,'             second order advection. If problems arise, the    '
     &,'             time step may have to be reduced. Refer to the    '
     &,'             manual                                            '
# endif

# if defined O_isopycmix
#  if defined O_consthmix && !defined O_biharmonic
      write (stdout,9000)
     & '==> Note: consthmix will only affect mixing of momentum        '
     &,'          since isopycmix was specified for tracer diffusion.  '
     &,'          kappa_h and Ah will be used as background mixing     '
     &,'          coefficients                                         '
#  endif
#  if defined O_biharmonic
      write (stdout,9000)
     & '==> Error: biharmonic is incompatable with isopycmix           '
     &,'           (use "smagnlmix" instead)                           '
      errorc = .true.
#  endif
# else
#  if defined O_tidal_kv
      write (stdout,9000)
     & '==> Error: "isopycmix" must be enabled for "tidal_kv" to work. '
      errorc = .true.
#  endif
#  if defined O_gent_mcwilliams
      write (stdout,9000)
     & '==> Error: "isopycmix" must be enabled for "gent_mcwilliams"   '
     &,'           to work.                                            '
      errorc = .true.
#  endif
# endif
# if defined O_implicit_free_surface && defined O_stream_function
      write (stdout,9000)
     & '==> Error: "implicit_free_surface" and "stream_function" are   '
     &,'           incompatible. choose one or the other option.       '
      errorc = .true.
# endif
# if defined O_rigid_lid_surface_pressure && defined O_stream_function
      write (stdout,9000)
     & '==> Error: "rigid_lid_surface_pressure" and "stream_function"  '
     &,'           are incompatible. choose one or the other option.   '
      errorc = .true.
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
#  if defined O_diagnostic_surf_height
      write (stdout,9000)
     & '==> Warning: "diagnostic_surf_height" does not apply under     '
     &,'             the "rigid_lid_surface_pressure" or               '
     &,'             "implicit_free_surface" options.                  '
#  endif
# endif
# if defined O_stream_function
#  if !defined O_conjugate_gradient && !defined O_oldrelax
#   if !defined O_hypergrid
      write (stdout,9000)
     & '==> Error: either "conjugate_gradient", "oldrelax", or         '
     &,'           "hypergrid" must be used with stream_function.      '
      errorc = .true.
#   endif
#  endif
#  if defined O_oldrelax && defined O_hypergrid
      write (stdout,9000)
     & '==> Error: "oldrelax" and "hypergrid" are incompatible         '
     &,'           choose one or the other option.                     '
      errorc = .true.
#  endif
#  if defined O_oldrelax && defined O_conjugate_gradient
      write (stdout,9000)
     & '==> Error: "oldrelax" and "conjugate_gradient" are             '
     &,'           incompatible. choose one or the other option.       '
      errorc = .true.
#  endif
#  if defined O_oldrelax && defined O_conjugate_gradient
      write (stdout,9000)
     & '==> Error: "oldrelax" and "conjugate_gradient" are             '
     &,'           incompatible. choose one or the other option.       '
      errorc = .true.
#  endif
#  if defined O_hypergrid && defined O_conjugate_gradient
      write (stdout,9000)
     & '==> Error: "hypergrid" and "conjugate_gradient" are            '
     &,'           incompatible. choose one or the other option.       '
      errorc = .true.
#  endif
#  if defined O_hypergrid && defined O_conjugate_gradient
      write (stdout,9000)
     & '==> Error: "hypergrid" and "conjugate_gradient" are            '
     &,'           incompatible. choose one or the other option.       '
      errorc = .true.
#  endif
#  if defined O_sf_9_point
      write (stdout,9000)
     & '==> Warning: although "sf_9_point" conserves total energy      '
     &,'             it has a 2 grid point null space.                 '
#  endif
#  if defined O_sf_5_point
      write (stdout,9000)
     & '==> Warning: although "sf_5_point" has no null space, it does  '
     &,'             not conserve total energy.                        '
#  endif
#  if defined O_sf_9_point && defined O_oldrelax
      write (stdout,9000)
     & '==> Error: "oldrelax" cannot be used with 9 point numerics     '
     &,'           use sf_5_point with oldrelax.                       '
      errorc = .true.
#  endif
# endif
# if defined O_biharmonic && !defined O_consthmix
      write (stdout,9000)
     & '==> Error: "consthmix" must be enabled when using "biharmonic" '
      errorc = .true.
# endif
# if defined O_isopycmix
#  if defined O_consthmix
      if ((ah+ahisop) .gt. 1.e11) then
        write (stdout,9000)
     & '==> Error: "ahisop"+"ah" is too large for the                  '
     &,'           "isopycmix" mixing option                           '
        errorc = .true.
      endif
#  endif
# endif

      if (dtsf .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the external mode time step "dtsf"      '
        errorc = .true.
      endif

      if (dtuv .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the internal mode time step "dtuv"      '
        errorc = .true.
      endif

      if (dtts .le. c0) then
        write (stdout,9000)
     & '==> Error: need to set the density time step "dtts"            '
        errorc = .true.
      endif

      critv = 1.e-6
      if (mod(rundays,dtts*secday) .gt. critv) then
        t1 = nint(rundays/(dtts*secday))*dtts*secday
        write (stdout,9000)
     & '==> Warning: there must be an integral number of density time  '
     &,'             steps within "rundays" (the integration time).    '
        write (stdout,*) '               (changed "rundays" from     '
     &,   rundays,' days to ', t1,' days to insure this condition)     '
          rundays = t1
      endif

      if (itmb) then
        write (stdout,9000)
     & '==> Warning: "itmb" is set to "true". set it to "false" in     '
     &,'             subsequent runs to prevent the time independent   '
     &,'             basin mask from being written more than once.     '
     &,'             This reduces the size of the diagnostic file.     '
      endif

      if (itrmb) then
        write (stdout,9000)
     & '==> Warning: "itrmb" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
      endif

      if (itavg) then
        write (stdout,9000)
     & '==> Warning: "itavg" is set to "true". set it to "false" in    '
     &,'             subsequent runs to prevent the time independent   '
     &,'             region masks from being written more than once.   '
     &,'             This reduces the size of the diagnostic file.     '
      endif
# if defined O_meridional_tracer_budget
      if (mod(tmbint,dtts*secday) .gt. critv .and. tmbint .gt. c0)
     &    then
        t1 = nint(tmbint/(dtts*secday))*dtts*secday
        write (stdout,9000)
     & '==> Warning: "tmbint" does not contain an integral number of   '
     &,'             density time steps "dtts".                        '
        write (stdout,*) '               (changed "tmbint" from '
     &,   tmbint,' days to ', t1,' days to insure this condition)'
        tmbint = t1
      endif
      if (tmbint .eq. c0) then
        write (stdout,9000)
     & '==> Error: averaging interval "tmbint"  must be specified when '
     &,'           the "meridional_tracer_budget" option is enabled    '
        errorc = .true.
      endif
      if (tmbint .gt. tmbper) then
        write (stdout,9000)
     & '==> Warning: the interval "tmbint" exceeds the averaging       '
     &,'             period "tmbper" for option                        '
     &,'             "meridional_tracer_budget"                        '
      endif
      if (tmbint .lt. tmbper) then
        write (stdout,9000)
     & '==> Warning: averaging period "tmbper" exceeds interval        '
     &,'             "tmbint". Setting tmbper = tmbint                 '
        tmbper = tmbint
      endif
      if (tmbper .eq. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "tmbper" is zero. The        '
     &,'             average will be over only one time step!          '
      endif
# else
      if (tmbint .gt. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "tmbint" is > zero but the '
     &,'             the "meridional_tracer_budget" option is not on.  '
      endif
# endif
# if defined O_xbts
      if (mod(xbtint,dtts*secday) .gt. critv .and. xbtint .gt. c0)
     &  then
        t1 = nint(xbtint/(dtts*secday))*dtts*secday
        write (stdout,9000)
     & '==> Warning: "xbtint" does not contain an integral number of   '
     &,'    density time steps "dtts".                                 '
        write (stdout,*) '               (changed "xbtint" from '
     &,   xbtint,' days to ', t1,' days to insure this condition)'
        xbtint = t1
      endif
      if (xbtint .eq. c0) then
        write (stdout,9000)
     & '==> Warning: averaging interval "xbtint" = 0.0 implies no      '
     &,'             averaging when option "xbts" is enabled           '
      endif
      if (numxbt .eq. 0) then
        write (stdout,9000)
     &   '==> Error: the number of XBT stations detected is zero       '
     &,  '           at least one station must be specified            '
          errorc = .true.
      endif
      if (xbtint .gt. xbtper) then
        write (stdout,9000)
     & '==> Warning: the interval "xbtint" exceeds the  averaging      '
     &,'             period "xbtper" for option "xbts"                 '
      endif
      if (xbtint .lt. xbtper) then
        write (stdout,9000)
     & '==> Warning: averaging period "xbtper" exceeds interval        '
     &,'             "xbtint". Setting xbtper = xbtint                 '
          xbtper = xbtint
      endif
      if (xbtper .eq. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "xbtper" is zero. The        '
     &,'             average will be over only one time step!          '
      endif
# else
      if (xbtint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "xbtint"  is > zero but    '
     &,'             the "xbts" option is not enabled                  '
      endif
# endif
# if defined O_diagnostic_surf_height
      if (mod(dspint,dtts*secday) .gt. critv .and. dspint .gt. c0)
     &  then
        t1 = nint(dspint/(dtts*secday))*dtts*secday
      write (stdout,9000)
     & '==> Warning: "dspint" does not contain an integral number of   '
     &,'             density time steps "dtts".                        '
      write (stdout,*) '               (changed "dspint" from '
     &,   dspint,' days to ', t1,' days to insure this condition)'
        dspint = t1
      endif
      if (dspint .eq. c0) then
        write (stdout,9000)
     & '==> Warning: averaging interval "dspint" = 0.0 implies no      '
     &,'             averaging when option "diagnostic_surf_height" is '
     &,'             enabled                                           '
      endif
      if (dspint .gt. dspper) then
        write (stdout,9000)
     & '==> Warning: the interval "dspint" exceeds the averaging       '
     &,'             period "dspper" for option                        '
     &,'             "diagnostic_surf_height"                          '
      endif
      if (dspint .lt. dspper) then
        write (stdout,9000)
     & '==> Warning: averaging period "dspper" exceeds interval        '
     &,'             "dspint". Setting dspper = dspint                 '
          dspper = dspint
      endif
      if (dspper .eq. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging period "dspper" is zero. The        '
     &,'             average will be over only one time step!          '
      endif
# else
      if (dspint .ne. c0) then
        write (stdout,9000)
     & '==> Warning: the averaging interval "dspint"  is > zero but    '
     &,'             option "diagnostic_surf_height" is not enabled    '
      endif
# endif

      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     & '==> Warning: use of unequal time steps implies the transient   '
     &,'             response is unimportant and multiple equilibria   '
     &,'             do not exist.                                     '
      endif

!     check for mixing coefficients larger than stability permits

      dymin  = dyt(2)
      dxmin  = dxt(2)
      do jrow=2,jmtm1
        dymin  = min(dymin,dyt(jrow))
      enddo
      do i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
      enddo
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
      do k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
      enddo

      if (xlmax .gt. c1) then
        write (stdout,9000)
     & '==> Warning: use of accelerated time steps implies the         '
     &,'             transient response is unimportant and multiple    '
     &,'             equilibria do not exist. stability tests will     '
     &,'             use "dtts" multiplied by the maximum "dtxcel"     '
      endif

# if defined O_consthmix
      dtxl = dtts*xlmax
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
#  if defined O_isopycmix
        fimax = 0.
        do k=1,km
          do i=2,imtm1
            fimax = max(fimax,fisop(i,j,k))
          enddo
        enddo
        if ((dtxl*(ah+ahisop*fimax))/dxdymn .ge. p25) then
#  else
        if ((dtxl*ah)/dxdymn .ge. p25) then
#  endif
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
#  if defined O_isopycmix
     & '==> Warning: lateral diffusive criteria exceeded for "ah" +    '
     &,'             "ahisop". use a smaller "dtts", "dtxcel", and/or  '
     &,'             "ah" + "ahisop"                                   '
#  else
     & '==> Warning: lateral diffusive criteria exceeded for "ah".     '
     &,'              use a smaller "dtts", "dtxcel", and/or "ah"      '
#  endif
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
      num = 0
      do j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
        if ((dtuv*am)/dxdymn .ge. p25) then
          num = num + 1
          if (num .eq. 1) write (stdout,9000)
     & '==> Warning: lateral diffusive criteria exceeded for "am".     '
     &,'             use a smaller "dtuv" and/or "am"                  '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     &,                                     ',  j=',j
        endif
      enddo
# endif
# if defined O_shortwave
#  if !defined O_source_term
      write (stdout,9000)
     & '==> Error: ifdef "shortwave" is enabled but "source_term"      '
     &,'           is not. "source_term" needs to be enabled           '
      errorc = .true.
#  endif
#  if !defined O_embm
      write (stdout,9000)
     & '==> Warning: ifdef "shortwave" is enabled but "coupled" is not.'
     &,'             Is the shortwave component of radiation defined?  '
#  endif
# else
      if (dzt(1) .lt. 20.0e2) then
        write (stdout,9000)
     & '==> Warning: if shallow mixed layers develop, then enabling    '
     &,'             ifdef "shortwave" may help to deepen them. note   '
     &,'             that either you or the atmosphere must provide    '
     &,'             the solar short wave as a boundary condition.     '
      endif
# endif
# if defined O_constvmix
      do k=1,km
        if ((dtts*dtxcel(k)*kappa_h)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_h" '
     &,'             use a smaller "dtts", "dtxcel", and/or "kappa_h"  '
       write (stdout,'(a48,i3)') ' at level =',k
        endif
      enddo
      if ((dtuv*kappa_m)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     & '==> Warning: vertical diffusive criteria exceeded on "kappa_m" '
     &,'             use a smaller "dtuv" and/or "kappa_m"             '
      endif
# endif
# if defined O_implicitvmix
      write (stdout,9000)
     & '==> Warning: "implicitvmix" turns convective adjustment off.   '
     &,'             It will be handled by the vertical mixing         '
#  if defined O_fullconvect
      write (stdout,9000)
     & '==> Error: the full convective scheme is incompatible with     '
     &,'           option "implicitvmix"                               '
      errorc = .true.
#  endif
# else
#  if defined O_fullconvect
      write (stdout,9000)
     & '==> Warning: the full convective scheme is enabled.            '
     &,'             it will ignore "ncon" and remove all instability  '
#  else
      write (stdout,9000)
     & '==> Warning: the standard convective scheme is enabled. Verify '
     &,'             that "ncon" is set to what is wanted:             '
      write (stdout,'(32x,a,i3)') '"ncon" is currently set to ',ncon
#  endif
# endif

!     check range of implicit factors

# if !defined O_implicitvmix && !defined O_isopycmix
      if (aidif .ne. c0) then
        write (stdout,9000)
     & '==> Warning: "aidif" was not = 0.0  but "implicivmix" or       '
     &,'             "isopycmix" is not defined, so vertical mixing    '
     &,'              will be fully explicit.                          '
      endif
# endif
# if defined O_damp_inertial_oscillation
      if (acor .lt. p5) then
        write (stdout,9000)
     & '==> Error: "acor" must be between 0.5 & 1.0 when option        '
     &,'           damp_inertial_oscillation is enabled.               '
        errorc = .true.
      endif
#  if defined O_fct
      write (stdout,9000)
     & '==> Error: non-linear instability will result if "fct" is used '
     &,'           "damp_inertial_oscillation"                         '
      errorc = .true.
#  endif
# else
      if (acor .ne. 0) then
        write (stdout,9000)
     & '==> Error: "acor" must=0 when option damp_inertial_oscillation '
     &,'           is not enabled.                                     '
        errorc = .true.
      else

!       check for marginally resolved inertial oscillation

        fmax = epsln
        do jrow=2,jmtm1
          do i=2,imtm1
            fmax = max(fmax,abs(cori(i,jrow,1)))
          enddo
        enddo
        if (dtuv .gt. (1.0/6.0)*(c2*pi)/fmax) then
          write (stdout,9000)
     & '==> Error: the inertial oscillation is not resolved. reduce    '
     &,'           "dtuv" or use option "damp_inertial_oscillation"    '
          errorc = .true.
        endif
      endif
# endif
# if defined O_consthmix

!-----------------------------------------------------------------------
!     search for topographic instabilities (based  on the  work of
!     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
!-----------------------------------------------------------------------

      num   = 50
      do j=2,jmtm1
        dysq = dyt(j)**2
        do i=2,imtm1
        if (kmu(i+1,j-1) .ne. 0 .and. kmu(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmu(i+1,j-1))
            h2   = zw(kmu(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write(stdout,*)
              write (stdout,'(a,a,i4,a,i4,a)')
     &        '==> Warning: Killworth topographic roughness condition'
     &,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     &         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
        enddo
      enddo
# endif

!     verify that the domain boundary is valid

# if defined O_cyclic && defined O_solid_walls
      write (stdout,9000)
     & '==> Error: Both "cyclic" and "solid_walls cannot be enabled    '
      errorc = .true.
# endif
# if !defined O_cyclic && !defined O_solid_walls
      write (stdout,9000)
     & '==> Error: Either "cyclic" or "solid_walls must be enabled     '
      errorc = .true.
# endif
# if defined O_solid_walls
      jw = 0
      je = 0
      do j=1,jmt
        if (kmt(1,j) .ne. 0) jw = j
        if (kmt(imt,j) .ne. 0) je = j
      enddo
      if (je .ne. 0 .or. jw .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the eastern or western boundary. This is not        '
     &,'           correct when enabeling option "solid_walls"         '
        if (jw .ne. 0) write (stdout,*) '           at i=1 and j=',jw
        if (je .ne. 0) write (stdout,*) '           at i=imt and j=',je
      endif
# endif
      in = 0
      is = 0
      do i=1,imt
        if (kmt(i,1) .ne. 0) is = i
        if (kmt(i,jmt) .ne. 0) in = i
      enddo
      if (is .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the southern boundary.                              '
        write (stdout,*) '           at j=1 and i=',is
      endif
# if !defined O_symmetry
      if (in .ne. 0) then
        errorc = .true.
        write (stdout,9000)
     & '==> Error: The basin is not closed. "kmt" is non zero along    '
     &,'           the northern boundary.                              '
        write (stdout,*) '           at j=jmt and i=',in
      endif
# endif

!     verify that each ocean point is at least 2 levels deep

      numk = 0
      do jrow=1,jmt
        do i=1,imt
          if (kmt(i,jrow) .eq. 1) then
            numk = numk + 1
            errorc = .true.
            write (stdout,*)
     &       ' Error: kmt(',i,',',jrow,') = 1 is not allowed    '
          endif
        enddo
      enddo
      if (numk .ne. 0) then
        write (stdout,9000)
     & '==> Error: "kmt" must be at least 2 levels deep at all ocean   '
     &,'           points.                                             '
      endif

# if !defined O_quicker && defined O_ncar_upwind3
      write (stdout,9000)
     & '==> Error: "ncar_upwind3" only works when option               '
     &,'           option "quicker"  is enabled                        '
      errorc = .true.
# endif

# if defined O_fct
#  if defined O_fourth_order_tracer_advection
      write (stdout,9000)
     & '==> Error: "fct" is not compatable with option                 '
     &,'           "fourth_order_tracer_advection"                     '
      errorc = .true.
#  endif
#  if  defined O_quicker
      write (stdout,9000)
     & '==> Error: "fct" is not compatable with option                 '
     &,'           "quicker"                                           '
      errorc = .true.
#  endif
#  if defined O_fct_dlm1 && defined O_fct_dlm2
      write (stdout,9000)
     & '==> Error: more than one fct delimiter type is specified       '
     &,'           specify either fct_dlm1 or fct_dlm2                 '
     &,'           Default: using fct_dlm1                             '
      errorc = .true.
#  endif
#  if !defined O_fct_dlm1 && !defined O_fct_dlm2
      write (stdout,9000)
     & '==> Warning: fct delimiter type is not specified               '
     &,'             specify either fct_dlm1 or fct_dlm2               '
     &,'    Default: using fct_dlm1                                    '
#  endif
# else
#  if defined O_fct_dlm1 || defined O_fct_dlm2 || defined O_fct_3d
      write (stdout,9000)
     & '==> Warning: main fct option is not specified                  '
     &,'             minor fct options will be ignored                 '
#  endif
# endif
      write (stdout,'(/,20x,a,/)')
     &         ' E N D    C O N S I S T E N C Y    C H E C K S'
      if (errorc) stop '=>checks'

9000  format (/,(1x,a))

#endif
      return
      end
