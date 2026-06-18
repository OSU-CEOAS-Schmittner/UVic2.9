      subroutine set_time_switches
!=======================================================================

!               S E T   T I M E   S W I T C H E S

!     Subroutine set_time_switches sets switches in "switch.h"
!     that trigger periodically recurring events in the model
!     such as diagnostics and end-of-run.

!     Currently implemented switches include end-of-day, end-of-
!     week, end-of-two-weeks, end-of-month, end-of-year, end-of-
!     run, mid-month, and switches active at prespecified intervals
!     from either start of run, initial conditions, or any other
!     reference time the user chooses.  It is relatively
!     easy to add additional switches by following the models of
!     switches already provided.

!     input:

!       a switch.h file with all switch intervals in units of days.

!     outputs:

!       a whole collection of useful logical switches in the
!       switch.h file.
!=======================================================================

      implicit none

      integer i

      logical alarm, avg_alarm, timeless

      real shortest_mon

      include "stdunits.h"
      include "tmngr.h"
      include "switch.h"
      include "calendar.h"

!=======================================================================
!     set all time dependent logical switches (except "first")
!=======================================================================

!-----------------------------------------------------------------------

!      C A L E N D A R     A N D      C L O C K     S W I T C H E S

!      here are some examples of setting logicals based on the calendar
!      and clock. "dayoyr" is time in days since the start of a year.

!      alarms go off when current time + dt/2 is later then the alarm
!      time.  arguments of alarm routine:
!            1) index to the alarm setting time
!            2) index to the model time+dt/2
!            3) interval of the alarm in real days
!            4) index to the reference time
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of the day ?
!-----------------------------------------------------------------------

      eoday = alarm (ieoday, ihalfstep, 1.0, isunday)

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of the week ?
!-----------------------------------------------------------------------

      eoweek = alarm (ieoweek, ihalfstep, 7.0, isunday)

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of two weeks ?
!-----------------------------------------------------------------------

      eo2wks = alarm (ieo2wks, ihalfstep, 14.0, isunday)

!-----------------------------------------------------------------------
!     set initial alarm times in the first time step for the
!     end of month, mid month, end of year, and end of run switches.
!-----------------------------------------------------------------------

      if (first) then

!-----------------------------------------------------------------------
!       initialize end of month alarm
!-----------------------------------------------------------------------

        call getswitch (ieomon)
        call getfulltime (ialarm(ieomon))
        call copyfulltime (itime, ialarm(ieomon))
        month (ialarm(ieomon)) = month (ialarm(ieomon)) + 1
        if (month(ialarm(ieomon)) .gt. 12) then
          month(ialarm(ieomon)) = 1
          year (ialarm(ieomon)) = year(ialarm(ieomon)) + 1
        endif
        day   (ialarm(ieomon)) = 1
        hour  (ialarm(ieomon)) = 0
        minute(ialarm(ieomon)) = 0
        second(ialarm(ieomon)) = 0
        call expandtime (ialarm(ieomon))

        call gettime (iinterval(ieomon))
        shortest_mon = daypm(1)
        do i=2,12
          if (daypm(i) .lt. shortest_mon) shortest_mon = daypm(i)
        enddo
        call settime3 (iinterval(ieomon), shortest_mon)

!-----------------------------------------------------------------------
!       initialize mid month alarm
!-----------------------------------------------------------------------

        call getswitch (imidmon)
        call getfulltime (ialarm(imidmon))
        call copyfulltime (itime, ialarm(imidmon))
        day   (ialarm(imidmon)) = daysinmon(ialarm(imidmon))/2 + 1
        hour  (ialarm(imidmon)) = 12*modulo(daysinmon(ialarm(imidmon))
     &                            ,2)
        minute(ialarm(imidmon)) = 0
        second(ialarm(imidmon)) = 0
        call expandtime (ialarm(imidmon))
        if (timeless (ialarm(imidmon), itime)) then
          month (ialarm(imidmon)) = month(ialarm(imidmon)) + 1
          if (month(ialarm(imidmon)) .gt. 12) then
            month(ialarm(imidmon)) = 1
            year (ialarm(imidmon)) = year(ialarm(imidmon)) + 1
          endif
          day   (ialarm(imidmon)) = daysinmon(ialarm(imidmon))/2 + 1
          hour  (ialarm(imidmon)) = 12*modulo(daysinmon(ialarm(imidmon))
     &                              ,2)
          minute(ialarm(imidmon)) = 0
          second(ialarm(imidmon)) = 0
          call expandtime (ialarm(imidmon))
        endif

        call gettime (iinterval(imidmon))
        call settime3 (iinterval(imidmon), shortest_mon)

!-----------------------------------------------------------------------
!       initialize end of year alarm
!-----------------------------------------------------------------------

        call getswitch (ieoyear)
        call getfulltime (ialarm(ieoyear))
        call copyfulltime (itime, ialarm(ieoyear))
        year  (ialarm(ieoyear)) = year(ialarm(ieoyear)) + 1
        month (ialarm(ieoyear)) = 1
        day   (ialarm(ieoyear)) = 1
        hour  (ialarm(ieoyear)) = 0
        minute(ialarm(ieoyear)) = 0
        second(ialarm(ieoyear)) = 0
        call expandtime (ialarm(ieoyear))

        call gettime (iinterval(ieoyear))
        call settime3 (iinterval(ieoyear), real(yrlen))

      endif

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of the month ?
!-----------------------------------------------------------------------

      eomon = timeless (ialarm(ieomon), ihalfstep)
      on(ieomon) = eomon
      if (eomon) then
        iday(ialarm(ieomon)) = iday(ialarm(ieomon)) +
     &                         daysinmon(ialarm(ieomon))
        call expandtime2 (ialarm(ieomon))
      endif

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the middle of the month ?
!-----------------------------------------------------------------------

      midmon = timeless (ialarm(imidmon), ihalfstep)
      on(imidmon) = midmon
      if (midmon) then
        day(ialarm(imidmon)) = 1
        month(ialarm(imidmon)) = month(ialarm(imidmon)) + 1
        if (month(ialarm(imidmon)) .gt. 12) then
          month(ialarm(imidmon)) = 1
          year(ialarm(imidmon))  = year(ialarm(imidmon)) + 1
        endif
        call expandtime (ialarm(imidmon))
        day   (ialarm(imidmon)) = daysinmon(ialarm(imidmon))/2 + 1
        hour  (ialarm(imidmon)) = 12*modulo(daysinmon(ialarm(imidmon))
     &                            , 2)
        minute(ialarm(imidmon)) = 0
        second(ialarm(imidmon)) = 0
        call expandtime (ialarm(imidmon))
      endif

!-----------------------------------------------------------------------
!     is it within 1/2 time step of the end of the year ?
!-----------------------------------------------------------------------

      eoyear = timeless (ialarm(ieoyear), ihalfstep)
      on(ieoyear) = eoyear
      if (eoyear) then
        iday(ialarm(ieoyear)) = iday(ialarm(ieoyear)) +
     &                          daysinyear(ialarm(ieoyear))
        call expandtime2 (ialarm(ieoyear))
      endif

!-----------------------------------------------------------------------
!     is it the last time step of the run ?
!-----------------------------------------------------------------------

      eorun = timeless (ialarm(ieorun), ihalfstep)
      on(ieorun) = eorun

!-----------------------------------------------------------------------
!     is time mixing to be done now ?
!-----------------------------------------------------------------------

      if (nmix .eq. 0 .or. nmix .eq. 1) then
        leapfrog = .true.
      else
        if (mod(itt,nmix) .eq. 1) then
          leapfrog = .false.
        else
          leapfrog = .true.
        endif
      endif

!-----------------------------------------------------------------------
!     is it the first time step of an ocean segment?
!     is it the last time step of an ocean segment?
!-----------------------------------------------------------------------

      if (first) then
        osegs = .true.
        call getswitch (iosegs)
        call gettime (iinterval(iosegs))
        call settime3 (iinterval(iosegs), segtim)
      else
        osegs = osege
      endif
      on(iosegs) = osegs

      osege = alarm (iosege, ihalfstep, segtim, irunstart)

!=======================================================================

!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L

!     each interval switch needs three variables in common. The
!     following naming convention is used.

!         1) an interval (real) for diagnostic output (e.g,.  glenint)
!         2) a switch (logical) for the interval (e.g.,  glents )

!     the third is an internal variable needed by the time manager
!     to support calculation of the logical switch

!         3) an index (integer)                       (e.g., iglenint)

!     the user must specify the interval [e.g., glenint] for diagnostic
!     output in units of days. set_time_switches sets the corresponding
!     logical switch [e.g., glents] every time step. It is set to true
!     when within half a time step of the requested interval, otherwise
!     it is false. All decisions relating to the interval [e.g., glenint]
!     are based on the logical switch [e.g., glents].

!     internal time structures

!     The switch index [e.g., iglenint] is used to subscript into
!     internal arrays maintained by tmngr.F. The switch index is
!     allocated on the first call to function "alarm".
!     The array entry [e.g., iinterval(iglenint)] is a time index to the
!     internal representation of the interval [e.g., glenint].
!     The array entry [e.g., ialarm(iglenint)] is a time index to the
!     next time the alarm will be true.
!=======================================================================

#if defined O_tracer_averages

!-----------------------------------------------------------------------
!     are regional tracer averages to be done at this time ?
!-----------------------------------------------------------------------

      tavgts = alarm (itavgint, ihalfstep, tavgint, iref)
#endif
#if defined O_save_cross_flow

!-----------------------------------------------------------------------
!     is saving of diapycnal and isopycnal components of flow to be done
!     at this time ?
!-----------------------------------------------------------------------

      crossts = alarm (icrossint, ihalfstep, crossint, iref)
#endif
#if defined O_save_mixing_coeff

!-----------------------------------------------------------------------
!     is saving of estimated mixing coefficients to be done at
!     this time ?
!-----------------------------------------------------------------------

      cmixts = alarm (icmixint, ihalfstep, cmixint, iref)
#endif
#if defined O_save_convection_full

!-----------------------------------------------------------------------
!     is saving of temperature change due to explicit convection
!     to be done at this time ?
!-----------------------------------------------------------------------

      exconvts = alarm (iexconvint, ihalfstep, exconvint, iref)
#endif
#if defined O_energy_analysis

!-----------------------------------------------------------------------
!     are global energetics to be done at this time ?
!-----------------------------------------------------------------------

      glents = alarm (iglenint, ihalfstep, glenint, iref)
#endif
#if defined O_term_balances

!-----------------------------------------------------------------------
!     are tracer & momentum term balances to be done at this time ?
!-----------------------------------------------------------------------

      trmbts = alarm (itrmbint, ihalfstep, trmbint, iref)
#endif
#if defined O_matrix_sections

!-----------------------------------------------------------------------
!     are x-z printouts to be done at this time ?
!-----------------------------------------------------------------------

      prxzts = alarm (iprxzint, ihalfstep, prxzint, iref)
#endif
#if defined O_show_external_mode

!-----------------------------------------------------------------------
!     is the external mode printout to be done at this time ?
!-----------------------------------------------------------------------

      extts = alarm (iextint, ihalfstep, extint, iref)
#endif
#if defined O_tracer_yz

!-----------------------------------------------------------------------
!     is time to write the zonally averaged tracer components?
!-----------------------------------------------------------------------

      tyzts = alarm (ityzint, ihalfstep, tyzint, iref)
#endif
#if defined O_meridional_overturning

!-----------------------------------------------------------------------
!     is time to write the vertical_meridional stream function?
!-----------------------------------------------------------------------

      vmsfts = alarm (ivmsfint, ihalfstep, vmsfint, iref)
#endif
#if defined O_gyre_components

!-----------------------------------------------------------------------
!     are the gyre components to be done at this time ?
!-----------------------------------------------------------------------

      gyrets = alarm (igyreint, ihalfstep, gyreint, iref)
#endif
#if defined O_stability_tests

!-----------------------------------------------------------------------
!     is it time to do a stability diagnosis ?
!-----------------------------------------------------------------------

      stabts = alarm (istabint, ihalfstep, stabint, iref)

#endif

!-----------------------------------------------------------------------
!     is it time to save a restart ?
!-----------------------------------------------------------------------

      restts = alarm (irestint, ihalfstep, restint, iref)
#if defined O_show_zonal_mean_of_sbc

!-----------------------------------------------------------------------
!      is it time to do a boundary condition diagnosis ?
!-----------------------------------------------------------------------

      zmbcts = alarm (izmbcint, ihalfstep, zmbcint, iref)
#endif

!=======================================================================

!     S W I T C H E S    B A S E D    O N    A N    I N T E R V A L

!              A N D   A V E R A G I N G   P E R I O D

!     each averaging period switch needs five variables in common. The
!     following naming convention is used.

!         1) an interval (real) for diagnostic output    (e.g. xbtint  )
!         2) a switch (logical) for the interval         (e.g. xbtts   )
!         3) an averaging period (real)                  (e.g. xbtper  )
!         4) a switch (logical) for accumulating         (e.g. xbtperts)

!     the third is an internal variable needed by the time manager
!     to support calculation of the logical switches

!         5) an index (integer)                         (e.g. ixbtint  )

!     The user must specify the interval [e.g., xbtint] for diagnostic
!     output in units of days and the averaging period [e.g., xbtper]
!     in units of days. The averaging period may be less than or equal
!     to the interval. For example, if the interval is 30.0 days and the
!     averaging period is 5.0 days, results will be averaged over all
!     time steps within days 26, 27, 28, 29, and 30.

!     The logical switch for writing output at the specified interval
!     [e.g., xbtts] is set to true on the last time step of the
!     averaging period. The logical switch for accumulating results
!     [e.g., xbtperts] is true for all time steps within the averaging
!     period, otherwise it is false.

!     internal time structures

!     The index [e.g., ixbtint] is allocated on the first call to
!     function "avg_alarm". The array element iperiod(ixbtint) is an
!     index to the time structure for the internal representation of
!     "xbtper", and ilastsw(ixbtint) is the index of the switch that
!     flags the last time step of the accumulation period.
!     Depending on use,  ilastsw(ixbtint) may either be the index
!     of another "named" switch or the index of a new switch
!     allocated on the first time step.
!     In the latter case, iinterval(ilastsw(ixbtint)) is the index of
!     the time structure where "xbtint" is stored in internal form,
!     and ialarm(ilastsw(ixbtint)) is the index of the time when an
!     accumulation period will next end.
!     The variable nextts(ixbtint) is true whenever the next
!     time step will begin the accumulation period.
!=======================================================================
#if defined O_diagnostic_surf_height

!-----------------------------------------------------------------------
!     is it time to accumulate the forcing or "average and write" the
!      diagnostic surface height?
!-----------------------------------------------------------------------

      dspts = avg_alarm(idspint, ihalfstep, dspint, dspper, iref, 0)
      dspperts = on(idspint)
#endif
#if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     is it time to accumulate or "average and write" time mean
!     integrated data?
!-----------------------------------------------------------------------

      tsits = avg_alarm(itsiint, ihalfstep, tsiint, tsiper, iref, 0)
      tsiperts = on(itsiint)
#endif
#if defined O_time_averages

!-----------------------------------------------------------------------
!     is it time to accumulate or "average and write" time mean data
!     on the "averaging" grid?
!-----------------------------------------------------------------------

# if defined O_time_averages_monthly
      timavgts    = eomon
      timavgperts = .true.
# else
      timavgts = avg_alarm(itimavgint, ihalfstep, timavgint
     &,                        timavgper, iref, 0)
      timavgperts = on(itimavgint)
# endif
#endif
#if defined O_xbts

!-----------------------------------------------------------------------
!     is it time to accumulate or "average and write" XBT data?
!-----------------------------------------------------------------------

      xbtts = avg_alarm(ixbtint, ihalfstep, xbtint, xbtper, iref, 0)
      xbtperts = on(ixbtint)
#endif
#if defined O_mom_tbt

!-----------------------------------------------------------------------
!     is it time to accumulate or "average and write" XBT data?
!-----------------------------------------------------------------------

      tbtts = avg_alarm(itbtint, ihalfstep, tbtint, tbtper, iref, 0)
      tbtperts = on(itbtint)
#endif
#if defined O_meridional_tracer_budget

!-----------------------------------------------------------------------
!     is it time to accumulate or "average and write" the meridional
!     tracer budget?
!-----------------------------------------------------------------------

      tmbts = avg_alarm(itmbint, ihalfstep, tmbint, tmbper, iref, 0)
      tmbperts = on(itmbint)
#endif

      return
      end

      function alarm (isw, ihalf, timeint, irefer)

!=======================================================================
!     inputs:

!     isw     = index to the switch number for this switch. isw is
!               allocated on the first call to avg_alarm.

!     ihalf   = index to the time one half time step ahead of current
!               time

!     timeint = specified time interval of interest (in days. eg: 1.0
!               day, 36 hours (1.5 days), a week (7.0 days), two weeks
!               (14.0 days), 365.0 days ... etc)

!     irefer  = index to the time to which the alarm is referenced

!     outputs:

!     alarm:
!       if timeint > 0 then
!         alarm = true when the current time + dt/2 is later than the
!                 alarm time.  The alarm is then incremented
!       if timeint = 0 then
!         alarm = true (i.e., do this option every step)
!       if timeint < 0 then
!         alarm = false (i.e., this option is disabled)

!     if first=.T. then alarm sets the initial alarm time

!     output in common:

!     on(isw)      = true whenever alarm is true

!     internal values in common:
!     the user need not be concerned with these:

!     iinterval(isw) = index of time field where timeint is stored in
!                      internal form.

!     ialarm(isw)    = index of time when alarm will next be true
!=======================================================================

      implicit none

      integer isw, irefer, i, ihalf

      logical ltemp, alarm, timeless

      real timeint, realintervals, realdays

      include "tmngr.h"
      include "switch.h"

      if (first) then

!-----------------------------------------------------------------------
!       initialize the alarm on first iteration
!-----------------------------------------------------------------------

        call getswitch (isw)
        call gettime (ialarm(isw))
        call gettime (iinterval(isw))
        call settime3 (iinterval(isw), timeint)

        if (timeint .gt. 0.0) then
          if (timeless (irunstart, irefer)) then

!-----------------------------------------------------------------------
!           if reference time is later than run start time then the
!           initial alarm is set to reference time.  (No alarms will
!           go off until the time reaches the reference time.)
!-----------------------------------------------------------------------

            call copytime (irefer, ialarm(isw))
          else

!-----------------------------------------------------------------------
!           set alarm to first time of the form:  reftime + i * timeint
!           that is at least dt/2 after runstart.
!-----------------------------------------------------------------------

            call subtime (irunstart, irefer, ireftime)
            call addtime (ireftime, idtd2, ireftime)
            realintervals = realdays(ireftime)/timeint

!-----------------------------------------------------------------------
!           i = least integer greater than realintervals (ceiling)
!-----------------------------------------------------------------------

            i = int(realintervals) + 1
            i = i - int(i - realintervals)
            call multime (i, iinterval(isw), ialarm(isw))
            call addtime (ialarm(isw), irefer, ialarm(isw))
          endif
        else
          call addtime (itime, idt, ialarm(isw))
        endif
      endif

!-----------------------------------------------------------------------
!     check alarm
!-----------------------------------------------------------------------

      if (timeint .gt. 0.0) then
        ltemp = timeless (ialarm(isw), ihalf)
        if (ltemp) then

!-----------------------------------------------------------------------
!         increment the alarm time
!-----------------------------------------------------------------------

 100      continue
          call addtime (ialarm(isw), iinterval(isw), ialarm(isw))
          if (timeless (ialarm(isw), ihalf)) goto 100
        endif
      elseif (timeint .lt. 0.0) then
        ltemp = .false.
        call addtime (itime, idt, ialarm(isw))
      else
        ltemp = .true.
        call addtime (itime, idt, ialarm(isw))
      endif
      alarm = ltemp
      on(isw) = ltemp

      return
      end

      function avg_alarm(isw, ihalf, swinterval, period, irefer
     &,                    ilastswitch)

!=======================================================================
!     inputs:

!     isw     = index to the switch number for this switch. isw is
!               allocated on the first call to avg_alarm.

!     ihalf   = index to the time one half time step ahead of current
!               time

!     swinterval = specified time interval of the switch (in days.
!               eg: 1.0 day, 36 hours (1.5 days), a week (7.0 days),
!               two weeks (14.0 days), 365.0 days ... etc)

!     period = period of the averaging or accumulation part of the
!                interval (in days).
!                It is permissible for period and interval to be equal,
!                in which case, the logical on(isw) is always true, and
!                the logical lastts(isw) and the function value are
!                true become true once every interval days.

!     irefer  = index to the time to which the alarm is referenced

!     ilastswitch = index of switch which signals the last timestep of
!                   the interval.  If ilastswitch = 0, then a new,
!                   unnamed switch is allocated which turns on every
!                   interval days after the reference time.  If
!                   ilastswitch is nonzero, the argument interval is
!                   ignored and the "named" switch indexed by
!                   ilastswitch is used to signal the last time step
!                   of the interval.  Note that ilastswitch need not
!                   be a simple interval switch; it may be a calendar
!                   switch such as imidmon (mid month) or ieoyear
!                   (end of year).

!     outputs:

!     avg_alarm = true only on the last time step of the interval.

!     isw       = index to switch is allocated on first time step

!     outputs in common:

!     on(isw)      = true whenever current time step is within accumula-
!                    tion part of the interval.

!     lastts(isw)  = true only on last time step of interval.
!                    The alarm is then incremented

!     internal values in common:
!     the user need not be concerned with these:

!     nextts(isw)  = true when next time step will be first of interval.
!                    (i.e., on the time step nearest to period days
!                    before lastts(isw) will be true.)

!     iperiod(isw) = index of time field where length of the averaging
!                    period is stored in internal form.

!     ilastsw(isw) = index of switch which signals the last time step
!                    of the interval.  if ilastsw(isw) is a "named"
!                    switch, be sure to update it before isw is updated
!                    in the subroutine set_time_switches.
!=======================================================================

      implicit none

      integer isw, irefer, ihalf, ilastswitch

      logical avg_alarm, timeless, alarm, named_lastswitch

      real realdays, swinterval, period

      include "tmngr.h"
      include "switch.h"

      named_lastswitch = (ilastswitch .ne. 0)

      if (first) then

!-----------------------------------------------------------------------
!       allocate a new switch
!-----------------------------------------------------------------------

        call getswitch(isw)

!-----------------------------------------------------------------------
!       check for invalid user input
!-----------------------------------------------------------------------

        if (named_lastswitch) then
          ilastsw(isw) = ilastswitch
          swinterval = realdays(iinterval(ilastswitch))
        endif
        if (period .gt. swinterval) then
          print *, 'ERROR: switch period exceeds its interval'
          print *, '       period = ',period
          print *, '       interval = ',swinterval
          stop '==>avgalarm'
        endif

!-----------------------------------------------------------------------
!       initialize internal switch variables
!-----------------------------------------------------------------------

        call gettime (iperiod(isw))
        call settime3 (iperiod(isw), period)

!-----------------------------------------------------------------------
!       initialize switches lastts and nextts
!-----------------------------------------------------------------------

        if (named_lastswitch) then
          lastts(isw) = on(ilastsw(isw))
        else
          call getswitch(ilastsw(isw))
          lastts(isw) = alarm (ilastsw(isw), ihalf, swinterval, irefer)
        endif

        call subtime (ialarm(ilastsw(isw)), iperiod(isw), itmptime)
        call subtime (ihalf, idt, itmptime2)
        nextts(isw) = timeless (itmptime2, itmptime) .and.
     &                timeless (itmptime, ihalf)

!-----------------------------------------------------------------------
!       set on(isw) only if start of averaging interval is within a half
!       time step of runstart.
!-----------------------------------------------------------------------

        call subtime (itmptime2, idt, itmptime3)
        on(isw) = timeless (itmptime3, itmptime) .and.
     &            timeless (itmptime, itmptime2)

!-----------------------------------------------------------------------
!       turn lastts(isw) off if on(isw) is false
!-----------------------------------------------------------------------

        lastts(isw) = lastts(isw) .and. on(isw)

!-----------------------------------------------------------------------

      else

!-----------------------------------------------------------------------
!       this is not the first time step.
!       reset on(isw) based on events in previous call to increment_time
!-----------------------------------------------------------------------

        if (lastts(isw)) on(isw) = .false.
        if (nextts(isw)) on(isw) = .true.

!-----------------------------------------------------------------------
!       set lastts based on current call to increment_time.
!       end of run turns lastts(isw) on; however lastts(isw) must
!       never be true unless on(isw) is also true.
!-----------------------------------------------------------------------

        if (named_lastswitch) then
          lastts(isw) = on(ilastsw(isw))
        else
          lastts(isw) = alarm (ilastsw(isw), ihalf, swinterval, irefer)
        endif
        lastts(isw) = lastts(isw) .or. eorun
        lastts(isw) = lastts(isw) .and. on(isw)

!-----------------------------------------------------------------------
!       set nextts based on current call to increment_time.
!       these tests must follow the setting of lastts so that in the
!       case that the averaging period is equal to the entire interval,
!       the alarm setting for lastts will already be pushed ahead.

!       a period of 0 really means average over one time step.
!-----------------------------------------------------------------------

        call subtime (ialarm(ilastsw(isw)), iperiod(isw), itmptime)
        if (period .eq. 0) then
          call subtime (itmptime, idt, itmptime)
        endif
        call subtime (ihalf, idt, itmptime2)
        nextts(isw) = timeless (itmptime2, itmptime) .and.
     &                timeless (itmptime, ihalf)

      endif

      avg_alarm = lastts(isw)

      return
      end

      subroutine getswitch (isw)

      implicit none

      integer isw

      include "switch.h"

      nsw = nsw + 1
      if (nsw .gt. maxsw) then
        print *, 'ERROR:  not enough switches.'
        print *, '        increase maxsw = ',maxsw,' in switch.h'
        stop '==>getsw'
      endif
      isw = nsw

      return
      end

      subroutine set_eorun (runlen0, rununits0, rundays0)
!=======================================================================
!     initialize end of run alarm.
!     place dummy arguments in corresponding variables in "switch.h"

!     inputs:
!       runlen0   = length of run [see rununits for units]
!       rununits0 = time units ('days', 'months', or 'years' for runlen0
!     output:
!       rundays0  = length of run converted to (real) days

!     at present, arbitrary real run lengths in 'days' are handled.
!     run lengths in real 'months' or 'years' are rounded to the
!     nearest integral number of months or years before use.
!=======================================================================

      implicit none

      character(*) :: rununits0

      real runlen0, rundays0, realdays

      include "switch.h"
      include "stdunits.h"
      include "tmngr.h"

!     copy inputs to common in "switch.h"

      runlen = runlen0
      rununits = rununits0

!     calculate end of run based on rununits0

      call getswitch (ieorun)
      call getfulltime (ialarm(ieorun))
      call copyfulltime (irunstart, ialarm(ieorun))
      if (rununits0 .eq. 'years') then
        year(ialarm(ieorun)) = year(ialarm(ieorun)) + nint(runlen0)
        call expandtime (ialarm(ieorun))
        day  (ialarm(ieorun)) = min (day(irunstart),
     &                               daysinmon(ialarm(ieorun)))
        call expandtime (ialarm(ieorun))
        if (real(nint(runlen0)) .ne. runlen0) then
          print '(/,a,1pg14.7,a,i3,a,/)',
     &          'WARNING:  run length ',runlen0,' years rounded to ',
     &          nint(runlen0), ' years'
          runlen0 = real(nint(runlen0))
          runlen = runlen0
        endif
      elseif (rununits0 .eq. 'months') then
        month(ialarm(ieorun)) = month(ialarm(ieorun)) + nint(runlen0)
     &                                                - 1
        year (ialarm(ieorun)) = year (ialarm(ieorun)) + month(ieorun)/12
        month(ialarm(ieorun)) = modulo (month(ialarm(ieorun)), 12) + 1
        day  (ialarm(ieorun)) = 1
        call expandtime (ialarm(ieorun))
        if (real(nint(runlen0)) .ne. runlen0) then
          print '(/,a,1pg14.7,a,i3,a,/)',
     &          'WARNING:  run length ',runlen0,' months rounded to ',
     &          nint(runlen0), ' months'
          runlen0 = real(nint(runlen0))
          runlen = runlen0
        endif
        day  (ialarm(ieorun)) = min (day(irunstart),
     &                               daysinmon(ialarm(ieorun)))
        call expandtime (ialarm(ieorun))
      elseif (rununits0 .eq. 'days') then
        call settime3 (itemptime, runlen0)
        call addtime (irunstart, itemptime, ialarm(ieorun))
        call expandtime2 (ialarm(ieorun))
      else
        write (stdout, *)
     &       'Warning: No units given for run length--days assumed'
        call settime3 (itemptime, runlen0)
        call addtime (irunstart, itemptime, ialarm(ieorun))
        call expandtime2 (ialarm(ieorun))
      endif

      call subtime (ialarm(ieorun), irunstart, itemptime)
      rundays0 = realdays (itemptime)
      rundays  = rundays0

      return
      end

      subroutine initswitch

!=======================================================================
!     initialize switch indices for getswitch.
!=======================================================================

      implicit none

      include "switch.h"

      nsw = 0

      return
      end
