!=======================================================================

!               T I M E      M A N A G E R   M O D U L E

!     The time manager does three basic things for the model
!           1) It contains subroutines to keep, manipulate, and convert
!              all times for the model in internal representations in
!              arrays in "tmngr.h"
!           2) Subroutine increment_time updates times using the chosen
!              time step "dt" and calendar type without any roundoff or
!              drift for arbitrarily long integrations, even when using
!              32 bit arithmetic.
!           3) Subroutine set_time_switches sets switches in "switch.h"
!              that trigger periodically recurring events in the model
!              such as diagnostics and end-of-run.

!     Times are kept internal to the time manager in the primary form
!     of integer days (in the array "iday") and nonnegative integer
!     milliseconds in a fractional day (in the array "msday").  A unique
!     subscript for each time accesses all the arrays.  (See "tmngr.h"
!     for a list of predeclared time subscripts.)

!     Times that have a calendar reference are called "full times"
!     because in addition to the primary time fields "iday" and
!     "msday", they also have integer secondary time fields "year",
!     "month", "day", "hour", "minute", and "second", and a 32
!     character time stamp, "tstamp", in the form written and read
!     by the model. They also have tertiary fields "dayofyear",
!     "dayofweek", "daysinmon", and "daysinyear".  For such times,
!     the primary fields "iday" and "msday" represent days and
!     fractional days since the calendar base, (0,0) =
!     December 31, 1899 at the start of the day, 00:00:00.

!     For differences of two calendar times, the secondary calendar
!     fields make no sense, so these times are kept as "short times"
!     with only the two primary fields "iday" and "msday".

!     tmngr can keep track of time using either
!           a) a leap-year corrected Gregorian calendar
!           b) a constant 365-day year calendar
!           c) a constant month calendar, 12 months/year, with an
!              arbitrary number of days per month

!     Although the model does not currently use this feature, the
!     subroutine "tmngr" is designed to handle time increments, "dt",
!     varying at every time step.

!     Currently implemented switches include end-of-day, end-of-
!     week, end-of-two-weeks, end-of-month, end-of-year, end-of-
!     run, mid-month, and switches active at prespecified intervals
!     from either start of run, initial conditions, or any other
!     reference time the user chooses.  It is relatively
!     easy to add additional switches by following the models of
!     switches already provided.

!     Finally, this module provides a collection of utility
!     subroutines to perform arithmetic on the internal representations
!     of time, to convert between primary and secondary internal
!     representations, and between internal representations and
!     external representations such as real seconds or real days or
!     calendar year, month, day, hour, minute, and second.
!     The user who needs to work with time quantities, say to
!     write time stamps on a data set prepared by a user-written
!     program, might find some of them useful.
!=======================================================================

      subroutine tmngri (icyear, icmonth, icday, ichour, icmin, icsec
     &,                     rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
     &,                     idayrestart, msrestart
     &,                     runlen0, rununits0, rundays0, timestep)

!=======================================================================
!     initialize internal time variables and enter initial values
!     for externally defined times (initial conditions, user reference
!     and model time at start of run).
!=======================================================================

      implicit none

      character(*) :: rununits0

      integer icyear, icmonth, icday, ichour, icmin, icsec
      integer rfyear, rfmonth, rfday, rfhour, rfmin, rfsec
      integer idayrestart, msrestart, irefs, iddt, id, msdt, msec

      logical error, timeless

      real realdays, timestep, runlen0, rundays0

      include "stdunits.h"
      include "switch.h"
      include "tmngr.h"
      include "calendar.h"

      write (stdout,'(//,10x,a,/)') 'Time manager initialization'

      call inittime
      call initswitch
      first = .true.
      error = .false.
      call calendari (eqyear, eqmon, monlen,
     &            yrlen, daypm, msum, dayname, monname, error)
      if (error) then
        stop 'badcal'
      endif

!-----------------------------------------------------------------------
!     enter and print initial date and time and check bounds.
!-----------------------------------------------------------------------

      call getfulltime (initial)
      call setfulltime (initial, icyear, icmonth, icday, ichour, icmin
     &,                 icsec)
      write (stdout,9000) 'Initial Conditions: ', tstamp(initial)
      write (stdout,*) ' '
      call ckdate (initial, error)
      if (error) then
        stop '=>tmngri'
      endif

!----------------------------------------------------------------------
!     set model time counter.
!----------------------------------------------------------------------

      call gettime  (imodeltime)
      call settime2 (imodeltime, idayrestart, msrestart)
      write (stdout,'(a,1pg14.7,a/)')
     & '  Time since initial conditions =',
     & realdays(imodeltime), ' days'

!----------------------------------------------------------------------
!     calculate y/m/d and h/m/s of start of run.
!----------------------------------------------------------------------

      call getfulltime (irunstart)
      call addtime (initial, imodeltime, irunstart)
      call expandtime2 (irunstart)
      call getfulltime (itime)
      call copyfulltime (irunstart, itime)
      stamp = tstamp(itime)
      pstamp = stamp

!-----------------------------------------------------------------------
!     calculate real output quantities:
!       relyr  = years of model time since initial conditions
!       dayoyr = days since start of current year
!-----------------------------------------------------------------------

      dayoyr = dayofyear(itime) - 1 + (msday(itime)/(daylen*1000.0))
      relyr  = year(itime) - year(initial) + dayoyr/daysinyear(itime)
     &       - (dayofyear(initial) - 1
     &       + (msday(initial)/(daylen*1000.0)))/daysinyear(initial)
      prelyr = relyr

      write (stdout,9000) 'Start of Run: ', tstamp(irunstart)
      write (stdout,'(a,i10/)') 'Corresponding to time step "itt" =',itt

!----------------------------------------------------------------------
!       select reference time for computing diagnostic switches
!----------------------------------------------------------------------

      irefs = 0
      call getfulltime (iref)
      if (refrun) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refrun = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  the beginning of each run.'
        call copyfulltime (irunstart, iref)
      endif
      if (refinit) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refinit = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  the initial conditions time.'
        call copyfulltime (initial, iref)
      endif
      if (refuser) then
        irefs = irefs + 1
        write (stdout,*)
     &  ' "refuser = .true." selected.'
        write (stdout,'(a,a)')
     &  '  intervals for diagnostic switches are referenced to ',
     &  '  user specified date and time.'
        call getfulltime (iuser)
        call setfulltime (iuser, rfyear, rfmonth, rfday, rfhour, rfmin
     &,                   rfsec)
        call ckdate (iuser, error)
        if (error) then
          stop '=>tmngri'
        endif
        write (stdout,9000) '    Reference time:     ', tstamp(iuser)
        write (stdout,*) ' '
        call copyfulltime (iuser, iref)
      endif

      if (irefs .ne. 1) then
        write (stdout, *) 'You must choose exactly one of the'
     &  // ' options: refrun, refinit, or refuser.'
        stop '=>tmngr'
      endif

      call gettime (idt)
      call gettime (idtd2)
      call gettime (iusertime)
      call gettime (ireftime)
      call gettime (iruntime)
      call gettime (ihalfstep)

      call getfulltime (itemptime)
      call getfulltime (itemptime2)
      call gettime (itmptime)
      call gettime (itmptime2)
      call gettime (itmptime3)

!-----------------------------------------------------------------------
!     set a reference to Sunday Jan 0, 1900, 0:00:00, the base date
!     or earier if necessary.
!-----------------------------------------------------------------------

      call gettime  (isunday)
      call settime2 (isunday, 0, 0)
      if (timeless (initial, isunday)) then
        iday(isunday) = 14*int(iday(initial)/14) - 14
      endif

!-----------------------------------------------------------------------
!     initialize idt as "timestep" rounded to the nearest millisecond
!-----------------------------------------------------------------------

      iddt = id (timestep)
      msdt = msec (timestep)
      call settime2 (idt, iddt, msdt)

      call set_eorun (runlen0, rununits0, rundays0)

      write (stdout,'(/,10x,a,//)') 'Initialization completed'

9000  format (a,a)
      return
      end

      subroutine increment_time (dt)
!=======================================================================

!               I N C R E M E N T   T I M E

!     increment_time keeps track of model time using either
!           a) leap-year corrected Gregorian calendar
!           b) constant 365-day year calendar
!           c) constant month calendar, 12 months/year, arbitrary
!              number of days per month

!     date and time are accurate to within one millisecond for arbitrary
!     length time steps (even on computers with 32 bit word lengths!).

!               julian base = Jan 0, 1900 at 00:00:00

!       For accuracy, all fundamental times are kept
!       in the form:  integer days (with Jan 1, 1900=1)
!                     non-negative integer milliseconds within the day

!     input:

!       dt     = length of time step in seconds. (need not be constant)

!     outputs:

!       dt     = dt rounded to nearest millisec if needed

!       updated time fields: year, month, day, hour, minute, second,
!                            tstamp, dayofyear, dayofweek, daysinmonth,
!                            daysinyear
!       times updated:
!         itime  = "absolute" (y/m/d...) time after adding dt
!         ihalfstep = dt/2 beyond itime
!         imodeltime = time since initial conditions

!       stamp   = 32 character time stamp (m/d/y h:m:s)
!       pstamp  = stamp returned in previous call to increment_time

!       relyr  = model time in years since initial conditions
!       dayoyr = days + fractional days since start of calendar year
!=======================================================================

      implicit none

      integer iddt, id, msdt, msec, iddtd2, msdtd2

      logical timeequal

      real dt, realsecs

      include "stdunits.h"
      include "tmngr.h"
      include "switch.h"
      include "calendar.h"

!=======================================================================
!     date and time calculations (done every time step)
!=======================================================================

!     set flag "first" if first iteration of a run.
!     this flag must be set BEFORE the time is incremented.

      first = timeequal (itime, irunstart)

!-----------------------------------------------------------------------
!     round dt to the nearest millisecond
!-----------------------------------------------------------------------

      iddt = id (dt)
      msdt = msec (dt)
      call settime2 (idt, iddt, msdt)
      dt   = realsecs(idt)

!-----------------------------------------------------------------------
!     calculate half time step
!-----------------------------------------------------------------------

      iddtd2 = iddt/2
      msdtd2 = msdt/2 + 43200000*modulo (iddt, 2)
      call settime2 (idtd2, iddtd2, msdtd2)

!-----------------------------------------------------------------------
!     save previous values of stamp and relyr
!-----------------------------------------------------------------------

      pstamp = tstamp(itime)
      prelyr = relyr

!-----------------------------------------------------------------------
!     increment time counters
!-----------------------------------------------------------------------

      call addtime (itime, idt, itime)
      call expandtime2 (itime)

!     set current time stamp for MOM

      stamp = tstamp(itime)

      call addtime (itime, idtd2, ihalfstep)

!-----------------------------------------------------------------------
!     calculate number of days since reference time and start of run
!     all times are of form: (integer days, fractional day in millisec)
!----------------------------------------------------------------------

      call subtime (itime, irunstart, iruntime)
      call subtime (itime, initial, imodeltime)
      if (refuser) then
        call subtime (itime, iuser, iusertime)
      endif

!-----------------------------------------------------------------------
!     calculate real output quantities:
!       relyr  = years of model time since initial conditions
!       dayoyr = days since start of current year
!-----------------------------------------------------------------------

      dayoyr = dayofyear(itime) - 1 + (msday(itime)/(daylen*1000.))
      relyr  = year(itime) - year(initial) + dayoyr/daysinyear(itime)
     &       - (dayofyear(initial) - 1
     &       + (msday(initial)/(daylen*1000.)))/daysinyear(initial)

      return
      end

      subroutine calendari (eqyear, eqmon, monlen,
     &                   yrlen, daypm, msum, dayname, monname, error)
!=======================================================================
!     set up the calendar by choosing one of the following:

!        a) fully leap-year corrected calendar     (eqyear=F eqmon=F)
!        b) equal 365-day years (variable months)  (eqyear=T eqmon=F)
!        c) 12 equal months of "monlen" days each  (eqyear=T eqmon=T)

!     inputs:
!        eqyear, eqmon = logicals degined as in "a", "b", and "c"
!        monlen = days per month when using option "c"

!     outputs:
!        yrlen   = length of year in days
!        daypm   = days per month
!        msum    = accumulated days per month
!        dayname = character day names
!        monname = character month names
!        error   = .t. if something went wrong
!=======================================================================

      implicit none

      integer i, monlen, n, nmonth, nday, yrlen
      parameter (nmonth = 12, nday = 7)

      logical eqyear, eqmon, error

      include "stdunits.h"

      character(10) :: daynamei(nday)
      character(12) :: monnamei(nmonth)
      character(*) :: dayname(nday), monname(nmonth)

      integer daypm(nmonth), daypmi(nmonth), msum(nmonth)

      data daypmi/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      data daynamei/'sunday', 'monday', 'tuesday', 'wednesday',
     &            'thursday', 'friday', 'saturday'/

      data monnamei /'january', 'febuary', 'march', 'april', 'may'
     &,           'june', 'july', 'august', 'september', 'october'
     &,           'november', 'december'/

!-----------------------------------------------------------------------
!     first check that a valid combination of eqyear and eqmon
!     has been specified
!-----------------------------------------------------------------------

      if (eqmon .and. (.not. eqyear)) then
        write (stdout,999) eqyear, eqmon
        stop 'badcal'
      endif

      if (.not. (eqyear .and. eqmon)) then
        do i = 1, nmonth
          daypm(i) = daypmi(i)
        enddo
      else
        error = error .or. (monlen .le. 0)
        do i = 1, nmonth
          daypm(i) = monlen
        enddo
      endif

      msum(1) = 0
      do i = 2,nmonth
        msum(i) = msum(i-1) + daypm(i-1)
      enddo
      yrlen = msum(nmonth) + daypm(nmonth)

!     initialize day and month names

      do n=1,nday
        dayname(n) = daynamei(n)
      enddo

      do n=1,nmonth
        monname(n) = monnamei(n)
      enddo

      write (stdout,'(a)') 'Calendar selection:'
      if (.not. eqyear) then
        print '(/,a)', '  leap year corrected calendar'
      else
        print '(/,a,i4,a)', '  equal years of ', yrlen, ' days'
      endif

      print '(a,12i3,/)', '  days per month =',daypm

      return

999   format(/' error => An inappropriate calendar type was selected  ',
     &   /'     eqyear was set to ',l1,'   eqmon was set to ',l1,
     &   /'     valid combinations are:',
     &   /'     fully leap-year corrected calendar  (eqyear=F eqmon=F)'
     &   /'     equal 365-day years                 (eqyear=T eqmon=F)'
     &   /'     12 equal months of monlen days each (eqyear=T eqmon=T)')

      end

      subroutine ckdate (i, error)
!=======================================================================
!     do bounds checking on clock parameters
!     year is not checked since all years are ok.
!     at present, one extra day is allowed in month 2, even if year is
!     not a leap year or if there are no leap years.
!=======================================================================

      implicit none

      integer i

      logical error

      include "stdunits.h"
      include "tmngr.h"
      include "calendar.h"

      error = .false.

      if (month(i) .lt. 1 .or. month(i) .gt. 12) then
        write (stdout,*) ' Error:  month is out of bounds'
        error = .true.
      endif

      if ((day(i) .lt. 1 .or. day(i) .gt. daypm(month(i))) .and. .not.
     &    (month(i) .eq. 2 .and. day(i) .eq. daypm(2)+1)) then
        write (stdout,*) ' Error:  day is out of bounds'
        error = .true.
      endif

      if (hour(i) .lt. 0 .or. hour(i) .gt. 23) then
        write (stdout,*) ' Error:  hour is out of bounds'
        error = .true.
      endif

      if (minute(i) .lt. 0 .or. minute(i) .gt. 59) then
        write (stdout,*) ' Error:  minute is out of bounds'
        error = .true.
      endif

      if (second(i) .lt. 0 .or. second(i) .gt. 59) then
        write (stdout,*) ' Error:  second is out of bounds'
        error = .true.
      endif

      return
      end

      subroutine d2ymd (inday,
     &                  iyear, imonth, iday, idoy, idow, ndim, ndiy)
!======================================================================
!     d2ymd takes the number of days since Dec 31, 1899 and converts
!     to year, month, day, and other defining quantities.
!     ie. inday=1   yields  Jan 1, 1900.
!     if eqyear=.true. then calculation is done in subroutine d2ymdc

!     inday  = input days since Dec. 31, 1899
!     eqyear = .false. ==> use leap year corrected calendar
!     eqyear = .true.  ==> use constant length year calendar

!     output:
!       iyear  = output year
!       imonth = output month
!       iday   = output day of the month
!       idoy   = output day of the current year
!       idow   = output day of the week   1=sun - 7=sat
!       ndim   = number of days in the current month
!       ndiy   = number of days in the current year
!======================================================================

      implicit none

      integer inday, iyear, imonth, iday, idoy, idow, ndim, ndiy
      integer nday, ibasyr, nfh, nhund, nfour, nex, m

      logical leap

      include "calendar.h"

!     define the number of days per month for a non leap year (daypmi)

      integer daypmi(12)
      data daypmi/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

!-----------------------------------------------------------------------
!     equal-year calculations are handled by d2ymdc
!-----------------------------------------------------------------------

      if (eqyear) then
        call d2ymdc (inday,
     &               iyear, imonth, iday, idoy, idow, ndim, ndiy)
        return
      endif

!-----------------------------------------------------------------------
!     calculate day of the week (idow)
!     inday = 0, i.e.,  Dec 31, 1899 was a Sunday.
!-----------------------------------------------------------------------

      idow = modulo (inday, 7) + 1

!-----------------------------------------------------------------------
!     convert to days after january 1, 1601 using a virtual calendar
!     consistent with current leap-year systems to extend backward.
!     (ie. 1700, 1800, 1900 are not considered leap-years but 2000
!     is a leap-year)
!-----------------------------------------------------------------------

      nday  = inday + 109206
      ibasyr = 1601

!-----------------------------------------------------------------------
!     make corrections for years before 1601 (assuming modern leap year
!     conventions)  note that 400 years does not affect the day of the
!     week.   146097 days = 20871 weeks exactly.
!-----------------------------------------------------------------------

      if (nday .lt. 0) then
        nfh     = 1 + (-nday - 1)/146097
        ibasyr  = ibasyr - 400*nfh
        nday    = nday + 146097*nfh
      endif

!-----------------------------------------------------------------------
!     peel off year.  Accurate from before Jan 1, 1900 until the
!     julian day overflows machine limitations.
!----------------------------------------------------------------------

      nfh   = nday/146097
      nday  = modulo (nday, 146097)
      nhund = nday/36524
      if (nhund .gt. 3) then
        nhund = 3
        nday  = 36524
      else
        nday  = modulo (nday, 36524)
      endif
      nfour = nday/1461
      nday  = modulo (nday, 1461)
      nex   = nday/365
      if (nex .gt. 3) then
        nex   = 3
        nday  = 365
      else
        nday  = modulo (nday, 365)
      endif
      leap  = (nex .eq. 3) .and. ((nfour .ne. 24) .or. (nhund .eq. 3))
      if (leap) then
        ndiy = 366
      else
        ndiy = 365
      endif
      iyear = ibasyr + 400*nfh + 100*nhund + 4*nfour + nex
      idoy  = nday + 1

!-----------------------------------------------------------------------
!     peel off month and day.
!-----------------------------------------------------------------------

      iday = idoy
      do m=1,12
        if (leap .and. (m .eq. 2)) then
          if (iday .le. (daypmi(2)+1)) then
            imonth = 2
            ndim = daypmi(2) + 1
            goto 200
          endif
          iday = iday - daypmi(2) - 1
        else
          if (iday .le. daypmi(m))  then
            imonth = m
            ndim = daypmi(m)
            goto 200
          endif
          iday = iday - daypmi(m)
        endif
      enddo

 200  continue

      return
      end

      subroutine d2ymdc (inday,
     &                   iyear, imonth, iday, idoy, idow, ndim, ndiy)
!======================================================================
!     inverse of "ymd2dc"
!     d2ymdc takes the number of days since Dec 31, 1899 (or the last
!     day of December 1899 in case equal months of length < 31 days are
!     used) and converts to year, month, day, and other defining
!     quantities.
!     eg: inday=1   yields  Jan 1, 1900.

!     equal length years are assumed (no leap years here)

!     input:
!       inday  = days since Dec. 31, 1899

!     output:
!       iyear  = output year
!       imonth = output month
!       iday   = output day of the month
!       idoy   = day of the current year
!       idow   = day of the week   1=sun - 7=sat
!       ndim   = number of days in the current month
!       ndiy   = number of days in the current year
!======================================================================

      implicit none

      integer idow, inday, iyr, iyear, idoy, iday, m, imonth, ndim, ndiy

      include "calendar.h"

!-----------------------------------------------------------------------
!     calculate day of the week (sunday=1)
!     inday = 0, i.e., Dec 31, 1899 was a Sunday.
!-----------------------------------------------------------------------

      idow = modulo (inday, 7) + 1

!-----------------------------------------------------------------------
!     peel off the year
!-----------------------------------------------------------------------

      if (inday .le. 0) then
        iyr = (inday / yrlen) - 1
      else
        iyr = (inday - 1) / yrlen
      endif
      iyear = iyr + 1900
      idoy  = inday - yrlen*iyr

!-----------------------------------------------------------------------
!     peel off the month and day.
!-----------------------------------------------------------------------

      iday = idoy
      do m=1,12
        if (iday .le. daypm(m))  then
          imonth = m
          ndim = daypm(m)
          goto 200
        endif
        iday = iday - daypm(m)
      enddo

 200  continue
      ndiy = yrlen

      return
      end

      subroutine ymd2d (iyear, imonth, iday,
     &                  nday, idoy, idow, ndim, ndiy)
!=======================================================================
!     inverse of "d2ymd"
!     ymd2d takes a date in the year, month, day form and converts it
!     to days since Dec 31. 1899 (nday).  It also needs the
!     cumulative sums of days from the beginning of the year to each
!     month "msum". It returns the number of days in the
!     specified year (ndiy).

!     equal year calculations are passed on to subroutine ymd2dc

!     input:
!       iyear  = input year
!       imonth = input month
!       iday   = input day of the month

!     output:
!       nday   = days since Dec. 31, 1899
!       idoy   = day of the current year
!       idow   = day of the week   1=sun - 7=sat
!       ndim   = number of days in the current month
!       ndiy   = number of days in the current year
!======================================================================

      implicit none

      integer iyear, imonth, iday, nday, idoy, idow, ndim, ndiy, iyr
      integer nfh

      logical leap

      include "calendar.h"

!-----------------------------------------------------------------------
!     have subroutine ymd2dc do equal-year calculations
!-----------------------------------------------------------------------

      if (eqyear) then
        call ymd2dc (iyear, imonth, iday, nday, idoy, idow, ndim, ndiy)
        return
      endif

      if (mod (iyear, 400) .eq. 0) then
        leap = .true.
      elseif ((mod(iyear,4) .eq. 0) .and. (mod(iyear,100) .ne. 0)) then
        leap = .true.
      else
        leap = .false.
      endif
      if (leap) then
        ndiy = 366
      else
        ndiy = 365
      endif
      nday = iday + msum(imonth)
      if (leap .and. (imonth .gt. 2)) then
        nday = nday + 1
      endif
      idoy = nday
      iyr = iyear - 1601

!-----------------------------------------------------------------------
!     make corrections for years before 1601.
!-----------------------------------------------------------------------

      if (iyr .lt. 0) then
        nfh  = 1 + (-iyr)/400
        nday = nday - 146097*nfh
        iyr  = iyr + 400*nfh
      endif
      nday = nday + 365*iyr

!-----------------------------------------------------------------------
!     correct for leap-years between 1601 and iyear-1.
!     A virtual calendar consistent with current leap-year systems
!     is used to extend backward.  (ie. 1700, 1800, 1900 are not
!     considered leap-year but 2000 is a leap-year)
!-----------------------------------------------------------------------

      iyr = iyr/4
      nday = nday + iyr
      iyr = iyr/25
      nday = nday - iyr
      iyr = iyr/4
      nday = nday + iyr

!-----------------------------------------------------------------------
!     nday is now in days since Dec. 31, 1600.  Convert to days since
!     Dec 31, 1899.
!-----------------------------------------------------------------------

      nday = nday - 109207
      idow = modulo (nday, 7) + 1
      ndim = daypm(imonth)
      if (leap .and. (imonth .eq. 2)) then
        ndim = ndim + 1
      endif

      return
      end

      subroutine ymd2dc (iyear, imonth, iday,
     &                   nday, idoy, idow, ndim, ndiy)
!=======================================================================
!     inverse of "d2ymdc"
!     ymd2dc takes a date in the year, month, day form and converts
!     to days since Dec 31. 1899 (nday).  It also returns the number
!     of days in the specified year (ndiy).

!     ymd2dc uses a year of constant length with no leap years

!     input:
!       iyear  = input year
!       imonth = input month
!       iday   = input day of the month

!     output:
!       nday   = days since Dec. 31, 1899
!       idoy   = day of the current year
!       idow   = day of the week   1=sun - 7=sat
!       ndim   = number of days in the current month
!       ndiy   = number of days in the current year
!=======================================================================

      implicit none

      integer idoy, iday, imonth, iyr, iyear, nday, idow, ndim, ndiy

      include "calendar.h"

      idoy = iday + msum(imonth)
      iyr = iyear - 1900
      nday = idoy + yrlen*iyr
      idow = modulo (nday, 7) + 1
      ndim = daypm(imonth)
      ndiy = yrlen

      return
      end

      subroutine mkstmp (stamp, year, month, day, hour, min, sec)
!=======================================================================
!     make a 32 character time stamp from day,month,year,sec,min,hour
!     writes an i5 or i9 year in y/m/d, d/m/y or m/d/y format
!=======================================================================

      implicit none

      character(32) :: stamp

      integer year, month, day, hour, min, sec, y, yr

      yr = year

!     eliminate the leading digits of year if it is too large
      if (yr .gt. 0) then
        y = mod(yr,1000000000)
      else
        y = mod(yr,100000000)
      endif
      if (y .ne. yr) write (*,*)
     &  '==>Warning: in mkstmp, year=',yr,'is modified to year=',y
      if (y .gt. -10000 .and. y .lt. 100000) then
#if !defined O_stamp_ymd && !defined O_stamp_mdy
        write (stamp,'(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)') 'd/m/y='
     &,   day,'/',month,'/',y,',h:m:s=', hour,':', min,':', sec
      else
        write (stamp,'(a4,i2,a1,i2,a1,i9,a5,i2,a1,i2,a1,i2)') 'dmy='
     &,   day,'/',month,'/',y,',hms=', hour,':', min,':', sec
#elif defined O_stamp_mdy
        write (stamp,'(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)') 'm/d/y='
     &,   month,'/',day,'/',y,',h:m:s=', hour,':', min,':', sec
      else
        write (stamp,'(a4,i2,a1,i2,a1,i9,a5,i2,a1,i2,a1,i2)') 'mdy='
     &,   month,'/',day,'/',y,',hms=', hour,':', min,':', sec
#else
        write (stamp,'(a6,i5,a1,i2,a1,i2,a7,i2,a1,i2,a1,i2)') 'y/m/d='
     &,   y,'/',month,'/',day,',h:m:s=', hour,':', min,':', sec
      else
        write (stamp,'(a4,i9,a1,i2,a1,i2,a5,i2,a1,i2,a1,i2)') 'ymd='
     &,   y,'/',month,'/',day,',hms=', hour,':', min,':', sec
#endif
      endif

      return
      end

      subroutine rdstmp (stamp, year, month, day, hour, min, sec)
!=======================================================================
!     convert 32 character time stamp into day,month,year,sec,min,hour
!     reads stamp of form y/m/d, d/m/y or m/d/y (for many year lengths)
!=======================================================================

      implicit none

      character(32) :: stamp

      integer year, month, day, hour, min, sec

      character(1) :: skip1
      character(4) :: skip4
      character(5) :: skip5
      character(6) :: skip6
      character(7) :: skip7
      character(8) :: skip8

!     old MOM2 and MOM3 short year stamp (year=i4)
      if (stamp(18:18) .eq. ' ') then
        if (stamp(1:1) .eq. 'y') then
          read (stamp, '(a6,i4,a1,i2,a1,i2,a8,i2,a1,i2,a1,i2)')
     &      skip6, year, skip1, month, skip1, day, skip8, hour
     &,     skip1, min, skip1, sec
        elseif (stamp(1:1) .eq. 'd') then
          read (stamp, '(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &      skip6, day, skip1, month, skip1, year, skip8, hour
     &,     skip1, min, skip1, sec
        else
          read (stamp, '(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip8, hour
     &,     skip1, min, skip1, sec
        endif

!     new MOM3 short year stamp (year=i5)
      elseif (stamp(18:18) .eq. ',') then
        if (stamp(1:1) .eq. 'y') then
          read (stamp, '(a6,i5,a1,i2,a1,i2,a7,i2,a1,i2,a1,i2)')
     &      skip6, year, skip1, month, skip1, day, skip7, hour
     &,     skip1, min, skip1, sec
        elseif (stamp(1:1) .eq. 'd') then
          read (stamp, '(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &      skip6, day, skip1, month, skip1, year, skip7, hour
     &,     skip1, min, skip1, sec
        else
          read (stamp, '(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip7, hour
     &,     skip1, min, skip1, sec
        endif

!     old MOM3 long year stamp (year=i6)
      elseif (stamp(2:2) .eq. "/") then
        if (stamp(1:1) .eq. 'y') then
          read (stamp, '(a6,i6,a1,i2,a1,i2,a6,i2,a1,i2,a1,i2)')
     &      skip6, year, skip1, month, skip1, day, skip6, hour
     &,     skip1, min, skip1, sec
        elseif (stamp(1:1) .eq. 'd') then
          read (stamp, '(a6,i2,a1,i2,a1,i6,a6,i2,a1,i2,a1,i2)')
     &      skip6, day, skip1, month, skip1, year, skip6, hour
     &,     skip1, min, skip1, sec
        else
          read (stamp, '(a6,i2,a1,i2,a1,i6,a6,i2,a1,i2,a1,i2)')
     &      skip6, month, skip1, day, skip1, year, skip6, hour
     &,     skip1, min, skip1, sec
        endif

!     new MOM3 long year stamp (year=i9)
      else
        if (stamp(1:1) .eq. 'y') then
          read (stamp, '(a4,i9,a1,i2,a1,i2,a5,i2,a1,i2,a1,i2)')
     &      skip4, year, skip1, month, skip1, day, skip5, hour
     &,     skip1, min, skip1, sec
        elseif (stamp(1:1) .eq. 'd') then
          read (stamp, '(a4,i2,a1,i2,a1,i9,a5,i2,a1,i2,a1,i2)')
     &      skip4, day, skip1, month, skip1, year, skip5, hour
     &,     skip1, min, skip1, sec
        else
          read (stamp, '(a4,i2,a1,i2,a1,i9,a5,i2,a1,i2,a1,i2)')
     &      skip4, month, skip1, day, skip1, year, skip5, hour
     &,     skip1, min, skip1, sec
        endif
      endif

      return
      end

      function file_stamp (fname, stamp, suffix)
!=======================================================================
!     add day, month and year information from "stamp" to "fname"
!     writes year as 5 or 9 characters padded with zeros.
!=======================================================================

      implicit none

      character(*) :: fname, stamp, suffix
      character(3) :: dd(1:31)
      character(3) :: mm(1:12)
      character(10) :: yy
      character(120) :: file_stamp

      integer year, month, day, hour, min, sec, i

      save dd, mm

      data (dd(i),i=1,7)   /'.01','.02','.03','.04','.05','.06','.07'/
      data (dd(i),i=8,14)  /'.08','.09','.10','.11','.12','.13','.14'/
      data (dd(i),i=15,21) /'.15','.16','.17','.18','.19','.20','.21'/
      data (dd(i),i=22,28) /'.22','.23','.24','.25','.26','.27','.28'/
      data (dd(i),i=29,31) /'.29','.30','.31'/
      data (mm(i),i=1,7)   /'.01','.02','.03','.04','.05','.06','.07'/
      data (mm(i),i=8,12)  /'.08','.09','.10','.11','.12'/

      call rdstmp (stamp, year, month, day, hour, min, sec)
      if (year .ge. 0) yy = '.00000    '
      if (year .ge. 100000) yy = '.000000000'
      if (year .lt. 0) yy = '.-0000    '
      if (year .le. -10000) yy = '.-00000000'

!     use long year form (9 characters)
      if (abs(year) .ge. 100000000) then
        write (yy(2:10),'(i9)') abs(year)
      elseif (abs(year) .ge. 10000000) then
        write (yy(3:10),'(i8)') abs(year)
      elseif (abs(year) .ge. 1000000) then
        write (yy(4:10),'(i7)') abs(year)
      elseif (abs(year) .ge. 100000) then
        write (yy(5:10),'(i6)') abs(year)
      elseif (year .le. -10000) then
        write (yy(6:10),'(i5)') abs(year)
!     switch to short year form (5 characters)
      elseif (year .ge. 10000) then
        write (yy(2:6),'(i5)') abs(year)
      elseif (abs(year) .ge. 1000) then
        write (yy(3:6),'(i4)') abs(year)
      elseif (abs(year) .ge. 100) then
        write (yy(4:6),'(i3)') abs(year)
      elseif (abs(year) .ge. 10) then
        write (yy(5:6),'(i2)') abs(year)
      else
        write (yy(6:6),'(i1)') abs(year)
      endif

      file_stamp=trim(fname)//trim(yy)//mm(month)//dd(day)//trim(suffix)

      return
      end

      subroutine incstamp (stamp1, days, stamp2)
!=======================================================================
!     increment a time stamp by a real number of days
!     input:
!      stamp1 =  character date and time stamp
!      days   =  days to increment stamp (may be negative)
!     output:
!      stamp2 = stamp1 incremented by days
!=======================================================================

      implicit none

      character(*) :: stamp1, stamp2

      integer id, msec

      real days

      include "tmngr.h"
      include "calendar.h"

      call settimefromstamp (stamp1, itemptime)
      call settime2 (itemptime2, id(days*daylen), msec(days*daylen))
      call addtime (itemptime, itemptime2, itemptime)
      call expandtime2 (itemptime)
      stamp2 = tstamp(itemptime)
      return
      end

      subroutine inctime (index1, days, index2)
!=======================================================================
!     increment a full time by a real number of days
!     input:
!      index1 =  subscript of a full time
!      days   =  days to increment full time (may be negative)
!      index2 = subscript of the answer
!     output:
!      time arrays for index2 (year, month, day, etc) in tmngr.h
!=======================================================================

      implicit none

      integer id, msec, index1, index2

      real days

      include "tmngr.h"
      include "calendar.h"

      call settime2 (itemptime, id(days*daylen), msec(days*daylen))
      call addtime (index1, itemptime, index2)
      call expandtime2 (index2)
      return
      end

      subroutine settimefromstamp (stamp, index)
!=======================================================================
!     input:
!      stamp =  character dat and time stamp
!      index = subscript of the answer
!     output:
!      time arrays for index (year, month, day, etc) in tmngr.h
!=======================================================================

      implicit none

      character(*) :: stamp

      integer iyear, imonth, iday, ihour, imin, isec, index

      call rdstmp (stamp, iyear, imonth, iday, ihour, imin, isec)
      call setfulltime (index, iyear, imonth, iday, ihour, imin, isec)
      return
      end

      subroutine addtime (index1, index2, index)
!=======================================================================
!     add two times given in (integer day, nonneg integer ms) form
!     input:
!      index1 =  subscript of the first time into time arrays in tmngr.h
!      index2 =  subscript of the second time
!      index  = subscript of the answer
!     output:
!      iday(index) = integer day number of answer in tmngr.h
!      msday(index)= millisec fractional day of answer in tmngr.h
!=======================================================================

      implicit none

      integer mst, index1, index2, it, index

      include "tmngr.h"
      include "calendar.h"

      mst = msday(index1) + msday(index2)
      it  = iday(index1) + iday(index2)
      if (mst .ge. int(daylen*1000.)) then
        mst = mst - int(daylen*1000.)
        it  = it + 1
      endif
      iday(index) = it
      msday(index) = mst

      return
      end

      subroutine subtime (index1, index2, index)
!=======================================================================
!     subtract two times given in (integer day, nonneg integer ms) form
!     input:
!      index1 =  subscript of the first time into time arrays in tmngr.h
!      index2 =  subscript of the second time
!      index  = subscript of the answer
!     output:
!      iday(index) = integer day number of answer in tmngr.h
!      msday(index)= millisec fractional day of answer in tmngr.h
!=======================================================================

      implicit none

      integer mst, index1, index2, it, index

      include "tmngr.h"
      include "calendar.h"

      mst = msday(index1) - msday(index2)
      it  = iday(index1) - iday(index2)
      if (mst .lt. 0) then
        mst = mst + int(daylen*1000.)
        it = it - 1
      endif
      iday(index)  = it
      msday(index) = mst

      return
      end

      subroutine multime (n, index1, index)
!=======================================================================
!     multiply time (integer day, nonneg integer ms) by an integer
!     input:
!      n      = integer multiple
!      index1 = subscript of the time
!      index  = subscript of the answer
!     output:
!      iday(index) = integer day number of answer in tmngr.h
!      msday(index)= millisec fractional day of answer in tmngr.h
!=======================================================================

      implicit none

      integer n, index1, it, mst, index

      double precision dmst, dayl

      include "tmngr.h"
      include "calendar.h"

      dayl = daylen*1000.
      dmst = n*dble(msday(index1))
      it   = n*iday(index1)
      it   = it + int(dmst/dayl)
      mst  = nint(mod (dmst, dayl))
      iday(index)  = it
      msday(index) = mst

      return
      end

      subroutine ms2hms (ms, hour, min, sec)
!=======================================================================
!     converts fraction of a day in milliseconds to hour/min/sec
!=======================================================================

      implicit none

      integer ms, hour, min, sec, msrem

      hour  = ms / 3600000
      msrem = ms - hour * 3600000
      min   = msrem / 60000
      sec   = (msrem - min * 60000) / 1000

      return
      end

      function hms2ms (hour, min, sec)
!=======================================================================
!     converts fraction of a day in hour/min/sec to milliseconds
!=======================================================================

      implicit none

      integer hms2ms, hour, min, sec

      hms2ms = hour*3600000 + min*60000 + sec*1000

      return
      end

      function ifloor (realvalue)
!=======================================================================
!     f90 intrinsic "floor"
!     largest integer less than or equal to realvalue
!=======================================================================

      implicit none

      integer iflr, ifloor

      real realvalue

      iflr = int(realvalue) -1
      ifloor = iflr + int(realvalue - iflr)

      return
      end

      function id (realsec)
!=======================================================================
!     converts time in real seconds to the integer day part
!     see also "msec"
!=======================================================================

      implicit none

      integer id, ifloor

      real realsec

      include "calendar.h"

      id = ifloor(realsec/daylen)

      return
      end

      function msec (realsec)
!=======================================================================
!     extracts integer milliseconds from time in real seconds
!     see also "id"
!=======================================================================

      implicit none

      integer id, msec

      real realsec

      include "calendar.h"

      msec = nint (1000.0*(realsec-daylen*id(realsec)))

      return
      end

      function realsecs (index)
!=======================================================================
!     converts integer days and milliseconds to real seconds
!     input:
!      index = index of any time in tmngr.h
!=======================================================================

      implicit none

      integer index

      real realsecs

      include "tmngr.h"
      include "calendar.h"

      realsecs = daylen*iday(index) + msday(index)/1000.0

      return
      end

      function realdays (index)
!=======================================================================
!     converts integer days and milliseconds to real days
!     input:
!      index = index of any time in tmngr.h
!=======================================================================

      implicit none

      integer index

      real realdays

      include "tmngr.h"
      include "calendar.h"

      realdays = iday(index) + msday(index)/(daylen*1000.)

      return
      end

      function modulo (m, n)
!=======================================================================
!     Fortran 90 intrinsic
!     similar to mod, but remainder has sign of n.
!     always gives positive remainder when n .gt. 0, even if m .lt. 0
!=======================================================================

      implicit none

      integer m, modulo, n

      if (m .ge. 0) then
        modulo = mod(m,n)
      else
        modulo = abs(n) - mod(-m-1,n) - 1
      endif

      return
      end

      subroutine inittime
!=======================================================================
!     initialize time indices for gettime and getfulltime.
!=======================================================================

      implicit none

      include "tmngr.h"

      nextfulltime = 1
      nexttime = nfulltimes + 1

      return
      end

      subroutine getfulltime (index)
!=======================================================================
!     allocates and returns an index for a full time and increments the
!     next available full time index counter.

!     output:
!      index = index of a full time in tmngr.h
!=======================================================================

      implicit none

      integer index

      include "stdunits.h"
      include "tmngr.h"

      if (nextfulltime .gt. nfulltimes) then
        write (stdout, "(a,a,i4,a)") '     Too many full times.',
     &                    '  Increase nfulltimes = ', nfulltimes,
     &                    ' in tmngr.h'
        stop 'getfulltime'
      endif
      index = nextfulltime
      nextfulltime = nextfulltime + 1

      return
      end

      subroutine gettime (index)
!=======================================================================
!     returns an index for a time and increments the next available
!     time index counter.
!     output:
!      index = index of a short time (iday, msday only) in tmngr.h
!=======================================================================

      implicit none

      integer index

      include "stdunits.h"
      include "tmngr.h"

      if (nexttime .gt. ntimes) then
        write (stdout, "(a,a,i4,a)") '     Too many times.',
     &                    '  Increase ntimes = ', ntimes,
     &                    ' in tmngr.h'
        stop 'gettime'
      endif
      index = nexttime
      nexttime = nexttime + 1

      return
      end

      subroutine copyfulltime (index1, index)
!=======================================================================
!     copy all fields of a full time structure from one index to
!     another.
!     input:
!      index1 = index of a full time in tmngr.h
!      index  = index of a full time in tmngr.h
!     output:
!      full time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, index1

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index1, 'full')) stop 'copyfulltime'
      if (badindex (index,  'full')) stop 'copyfulltime'
      iday  (index) = iday  (index1)
      msday (index) = msday (index1)
      year  (index) = year  (index1)
      month (index) = month (index1)
      day   (index) = day   (index1)
      hour  (index) = hour  (index1)
      minute(index) = minute(index1)
      second(index) = second(index1)
      tstamp (index) = tstamp (index1)
      dayofyear (index) = dayofyear (index1)
      dayofweek (index) = dayofweek (index1)
      daysinmon (index) = daysinmon (index1)
      daysinyear(index) = daysinyear(index1)

      return
      end

      subroutine copytime (index1, index)
!=======================================================================
!     copy both fields of a short time structure from one index to
!     another.
!     input:
!      index1 = index of a time in tmngr.h
!      index  = index of a short time in tmngr.h
!     output:
!      short time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, index1

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index1, 'short')) stop 'copytime'
      if (badindex (index,  'short')) stop 'copytime'
      iday (index) = iday (index1)
      msday(index) = msday(index1)

      return
      end

      subroutine settime2 (index, it, mst)
!=======================================================================
!     set a time by placing given day and millisecond information in
!     a time structure array indexed by index.
!     input:
!      index  = index of a short time in tmngr.h
!      it     = integer day part
!      mst    = fractional part of day in millisec
!     output:
!      short time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, it, mst

      logical badindex

      include "stdunits.h"
      include "tmngr.h"
      include "calendar.h"

      if (mst .ge. int(daylen*1000.) .or. mst .lt. 0) then
        write (stdout, *) '     Invalid millisecond designation'
        stop 'settime2'
      endif
      if (badindex (index,  'short')) stop 'settime2'
      iday(index)  = it
      msday(index) = mst

      return
      end

      subroutine settime3 (index, rd)
!=======================================================================
!     set a time structure given a time in real days
!     input:
!      index  = index of a short time in tmngr.h
!      rd     = time in real days
!     output:
!      short time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, id, msec

      logical badindex

      real rd

      include "stdunits.h"
      include "tmngr.h"
      include "calendar.h"

      if (badindex (index,  'short')) stop 'settime3'

      iday(index)  = id(rd*daylen)
      msday(index) = msec(rd*daylen)

      return
      end

      subroutine setfulltime2 (index, it, mst)
!=======================================================================
!     set a full time by placing given day and millisecond information
!     in a time structure array indexed by index and expanding it
!     using the d2ymd conversion and making a stamp.
!     input:
!      index  = index of a full time in tmngr.h
!      it     = integer day part
!      mst    = fractional part of day in millisec
!     output:
!      full time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, it, mst

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index, 'full')) stop 'setfulltime2'
      call settime2 (index, it, mst)
      call d2ymd (it, year(index), month(index), day(index),
     &            dayofyear(index), dayofweek(index), daysinmon(index),
     &            daysinyear(index))
      call ms2hms(mst, hour(index), minute(index), second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))

      return
      end

      subroutine setfulltime (index, yr, mon, dy, hr, min, sec)
!=======================================================================
!     set a full time by specifying the y/m/d, h/m/s and expanding it
!     using the ymd2d conversion and making a stamp.
!     input:
!      index  = index of a full time in tmngr.h
!      yr     = year to set
!      mon    = month to set
!      dy     = day to set
!      hr     = hour to set
!      min    = minute to set
!      sec    = second to set
!     output:
!      full time fields set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, yr, mon, dy, hr, min, sec, hms2ms

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index, 'full')) stop 'setfulltime'
      year  (index) = yr
      month (index) = mon
      day   (index) = dy
      hour  (index) = hr
      minute(index) = min
      second(index) = sec
      call ymd2d (yr, mon, dy, iday(index),
     &            dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))
      msday(index) = hms2ms(hr, min, sec)
      call mkstmp (tstamp(index), yr, mon, dy, hr, min, sec)

      return
      end

      subroutine expandtime (index)
!=======================================================================
!     expand a full time specified by the y/m/d, h/m/s
!     using the ymd2d conversion and making a stamp.
!     input:
!      index  = index of a full time in tmngr.h
!      y/m/d, h/m/s must be set for time index
!     output:
!      full time fields (iday,msday, tstamp) set for index in tmngr.h
!=======================================================================

      implicit none

      integer index, hms2ms

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index, 'full')) stop 'expandtime'
      call ymd2d (year(index), month(index), day(index),
     &            iday(index), dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))
      msday(index) = hms2ms(hour(index), minute(index), second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))

      return
      end

      subroutine expandtime2 (index)
!=======================================================================
!     expand a full time specified by absolute day and millisecond
!     information using d2ymd conversion and making a stamp.
!     input:
!      index  = index of a full time in tmngr.h
!      iday, msday must be set for time index
!     output:
!      full time fields (year,month,day,hour,minute,second,tstamp)
!      set for index in tmngr.h
!=======================================================================

      implicit none

      integer index

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index, 'full')) stop 'expandtime2'
      call d2ymd (iday(index), year(index), month(index),
     &            day(index), dayofyear(index), dayofweek(index),
     &            daysinmon(index), daysinyear(index))
      call ms2hms(msday(index), hour(index), minute(index),
     &            second(index))
      call mkstmp (tstamp(index), year(index), month(index), day(index),
     &             hour(index), minute(index), second(index))

      return
      end

      function timeless (index1, index2)
!=======================================================================
!     compare times referenced by index1 and index2.  timeless is
!     true if time1 is less than time2.
!=======================================================================

      implicit none

      integer index1, index2

      logical timeless, badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index1, 'short')) stop 'timeless'
      if (badindex (index2, 'short')) stop 'timeless'
      if (iday(index1) .lt. iday(index2)) then
        timeless = .true.
      elseif ((iday (index1) .eq. iday (index2)) .and.
     &        (msday(index1) .lt. msday(index2))) then
        timeless = .true.
      else
        timeless = .false.
      endif

      return
      end

      function timeequal (index1, index2)
!=======================================================================
!     compare times referenced by index1 and index2.  timeequal is
!     true if time1 is equal to  time2.
!=======================================================================

      implicit none

      integer index1, index2

      logical timeequal, badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index1, 'short')) stop 'timeequal'
      if (badindex (index2, 'short')) stop 'timeequal'
      if (iday(index1) .eq. iday(index2) .and.
     &    msday(index1) .eq. msday(index2)) then
        timeequal = .true.
      else
        timeequal = .false.
      endif

      return
      end

      function timemore (index1, index2)
!=======================================================================
!     compare times referenced by index1 and index2.  timemore is
!     true if time1 is more than time2.
!=======================================================================

      implicit none

      integer index1, index2

      logical timemore, badindex

      include "stdunits.h"
      include "tmngr.h"

      if (badindex (index1, 'short')) stop 'timemore'
      if (badindex (index2, 'short')) stop 'timemore'
      if (iday(index1) .gt. iday(index2)) then
        timemore = .true.
      elseif ((iday (index1) .eq. iday (index2)) .and.
     &        (msday(index1) .gt. msday(index2))) then
        timemore = .true.
      else
        timemore = .false.
      endif

      return
      end

      function badindex (index, timetype)
!=======================================================================
!     check to see if a given index is in bounds in tmngr.h
!     type = 'full'  means a full time
!     type = 'short' means a short time
!=======================================================================

      implicit none

      character(*) :: timetype

      integer index

      logical badindex

      include "stdunits.h"
      include "tmngr.h"

      if (timetype .eq. 'full') then
        if (index .lt. 1 .or. index .gt. nfulltimes) then
          badindex = .true.
          print *, '     Invalid full time index = ',index
        else
          badindex = .false.
        endif
      else
        if (index .lt. 1 .or. index .gt. ntimes) then
          badindex = .true.
          print *, '     Invalid time index = ',index
        else
          badindex = .false.
        endif
      endif

      return
      end
