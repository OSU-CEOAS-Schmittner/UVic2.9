      subroutine timeinterpi (nrec, stamp1, aprec, tdrec, isbcstart
     &,                       period)

!=======================================================================
!     initializes time centre of each data record based on the time
!     stamp and average period.
!=======================================================================

      implicit none

      character(*) :: stamp1

      integer nrec, isbcend, m, isbcstart, isbcyear, isbcmon, isbcday
      integer isbchour, isbcmin, isbcsec

      logical period

      real sum

      include "tmngr.h"

      real aprec(nrec), tdrec(nrec)

      data isbcend /0/
      save isbcend

!     define each climatological data record to be at the centre of the
!     month starting with month "isbcmonth"

      sum = 0.0
      do m=1,nrec
        sum      = sum + aprec(m)
        tdrec(m) = sum - 0.5*aprec(m)
      enddo

!     calculate time at start of first record: "isbcstart"

      if (isbcend .eq. 0) call getfulltime (isbcend)
      call getfulltime (isbcstart)

      call rdstmp (stamp1, isbcyear, isbcmon, isbcday, isbchour
     &,            isbcmin, isbcsec)
      call setfulltime (isbcend, isbcyear, isbcmon, isbcday, isbchour
     &,            isbcmin, isbcsec)
      call inctime (isbcend, -aprec(1), isbcstart)

!     check integrity of data record times. also when using datasets
!     as periodic, add 0.2425 days to febuary and adjust subsequent
!     months to account for this change when using the time manager
!     with a leap year calendar.

      if (period) call checkinterp (nrec, tdrec, aprec)
      return
      end

      subroutine timeinterp (tm, n, tdrec, aprec, ndr, period, method
     &,                     ia, ib, wb, change, inext, iprev)

!=======================================================================

!     time interpolator ... constructs indices & weight needed for
!     linearly interpolating data defined at arbitrary time intervals
!     (midpoints of years, months, days or  random intervals) to
!     the time of the current model time step.

!     inputs:

!     tm     = the time at which the data is desired (units of "tdrec")

!     tdrec  = the times at which the data records in the dataset are
!              defined. times must be monotonically increasing and are
!              assumed to be at the centres of the averaging periods.
!              (eg: the centres of the months if using monthly averaged
!               climatology. units are arbitrary)

!     aprec  = array of averaging periods for the data records
!              (eg: the number of days per month)

!     ndr    = number of data records in the dataset. (eg: 12 if using
!              monthly climatology)

!     period = (true,false) if the dataset is to be treated as
!              (perodic, not periodic). if periodic, then the model
!               time is always mapped into the dataset. if not, then
!               record 1 is used for all model time before the
!               beginning of the dataset and record "ndr" is used for
!               all model time after the end of the dataset.

!     method = interpolation scheme desired.  (0..3)
!                0 = no interpolation; the average value is used
!                    for all times in the entire averaging period.
!                    (preserves the integral over averaging periods,
!                    but is discontinuous at period boundaries.)
!                1 = linear interpolation between the middles of
!                    two adjacent averaging periods.
!                    (continuous but does not preserve integral for
!                    unequal periods.)
!                2 = equal linear interpolation.  Assumes that the
!                    value on the boundary between two adjacent
!                    averaging periods is the unweighted average of
!                    the two average values.  Linearly interpolates
!                    between the midperiod and period boundary.
!                    (continuous but does not preserve integral for
!                    unequal periods.)
!                3 = equal area (midperiod to midperiod) interpolation
!                    chooses a value for the boundary between two
!                    adjacent periods such that linear interpolation
!                    between the two midperiods and this value will
!                    preserve the integral midperiod to midperiod.
!                Note that methods 1,2, and 3 are equivalent if
!                all periods lengths are equal.

!     n      = a number denoting which dataset is being interpolated
!              (each dataset should be referenced by a unique number
!               starting with 1 for the 1st, 2 for the 2nd, ...etc)

!     outputs:

!     ia     = index for pointing to the next data record which will be
!              reached by the model. (eg: ahead of the model. "ia" would
!              be 3 if "tm" was beyond the  middle of {but still within}
!              february)
!     ib     = index for pointing to the data record which was just
!              passed by the model. (eg: behind the model. "ib" would
!              be 2 if "tm" was beyond the middle of {but still within}
!              february)
!     inext  = index to memory buffer containing data from "ia"
!     iprev  = index to memory buffer containing data from "ib"
!     wb     = interpolation weight for defining data at "tm"
!              schematically the interpolation is defined by:

!              data(iprev) <== disk data "ib"
!              data(inext) <== disk data "ia"
!              data(tm) = wb*data(iprev) + (1-wb)*data(inext)

!     change = logical for sensing when "ia" and "ib" change.
!              when change = T then it is time to read the disk
!              and update "inext" and "iprev"
!=======================================================================

      implicit none

      integer maxsets, iflag
      parameter (maxsets=15, iflag=-99999)

      integer ndr, n, method, ib, indp, ia, ic, io, itemp, iprev, inext
      integer iaold(maxsets), imethod(maxsets)

      logical change, period

      real dstart, dend, dlen, tm, d, f, time, startaft, dtmid, dtbnd
      real dtomid, wc, wb, tdrec(ndr), aprec(ndr)

      data iaold /maxsets*iflag/
      save iaold, imethod

!-----------------------------------------------------------------------
!     statement function
!-----------------------------------------------------------------------

      if (n .gt. maxsets) then
        write (*,'(a,i10,a,i10)') 'Error: n=', n, ' maxsets=',maxsets
        stop '=>timeinterp'
      endif

      if (iaold(n) .eq. iflag) then
        write (*,'(/1x,a,i2,a,i3/)')
     &      'Assigning interpolation method ',method, ' to dataset # ',n
        imethod(n) = method
      endif

      if (method .ne. imethod(n)) then
        write (*,'(/a,i2,a,i3/a,i2,a/)')
     &   'Error: trying to use method ',method, ' on dataset # ',n
     &,  'originally, method ',imethod(n),' was used in timeinterp'
        stop
      endif

      if (period) then

!       define the position of the dataset in time

        dstart = tdrec(1) - 0.5*aprec(1)
        dend   = tdrec(ndr) + 0.5*aprec(ndr)
        dlen   = dend - dstart

!       map the model time into the dataset assuming dataset periodicity

        if (tm .lt. dstart) then
          d = dstart - tm
          f = d/dlen - int(d/dlen)
          time = dend - f*dlen
        elseif (tm .gt. dend) then
          d = tm - dend
          f = d/dlen - int(d/dlen)
          time = dstart + f*dlen
        else
          time = tm
        endif
      else

!       define the position of the dataset in time. no periodicity

        dstart = tdrec(1)
        dend   = tdrec(ndr)
        dlen   = dend - dstart

!       map the model time into the dataset. assume data is constant
!       before the beginning and after the end of the dataset

        if (tm .lt. dstart) then
          time = dstart
        elseif (tm .gt. dend) then
          time = dend
        else
          time = tm
        endif
      endif

!     calculate record pointers and weighting for interpolation of
!     dataset records to the model time step.

      ib = indp (time, tdrec, ndr)
      if (tdrec(ib) .gt. time) ib = ib - 1
      if (period) then
        ia = mod(ib, ndr) + 1
        if (ib .lt. 1) ib = ndr
      else
        ia = ib + 1
        if (ia .gt. ndr) ia = ib
        if (ib .lt. 1)   ib = ia
      endif

!     find whether "time" is closer to midpoint of record "ia" or ib"
!     ic is the index of the closest midpoint
!     io is the index of the other midpoint

      startaft = tdrec(ia) - 0.5*aprec(ia)
      if (time .ge. startaft .and. time .le. tdrec(ia)) then
        ic = ia
        io = ib
      else
        ic = ib
        io = ia
      endif

!     dtmid = distance from "time" to midpoint of closer record
!     dtbnd = distance from "time" to boundary of closer record
!     dtomid = distance from "time" to midpoint of other record

      dtmid  = abs(time - tdrec(ic))
      dtbnd  = 0.5*aprec(ic) - dtmid
      dtomid = 0.5*aprec(io) + dtbnd

!-----------------------------------------------------------------------
!     3) equal area (midperiod to midperiod) interpolation formula
!-----------------------------------------------------------------------

      if (method .eq. 3) then
        wc = 2.0*dtbnd/aprec(ic) + 2.0*dtmid/(aprec(ic) + aprec(io))

!-----------------------------------------------------------------------
!     2) equal linear interpolation
!             value on period boundary assumed to be average of values
!             on the two adjacent periods.
!-----------------------------------------------------------------------

      elseif (method .eq. 2) then
        wc = (2.0*dtbnd + dtmid)/aprec(ic)

!-----------------------------------------------------------------------
!     1) linear interpolation
!-----------------------------------------------------------------------

      elseif (method .eq. 1) then
        wc = dtomid/(dtmid + dtomid)

!-----------------------------------------------------------------------
!     0) no interpolation
!-----------------------------------------------------------------------

      elseif (method .eq. 0) then
        wc = 1.0
      else

!-----------------------------------------------------------------------
!     anything else is not allowed for (unless you want to add one!)
!-----------------------------------------------------------------------

        print *,'=>Error: method = ',method,' not allowed in timeinterp'
        stop
      endif

      if (ib .eq. ic) then
        wb = wc
      else
        wb = 1.0 - wc
      endif
      if (wc .lt. 0.0 .or. wc .gt. 1.0) then
        print *,' ic=',ic,' io=',io, ' dtmid=',dtmid,' dtbnd=',dtbnd
     &,' dtomid=',dtomid, ' time=',time, ' ia=',ia,' ib=',ib
     &, ' wc=',wc
        print *,' =>Error: bad interpolation weight in timeinterp'
        stop
      endif

!     refresh pointers to memory buffers when reading disk data

      if (iaold(n) .ne. ia) then
        change = .true.
        itemp = iprev
        iprev = inext
        inext = itemp
      else
        change = .false.
      endif
      iaold(n) = ia

      return
      end

      subroutine checkinterp (ntdrec, tdrec, aprec)

!=======================================================================
!     check for consistency between interpolation period centres "tdrec"
!     and period lengths "aprec".
!     adjust tdrec and aprec for leap years
!     check for and compensate for some mismatches between data and
!     calendar
!=======================================================================

      implicit none

      integer ntdrec, m

      logical febdone, monthly, error

      real dlen, sum, time

      include "calendar.h"
      include "stdunits.h"

      real tdrec(ntdrec), aprec(ntdrec)

!     test for consistency of tdrec and aprec times

      monthly = .true.
      error   = .false.
      sum     = 0.5*aprec(1)
      do m=2,ntdrec
        sum = sum + 0.5*(aprec(m) + aprec(m-1))
        if (abs(tdrec(m) - sum) .gt. 0.01*tdrec(m)) then
          error = .true.
          write (stdout,*) 'Error in time interpolation data'
          write (stdout,*) 'Date for middle of record ',m
     &,                   ' is not centred'
        endif
        if (.not.(28.0 .le. aprec(m) .and. aprec(m) .le. 32)) then
          monthly = .false.
        endif
      enddo

      if (error) then
        write (stdout,*) 'STOP in checkinterp'
!        stop
      endif

      dlen = (tdrec(ntdrec) + 0.5*aprec(ntdrec)) -
     &       (tdrec(1) - 0.5*aprec(1))

!     if using leap years, add 1/4 day to feburary (or last record in
!     feb if data is other than monthly. eg: daily)

      if (.not.eqyear) then
        if (mod(dlen, real(yrlen)) .lt. 0.01) then

!         calendar has leap years but data does not, add 1/4 day to
!         feburary (or last record in feb if data is other than monthly.
!         eg: daily)

          write (stdout, '(/,a,a)')
     &           'Checkinterp: Modifying equal year interpolation'
     &,          ' data for use with leap year calendar'
          febdone = .false.
          time    = 0.0
          do m=1,ntdrec
            time = time + aprec(m)
            if (time .ge. yrlen) then
              time = time - yrlen
              febdone = .false.
            endif
            if (time .ge. msum(3)) then
              if (.not. febdone) then
                aprec(m)   = aprec(m) + 0.2425
                write (stdout, '(a,i4)')
     &                'Checkinterp: Adding 0.2425 days to record ',m
                febdone = .true.
              endif
            endif
          enddo
          sum = tdrec(1) - 0.5*aprec(1)
          do m=1,ntdrec
            sum = sum + aprec(m)
            tdrec(m) = sum - 0.5*aprec(m)
!            print *,' m=',m,' tdrec=',tdrec(m), ' aprec=',aprec(m)
          enddo

        elseif (mod(dlen, real(yrlen)) - 0.2425*nint(dlen/yrlen)
     &          .le. 0.01) then

!         calendar has leap years and data is leap year corrected by adding
!         0.2425 days per year.  interpolation data is consistent.

        else

!         calendar has leap years but data is neither leap year
!         compensated nor an exact number of years, it is not clear
!         what user wants.

          write (stdout,*) 'Problem in checkinterp'
          write (stdout,*) 'Calendar uses leap years, but interpolation'
     &,     ' data is neither an integer number of years or leap year'
     &,     ' corrected by adding 0.2425 days per year.'
          stop
        endif

      else
        if (mod(dlen, real(yrlen)) .lt. 0.01) then

!         calendar uses equal years and data is an integral number of these
!         years.  interpolation data is consistent.

        elseif (mod(dlen, real(yrlen)) - 0.2425*nint(dlen/yrlen)
     &          .le. 0.01) then

!         calendar uses equal years, but data is leap year corrected.
!         subtract 1/4 day from feburary (or last record in feb if data is other
!         than monthly. eg: daily)

          write (stdout, '(/,a,a)')
     &      'Checkinterp: Modifying leap year corrected'
     &,     ' interpolation data for use with equal years'
          febdone = .false.
          time    = 0.0
          do m=1,ntdrec
            time = time + aprec(m)
            if (time .ge. yrlen + 0.2425) then
              time = time - yrlen - 0.2425
              febdone = .false.
            endif
            if (time .ge. msum(3)+0.2425) then
              if (.not. febdone) then
                aprec(m)   = aprec(m) - 0.2425
                write (stdout, '(a,i4)')
     &           'Checkinterp: Subtracting 0.2425 days from record ',m
                febdone = .true.
              endif
            endif
          enddo
          sum = tdrec(1) - 0.5*aprec(1)
          do m=1,ntdrec
            sum = sum + aprec(m)
            tdrec(m) = sum - 0.5*aprec(m)
!            print *,' m=',m,' tdrec=',tdrec(m), ' aprec=',aprec(m)
          enddo

        else

!         calendar has equal years but data is neither leap year
!         compensated nor an exact number of years, it is not clear
!         what user wants.

          write (stdout,*) 'Problem in checkinterp'
          write (stdout,*) 'Calendar uses equal years, but'
     &,     ' interpolation data is neither an integer number of years'
     &,     '  or leap year corrected by adding 0.2425 days per year.'
          stop

        endif

      endif
      return
      end
