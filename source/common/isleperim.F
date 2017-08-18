      subroutine isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,                    imt, jmt, km, mnisle, maxipp
     &,                    xu, yu, zw)

!=======================================================================
!         Island and Island Perimeter Mapping Routines

!          The main computational subroutine, expand, uses a "floodfill"
!          algorithm to expand one previously unmarked land
!          point to its entire connected land mass and its perimeter
!          ocean points.   Diagonally adjacent land points are
!          considered connected.  Perimeter "collisions" (i.e.,
!          ocean points that are adjacent to two unconnected
!          land masses) are detected and error messages generated.

!          Perimeter collisions must be removed because the stream
!          function cannot be unambiguously defined on such cells.
!          The subroutine isleperim will not return until all perimeter
!          collisions are removed.

!          The subroutine expand uses a queue of size maxq of
!          coordinate pairs of candidate points.  Suggested
!          size for maxq is 4*(imt+jmt).  Queue overflow stops
!          execution with a message to increase the size of maxq.
!          Similarly a map with more that maxipp island perimeter
!          points or more than mnisle land masses stops execution
!          with an appropriate error message.

!          Computes map of land masses and island perimeters

!          Input:
!                  kmt = array of depths.  (0 for land) (>0 for ocean)
!          Outputs:
!                  map = map of land masses and their ocean perimeters
!                           mid-ocean cells are labelled 0
!                           land masses are labelled 1, 2, 3, ...,
!                           their perimeter ocean cells -1, -2, -3, ...,
!                  iperm = i coordinates of perimeter points
!                  jperm = j coordinates of perimeter points
!                  iofs = offset of each land mass in iperm, jperm
!                  nippts = number of island perimeter points by isle
!                  nisle = number of land masses
!          Array size inputs:
!                  imt = east/west array extents
!                  jmt = north/south array extents
!                  mnisle = maximum number of land masses
!                  maxipp = maximum number of island perimeter points

!          Arguments used in editing kmt field
!                  kmt_opt= selectable options for kmt changes
!                  kmt_changes = changes to kmt field
!                  nchanges = number of changes to kmt field
!                  i_del_kmt = i/o unit for changes in kmt
!                  xu = longitude (degrees) at "u" points
!                  yu = latitude (degrees) at "u" points
!                  zw = depth at bottom of cells
!=======================================================================

      implicit none

      integer imt, jmt, maxipp, mnisle, km, maxq, land, kmt_land
      integer kmt_ocean, i, j, maxqsize, label, nerror, jnorth
      integer iwest, ieast, nisle, isle, linewidth

      common /qsize/ maxqsize

      integer kmt(imt,jmt), map(imt,jmt), iperm(maxipp), jperm(maxipp)
      integer nippts(mnisle), iofs(mnisle)

      include "isleperim.h"

      real xu(imt), yu(jmt), zw(km)

      parameter (maxq=10000)

      integer iq(maxq)
      integer jq(maxq)
      integer qfront, qback
      integer ocean

      parameter (land=1, ocean=0)
      parameter (kmt_land=0, kmt_ocean=1)

      print '(/,a,/)','Finding perimeters of all land masses'

!     initialize number of changes to kmt

      nchanges = 0

    1 continue

!-----------------------------------------------------------------------
!     copy kmt to map changing notation
!     initially, 0 means ocean and 1 means unassigned land in map
!     as land masses are found, they are labelled 2, 3, 4, ...,
!     and their perimeter ocean cells -2, -3, -4, ...,
!     when no land points remain unassigned, land mass numbers are
!     reduced by 1 and their perimeter ocean points relabelled accordingly
!-----------------------------------------------------------------------

      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .gt. 0) then
            map(i,j) = ocean
          else
            map(i,j) = land
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     find unassigned land points and expand them to continents
!-----------------------------------------------------------------------

      maxqsize = 0
      call qinit (iq, jq, qfront, qback)
      label = 2
      iofs(label) = 0
      nippts(label) = 0
      nerror = 0
#if defined O_symmetry
      jnorth = jmt-1
#else
      jnorth = jmt
#endif
#if defined O_cyclic
      iwest = 2
      ieast = imt-1
#else
      iwest = 1
      ieast = imt
#endif
      do j=jnorth,1,-1
        do i=iwest,ieast
          if (map(i,j) .eq. land) then
            call qpush (i, j, iq, jq, qfront, qback)
            call expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt
     &,                  xu, yu, zw)
            print '(a,i2,a,i4)',
     &        'number of island perimeter points: nippts(',label-1,')=',
     &         nippts(label)
            label = label + 1
            if (label .gt. mnisle) then
              print '(a,i3,a)','ERROR==> mnisle=',mnisle,' is too small'
              stop '==> expand'
            endif
            iofs(label) = iofs(label-1) + nippts(label-1)
            nippts(label) = 0
          endif
        enddo
      enddo
      nisle = label - 1

!-----------------------------------------------------------------------
!     relabel land masses and their ocean perimeters
!------------------------------------------------------------------------

      do i=iwest,ieast
        do j=1,jnorth
          if (map(i,j) .ne. 0) then
            map(i,j) = map(i,j) - sign(1, map(i,j))
          endif
        enddo
      enddo
      do isle=2,nisle
        iofs(isle-1) = iofs(isle)
        nippts(isle-1) = nippts(isle)
      enddo
      nisle = nisle - 1

#if defined O_symmetry
      do i=iwest,ieast
        map(i,jmt) = map(i,jmt-1)
      enddo
#endif
#if defined O_cyclic
      do j=1,jmt
        map(1,j) = map(imt-1,j)
        map(imt,j) = map(2,j)
      enddo
#endif

!      print '(/a)',
!     &       'The following changes fix "PERIMETER VIOLATIONS"'

      call enter_kmt_changes()

      if (nerror .gt. 0) then
        print *,' Island perimeter statistics:'
        print *,'maximum queue size was ',maxqsize
        print *,'number of land masses is ', nisle
        print *,'number of island perimeter points is ',
     &          nippts(nisle) + iofs(nisle)
        print *, ' '
        print *
     &, '==>Remapping land masses to see if PERIMETER VIOLATIONS remain'
        print *, ' '
        goto 1
      endif

      print *,' Island perimeter statistics:'
      print *,'maximum queue size was ',maxqsize
      print *,'number of land masses is ', nisle
      print *,'number of island perimeter points is ',
     &        nippts(nisle) + iofs(nisle)

      return
      end

      subroutine showmap (map, imt, jmt, linewidth)

      implicit none

      integer iline, imt, iremain, istart, isweep, j, jmt, linewidth
      integer map(imt,jmt), l, i, mmm

      linewidth = 125
      print '(/,132a)',(' ',l=1,5+min(linewidth,imt)/2-13)
     &,     'Land Masses and Perimeters'
      istart = 0
      iremain = imt
      do isweep=1,imt/linewidth + 1
        iline = min(iremain, linewidth)
        iremain = iremain - iline
        if (iline .gt. 0) then
          print *, ' '
          print '(t6,32i5)', (istart+i+4,i=1,iline,5)
          do j=jmt,1,-1
            print '(i4,t6,160i1)', j,(mmm(map(istart+i,j)),i=1,iline)
          enddo
          print '(t6,32i5)', (istart+i+4,i=1,iline,5)
          istart = istart + iline
        endif
      enddo
      print *, ' '
      return
      end

      function mmm(m)

      implicit none

      integer mmm, m

      mmm = 0
      if (m .gt. 0) mmm = mod(m,10)

      return
      end

      subroutine expand (map, label, iq, jq, qfront, qback, nerror
     &,                  iperm, jperm, iofs, nippts
     &,                  imt, jmt, km, mnisle, maxipp, kmt
     &,                  xu, yu, zw)

!-----------------------------------------------------------------------
!          The subroutine expand uses a "flood fill" algorithm
!          to expand one previously unmarked land
!          point to its entire connected land mass and its perimeter
!          ocean points.   Diagonally adjacent land points are
!          considered connected.  Perimeter "collisions" (i.e.,
!          ocean points that are adjacent to two unconnected
!          land masses) are detected and error messages generated.

!          The subroutine expand uses a queue of size maxq of
!          coordinate pairs of candidate points.  Suggested
!          size for maxq is 4*(imt+jmt).  Queue overflow stops
!          execution with a message to increase the size of maxq.
!          Similarly a map with more that maxipp island perimeter
!          points or more than mnisle land masses stops execution
!          with an appropriate error message.
!-----------------------------------------------------------------------

      implicit none

      integer imt ,jmt ,maxipp ,mnisle ,km ,maxq ,land ,mnisle2
      integer label ,isle ,i ,j ,jn ,ie ,js ,iw ,nerror ,i1 ,j1 ,n

      integer map(imt,jmt), kmt(imt,jmt)

      integer iperm(maxipp)
      integer jperm(maxipp)
      integer nippts(mnisle)
      integer iofs(mnisle)

      include "isleperim.h"
      real xu(imt), yu(jmt), zw(km)
      character(32) :: problem

      parameter (maxq=10000)
      integer iq(maxq)
      integer jq(maxq)
      integer qfront, qback
      logical qempty

      integer offmap, ocean
      parameter (offmap = -1)
      parameter (land = 1, ocean = 0)

      parameter (mnisle2=100)
      logical bridge_to(1:mnisle2)

      print '(a,i3)', 'Exploring land mass ',label-1

      if (mnisle2 .lt. mnisle) then
        print '(a,i4,a)',
     & 'ERROR:  change parameter (mnisle2=',mnisle,') in isleperim.F'
        stop '==>isleperim'
      endif

      do isle=1,mnisle
        bridge_to(isle) = .false.
      enddo

!-----------------------------------------------------------------------
!     main loop:
!        Pop a candidate point off the queue and process it.
!-----------------------------------------------------------------------

 1000 continue

      if (qempty (iq, jq, qfront, qback)) then
        call qinit (iq, jq, qfront, qback)
        return
      else
        call qpop (i, j, iq, jq, qfront, qback)

!       case: (i,j) is off the map
        if (i .eq. offmap .or. j .eq. offmap) then
          goto 1000

!       case: map(i,j) is already labelled for this land mass
        elseif (map(i,j) .eq. label) then
          goto 1000

!       case: map(i,j) is an ocean perimeter point of this land mass
        elseif (map(i,j) .eq. -label) then
          goto 1000

!       case: map(i,j) is an unassigned land point
        elseif (map(i,j) .eq. land) then
          map(i,j) = label
!         print *, 'labeling ',i,j,' as ',label
          call qpush (i,         jn(j,jmt), iq, jq, qfront, qback)
          call qpush (ie(i,imt), jn(j,jmt), iq, jq, qfront, qback)
          call qpush (ie(i,imt), j,         iq, jq, qfront, qback)
          call qpush (ie(i,imt), js(j,jmt), iq, jq, qfront, qback)
          call qpush (i,         js(j,jmt), iq, jq, qfront, qback)
          call qpush (iw(i,imt), js(j,jmt), iq, jq, qfront, qback)
          call qpush (iw(i,imt), j,         iq, jq, qfront, qback)
          call qpush (iw(i,imt), jn(j,jmt), iq, jq, qfront, qback)
          goto 1000

!       case: map(i,j) is an ocean point adjacent to this land mass
        elseif (map(i,j) .eq. ocean .or. map(i,j) .lt. 0) then

!         subcase: map(i,j) is a perimeter ocean point of another mass
          if (map(i,j) .lt. 0) then
            nerror = nerror + 1
            print '(a,a,i3,a,i3,a,a,i3,a,i3)',
     &            'PERIMETER VIOLATION==> ',
     &            'map(',i,',',j,') is in the perimeter of both ',
     &            'land masses ', -map(i,j)-1, ' and ', label-1
!           if we just quit processing this point here, problem points
!           will be flagged several times.
!           if we relabel them, then they are only flagged once, but
!           appear in both island perimeters, which causes problems in
!           island integrals.  current choice is quit processing.

!           only fill first common perimeter point detected.
!           after the first land bridge is built, subsequent collisions
!           are not problems.

            if (.not. bridge_to(-map(i,j)-1)) then
              call clear_kmt_options ()

!             OK to set kmt to 0, but don`t change map because that would
!             require complete relabelling of land masses

!             option 1: fill common perimeter point to make land bridge

              call kmt_option (1, i, j, kmt(i,j), 0, kmt)

!             option 2: user selected interactive changes

!             option 3: change nearby land points on older land mass to ocean.
!                       we do not want to change points on the newer land
!                       mass because it may be incompletely explored at the
!                       current moment and some conflicting land points may
!                       be missed.

              do i1=-1,1
                do j1=-1,1
                  if (map(i+i1,j+j1) .eq. -map(i,j) .and.
     &                kmt(i+i1,j+j1) .eq. 0) then
                    call kmt_option (3, i+i1, j+j1
     &,                              kmt(i+i1,j+j1), kmt(i,j), kmt)
                  endif
                enddo
              enddo

              problem = 'perim'
              call select_option ('perim'
     &,                            i, j, kmt, xu, yu, zw)

!             see if option selected builds a land bridge

              do n=1,nchanges
                if (kmt_changes(n,1) .eq.i .and.
     &              kmt_changes(n,2) .eq.j .and.
     &              kmt_changes(n,4) .eq.0) then
                  bridge_to(-map(i,j)-1) = .true.
                endif
              enddo
            endif

            goto 1000
          endif

!         case: map(i,j) is a ocean point--label it for current mass
          map(i,j) = -label
          nippts(label) = nippts(label) + 1
!         print *, 'iofs(label)=',iofs(label)
!         print *, 'nippts(label)=',nippts(label)
          if (iofs(label) + nippts(label) .gt. maxipp) then
            print *, 'ERROR==>  maxipp=',maxipp,' is not large enough'
            stop '==>expand'
          endif
          iperm(iofs(label) + nippts(label)) = i
          jperm(iofs(label) + nippts(label)) = j
          goto 1000

!       case: map(i,j) is probably labelled for another land mass
!       ************* this case should not happen **************
        else
          nerror = nerror + 1
          print '(a,a,i3,a,i3,a,a,i3,a,i3)',
     &          'ERROR ==>  ',
     &          'map(',i,',',j,') is labelled for both ',
     &          'land masses ', map(i,j)-1,' and ',label-1
        endif
        goto 1000

      endif
      return
      end

      subroutine qinit (iq, jq, qfront, qback)

      implicit none

      integer maxq
      parameter (maxq=10000)
      integer qfront, qback, iq(maxq), jq(maxq)

      qfront = 1
      qback = 0

!     fake assignments to iq and jq to avoid "flint" warning
      iq(qfront) = 0
      jq(qfront) = 0

      return
      end

      subroutine qpush (i, j, iq, jq, qfront, qback)

      implicit none

      integer maxq, ishift, ip, i, j
      parameter (maxq=10000)
      integer qfront, qback, iq(maxq), jq(maxq), maxqsize

      common /qsize/ maxqsize

      qback = qback + 1
      if (qback .gt. maxq) then
        if (qfront .ne. 1) then
!         shift queue left to make room
          ishift = qfront - 1
          do ip=qfront,qback-1
            iq(ip-ishift) = iq(ip)
            jq(ip-ishift) = jq(ip)
          enddo
          qfront = 1
          qback = qback - ishift
        else
          print *, 'queue fault in qpush'
          stop '==>qpush'
        endif
      endif
      iq(qback) = i
      jq(qback) = j

      maxqsize = max(maxqsize, (qback-qfront))

      return
      end

      subroutine qpop (i, j, iq, jq, qfront, qback)

      implicit none

      integer maxq, i, j
      parameter (maxq=10000)
      integer qfront, qback, iq(maxq), jq(maxq)

      i = iq(qfront)
      j = jq(qfront)
      qfront = qfront + 1

      return
      end

      function qempty (iq, jq, qfront, qback)

      implicit none

      integer maxq
      parameter (maxq=10000)
      integer qfront, qback, iq(maxq), jq(maxq)

      logical qempty

      qempty = (qfront .gt. qback)

      return
      end

      function jn(j, jmt)

      implicit none

!     j coordinate to the north of j

      integer j, jmt, jn, offmap
      parameter (offmap = -1)

#if defined O_symmetry
      if (j .lt. jmt-1) then
        jn = j + 1
      elseif (j .eq. jmt-1) then
        jn = jmt-2
      else
        jn = offmap
      endif
#else
      if (j .lt. jmt) then
        jn = j + 1
      else
        jn = offmap
      endif
#endif
      return
      end

      function js(j, jmt)

      implicit none

!     j coordinate to the south of j

      integer j, jmt, js, offmap
      parameter (offmap = -1)

      if (j .gt. 1) then
        js = j - 1
      else
        js = offmap
      endif

      return
      end

      function ie(i,imt)

!     i coordinate to the east of i

      implicit none

      integer i, imt, ie, offmap
      parameter (offmap = -1)

#if defined O_cyclic
      if (i .lt. imt-1) then
        ie = i + 1
      else
        ie = (i+1) - imt + 2
      endif
#else
      if (i .lt. imt) then
        ie = i + 1
      else
        ie = offmap
      endif
#endif
      return
      end

      function iw(i,imt)

      implicit none

!     i coordinate to the west of i

      integer i, imt, iw, offmap
      parameter (offmap = -1)

#if defined O_cyclic
      if (i .gt. 2) then
        iw = i - 1
      else
        iw = (i-1) + imt - 2
      endif
#else
      if (i .gt. 1) then
        iw = i - 1
      else
        iw = offmap
      endif
#endif
      return
      end

      subroutine enter_kmt_changes()

!-----------------------------------------------------------------------
!     copy accumulated changes from the array kmt_changes
!     to kmt
!-----------------------------------------------------------------------

      implicit none

      integer n

      include "size.h"
      include "isleperim.h"

      integer kmt(imt,jmt)

      if (nchanges .eq. 0) then
        print '(t4,a,/)', '<none>'
      else
        auto_kmt_changes = .true.
      endif
      do n=1,nchanges
        kmt(kmt_changes(n,1), kmt_changes(n,2)) = kmt_changes(n,4)
      enddo
      print '(/)'
      nchanges = 0
      return
      end

      subroutine clear_kmt_options ()

!-----------------------------------------------------------------------
!     clear all potential options for kmt changes
!-----------------------------------------------------------------------

      implicit none

      integer i_opt, j_opt, k_opt

      include "isleperim.h"

      do i_opt=1,max_opt
        do j_opt=1,len_opt
          do k_opt=1,4
            kmt_opt(i_opt, j_opt, k_opt) = 0
          enddo
        enddo
      enddo
      return
      end

      subroutine kmt_option (i_opt, i, jrow, kmt_old, kmt_new, kmt)

!-----------------------------------------------------------------------
!     add a change for kmt(i,jrow) from kmt_old to kmt_new to the
!     options list for option i_opt
!-----------------------------------------------------------------------

      implicit none

      integer kmt_new, kmt_old, i_opt, n, k_opt, j_opt, nn, i, jrow

      include "isleperim.h"
      include "size.h"

      integer kmt(imt,jmt)

      if (kmt_new .ne. kmt_old) then

!       if current option is full, accept it and clear it

        if (kmt_opt(i_opt,len_opt,1) .ne. 0) then
          call accept_option(i_opt, kmt)
          do n=1,len_opt
            do k_opt=1,4
              kmt_opt(i_opt,n,k_opt) = 0
            enddo
          enddo
        endif

!      enter the change in kmt_opt

       j_opt = len_opt
       do nn=len_opt,1,-1
          if (kmt_opt(i_opt,nn,1) .eq. 0) then
            j_opt = nn
          endif
        enddo
        kmt_opt(i_opt,j_opt,1) = i
        kmt_opt(i_opt,j_opt,2) = jrow
        kmt_opt(i_opt,j_opt,3) = kmt_old
        kmt_opt(i_opt,j_opt,4) = kmt_new
      endif
      return
      end

      subroutine select_option (problem
     &,                         i, jrow, kmt, xu, yu, zw)

!-----------------------------------------------------------------------
!     use define options to select a default option for changing kmt,
!     then accept and enter the changes.
!-----------------------------------------------------------------------

      implicit none

      character(*) :: problem

      integer idefault, i_opt, j_opt, k_opt, jrow, i

      include "size.h"
      include "isleperim.h"

      integer kmt(imt, jmt)

      real xu(imt), yu(jmt), zw(km)

!     select default option based on ifdef options

      if (problem .eq. 'perim') then
        idefault = 1
        i_opt = idefault
        call accept_option(i_opt, kmt)

      elseif (problem .eq. 'bay') then
        idefault = 1
        if (n_del_kmt .ne. 0) then
          idefault = 3
        endif
        i_opt = idefault
        call accept_option(i_opt, kmt)

       elseif (problem .eq. 'pothole') then
        idefault = 1
        if (n_del_kmt .ne. 0) then
          idefault = 3
        endif
        i_opt = idefault
        call accept_option(i_opt, kmt)

       elseif (problem .eq. 'shallow') then
        idefault = 2
        if (n_del_kmt .ne. 0) then

!         clear option 2

          i_opt = 2
          do j_opt=1,len_opt
            do k_opt=1,4
              kmt_opt(i_opt,j_opt,k_opt) = 0
            enddo
          enddo
          idefault = 2
        endif
        i_opt = idefault
        call accept_option(i_opt, kmt)
      endif

      return
      end

      subroutine accept_option(i_opt,kmt)

!-----------------------------------------------------------------------
!     copy selected option to kmt and kmt_changes
!-----------------------------------------------------------------------

      implicit none

      integer j_opt, i_opt, ii, jj, kk, nn, jch

      include "size.h"
      include "isleperim.h"

      integer kmt(imt, jmt)

      do j_opt=1,len_opt
        if (kmt_opt(i_opt, j_opt, 1) .ne. 0 .and.
     &      kmt_opt(i_opt, j_opt, 2) .ne. 0) then
          ii = kmt_opt(i_opt,j_opt,1)
          jj = kmt_opt(i_opt,j_opt,2)
          kk = kmt_opt(i_opt,j_opt,4)
          kmt(ii,jj) = kk
          nn = nchanges + 1
          if (nn .gt. max_change) then
            print '(/,a,a,/,a,a,i6,a)'
     &,           'WARNING: kmt_changes buffer full.  '
     &,           'Changes written early to delta.kmt file.'
     &,           'To avoid this message, '
     &,          'increase max_change = ', max_change,' in isleperim.h'
            call enter_kmt_changes()
          endif
          nchanges = nchanges + 1
          do jch=1,4
            kmt_changes(nchanges,jch) = kmt_opt(i_opt, j_opt, jch)
          enddo
          kmt(kmt_changes(nchanges,1), kmt_changes(nchanges,2)) =
     &      kmt_changes(nchanges,4)
        endif
      enddo

      return
      end
