      subroutine findex (kxx, jjmax, kmax, jf1, jf2, imax, isf, ief)

#if defined O_fourfil || defined O_firfil
!=======================================================================

!  findex finds and prints starting and ending indices
!         for filtering, where:
!             kxx  = field of maximum levels for the quantity
!                    being filtered
!             jjmax = number of rows to be filtered
!             kmax  = maximum number of levels to be filtered
!             jf1   = last row in the south to be filtered
!             jf2   = first row in the north to be filtered
!             imax  = maximum i index to be filtered
!             isf   = returned values of starting indices
!             ief   = returned values of ending indices
!=======================================================================

      implicit none

      integer jjmax, kmax, jj, jrow, jf1, jf2, k, l, i, imax, lm, llast

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "index.h"

      integer kxx(imt,jmt)
      integer isf(jjmax,lsegf,kmax), ief(jjmax,lsegf,kmax)
      integer iis(lsegf+1), iie(lsegf+1)

!=======================================================================
!     begin executable code to find start and end indices
!=======================================================================

      jj = 0
      do jrow = jfrst,jmtm1
        if (jrow.le.jf1 .or. jrow.ge.jf2) then
          jj = jj+1
          do k=1,kmax
            do l=1,lsegf+1
              iis(l) = 0
              iie(l) = 0
            enddo
            l = 1
            if (kxx(2,jrow) .ge. k) then
              iis(1) = 2
            endif
            do i=2,imax-1
              if (kxx(i-1,jrow).lt.k .and. kxx(i,jrow).ge.k) then
                iis(l) = i
              endif
              if (kxx(i,jrow).ge.k .and. kxx(i+1,jrow).lt.k) then
                if (i.ne.iis(l) .or. (i.eq.2 .and. kxx(1,jrow).ge.k))
     &            then
                  iie(l) = i
                  l = l+1
                else
                  iis(l) = 0
                endif
              endif
            enddo
            if (kxx(imax-1,jrow).ge.k .and. kxx(imax,jrow).ge.k) then
              iie(l) = imax-1
              l = l+1
            endif
            lm = l-1
# if defined O_cyclic
            if (lm.gt.1) then
              if (iis(1).eq.2 .and. iie(lm).eq.imax-1
     &                        .and. kxx(1,jrow).ge.k) then
                iis(1) = iis(lm)
                iie(1) = iie(1) + imax-2
                iis(lm) = 0
                iie(lm) = 0
                lm = lm-1
              endif
            endif
# endif
            if (lm .gt. lsegf) then
              write (stdout,99) lsegf, jrow, k
              write (stderr,99) lsegf, jrow, k
              stop ' findex '
            endif
            do l = 1,lsegf
              isf(jj,l,k) = iis(l)
              ief(jj,l,k) = iie(l)
            enddo
          enddo
        endif
      enddo

99    format (/' Error: => Increase lsegf = ',i4
     &,'(the maximum number of ocean strips at any latitude and depth)'
     &,/'           The problem occurs at jrow=',i4 ,'  k=',i3
     &,/'           use "grep lsegf *.h" to find the parameter "lseg"')
#endif

      return
      end
