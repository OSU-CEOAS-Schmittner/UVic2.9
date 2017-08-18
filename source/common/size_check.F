      subroutine size_check(imt2, jmt2, km2, sub, option)

!-----------------------------------------------------------------------
!     check that array bounds (imt2, jmt2, km2) = (imt, jmt, km)

!     inputs:
!       imt2 = input value for imt
!       jmt2 = input value for jmt
!       km2  = input value for km
!       sub  = name of the subroutine requesting the size check
!       option = what to do if size check fails:
!              'stop' is the only option
!-----------------------------------------------------------------------

      implicit none

      character(*) :: sub, option
      character(60) :: msg

      integer imt2, jmt2, km2

      include "size.h"

      if (imt .ne. imt2 .or. jmt .ne. jmt2 .or. km .ne. km2) then
        print '(a/2(a,i4,a,i4,a,i4,a,a,/))', '==>Error:  size_check '
     &,       'imt = ', imt2, '  jmt = ',jmt2,'  km = ',km2
     &,       ' in ', sub
     &,       'imt = ', imt, '  jmt = ',jmt,'  km = ',km
     &,       ' in "size.h"'
        print '(/,a,a,a)'
     &,       'Sizes in ',sub,' are incompatible with "size.h"'
        stop
      endif
      return
      end
