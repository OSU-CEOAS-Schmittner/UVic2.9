      subroutine blmixi

#if defined O_mom
# if defined O_bryan_lewis_vertical || defined O_bryan_lewis_horizontal

!-----------------------------------------------------------------------
!     Reference:
!     A Water Mass Model of the World Ocean  K. Bryan, L.J. Lewis
!     JGR, vol 84, No. C5, May 20, 1979
!-----------------------------------------------------------------------

      implicit none

      integer k, ioun

      real afkph, dfkph, sfkph, zfkph, pi, ahs, ahb

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "hmixc.h"
      include "vmixc.h"

      data afkph, dfkph, sfkph, zfkph /0.8, 1.05, 4.5e-5, 2500.0e2/

      namelist /blmix/  Ahv, Ahh
#  if defined O_bryan_lewis_vertical

!------------------------------------------------------------------------
!     Use Bryan & Lewis values for vertical tracer diffusion
!     Ahv range of 0.3 to 1.3, crossover at 2500m.
!------------------------------------------------------------------------

!     compute depth dependent vertical diffusion coefficients for
!     tracers using the relationship of Bryan and Lewis

      pi = 4.0 * atan(1.0)

      do k=1,km
        Ahv(k) = (afkph + (dfkph/pi)*(atan(sfkph*(zw(k) - zfkph))))
      enddo

#  endif
#  if defined O_bryan_lewis_horizontal

!------------------------------------------------------------------------
!     compute depth dependent horizontal diffusion coefficients for
!     tracers using the relationship of Bryan and Lewis
!------------------------------------------------------------------------

      ahs = 5.0e+3
      ahb = 1.0e+3
      do k=1,km
        Ahh(k) = (ahb + (ahs - ahb)*exp(-zt(k)/50000.0))*1.0e4
      enddo
#  else
      Ahh = 0.
#  endif

      write (stdout,'(/a/)')
     &  'B R Y A N - L E W I S   M I X I N G   C O E F F S'

      call getunit (ioun, 'control.in'
     &,               'formatted sequential rewind')
      read  (ioun, blmix, end=100)
100   continue
      write (stdout,blmix)
      call relunit (ioun)

      write (stdout,'(/)')
# endif
#endif

      return
      end
