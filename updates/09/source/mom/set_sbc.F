      subroutine set_sbc(joff, js, je, is, ie, isbc, itr, 
     &                   taut, doAccum)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie, taut
      logical doAccum

      real rts

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"

      if (isbc .le. 0 .or. itr .le. 0) return
      
      if (doAccum) then
!       initialize the Atmos S.B.C. at the start of each ocean segment
!       (do not alter values in land)		

		if (eots .and. osegs) then
		  do j=js,je
			jrow   = j + joff
			do i=is,ie
			  if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
			enddo
		  enddo
		endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

		if (eots) then
		  do j=js,je
			jrow = j + joff
			do i=is,ie
			  sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taut)
			enddo
		  enddo
		endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

		if (eots .and. osege) then
		  rts = c1/ntspos
		  do j=js,je
			jrow   = j + joff
			do i=is,ie
			  if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
			enddo
		  enddo
		endif
      
      else

		do j=js,je
		  jrow = j + joff
		  do i=is,ie
			sbc(i,jrow,isbc) = t(i,1,j,itr,taut)
		  enddo
		enddo
      
      endif

      return
      end
