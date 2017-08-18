      subroutine setkmp (alat1, slon1, elon1, alat2, slon2, elon2, num)

#if defined O_mom
!-----------------------------------------------------------------------
!     set the topography mask "kmt(i,j)" = "num" within the area of the
!     parallelogram bounded by vertices:
!     (alat1,slon1), (alat1,elon1), (alat2,slon1), & (alat2,elon2)

!-----------------------------------------------------------------------

      implicit none

      integer j1, indp, j2, js, je, i1, i2, is1, ie1, is2, ie2, is, ie
      integer j, i, num

      real alat1, slon1, alat2, elon1, slon2, elon2, rdj

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "levind.h"

!     convert the four vertices into model indices
!     (js,is1), (js,ie1), (je,is2), (je,ie2)

      j1 = indp (alat1, yt, jmt)
      j2 = indp (alat2, yt, jmt)
      js = min (j1,j2)
      je = max (j1,j2)

      i1  = indp (slon1, xt, imt)
      i2  = indp (elon1, xt, imt)
      is1 = min (i1,i2)
      ie1 = max (i1,i2)

      i1  = indp (slon2, xt, imt)
      i2  = indp (elon2, xt, imt)
      is2 = min (i1,i2)
      ie2 = max (i1,i2)

      is = is1
      ie = ie1

!     fill in the area bounded by (js,is1), (js,ie1), (je,is2), (je,ie2)

      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do 100 j=js,je
        do 90 i=is,ie
          kmt(i,j) = num
90      continue
        is = nint(rdj*((j-js)*is2 + (je-j)*is1))
        ie = nint(rdj*((j-js)*ie2 + (je-j)*ie1))
100   continue
#endif

      return
      end
