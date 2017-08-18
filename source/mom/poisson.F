      subroutine border (v, sym)

#if defined O_mom
!-----------------------------------------------------------------------
!     adjust borders of an array for cyclic and symmetry settings

!     symmetry conditions are:
!     't odd':  "t" grid variable: asymmetric reflection at north
!     't even':  "t" grid variable: symmetric reflection at north
!     'u odd':  "u" grid variable: asymmetric reflection at north
!     'u even':  "u" grid variable: symmetric reflection at north
!-----------------------------------------------------------------------

      implicit none

      character(*) :: sym

      integer i, j

      include "size.h"
      include "stdunits.h"

      real v(imt,jmt)

!     set southern border

      do i=2,imt-1
        v(i,1) = 0.0
      enddo

# if defined O_symmetry
      if     (sym .eq. 't even') then
        do i=2,imt-1
          v(i,jmt) = v(i,jmt-1)
        enddo
      elseif (sym .eq. 't odd') then
        do i=2,imt-1
          v(i,jmt) = -v(i,jmt-1)
        enddo
      elseif (sym .eq. 'u even') then
        do i=2,imt-1
          v(i,jmt)   = v(i,jmt-2)
        enddo
      elseif (sym .eq. 'u odd') then
        do i=2,imt-1
          v(i,jmt)   = -v(i,jmt-2)
          v(i,jmt-1) = 0.0
        enddo
      else
        write (stdout,*) '=> Error: illegal sym =',sym
        stop '=>border'
      endif
# else
      do i=2,imt-1
        v(i,jmt) = 0.0
      enddo
# endif

# if defined O_cyclic
      do j=1,jmt
        v(1,j) = v(imt-1,j)
        v(imt,j) = v(2,j)
      enddo
# else
      do j=1,jmt
        v(1,j) = 0.0
        v(imt,j) = 0.0
      enddo
# endif
      return
      end

      subroutine iborder (iv, sym)

!-----------------------------------------------------------------------
!     adjust borders of an array for cyclic and symmetry settings

!     symmetry conditions are:
!     't odd':  "t" grid variable: asymmetric reflection at north
!     't even':  "t" grid variable: symmetric reflection at north
!     'u odd':  "u" grid variable: asymmetric reflection at north
!     'u even':  "u" grid variable: symmetric reflection at north
!-----------------------------------------------------------------------

      implicit none

      character(*) :: sym

      integer i, j

      include "size.h"
      include "stdunits.h"

      real iv(imt,jmt)

!     set southern border

      do i=2,imt-1
        iv(i,1) = 0
      enddo

# if defined O_symmetry
      if     (sym .eq. 't even') then
        do i=2,imt-1
          iv(i,jmt) = iv(i,jmt-1)
        enddo
      elseif (sym .eq. 't odd') then
        write (stdout,*) '=> Error: illegal sym =',sym
        stop '=>iborder'
      elseif (sym .eq. 'u even') then
        do i=2,imt-1
          iv(i,jmt)   = iv(i,jmt-2)
        enddo
      elseif (sym .eq. 'u odd') then
        write (stdout,*) '=> Error: illegal sym =',sym
        stop '=>iborder'
      else
        write (stdout,*) '=> Error: illegal sym =',sym
        stop '=>iborder'
      endif
# else
      do i=2,imt-1
        iv(i,jmt) = 0
      enddo
# endif

# if defined O_cyclic
      do j=1,jmt
        iv(1,j)   = iv(imt-1,j)
        iv(imt,j) = iv(2,j)
      enddo
# else
      do j=1,jmt
        iv(1,j)   = 0
        iv(imt,j) = 0
      enddo
# endif
      return
      end

      subroutine checkerboard (solution, map)

!-----------------------------------------------------------------------
!     removes "checkboard" null space from an array "solution"
!-----------------------------------------------------------------------

      implicit none

      integer noceansum(0:1,0:1), i1, j1, j, i

      real redsum, blacksum, nred, nblack, diff, c
      real sum(0:1,0:1), correction(0:1,0:1)

      include "size.h"

      integer map(imt,jmt)

      real solution(imt,jmt)

      do i1=0,1
        do j1=0,1
          sum(i1,j1) = 0.0
          noceansum(i1,j1) = 0
        enddo
      enddo

      do i1=0,1
        do j1=0,1
          do j=2+j1,jmt-1,2
            do i=2+i1,imt-1,2
              sum(i1,j1) = sum(i1,j1) + solution(i,j)
            enddo
          enddo
        enddo
      enddo

      do i1=0,1
        do j1=0,1
          do j=2+j1,jmt-1,2
            do i=2+i1,imt-1,2
              if (map(i,j) .le. 0) then
                noceansum(i1,j1) = noceansum(i1,j1) + 1
              endif
            enddo
          enddo
        enddo
      enddo

      redsum   = sum(0,0) + sum(1,1)
      blacksum = sum(1,0) + sum(0,1)
      nred     = noceansum(0,0) + noceansum(1,1)
      nblack   = noceansum(1,0) + noceansum(0,1)
      diff = redsum/nred - blacksum/nblack
      c    = diff / 2.0

      print *, ' '
      print '(a,i6,a,i6,a,e14.7)'
     &,         '=> checkerboard: nred = ',nred, ', nblack = ',nblack
     &,         ', removing a checkerboard correction of ', c

      correction (0,0) = -c
      correction (1,1) = -c
      correction (1,0) =  c
      correction (0,1) =  c

      do i1=0,1
        do j1=0,1
          do j=2+j1,jmt-1,2
            do i=2+i1,imt-1,2
              if (map(i,j) .le. 0) then
                solution(i,j) = solution(i,j) + correction(i1,j1)
              endif
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine fill_land (solution, map, noslip,
     &                      nisle, iperm, jperm, iofs, nippts)

      implicit none

      integer nisle

      logical noslip

      include "size.h"

      integer map(imt,jmt), iperm(maxipp), jperm(maxipp)
      integer nippts(mnisle), iofs(mnisle)

      real solution(imt,jmt)

      if (noslip) then
        call fill_land1 (solution, map,
     &                   nisle, iperm, jperm, iofs, nippts)
      else
        call fill_land2 (solution, map,
     &                   nisle, iperm, jperm, iofs, nippts)
      endif

      return
      end

      subroutine fill_land1 (solution, map,
     &                      nisle, iperm, jperm, iofs, nippts)

!=======================================================================

!     fills each land area with the [presumed constant] value
!     of solution along its ocean perimeter.
!=======================================================================

      implicit none

      integer i, j, isle, nisle

      real fill

      include "size.h"

      integer map(imt,jmt), iperm(maxipp), jperm(maxipp)
      integer nippts(mnisle), iofs(mnisle)

      real solution(imt,jmt)

      do i=2,imt-1
        do j=2,jmt-1
          if (map(i,j) .gt. 0) then
            isle = map(i,j)
            fill = solution(iperm(iofs(isle)+1),jperm(iofs(isle)+1))
            solution(i,j) = fill
          endif
        enddo
      enddo

      call mirror_adjust (solution)

      return
      end

      subroutine fill_land2 (solution, map,
     &                      nisle, iperm, jperm, iofs, nippts)

!=======================================================================

!     fills the boundary cells of each land area with the value
!     of solution at the adjacent ocean perimeter point.
!     only NSEW directions are searched [no diagonal directions]
!     in case of multiple ocean perimeter points, their average
!     is used.
!=======================================================================

      implicit none

      integer i, j, isle, nbrs, nisle

      logical last_pass

      real sum

      include "size.h"

      integer map(imt,jmt), iperm(maxipp), jperm(maxipp), nippts(mnisle)
      integer iofs(mnisle)

      real solution(imt,jmt)

1000  continue
      last_pass = .true.

      do i=2,imt-1
        do j=2,jmt-1
          if (map(i,j) .gt. 0) then
            isle = map(i,j)
            nbrs = 0
            sum = 0.0
            if (map(i,j+1) .eq. -isle) then
              nbrs = nbrs + 1
              sum = sum + solution(i,j+1)
            endif
            if (map(i+1,j) .eq. -isle) then
              nbrs = nbrs + 1
              sum = sum + solution(i+1,j)
            endif
            if (map(i,j-1) .eq. -isle) then
              nbrs = nbrs + 1
              sum = sum + solution(i,j-1)
            endif
            if (map(i-1,j) .eq. -isle) then
              nbrs = nbrs + 1
              sum = sum + solution(i-1,j)
            endif
            if (nbrs .gt. 0) then
              solution(i,j) = sum / nbrs
!             last_pass = .false.
            endif
          endif
        enddo
      enddo

      call mirror_adjust (solution)

      if (.not. last_pass) goto 1000

      return
      end

      subroutine mirror_adjust (solution)

!=======================================================================
!     fills each border cell with the value adjacent to it
!=======================================================================

      implicit none

      integer i

      include "size.h"

      real solution(imt,jmt)

      call border(solution, 't odd')
      do i=1,imt
        solution(i,1)   = solution(i,2)
        solution(i,jmt) = solution(i,jmt-1)
      enddo
# if !defined O_cyclic
      do j=2,jmt-1
        solution(1,j)   = solution(2,j)
        solution(imt,j) = solution(imt-1,j)
      enddo
      solution(  1,  1) = solution(    2,    2)
      solution(  1,jmt) = solution(    2,jmt-1)
      solution(imt,  1) = solution(imt-1,    2)
      solution(imt,jmt) = solution(imt-1,jmt-1)
# endif
      return
      end

      subroutine zero_level (surfpres, variable, map, dxt, dyt, cst)

      implicit none

      character(*) :: variable

      integer i, j

      real sum, area_ocean, area, surfpres0

      include "size.h"

      integer map(imt,jmt)

      real surfpres(imt,jmt), dxt(imt), dyt(jmt), cst(jmt)

!     this does not correctly handle multiple basins

      sum = 0.0
      area_ocean = 0.0
      do i=2,imt-1
        do j=2,jmt-1
          if (map(i,j) .le. 0) then
            area = dxt(i)*cst(j)*dyt(j)
            sum = sum + surfpres(i,j)*area
            area_ocean = area_ocean + area
          endif
        enddo
      enddo
      surfpres0 = sum / area_ocean
      call con_adjust (surfpres, surfpres0, map)
      print '(a,e14.7,a,a/)'
     &, '=> zero_level: removing a mean of ', surfpres0, ' from '
     &, variable
      return
      end

      subroutine ddxu (tquant, uquant, dxu, cosu)

!=======================================================================
!     Calculates x partial derivative of field tquant
!     Answer is centered at u/v points
!=======================================================================

      implicit none

      integer i, j

      include "size.h"

      real tquant(imt,jmt), uquant(imt,jmt)
      real dxu(imt), cosu(jmt)

!     calculate partial derivative = ddx (tquant)

      call diffdxu (tquant, uquant)
      do i=1,imt-1
        do j=1,jmt-1
          uquant(i,j) = uquant(i,j) / (dxu(i)*cosu(j))
        enddo
      enddo

      return
      end

      subroutine diffdxu (tquant, uquant)

!=======================================================================
!     Calculates x partial difference of field tquant
!     Answer is centered at u/v points
!=======================================================================

      implicit none

      integer i, j, i1, j1

      real c0, p5, cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real tquant(imt,jmt), uquant(imt,jmt)

!-----------------------------------------------------------------------
!     set locally needed constants
!-----------------------------------------------------------------------

      c0    = 0.0
      p5    = 0.5

!-----------------------------------------------------------------------
!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     calculate partial difference
!        diffdx (tquant) = deltax (bary (tquant))
!-----------------------------------------------------------------------

      do i=1,imt
        do j=1,jmt
          uquant(i,j) = 0.0
        enddo
      enddo

      do i1=0,1
        do j1=0,1
          do i=1,imt-1
            do j=1,jmt-1
              uquant(i,j) = uquant(i,j) + cddxu(i1,j1)*tquant(i+i1,j+j1)
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine ddyu (tquant, uquant, dyu)

!=======================================================================
!     Calculates y partial derivative of field tquant
!     Answer is centered at u/v points
!=======================================================================

      implicit none

      integer i, j

      include "size.h"

      real tquant(imt,jmt), uquant(imt,jmt), dyu(jmt)

      call diffdyu (tquant, uquant)
      do i=1,imt-1
        do j=1,jmt-1
          uquant(i,j) = uquant(i,j) / dyu(j)
        enddo
      enddo

      return
      end

      subroutine diffdyu (tquant, uquant)

!=======================================================================
!     Calculates y partial difference of field tquant
!     Answer is centered at u/v points
!=======================================================================

      implicit none

      integer i, j, i1, j1

      real c0, p5, cddxu(0:1,0:1), cddyu(0:1,0:1)
      real cddxt(-1:0,-1:0), cddyt(-1:0,-1:0)

      include "size.h"

      real tquant(imt,jmt), uquant(imt,jmt)

!-----------------------------------------------------------------------
!     set locally needed constants
!-----------------------------------------------------------------------

      c0    = 0.0
      p5    = 0.5

!-----------------------------------------------------------------------
!     construct coefficients for partial differences. a partial
!     difference in "x" is defined as an "x" difference of a quantity
!     which is averaged in "y". (and symmetrically for "y" differences).
!     Note that this is an x difference and NOT an x derivitive.
!     partial differences of quantities on the "t" grid are defined on
!     the "u" grid and visa versa.
!     therefore partial differences at:
!     u/v points (i,j), involve nearby t/s points with subscripts:
!        (i  ,j+1)    (i+1,j+1)
!        (i  ,j  )    (i+1,j  )
!     t/s points (i,j), involve nearby u/v points with subscripts:
!        (i-1,j  )    (i  ,j  )
!        (i-1,j-1)    (i  ,j-1)
!     thus if qu(i,j) is defined on u/v points, its partial
!     difference ddxqt = ddxt(qu) is defined on t/s points and has the
!     value
!     ddxqt(i,j) = cddxt(-1,-1)*qu(i-1,j-1) + cddxt(-1,0)*qu(i-1,j+0)
!                + cddxt( 0,-1)*qu(i+0,j-1) + cddxt( 0,0)*qu(i+0,j+0)
!-----------------------------------------------------------------------

      cddxu( 0, 0) = -p5
      cddxu( 0, 1) = -p5
      cddxu( 1, 0) =  p5
      cddxu( 1, 1) =  p5

      cddxt(-1,-1) = -p5
      cddxt(-1, 0) = -p5
      cddxt( 0,-1) =  p5
      cddxt( 0, 0) =  p5

      cddyu( 0, 0) = -p5
      cddyu( 0, 1) =  p5
      cddyu( 1, 0) = -p5
      cddyu( 1, 1) =  p5

      cddyt(-1,-1) = -p5
      cddyt(-1, 0) =  p5
      cddyt( 0,-1) = -p5
      cddyt( 0, 0) =  p5

!-----------------------------------------------------------------------
!     calculate partial difference
!        diffdy (tquant) = deltay (barx (tquant))
!-----------------------------------------------------------------------

      do i=1,imt
        do j=1,jmt
          uquant(i,j) = 0.0
        enddo
      enddo

      do i1=0,1
        do j1=0,1
          do i=1,imt-1
            do j=1,jmt-1
              uquant(i,j) = uquant(i,j) + cddyu(i1,j1)*tquant(i+i1,j+j1)
            enddo
          enddo
        enddo
      enddo

      return
      end

      subroutine con_adjust (dpsi, dpsi1, map)

!-----------------------------------------------------------------------
!     the constant dpsi1 is subtracted from dpsi(i,j) at all
!     ocean points (i.e., where map(i,j) .le. 0)
!-----------------------------------------------------------------------

      implicit none

      integer i, j

      real dpsi1

      include "size.h"

      integer map(imt,jmt)

      real dpsi(imt,jmt)

      do i=1,imt
        do j=1,jmt
          if (map(i,j) .le. 0) then
            dpsi(i,j) = dpsi(i,j) - dpsi1
          endif
        enddo
      enddo
      return
      end

      function conv (converged)

!-----------------------------------------------------------------------
!     converts logical to character form for printing
!-----------------------------------------------------------------------

      implicit none

      character(11) :: conv

      logical converged

      if (converged) then
        conv = '[converged]'
      else
        conv = '[diverged] '
      endif

      return
      end

      subroutine compare2(x, y, ax, ay, imax, jmax)

!-----------------------------------------------------------------------
!     compare arrays "x" and "y"
!     ax = alphabetical identifier for "x"
!     ay = alphabetical identifier for "y"
!     prints count of relative differences greater than threshhold
!     which is set in parameter statement below
!-----------------------------------------------------------------------

      implicit none

      character(*) :: ax, ay

      integer imax, jmax, numneq, numbadrel, i, j

      real threshhold, relerrormax, badabserrormax, abserrormax
      real relerror, relerr, x(imax,jmax), y(imax,jmax)
      parameter (threshhold= 1.0e-6)

      write (6,*) 'comparing ',ax,' with ',ay

      numneq = 0
      numbadrel = 0
      relerrormax = 0.0
      badabserrormax = 0.0
      abserrormax = 0.0

      do i=1,imax
        do j=1,jmax
          if (x(i,j) .ne. y(i,j) ) then
            numneq = numneq + 1
            relerror = relerr (x(i,j), y(i,j))
            relerrormax = max (relerror, relerrormax)
            if (numneq .le. 20 .and.
     &          relerror .gt. threshhold) then
              write(6,'(3a10,2i4,tr1,2e25.18,a,e25.18)')
     &            '-->',ax,ay,i,j,x(i,j),y(i,j),
     &            ' rel error = ', relerror
            endif
            if (relerror .gt. threshhold) then
              numbadrel = numbadrel + 1
              badabserrormax = max(badabserrormax, abs(x(i,j)-y(i,j)))
            endif
            abserrormax = max(abserrormax, abs(x(i,j)-y(i,j)))
          endif
        enddo
      enddo

      if (numneq .ne. 0) then
        write(6,*) numneq, ' entries differ'
        write(6,*) numbadrel, ' have relative error >', threshhold
        write(6,*) '  bad max absolute error = ',badabserrormax
        write(6,*) '  maximum relative error = ',relerrormax
        write(6,*) '  maximum absolute error = ',abserrormax
      else
        write(6,*) 'no differences detected'
      endif

      return
      end

      function relerr (x, y)

      implicit none

      real x, y, relerr, xymin, xymax

      if (x .eq. y) then
        relerr = 0.0
      else
        xymin = min(abs(x), abs(y))
        xymax = max(abs(x), abs(y))
        if (xymin .gt. 1.0e-23 * xymax) then
          relerr = abs(x-y)/xymin
        else
          relerr = xymax
        endif
      endif
#endif

      return
      end
