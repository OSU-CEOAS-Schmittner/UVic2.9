      function indp (value, array, ia)

!=======================================================================

!     indp = index of nearest data point within "array" corresponding to
!            "value".

!     inputs:

!     value  = arbitrary data...same units as elements in "array"
!     array  = array of data points  (must be monotonically increasing)
!     ia     = dimension of "array"

!     output:

!     indp =  index of nearest data point to "value"
!             if "value" is outside the domain of "array" then indp = 1
!             or "ia" depending on whether array(1) or array(ia) is
!             closest to "value"

!             note: if "array" is dimensioned array(0:ia) in the calling
!                   program, then the returned index should be reduced
!                   by one to account for the zero base.

!     example:

!     let model depths be defined by the following:
!     parameter (km=5)
!     dimension z(km)
!     data z /5.0, 10.0, 50.0, 100.0, 250.0/

!     k1 = indp (12.5, z, km)
!     k2 = indp (0.0, z, km)

!     k1 would be set to 2, & k2 would be set to 1 so that
!     z(k1) would be the nearest data point to 12.5 and z(k2) would
!     be the nearest data point to 0.0

!=======================================================================

      implicit none

      integer ia, i, ii, indp

      real value

      include "stdunits.h"

      real array(ia)

      do i=2,ia
        if (array(i) .lt. array(i-1)) then
         write (stdout,*)
     &   ' => Error: array must be monotonically increasing in "indp"'
     &,  '           when searching for nearest element to value=',value
          write (stdout,*) '           array(i) < array(i-1) for i=',i
          write (stdout,*) '           array(i) for i=1..ia follows:'
          do ii=1,ia
            write (stdout,*) 'i=',ii, ' array(i)=',array(ii)
          enddo
          stop '=>indp'
        endif
      enddo
      if (value .lt. array(1) .or. value .gt. array(ia)) then
        if (value .lt. array(1))  indp = 1
        if (value .gt. array(ia)) indp = ia
        return
      else
        do i=2,ia
          if (value .le. array(i)) then
            indp = i
            if (array(i)-value .gt. value-array(i-1)) indp = i-1
            go to 101
          endif
        enddo
101     continue
      endif
      return
      end

      subroutine ftc (f, if, jf, xf, yf, c, ic, jc, istart, iend
     &,               jstart, jend, xc, yc, init, work, lenw)

!=======================================================================

!     "ftc" is a mnemonic for "fine to coarse".

!     obtain a coarse grid representation of a fine grid dataset by area
!     averaging grid box values on the fine grid which overlay coarse
!     grid boxes. note: the coarse grid boxes do not have to contain an
!     integral number of fine grid boxes.

!     inputs:

!     f      = data on fine grid
!     if     = inner dimension of "f"
!     jf     = outer dimension of "f"
!     xf     = coordinates for inner dimension of "f" (eg: longitudes)
!     yf     = coordinates for outer dimension of "f" (eg: latitudes)

!     ic     = inner dimension of coarse grid "c"
!     jc     = outer dimension of coarse grid "c"
!     istart = starting index along inner dimension of "c" for which
!              averaged values are desired
!     iend   = ending index along inner dimension of "c" for which
!              averaged values are desired
!     jstart = starting index along outer dimension of "c" for which
!              averaged values are desired
!     jend   = ending index along outer dimension of "c" for which
!              averaged values are desired
!     xc     = coordinates for inner dimension of "c" (eg: longitudes)
!     yc     = coordinates for outer dimension of "c" (eg: latitudes)
!     init   = initialize the averaging factors
!              "init" should be set = 1 on the first call.
!              "init" <> 1 uses the previously computed factors stored
!              in "work" array.
!     work   = work array of averaging factors when "init" <> 1
!              (previously calculated by "ftc" when "init" = 1)
!     lenw   = size of work array. lenw should be >= 9*max(if,jf)

!     output:

!     c  = coarse grid average of "f" defined over
!          ((c(i,j),i=istart,iend),j=jstart,jend)
!     work   = work array of averaging factors when "init" = 1

!     restrictions:

!     fine and coarse grids are assumed rectangular with "xf" and "xc"
!     having the same units. "yf" and "yc" must also have the same units
!     the coarse domain xc(istart)...xc(iend) must be within
!     the fine domain xf(1)...xf(if). similarly,
!     yc(jstart)...yc(jend) must be within yf(1)...yf(jf). all
!     coordinates must be strictly monotonically increasing.
!=======================================================================

      include "stdunits.h"
      logical error, show_coord
      parameter (len=10000, p5=0.5, c0=0.0)
      dimension iso(0:len), ieo(0:len), jso(0:len), jeo(0:len)
     &,         dx(0:len,2), dy(0:len,2), edgecx(0:len), edgecy(0:len)
     &,         edgefx(0:len), edgefy(0:len)
      dimension f(if,jf), xf(if), yf(jf)
      dimension c(ic,jc), xc(ic), yc(jc)
      dimension work(lenw)

        write (stdout,*) ' '
        write (stdout,*)
     & '            Averaging data from "fine" to "coarse" grid'

!-----------------------------------------------------------------------
!     initialize weights or use previously calculated weights
!-----------------------------------------------------------------------

      if (init .eq. 1) then
        error = .false.
        write (stdout,*)
     & '              (initializing the averaging weights)'
        write (stdout,*) ' '

!       test to verify that array sizes do not exceed limits

        if (if .gt. len .or. jf .gt. len) then
          i = max(if,jf)
          write (stdout,*) '=>Error: increase "len" in "ftc" to ',i
          stop '=>ftc'
        endif
        if (lenw .lt. 9*max(if,jf)) then
          write (stdout,*) '=>Error: increase size of "work" array',
     &      ' to at least ',9*max(if,jf),' for calls to "ftc"'
          error = .true.
        endif

!       verify that the "coarse" grid lies within the "fine" grid

        if (xf(1) .gt. xc(istart) .or. xf(if) .lt. xc(iend)) then
          write (stdout,*)
     &     '=>Warning: Coarse grid "xc" is outside "fine" grid "xf".'
          if (xf(1) .gt. xc(istart)) then
            write (stdout,*) '  xc(',istart,')  .lt.  xf(1)'
          endif
          if (xf(if) .lt. xc(iend)) then
            write (stdout,*) '  xc(',iend,')  .gt.  xf(',if,')'
          endif
        endif
        if (yf(1) .gt. yc(jstart) .or. yf(jf) .lt. yc(jend)) then
          write (stdout,*)
     &     '=>Warning: Coarse grid "yc" is outside "fine" grid "yf".'
          if (yf(1) .gt. yc(jstart)) then
            write (stdout,*) '  yc(',jstart,') .lt. yf(1)'
          endif
          if (yf(jf) .lt. yc(jend)) then
            write (stdout,*) '  yc(',jend,')  .gt.  yf(',jf,')'
          endif
        endif

!       construct edges of "coarse" grid boxes

        do i=1,ic-1
          edgecx(i) = p5*(xc(i) + xc(i+1))
        enddo
        edgecx(0)  = xc(1) - (edgecx(1) - xc(1))
        edgecx(ic) = xc(ic) + (xc(ic) - edgecx(ic-1))

        do j=1,jc-1
          edgecy(j) = p5*(yc(j) + yc(j+1))
        enddo
        edgecy(0)  = yc(1) - (edgecy(1) - yc(1))
        edgecy(jc) = yc(jc) + (yc(jc) - edgecy(jc-1))

!       construct edges of "fine" grid boxes

        do i=1,if-1
          edgefx(i) = p5*(xf(i) + xf(i+1))
        enddo
        edgefx(0)  = xf(1) - (edgefx(1) - xf(1))
        edgefx(if) = xf(if) + (xf(if) - edgefx(if-1))

        do j=1,jf-1
          edgefy(j) = p5*(yf(j) + yf(j+1))
        enddo
        edgefy(0)  = yf(1) - (edgefy(1) - yf(1))
        edgefy(jf) = yf(jf) + (yf(jf) - edgefy(jf-1))

!       calculate "dx" and "dy" for the "fine" grid boxes

        do i=1,if
          dx(i,1) = edgefx(i) - edgefx(i-1)
          dx(i,2) = dx(i,1)
        enddo
        dx(0,1) = dx(1,1)
        dx(0,2) = dx(1,2)

        do j=1,jf
          dy(j,1) = edgefy(j) - edgefy(j-1)
          dy(j,2) = dy(j,1)
        enddo
        dy(0,1) = dy(1,1)
        dy(0,2) = dy(1,2)

!       modify "dx" and "dy" for possibly partial "fine" grid boxes
!       near the edges of each coarse grid box.
!       "ii" is the index of the fine grid box which contains the
!       eastern edge of coarse grid box with index "i".
!       dx(ii,1) is the portion of the fine grid box to the west of the
!       edge and dx(ii,2) is the portion to the east. similarly,
!       dy(jj,1) is to the south and dy(jj,2) is to the north of the
!       northern edge of coarse box with index "j".

!       note: edgefx and edgefy are zero based and need the -1 when
!       using "indp"

        do i=0,ic
          ii   = indp (edgecx(i), edgefx, if+1) - 1
          frac = abs(edgecx(i) - edgefx(ii))
          if (edgefx(ii) .lt. edgecx(i)) then
            ii = ii + 1
            dx(ii,2) = (edgefx(min(if,ii)) - edgefx(ii-1)) - frac
            dx(ii,1) = frac
          else
            dx(ii,2) = frac
            dx(ii,1) = (edgefx(ii) - edgefx(max(ii-1,0))) - frac
          endif
          ieo(i) = min(if,max(1,ii))
        enddo
        do i=1,ic
          iso(i) = max(1,ieo(i-1))
        enddo
        iso(0) = ieo(0)

        do j=0,jc
          jj   = indp (edgecy(j), edgefy, jf+1) - 1
          frac = abs(edgecy(j) - edgefy(jj))
          if (edgefy(jj) .lt. edgecy(j)) then
            jj = jj + 1
            dy(jj,2) = (edgefy(min(jf,jj)) - edgefy(jj-1)) - frac
            dy(jj,1) = frac
          else
            dy(jj,2) = frac
            dy(jj,1) = (edgefy(jj) - edgefy(max(jj-1,0))) - frac
          endif
          jeo(j) = min(jf,max(1,jj))
        enddo
        do j=1,jc
          jso(j) = max(1,jeo(j-1))
        enddo
        jso(0) = jeo(0)

!       store the weights into the "work" array

        indx = 1
        do j=0,jc
          work(indx)   = jso(j)
          work(indx+1) = jeo(j)
          indx         = indx + 2
        enddo

        do i=0,ic
          work(indx)   = iso(i)
          work(indx+1) = ieo(i)
          indx         = indx + 2
        enddo

        do j=0,jf
          work(indx)   = dy(j,1)
          work(indx+1) = dy(j,2)
          indx         = indx + 2
        enddo

        do i=0,if
          work(indx)   = dx(i,1)
          work(indx+1) = dx(i,2)
          indx         = indx + 2
        enddo

!       verify that coarse grid is coarser than the fine grid

        do j=jstart,jend
          if ((jso(j) .eq. jso(j+1)) .and. (yc(j) .ge. yf(1))
     &      .and. (yc(j) .le. yf(jf))) then
            write (stdout,*)
     &         '=>Warning: "Coarse" grid is finer than "fine" grid'
     &,     ' near yf(',jso(j),') =',yf(jso(j))
     &,     ' (average may not be accurate)'
          endif
        enddo

        do i=istart,iend
          if ((iso(i) .eq. iso(i+1)) .and. (xc(i) .ge. xf(1))
     &      .and. (xc(i) .le. xf(if))) then
            write (stdout,*)
     &            '=>Warning: "Coarse" grid is finer than "fine" grid'
     &,     ' near xf(',iso(i),') = ',xf(iso(i))
     &,     ' (average may not be accurate)'
          endif
        enddo
        show_coord = .false.
        if (error .or. show_coord) then
          write (stdout,*)
     & ' Indices for averaging fine grid to coarse grid:'
          write (stdout,*)
     & ' (fractional grid boxes are accounted for)'
          write (stdout,8700)
          write (stdout,9000) (m,iso(m),ieo(m),m=istart,iend)
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "xc" ='
          write (stdout,8500) xc
          write (stdout,*) ' Coordinates for fine grid points "xf" ='
          write (stdout,8500) xf

          write (stdout,8800)
          write (stdout,9000) (m,jso(m),jeo(m),m=jstart,jend)
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "yc" ='
          write (stdout,8500) yc
          write (stdout,*) ' Coordinates for fine grid points "yf" ='
          write (stdout,8500) yf
        endif
        if (error) stop '=>ftc'
      else
        write (stdout,*)
     & '              (using previously initialized averaging weights)'
        write (stdout,*) ' '

!       extract the weights from the "work" array

        indx = 1
        do j=0,jc
          jso(j) = nint(work(indx))
          jeo(j) = nint(work(indx+1))
          indx   = indx + 2
        enddo

        do i=0,ic
          iso(i) = nint(work(indx))
          ieo(i) = nint(work(indx+1))
          indx   = indx + 2
        enddo

        do j=0,jf
          dy(j,1) = work(indx)
          dy(j,2) = work(indx+1)
          indx    = indx + 2
        enddo

        do i=0,if
          dx(i,1) = work(indx)
          dx(i,2) = work(indx+1)
          indx    = indx + 2
        enddo
      endif

!-----------------------------------------------------------------------
!     average the "fine" grid to the "coarse" grid
!-----------------------------------------------------------------------

      do m=jstart,jend
        do i=istart,iend
          weight = c0
          sum    = c0
          do j=jso(m),jeo(m)
            indy = 2
            if (j .eq. jeo(m)) indy = 1
            wty = dy(j,indy)
            do ii=iso(i),ieo(i)
              indx = 2
              if (ii .eq. ieo(i)) indx = 1
              area   = dx(ii,indx)*wty
              weight = weight + area
              sum    = sum + f(ii,j)*area
            enddo
          enddo
          c(i,m) = sum/weight
        enddo
      enddo
      return
8500  format (1x,10g11.4)
8700  format (/' Along the 1st dimension, the form is (Coarse grid'
     &,' point "xc": range of fine grid points "xf" to average)'/)
8800  format (/' Along the 2nd dimension, the form is (Coarse grid'
     &,' point "yc": range of fine grid points "yf" to average)'/)
9000  format (5(1x,'(',i4,': ',i4,' to ',i4,')'),/)
      end

      subroutine ctf (c, ic, jc, xc, yc, f, if, jf, istart, iend
     &,               jstart, jend, xf, yf, init, work, lenw)

!=======================================================================

!     "ctf" is a mnemonic for "coarse to fine".

!     obtain a fine grid representation of a coarse grid dataset by
!     linear interpolation of grid box values on the coarse grid to grid
!     boxes on the fine grid.

!     inputs:

!     c      = coarse grid data
!     ic     = inner dimension of coarse grid "c"
!     jc     = outer dimension of coarse grid "c"
!     xc     = coordinates for inner dimension of "c" (eg: longitudes)
!     yc     = coordinates for outer dimension of "c" (eg: latitudes)

!     if     = inner dimension of "f"
!     jf     = outer dimension of "f"
!     xf     = coordinates for inner dimension of "f" (eg: longitudes)
!     yf     = coordinates for outer dimension of "f" (eg: latitudes)

!     istart = starting index along inner dimension of "f" for which
!              interpolated values are desired
!     iend   = ending index along inner dimension of "f" for which
!              interpolated values are desired
!     jstart = starting index along outer dimension of "f" for which
!              interpolated values are desired
!     jend   = ending index along outer dimension of "f" for which
!              interpolated values are desired
!     init   = initialize the interpolation factors
!              "init" should be set = 1 on the first call.
!              "init" <> 1 uses the previously computed factors stored
!              in "work" array.
!     work   = work array of interpolation factors when "init" <> 1
!              (previously calculated by "ctf" when "init" = 1)
!     lenw   = size of work array. lenw should be >= 8*max(if,jf)

!     output:

!     f      = interpolated data on fine grid defined over
!              ((f(i,j),i=istart,iend),j=jstart,jend)
!     work   = work array of interpolation factors when "init" = 1

!     restrictions:

!     fine and coarse grids are assumed rectangular with "xf" and "xc"
!     having the same units. "yf" and "yc" must also have the same units
!     the fine domain xf(istart)...xf(iend) must be within
!     the coarse domain xc(1)...xc(ic). Similarly,
!     and yc(js)...yc(je) must be within yf(1)...yf(jf). all coordinates
!     must be strictly monotonically increasing.
!=======================================================================

      include "stdunits.h"
      logical error, show_coord
      parameter (len=10000, p5=0.5, c0=0.0)
      dimension indxi(len), indxj(len), dnorth(len), dsouth(len)
     &,         deast(len), dwest(len), width(len), height(len)
      dimension f(if,jf), xf(if), yf(jf)
      dimension c(ic,jc), xc(ic), yc(jc)
      dimension work(lenw)

        write (stdout,*) ' '
        write (stdout,*)
     & '            Interpolating data from "coarse" to "fine" grid'

!-----------------------------------------------------------------------
!     initialize weights or use previously calculated weights
!-----------------------------------------------------------------------

      if (init .eq. 1) then
        error = .false.
        write (stdout,*)
     & '              (initializing interpolation weights)'
        write (stdout,*) ' '

!       test to verify that array sizes do not exceed limits

        if (if .gt. len .or. jf .gt. len) then
          i = max(if,jf)
          write (stdout,*) '=>Error: increase "len" in "ctf" to ',i
          error = .true.
        endif
        if (lenw .lt. 8*max(if,jf)) then
          write (stdout,*) '=>Error: increase size of "work" array',
     &      ' to at least ',8*max(if,jf),' for calls to "ctf"'
          error = .true.
        endif

!       verify that the "fine" grid lies within the "coarse" grid

        epsilon = 1.e-5
        xcminus = xc(1) - epsilon*(xc(2)-xc(1))
        xcplus = xc(ic) + epsilon*(xc(ic)-xc(ic-1))
        if (xf(istart) .lt. xcminus .or. xf(iend) .gt. xcplus) then
          error = .true.
          write (stdout,*)
     &          '=>Warning: "fine" grid outside "coarse" grid in "ctf".'
          if (xf(istart) .lt. xc(1)) then
            write (stdout,*) '  xf(',istart,')  .lt.  xc(1)'
          endif
          if (xc(ic) .lt. xf(iend)) then
            write (stdout,*) '  xf(',iend,')  .gt.  xc(',ic,')'
          endif
        endif
        ycminus = yc(1) - epsilon*(yc(2)-yc(1))
        ycplus = yc(jc) + epsilon*(yc(jc)-yc(jc-1))
        if (ycminus .gt. yf(jstart) .or. ycplus .lt. yf(jend)) then
          error = .true.
          write (stdout,*)
     &          '=>Warning: "fine" grid outside "coarse" grid in "ctf".'
          if (yc(1) .gt. yf(jstart)) then
            write (stdout,*) '  yf(',jstart,') .lt. yc(1)'
          endif
          if (yc(jc) .lt. yf(jend)) then
            write (stdout,*) '  yf(',jend,')  .gt.  yc(',jc,')'
          endif
        endif

!       find interpolation factors

        indx = 1
        do j=jstart,jend
          jj = indp (yf(j), yc, jc)
          if (yc(jj) .gt. yf(j) .or. jj .eq. jc) jj = jj - 1
          indxj(j) = jj
          dnorth(j) = yc(jj+1) - yf(j)
          dsouth(j) = yf(j) - yc(jj)
          height(j) = yc(jj+1) - yc(jj)

!         store into "work" array for future use (when "init" <> 1)

          work(indx)   = indxj(j)
          work(indx+1) = dnorth(j)
          work(indx+2) = dsouth(j)
          work(indx+3) = height(j)
          indx         = indx + 4
        enddo

        do i=istart,iend
          ii = indp (xf(i), xc, ic)
          if (xc(ii) .gt. xf(i) .or. ii .eq. ic) ii = ii - 1
          indxi(i) = ii
          deast(i) = xc(ii+1) - xf(i)
          dwest(i) = xf(i) - xc(ii)
          width(i) = xc(ii+1) - xc(ii)

!         store into "work" array for future use (when "init" <> 1)

          work(indx)   = indxi(i)
          work(indx+1) = deast(i)
          work(indx+2) = dwest(i)
          work(indx+3) = width(i)
          indx         = indx + 4
        enddo
        show_coord = .false.
        if (error .or. show_coord) then
          write (stdout,*) ' '
          write (stdout,*) ' Coordinates for coarse grid points "xc" ='
          write (stdout,8500) xc
          write (stdout,*) ' Coordinates for fine grid points "xf" ='
          write (stdout,8500) xf

          write (stdout,*) ' Coordinates for coarse grid points "yc" ='
          write (stdout,8500) yc
          write (stdout,*) ' Coordinates for fine grid points "yf" ='
          write (stdout,8500) yf
        endif
        if (error) stop '=>ctf'
      else

!       extract previously calculated interpolation weights from "work"

        write (stdout,*)
     &'      (using previously initialized interpolation weights)'

        indx = 1
        do j=jstart,jend
          indxj(j)  = nint(work(indx))
          dnorth(j) = work(indx+1)
          dsouth(j) = work(indx+2)
          height(j) = work(indx+3)
          indx      = indx + 4
        enddo

        do i=istart,iend
          indxi(i) = nint(work(indx))
          deast(i) = work(indx+1)
          dwest(i) = work(indx+2)
          width(i) = work(indx+3)
          indx     = indx + 4
        enddo
      endif

!-----------------------------------------------------------------------
!     interpolate data from "coarse" to "fine" grid
!-----------------------------------------------------------------------

      do jj=jstart,jend
        j = indxj(jj)
        do ii=istart,iend
          i = indxi(ii)
          f(ii,jj) = (c(i,j)    *deast(ii)*dnorth(jj)
     &              + c(i+1,j)  *dwest(ii)*dnorth(jj)
     &              + c(i,j+1)  *deast(ii)*dsouth(jj)
     &              + c(i+1,j+1)*dwest(ii)*dsouth(jj)) /
     &               (width(ii)*height(jj))
        enddo
      enddo

      return
8500  format (1x,10g11.3)
      end

      subroutine extrap (a, land, sor, res, il, jl, maxscn, crit, text
     &,                  gtype)

!=======================================================================

!     utility to extrapolate values into land areas neglecting
!     non-uniformity or asymmetry in the grid by solving a simple
!     heat eqn: del**2(a) = 0 over land areas using values over ocean
!     areas as boundary conditions.
!     this alleviates the problem of mismatched land/sea areas due to
!     different geometries or resolutions when interpolating between
!     atmospheric and ocean model grids.
!     the intent is to force reasonable values into land areas
!     near coastlines. far from coasts, the extrapolations may not be
!     reasonable.

!     note: the values over land are used as an initial guess field
!           and need to be specified

!     inputs:

!     a       = array with land areas to be filled. land areas contain
!               initial guess field.
!     land    = mask = (0, non zero) to indicate (land, non land) area
!     il      = number of points along 1st dimension to be filled
!     jl      = number of points along 2nd dimension to be filled
!     maxscn  = maximum number of passes allowed in relaxation
!     crit    = criterion for ending relaxation before "maxscn" limit
!     text    = character string (up to 15 chars) to identify data
!     gtype   = grid type = (1,2) to identify (ocean, atmosphere) grid
!     sor     = scratch area
!     res     = scratch area

!     outputs:

!     a       = array with extrapolated values in land areas.
!               non land areas remain unchanged.
!=======================================================================

      logical done
      include "stdunits.h"
      integer gtype
      character(*) :: text
      parameter (c0=0.0, p25=0.25)
      dimension a(il,jl), land(il,jl), res(il,jl), sor(il,jl)

!-----------------------------------------------------------------------

!     solve a simple Poisson eqn by relaxation to extrapolate data into
!     land areas using values over non land areas as boundary values.

!     note: successive calls to extrap will require fewer scans because
!           the initial guess field over land areas gets better with
!           each call.
!-----------------------------------------------------------------------

!     check on the grid type: atmosphere or ocean

      if (gtype .ne. 1 .and. gtype .ne. 2) then
        write (stdout,98) gtype
        stop '=>extrap'
      endif

!-----------------------------------------------------------------------
!     set the relaxation coefficient to zero over ocean or air
!     relc is somewhat arbitrary
!-----------------------------------------------------------------------

      relc = 0.6
      do j=1,jl
        do i=1,il
          if (land(i,j) .eq. 0) then
            sor(i,j) = relc
          else
            sor(i,j) = c0
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     iterate until errors are acceptable.
!-----------------------------------------------------------------------

      n = 0
100   continue
        resmax = c0
        done   = .true.
        n    = n + 1
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = p25*(a(i-1,j) + a(i+1,j) + a(i,j-1) + a(i,j+1))
     &                 - a(i,j)
          enddo
        enddo
        do j=2,jl-1
          do i=2,il-1
            res(i,j) = res(i,j)*sor(i,j)
            a(i,j) = a(i,j) + res(i,j)
            absres = abs(res(i,j))
            if (absres .gt. crit) done = .false.
            resmax = max(absres,resmax)
          enddo
        enddo

!-----------------------------------------------------------------------
!       set conditions at edge of grid
!-----------------------------------------------------------------------

        if (gtype .eq. 1) then

!         use cyclic or no flux conditions on ocean grids

          do j=1,jl
#if defined O_cyclic
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
#else
            a(1,j)  = a(2,j)
            a(il,j) = a(il-1,j)
#endif
          enddo
        elseif (gtype .eq. 2) then

!         always put cyclic conditions on atmosphere grids

          do j=1,jl
            a(1,j)  = a(il-1,j)
            a(il,j) = a(2,j)
          enddo
        endif

!       no flux condition at northern and southern boundaries

        do i=1,il
          a(i,1)  = a(i,2)
          a(i,jl) = a(i,jl-1)
          enddo

      if (.not. done .and. n .le. maxscn) go to 100

      write (stdout,99) text, n, resmax
99    format (1x,'==> Extrapolated ',a15,' into land using ',i4
     &,       ' scans.  max residual=', g14.7)
98    format (1x,'==> Error:   gtype =',i6,' in extrap')
      return
      end

      subroutine setbcx (a, imt, jmtorkm)

!=======================================================================
!     set zonal boundary condition on the first index of array "a" for
!     every second index. the first index corresponds to the "x"
!     or longitude direction.

!     input:
!      a = array in need of setting the zonal b.c.
!     output
!      a = array with zonal b.c. set
!=======================================================================

      dimension a(imt,jmtorkm)
      do k=1,jmtorkm
#if defined O_cyclic
        a(1,k)   = a(imt-1,k)
        a(imt,k) = a(2,k)
#else
        a(1,k)   = 0.0
        a(imt,k) = 0.0
#endif
      enddo

      return
      end

      subroutine iplot (iarray, im, il, jl)

!=======================================================================
!      map integer array "iarray" into characters for printing with
!      format (a1) to provide a contour map of the integer field.
!      note: max number of unique characters = 120

!     inputs:

!     iarray = integer array to be plotted
!     im     = inner dimension of "iarray"
!     il     = number of points along inner dimension to plot (along x)
!     jl     = number of points along outer dimension to plot (along y)

!     output: prints contour map of "iarray"
!=======================================================================

      include "stdunits.h"
      dimension iarray(im,jl)
      character(120) :: levels, lev1
      save levels
      write (stdout,*) ' '

!     set character markers

      lev1(1:51) = '.abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWX'
      levels = lev1(1:51)//'YZ0123456789+*-=!@#$%<>[]{}()'

!     find range of integers

      maxint = iarray(1,1)
      minint = iarray(1,1)
      do j=1,jl
        do i=1,il
          maxint = max(maxint,iarray(i,j))
          minint = min(minint,iarray(i,j))
        enddo
      enddo

!     show mapping of integers into characters

      write (stdout,*) ' '
      write (stdout,*)
     & ' "iplot" mapping of integers to characters is as follows:'
      inc  = 3
      last = min(minint+80-1,maxint)
      do i=minint,last,inc
        ii = i-minint+1
        if (i+inc .le. last) then
          jinc = inc
        else
          jinc = last-i+1
        endif
        write (stdout,'(6(1x,i6,a,a,3x))')
     &  (j+minint-1," is printed as ",levels(j:j),j=ii,ii+jinc-1)
      enddo
      write (stdout,*) ' '

      if (maxint - minint + 1 .gt. 80) then
        write (stdout,*)
     & ' => Note: there are ',maxint-minint+1,' integers in the field'
        write (stdout,*) '          "iplot" cannot uniquely assign '
     &,' more than 120 characters for plotting symbols.'
        write (stdout,*) '          therefore integers are '
     &, 'represented by cyclically reusing the list of plotting symbols'
        write (stdout,*) '   '
      endif

!     print character representation of integers

      inc=124
      do l=0,il,inc
        incr = min(inc,il-l)
        write (stdout,8800) (l+i,i=1,incr,4)
        do  jj=1,jl
          j  = jl+1-jj
!          write (stdout,8900) j, (levels(min(80,iarray(l+i,j)-minint+1):
!     &                        min(80,iarray(l+i,j)-minint+1)),i=1,incr)
          write (stdout,8900) j,
     &  (levels(mod(iarray(l+i,j)-minint+1-1,80)+1:
     &          mod(iarray(l+i,j)-minint+1-1,80)+1),i=1,incr)
        enddo
      enddo
      return
8800  format (/, 2x, 31i4)
8900  format (1x,i3,1x, 124a1)
      end

      subroutine imatrx (iarray, im, istrt, iend, jstrt, jend, iform)

!=======================================================================
!      integer matrix print with various formats

!     inputs:

!     iarray = integer array to be printed
!     im     = the 1st dimension of array
!     istrt  = starting index along 1st dimension for plot
!     iend   = ending index along 1st dimension for plot
!     jstrt  = starting index along 2nd dimension for plot
!     jend   = ending index along 2nd dimension for plot
!              note: if jstrt and jend are negative then the vertical
!                   y-axis is inverted.
!     iform  = format designator: 1 => format (i1)
!                                 2 => format (i2)
!                                 3 => format (i3)
!                                 4 => format (i4)

!     output: print integer "iarray" with user format control
!=======================================================================

      include "stdunits.h"
      dimension iarray(im,1000)

!     choose the plotting domain

      js = min(abs(jstrt),abs(jend))
      je = max(abs(jstrt),abs(jend))
      is = min(abs(istrt),abs(iend))
      ie = max(abs(istrt),abs(iend))
      il = ie-is+1

      write (stdout,*) ' '
      if (iform .eq. 1) then

!       use I1 format to print the integer array

        inc=120
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,8800) (l+i+is-1,i=1,incr,4)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,8900) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 2) then

!       use I2 format to print the integer array

        inc=60
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9000) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9200) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 3) then

!       use I3 format to print the integer array

        inc=40
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9400) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9500) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      elseif (iform .eq. 4) then

!       use I4 format to print the integer array

        inc=30
        do l=0,il,inc
          incr = min(inc,il-l)
          write (stdout,9600) (l+i+is-1,i=1,incr,2)
          do jj=js,je
            if (jstrt .lt. 0 .or. jend .lt. 0) then
              j = je - (jj-js)
            else
              j = jj
            endif
            write (stdout,9700) j, (iarray(l+i+is-1,j),i=1,incr)
          enddo
        enddo
      endif
      return
8800  format (/, 2x, 30i4)
8900  format (1x,i3,1x, 120i1)
9000  format (/, 3x, 30i4)
9200  format (1x,i3,1x, 60i2)
9400  format (/,/,/,2x,20i6/)
9500  format (1x,i3,1x,40i3)
9600  format (/,/,/,1x,20i8/)
9700  format (1x,i3,1x,30i4)
      end

      subroutine matrix (array, irdim, istrt, im, jstrt, jm, scale)

!=======================================================================

!     matrix is a general two-dimensional array printing routine,
!     input:
!     array = the array to be printed
!     irdim = the 1st dimension of array
!     istrt = the 1st element of the 1st dimension to be printed
!     im    = the last element of the 1st dimension to be printed
!     jstrt = the 1st element of the 2nd dimension to be printed
!     jm    = the last element of the 2nd dimension to be printed
!             the 2nd dimension is printed in reverse order if both
!             jstrt & jm are negative
!     scale = a scaling factor by which array is divided before
!             printing.  (if this is zero, no scaling is done.)
!             if scale=0, 10 columns are printed across in e format
!             if scale>0, 20 columns are printed across in f format

!     output: print "array" as a matrix
!=======================================================================

      include "stdunits.h"
      parameter (c0=0.0, c1=1.0)
      dimension array(irdim,1000)

      if (jstrt*jm .lt. 0) then
        write (stdout,999)  jstrt, jm
        stop '=>matrix'
      endif

!     allow for inversion of 2nd dimension

      if (jm .lt. 0) then
        js   = -jm
        je   = -jstrt
        jinc = -1
      else
        js   = jstrt
        je   = jm
        jinc = 1
      endif

      if (scale .eq. c0) then

        do is=istrt,im,10
          ie = min(is + 9,im)
          write (stdout,9001) (i, i=is,ie)
          do l=js,je,jinc
            write (stdout,9002) l, (array(i,l),i=is,ie)
          enddo
          write (stdout,'(/)')
        enddo
      else
        scaler = c1/scale
        do is=istrt,im,20
          ie = min(is + 19,im)
          write (stdout,9003) (i, i=is,ie)
          do l=js,je,jinc
            write (stdout,9004) l, (array(i,l)*scaler,i=is,ie)
          enddo
          write (stdout,'(/)')
        enddo
      endif
      return
999   format (1x,'jstrt=',i5,' jm=',i5,' in matrix')
9001  format(10i13)
9002  format(1x,i2,10(1pe13.5))
9003  format(3x,20i6)
9004  format(1x,i3,1x,20f6.2)
      end

      subroutine scope (array, im, il, jl, text)

!=======================================================================
!     scope interrogates "array" for the min, max (with respective
!            locations), and simple unweighted average

!     inputs:
!     array = the array to be interrogated
!     im    = the inner dimension of "array"
!     il    = the number of points along the inner dimension to consider
!     jl    = the number of points along the inner dimension to consider
!     text  = descriptive text (up to 15 chars) to be printed

!     output: prints min, max (with locations), and average of "array"
!=======================================================================

      character(*) :: text
      dimension array(im,jl)
      umax  = array(1,1)
      umin  = array(1,1)
      iumax = 1
      jumax = 1
      iumin = 1
      jumin = 1
      sum    = 0.0
      do j=1,jl
        do i=1,il
          sum = sum + array(i,j)
          if (array(i,j) .gt. umax) then
            umax = array(i,j)
            iumax = i
            jumax = j
          endif
          if (array(i,j) .lt. umin) then
            umin = array(i,j)
            iumin = i
            jumin = j
          endif
        enddo
      enddo
      avg = sum/(il*jl)
      write (*,9200) text, il, jl, iumin, jumin, umin, iumax, jumax
     &,              umax, avg
      return
9200  format (1x,'Scope: ',a15,': il=',i4,', jl=',i4
     &,', min(',i4,',',i4,')=',g10.3,', max(',i4,',',i4,')=',g10.3
     &,', avg=',g10.3/)
      end

      subroutine sum1st (a, imt, jmt, text)

!=======================================================================
!     inputs:
!     a     = the array to be interrogated
!     imt   = the 1st dimension of "a"
!     imt   = the 2nd dimension of "a"
!     text  = descriptive text to be printed

!     output: sums the first index of array "a" for each 2nd index
!=======================================================================

      character(*) :: text
      dimension a(imt,jmt)
      print *,' '
      print *,text
      big = abs(a(1,1))
      imax = 0
      jmax = 0
      do j=1,jmt
        sum = 0.0
        do i=1,imt
          sum = sum + a(i,j)
          if (abs(a(i,j)) .gt. big) then
            big = abs(a(i,j))
            imax = i
            jmax = j
          endif
        enddo
        if (sum .ne. 0.0) then
          write (*,'(a,i4,a,e14.7)')
     &    '2nd index=',j,'. sum over 1st index =',sum
        endif
      enddo
      if (imax .ne. 0 .and. jmax .ne. 0) then
        write (*,*) ' biggest=',a(imax,jmax),' at i=',imax,' j=',jmax
      else
        write (*,*) ' ---field is a constant =', a(1,1)
      endif
      return
      end

      subroutine plot (array, im, istrt, iend, jstrt,  jend
     &,                  zmin, zmax, nbin, title)

!=======================================================================
!     "plot" contours "array" by dividing the array values into bins,
!     assigning a character to each bin, and printing the characters to
!     produce a contour map.
!     The inner (1st) dimension is plotted horizontally (x-axis) and the
!     outer (2nd) dimension is plotted vertically (y-axis) and can be
!     inverted.

!     inputs:
!     array = the array to be contoured
!     im    = the inner dimension of "array"
!     istrt = starting index along 1st dimension for plot
!     iend  = ending index along 1st dimension for plot
!     jstrt = starting index along 2nd dimension for plot
!     jend  = ending index along 2nd dimension for plot
!             note: if jstrt and jend are negative then the vertical
!                   y-axis is inverted.
!     zmin  = the minimum value to plot
!     zmax  = the maximum value to plot
!             note: if zmin=zmax then then min and max of "array"
!                   is used.
!     nbin  = the number of bins between zmin and zmax
!     title = descriptive text string

!     output: contours "array"
!=======================================================================

      parameter (maxbin=52, ncols=124)
      character(*) :: title
      dimension array(im,1000), bin(maxbin,3)
      character(maxbin) :: c, cc
      character(1) :: line(ncols)
      data c/'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'/
      save c

!     choose the plotting domain

      js = min(abs(jstrt),abs(jend))
      je = max(abs(jstrt),abs(jend))
      is = min(abs(istrt),abs(iend))
      ie = max(abs(istrt),abs(iend))
      il = ie-is+1

!     set the max and min

      if (zmax .eq. zmin) then
        big = array(is,js)
        sml = array(is,js)
        do j=js,je
          do i=is,ie
            if (array(i,j) .gt. big) big = array(i,j)
            if (array(i,j) .lt. sml) sml = array(i,j)
          enddo
        enddo
      else
        big = max(zmin,zmax)
        sml = min(zmin,zmax)
      endif

!     set up the discretization into bins

      write (*,*) ' '
      write (*,*) ' => contour bins for: ', title
      delta = (big-sml)/max(1,abs(nbin))
      limit = min(maxbin,max(1,abs(nbin)))
      do n=1,limit
        bin(n,1) = sml + (n-1)*delta
        bin(n,2) = sml + n*delta
        bin(n,3) = 0.5*(bin(n,1)+bin(n,2))
        cc(n:n)  = c(n:n)
        if (n .eq. 1) cc(n:n) = '-'
        if (n .ge. max(1,abs(nbin)) .or. n .eq. maxbin) cc(n:n) = '+'
        if (bin(n,1) .le. 0.0 .and. bin(n,2) .gt. 0.0) cc(n:n) = '.'
        if (n .lt. limit) then
          write (*,'(1x,a,a,g10.3,a,g10.3,a,g10.3,a)')  cc(n:n)
     &,   ' = ', bin(n,3)
     &,   ', (', bin(n,1),' <= x < ', bin(n,2),')'
        else
          write (*,'(1x,a,a,g10.3,a,g10.3,a,g10.3,a)') cc(n:n)
     &,   ' = ', bin(n,3)
     &,   ', (', bin(n,1),' <= x <=', bin(n,2),')'
        endif
      enddo
      if (abs(nbin) .gt. maxbin) then
        write (*,*)
     &    ' => Note: numbers larger than this are also plotted as +'
      endif

!     plot the character representation of the bins

      write (*,*) ' '
      do l=0,il,ncols
        incr = min(ncols,il-l)
        write (*,8800) (l+i+is-1,i=1,incr,4)
        do jj=js,je
          if (jstrt .lt. 0 .or. jend .lt. 0) then
            j = jj
          else
            j = je - (jj-js)
          endif
          do i=1,incr
            k = min(ifix((array(l+i+is-1,j) - sml)/delta) + 1, limit)
            line(i) = cc(k:k)
          enddo
          write (*,8900) j, (line(i),i=1,incr)
        enddo
      enddo
      return
8800  format (/, 2x, 31i4)
8900  format (1x,i3,1x, 124a1)
      end

      subroutine print_checksum (a, im, jm, text)
      dimension a(im,jm)
      character(*) :: text
      sum = checksum (a, im, jm)
      print *, text, sum
      return
      end

      function checksum (a, im, jm)

      implicit none

      integer i, im, j, jm
      real checksum, sum, a(im,jm)

      sum = 0.0
      do j=1,jm
        do i=1,im
          sum = sum + abs(a(i,j))
        enddo
      enddo
      checksum = sum
      return
      end

      subroutine wrufio (iounit, array, len)

!=======================================================================
!     write unformatted fortran i/o

!     input:

!     iounit = fortran unit number
!     array  = array to be written to "iounit"
!     len    = length of array

!     output: writes to unit "iounit"
!=======================================================================

      implicit none

      integer iounit, l, len
      real array(len)
# if defined O_io_single
      real (kind=4) :: array4(len)

      do l=1,len
        array4(l) = array(l)
      enddo
      write (iounit) array4
# else

      write (iounit) array
# endif

      return
      end

      subroutine rdufio (iounit, array, len)

!=======================================================================
!     read unformatted fortran i/o

!     input:

!     iounit = fortran unit number
!     array  = array to be read from "iounit"
!     len    = length of array

!     output: none
!=======================================================================

      implicit none

      integer iounit, len
      real array(len)

      read (iounit) array

      return
      end

      subroutine tranlon (c, ic, il, jl, t, cx, fx, ifl, tx)

!-----------------------------------------------------------------------
!     grid interpolators in MOM require the grid coordinates to be
!     monotonic. when the prime meridian lies within the model grid,
!     global datasets (e.g. Scripps topography) must be translated in
!     longitude to remove the jump in
!     longitudes (358.5 359.5, 0.5, 1.5) across the meridian before
!     interpolating to the model grid. This is only of concern in
!     limited domain grids (e.g. Atlantic basin) that contain the
!     prime meridian.

!     translate longitudes  "cx" to "tx" so that tx(i) i=1..ic
!     completely encloses model longitudes  fx(i) i=1..ifl
!     note that "tx" may extend beyond 360 degrees to contain "fx".
!     use same mapping to translate data in "c"

!     input:
!     c  = original data array
!     t  = temp array for translating data
!     cx = original data longitudes
!     tx = translated data longitudes
!     fx = model longitudes

!     output
!     c  = translated data array
!-----------------------------------------------------------------------

      implicit none

      integer iw, ic, il, jl, ifl, indp, i, im1, j

      real c(ic,jl), t(ic), tx(ic), cx(ic), fx(ifl)

!-----------------------------------------------------------------------
!     find the index of the 1st model grid point on the data grid
!-----------------------------------------------------------------------

      iw = indp (fx(1), cx, ic)
      if (cx(iw) .gt. fx(1)) iw = max(1,iw-1)

!-----------------------------------------------------------------------
!     translate data longitudes so that tx(1) = cx(iw), tx(2) = cx(iw+1)
!-----------------------------------------------------------------------

      do i=1,ic
        tx(i) = cx(mod(i+iw-2,il) + 1)
        im1   = max(1,i-1)
        if (tx(i) .lt. tx(im1)) tx(i) = tx(i) + 360.0
      enddo
      if (fx(ifl) .gt. tx(ic)) then
        write (6,997) iw, ic, ifl
        write (6,998) 'tx= ',(tx(i),i=1,ic)
        write (6,998) 'fx= ',(fx(i),i=1,ifl)
        stop
      endif

!-----------------------------------------------------------------------
!     translate data to match translated longitudes
!-----------------------------------------------------------------------

      do j=1,jl
        do i=1,ic
          t(i) = c(mod(i+iw-2,il) + 1,j)
        enddo
        do i=1,ic
          c(i,j) = t(i)
        enddo
      enddo

      return
997   format (1x, ' ===>  tx(ic) < fx(ifl) in tranlon. iw=',i6,
     1         ' il=',i6,' ifl=',i6)
998   format (1x,a4,(5x,10e11.4))
      end

      subroutine areatot (data, dmsk, tot)

!=======================================================================
!     calculate the area weighted total of data

!     input:
!       data = data to be totalled
!       dmsk = data mask
!     output:
!       tot = area weighted total of the data
!=======================================================================

      implicit none

      integer i, j, land
      real sl, el, fx, tot

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"

      real data(imt,jmt), dmsk(imt,jmt)

      tot = 0.0

      do j=2,jmtm1
        fx = cst(j)*dyt(j)
        do i=2,imtm1
          tot = tot + data(i,j)*dxt(i)*fx*dmsk(i,j)
        enddo
      enddo

      return
      end

      subroutine areaavg (data, dmsk, avg)

!=======================================================================
!     calculate the area weighted average of data

!     input:
!       data = data to be averaged
!       dmsk = data mask
!     output:
!       avg = area weighted average of the data
!=======================================================================

      implicit none

      integer i, j
      real fx, tarea, tdata, area, avg

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"

      real data(imt,jmt), dmsk(imt,jmt)

      tarea = 0.0
      tdata = 0.0
      avg = 0.0

      do j=2,jmtm1
        fx = cst(j)*dyt(j)
        do i=2,imtm1
          area = dxt(i)*fx*dmsk(i,j)
          tarea = tarea + area
          tdata = tdata + data(i,j)*area
        enddo
      enddo

      if (tarea .ne. 0.0) avg = tdata/tarea

      return
      end
