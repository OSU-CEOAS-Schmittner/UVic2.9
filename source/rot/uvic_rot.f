      subroutine set_land (data, value, mask, imt, jmt, k)

!=======================================================================
!     set land to value
!=======================================================================

      implicit none

      real data(imt,jmt), value
      integer mask(imt,jmt), i, j, k, imt, jmt

      if (k .ge. 0) then
        do j=1,jmt
          do i=1,imt
            if (mask(i,j) .lt. k) data(i,j) = value
          enddo
        enddo
      else
        do j=1,jmt
          do i=1,imt
            if (mask(i,j) .ge. -k) data(i,j) = value
          enddo
        enddo
      endif

      return
      end

      subroutine load_vctr (u, v, data, imt, jmt)

!=======================================================================
!     load vector components into a single array
!=======================================================================

      implicit none

      real u(imt,jmt), v(imt,jmt), data(imt,jmt,2)
      integer i, j, imt, jmt

      do j=1,jmt
        do i=1,imt
          data(i,j,1) = u(i,j)
          data(i,j,2) = v(i,j)
        enddo
      enddo

      return
      end

      subroutine unload_vctr (u, v, data, imt, jmt)

!=======================================================================
!     unload vector components to separate arrays
!=======================================================================

      implicit none

      real u(imt,jmt), v(imt,jmt), data(imt,jmt,2)
      integer i, j, imt, jmt

      do j=1,jmt
        do i=1,imt
          u(i,j) = data(i,j,1)
          v(i,j) = data(i,j,2)
        enddo
      enddo

      return
      end

      subroutine intrp_vert (wrka, wrkb, wtb, bad, data, imt, jmt)

!=======================================================================
!     interpolate in the vertical
!=======================================================================

      implicit none

      real wrka(imt,jmt), wrkb(imt,jmt), data(imt,jmt), bad, wt, wtb
      integer i, j, imt, jmt

      do j=1,jmt
        do i=1,imt
          wt = wtb
          if (abs(wrkb(i,j)) .ge. abs(bad)) wt = 0
          data(i,j) = wrka(i,j)*(1.0 - wt) + wrkb(i,j)*wt
        enddo
      enddo

      return
      end

      subroutine extrap (data, bad, xt, mask, imt, jmt, km)

!=======================================================================
!     replace missing values with extrapolated ones. looks along a
!     level and if land locked looks up.
!=======================================================================

      implicit none

      real data(imt,jmt,km), wrk(imt,jmt,km), xt(imt), sum, bad, abad, d
      integer mask(imt,jmt), ir(imt*jmt*km), jr(imt*jmt*km)
     &,       kr(imt*jmt*km), i, j, k, ngood, l, m, n, nrp, nr
     &,       imt, jmt, km, km1

!     find indices of missing data

      abad = abs(bad)
      nr = 0
      nrp = 0
      do k=1,km
        do j=1,jmt
          do i=1,imt
            wrk(i,j,k) = data(i,j,k)
            if (abs(data(i,j,k)) .ge. abad) then
              if (mask(i,j) .ge. k .or. km .eq. 1) then
                nr = nr + 1
                ir(nr) = i
                jr(nr) = j
                kr(nr) = k
              endif
            endif
          enddo
        enddo
      enddo

      do while (nrp .ne. nr)

        nrp = nr
        nr = 0
        do n=1,nrp
          i = ir(n)
          j = jr(n)
          k = kr(n)
          sum = 0.0
          ngood = 0
          do l=-1,1
            do m=-1,1
              d = data(min(imt,max(1,i+m)), min(jmt,max(1,j+l)), k)
              if (abs(d) .lt. abad) then
                ngood = ngood + 1
                sum = sum + d
              endif
            enddo
          enddo
          if (ngood .ne. 0) then
            wrk(i,j,k) = sum/ngood
            if (xt(imt-1) - 359.999 .gt. xt(1)) then
              wrk(1,j,k) = wrk(imt-1,j,k)
              wrk(imt,j,k) = wrk(2,j,k)
            endif
          else
            nr = nr + 1
            ir(nr) = i
            jr(nr) = j
            kr(nr) = k
          endif
        enddo

        do k=1,km
          do j=1,jmt
            do i=1,imt
              data(i,j,k) = wrk(i,j,k)
            enddo
          enddo
        enddo

        if (xt(imt-1) - 359.999 .gt. xt(1)) then
          do k=1,km
            do j=1,jmt
              data(1,j,k) = data(imt-1,j,k)
              data(imt,j,k) = data(2,j,k)
            enddo
          enddo
        endif

      enddo

      do k=1,km
        km1 = max(1, k-1)
        do j=1,jmt
          do i=1,imt
            if (abs(data(i,j,k)) .gt. abad)
     &        data(i,j,k) = data(i,j,km1)
          enddo
        enddo
      enddo

      return
      end

      subroutine extrap2 (data, bad, xt, imt, jmt)

!=======================================================================
!     replace missing values with extrapolated ones.
!=======================================================================

      implicit none

      real data(imt,jmt), wrk(imt,jmt), xt(imt), sum, bad, abad, d
      integer ir(imt*jmt), jr(imt*jmt)
      integer i, j, ngood, l, m, n, nrp, nr, imt, jmt

!     find indices of missing data

      abad = abs(bad)
      nr = 0
      nrp = 0
      do j=1,jmt
        do i=1,imt
          wrk(i,j) = data(i,j)
          if (abs(data(i,j)) .ge. abad) then
            nr = nr + 1
            ir(nr) = i
            jr(nr) = j
          endif
        enddo
      enddo

      do while (nrp .ne. nr)

        nrp = nr
        nr = 0
        do n=1,nrp
          i = ir(n)
          j = jr(n)
          sum = 0.0
          ngood = 0
          do l=-1,1
            do m=-1,1
              d = data(min(imt,max(1,i+m)), min(jmt,max(1,j+l)))
              if (abs(d) .lt. abad) then
                ngood = ngood + 1
                sum = sum + d
              endif
            enddo
          enddo
          if (ngood .ne. 0) then
            wrk(i,j) = sum/ngood
            if (xt(imt-1) - 359.999 .gt. xt(1)) then
              wrk(1,j) = wrk(imt-1,j)
              wrk(imt,j) = wrk(2,j)
            endif
          else
            nr = nr + 1
            ir(nr) = i
            jr(nr) = j
          endif
        enddo

        do j=1,jmt
          do i=1,imt
            data(i,j) = wrk(i,j)
          enddo
        enddo

        if (xt(imt-1) - 359.999 .gt. xt(1)) then
          do j=1,jmt
            data(1,j) = data(imt-1,j)
            data(imt,j) = data(2,j)
          enddo
        endif

      enddo

      return
      end

      subroutine rot_intrp_vctr (g, xg, yg, ig, jg, r, xr, yr, ir, jr
     &,                          phir, thetar, psir, bad, near)

!=======================================================================
!     interpolate vector data from an geographic data grid to a
!     rotated model grid

!     input
!     g    = vector on geographic data grid
!     xg   = longitude of data points on geographic data grid
!     yg   = latitude of data points on geographic data grid
!     ig   = number of longitudes in geographic data grid
!     jg   = number of latitudes in geographic data grid
!     xr   = longitude of points on rotated model grid
!     yr   = latitude of points on rotated model grid
!     ir   = number of longitudes in rotated model grid
!     jr   = number of latitudes in rotated model grid
!     phir, thetar, psir = Euler angles defining rotation
!     bad  = bad data value (make negative to suppress warnings)
!     near = flag to find nearest neighbour instead of interpolation

!     output
!     r   = vector on rotated model grid

!     internal
!     (rln,rlt) = (longitude,latitude) in rotated coordinates
!     (gln,glt) = (longitude,latitude) in geographic coordinates
!     vm        = vector magnitude
!     x         = x component at (rln,rlt)
!     y         = y component at (rln,rlt)
!     a         = vector angle
!=======================================================================

      implicit none

      real g(ig,jg,2), xg(ig), yg(jg), r(ir,jr,2), xr(ir), yr(jr)
      real phir, thetar, psir, bad

      integer i, j, ig, jg, ir, jr, near

      real (kind=8) :: rad, abad, dphir, dthetar, dpsir, rlt, rln
      real (kind=8) :: glt, gln, x, y, vm, a, d0, d1, d90, d180, d360

      d0   = 0.D+00
      d1   = 1.D+00
      d90  = 90.D+00
      d180 = 180.D+00
      d360 = 360.D+00
      rad = acos(-d1)/d180

!     interpolate vector components as scalers on rotated model grid

      call rot_intrp_sclr (g(1,1,1), xg, yg, ig, jg, r(1,1,1), xr, yr
     &,                          ir, jr, phir, thetar, psir, bad, near)
      call rot_intrp_sclr (g(1,1,2), xg, yg, ig, jg, r(1,1,2), xr, yr
     &,                          ir, jr, phir, thetar, psir, bad, near)

!     correct vector direction

!     convert some variables to double precision
      abad = abs(bad)
      dphir = phir
      dthetar = thetar
      dpsir = psir
      do j=1,jr
        do i=1,ir
          x = r(i,j,1)
          y = r(i,j,2)
          vm = sqrt(x**2 + y**2)
          if (abs(vm) .lt. abad) then
            rlt = yr(j)
            if (vm .gt. d0 .and. abs(rlt) .lt. d90) then
              rln = xr(i)
!             keep the range of rln between -180 and 180 degrees
              rln = mod(rln,d360)
              if (rln .lt. -d180) rln = rln + d360
              if (rln .gt. d180) rln = rln - d360
              call drot_angle (rlt, rln, dphir, dthetar, dpsir, a)
!             add correction to original vector angle
              vm = sqrt(x**2 + y**2)
              if (abs(x) .lt. abs(y)) then
                if (y .lt. d0) then
                  a = a - acos(max(min(x/vm, d1),-d1))
                else
                  a = a + acos(max(min(x/vm, d1), -d1))
                endif
              else
                if (x .lt. d0) then
                  a = a - asin(max(min(y/vm, d1), -d1)) + d180*rad
                else
                  a = a + asin(max(min(y/vm, d1), -d1))
                endif
              endif
              x = vm*cos(a)
              y = vm*sin(a)
            else
              x = d0
              y = d0
            endif
          else
            x = abad
            y = abad
          endif
          r(i,j,1) = x
          r(i,j,2) = y
        enddo
      enddo

      return
      end

      subroutine drot_angle (rlt, rln, dphir, dthetar, dpsir, a)

!=======================================================================
!     find local grid rotation angle

!     input
!     (rln,rlt) = (longitude,latitude) in rotated coordinates
!     dphir, dthetar, dpsir = Euler angles defining rotation

!     output
!     a   = angle of local grid rotation

!     internal
!     (rln,rlt) = (longitude,latitude) in rotated coordinates
!     (gln,glt) = (longitude,latitude) in geographic coordinates
!     rltp      = rotated grid latitude of rotated pole
!     rlnp      = rotated grid longitude of rotated pole
!     gltp      = geographic grid latitude of rotated pole
!     glnp      = geographic grid longitude of rotated pole
!     dgp       = distance to geographic pole
!     drp       = distance to rotated pole
!     dbp       = distance between geographic and rotated pole
!=======================================================================

      implicit none

      real (kind=8) :: rad, dphir, dthetar, dpsir, rlt, rln, glt, gln
      real (kind=8) :: rltp, rlnp, gltp, glnp, dgp, drp, dbp, a
      real (kind=8) :: d0, d1, d90, d180, d360

      d0   = 0.D+00
      d1   = 1.D+00
      d90  = 90.D+00
      d180 = 180.D+00
      d360 = 360.D+00
      rad = acos(-d1)/d180

      call drotate (rlt, rln, -dpsir, -dthetar, -dphir, glt, gln)

!     find pointers to the furthest pole on the geographic grid.
      rltp = d90
      rlnp = rln
!     if the point is in the north use the south pole instead.
      if (glt .gt. d0) rltp = -d90
      call drotate (rltp, rlnp, dphir, dthetar, dpsir, gltp, glnp)
!     distance between point and the geographic pole
      dgp = abs(rltp - glt)*rad
!     distance between point and rotated pole
      drp = abs(rltp - rlt)*rad
!     distance between the rotated and geographic poles
      dbp = abs(rltp - gltp)*rad

!     find the angle between pointers with the law of cosines.
      a = (cos(dbp) - cos(dgp)*cos(drp))/(sin(dgp)*sin(drp))
      a = acos(max(min(a, d1), -d1))

!     determine the sign by checking the offset longitudes.
      if (glnp - rlnp .gt. d180) rlnp = rlnp + d360
      if (rlnp - glnp .gt. d180) glnp = glnp + d360
      if (glnp + d360 .gt. rlnp + d360) a = -a

!     change sign if we used the south pole.
      if (glt .ge. 0.) a = -a

      return
      end

      subroutine rot_intrp_sclr (g, xg, yg, ig, jg, r, xr, yr, ir, jr
     &,                          phir, thetar, psir, bad, near)

!=======================================================================
!     interpolate scaler data from an geographic data grid to a
!     rotated model grid

!     input
!     psir, thetar, phir = Euler angles defining rotation
!     g    = scaler on geographic data grid
!     xg   = longitude of data points on geographic data grid
!     yg   = latitude of data points on geographic data grid
!     ig   = number of longitudes in on geographic data grid
!     jg   = number of latitudes in on geographic data grid
!     xr   = longitude of points on rotated model grid
!     yr   = latitude of points on rotated model grid
!     ir   = number of longitudes in rotated model grid
!     jr   = number of latitudes in rotated model grid
!     bad  = bad data value (make negative to suppress warnings)
!     near = flag to find nearest neighbour instead of interpolation

!     output
!     r   = scaler on rotated model grid

!     internal
!     (rln,rlt) = (longitude,latitude) in rotated coordinates
!     (gln,glt) = (longitude,latitude) in geographic coordinates
!     xg(iw) = point on the geographic grid to the west of (gln,glt)
!     xg(ie) = point on the geographic grid to the east of (gln,glt)
!     yg(js) = point on the geographic grid to the south of (gln,glt)
!     yg(jn) = point on the geographic grid to the north of (gln,glt)
!=======================================================================

      implicit none

      real g(ig,jg), xg(ig), yg(jg), r(ir,jr), xr(ir), yr(jr), gln, glt
     &,    gln_max, glt_max, gln_min, glt_min, wt, wtn, wts, wte, wtw
     &,    wtne, wtse, wtnw, wtsw, phir, thetar, psir, bad, abad, del
     &,    epsln, tg(ig,jg), txg(ig), onep

      integer i, j, ig, jg, ir, jr, iw, ie, jn, js, lt_err, ln_err, near
     &,       jstrt, jend, istrt, iend, indp, ii, izero

      epsln = 1.e-10
      glt_min = 90.
      glt_max = -90.
      gln_min = 360.
      gln_max = -360.
      ln_err = 0
      lt_err = 0
      abad = abs(bad)
      iend = ig
      istrt = 1
      onep = 1.000001

!     look for extra points (up to 2 for cyclic conditions)
      if (abs(xg(ig) - 360.0 - xg(1)) .lt. 0.0001) iend = ig-1
      if (abs(xg(ig) - 360.0 - xg(2)) .lt. 0.0001) then
        istrt = 2
        iend = ig-1
      endif

!     find longitude points of data within interval [0., 360.]
      izero = istrt
      do i=1,ig
        txg(i) = mod(xg(i),360.0)
        if (txg(i) .lt. 0.0) txg(i) = txg(i) + 360.0
        if (i .gt. istrt .and. i .le. iend) then
          if (txg(i-1) .gt. txg(i)) izero = i
        endif
      enddo

!     shift the data if needed
      do i=istrt,iend
        ii = i + izero - istrt
        if ( ii .gt. iend) ii = ii - iend + istrt - 1
        txg(i) = mod(xg(ii),360.0)
        if (txg(i) .lt. 0.0) txg(i) = txg(i) + 360.0
        do j=1,jg
          tg(i,j) = g(ii,j)
        enddo
      enddo

!     set cyclic conditions if needed
      if (iend .lt. ig) then
        txg(ig) = txg(istrt) + 360.0
        do j=1,jg
          tg(ig,j) = g(istrt,j)
        enddo
      endif
      if (istrt .gt. 1) then
        txg(1) = txg(iend) - 360.0
        do j=1,jg
          tg(1,j) = g(iend,j)
        enddo
      endif

!     find latitude points of data within interval [-90., 90.]
      jstrt = 1
      do j=2,jg
        if (yg(j-1) .lt. -90. .and. yg(j) .ge. -90.) jstrt = j
      enddo
      jend = jg
      do j=2,jg
        if (yg(j) .le. 90.) jend = j
      enddo

!     interpolate data to model grid

      do j=1,jr
        do i=1,ir
          call rotate (yr(j), xr(i), -psir, -thetar, -phir, glt, gln)
          if (gln .lt. 0.) gln = gln + 360.
          if (gln .ge. 360.) gln = gln - 360.
          glt_min = min(glt,glt_min)
          glt_max = max(glt,glt_max)
          gln_min = min(gln,gln_min)
          gln_max = max(gln,gln_max)

          iw = indp (gln, txg(istrt), iend - istrt + 1) + istrt - 1
          if (txg(iw) .gt. gln) iw = iw - 1
          ie = iw + 1
          if (iw .ge. istrt .and. ie .le. iend) then
            del = txg(ie) - txg(iw)
              wtw = (txg(ie) - gln)/del
          else
!     east or west of the last data value. this could be because a
!     cyclic condition is needed or the dataset is too small. in either
!     case apply a cyclic condition
            ln_err = 1
            iw = iend
            ie = istrt
            del = txg(ie) + 360. + epsln - txg(iw)
            if (txg(ie) .ge. gln) then
                wtw = (txg(ie) - gln)/del
            else
              wtw = (txg(ie) + 360. + epsln - gln)/del
            endif
          endif
          wte = 1. - wtw

          js = indp (glt, yg(jstrt), jend - jstrt + 1) + jstrt - 1
          if (yg(js) .gt. glt) js = max(js - 1,jstrt)
          jn = min(js + 1,jend)
          if (yg(jn) .ne. yg(js) .and. yg(js) .le. glt) then
            wts = (yg(jn) - glt)/(yg(jn) - yg(js))
          else
!     north or south of the last data value. this could be because a
!     pole is not included in the data set or the dataset is too small.
!     in either case extrapolate north or south
            lt_err = 1
            wts = 1.
          endif
          wtn = 1. - wts

!     check for weighting error
          if (abs(wtn) .gt. onep .or. abs(wts) .gt. onep .or.
     &        abs(wte) .gt. onep .or. abs(wtw) .gt. onep) then
            print*, 'Weighting error: i, j, wtn, wts, wte, wtw '
     &,                 i,j,wtn,wts,wte,wtw
            stop 'rot_intrp_sclr'
          endif

!     set weights to zero for bad data values and normalize the good

          wtne = wtn*wte
          if (abs(tg(ie,jn)) .ge. abad) wtne = 0.
          wtse = wts*wte
          if (abs(tg(ie,js)) .ge. abad) wtse = 0.
          wtnw = wtn*wtw
          if (abs(tg(iw,jn)) .ge. abad) wtnw = 0.
          wtsw = wts*wtw
          if (abs(tg(iw,js)) .ge. abad) wtsw = 0.
          if (near .ne. 0) then
            wt = max(wtne, wtse, wtnw, wtsw)
            if (wtne .eq. wt) then
              wtne = 1.
              wtse = 0.
              wtnw = 0.
              wtsw = 0.
            endif
            if (wtse .eq. wt) then
              wtne = 0.
              wtse = 1.
              wtnw = 0.
              wtsw = 0.
            endif
            if (wtnw .eq. wt) then
              wtne = 0.
              wtse = 0.
              wtnw = 1.
              wtsw = 0.
            endif
            if (wtsw .eq. wt) then
              wtne = 0.
              wtse = 0.
              wtnw = 0.
              wtsw = 1.
            endif
          endif
          wt = wtne + wtse + wtnw + wtsw

          if (wt .gt. 0) then
            r(i,j) = (tg(ie,jn)*wtne + tg(ie,js)*wtse
     &             + tg(iw,jn)*wtnw + tg(iw,js)*wtsw)/wt
          else
            r(i,j) = abad*2.
          endif

        enddo
      enddo

      if (ln_err .eq. 1 .and. bad .gt. 0.0) then
        write (*,'(/,(1x,a))')
     &    '==> Warning: the geographic data set does not extend far   '
     &,   '             enough east or west - a cyclic boundary       '
     &,   '             condition was applied. check if appropriate   '
        write (*,'(/,(1x,a,2f8.2))')
     &    '    data required between longitudes:', gln_min, gln_max
     &,   '      data set is between longitudes:', txg(istrt), txg(iend)
      endif

      if (lt_err .eq. 1 .and. bad .gt. 0.0) then
        write (*,'(/,(1x,a))')
     &    '==> Warning: the geographic data set does not extend far   '
     &,   '             enough north or south - extrapolation from    '
     &,   '             the nearest data was applied. this may create '
     &,   '             artificial gradients near a geographic pole   '
        write (*,'(/,(1x,a,2f8.2))')
     &    '    data required between latitudes:', glt_min, glt_max
     &,   '      data set is between latitudes:', yg(jstrt), yg(jend)
      endif

      if (abs(xr(ir) - 360.0 - xr(1)) .lt. 0.0001) then
        do j=1,jr
          r(ir,j) = r(1,j)
        enddo
      endif
      if (abs(xr(ir) - 360.0 - xr(2)) .lt. 0.0001) then
        do j=1,jr
          r(1,j) = r(ir-1,j)
          r(ir,j) = r(2,j)
        enddo
      endif

      return
      end

      subroutine rotate (glt, gln, phir, thetar, psir, rlt, rln)

!     double precision translator for drotate

      implicit none

      real glt, gln, phir, thetar, psir, rlt, rln
      real (kind=8) :: dglt, dgln, dphir, dthetar, dpsir, drlt, drln

      dglt = glt
      dgln = gln
      dphir = phir
      dthetar = thetar
      dpsir = psir
      drlt = rlt
      drln = rln
      call drotate (dglt, dgln, dphir, dthetar, dpsir, drlt, drln)
      rlt = drlt
      rln = drln

      return
      end

      subroutine drotate (glt, gln, phir, thetar, psir, rlt, rln)

!=======================================================================
!     subroutine drotate takes a geographic latitude and longitude and
!     finds the the equivalent latitude and longitude on a rotated grid.
!     when going from a geographic grid to a rotated grid, all of the
!     defined rotation angles given to drotate by the calling program
!     are positive, but when going from a rotated grid back to the
!     geographic grid, the calling program must reverse the angle order
!     (phir and psir are switched) and all of the angles made negative.

!     the first rotation angle phir is defined as a rotation about the
!     original z axis. the second rotation angle thetar is defined as a
!     rotation about the new x axis. the final rotation angle psir is
!     defined as a rotation about the new z axis. these rotation angles
!     are just the Euler angles as defined in "classical mechanics"
!     Goldstein (1951).
!=======================================================================

      implicit none

      real (kind=8) :: glt, gln, phir, thetar, psir, phis, thetas, rlt
      real (kind=8) :: rln, gy, gx, gz, rx, ry, rz, d0, d1, d90, d180
      real (kind=8) :: d360, rad

!     g...  = geographic value
!     r...  = rotated value
!     ...lt = latitude (or equivalent spherical coordinate)
!     ...ln = longitude (or equivalent spherical coordinate)
!     ...x  = x coordinate
!     ...y  = y coordinate
!     ...z  = z coordinate
!     psir, thetar, phir = Euler angles defining rotation

      d0   = 0.D+00
      d1   = 1.D+00
      d90  = 90.D+00
      d180 = 180.D+00
      d360 = 360.D+00

!     define rad for conversion to radians.
      rad = acos(-d1)/d180

!     convert latitude and longitude to spherical coordinates
      thetas = gln
      if (thetas .gt. d180) thetas = thetas - d360
      if (thetas .lt. -d180) thetas = thetas + d360
      thetas = thetas*rad
      phis = (d90 - glt)*rad

!     translate point into Cartesian coordinates for rotation.
      gx = sin(phis)*cos(thetas)
      gy = sin(phis)*sin(thetas)
      gz = cos(phis)

!     rotate the point (gx, gy, gz) about the z axis by phir then the x
!     axis by thetar and finally about the z axis by psir.

      rx = gx*(cos(psir)*cos(phir) - cos(thetar)*sin(phir)*sin(psir)) +
     &     gy*(cos(psir)*sin(phir) + cos(thetar)*cos(phir)*sin(psir)) +
     &     gz*sin(psir)*sin(thetar)

      ry = gx*(-sin(psir)*cos(phir) - cos(thetar)*sin(phir)*cos(psir)) +
     &     gy*(-sin(psir)*sin(phir) + cos(thetar)*cos(phir)*cos(psir)) +
     &     gz*(cos(psir)*sin(thetar))

      rz = gx*(sin(thetar)*sin(phir)) + gy*(-sin(thetar)*cos(phir)) +
     &     gz*(cos(thetar))

!     convert rotated point back to spherical coordinates

!     check for rounding error (arccos(x): abs(x) must be .le. 1)
      rz = min(rz, d1)
      rz = max(rz, -d1)
      rlt = acos(rz)
!     if point is at a pole set rotated longitude equal to initial.
      if (rlt .le. d0 .or. rlt .ge. d180*rad) then
        rln = thetas
      else
!       if rln lies between -135 and -45 or between 45 and 135 degrees
!       it is more accurate to use an arccos calculation.
        if (abs(rx/sin(rlt)) .lt. cos(45.D+00*rad)) then
          rln = acos(max(min(rx/sin(rlt), d1), -d1))
!         arccos will give rln between 0 and 180 degrees.  if the point
!         is negative in y, rln must be equal to negative rln.
          if (ry .lt. d0) rln = -rln
        else
!         if rln lies between -45 and 45 or between 135 and -135 degrees
!         it is more accurate to use an arcsin calculation.
          rln = asin(max(min(ry/sin(rlt), d1), -d1))
!         arcsin will give rln between -90 and 90 degrees. if the point
!         is negative in x, rln must be equal to 180 degrees minus rln.
          if (rx .lt. d0) rln = d180*rad - rln
        endif
      endif

!     convert back to degrees of latitude and longitude.
      rlt = d90 - rlt/rad
      rln = rln/rad
      if (rln .gt. d180) rln = rln - d360
      if (rln .le. -d180) rln = rln + d360

      return
      end

      integer function indp (value, array, ia)

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

      integer stdin, stdout, stderr
      parameter (stdin = 5, stdout = 6, stderr = 6)

      integer i, l, ia
      real array(ia), value

      do i=2,ia
        if (array(i) .lt. array(i-1)) then
           print*,
     &   ' => Error: array must be monotonically increasing in "indp"'
     &,  '           when searching for nearest element to value=',value
          print*, '           array(i) < array(i-1) for i=',i
          print*, '           array(i) for i=1..ia follows:'
          do l=1,ia
            print*, 'i=',l, ' array(i)=',array(l)
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

      dimension array(len)
      integer l, len, iounit
      real array
      real (kind=4) :: array4(len)

      do l=1,len
        array4(l) = array(l)
      enddo
      write (iounit) array4
      return
      end
