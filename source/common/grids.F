      subroutine gcoord (maxlen, imt, jmt, km, dxtdeg, dytdeg, dxudeg
     &,                  dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)

!=======================================================================

!                      G R I D    C O O R D I N A T E S

!     Construct grid point coordinates and resolution

!     input:

!     maxlen = maximum number of grid cells in latitude, longitude,
!              and depth

!              set grid specifications in USER INPUT section.

!     output:

!     imt    = number of longitudes
!     jmt    = number of latitudes
!     km     = number of depths
!     dxtdeg = width of "t" grid cells (degrees)
!     dytdeg = height of "t" grid cells (degrees)
!     dxudeg = width of "u" grid cells (degrees)
!     dyudeg = height of "u" grid cells (degrees)
!     dzt    = thickness of "t" grid cells (cm)
!     dzw    = thickness of "w" grid cells (cm)
!     xt     = longitude at centres of "t" grid cells (degrees)
!     xu     = longitude at centres of "u" grid cells (degrees)
!     yt     = latitude at centres of "t" grid cells (degrees)
!     yu     = latitude at centres of "u" grid cells (degrees)
!     zt     = depth at centres of "t" grid cells (centimetres)
!     zw     = depth at centres of "u" grid cells (centimetres)
!=======================================================================

      implicit none

      integer ib(10), ic(10), maxbounds, imt, jmt, jrow, km
      integer i, k, ncase, num, maxlen, iou
      parameter (maxbounds=11)

      character(120) :: fname, new_file_name

      logical exists

      real p5, tolr, dxubar, dyubar, dzwbar, cksum, checksum, c0, c1
      real c100

      include "stdunits.h"

      real xt(imt), yt(jmt), xu(imt), yu(jmt), zw(km), zt(km)
      real dxtdeg(imt), dytdeg(jmt), dzt(km), dxudeg(imt), dyudeg(jmt)
      real dzw(0:km)

!     set some constants

      p5 = 0.5

      ncase = 0
      ncase = ncase + 1
      c0 = 0.
      c1 = 1.
      c100 = 100.
      fname = new_file_name ("G_grid.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'gcoord in grids.f'
      endif
      dxtdeg(:) = c0
      dytdeg(:) = c0
      dxudeg(:) = c0
      dyudeg(:) = c0
      dzt(:) = c0
      dzw(:) = c0
      xt(:) = c0
      yt(:) = c0
      xu(:) = c0
      yu(:) = c0
      zt(:) = c0
      zw(:) = c0
      ib(:) = 1
      call openfile (fname, iou)
      ic(:) = imt
      call getvara ('G_dxt', iou, imt, ib, ic, dxtdeg, c1, c0)
      call getvara ('G_dxu', iou, imt, ib, ic, dxudeg, c1, c0)
      call getvara ('longitude', iou, imt, ib, ic, xt, c1, c0)
      call getvara ('longitude_V', iou, imt, ib, ic, xu, c1, c0)
      ic(1) = jmt
      call getvara ('G_dyt', iou, jmt, ib, ic, dytdeg, c1, c0)
      call getvara ('G_dyu', iou, jmt, ib, ic, dyudeg, c1, c0)
      call getvara ('latitude', iou, jmt, ib, ic, yt, c1, c0)
      call getvara ('latitude_V', iou, jmt, ib, ic, yu, c1, c0)
      ic(1) = km
      call getvara ('G_dzt', iou, km, ib, ic, dzt, c100, c0)
      call getvara ('G_dzw', iou, km, ib, ic, dzw(1:km), c100, c0)
      call getvara ('depth', iou, km, ib, ic, zt, c100, c0)
      call getvara ('depth_W', iou, km, ib, ic, zw, c100, c0)

      if (zt(km) .gt. 8000.e2) then
!       assume depths are in cm
        print*, "==> Warning: converting depths to cm in gcoord"
        zt(1:km) = zt(1:km)/100.
        dzt(1:km) = dzt(1:km)/100.
        zw(1:km) = zw(1:km)/100.
        dzw(1:km) = dzw(1:km)/100.
      endif
      dzw(0)  = zt(1)
      dzw(km) = zw(km) - zt(km)

!-----------------------------------------------------------------------
!     Check if the "t" grid resolution is an average of the
!     "u" cell resolution. This insures more accurate advection of
!     tracers within a stretched grid.
!-----------------------------------------------------------------------

      num  = 0
      tolr = 1.e-5
      write (stdout,'(/)')
      do i=2,imt-1
        dxubar = p5*(dxudeg(i) + dxudeg(i-1))
        if (abs(dxubar-dxtdeg(i)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta x at i=',i
     &,   ' is not an average of adjacent "u" cell delta x`s'
        endif
      enddo

      do jrow=2,jmt-1
        dyubar = p5*(dyudeg(jrow) + dyudeg(jrow-1))
        if (abs(dyubar-dytdeg(jrow)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta y at jrow=',jrow
     &,   ' is not an average of adjacent "u" cell delta y`s'
        endif
      enddo

      tolr = 1.e0
      do k=2,km-1
        dzwbar = p5*(dzw(k) + dzw(k-1))
        if (abs(dzwbar-dzt(k)) .gt. tolr) then
          num = num + 1
          write (stdout,'(a,i5,a)')
     &    '=>Warning: "t" cell delta z at k=',k
     &,   ' is not an average of adjacent "w" cell delta z`s'
        endif
      enddo

      if (num .ne. 0) then
        write (stdout,'(/a/a/a/a//a,a/)')
     &  '==>Warning, At the above locations, advection of tracers is'
     &, 'not as accurate as it could be. If you are reading in your own'
     &, 'grid or constructing a grid as in MOM 1, we assume you want to'
     &, 'define the grid this way and we let you proceed from here...'
     &, 'Please read ALL the information in the USER INPUT section to '
     &, 'understand what this means'
      endif

!-----------------------------------------------------------------------
!     Print all grid coordinates
!-----------------------------------------------------------------------

      write (stdout
     &,'(//,40x,a,//,a,g14.7,a,/a/,a,g14.7,a/a,/,a,g14.7,a)')
     &  ' Grid Point Coordinate details: '
     &, ' The western edge of the 2nd "t" cell is at longitude:'
     &,   xu(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
     &, ' The southern edge of the 2nd "t" cell is at latitude:'
     &,   yu(1),' (deg)',' (the 1st "t" cell is a boundary cell)'
      write (stdout,'(/,a,g14.7,a/a/,a,g14.7,a/a/,a,g14.7,a/)')
     &  ' The western edge of the 1st "u" cell is at longitude:', xt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
     &, ' The southern edge of the 1st "u" cell is at latitude:', yt(1)
     &, ' (deg)',' (the 1st "u" cell is a boundary point)'
      write (stdout,9103) km
      write (stdout,9002) (zt(k),k=1,km)
      write (stdout,9104) km
      write (stdout,9002) (zw(k),k=1,km)
      write (stdout,9105) jmt
      write (stdout,9001) (yt(jrow),jrow=1,jmt)
      write (stdout,9106) jmt
      write (stdout,9001) (yu(jrow),jrow=1,jmt)
      write (stdout,9107) imt
      write (stdout,9001) (xt(i),i=1,imt)
      write (stdout,9108) imt
      write (stdout,9001) (xu(i),i=1,imt)

#if defined O_symmetry

!-----------------------------------------------------------------------
!     insure that yu(jmt-1) = 0.0 (equator) when using symmetry
!-----------------------------------------------------------------------

      if (yu(jmt-1) .ne. 0.0) then
        write (stdout,*) '=> Error: yu(jmt-1) must = 0.0 for symmetry'
        stop '=>gcoord'
      endif
#endif

!---------------------------------------------------------------------
!     compute a grid checksum
!---------------------------------------------------------------------

      cksum = 0.0
      cksum = cksum + checksum (xt, imt, 1)
      cksum = cksum + checksum (yt, jmt, 1)
      cksum = cksum + checksum (zt, km, 1)
      cksum = cksum + checksum (xu, imt, 1)
      cksum = cksum + checksum (yu, jmt, 1)
      cksum = cksum + checksum (zw, km, 1)
      cksum = cksum + checksum (dxtdeg, imt, 1)
      cksum = cksum + checksum (dytdeg, jmt, 1)
      cksum = cksum + checksum (dxudeg, imt, 1)
      cksum = cksum + checksum (dyudeg, jmt, 1)
      cksum = cksum + checksum (dzt, km, 1)
      cksum = cksum + checksum (dzw, km+1, 1)
      write (stdout,'(/)')
      write (stdout,*) 'Grid checksum = ',cksum
      write (stdout,'(/)')
      return
9001  format (1x,10f10.4)
9002  format (1x,10f10.2)
9103  format (/,' Depth to "t" & "u" grid points (cm): zt(k) k=1,',i3)
9104  format (/,' Depth to "w" grid points (cm): zw(k) k=1,',i3)
9105  format (/,' Latitude of "t" points (deg): yt(j) j=1,',i4)
9106  format (/,' Latitude of "u" points (deg): yu(j) j=1,',i4)
9107  format (/,' Longitude of "t" points (deg): xt(i) i=1,',i4)
9108  format (/,' Longitude of "u" points (deg): xu(i) i=1,',i4)
      end

      subroutine gcell (maxlen, n_bounds, bounds, d_bounds, nbpts
     &,                 num, deltat, deltau, stretch)

!=======================================================================

!              G R I D   C E L L   C O N S T R U C T I O N

!     A domain is composed of one or more regions:
!     Build "num" "t"  cells with resolution "deltat(n) n=1,num"
!     within the domain composed of regions bounded by "bounds".
!     Also construct "num" "u"  cells of resolution "deltau(n) n=1,num"
!     with the relation between "t" and "u" cells given by:
!     deltat(n) = 0.5*(deltau(n-1) + deltau(n))
!     Resolution may be constant or smoothly varying within each
!     region AND there must be an integral number of grid cells within
!     each region. The domain is the sum of all regions.

!     inputs:

!     maxlen   = maximum length of "deltat" and "deltau"
!     n_bounds = number of bounds needed to define the regions
!     bounds   = latitude, longitude, or depth at each bound
!     d_bounds = delta (resolution) at each of the "bounds"
!     nbpts    = number of extra boundary cells to add to the domain.
!                (usually one at the beginning and end)
!     stretch  = stretching factor for last region (should only be used
!                in the vertical to provide increased stretching of grid
!                points. "stretch" = 1.0 gives no increased stretching.
!                "stretch" = 1.2 gives increased stretching...etc

!     outputs:

!     num    = total number of grid cells within the domain
!     deltau = resolution of "u" grid cells: n=1,num
!     deltat = resolution of "t" grid cells: n=1,num
!=======================================================================

      implicit none

      integer maxlen, n_bounds, num, l, m, n, i, nbpts

      real p5, pi, avg_res, stretch, chg_res, tol, wid, an, sum
      real del

      include "stdunits.h"

      real deltat(maxlen), deltau(maxlen), d_bounds(n_bounds)
      real bounds(n_bounds)

!     Set some constants

      p5 = 0.5
      pi = 4.0*atan(1.0)

!     Do all regions, one at a time, to construct the domain

      num  = 1
      do l=1,n_bounds-1

        write (stdout,'(2x,a,i2,a,g14.7,a,g14.7,a,g14.7,a,g14.7,a)')
     & ' region # ',l,'  going from ',bounds(l),' (res=',d_bounds(l)
     &,') to ',  bounds(l+1),' (res=',d_bounds(l+1),')'

!       avg_res = average resolution of "t" cells within region
!       chg_res = change in resolution across the region
!       wid     = width of region
!       tol     = tolerance for fitting "t" cells within region width

!       provide for stretching last region if needed

        if (l .eq. n_bounds-1) then
          avg_res = p5*(d_bounds(l) + stretch*d_bounds(l+1))
          chg_res = (stretch*d_bounds(l+1) - d_bounds(l))
        else
          avg_res = p5*(d_bounds(l) + d_bounds(l+1))
          chg_res = (d_bounds(l+1) - d_bounds(l))
        endif

        tol = 1.e-5
        wid = abs(bounds(l+1) - bounds(l))
        an  = wid/avg_res
        m   = nint(an)

!       Calculate resolution of "u" cells: "deltau"
!       "u" grid points will be centered in these cells
!       n = number of "t" cells fitting within the region boundaries
!       note: "sum" initially discounts half of the "u" cells widths
!       at the boundaries

        sum = 0.5*d_bounds(l) - 0.5*d_bounds(l+1)
        n   = 0
        do i = 1,100000
          del = avg_res - p5*chg_res*cos((pi/m)*i)
          if (sum + del .le. wid*(1.0 + tol)) then
            sum = sum + del
            if (num+i-1 .gt. maxlen) then
              write (stdout,*) "=>Error: maxlen exceeded in gcell. "
     &,                        " ...increase size of maxlen"
              stop
            endif
            deltau(num+i-1) = del
            n = n + 1
          else
            go to 100
          endif
        enddo

100     continue
        num = num + n
      enddo

!     adjust "num" to reflect the total number of cells contained in
!     all regions

      num = num - 1

      do i=1,num

!       build resolution for "T" cells: "deltat". Note that
!       variable resolution (stretched grid) implies "T" points are
!       off centre

        if (i .eq. 1) then
          deltat(i) = p5*(d_bounds(1) + deltau(i))
        else
          deltat(i) = p5*(deltau(i) + deltau(i-1))
        endif

      enddo

!     add boundary points if needed

      if (nbpts .ne. 0) then
        do i=num,1,-1
          deltat(i+1) = deltat(i)
          deltau(i+1) = deltau(i)
        enddo
        deltat(1)     = deltat(2)
        deltau(1)     = d_bounds(1)
        deltat(num+2) = deltat(num+1)
        deltau(num+2) = deltau(num+1)
        num           = num + 2
      endif
      return
      end

      subroutine grids

!=======================================================================
!     set up a staggered "B" grid for MOM and compute grid related
!     variables
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer ib(10), ic(10), maxlen, imt2, jmt2, km2, jrow, i, k, iou
      integer ip2, j, jp2, jp1, jm1, kp2, kp1, km1

      real degtcm, tiny

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
#if !defined O_implicitvmix || defined O_isopycmix
      include "accel.h"
#endif
      include "scalar.h"
#if defined O_mom
      include "hmixc.h"
      include "vmixc.h"
#endif

!-----------------------------------------------------------------------
!     set some constants
!-----------------------------------------------------------------------

      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      degtcm = radius/radian

!-----------------------------------------------------------------------
!     calculate coordinates for "t" and "u" grid cells.
!-----------------------------------------------------------------------

      maxlen = max(imt,jmt,km)
      imt2 = imt
      jmt2 = jmt
      km2 = km
      call gcoord (maxlen, imt2, jmt2, km2, dxtdeg, dytdeg, dxudeg
     &,            dyudeg, dzt, dzw, xt, xu, yt, yu, zt, zw)

!-----------------------------------------------------------------------
!     verify that the number of grid points match the number set in
!     the parameter statement in "size.h".
!-----------------------------------------------------------------------

      call size_check (imt2, jmt2, km2, 'sub grids', 'stop')

!-----------------------------------------------------------------------
!     convert grid resolution to cm
!-----------------------------------------------------------------------

      do jrow=1,jmt
        dyt(jrow) = dytdeg(jrow)*degtcm
        dyu(jrow) = dyudeg(jrow)*degtcm
      enddo

      do i=1,imt
        dxt(i) = dxtdeg(i)*degtcm
        dxu(i) = dxudeg(i)*degtcm
      enddo
#if defined O_cyclic
      dxt(1)   = dxt(imt-1)
      dxt(imt) = dxt(2)
      dxu(1)   = dxu(imt-1)
      dxu(imt) = dxu(2)
#endif

!-----------------------------------------------------------------------
!     compute all quantities derived from the grid spacings
!-----------------------------------------------------------------------

      do k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo

      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)

      do k=1,km
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo

      do k=1,km
#if defined O_implicitvmix || defined O_isopycmix
        dztur(k) = c1/(dzw(k-1)*dzt(k))
        dztlr(k) = c1/(dzw(k)*dzt(k))
#endif
        dztr(k)  = c1/dzt(k)
      enddo

      tiny = 1.e-20
      do jrow=1,jmt
        dytr(jrow)  = c1/dyt(jrow)
        dyt2r(jrow) = p5/dyt(jrow)
        dyt4r(jrow) = p25/dyt(jrow)
        dyur(jrow)  = c1/dyu(jrow)
        dyu2r(jrow) = p5/dyu(jrow)
        dyu4r(jrow) = p25/dyu(jrow)
        phi(jrow)   = yu(jrow)/radian
        phit(jrow)  = yt(jrow)/radian
        cst(jrow)   = cos(phit(jrow))
        csu(jrow)   = cos(phi (jrow))
        sine(jrow)  = sin(phi(jrow))
        if (cst(jrow)  .eq. 0.0) then
          print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          cst(jrow)  = tiny
        endif
        if (csu(jrow)  .eq. 0.0) then
          print '(/a,e14.7,a,i4,/a)'
     &  ,' Warning: setting cst(jrow) = ',tiny, ' for jrow =',jrow
     &,  '          to prevent division by zero at the pole'
          csu(jrow)  = tiny
        endif
        cstr(jrow)     = c1/cst(jrow)
        csur(jrow)     = c1/csu(jrow)
        tng(jrow)      = sine(jrow)/csu(jrow)
        cstdytr(jrow)  = c1/(cst(jrow)*dyt(jrow))
        cstdyt2r(jrow) = cstdytr(jrow)*p5
        csudyur(jrow)  = c1/(csu(jrow)*dyu(jrow))
        csudyu2r(jrow) = p5/(csu(jrow)*dyu(jrow))
        cst_dytr(jrow) = cst(jrow)/dyt(jrow)
        csu_dyur(jrow) = csu(jrow)/dyu(jrow)

      enddo

      do i=1,imt
        dxtr(i)  = c1/dxt(i)
        dxt2r(i) = p5/dxt(i)
        dxt4r(i) = p25/dxt(i)
        dxur(i)  = c1/dxu(i)
        dxu2r(i) = p5/dxu(i)
        dxu4r(i) = p25/dxu(i)
      enddo

      do i=2,imtm1
        dxmetr(i) = c1/(dxt(i) + dxt(i+1))
      enddo

      do i=1,imt
        duw(i) = (xu(i) - xt(i))*degtcm
      enddo
      do i=1,imtm1
        due(i) = (xt(i+1) - xu(i))*degtcm
      enddo
#if defined O_cyclic
      due(imt) = due(2)
#else
      due(imt) = due(imtm1)
#endif

      do jrow=1,jmt
        dus(jrow) = (yu(jrow) - yt(jrow))*degtcm
      enddo

      do jrow=1,jmtm1
        dun(jrow) = (yt(jrow+1) - yu(jrow))*degtcm
      enddo
      dun(jmt) = dun(jmtm1)

#if defined O_mom
# if !defined O_implicitvmix || defined O_isopycmix

!     for convection code, compute values needed to include effects
!     of tracer timestep acceleration on effective layer thicknesses

      do k=1,km
        dztxcl(k) = dzt(k)/dtxcel(k)
      enddo

      do k=1,kmm1
        dzwxcl(k) = c1/(dztxcl(k)+dztxcl(k+1))
      enddo
      dzwxcl(km) = c0

# endif
# if defined O_quicker

!     coefficients for 3rd order advective scheme

      do i=2,imt-1
        if (i .lt. imt-1) then
          ip2 = i+2
        else
          ip2 = 3
        endif
        quick_x(i,1) = c2*dxt(i+1)/(dxt(i+1)+dxt(i))
        quick_x(i,2) = c2*dxt(i  )/(dxt(i+1)+dxt(i))

        curv_xp(i,1) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i)+dxt(i+1)))
        curv_xp(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i)+dxt(i+1))
     &                                  *(dxt(i-1)+dxt(i)))
        curv_xp(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i-1)+2.0*dxt(i)+dxt(i+1))*(dxt(i-1)+dxt(i)))

        curv_xn(i,1) = c2*dxt(i)*dxt(i+1)/
     &             ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i+1)+dxt(ip2)))
        curv_xn(i,2) =-c2*dxt(i)*dxt(i+1)/((dxt(i+1)+dxt(ip2))
     &                                  *(dxt(i)+dxt(i+1)))
        curv_xn(i,3) = c2*dxt(i)*dxt(i+1)/
     &               ((dxt(i)+2.0*dxt(i+1)+dxt(ip2))*(dxt(i)+dxt(i+1)))
      enddo
      call setbcx (quick_x, imt, 2)
      call setbcx (curv_xp, imt, 3)
      call setbcx (curv_xn, imt, 3)

      do j=1,jmt
        jp2 = min(j+2,jmt)
        jp1 = min(j+1,jmt)
        jm1 = max(j-1,1)
        quick_y(j,1) = c2*dyt(jp1)/(dyt(jp1)+dyt(j))
        quick_y(j,2) = c2*dyt(j  )/(dyt(jp1)+dyt(j))

        curv_yp(j,1) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(j)+dyt(jp1)))
        curv_yp(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(j)+dyt(jp1))
     &                                  *(dyt(jm1)+dyt(j)))
        curv_yp(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(jm1)+2.0*dyt(j)+dyt(jp1))*(dyt(jm1)+dyt(j)))

        curv_yn(j,1) = c2*dyt(j)*dyt(jp1)/
     &             ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(jp1)+dyt(jp2)))
        curv_yn(j,2) =-c2*dyt(j)*dyt(jp1)/((dyt(jp1)+dyt(jp2))
     &                                  *(dyt(j)+dyt(jp1)))
        curv_yn(j,3) = c2*dyt(j)*dyt(jp1)/
     &               ((dyt(j)+2.0*dyt(jp1)+dyt(jp2))*(dyt(j)+dyt(jp1)))
      enddo

      do k=1,km
        kp2 = min(k+2,km)
        kp1 = min(k+1,km)
        km1 = max(k-1,1)
        quick_z(k,1) = c2*dzt(kp1)/(dzt(kp1)+dzt(k))
        quick_z(k,2) = c2*dzt(k  )/(dzt(kp1)+dzt(k))

        curv_zp(k,1) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(k)+dzt(kp1)))
        curv_zp(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(k)+dzt(kp1))
     &                                  *(dzt(km1)+dzt(k)))
        curv_zp(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(km1)+2.0*dzt(k)+dzt(kp1))*(dzt(km1)+dzt(k)))

        curv_zn(k,1) = c2*dzt(k)*dzt(kp1)/
     &             ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(kp1)+dzt(kp2)))
        curv_zn(k,2) =-c2*dzt(k)*dzt(kp1)/((dzt(kp1)+dzt(kp2))
     &                                  *(dzt(k)+dzt(kp1)))
        curv_zn(k,3) = c2*dzt(k)*dzt(kp1)/
     &               ((dzt(k)+2.0*dzt(kp1)+dzt(kp2))*(dzt(k)+dzt(kp1)))
      enddo

# endif
#endif
      fname = new_file_name ("G_grid.nc")
      ib(:) = 1
      ic(:) = imt
      ic(2) = jmt
      call openfile (fname, iou)
      call getvara ('G_latT', iou, imt*jmt, ib, ic, tlat, c1, c0)
      call getvara ('G_lonT', iou, imt*jmt, ib, ic, tlon, c1, c0)
      call getvara ('G_latU', iou, imt*jmt, ib, ic, ulat, c1, c0)
      call getvara ('G_lonU', iou, imt*jmt, ib, ic, ulon, c1, c0)

      return
      end
