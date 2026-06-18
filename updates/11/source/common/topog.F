      subroutine topog (kmt, kmu, map, xt, yt, zt, xu, yu, zw, imt2
     &,                 jmt2, km2, sg_bathy)

!-----------------------------------------------------------------------
!     construct the topographic mask "kmt" which determines the geometry
!     and topography by defining the number of vertical levels on model
!     "t" points.

!-----------------------------------------------------------------------
!     input:

!     xt  = longitudes of model grid "t" points (degrees)
!     yt  = latitudes of model grid "t" points (degrees)
!     zt  = depths of box "u/v" points on model grid (cm)
!     xu  = longitudes of model grid "u" points (degrees)
!     yu  = latitudes of model grid "u" points (degrees)
!     zw  = depths of box bottoms on model grid "t" boxes (cm)
!     imt = number of model longitudes
!     jmt = number of model latitudes
!     km  = maximum number of levels

!     output:

!     kmt = number of vertical levels on "t" points
!           0 is a land point.
!           n implies n boxes between the surface and bottom
!             where (n >= kmt_min) & (n <= km)
!     map = array indicating land masses and their perimeters
!           0 is a mid ocean point
!           n > 0 is a land mass
!           n < 0 is the ocean perimeter of land mass |n|

!-----------------------------------------------------------------------

      implicit none

      character(40) :: cktext
      character(120) :: iotext
      character(60) :: expnam
      character(32) :: stamp
      character(120) :: fname, new_file_name

!     kmt_min is the minimum permitted number of levels at an ocean
!     "T" cell.  kmt_min must be at least 2.
      integer kmt_min, imt2, jmt2, km2, ncase, iou, j, i, n
      integer nerror, jrow, npass, k, nisle, linewidth
      parameter (kmt_min = 2)

      real c0, c1, cksum1, checksumi, cksum2

      logical exists

      include "stdunits.h"
      include "size.h"
      include "isleperim.h"

      integer map(imt,jmt), iperm(maxipp), jperm(maxipp),nippts(mnisle)
      integer iofs(mnisle), kmt(imt,jmt), kmu(imt,jmt), ib(10), ic(10)

      real tmpij(imt-2,jmt-2), tmpijk(imt-2,jmt-2,km)
      real xt(imt), yt(jmt), zt(km), xu(imt), yu(jmt), zw(km)
      real sg_bathy(imt,jmt,km)

!     the variable auto_kmt_changes tells whether any changes have been
!     made to the kmt field as a result of define options.
!     nchanges counts how many changes in each group.

      auto_kmt_changes = .false.
      nchanges = 0

!-----------------------------------------------------------------------
!     check that grid sizes in argument list and file "size.h" agree
!-----------------------------------------------------------------------

      call size_check (imt2, jmt2, km2, 'topog', 'stop')

      ncase = 0

      write (stdout,'(//,t37,a,/)')
     & 'T O P O G R A P H Y    G E N E R A T I O N'
      write (stdout,'(a,i4,/)') 'kmt_min = ', kmt_min

      ncase = ncase + 1
      c0 = 0.
      c1 = 1.
      fname = new_file_name ("G_kmt.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'topog in topog.f'
      endif
      ib(:) = 1
      ic(:) = imt-2
      ic(2) = jmt-2
      call openfile (fname, iou)
      if (.not. exists) then
        print*, "==>  Warning: G_kmt data does not exist."
      else
        call getvara ('G_kmt', iou, ic(1)*ic(2), ib, ic, tmpij, c1, c0)
      endif
      kmt(2:imt-1,2:jmt-1) = tmpij(1:imt-2,1:jmt-2)

      nerror = 0
      do i=1,imt
        do j=1,jmt
          if (kmt(i,j) .lt. 0 .or. kmt(i,j) .gt. km) then
            nerror = nerror + 1
          endif
        enddo
      enddo
      if (nerror .gt. 0) then
        write (stdout,'(a,i4,a,a,i4,/)')
     &    '==>ERROR: There are ', nerror, ' depth values in the file '
     &,   '"kmt.dta" that are outside the range 0 <= kmt(i,j) <= km = '
     &,   km
        stop '==>topog'
      endif
      write (stdout,'(/a/)') '==> The kmt field has been imported.'

!-----------------------------------------------------------------------
!     apply boundary conditions
!-----------------------------------------------------------------------

      call kmtbc (kmt, imt, jmt)
      call area_volume (kmt, xu, yu, zw)
      cksum1 = checksumi (kmt, imt, jmt, 'Original "kmt" checksum =')

!-----------------------------------------------------------------------
!     iterate on changes to "kmt" until all are consistent
!-----------------------------------------------------------------------

      npass = 0
      n_del_kmt = 0
1000  continue
      npass = npass + 1

!-----------------------------------------------------------------------
!     fill isolated cells or potholes at all levels in kmt
!     filled kmt array is just the maximum of the surrounding kmu levels
!     symmetry conditions are automatic, if set on kmu
!-----------------------------------------------------------------------

      do i=1,imt-1
        do j=1,jmt-1
          kmu(i,j) = min(kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
        enddo
      enddo
      n = 0
      write (*,'(/a)') 'Searching for isolated T cells...'
      do i=2,imt-1
        do j=2,jmt-1
          k = max (kmu(i,j), kmu(i-1,j), kmu(i,j-1), kmu(i-1,j-1))
          if (kmt(i,j) .ne. k) n = n + 1
          kmt(i,j) = k
        enddo
      enddo
      if (n .gt. 0) then
        write (*,*) '-> Found ',n,' and filled them in.'
      else
        write (*,*) '-> None Found.'
      endif

      call kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     USER INPUT
!     other changes to "kmt" may be made here.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     ALL USER MODIFICATIONS MUST PRECEDE THESE FINAL CHECKS
!-----------------------------------------------------------------------

      call kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     limit the minimum number of levels. kmt_min should be >= 2
!-----------------------------------------------------------------------

      call min_depth (kmt, kmt_min, xu, yu, zw)
      call kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     test for "island perimeter violations"
!     channels between different land masses must be at least 2 ocean
!     cells wide to prevent conflicting values of the stream function
!     on these island perimeter points.
!-----------------------------------------------------------------------

      write (stdout,'(/a/)')
     & 'Searching for and correcting island PERIMETER VIOLATIONS...'
      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle
     &,               imt, jmt, km, mnisle, maxipp, xu, yu, zw)
      call kmtbc (kmt, imt, jmt)
      call area_volume (kmt, xu, yu, zw)

      write (cktext,'(a,i2,a)') 'Pass # ',npass
     &,                         ' complete.  kmt checksum ='
      cksum2 = checksumi (kmt, imt, jmt, cktext)
      if (cksum1 .ne. cksum2) then
        cksum1 = cksum2
        if (npass .gt. 15) then
          print *,'=>Error: changes to "kmt" may be incompatible'
          print *,'         arbitrarily stopping after 15 passes'
          stop
        endif
        go to 1000
      endif
      print '(/a/)', 'All changes have now been made to kmt'

!-----------------------------------------------------------------------
!     end of iteration on "kmt" changes
!-----------------------------------------------------------------------

      call showmap (map, imt, jmt, linewidth)

!-----------------------------------------------------------------------
!     show the kmt field
!-----------------------------------------------------------------------

      write (stdout,'(/,20x,a/)') ' The "kmt" field follows:'
      call iplot (kmt, imt, imt, jmt)
      call area_volume (kmt, xu, yu, zw)

!---------------------------------------------------------------------
!     compute a topography checksum
!---------------------------------------------------------------------

      cksum1 = checksumi (kmt, imt, jmt, 'Final "kmt" checksum =')

!     construct depth arrays associated with "u" cells

      do j=1,jmt
        kmu(imt,j) = 0
      enddo
      do i=1,imt
        kmu(i,jmt) = 0
      enddo
      do j=1,jmt-1
        do i=1,imt-1
           kmu(i,j) = min (kmt(i,j), kmt(i+1,j), kmt(i,j+1)
     &,                    kmt(i+1,j+1))
        enddo
      enddo
#if defined O_cyclic
      do j=1,jmt
        kmu(imt,j) = kmu(2,j)
      enddo
#endif
#if defined O_symmetry
      do i=1,imt
        kmu(i,jmt) = kmu(i,jmt-2)
      enddo
#endif
      sg_bathy(:,:,:) = 0.
      fname = new_file_name ("G_subgrid_bathy.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
         print*, "Error => ", trim(fname), " does not exist."
         stop 'topog in topog.f'
      endif
      ib(:) = 1
      ic(:) = imt-2
      ic(2) = jmt-2
      ic(3) = km
      call openfile (fname, iou)
      if (.not. exists) then
         print*, "==> Warning: G_subgrid_bathy data does not exist."
      else
         print*, "read G_subgrid_bathy"
         call getvara ('SUBGRID_BATHY', iou, ic(1)*ic(2)*ic(3), ib, ic
     &,             tmpijk, c1, c0)
         print*, "finish read G_subgrid_bathy"
      endif
      sg_bathy(2:imt-1,2:jmt-1,:) = tmpijk(1:imt-2,1:jmt-2,:)
c     make sure sg_bathy = 1 at the bottom (required for conservation)
      do i=2,imt-1
         do j=2,jmt-1
            if (kmt(i,j) .gt. 0) sg_bathy(i,j,kmt(i,j)) = 1.
         enddo
      enddo
      print *,'sg_bathy(kmt)=1'
#if defined O_cyclic
      do j=2,jmt-1
         sg_bathy(1,j,:) = sg_bathy(imt-1,j,:)
         sg_bathy(imt,j,:) = sg_bathy(2,j,:)
      enddo
      print *,'finished setting up sg_bathy'
#else
      print *,'stop in topog.F:'
      print *,'boundary condition not implemented for sg_bathy'
      stop
#endif

      return
      end

      function checksumi (kmt, im, jm, text)

      implicit none

      integer i, im, jm, jrow, kmt(im,jm)
      character(*) :: text
      real checksumi, cksum

      cksum = 0.0
      do jrow=1,jm
        do i=1,im
          cksum = cksum + i*jrow*kmt(i,jrow)
        enddo
      enddo
      print *, text, cksum
      checksumi = cksum
      return
      end

      subroutine area_volume(kmt, xu, yu, zw)

!-----------------------------------------------------------------------
!     compute surface area and volume of ocean ("t" cells and "u" cells)
!     (note that areas are defined at each level)
!-----------------------------------------------------------------------

      implicit none

      integer k, i, jrow
      real ocnp, tcellv, degrad

      include "size.h"
      include "scalar.h"

      integer kmt(imt,jmt)

      real xu(imt), yu(jmt), zw(km), tcella(km), cst(jmt), dxt(imt)
      real dyt(jmt)

      do k=1,km
        tcella(k) = 0
      enddo
      ocnp   = 0
      tcellv = 0
      pi = 4.0 * atan(1.0)
      degrad = pi / 180.0
      do i=2,imt-1
        dxt(i) = radius * degrad * (xu(i) - xu(i-1))
      enddo
      do jrow = 2,jmt-1
        cst(jrow) = cos(0.5*(yu(jrow)+yu(jrow-1))*degrad)
        dyt(jrow) = radius * degrad * (yu(jrow) - yu(jrow-1))
      enddo

      do jrow=2,jmt-1
        do i=2,imt-1
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
              tcella(k) = tcella(k) +
     &           cst(jrow)*dxt(i)*dyt(jrow)
            enddo
            tcellv = tcellv + cst(jrow)*dxt(i)*dyt(jrow)*zw(kmt(i,jrow))
            ocnp   = ocnp + float(kmt(i,jrow))
          endif
        enddo
      enddo

      print '(/)'
      print '(a,f8.0)',   'number of ocean T cells =', ocnp
      print '(a,e15.6,a)','surface area (T cells)  =',tcella(1),' cm**2'
      print '(a,e15.6,a)', 'ocean volume (T cells)  =', tcellv,' cm**3'
      print '(/)'

      return
      end

      subroutine kmtbc (kmt, imt, jmt)

!-----------------------------------------------------------------------
!     set lateral boundary conditions on kmt
!     set cyclic, solid wall and symmetry conditions on "t" grid
!-----------------------------------------------------------------------

      implicit none

      integer kmt(imt,jmt), imt, jmt, jrow, i

#if defined O_cyclic

!     set cyclic conditions on eastern and western boundary

      do jrow=1,jmt
        kmt(1,jrow)   = kmt(imt-1,jrow)
        kmt(imt,jrow) = kmt(2,jrow)
      enddo
#else

!     set solid wall conditions on eastern and western boundary
!     (if "cyclic" is not enabled then option "solid_walls" is assumed)

      do jrow=1,jmt
        kmt(1,jrow)   = 0
        kmt(imt,jrow) = 0
      enddo
#endif
#if defined O_symmetry

!     set symmetry conditions at "u" row jmt-1

      do i=1,imt
        kmt(i,jmt) = kmt(i,jmt-1)
      enddo
#else

!     close the domain at the northernmost row

      do i=1,imt
        kmt(i,jmt) = 0
      enddo
#endif

!     close the domain at the southernmost row

      do i=1,imt
        kmt(i,1)   = 0
      enddo

      return
      end

      subroutine setkmt (kmt, xt, yt, imt, jmt, alat1, slon1, elon1
     &,                  alat2, slon2, elon2, num)

!-----------------------------------------------------------------------
!     set the topography mask "kmt(i,j)" = "num" within the area of
!     the parallelogram bounded by vertices:
!     (alat1,slon1), (alat1,elon1), (alat2,slon2), & (alat2,elon2)

!     inputs:

!     xt = longitudes of "t" points in degrees
!     yt = latitudes of "t" points in degrees
!     imt = number of model longitudes
!     jmt = number of model latitudes
!     alat1 = southern latitude of parallelogram (degrees)
!     slon1 = starting longitude of southern edge of parallelogram (deg)
!     elon1 = ending longitude of southern edge of parallelogram (deg)
!     alat2 = northern latitude of parallelogram (degrees)
!     slon2 = starting longitude of northern edge of parallelogram (deg)
!     elon2 = ending longitude of northern edge of parallelogram (deg)
!     num   = number of vertical levels

!     outputs:

!     kmt   = mask of vertical levels on model "t" points
!-----------------------------------------------------------------------

      implicit none

      integer imt, jmt, j1, indp, j2, js, je, i1, i2, is1, ie1
      integer is2, ie2, is, ie, jrow, i, num, kmt(imt,jmt)

      real alat1, alat2, slon1, elon1, slon2, elon2, c1, rdj
      real xt(imt), yt(jmt)

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
!     the nudging of 1.e-5 is to insure that the test case resolution
!     generates the same topography and geometry on various computers.

      c1 = 1.0
      if (js .eq. je) then
        rdj = c1
      else
        rdj = c1/(je-js)
      endif
      do jrow=js,je
        do i=is,ie
          kmt(i,jrow) = num
        enddo
        is = nint(rdj*((jrow-js)*is2 + (je-jrow)*is1) + 1.0e-5)
        ie = nint(rdj*((jrow-js)*ie2 + (je-jrow)*ie1) + 1.0e-5)
      enddo
      return
      end

      subroutine min_depth (kmt, kmt_min, xu, yu, zw)

!     limit the minimum number of levels. kmt_min should be >= 2

      implicit none

      integer n, i, jrow, kmt_min, kmt_shallow

      include "size.h"
      include "stdunits.h"

      integer kmt(imt,jmt)
      real xu(imt), yu(jmt), zw(km)

      write (stdout,'(/a/)')
     & 'Searching for and correcting minimum level violations'
      n = 0
      do i=2,imt-1
        do jrow=jmt-1,2,-1
          if (kmt(i,jrow) .ne. 0 .and. kmt(i,jrow) .lt. kmt_min) then
            n = n + 1
            kmt_shallow = kmt(i,jrow)
            if (zw(kmt(i,jrow)) .lt. 0.5*zw(kmt_min)) then
              kmt(i,jrow) = 0
            else
              kmt(i,jrow) = kmt_min
            endif
          endif
        enddo
      enddo
      if (n .gt. 0) then
        write (stdout,'(a,i5,a/)')  '->Modified', n,' shallow cells'
      else
        write (stdout,'(a/)') '->No modifications needed'
      endif
      return
      end
