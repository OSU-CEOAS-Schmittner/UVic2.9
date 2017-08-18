      subroutine ppmixi (vmixset)

#if defined O_mom && defined O_ppvmix
      implicit none

      integer io, j, i, k, jrow

      logical vmixset, error

      real dzmin, p25, c0, extlat

      include "size.h"
      include "accel.h"
      include "coord.h"
      include "iounit.h"
      include "scalar.h"
      include "stdunits.h"
      include "vmixc.h"

!=======================================================================
!     Initialization for the Pacanowski/Philander vertical mixing scheme
!     Pacanowski & Philander (JPO vol 11, #11, 1981).

!     input:
!       dzt    = thickness of vertical levels (cm)
!       km     = number of vertical levels
!       yt     = latitude of grid points (deg)
!       jmt    = number of latitudes
!       dtxcel = time step accelerator as a function of level
!       dtts   = density time step (sec)
!       dtuv   = internal mode time step (sec)
!       vmixset= logical to determine if a vertical mixing scheme was
!                chosen

!     output:
!       wndmix = min value for mixing at surface to simulate high freq
!                wind mixing (if absent in forcing). (cm**2/sec)
!       fricmx = maximum mixing (cm**2/sec)
!       diff_cbt_back = background "diff_cbt" (cm**2/sec)
!       visc_cbu_back = background"visc_cbu" t(cm**2/sec)
!       diff_cbt_limit = largest "diff_cbt" (cm**2/sec)
!       visc_cbu_limit = largest "visc_cbu" (cm**2/sec)
!       vmixset= true
!=======================================================================

      namelist /ppmix/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
     &,                visc_cbu_limit, diff_cbt_limit

      write (stdout,'(/,20x,a,/)')
     & 'P P V M I X    I N I T I A L I Z A T I O N'

!-----------------------------------------------------------------------
!     initialize variables (all mixing units are cm**2/sec.)
!-----------------------------------------------------------------------

      wndmix    = 10.0
      fricmx    = 50.0
      diff_cbt_back =  0.1
      visc_cbu_back =  1.0
      dzmin     = 1.e10
      p25       = 0.25
      c0        = 0.0

# if defined O_implicitvmix

!     simulate convective adjustment with large mixing coefficient
!     limits

      visc_cbu_limit = fricmx
      diff_cbt_limit = 1.0e6
# else

!     in regions of gravitational instability set mixing limits to the
!     maximum consistant with the "cfl" criterion. convective adjustment
!     will also act on the instability.

      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx
# endif

!-----------------------------------------------------------------------
!     provide for namelist over-ride of above settings + documentation
!-----------------------------------------------------------------------

      call getunit (io, 'control.in'
     &,             'formatted sequential rewind')
      read (io,ppmix,end=100)
100   continue
      call relunit (io)

!-----------------------------------------------------------------------
!     set no-flux condition on density difference across bottom level
!-----------------------------------------------------------------------

      do j=1,jmw
        do i=1,imt
          rhom1z(i,km,j) = c0
        enddo
      enddo

!-----------------------------------------------------------------------
!     check for problems
!-----------------------------------------------------------------------

# if defined O_ppvmix && !defined O_implicitvmix
#  if defined O_isopycmix
        write (stdout,'(/,(1x,a))')
     & '==> Error:  "ppvmix" must use "implicitvmix" when "isopycmix" '
     &,'            is also enabled. Also "aidif" should = 0.5        '
        error = .true.
#  endif
# endif
      if (vmixset) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: "ppvmix" cannot be enabled because another   '
     &,'            vertical mixing scheme has been enabled           '
        error = .true.
      else
        vmixset = .true.
      endif

      do k=1,km
        dzmin = min(dzmin,dzt(k))
      enddo
      if (dzmin .ge. 25.e2) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well with coarse vertical  '
     &,'              resolution                                      '
      endif

      extlat = c0
      do jrow=1,jmt
        extlat = max(abs(yt(jrow)),extlat)
      enddo
      if (extlat .gt. 10.0) then
        write (stdout,'(/,(1x,a))')
     & '==> Warning: "ppvmix" may not work well outside the tropics   '
     &,'              where vertical shear is small unless solar      '
     &,'              shortwave penetration into the ocean is         '
     &,'              accounted for by enabeling  "shortwave"         '
      endif

# if !defined O_implicitvmix
      do k=1,km
        if ((dtts*dtxcel(k)*fricmx)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'           "fricmx".  use a smaller "dtts", "dtxcel", and/or  '
     &,'           "fricmx" .... or enable "implicitvmix"             '
          write (stdout,'(a48,i3)') ' at level =',k
          error = .true.
        endif
        if ((dtts*dtxcel(k)*diff_cbt_limit)/dzt(k)**2 .ge. p25) then
          write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "diff_cbt_limit". use a smaller "dtts", "dtxcel"  '
     &,'         ,and/or  "diff_cbt_limit" ...or enable "implicitvmix"'
          write (stdout,'(a48,i3)') ' at level =',k
          error = .true.
        endif
      enddo

      if ((dtuv*fricmx)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "fricmx". use a smaller "dtuv" and/or "fricmx"    '
     &,'             or enable "implicitvmix"                         '
        error = .true.
      endif

      if ((dtuv*visc_cbu_limit)/dzmin**2 .ge. p25) then
        write (stdout,'(/,(1x,a))')
     & '==> Error: vertical diffusive criteria exceeded for           '
     &,'            "visc_cbu_limit". use a smaller "dtuv" or         '
     &,'            "visc_cbu_limit" or enable "implicitvmix"         '
        error = .true.
      endif
# else
        write (stdout,'(/,(1x,a))')
     & '==> Warning: enabeling "implicitvmix" with "ppvmix" uses      '
     &,'             variables defined at "tau" rather than at "tau-1"'
     &,'             as was done in MOM 1.x'
# endif
# if defined O_bryan_lewis_vertical
        write (stdout,'(/,(1x,a/1x,a/1x,a/1x,a))')
     & '==> Warning: "bryan_lewis_vertical" tracer diffus coefficients'
     &,'              will  be added to "ppvmix" diffus coefficients  '
     &,'              Note that diff_cbt_back is being reset to zero  '
     &,'              while diff_cbu_back is unchanged                '
      diff_cbt_back =  0.0
# endif

!     write out namelist values

      write (stdout,ppmix)
      call relunit (io)

      if (error) stop '=> pmixi'

      return
      end

      subroutine ppmix (joff, js, je, is, ie)

!=======================================================================
!     Compute vertical mixing coefficients based on...
!     Pacanowski & Philander (JPO vol 11, #11, 1981).

!     Note: this parameterization was designed for equatorial models
!     and may not do a good job in mid or high latitudes. Simulations
!     in these regions (where vertical shear is small) are improved with
!     the addition of solar short wave penetration into the ocean which
!     reduces buoyancy and enhances vertical mixing.

!     inputs:

!      joff   = offset between rows in the MW and latitude rows
!      js     = starting row for loading variables to calculate
!               coefficients. calculations start at jstrt=max(js-1,jsmw)
!      je     = ending row for loading variables to calculate
!               coefficients. calculations end at je-1
!      is     = starting index for calculating coefficients in the
!               longitude direction
!      ie     = ending index for calculating coefficients in the
!               longitude direction
!      km     = number of vertical levels
!      grav   = gravity (cm/sec**2)
!      umask  = land/sea mask on "u" grid (land=0.0, sea=1.0)
!      tmask  = land/sea mask on "t" grid (land=0.0, sea=1.0)
!      fricmx = max viscosity (cm**2/sec)
!      wndmix = min viscosity at bottom of 1st level to simulate
!               missing high frequency windstress components (cm**2/sec)
!      visc_cbu_back = background "visc_cbu" (cm**2/sec)
!      diff_cbt_back = background "diff_cbt" (cm**2/sec)
!      visc_cbu_limit = largest "visc_cbu" in regions of gravitational
!                      instability (cm**2/sec)
!      diff_cbt_limit = largest "diff_cbt" in regions of gravitational
!                       instability (cm**2/sec)

!     outputs:

!      riu     = richardson number at bottom of "u" cells
!      rit     = richardson number at bottom of "t" cells
!      visc_cbu = viscosity coefficient at bottom of "u" cells (cm**2/s)
!      diff_cbt = diffusion coefficient at bottom of "t" cells (cm**2/s)
!=======================================================================

      implicit none

      integer tlev, istrt, is, iend, ie, joff, js, k, i, ks, je, j
      integer jstrt

      real c0, c1, c5, p25, epsln, fx, t1, rit, t2

      include "size.h"
      include "param.h"
      include "coord.h"
      include "grdvar.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

      real ro(imt,km,1:jmw)

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      c0    = 0.0
      c1    = 1.0
      c5    = 5.0
      p25   = 0.25
      epsln = 1.e-25
      fx    = -4.0*grav
      istrt = max(2,is)
      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     set time level
!-----------------------------------------------------------------------

# if defined O_implicitvmix
      tlev = tau
# else
      tlev = taum1
# endif

!-----------------------------------------------------------------------
!     set vertical density difference for jrow=1 to zero
!-----------------------------------------------------------------------

      if (joff + js .eq. 1) then
        do k=1,km
          do i=istrt,iend
            rhom1z(i,k,1) = c0
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     compute density difference across bottom of "t" cells at tau-1
!-----------------------------------------------------------------------

      do ks=1,2
        call statec (t(1,1,1,1,tlev), t(1,1,1,2,tlev), ro(1,1,jsmw)
     &,              max(js,jsmw), je, istrt, iend, ks)
        do j=max(js,jsmw),je
          do k=ks,km-1,2
            do i=istrt,iend
              rhom1z(i,k,j) = (ro(i,k,j) - ro(i,k+1,j))*tmask(i,k+1,j)
            enddo
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute vertical difference of velocity squared
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km-1
          do i=istrt-1,iend+1
            uzsq(i,k,j) = (u(i,k,j,1,tlev) - u(i,k+1,j,1,tlev))**2 +
     &                    (u(i,k,j,2,tlev) - u(i,k+1,j,2,tlev))**2
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     diffusion and viscosity coefficients on bottom of "T" cells
!-----------------------------------------------------------------------

      jstrt = max(js,jsmw)
      do j=jstrt,je
        do k=1,km-1
          t1 = fx*dzw(k)
          do i=istrt,iend
            rit = t1*rhom1z(i,k,j)/(uzsq(i,k,j) + uzsq(i-1,k,j)
     &                   + uzsq(i,k,j-1) + uzsq(i-1,k,j-1)+ epsln)
            t2  = c1/(c1 + c5*rit)
            diff_cbt(i,k,j) = (fricmx*t2**3 + diff_cbt_back)
     &                        *tmask(i,k+1,j)
            visc_cbt(i,k,j) = (fricmx*t2**2 + visc_cbu_back)
     &                        *tmask(i,k+1,j)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     limit coeffs on bottom of "T" cells
!-----------------------------------------------------------------------

      do j=jstrt,je
        do k=1,km-1
          do i=istrt,iend
            if (rhom1z(i,k,j) .gt. c0) then
              diff_cbt(i,k,j) = diff_cbt_limit
              visc_cbt(i,k,j) = visc_cbu_limit
            endif
          enddo
        enddo
        call setbcx (visc_cbt(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     compute vertical viscosity coeff on "U" cell bottoms
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        do k=1,km-1
          do i=istrt,iend
            visc_cbu(i,k,j) = p25*(visc_cbt(i,k,j) + visc_cbt(i+1,k,j)
     &                        + visc_cbt(i,k,j+1) + visc_cbt(i+1,k,j+1))
     &                        *umask(i,k+1,j)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     approximation for high freq wind mixing near the surface
!     set no flux through bottom of bottom level "km"
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        do i=istrt,iend
          if (diff_cbt(i,1,j) .lt. wndmix) diff_cbt(i,1,j) = wndmix
          if (visc_cbu(i,1,j) .lt. wndmix) visc_cbu(i,1,j) = wndmix
          diff_cbt(i,km,j) = c0
          visc_cbu(i,km,j) = c0
        enddo
      enddo

#if defined O_bryan_lewis_vertical

!-----------------------------------------------------------------------
!     add Bryan-Lewis mixing if wanted
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        do k=1,km-1
          do i=istrt,iend
            diff_cbt(i,k,j) = diff_cbt(i,k,j) + Ahv(k)
          enddo
        enddo
      enddo
#endif

!-----------------------------------------------------------------------
!       set lateral bc
!-----------------------------------------------------------------------

      do j=jsmw,je-1
        call setbcx (visc_cbu(1,1,j), imt, km)
        call setbcx (diff_cbt(1,1,j), imt, km)
      enddo

# if defined O_matrix_sections
      if (prxzts .and. eots) then
        call diagpp (joff, jstrt, je-1)
      endif
# endif

      return
      end

# if defined O_matrix_sections
      subroutine diagpp (joff, js, je)

      implicit none

      integer j, js, je, jrow, joff, reltim, jlat, jj, indp, is, ie
      integer ks, ke, io

      real fx, scl

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cprnts.h"
      include "iounit.h"
      include "switch.h"
      include "tmngr.h"
      include "vmixc.h"

      do j=js,je
        jrow = j + joff
        reltim = relyr
        do jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt)) then
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2

            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'diff_cbt', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (diff_cbt(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                     'unformatted sequential append ieee')
              write (stdout,*) ' => diff_cbt ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((diff_cbt(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, diff_cbt(1,1,j), imt*km)
              call relunit (io)
            endif

            scl = c1
            if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
              write (stdout,9100) 'visc_cbu', itt, jrow
     &,       yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (visc_cbu(1,1,j), imt, is, ie, ks, ke, scl)
            endif
            if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then
              call getunit (io, 'sections.dta'
     &,                     'unformatted sequential append ieee')
              write (stdout,*) ' => visc_cbu ', ' slice: lat='
     &,       yt(jrow), ' written unformatted to file sections.dta'
     &,       ' on ts=', itt, stamp
              iotext = ' read (ioprxz) imt, km, reltim'
              write (io) stamp, iotext, expnam
              write (io) imt, km, reltim
              write(iotext,'(a10,i4)') ' for jrow=',jrow
              iotext(15:)=
     &        ':read(ioprxz)((visc_cbu(i,k),i=1,imt),k=1,km)'
              write (io) stamp, iotext, expnam
              call wrufio (io, visc_cbu(1,1,j), imt*km)
              call relunit (io)
            endif
          endif
        enddo
      enddo

9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
# endif
#endif

      return
      end
