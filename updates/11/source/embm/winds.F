      subroutine add_awind (is, ie, js, je)

#if defined O_embm && defined O_embm_awind
!=======================================================================
!     update winds
!=======================================================================

      implicit none

      integer i, ie, is, j, je, js

      real angle, ax, ay, contr, cosa, drag, f, s, sina, x, y

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "scalar.h"
      include "atm.h"
      include "cembm.h"
      include "csbc.h"

!     calculate advecting wind anomaly

      contr = 0.4           ! (0.8x0.5) contraction x reduction factor
      angle = 10.0/radian   ! half turning angle (from Gill 1982 p.328)
      cosa = cos(angle)
      sina = sin(angle)
      do j=js+1,je
        do i=is+1,ie
          sina = sign(sina, ulat(i,j))
!         damping term for the poles
          f = 1. - (1. - cos(ulat(i,j)/radian))**5.
          x = (awx(i,j)*cosa - awy(i,j)*sina)*f
          y = (awx(i,j)*sina + awy(i,j)*cosa)*f
          sbc(i,j,iwxq) = sbc(i,j,iwxq) + contr*x
          sbc(i,j,iwyq) = sbc(i,j,iwyq) + contr*y
        enddo
      enddo
      call embmbc (sbc(1,1,iwxq))
      call embmbc (sbc(1,1,iwyq))

      drag = cdatm*rhoatm
      contr = 0.8           ! contraction
      angle = 20.0/radian   ! turning angle (from Gill 1982 p.328)
      cosa = cos(angle)
      sina = sin(angle)
      do j=js+1,je
        do i=is+1,ie
          sina = sign(sina, ulat(i,j))
          x = awx(i,j)*cosa - awy(i,j)*sina
          y = awx(i,j)*sina + awy(i,j)*cosa

!         add surface anomaly to wind stress
          f = c1/drag/(sqrt(sqrt(sbc(i,j,itaux)**2
     &      + sbc(i,j,itauy)**2)/drag) + epsln)
          x = contr*x + f*sbc(i,j,itaux)
          y = contr*y + f*sbc(i,j,itauy)
          s = sqrt(x**2 + y**2)
          sbc(i,j,itaux) = drag*x*s
          sbc(i,j,itauy) = drag*y*s

!         add surface anomaly to wind speed
          ax = p25*(awx(i,j) + awx(i-1,j) + awx(i,j-1) + awx(i-1,j-1))
          ay = p25*(awy(i,j) + awy(i-1,j) + awy(i,j-1) + awy(i-1,j-1))
          x = contr*(ax*cosa - ay*sina) + cos(sbc(i,j,iwa))*sbc(i,j,iws)
          y = contr*(ax*sina + ay*cosa) + sin(sbc(i,j,iwa))*sbc(i,j,iws)
          sbc(i,j,iws) = sqrt(x**2 + y**2)
        enddo
      enddo
      call embmbc (sbc(1,1,itaux))
      call embmbc (sbc(1,1,itauy))
      call embmbc (sbc(1,1,iws))

      return
      end

      subroutine calc_awind (is, ie, js, je)

!=======================================================================
!     calculate anomalous pressure and geostrophic wind
!=======================================================================

      implicit none

      integer i, ie, is, j, je, js, n

      real b, rd, s, tclm, tmdl, tmp, adpdx, adpdy, const, diag0, diag1
      real dlat, rlat, slat, rnot, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "scalar.h"
      include "cembm.h"
      include "atm.h"

      real dmsk(imt,jmt)

!-----------------------------------------------------------------------
!     calculate sea level pressure (slp) from temperature. use the
!     equation of state for an ideal gas, P=rho*R*T, and assume a linear
!     relationship between density and temperature rho = s*T + b, where
!     s in g/cm3/K and b in g/cm3 are derived from NCEP and ECMWF data.
!     for the pressure anomaly use dP = R*(s*d(T**2) + b*dT)
!-----------------------------------------------------------------------

      rd = 287.0e4       ! ideal gas constant in cm^2/K/s^2
      s = -4.67e-6       ! slope of rho, T_c relationship in g/cm3/K
      b = 2.58e-3        ! intercept of rho, T_c relationship in g/cm3
      rnot = c1/3600.0   ! time scale for equatorial damping
      dlat = 22.5        ! latitude for equatorial damping
      slat = 30.         ! latitude for slope roll off
      C2K = 273.15       ! convert C to K

      const = 180./(90.-slat)/radian
      do j=js,je
        do i=is,ie
          s = -4.67e-6
!         roll off the slope at high latitudes (more in the south)
          if (tlat(i,j) .lt. -slat) then
            s = s + 1.8e-6*(cos((tlat(i,j) + slat)*const)*0.5 - 0.5)
          elseif (tlat(i,j) .gt. slat) then
            s = s + 0.9e-6*(cos((tlat(i,j) - slat)*const)*0.5 - 0.5)
          endif
          tmdl = rtbar(i,j) + C2K
          tclm = tbar(i,j) + C2K
          apress(i,j) = rd*(s*(tmdl**2 - tclm**2) + b*(tmdl - tclm))
        enddo
      enddo
!     normalise pressure anomaly for output
      dmsk(:,:) = 1.
      call areaavg (apress, dmsk, tmp)
      apress(:,:) = apress(:,:) - tmp

      call embmbc (apress)

      do j=js,je-1
        do i=is,ie-1
          diag1 = apress(i+1,j+1) - apress(i,j)
          diag0 = apress(i,j+1) - apress(i+1,j)
          adpdy  = (diag1 + diag0)*dyu2r(j)
          adpdx  = (diag1 - diag0)*dxu2r(i)*cstr(j)
          rlat = rnot*exp(-abs(ulat(i,j)/dlat))
          const = c1/(rhoatm*(rlat**2 + fcor(i,j)**2))
          awy(i,j) =  const*(fcor(i,j)*adpdx - rlat*adpdy)
          awx(i,j) = -const*(rlat*adpdx + fcor(i,j)*adpdy)
        enddo
      enddo

      call embmbc (awx)
      call embmbc (awy)
#endif

      return
      end
