      subroutine state (t, s, rho, js, je, is, ie)

#if defined O_mom
!=======================================================================

!     state computes normalized densities by using a 3rd
!     order polynomial fit to the knudsen formula, for each level
!     subtract normalizing constants from temperature and salinity
!     and compute polynomial approximation of knudsen density.

!     note.. for precision purposes, there is a depth dependent
!     constant subtracted from the density returned by this routine.
!     so... this routine should be used only for horizontal gradients
!     of density.

!     inputs:

!     t  = the input row of temperatures (potential deg C)
!     s  = the input row of salinities (units: (ppt-35)/1000)
!     js = starting row for computing density within the MW
!     je = ending row for computing density within the MW
!     is = starting longitude index for computing density within the MW
!     ie = ending longitude index for computing density within the MW

!     output:

!     rho = normalized densities
!     These densities are in cgs units(g/cm3) and represent
!     the in situ density at a level minus a depth dependent
!     normalization. The complete in situ density is given by
!     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k, c)
!                           + rho_norm(k)*10-3,
!     where rho_norm(k) are the depth dependent normalization densities
!     [in sigma units (density-1)*1000] given at the bottom of dncoef.h

!=======================================================================

      implicit none

      integer k, j, js, je, i, is, ie, ind, l

      real dens, tq, sq, drodt, drods

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "state.h"

      real t(imt,km,jmw), s(imt,km,jmw), rho(imt,km,jsmw:jmw)

      include "dens.h"

      do j=js,je
        do k=1,km
          do i=is,ie
            rho(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k)
          enddo
        enddo
      enddo

      return

      entry statec (t, s, rho, js, je, is, ie, ind)

!=======================================================================

!     statec computes, for one row, the normalized densities by using
!     a 3rd order polynomial fit to the knudsen formula. For
!     purposes of checking vertical stability between adjacent
!     levels, the reference depth for pressure dependence in
!     the knudsen formula must be held constant. that level is
!     determined by "ind".

!     inputs:

!     t   = the input row of temperatures (potential deg C)
!     s   = the input row of salinities (units: (ppt-35)/1000)
!     js  = starting row for computing density within the MW
!     je  = ending row for computing density within the MW
!     is  = starting longitude index for computing density within the MW
!     ie  = ending longitude index for computing density within the MW
!     ind = 1 for comparing levels 1 to 2, 3 to 4, etc.
!           (coefficients for the lower of the 2 levels are used)
!           2 for comparing levels 2 to 3, 4 to 5, etc.
!           (coefficients for the lower of the 2 levels are used)

!     output:

!     rho = normalized densities
!     These densities are in cgs units(g/cm3) and represent
!     the in situ density at a level minus a depth dependent
!     normalization. The complete in situ density is given by
!     rho_complete(i,k,j) = dens (t(i,k,j)-to(k), s(i,k,j)-so(k), k, c)
!                           + rho_norm(k)*10-3,
!     where rho_norm(k) are the depth dependent normalization densities
!     [in sigma units (density-1)*1000] given at the bottom of dncoef.h

!=======================================================================

      if (ind .lt. 1 .or. ind .gt. 2) then
        write (stderr,99) ind
        stop '=>statec'
      endif

      do j=js,je
        do l=1,km,2
          if (ind .eq. 1) then
            k = min(l+1,km)
          else
            k = l
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo

      do j=js,je
        do l=2,km,2
          if (ind .eq. 1) then
            k = l
          else
            k = min(l+1,km)
          endif
          do i=is,ie
            rho(i,l,j) = dens (t(i,l,j)-to(k), s(i,l,j)-so(k), k)
          enddo
        enddo
      enddo

      return
   99 format(/' error => bad "ind" in statec: ind =',i10)

      end

      subroutine state_ref (t, s, rho, js, je, is, ie, kr)

!=======================================================================

!     Construct potential density referenced to the level kr.

!=======================================================================

      implicit none

      integer j, js, je, i, is, ie, kr, k

      real dens, tq, sq, drodt, drods

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "state.h"

      real t(imt,km,jmw), s(imt,km,jmw), rho(imt,km,jsmw:jmw)

      include "dens.h"

      do j=js,je
        do k=1,km
          do i=is,ie
            rho(i,k,j) = dens (t(i,k,j)-to(kr), s(i,k,j)-so(kr), kr)
          enddo
        enddo
      enddo
#endif

      return
      end
