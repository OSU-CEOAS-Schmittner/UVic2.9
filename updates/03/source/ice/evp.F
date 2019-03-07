      subroutine evp (is, ie, js, je)

#if defined O_ice && defined O_ice_evp && defined O_embm
!=======================================================================
!     calculate velocities using an "elastic-viscous-plastic" rheology

!     see E. C. Hunke and J. K. Dukowicz. An elastic-viscous-plastic
!       model for sea ice dynamics. J. Phys. Oceanogr., 1997.
!=======================================================================

      implicit none

      integer i, j, k, ie, je, is, js

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"

      do j=1,jmt
        do i=1,imt
          xint(i,j) = 0.0
          yint(i,j) = 0.0
        enddo
      enddo

      call mass_prss
      call viscevp
      call stressprep

      do k=1,ndte
        call stressevp
        call stepu
        call embmbc (uice)
        call embmbc (vice)
      enddo

# if defined O_ice_cpts
!  compute the strain rate tensor
      call strain
# endif

      return
      end

      subroutine viscevp

! $Id: dyn.F,v 1.3 1997/02/11 18:04:27 eclare Exp $

!.. Elastic-viscous-plastic sea ice dynamics model
!.. Computes ice velocity

!.. author Elizabeth C. Hunke
!..        Fluid Dynamics Group, Los Alamos National Laboratory

!.. See E. C. Hunke and J. K. Dukowicz. An elastic-viscous-plastic model
!..     for sea ice dynamics. J. Phys. Oceanogr., 1997.

!.. Copyright, 1997.  The Regents of the University of California.
!.. This software was produced under a U.S. Government contract
!.. (W-7405-ENG-36) by Los Alamos National Laboratory, which is operated
!.. by the University of California for the U.S. Department of Energy.
!.. The U.S. Government is licensed to use, reproduce, and distribute this
!.. software.  Permission is granted to the public to copy and use this
!.. software without charge, provided that this Notice and any statement
!.. of authorship are reproduced on all copies.  Neither the Government
!.. nor the University makes any warranty, express or implied, or assumes
!.. any liability or responsibility for the use of this software.

!=======================================================================
!.. Computes the rates of strain, and the bulk and shear viscosities
!.. zeta and eta.  The rates of strain (xi*) and viscosities are
!.. calculated for each of the four triangles in each grid cell
!.. (north, south, east, west).
!=======================================================================

      implicit none

      integer n, j, i

      real prs, pf, cc, dd, xi11n, xi12n, xi22n, xi11e, xi12e, xi22e
      real xi11s, xi12s, xi22s, xi11w, xi12w, xi22w, deltan, deltae
      real deltas, deltaw, zetamax

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"

      do n=1,nseg
        do j=jsi(n),jei(n)+1
          do i=2,imt
            prs = 0.5*pice(i,j)     ! g/s^2  (P/2, varies with c*H)
            pf = prs/floor            ! initializes zeta to large values

!           initialize zeta                   ! g/s
            zetan(i,j) = pf
            zetae(i,j) = pf
            zetas(i,j) = pf
            zetaw(i,j) = pf

!           rates of strain                    ! 1/s
            cc = (uice(i,j) + uice(i-1,j) - uice(i,j-1)
     &         - uice(i-1,j-1))*dyt2r(j)
            dd = (vice(i,j) + vice(i,j-1) - vice(i-1,j)
     &         - vice(i-1,j-1))*cstr(j)*dxt2r(i)

            xi11n = (uice(i,j) - uice(i-1,j))*csur(j)*dxur(i)
            xi12n = ((vice(i,j) - vice(i-1,j))*csur(j)*dxur(i)
     &            + cc)*0.5
            xi22n = (vice(i,j) + vice(i-1,j) - vice(i,j-1)
     &                  - vice(i-1,j-1))*dyt2r(j)

            xi11e = (uice(i,j) + uice(i,j-1) - uice(i-1,j)
     &            - uice(i-1,j-1))*cstr(j)*dxt2r(i)
            xi12e = ((uice(i,j) - uice(i,j-1))*dyur(j) + dd)*0.5
            xi22e = (vice(i,j) - vice(i,j-1))*dyur(j)

            xi11s = (uice(i,j-1) - uice(i-1,j-1))*csur(j-1)*dxur(i)
            xi12s = ((vice(i,j-1) - vice(i-1,j-1))*
     &              csur(j-1)*dxur(i) + cc)*0.5
            xi22s = xi22n

            xi11w = xi11e
            xi12w = ((uice(i-1,j) - uice(i-1,j-1))*dyur(j)
     &            + dd)*0.5
            xi22w = (vice(i-1,j) - vice(i-1,j-1))*dyur(j)

!           Delta (in the denominator of zeta, eta)        ! 1/s
            deltan = sqrt( (xi11n**2 + xi22n**2)*ecc2p
     &        + 4.0*xi12n**2*ecc2 + xi11n*xi22n*ecc2m)
            deltae = sqrt( (xi11e**2 + xi22e**2)*ecc2p
     &        + 4.0*xi12e**2*ecc2 + xi11e*xi22e*ecc2m)
            deltas = sqrt( (xi11s**2 + xi22s**2)*ecc2p
     &        + 4.0*xi12s**2*ecc2 + xi11s*xi22s*ecc2m)
            deltaw = sqrt( (xi11w**2 + xi22w**2)*ecc2p
     &        + 4.0*xi12w**2*ecc2 + xi11w*xi22w*ecc2m)

            deltan = amax1(1.e-20,deltan)
            deltae = amax1(1.e-20,deltae)
            deltas = amax1(1.e-20,deltas)
            deltaw = amax1(1.e-20,deltaw)

!           bulk viscosity zeta, bounded by maximum, minimum values
            zetamax = 2.5e8*pice(i,j)  !Hibler pg 819    ! g/s

            zetan(i,j) = prs/deltan
            zetan(i,j) = min(zetamax,zetan(i,j))
            zetan(i,j) = max(zetamin,zetan(i,j))

            zetae(i,j) = prs/deltae
            zetae(i,j) = min(zetamax,zetae(i,j))
            zetae(i,j) = max(zetamin,zetae(i,j))

            zetas(i,j) = prs/deltas
            zetas(i,j) = min(zetamax,zetas(i,j))
            zetas(i,j) = max(zetamin,zetas(i,j))

            zetaw(i,j) = prs/deltaw
            zetaw(i,j) = min(zetamax,zetaw(i,j))
            zetaw(i,j) = max(zetamin,zetaw(i,j))

!           mask zeta (and therefore eta)
            zetan(i,j) = zetan(i,j)*tmsk(i,j)
            zetae(i,j) = zetae(i,j)*tmsk(i,j)
            zetas(i,j) = zetas(i,j)*tmsk(i,j)
            zetaw(i,j) = zetaw(i,j)*tmsk(i,j)

!           shear viscosity eta                        ! g/s
            etan(i,j) = zetan(i,j)*ecc2
            etae(i,j) = zetae(i,j)*ecc2
            etas(i,j) = zetas(i,j)*ecc2
            etaw(i,j) = zetaw(i,j)*ecc2

          enddo
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine stressprep
!=======================================================================
!.. Computes quantities used in the subroutines stress and stepu.
!.. This subroutine is cryptic - I apologize - but it made the code
!.. faster by almost a factor of 2 on a Cray-YMP.  Here we compute
!.. quantities needed in the stress tensor (sigma) and momentum (u)
!.. equations, but which don't change during the subcycling. Many of
!.. these variables are grouped in parties of four, one for each
!.. triangle (north, east, south, west).
!=======================================================================

      implicit none

      integer n, j, i, nc

      real costh, sinth, econst, ey, e2, en, ee, es, ew, zn, ze, zs
      real zw, c2n, c2e, c2s, c2w, c3n, c3e, c3s, c3w, d2n, d2e
      real d2s, d2w, a2n, a2e, a2s, a2w, prs

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "levind.h"
      include "grdvar.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"
      include "csbc.h"

      costh = 0.9063
      sinth = 0.4226
!     used to compute E = Econst*c*H
      Econst = 2.*eyc*rhoice*xyminevp*dtei**2  ! g/cm s^2

      do n=1,nseg
        do j=jsi(n),jei(n)+1
          do i=2,imt

            if (tmsk(i,j) .gt. floor) then

# if defined O_ice_cpts
              ey = heff(i,j,idx,1)
              do nc=2,ncat
                 ey = ey + heff(i,j,idx,nc)
              enddo
              ey = Econst*ey                   ! E, g/s^2
# else
              ey = Econst*hice(i,j,2)          ! E, g/s^2
# endif
              ey = max(ey,floor)

              e2   = 0.5*ey                    ! E/2
              edy(i,j)  = e2*dytr(j)           ! E/(2*dy)
              edx(i,j)  = e2*cstr(j)*dxtr(i)   ! E/(2*dx)
              eHN(i,j)  = e2/(csu(j)*dxu(i))   ! E/(2*HTN)
              eHE(i,j)  = e2/dyu(j)            !   etc
              eHNm(i,j) = e2/(csu(j-1)*dxu(i))
              eHEm(i,j) = e2/dyu(j)

              en = e2/etan(i,j)         !     E
              ee = e2/etae(i,j)         !  -------
              es = e2/etas(i,j)         !   2*eta
              ew = e2/etaw(i,j)

              zn = e2/zetan(i,j)        !     E
              ze = e2/zetae(i,j)        !  --------
              zs = e2/zetas(i,j)        !   2*zeta
              zw = e2/zetaw(i,j)

              c2n = dtei + en           !   1      E
              c2e = dtei + ee           !  --- + -----
              c2s = dtei + es           !  dte   2*eta
              c2w = dtei + ew

              c3n = 0.5*(en - zn)       !  E      1     1
              c3e = 0.5*(ee - ze)       !  - * ( --- - ---- )
              c3s = 0.5*(es - zs)       !  4     eta   zeta
              c3w = 0.5*(ew - zw)

              d2n = c2n - c3n           !  1      E     E     1     1
              d2e = c2e - c3e           ! --- + ----- - - * (--- - ----)
              d2s = c2s - c3s           ! dte   2*eta   4    eta   zeta
              d2w = c2w - c3w

              h2n(i,j) = 1./c2n         ! this rapidly gets out of hand
              h2e(i,j) = 1./c2e
              h2s(i,j) = 1./c2s
              h2w(i,j) = 1./c2w

              a2n = h2n(i,j)/(d2n - c3n)
              a2e = h2e(i,j)/(d2e - c3e)
              a2s = h2s(i,j)/(d2s - c3s)
              a2w = h2w(i,j)/(d2w - c3w)

              b2n(i,j) = a2n*d2n
              b2e(i,j) = a2e*d2e
              b2s(i,j) = a2s*d2s
              b2w(i,j) = a2w*d2w

              a2na(i,j) = a2n*c3n
              a2ea(i,j) = a2e*c3e
              a2sa(i,j) = a2s*c3s
              a2wa(i,j) = a2w*c3w

              prs  = 0.5*pice(i,j)
              prssn(i,j) = prs*zn       !    P*E
              prsse(i,j) = prs*ze       !  --------
              prsss(i,j) = prs*zs       !   4*zeta
              prssw(i,j) = prs*zw

            endif

          enddo
        enddo
      enddo

!     for subroutine stepu
      do n=1,nseg
        do j=jsi(n),jei(n)+1
          do i=2,imt
            HTN4(i,j) = 0.25/(csu(j)*dxu(i))
            HTE4(i,j) = 0.25/dyu(j)
            dxt8(i,j) = 0.125/(cst(j)*dxt(i))
            dyt8(i,j) = 0.125/dyt(j)
          enddo
        enddo
      enddo

      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=2,imtm1
            fmass(i,j) = fcor(i,j)*umass(i,j)      ! Coriolis * mass
            sinth = sign(sinth,fmass(i,j))
# if defined O_correct_ice_to_ocean_stress
!           for water stress
            waterx(i,j) = umsk(i,j)*(sbc(i,j,isu)*costh -
     &                sbc(i,j,igv)*sinth)
            watery(i,j) = umsk(i,j)*(sbc(i,j,isv)*costh +
     &                sbc(i,j,igu)*sinth)
# else
!           for water stress
            waterx(i,j) = umsk(i,j)*(sbc(i,j,igu)*costh -
     &                sbc(i,j,igv)*sinth)
            watery(i,j) = umsk(i,j)*(sbc(i,j,igv)*costh +
     &                sbc(i,j,igu)*sinth)
# endif
!           air stress and -mg*gradH_o term (tilt)
            strairx(i,j) = umsk(i,j)*(sbc(i,j,itaux) -
     &                fmass(i,j)*sbc(i,j,igv))
            strairy(i,j) = umsk(i,j)*(sbc(i,j,itauy) +
     &                fmass(i,j)*sbc(i,j,igu))
          enddo
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine stressevp
!=======================================================================
!.. Calculates the internal stress components, sigma_ij, in the four
!.. triangles of each cell.
!=======================================================================

      implicit none

      integer n, j, i

      real dun, dus, due, duw, dvn, dvs, dve, dvw, cc, dd, xi11n, xi12n
      real xi22n, xi11e, xi12e, xi22e, xi11s, xi12s, xi22s, xi11w
      real xi12w, xi22w, c4n, c4e, c4s, c4w, c5n, c5e, c5s, c5w

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"

      do n=1,nseg
        do j=jsi(n), jei(n) + 1
          do i=2,imt

            if (tmsk(i,j) .gt. floor) then

              dun = uice(i,j) - uice(i-1,j)
              dus = uice(i,j-1) - uice(i-1,j-1)
              due = uice(i,j) - uice(i,j-1)
              duw = uice(i-1,j) - uice(i-1,j-1)

              dvn = vice(i,j) - vice(i-1,j)
              dvs = vice(i,j-1) - vice(i-1,j-1)
              dve = vice(i,j) - vice(i,j-1)
              dvw = vice(i-1,j) - vice(i-1,j-1)

              cc = 0.5*edy(i,j)*(due + duw)
              dd = 0.5*edx(i,j)*(dvn + dvs)

!             NOTE these are rates of strain * E
              xi11n = 2.0*dun*eHN(i,j)
              xi12n = dvn*eHN(i,j) + cc
              xi22n = edy(i,j)*(dve + dvw)

              xi11e = edx(i,j)*(dun + dus)
              xi12e = due*eHE(i,j) + dd
              xi22e = 2.0*dve*eHE(i,j)

              xi11s = 2.0*dus*eHNm(i,j)
              xi12s = dvs*eHNm(i,j) + cc
              xi22s = xi22n

              xi11w = xi11e
              xi12w = duw*eHEm(i,j) + dd
              xi22w = 2.0*dvw*eHEm(i,j)

!             solve for the three components of sigma in each triangle
              c4n = dtei*sig11n(i,j) + xi11n - prssn(i,j)
              c4e = dtei*sig11e(i,j) + xi11e - prsse(i,j)
              c4s = dtei*sig11s(i,j) + xi11s - prsss(i,j)
              c4w = dtei*sig11w(i,j) + xi11w - prssw(i,j)

              c5n = dtei*sig22n(i,j) + xi22n - prssn(i,j)
              c5e = dtei*sig22e(i,j) + xi22e - prsse(i,j)
              c5s = dtei*sig22s(i,j) + xi22s - prsss(i,j)
              c5w = dtei*sig22w(i,j) + xi22w - prssw(i,j)

              sig11n(i,j) = a2na(i,j)*c5n + c4n*b2n(i,j)        ! g/s^2
              sig11e(i,j) = a2ea(i,j)*c5e + c4e*b2e(i,j)
              sig11s(i,j) = a2sa(i,j)*c5s + c4s*b2s(i,j)
              sig11w(i,j) = a2wa(i,j)*c5w + c4w*b2w(i,j)

              sig22n(i,j) = a2na(i,j)*c4n + c5n*b2n(i,j)
              sig22e(i,j) = a2ea(i,j)*c4e + c5e*b2e(i,j)
              sig22s(i,j) = a2sa(i,j)*c4s + c5s*b2s(i,j)
              sig22w(i,j) = a2wa(i,j)*c4w + c5w*b2w(i,j)

              sig12n(i,j) = h2n(i,j)*(xi12n + dtei*sig12n(i,j))
              sig12e(i,j) = h2e(i,j)*(xi12e + dtei*sig12e(i,j))
              sig12s(i,j) = h2s(i,j)*(xi12s + dtei*sig12s(i,j))
              sig12w(i,j) = h2w(i,j)*(xi12w + dtei*sig12w(i,j))

            endif

          enddo
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------

      subroutine mass_prss
!=======================================================================
!.. Computes ice mass and pressure (strength)
!=======================================================================

      implicit none

      integer n, j, i, nc

      real pstar, volm, area

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"

      real tmass(imt,jmt)

      pstar = 2.75e5
!     total mass of ice and snow, centred in T-cell
      do n=1,nseg
        do j=jsi(n), jei(n) + 1
          do i=1,imt
            if (tmsk(i,j) .ne. 0.0) then
!             if you, like Hibler (1979), have only thick ice
# if defined O_ice_cpts
              tmass(i,j) = heff(i,j,idx,1)
              do nc=2,ncat
                tmass(i,j) = tmass(i,j) + heff(i,j,idx,nc)
              enddo
              tmass(i,j) = rhoice*tmass(i,j)
# else
              tmass(i,j) = rhoice*hice(i,j,2)    ! g/cm^2
# endif
            else
                tmass(i,j) = 0.0
            endif
          enddo
        enddo
      enddo

!     mass centred at velocity nodes (U-cells)
      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=2,imtm1
            umass(i,j) = 0.25*(tmass(i,j) + tmass(i+1,j)
     &         + tmass(i,j+1) + tmass(i+1,j+1))              ! g/cm^2
          enddo
        enddo
      enddo

!     pressure P
      do n=1,nseg
        do j=jsi(n),jei(n)+1
          do i=2,imt
#  if defined O_ice_cpts
#   if defined O_ice_cpts_roth_press
            pice(i,j) = strength(i,j,idx)
#   else
            volm = heff(i,j,idx,1)
            area = A(i,j,idx,1)
            do nc=2,ncat
               volm = volm + heff(i,j,idx,nc)
               area = area + A(i,j,idx,nc)
            enddo
            pice(i,j) = pstar*volm*exp(-20.0*(1.0 - area))
#   endif
#  else
            pice(i,j) = pstar*hice(i,j,2)
     &                    *exp(-20.0*(1.0 - aice(i,j,2)))      ! g/s^2
#  endif
          enddo
        enddo
      enddo

      call embmbc (pice)

      return
      end
!-----------------------------------------------------------------------

      subroutine stepu
!=======================================================================
!.. Calculation of the surface stresses
!.. Integration of the momentum equation to find velocity (u,v)
!=======================================================================

      implicit none

      integer n, j, i

      real costh, sinth, uorel, vorel, vrel, umassdtei, cca, ccb, ab2
      real s11, s12, s21, s22, c1evp, c2evp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"
      include "csbc.h"

      real s11ns(imt,jmt), s11ew(imt,jmt)
      real s22ew(imt,jmt), s22ns(imt,jmt)
      real s12ewi(imt,jmt), s12ns(imt,jmt)
      real s12nsj(imt,jmt), s12ew(imt,jmt)

      costh = 0.9063
      sinth = 0.4226
!     some cryptic but useful arrays
      do n=1,nseg
        do j=jsi(n), jei(n) + 1
          do i=2,imt
            s11ew(i,j) = dxt8(i,j)*(sig11e(i,j) + sig11w(i,j))
            s22ns(i,j) = dyt8(i,j)*(sig22n(i,j) + sig22s(i,j))
            s12ns(i,j) = dyt8(i,j)*(sig12n(i,j) + sig12s(i,j))
            s12ew(i,j) = dxt8(i,j)*(sig12e(i,j) + sig12w(i,j))
          enddo
        enddo
      enddo

      do n=1,nseg
        do j=jsi(n), jei(n) + 1
          do i=2,imtm1
            s22ew(i,j) = HTE4(i,j)*(sig22e(i,j) + sig22w(i+1,j))
            s12ewi(i,j) = HTE4(i,j)*(sig12e(i,j) + sig12w(i+1,j))
          enddo
        enddo
      enddo

      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=2,imt
            s11ns(i,j) = HTN4(i,j)*(sig11s(i,j+1) + sig11n(i,j))
            s12nsj(i,j) = HTN4(i,j)*(sig12s(i,j+1) + sig12n(i,j))
          enddo
        enddo
      enddo

!     integrate the momentum equation
      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=2,imtm1

            if (umsk(i,j) .gt. floor .and. umass(i,j) .gt. 0.01) then
# if defined O_correct_ice_to_ocean_stress
!             ocean currents relative to ice velocity
!             (isu -> surface layer, isv -> second layer)
              uorel = sbc(i,j,isu) - uice(i,j)
              vorel = sbc(i,j,isv) - vice(i,j)
!             (magnitude of relative geostrophic ocean current)*rhow*drag
              vrel = 0.0055*1.03*sqrt(uorel**2 + vorel**2)  ! cm/s
!              vrel = dragw*sqrt(uorel**2 + vorel**2)  ! cm/s
# else
!             geostrophic ocean currents relative to ice velocity
              uorel = sbc(i,j,igu) - uice(i,j)
              vorel = sbc(i,j,igv) - vice(i,j)
!             (magnitude of relative geostrophic ocean current)*rhow*drag
              vrel = 0.0055*1.03*sqrt(uorel**2 + vorel**2)  ! cm/s
!              vrel = dragw*sqrt(uorel**2 + vorel**2)  ! cm/s
# endif
!             alpha, beta are defined in Hunke and Dukowicz sec 3.2
              umassdtei = umass(i,j)*dtei  ! m/dte,alpha, beta, g/cm^2 s
              cca = umassdtei + vrel*costh
              ccb = fmass(i,j) + vrel*sign(sinth,fmass(i,j))
              ab2 = cca**2 + ccb**2

!             more cryptic stuff
              s11 = - s11ns(i,j) + s11ns(i+1,j) + s11ew(i+1,j+1)
     &           + s11ew(i+1,j) - s11ew(i,j+1) - s11ew(i,j)

              s12 = - s12ewi(i,j) + s12ewi(i,j+1) + s12ns(i+1,j+1)
     &           + s12ns(i,j+1) - s12ns(i+1,j) - s12ns(i,j)

              s21 = - s12nsj(i,j) + s12nsj(i+1,j) + s12ew(i+1,j+1)
     &           + s12ew(i+1,j) - s12ew(i,j+1) - s12ew(i,j)

              s22 = - s22ew(i,j) + s22ew(i,j+1) + s22ns(i+1,j+1)
     &           + s22ns(i,j+1) - s22ns(i+1,j) - s22ns(i,j)
# if defined O_correct_ice_to_ocean_stress
!             finally, the velocity components (in cm/s)
              c1evp = s11 + s12 + strairx(i,j) + vrel*waterx(i,j)
     &              + umassdtei*uice(i,j)
              c2evp = s21 + s22 + strairy(i,j) + vrel*watery(i,j)
     &              + umassdtei*vice(i,j)
              uice(i,j) = (cca*c1evp + ccb*c2evp)/ab2
              vice(i,j) = (cca*c2evp - ccb*c1evp)/ab2
              sinth = sign(sinth,fmass(i,j))
              xint(i,j) = - (vrel*waterx(i,j)-vrel*umsk(i,j)
     &           *(uice(i,j)*costh - vice(i,j)*sinth))
              yint(i,j) = - (vrel*watery(i,j)-vrel*umsk(i,j)
     &           *(vice(i,j)*costh + uice(i,j)*sinth))
# else
              xint(i,j) = s11 + s12
              yint(i,j) = s21 + s22

!             finally, the velocity components (in cm/s)
              c1evp = xint(i,j) + strairx(i,j) + vrel*waterx(i,j)
     &              + umassdtei*uice(i,j)
              c2evp = yint(i,j) + strairy(i,j) + vrel*watery(i,j)
     &              + umassdtei*vice(i,j)
              uice(i,j) = (cca*c1evp + ccb*c2evp)/ab2
              vice(i,j) = (cca*c2evp - ccb*c1evp)/ab2
# endif
            else
!             set velocity to zero on land and (nearly) open water
              uice(i,j) = 0.
              vice(i,j) = 0.
            endif

          enddo
        enddo
      enddo

      return
      end

      subroutine strain

!     B-grid strain rate tensor for EVP model (lacks metric
!     terms, see subroutine strain in vpadi.F for full calc)

      implicit none

      integer n, j, i

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "evp.h"
      include "grdvar.h"

      real e11(imt,jmt),e22(imt,jmt),e12(imt,jmt)

      do n=1,nseg
        do j=jsi(n),jei(n)
          do i=2,imtm1
            e11(i,j) = 0.5*cstr(j)*dxtr(i)*( uice(i,j) + uice(i,j-1)
     &               - uice(i-1,j) - uice(i-1,j-1) )
            e22(i,j) = 0.5*cstr(j)*dytr(j)*(
     &                 csu(j)*( vice(i,j) + vice(i-1,j) )
     &               - csu(j-1)*( vice(i,j-1) + vice(i-1,j-1) ) )
            e12(i,j) = 0.25*cstr(j)*( dytr(j)*(
     &                 csu(j)*( uice(i,j) + uice(i-1,j) )
     &               - csu(j-1)*( uice(i,j-1) - uice(i-1,j-1) ) )
     &             +  dxtr(i)*( vice(i,j) +  vice(i,j-1)
     &               - vice(i-1,j) - vice(i-1,j-1) ) )
            eI(i,j)=e11(i,j)+e22(i,j)
            del(i,j) = (e11(i,j)**2 + e22(i,j)**2)*(1.0 + ecc2)
     &        + 4.0*ecc2*e12(i,j)**2 + 2.0*e11(i,j)*e22(i,j)*
     &        (1.0 - ecc2)
            del(i,j) = sqrt(del(i,j))
          enddo
        enddo
      enddo

      call embmbc (del)
      call embmbc (eI)

#endif
      return
      end
