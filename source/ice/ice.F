      subroutine ice (is, ie, js, je)

#if defined O_ice && defined O_embm
!=======================================================================
!     calling routine for ice model subroutines
!=======================================================================

      implicit none

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

# if defined O_ice_evp
      if (nivc .eq. 1) then

!-----------------------------------------------------------------------
!       find latitudes with ice for dynamic calculations
!-----------------------------------------------------------------------

        call icelats (is, ie, js, je)
#  if defined O_ice_evp

!-----------------------------------------------------------------------
!       calculate velocities with elastic-viscous-plastic rheology
!-----------------------------------------------------------------------

        if (nseg .ne. 0) then
          call evp (is, ie, js, je)

#   if defined O_mom
#    if defined O_ice_fourfil || defined O_ice_firfil
          call filuvice (js, je)
#    endif
#   endif
        endif
#  endif

      endif
      if (nivc .eq. nivts) nivc = 0
      nivc = nivc + 1

#  if !defined O_ice_cpts
!-----------------------------------------------------------------------
!     advect ice and snow
!-----------------------------------------------------------------------

      if (lf. eq. 2) then
        do j=jsp1,jem1
          do i=isp1,iem1
            hice(i,j,1) = hice(i,j,2)
            aice(i,j,1) = aice(i,j,2)
            hsno(i,j,1) = hsno(i,j,2)
          enddo
        enddo
      endif
      call advupb (uice, vice, hice(1,1,1), is, ie, js, je)
      call advupb (uice, vice, aice(1,1,1), is, ie, js, je)
      call advupb (uice, vice, hsno(1,1,1), is, ie, js, je)
#  endif
#  if defined O_ice_cpts

!-----------------------------------------------------------------------
!     calculate advection and mechanical redistribution for cpts ice
!-----------------------------------------------------------------------

      call adv_ridge_cpts (is, ie, js, je)
#  endif
# endif
# if defined O_ice_cpts

!-----------------------------------------------------------------------
!     calculate thermodynamics for cpts ice
!-----------------------------------------------------------------------

      call thermo (is, ie, js, je)
# endif
!-----------------------------------------------------------------------
!     calculate simple thermodynamics
!-----------------------------------------------------------------------

      call therm (is, ie, js, je)
#endif

      return
      end

#if defined O_ice && defined O_ice_evp && defined O_embm
      subroutine icelats (is, ie, js, je)

!=======================================================================
!     find latitudes for ice dynamics calculations
!     calculate starting and ending latitudes for velocity calculations
!=======================================================================

      implicit none

      integer i, ice, ie, iem1, inc, is, isp1, j, je, jem1, jm, jp
      integer js, jsp1, k, kmax, nc

      real hi

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      nseg = 0
      inc = 1
      kmax = 3
      do j=jsp1,jem1
        ice = 0
        do i=isp1,iem1
#   if defined O_ice_cpts
          hi = heff(i,j,1,1) + heff(i,j,2,1)
          do nc=2,ncat
             hi = hi + heff(i,j,1,nc) + heff(i,j,2,nc)
          enddo
#   else
          hi = hice(i,j,1) + hice(i,j,2)
#   endif
          do k=1,kmax
            jm = max(1,j-k)
            jp = min(jmt,j+k)
#   if defined O_ice_cpts
            hi = hi + heff(i,jm,1,1) + heff(i,jp,1,1) + heff(i,jm,2,1)
     &         + heff(i,jp,2,1)
            do nc=2,ncat
               hi = hi + heff(i,jm,1,nc) + heff(i,jp,1,nc)
     &              + heff(i,jm,2,nc) + heff(i,jp,2,nc)
            enddo
#   else
            hi = hi + hice(i,jm,1) + hice(i,jp,1) + hice(i,jm,2)
     &         + hice(i,jp,2)
#   endif
            if (hi .ne. 0.0) ice = 1
          enddo
        enddo
        if (ice .eq. 1) then
          nseg = nseg + inc
          if (inc .eq. 1) jsi(nseg) = j
          if (j .eq. jmtm1) jei(nseg) = jmtm1
          inc = 0
        else
          do i=is,ie
            uice(i,j) = 0.0
            vice(i,j) = 0.0
          enddo
          if (inc .eq. 0) jei(nseg) = j
          inc = 1
        endif
      enddo

      return
      end

      subroutine filuvice (js, je)
# if defined O_fourfil || defined O_firfil

!=====================================================================
!     filuvice sets up input needed for Fourier filtering
!     (when the "fourfil" option is defined) -or- symmetric finite
!     impulse response filtering (when ifdef "firfil" is defined) of
!     baroclinic velocities at the specified high latitude row "j".
!=====================================================================

      implicit none

      integer j, js, je, jj, isave, ieave, l, is, ie, iredo, im, m
      integer n, ism1, iea, i, ieb, ii, icyc

      real fx, avga, avgb

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cpolar.h"
      include "emode.h"
      include "grdvar.h"
      include "index.h"
      include "atm.h"
# if defined O_ice_cpts
      include "cpts.h"
# endif
      include "ice.h"
      include "scalar.h"
      include "switch.h"

      real tempi(imt,2)

      call embmbc (uice)
      call embmbc (vice)

!---------------------------------------------------------------------
!     fourier filter u and v at high latitudes
!---------------------------------------------------------------------

      do j=js,je
        if ((j .gt. jfu1 .and. j .lt. jfu2) .or. j .lt. jfrst) goto 701
        jj = j - jfrst + 1
        if (j .ge. jfu2) jj = jj - jskpu + 1
        fx = -c1
        if (phi(j) .gt. c0) fx = c1
        isave = 0
        ieave = 0

        do l=1,lsegf
          if (isuf(jj,l,1) .ne. 0) then
            is = isuf(jj,l,1)
            ie = ieuf(jj,l,1)
            iredo = 1
            if (is .ne. isave .or. ie .ne. ieave) then
              iredo = 0
              im = ie - is + 1
              isave = is
              ieave = ie
#  if defined O_cyclic
              if (im .ne. imtm2) then
                m = 2
                n = nint(im*csu(j)*csur(jfu0))
              else
                m = 3
                n = nint(im*csu(j)*csur(jfu0)*p5)
              endif
#  else
              m = 2
              n = nint(im*csu(j)*csur(jfu0))
#  endif
            endif
            ism1 = is - 1
            iea = ie
            if (ie .ge. imt) iea = imtm1
            do i=is,iea
              tempi(i-ism1,1) = -fx*uice(i,j)*spsin(i)
     &                        - vice(i,j)*spcos(i)
              tempi(i-ism1,2) =  fx*uice(i,j)*spcos(i)
     &                        - vice(i,j)*spsin(i)
            enddo
            if (ie .ge. imt) then
              ieb = ie - imtm2
              ii  = imtm1 - is
              do i=2,ieb
                tempi(i+ii,1) = -fx*uice(i,j)*spsin(i)
     &                        - vice(i,j)*spcos(i)
                tempi(i+ii,2) =  fx*uice(i,j)*spcos(i)
     &                        - vice(i,j)*spsin(i)
              enddo
            endif
#  if defined O_fourfil
            call filtr (tempi(1,1), im, m, n, iredo)
            call filtr (tempi(1,2), im, m, n, 1)
#  endif
#  if defined O_firfil
            if (im .eq. imt-2) then
#   if defined O_cyclic
              icyc = 1
#   else
              icyc = 0
#   endif
            else
              icyc = 0
            endif

!           filter then restore the original mean

            do n=1,2
              avgb = c0
              do i=1,im
                avgb = avgb + tempi(i,n)
              enddo
              avgb = avgb / im

              call filtra (tempi(1,n), im, numflu(jj), icyc)

              avga = c0
              do i=1,im
                avga = avga + tempi(i,n)
              enddo
              avga = avga / im

              do i=1,im
                tempi(i,n) = tempi(i,n) - avga + avgb
              enddo
            enddo

#  endif
            do i=is,iea
              uice(i,j) = fx*(-tempi(i-ism1,1)*spsin(i)
     &                   + tempi(i-ism1,2)*spcos(i))
              vice(i,j) = -tempi(i-ism1,1)*spcos(i)
     &                   - tempi(i-ism1,2)*spsin(i)
            enddo
            if (ie .ge. imt) then
              do i=2,ieb
                uice(i,j) = fx*(-tempi(i+ii,1)*spsin(i)
     &                     + tempi(i+ii,2)*spcos(i))
                vice(i,j) = -tempi(i+ii,1)*spcos(i)
     &                     - tempi(i+ii,2)*spsin(i)
              enddo
            endif
          endif
        enddo

        if (isave .ne. 0 .and. ieave .ne. 0) then

          do i=1,imt
            uice(i,j) = uice(i,j)*umsk(i,j)
            vice(i,j) = vice(i,j)*umsk(i,j)
          enddo

        endif

701     continue
      enddo

      call embmbc (uice)
      call embmbc (vice)
# endif

      return
      end
#endif
