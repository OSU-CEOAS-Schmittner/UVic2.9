      subroutine filuv (joff, js, je)

#if defined O_mom
# if defined O_fourfil || defined O_firfil
!=====================================================================
!     filuv sets up input needed for fourier filtering
!     (when the "fourfil" option is defined) -or- symmetric finite
!     impulse response filtering (when ifdef "firfil" is defined) of
!     baroclinic velocities at the specifiied high latitude row "jrow".
!=====================================================================

      implicit none

      integer n, j, js, je, jrow, joff, jj, isave, ieave, l, k, is, ie
      integer iredo, im, m, ism1, iea, i, ieb, ii, icyc

      real fx, avgb, avga

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cpolar.h"
      include "emode.h"
      include "grdvar.h"
      include "index.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

      real tempik(imt,km,2)

      do n=1,2
        do j=js,je
          call setbcx (u(1,1,j,n,taup1), imt, km)
        enddo
      enddo

#  if defined O_fourfil || defined O_firfil

!---------------------------------------------------------------------
!     fourier filter u and v at high latitudes
!---------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
      if ((jrow.gt.jfu1 .and. jrow.lt.jfu2) .or. jrow.lt.jfrst) goto 701
      jj = jrow - jfrst + 1
      if (jrow .ge. jfu2) jj = jj - jskpu + 1
      fx = -c1
      if (phi(jrow) .gt. c0) fx = c1
      isave = 0
      ieave = 0

      do l=1,lsegf
        do k=1,km
          if (isuf(jj,l,k) .ne. 0) then
            is = isuf(jj,l,k)
            ie = ieuf(jj,l,k)
            iredo = 1
            if (is.ne.isave .or. ie.ne.ieave) then
              iredo = 0
              im = ie - is + 1
              isave = is
              ieave = ie
#   if defined O_cyclic
              if (im .ne. imtm2) then
                m = 2
                n = nint(im*csu(jrow)*csur(jfu0))
              else
                m = 3
                n = nint(im*csu(jrow)*csur(jfu0)*p5)
              endif
#   else
              m = 2
              n = nint(im*csu(jrow)*csur(jfu0))
#   endif
            endif
            ism1 = is - 1
            iea = ie
            if (ie .ge. imt) iea = imtm1
            do i=is,iea
              tempik(i-ism1,k,1) = -fx*u(i,k,j,1,taup1)*spsin(i)
     &                             - u(i,k,j,2,taup1)*spcos(i)
              tempik(i-ism1,k,2) =  fx*u(i,k,j,1,taup1)*spcos(i)
     &                             - u(i,k,j,2,taup1)*spsin(i)
            enddo
            if (ie .ge. imt) then
              ieb = ie - imtm2
              ii  = imtm1 - is
              do i=2,ieb
                tempik(i+ii,k,1) = -fx*u(i,k,j,1,taup1)*spsin(i)
     &                             - u(i,k,j,2,taup1)*spcos(i)
                tempik(i+ii,k,2) =  fx*u(i,k,j,1,taup1)*spcos(i)
     &                            - u(i,k,j,2,taup1)*spsin(i)
              enddo
            endif
#  if defined O_fourfil
            call filtr (tempik(1,k,1), im, m, n, iredo)
            call filtr (tempik(1,k,2), im, m, n, 1)
#  endif
#  if defined O_firfil
            if (im .eq. imt-2) then
#    if defined O_cyclic
              icyc = 1
#    else
              icyc = 0
#    endif
            else
              icyc = 0
            endif

!           filter then restore the original mean

            do n=1,2
              avgb = c0
              do i=1,im
                avgb = avgb + tempik(i,k,n)
              enddo
              avgb = avgb / im

              call filtra (tempik(1,k,n), im, numflu(jj), icyc)

              avga = c0
              do i=1,im
                avga = avga + tempik(i,k,n)
              enddo
              avga = avga / im

              do i=1,im
                tempik(i,k,n) = tempik(i,k,n) - avga + avgb
              enddo
            enddo

#  endif
            do i=is,iea
              u(i,k,j,1,taup1) = fx*(-tempik(i-ism1,k,1)*spsin(i)
     &                   + tempik(i-ism1,k,2)*spcos(i))
              u(i,k,j,2,taup1) = -tempik(i-ism1,k,1)*spcos(i)
     &                   - tempik(i-ism1,k,2)*spsin(i)
            enddo
            if (ie .ge. imt) then
              do i=2,ieb
                u(i,k,j,1,taup1) = fx*(-tempik(i+ii,k,1)*spsin(i)
     &                     + tempik(i+ii,k,2)*spcos(i))
                u(i,k,j,2,taup1) = -tempik(i+ii,k,1)*spcos(i)
     &                     - tempik(i+ii,k,2)*spsin(i)
              enddo
            endif
          endif
        enddo
      enddo

      if (isave .ne. 0 .and. ieave .ne. 0) then
      do i=1,imt
        tempik(i,1,1) = c0
        tempik(i,1,2) = c0
      enddo

      do k=1,km
        do i=1,imt
          tempik(i,1,1) = tempik(i,1,1) + u(i,k,j,1,taup1)*dzt(k)
          tempik(i,1,2) = tempik(i,1,2) + u(i,k,j,2,taup1)*dzt(k)
        enddo
      enddo

      do i=1,imt
        tempik(i,1,1) = tempik(i,1,1)*hr(i,jrow)
        tempik(i,1,2) = tempik(i,1,2)*hr(i,jrow)
      enddo

      do k=1,km
        do i=1,imt
          u(i,k,j,1,taup1) = u(i,k,j,1,taup1) - tempik(i,1,1)
          u(i,k,j,2,taup1) = u(i,k,j,2,taup1) - tempik(i,1,2)
        enddo
      enddo

      do k=1,km
        do i=1,imt
          u(i,k,j,1,taup1) = u(i,k,j,1,taup1)*umask(i,k,j)
          u(i,k,j,2,taup1) = u(i,k,j,2,taup1)*umask(i,k,j)
        enddo
      enddo
      endif

701   continue
      enddo
#  endif
# endif
#endif

      return
      end

#if defined O_firfil
      subroutine filtra (a, im, num, icyc)

!-----------------------------------------------------------------------
!     finite impulse response filter used for velocity components.
!     boundary condition is zero velocity on walls

!     input:

!     a    = velocity component
!     im   = length of "a" (may be less than "imt")
!     num  = number of filterings (each is a double pass)
!     icyc = (0,1) = (do not apply, apply) cyclic condition
!            if im = "imt-2" then icyc=1 when cyclic option is used

!     output:

!     a = filtered velocity component

!     based on code by: R. C. Pacanowski
!-----------------------------------------------------------------------

      implicit none

      integer im, num, n, icyc, i

      include "size.h"

      real a(im), s(imt)

      if (num .eq. 0) return

      do n=1,num
        if (icyc .eq. 1) then
          s(1)  = 0.25*a(im)   + 0.5*a(1)  + 0.25*a(2)
          s(im) = 0.25*a(im-1) + 0.5*a(im) + 0.25*a(1)
        else
          s(1)  = 0.5*a(1) + 0.25*a(2)
          s(im) = 0.5*a(im) + 0.25*a(im-1)
        endif
        do i=2,im-1
          s(i) = 0.25*a(i-1) + 0.5*a(i) + 0.25*a(i+1)
        enddo

        if (icyc .eq. 1) then
          a(1)  = 0.25*s(im)   + 0.5*s(1)  + 0.25*s(2)
          a(im) = 0.25*s(im-1) + 0.5*s(im) + 0.25*s(1)
        else
          a(1)  = 0.5*s(1) + 0.25*s(2)
          a(im) = 0.5*s(im) + 0.25*s(im-1)
        endif
        do i=2,im-1
          a(i) = 0.25*s(i-1) + 0.5*s(i) + 0.25*s(i+1)
        enddo
      enddo

      return
      end
#endif
