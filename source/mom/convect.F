      subroutine convct (ts, ncon, joff, js, je, istrt, iend, kmt)

#if defined O_mom
# if !defined O_implicitvmix || defined O_isopycmix

      implicit none

      integer is, ie, nn, ncon, ks, js, je, j, jrow, joff, i, k, n
      integer iend, istrt

      real dense

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"

      integer kmt(imt,jmt)
      parameter (is=2, ie=imt-1)

      real ts(imt,km,1:jmw,nt), temp(imt,km,jsmw:jmw)

!-----------------------------------------------------------------------
!     standard explicit convection scheme
!     convectively adjust water column if gravitationally unstable

!     inputs:

!     ncon  = number of passes through convection routine
!     joff  = offset between "j" in MW and "jrow" latitude on disk
!     js    = starting row in MW
!     je    = ending row in MW
!     is    = starting longitude index
!     ie    = ending longitude index

!     Note: istrt,iend are currently bypassed. instead, is and ie are
!           set as parameters to optimize performance
!     kmt   = number of ocean "t" boxes in the vertical
!     ts    = temperature and salinity before convection

!     outputs:

!     ts    = tracers after convection

!-----------------------------------------------------------------------

!     ks=1: compare lev. 1 to 2; 3 to 4; etc.
!     ks=2: compare lev. 2 to 3; 4 to 5; etc.

      do nn=1,ncon
        do ks=1,2

!         find density for rows

          call statec (ts(1,1,1,1), ts(1,1,1,2), temp(1,1,jsmw)
     &,                max(js,jsmw), je, is, ie, ks)

!         set "heavy water" in land to stop convection

          dense = 1.e15
          do j=js,je
            jrow = j + joff
            do i=is,ie
              k = kmt(i,jrow) + 1
              if (k .le. km) then
                temp(i,k,j) = dense
              endif
            enddo
          enddo

!         if unstable,  mix tracers on adjoining levels

          do n=1,nt
            do j=js,je
              do k=ks,kmm1,2
                do i=is,ie
                  if (temp(i,k,j) .gt. temp(i,k+1,j)) then
                    ts(i,k,j,n)   = (dztxcl(k)*ts(i,k,j,n) +
     &                          dztxcl(k+1)*ts(i,k+1,j,n))*dzwxcl(k)
                    ts(i,k+1,j,n) = ts(i,k,j,n)
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      do n=1,nt
        do j=js,je
          call setbcx (ts(1,1,j,n), imt, km)
        enddo
      enddo

      return
      end

      subroutine convct2 (ts, joff, js, je, is, ie, kmt)

!=======================================================================
!     The following convection scheme is an alternative to the standard
!     scheme. In contrast to the standard scheme, it totally removes
!     all gravitational instability in the water column. It does that
!     in one pass, so the parameter ncon becomes irrelevant if this
!     option is selected. Since most convection propagates downward the
!     scheme looks downward first and follows any instability (upward or
!     downward) before checking the other direction. The routine mixes
!     passive tracers only after the entire instability is found. The
!     scheme is similar to that used by Rahmstorf (jgr 96,6951-6963) and
!     by Marotzke (jpo 21,903-907). It is discussed in a note to Ocean
!     Modelling (101). It uses as much cpu time as 1-3 passes of the
!     standard scheme, depending on the amount of static instability
!     found in the model, and is much faster than using "implicitvmix".

!     inputs:

!     joff   = offset between "j" in MW and "jrow" latitude on disk
!     js     = starting row in MW
!     je     = ending row in MW
!     is     = starting longitude index
!     ie     = ending longitude index
!     kmt    = number of ocean "t" boxes in the vertical
!     ts     = tracers before convection

!     outputs:

!     ts     = tracers after convection

!     other previously undefined variables:

!     chk_la = logical flag to check level above kt
!     chk_lb = logical flag to check level below kb
!     kb     = bottom level of (potential) instability
!     kbo    = bottom level of ocean
!     kt     = top level of (potential) instability
!     ktot   = total number of levels convecting in the column
!     kven   = number of levels that ventilated in the column
!     la     = test level above kt
!     lb     = test level below kb
!     rl     = lower level density referenced to lower level
!     ru     = upper level density referenced to lower level
!     tmx    = mixed tracer (1=temp, 2=salt, 3=other)
!     tsm    = sum of tracers (weighted by thickness) in the instability
!     zsm    = total thickness of the instability
!=======================================================================

      implicit none

      integer i, ie, is, j, je, joff, jrow, js, k, kb, kbo
      integer kt, ktot, kven, la, lb, n

      logical chk_la, chk_lb

      real rl, ru, tmx(3), tsm(3), zsm, dens, tq, sq, drodt, drods
      real drhodt, drhods, ddensdtdt, ddensdtds, ddensdsds

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
#  if defined O_save_convection
      include "coord.h"
      include "scalar.h"
      include "switch.h"
      include "diaga.h"
#  endif
      include "state.h"
      include "dens.h"

      integer kmt(imt,jmt)

      real ts(imt,km,1:jmw,nt)

      do j=js,je
        jrow = j + joff
        do i=is,ie
!         ktot = 0
          kbo = kmt(i,jrow)

#  if defined O_save_convection
          if (timavgperts) then
            totalk(i,j) = 0.0
            vdepth(i,j) = 0.0
            pe(i,j) = 0.0
            do k=1,km
              ru = dens (ts(i,k,j,1)-to(k), ts(i,k,j,2)-so(k), k)
              pe(i,j) = pe(i,j) + grav*zt(k)*ru*dztxcl(k)
            enddo
          endif

#  endif
!         search for unstable regions starting from the top
          kt = 1
          kb = 2
          do while (kt .lt. kbo)
            ru = dens (ts(i,kt,j,1)-to(kb), ts(i,kt,j,2)-so(kb), kb)
            rl = dens (ts(i,kb,j,1)-to(kb), ts(i,kb,j,2)-so(kb), kb)

!           sum the first pair found in an unstable region
            if (ru .gt. rl) then
              chk_la = .true.
              chk_lb = .true.
              zsm = dztxcl(kt) + dztxcl(kb)
              tsm(1) = ts(i,kt,j,1)*dztxcl(kt) + ts(i,kb,j,1)*dztxcl(kb)
              tmx(1) = tsm(1)/zsm
              tsm(2) = ts(i,kt,j,2)*dztxcl(kt) + ts(i,kb,j,2)*dztxcl(kb)
              tmx(2) = tsm(2)/zsm

              do while (chk_lb .or. chk_la)

!               check for an unstable level (lb) below kb
                if (kb .ge. kbo) chk_lb = .false.
                do while (chk_lb)
                  chk_lb = .false.
                  lb = kb + 1
                  ru = dens (tmx(1)-to(lb),      tmx(2)-so(lb),      lb)
                  rl = dens (ts(i,lb,j,1)-to(lb),ts(i,lb,j,2)-so(lb),lb)
                  if (ru .gt. rl) then
!                   add new level to sums
                    kb = lb
                   zsm = zsm + dztxcl(kb)
                    tsm(1) = tsm(1) + ts(i,kb,j,1)*dztxcl(kb)
                    tmx(1) = tsm(1)/zsm
                    tsm(2) = tsm(2) + ts(i,kb,j,2)*dztxcl(kb)
                    tmx(2) = tsm(2)/zsm
                    chk_la = .true.
                    if (kb .lt. kbo) chk_lb = .true.
                  endif
                enddo

!               check for an unstable level (la) above kt
! to get the equivalent of Rahmstorf's scheme uncomment the next line
                chk_la = .true.
                if (kt .le. 1) chk_la = .false.
                do while (chk_la)
                  chk_la = .false.
                  la = kt - 1
                  ru = dens (ts(i,la,j,1)-to(kt),ts(i,la,j,2)-so(kt),kt)
                  rl = dens (tmx(1)-to(kt),      tmx(2)-so(kt),      kt)
                  if (ru .gt. rl) then
!                   add new level to sums
                    kt = la
                    zsm = zsm + dztxcl(kt)
                    tsm(1) = tsm(1) + ts(i,kt,j,1)*dztxcl(kt)
                    tmx(1) = tsm(1)/zsm
                    tsm(2) = tsm(2) + ts(i,kt,j,2)*dztxcl(kt)
                    tmx(2) = tsm(2)/zsm
                    chk_lb = .true.
! to get the equivalent of Rahmstorf's scheme comment out the next line
!                    if (kt .gt. 1) chk_la = .true.
                  endif
                enddo
              enddo

!             mix all tracers from kt to kb
              do k=kt,kb
                ts(i,k,j,1) = tmx(1)
                ts(i,k,j,2) = tmx(2)
              enddo
              do n=3,nt
                tsm(3) = c0
                do k=kt,kb
                  tsm(3) = tsm(3) + ts(i,k,j,n)*dztxcl(k)
                enddo
                tmx(3) = tsm(3)/zsm
                do k=kt,kb
                  ts(i,k,j,n) = tmx(3)
                enddo
              enddo

!             some possible diagnostics
!              ktot = ktot + kb - kt + 1
!              if (kt .eq. 1) kven = kb
#  if defined O_save_convection

!             finished this region, save convection diagnostics
              if (timavgperts) then
                totalk(i,j) = totalk(i,j) + float(kb - kt + 1)
                if (kt .eq. 1) vdepth(i,j) = zw(kb)
              endif
#  endif

              kt = kb + 1
            else
              kt = kb
            endif

!           continue the search for other unstable regions
            kb = kt + 1
          enddo
#  if defined O_save_convection

          if (timavgperts) then
            do k=1,km
              ru = dens (ts(i,k,j,1)-to(k), ts(i,k,j,2)-so(k), k)
              pe(i,j) = pe(i,j) - grav*zt(k)*ru*dztxcl(k)
            enddo
            pe(i,j) = pe(i,j)/c2dtts
          endif

#  endif
        enddo

      enddo
# endif
#endif

      return
      end
