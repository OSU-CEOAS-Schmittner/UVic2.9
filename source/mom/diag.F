      subroutine diag (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     calculate diagnostics

!     input:

!      joff   = offset between row j in the MW and latitude jrow on disk
!      js     = starting row for calculations
!      je     = ending row for calculations
!      is     = starting longitude index for calculations
!      ie     = ending longitude index for calculations
!=======================================================================

      implicit none

      character(120) :: fname
      character(32) :: nstamp

      integer ntrec, nyear, nmonth, nday, nhour, nmin, nsec
      integer i, k, j, ip, kr, jq, js, je, istrt, is, iend, ie, joff
      integer jrow, n, jlat, jj, indp, ks, ke, m, io, i1, i2, iocm

      real zmau, zmat, zma1, zmsmf, zmsm, zma2, zmstf, zmst, reltim
      real fx, scl, period, ce, cn, cb, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
# if defined O_matrix_sections
      include "cprnts.h"
# endif
      include "cregin.h"
      include "diag.h"
      include "diaga.h"
      include "docnam.h"
      include "grdvar.h"
      include "iounit.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "vmixc.h"
      include "levind.h"
      include "emode.h"
# if defined O_npzd_out
      include "npzd.h"
# endif
# if defined O_embm
      include "cembm.h"
# endif

      real tmp_t(imt,km,nt), tmp_stf(imt,nt)
# if defined O_meridional_overturning
      real vbarx(km)
# endif
# if defined O_isopycmix
      real aibuf(imt,km)
# endif

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit longitudes
!-----------------------------------------------------------------------

      istrt  = max(2,is)
      iend   = min(imt-1,ie)

# if defined O_tai_otsf
      if (tsiperts .and. .not. euler2 .and. joff .eq. 0)
     &  nv_otsf = nv_otsf + 1
# endif
# if defined O_tai_slh
      if (tsiperts .and. .not. euler2 .and. joff .eq. 0)
     &  nt_slh = nt_slh + 1
# endif
      do j=js,je
        jrow = joff + j
# if defined O_time_step_monitor
!-----------------------------------------------------------------------
!       diagnostic: accumulate "tau" data for time step integrals
!-----------------------------------------------------------------------

        if (tsiperts .and. .not. euler2) then

#  if defined O_tai_otsf
          if (jrow .ge. jsot .and. jrow .le. jeot) then
            if (mrot .gt. 0 .and. mrot .le. nhreg) then
              do i=2,imtm1
                if (mskhr(i,jrow) .eq. mrot) then
                  do k=1,kmu(i,jrow)
                    v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                               dxu(i)
                  enddo
                endif
              enddo
            else
              do i=isot1,ieot1
                do k=1,kmu(i,jrow)
                  v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                             dxu(i)
                enddo
              enddo
              do i=isot2,ieot2
                do k=1,kmu(i,jrow)
                  v_otsf(jrow,k) = v_otsf(jrow,k) + u(i,k,j,2,tau)*
     &                             dxu(i)
                enddo
              enddo
            endif
          endif
#  endif
#  if defined O_tai_slh
          do i=1,imt
            do k=1,kmt(i,jrow)
              t_slh(i,jrow,k,1) = t_slh(i,jrow,k,1) + t(i,k,j,1,tau)
              t_slh(i,jrow,k,2) = t_slh(i,jrow,k,2) + t(i,k,j,2,tau)
            enddo
          enddo
#  endif

        endif
# endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!       diagnostic: accumulate "tau" data for time means
!-----------------------------------------------------------------------

        if (timavgperts .and. .not. euler2) then
          if (istrt .ne. 2 .and. iend .ne. imt-1) then
            write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "avgvar"'
            stop '=>diag'
          else
            call avgvar (j, jrow, adv_vbt(1,1,j), u(1,1,1,1,tau)
     &,                  t(1,1,1,1,tau), stf, smf, mapt)
          endif
        endif
# endif
# if defined O_stability_tests

!-----------------------------------------------------------------------
!       diagnostic: compute stability diagnostics
!-----------------------------------------------------------------------

        if (stabts .and. eots) then
          if (istrt .ne. 2 .and. iend .ne. imt-1) then
            write (stdout,*) '=>Error: istrt = ',istrt,' and iend ='
     &,     iend,' are not allowed when calling "stab"'
            stop '=>diag'
          else
            call stab (j, jrow)
          endif
        endif
# endif
# if defined O_meridional_overturning

!-----------------------------------------------------------------------
!       construct meridional overturning of mass
!-----------------------------------------------------------------------

        if (jrow .lt. jmtm1 .and. vmsfts .and. eots) then
          do k=1,km
            vbarx(k) = c0
          enddo

          do k=1,km
            do i=istrt,iend
              vbarx(k) = vbarx(k) + u(i,k,j,2,tau)*csu(jrow)*dxu(i)
            enddo
            if (k .eq. 1) then
              vmsf(jrow,k) = vbarx(k)*dzt(k)
            else
              vmsf(jrow,k) = vmsf(jrow,k-1) + vbarx(k)*dzt(k)
            endif
          enddo
        endif
# endif
# if defined O_show_zonal_mean_of_sbc

!-----------------------------------------------------------------------
!       construct zonal mean of surface b.c. and related items
!-----------------------------------------------------------------------

        if (zmbcts .and. eots) then
          zmau(jrow) = c0
          zmat(jrow) = c0
          do i=istrt,iend
            zma1            = umask(i,1,j)*csu(jrow)*dxu(i)*dyu(jrow)
            zmau(jrow)      = zmau(jrow)    + zma1
            zmsmf(jrow,1)   = zmsmf(jrow,1) + zma1*smf(i,j,1)
            zmsmf(jrow,2)   = zmsmf(jrow,2) + zma1*smf(i,j,2)
            zmsm(jrow,1)    = zmsm(jrow,1)  + zma1*u(i,1,j,1,tau)
            zmsm(jrow,2)    = zmsm(jrow,2)  + zma1*u(i,1,j,2,tau)
            zma2            = tmask(i,1,j)*cst(jrow)*dxt(i)*dyt(jrow)
            zmat(jrow)      = zmat(jrow) + zma2
            do n=1,nt
              zmstf(jrow,n) = zmstf(jrow,n) + zma2*stf(i,j,n)
              zmst(jrow,n)  = zmst(jrow,n)  + zma2*t(i,1,j,n,tau)
            enddo
          enddo
        endif
# endif
# if defined O_matrix_sections

!-----------------------------------------------------------------------
!       print "tau" (not "tau+1") variables  at specified latitudes
!-----------------------------------------------------------------------

        if (prxzts .and. eots) then
          reltim = relyr
          do jlat=1,nlatpr
            jj = indp (prlat(jlat), yt, jmt)
            if (jj .eq. jrow .and. prlat(jlat) .le. yt(jmt))
     &        then
              is = indp (prslon(jlat), xt, imt)
              ie = indp (prelon(jlat), xt, imt)
              ks = indp (prsdpt(jlat), zt, km)
              ke = indp (predpt(jlat), zt, km)
              fx = 1.0e-2
              if (jlat .eq. 1) write(stdout,9000)
              do m=1,nt
                scl = c1
                if (m .eq. 2) scl=1.e-3
                if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                  write (stdout,9100) trname(m), itt, jrow
     &,           yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                  call matrix (t(1,1,j,m,tau), imt, is, ie, ks, ke, scl)
                endif
                if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then

                  call getunit (io, 'sections.dta'
     &,                        'unformatted sequential append ieee')

                  write (stdout,*) ' =>Zonal ',trname(m), ' slice: lat='
     &,           yt(jrow), ' written unformatted to file sections.dta'
     &,           ' on ts=', itt, stamp
                  write (stdout,'(///)')
                  iotext = ' read (ioprxz) imt, km, m, nt, reltim'
                  write (io) stamp, iotext, expnam
                  write (io) imt, km, m, nt, reltim
                  write(iotext,'(a10,i4,a4,i2)') ' for jrow=',jrow
     &,           '  m=',m
                  iotext(18:)=':read(ioprxz)((t(i,k,m),i=1,imt),k=1,km)'
                  write (io) stamp, iotext, expnam
                  call wrufio (io, t(1,1,j,m,tau), imt*km)

                  call relunit (io)

                endif
              enddo
              scl = 1.e-3
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                i1 = max(is,2)
                i2 = min(ie,imtm1)
                write (stdout,9100)  'adv_vbu ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbu(1,1,j), imt, i1, i2, ks, ke, scl)
                write (stdout,9100)  'adv_vbt ', itt, jrow
     &,         yt(jrow), xt(i1), xt(i2), fx*zw(ks), fx*zw(ke), scl
                call matrix (adv_vbt(1,1,j), imt, i1, i2, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then

                call getunit (io, 'sections.dta'
     &,                      'unformatted sequential append ieee')

                write (stdout,*) ' => Zonal adv_vbu slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbu(i,k),i=1,imt),k=0,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbu(1,0,j), imt*(km+1))

                write (stdout,*) ' => Zonal adv_vbt slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:)=
     &          ': read(ioprxz)((adv_vbt(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, adv_vbt(1,0,j), imt*(km+1))

                call relunit (io)

              endif

              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'u velocity', itt
     &,        jrow, yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,1,tau), imt, is, ie, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then

                call getunit (io, 'sections.dta'
     &,                      'unformatted sequential append ieee')

                write (stdout,*) ' => Zonal u velocity slice: lat='
     &,         yt(jrow), ' written unformatted to file sections.dta'
     &,         ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((u(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,1,tau), imt*km)

                call relunit (io)

              endif

              scl = c1
              if (ioprxz .eq. stdout .or. ioprxz .lt. 0) then
                write (stdout,9100) 'v velocity', itt, jrow
     &,         yt(jrow), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
                call matrix (u(1,1,j,2,tau), imt, is, ie, ks, ke, scl)
              endif
              if (ioprxz .ne. stdout .or. ioprxz .lt. 0) then

                call getunit (io, 'sections.dta'
     &,                      'unformatted sequential append ieee')

                write (stdout,*) ' => Meridional v velocity slice: lat='
     &,         yt(jrow),' written unformatted to file sections.dta'
     &,          ' on ts=', itt, stamp
                write (stdout,'(///)')
                iotext = ' read (ioprxz) imt, km, reltim'
                write (io) stamp, iotext, expnam
                write (io) imt, km, reltim
                write(iotext,'(a10,i4)') ' for jrow=',jrow
                iotext(12:) = ': read (ioprxz)((v(i,k),i=1,imt),k=1,km)'
                write (io) stamp, iotext, expnam
                call wrufio (io, u(1,1,j,2,tau), imt*km)

                call relunit (io)

              endif
            endif
          enddo
        endif
9000    format(/' Zonal section printouts at specified latitudes:'/)
9100    format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     &,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     &,', scaling=',1pg10.3)
# endif
      enddo

# if defined O_save_mixing_coeff
!-----------------------------------------------------------------------
!     diagnostic: save estimated mixing coefficients on east, north, and
!                 bottom face of T and U cells
!-----------------------------------------------------------------------

      if (cmixts .and. eots) then

        reltim = relyr
        if (joff + js .eq. 2) then
          write (stdout,*) ' =>Writing mixing coefficients at ts=',itt
     &  , ' ',stamp
          call getunit (iocm, 'cmix.dta'
     &,                'unformatted sequential append ieee')

          period = 0.0
          iotext = 'read(iocm) reltim, period, imt, jmt, km'
          write (iocm) stamp, iotext, expnam
          write (iocm) reltim, period, imt, jmt, km

          iotext = 'read(iocm) (xt(i),i=1,imt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, xt, imt)

          iotext = 'read(iocm) (yt(j),j=1,jmt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, yt, jmt)

          iotext = 'read(iocm) (zt(k),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, zt, km)

          iotext = 'read(iocm) (xu(i),i=1,imt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, xu, imt)

          iotext = 'read(iocm) (yu(j),j=1,jmt)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, yu, jmt)

          iotext = 'read(iocm) (zw(k),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, zw, km)

          call relunit (iocm)
        endif

        call getunit (iocm, 'cmix.dta'
     &,               'unformatted sequential append ieee')

        do j=js,je
          jrow = j+joff

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_ceu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, ce(1,1,j,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cnu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cn(1,1,j,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cbu(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cb(1,1,j,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cet(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, ce(1,1,j,2), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cnt(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cn(1,1,j,2), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (diff_cbt(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          call wrufio (iocm, cb(1,1,j,2), imt*km)
#  if defined O_isopycmix

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (Ai_ez(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          do k=1,km
            do i=1,imt
              aibuf(i,k) = 0.25*(Ai_ez(i,k,j,0,1) + Ai_ez(i,k,j,0,0)
     &                         + Ai_ez(i,k,j,1,1) + Ai_ez(i,k,j,1,0))
            enddo
          enddo
          call wrufio (iocm, aibuf(1,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (Ai_nz(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          do k=1,km
            do i=1,imt
              aibuf(i,k) = 0.25*(Ai_nz(i,k,j,0,1) + Ai_nz(i,k,j,0,0)
     &                         + Ai_nz(i,k,j,1,1) + Ai_nz(i,k,j,1,0))
            enddo
          enddo
          call wrufio (iocm, aibuf(1,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (Ai_bx(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          do k=1,km
            do i=1,imt
              aibuf(i,k) = 0.25*(Ai_bx(i,k,j,0,0) + Ai_bx(i,k,j,1,0)
     &                         + Ai_bx(i,k,j,0,1) + Ai_bx(i,k,j,1,1))
            enddo
          enddo
          call wrufio (iocm, aibuf(1,1), imt*km)

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocm) (Ai_by(i,k),i=1,imt),k=1,km)'
          write (iocm) stamp, iotext, expnam
          do k=1,km
            do i=1,imt
              aibuf(i,k) = 0.25*(Ai_by(i,k,j,0,0) + Ai_by(i,k,j,1,0)
     &                         + Ai_by(i,k,j,0,1) + Ai_by(i,k,j,1,1))
            enddo
          enddo
          call wrufio (iocm, aibuf(1,1), imt*km)
#  endif
        enddo
        call relunit(iocm)
      endif
# endif
#endif

      return
      end
