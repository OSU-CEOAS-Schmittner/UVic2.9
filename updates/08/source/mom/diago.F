      subroutine diago

#if defined O_mom
!=======================================================================
!     write out diagnostic output
!=======================================================================

      implicit none

      integer numdsp, jrow, i, nislsp, mdscan, is, indp, ie, js, je
      integer io, m, mask, k, msk, num, n, iobadt, iobads, jte, jue
      integer jwte, jwue, ll, nv, ks, ke, iv, ih, lll, maxm, mloop
      integer ms, me, l, len, j, jj, iou

      logical defined

      real uext, vext, rnum, npt, dspcrt, rgrav, scl, reltim, rvolgk
      real rvolgt, rareag, avgper, fx, tnew, tsml, tbig, plicin, plicex
      real buoerr, enleak, dtconv, dtfilt, taux, tauy, erru, errv
      real contu, dchg, cont, tperiod, pwatts, csalt
      real tmbavg, convrt, sumsto, sumdiv, sumflx, sumdif, sumsor
      real sumvol, terror, serror, cmmday, cwatts, zmau, zmsmf,  zmsm
      real zmat, zmstf, zmst, tconv, tfilt, rnavg, timunit, tmp

      character(120) :: fname, ftbt, file_stamp, new_file_name
      character(32) :: nstamp
      save ftbt
      data ftbt  /' '/

      integer ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "ctmb.h"
      include "coord.h"
      include "cprnts.h"
      include "cregin.h"
      include "diag.h"
      include "docnam.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "stab.h"
      include "state.h"
      include "switch.h"
      include "tmngr.h"
      include "vmixc.h"
# if defined O_embm
      include "cembm.h"
# endif
# if defined O_tracer_averages
      include "ctavg.h"
# endif
# if defined O_save_carbon_carbonate_chem
      include "timeavgs.h"
      include "diaga.h"
# endif
      include "npzd.h"
# if defined O_mom_tbt
      include "tbt.h"

      real errt(nt), tm(imt,jmt,km)
# endif
# if defined O_diagnostic_surf_height && defined O_stream_function
      character(8) :: bc_symm

!     psgrad= gradient of sea surface pressure (1,2) = (u,v) component
!     dsp   = diagnostic sea surface pressure converted to height
!     divf  = divergence of sea surface pressure gradients

#  if !defined O_rigid_lid_surface_pressure && !defined O_implicit_free_surface
      real divf(imt,jmt)
      save divf
#  endif
      real psgrad(imt,jmt,2), dsp(imt,jmt)
      save dsp, numdsp

!-----------------------------------------------------------------------
!     reconstruct the diagnostic surface pressure
!-----------------------------------------------------------------------

!     initialize the surface pressure fields and averaging counter

      if (first) then
        numdsp = 0
        do jrow=1,jmt
          do i=1,imt
            divf(i,jrow) = c0
            dsp(i,jrow)  = c0
            res(i,jrow)  = c0
          enddo
        enddo
      endif

      if (dspperts) then

!       accumulate forcing for diagnostic surface height over the
!       prescribed averaging period "dspper"

!       increment the counter for averaging the divergence

        numdsp = numdsp + 1

!       construct surface pressure gradients
!       (note: uext & vext are returned but not needed)

        call calc_psgrad (psgrad, uext, vext, 1, jmt, 1, imt)

        do i=1,imt
          psgrad(i,jmt,1) = c0
          psgrad(i,jmt,2) = c0
          psgrad(i,1,1)   = c0
          psgrad(i,1,2)   = c0
        enddo
        call setbcx (psgrad(1,1,1), imt, jmt)
        call setbcx (psgrad(1,1,2), imt, jmt)

!       compute the divergence of the sea surface pressure gradients

        call spforc (psgrad, dxu, dyu, csu, h, res)

!       accumulate the divergence

        do jrow=1,jmt
          do i=1,imt
            divf(i,jrow) = divf(i,jrow) + res(i,jrow)
          enddo
        enddo
      endif
      if (dspts) then

        write (stdout,'(///,40x,a,/)')
     & 'D I A G O N O S T I C    S U R F A C E    H E I G H T'

!       average the divergence, zero the guess

        rnum = c1/numdsp
        do jrow=1,jmt
          do i=1,imt
            divf(i,jrow) = rnum*divf(i,jrow)
            dsp(i,jrow)  = c0
          enddo
        enddo

!       initialize coefficients for the conjugate gradient solver

        call spc9pt (dxu, dyu, csu, h, cf)
        npt        = 9
        variable   = 'surfpres'
        bc_symm    = 't even'

!       number of islands must be zero for surface pressure
!       dspcrt = tolerence of solution for dsp

        nislsp = 0
        dspcrt = tolrsp

        call congr (npt, variable, bc_symm, dsp, dsp, divf, res
     &,           cf
     &,           mxscan, mdscan, dspcrt
     &,           imask, iperm, jperm, iofs, nislsp, nippts
     &,           converged, esterr)

        write (stdout,'(/a,i5,a,e14.7)')
     &  '=> diagnostic surface pressure scans=',mdscan,' estimated err='
     &, esterr

        if (.not.converged) then
          write (stdout,'(/,a,/)')
     & '=> Warning: convergence not reached for diag surface pressure'
          converged = .true.
        endif

!       subtract the null space from the surface pressure

        call checkerboard (dsp, map)
        call border (dsp, 't even')

!       remove mean from the surface pressure

        call zero_level (dsp, 'surf press', map, dxt, dyt, cst)
        call border (dsp, 't even')

!       convert surface pressure to surface elevation and write out

        rgrav = c1/(rho0*grav)
        do jrow=1,jmt
          do i=1,imt
            dsp(i,jrow) = rgrav*dsp(i,jrow)
          enddo
        enddo

        if (iodsp .eq. stdout .or. iodsp .lt. 0) then
          is  = indp (slonxy, xt, imt)
          ie  = indp (elonxy, xt, imt)
          js  = indp (slatxy, yt, jmt)
          je  = indp (elatxy, yt, jmt)
          scl = c1
          write (stdout,'(/1x,a,i5,a)')
     &   'Diagnostic surf height is averaged over ',numdsp,' time steps'
          write (stdout,7100) 'Diagnostic surface height (cm)'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl

          call matrix (dsp(1,1), imt, is, ie, -js, -je, scl)

7100      format(1x,a30,1x,'ts=',i7
     &,    ', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,    ', scaling=',1pg10.3)

        endif
        if (iodsp .ne. stdout .or. iodsp .lt. 0) then

          reltim = prelyr

          call getunit (io, 'diag_surf.dta'
     &,                'unformatted sequential append ieee')

          iotext = 'read (iodsp) imt, jmt, reltim, dspper, xt, yt'
          write (io) pstamp, iotext, expnam
          write (io) imt, jmt, reltim, dspper, xt, yt

          iotext = 'read (iodsp) ((dsp(i,j),i=1,imt),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, dsp, imt*jmt)

          call relunit (io)

          write (stdout,'(a,a,f8.3,a,a,i10,a,a)')
     &    ' => Diagnostic surface height (averaged '
     &,   'over ',dspper,' days) written '
     &,  'unformatted to file diag_surf.dta on ts=', itt, ' ',stamp
        endif

!       reset counter for next average. zero the averaged divergence

        numdsp = 0
        do jrow=1,jmt
          do i=1,imt
            divf(i,jrow) = c0
          enddo
        enddo
      endif
# endif
# if defined O_tracer_averages

!-----------------------------------------------------------------------
!     compute tracer averages under horizontal regions and write them
!     out
!-----------------------------------------------------------------------

      if (tavgts) then

!       initialize sums for horizontal regions & tracer averages

        do m=1,nt
          sumgf(m) = c0
          avggf(m) = c0
          sumgt(m) = c0
          avggt(m) = c0
          do mask=1,nhreg
            sumbt(mask,m) = c0
            avgbt(mask,m) = c0
            avgbf(mask,m) = c0
          enddo
          do k=1,km
            sumgk(k,m) = c0
            avggk(k,m) = c0
          enddo
        enddo

!       compute sums for tracer averages over horizontal regions

        do m=1,nt
          do mask=1,nhreg
            sumgf(m) = sumgf(m) + sumbf(mask,m)
            do k=1,km
              sumbt(mask,m) = sumbt(mask,m) + sumbk(mask,k,m)
              sumgk(k,m) = sumgk(k,m) + sumbk(mask,k,m)
            enddo
            sumgt(m) = sumgt(m) + sumbt(mask,m)
          enddo
        enddo

        do k=1,km
          if (volgk(k) .gt. c0) then
            rvolgk = c1 / volgk(k)
            do m=1,nt
              avggk(k,m) = sumgk(k,m) * rvolgk
              do mask=1,nhreg
                if (volbk(mask,k) .gt. c0) then
                  avgbk(mask,k,m) = sumbk(mask,k,m) / volbk(mask,k)
                endif
              enddo
            enddo
          endif
        enddo

        rvolgt = c1 / volgt
        rareag = c1 / areag
        do m=1,nt
          avggt(m) = sumgt(m) * rvolgt
          avggf(m) = sumgf(m) * rareag
          do mask=1,nhreg
            if (volbt(mask) .gt. c0) then
              avgbt(mask,m) = sumbt(mask,m) / volbt(mask)
            endif
            if (areab(mask) .gt. c0) then
              avgbf(mask,m) = sumbf(mask,m) / areab(mask)
            endif
          enddo
        enddo

!       write out regional tracer means

        if (iotavg .eq. stdout .or. iotavg .lt. 0) then
          write (stdout,'(//,30x,a,/)') 'T R A C E R    A V E R A G E S'
          do m=1,nt
            write(stdout,9004) trname(m), itt, stamp
            write(stdout,9001) (mask,mask=1,nhreg)
            do k=1,km
              write(stdout,9002) k, avggk(k,m),
     &                             (avgbk(mask,k,m),mask=1,nhreg)
            enddo
            write(stdout,9003) avggt(m), (avgbt(mask,m),mask=1,nhreg)
            write(stdout,9014) m, avggf(m), (avgbf(msk,m),msk=1,nhreg)
          enddo
        endif
        if (iotavg .ne. stdout .or. iotavg .lt. 0) then
          reltim = prelyr

          call getunit (io, 'tracer_avg.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Regional tracer averages written'
     &,      ' unformatted to file tracer_avg.dta on ts=',itt, ' ',stamp

          iotext = 'read (iotavg) reltim, nt, nhreg, km'
          write (io) pstamp, iotext, expnam
          write (io) reltim, nt, nhreg, km

          iotext = 'read (iotavg) ((avggk(k,n),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggk, km*nt)

          iotext =
     &    'read (iotavg) (((avgbk(l,k,n),l=1,nhreg),k=1,km),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbk, nhreg*km*nt)

          iotext = 'read (iotavg) (avggt(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggt, nt)

          iotext = 'read (iotavg) (avggf(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avggf, nt)

          iotext = 'read (iotavg) ((avgbt(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbt, nhreg*nt)

          iotext = 'read (iotavg) ((avgbf(l,n),l=1,nhreg),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgbf, nhreg*nt)

          call relunit (io)
        endif
      endif
9001  format('    k','  All Regions ',9(1x,i7,5x))
9002  format(1x,i4,10(1x,e12.5))
9003  format('  AVG',10(1x,e12.5))
9004  format(/' Volume Weighted Averages for ',a12,' on ts =',i10,a33)
9014  format(/' FLX',i1,10(1x,e12.5),/)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     print integrals for monitoring the time step
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) call ta_mom_tsi (1)

      if (tsits .and. ntatio .gt. 0) then
        call ta_mom_tsi (2)

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = time
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (stdout,9601) itt, nstamp, tai_ek, tai_dt, tai_ds, tai_t
     &,                       tai_s, tai_tvar, tai_svar, tai_scan
        endif

        call def_tsi
        call def_tsi_mom (fname)
#  if defined O_save_carbon_totals
!       convert from umol cm-3 to Pg (redctn contained 1e-3)
        tai_cocn = (tai_dic 
#   if defined O_mobi
     &       + (tai_phyt + tai_zoop + tai_detr
#    if defined O_mobi_nitrogen
     &       + tai_diaz
#    endif
     &       )*redctn
#   endif
     &       )*tcellv*12.e-21
#   if defined O_mobi_caco3
        tai_cocn = tai_cocn + tai_caco3*tcellv*12.e-21
#   endif
#   if defined O_kk_ballast
        tai_cocn = tai_cocn + tai_d_B*redctn*tcellv*12.e-21
#   endif
#   if defined O_mobi_silicon
        tai_cocn = tai_cocn + tai_diat*redctn*tcellv*12.e-21
#   endif        
!       convert from umol cm-2 s-1 to Pg year-1
        tai_cfa2o = tai_dicflx*tcella(1)*yrlen*daylen*12.e-21
#   if defined O_sed && !defined O_sed_uncoupled
!       subtract out any land to ocean flux
        tai_cfa2o = tai_cfa2o - tai_dicwflx*yrlen*daylen*12.e-21
#   endif
#  else
        tai_cocn = 0.
        tai_cfa2o = 0.
#  endif

#  if !defined O_save_time_relyear0
!       make output time relative to year 1
        time = time - 1.
#  endif
        avgper = tsiper*accel
        if (avgper .le. 1e-6) avgper = 0.
#  if defined O_save_time_endper
        tmp = 0.
#  elif defined O_save_time_startper
        tmp = 1.
#  else
        tmp = 0.5
#  endif
#  if defined O_units_time_years
#   if defined O_calendar_360_day
        time = time - tmp*avgper/360.
#   elif defined O_calendar_gregorian
        time = time - tmp*avgper/365.25
#   else
        time = time - tmp*avgper/365.
#   endif
#  else
#   if defined O_calendar_360_day
        time = time*360. - tmp*avgper
#   elif defined O_calendar_gregorian
        time = time*365.25 - tmp*avgper
#   else
        time = time*365. - tmp*avgper
#   endif
#  endif

        call mom_tsi_out (fname, avgper, time, nstamp, tai_ek, tai_t
     &,                   tai_s, tai_tvar, tai_svar, tai_dt, tai_ds
     &,                   tai_scan, tai_hflx, tai_sflx, tai_dic
     &,                   tai_dicflx, tai_dic13, tai_dic13flx 
     &,                   tai_alk, tai_o2, tai_o2flx
     &,                   tai_po4, tai_dop, tai_phyt, tai_zoop, tai_detr
     &,                   tai_no3, tai_don, tai_diaz, tai_din15
     &,                   tai_don15, tai_phytn15, tai_zoopn15
     &,                   tai_detrn15, tai_diazn15, tai_doc13 
     &,                   tai_phytc13, tai_zoopc13, tai_detrc13
     &,                   tai_diazc13, tai_c14, tai_dc14
     &,                   tai_c14flx, tai_cfc11, tai_cfc11flx, tai_cfc12
     &,                   tai_cfc12flx, tai_otmax, tai_otmin
     &,                   tai_slh, tai_sspH, tai_ssCO3, tai_ssOc
     &,                   tai_ssOa, tai_sspCO2, tai_cocn, tai_cfa2o
     &,                   ntrec, tai_dfe, tai_ddfe
#  if defined O_kk_ballast
     &,                   tai_d_B
#  endif
#  if defined O_mobi_caco3
     &,                   tai_caco3, tai_caco3c13
#  endif
#  if defined O_mobi_silicon
     &,                   tai_diat
     &,                   tai_diatn15, tai_diatc13 
     &,                   tai_sil, tai_silflx, tai_opl
#  endif
     &       )

        call ta_mom_tsi (0)
      endif
9601  format (1x,'ts=',i7, 1x, a32, ', ke=', 1pe13.6,' |dT|=',1pe13.6
     &,      ' |dS|=',1pe13.6,' Tbar=',1pe13.6,' Sbar=',1pe13.6
     &,      ' Tvar=',1pe13.6,' Svar=',1pe13.6, ' scans=',1pe13.6)
# endif
# if defined O_stability_tests

!-----------------------------------------------------------------------
!     show stability and CFL conditions, reynolds and peclet numbers
!-----------------------------------------------------------------------

      if (stabts) then
        write (stdout,'(///20x,a/15x,a,/)')
     & ' S T A B I L I T Y     A N A L Y S I S'
     &, '(The indicated locations are the most unstable ones)'

        write (stdout,*) ' longitudinal domain: ',cflons, ' to ',cflone
        write (stdout,*) ' latitudinal  domain: ',cflats, ' to ',cflate
        write (stdout,*) ' depth domain (m)   : ',cfldps*0.01
     &,                  ' to ',cfldpe*0.01

        write (stdout,'(/60x,a/)') ' CFL summary'
        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local U velocity (',cflum,') is ',cflup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflu,',',jcflu,',',kcflu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflu),',',yu(jcflu),','
     &,  0.01*zt(kcflu),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     &   ' Local V velocity (',cflvm,') is ',cflvp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflv,',',jcflv,',',kcflv,')'
     &,  ' (lon,lat,dpt) = (',xu(icflv),',',yu(jcflv),','
     &,  0.01*zt(kcflv),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbu    (',cflwum,') is ',cflwup,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwu,',',jcflwu,',',kcflwu,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwu),',',yu(jcflwu),','
     &,  0.01*zw(kcflwu),')'

        write (stdout
     &,'(a,g10.3,a,f7.2,a,/a,i4,a,i4,a,i3,a,3x,a,f7.2,a,f7.2,a,f7.0,a)')
     & ' Local adv_vbt    (',cflwtm,') is ',cflwtp,' % of the CFL limit'
     &,  ' at location: (i,j,k) = (',icflwt,',',jcflwt,',',kcflwt,')'
     &,  ' (lon,lat,dpt) = (',xu(icflwt),',',yu(jcflwt),','
     &,  0.01*zw(kcflwt),')'

        fx = 100.0
        if (cflup .gt. fx .or. cflvp .gt. fx .or. cflwup .gt. fx .or.
     &      cflwtp .gt. fx) then
          write (stdout,*)
     &      ' => Warning. CFL exceeded... computational mode exists!'
        endif

        write (stdout,'(/60x,a24/)') ' Reynolds number summary'

        write (stdout,10300) reynx, ireynx, jreynx, kreynx
     &,                      xu(ireynx), yu(jreynx), 0.01*zt(kreynx)
        write (stdout,10310) reynu, reynmu

        write (stdout,10400) reyny, ireyny, jreyny, kreyny
     &,                      xu(ireyny), yu(jreyny), 0.01*zt(kreyny)
        write (stdout,10410) reynv, reynmv

        write (stdout,10500) reynz, ireynz, jreynz, kreynz
     &,                      xu(ireynz), yu(jreynz), 0.01*zt(kreynz)
        write (stdout,10510) reynw, reynmw

        write (stdout,'(/60x,a22/)') ' Peclet number summary'

        write (stdout,10600) peclx, ipeclx, jpeclx, kpeclx
     &,                      xu(ipeclx), yu(jpeclx), 0.01*zt(kpeclx)
        write (stdout,10610) peclu, peclmu

        write (stdout,10700) pecly, ipecly, jpecly, kpecly
     &,                      xu(ipecly), yu(jpecly), 0.01*zt(kpecly)
        write (stdout,10710) peclv, peclmv

        write (stdout,10800) peclz, ipeclz, jpeclz, kpeclz
     &,                      xu(ipeclz), yu(jpeclz), 0.01*zt(kpeclz)
        write (stdout,10810) peclw, peclmw

!       show ficticious tracer extremums

        call getunit (iostab, 'iostab', 'formatted sequential rewind')
        rewind iostab
        write (stdout,'(/,10x,a/,11x,a,1pe10.3,a/)')
     &  'Spurious creation of local tracer extremum (if any) follow:'
     &, '(where tracer exceeds local extremum by ',tdig,'*tracer)'
        do num=1,1000
          read (iostab,'(i4, i4, i4, i2, 3g14.7)', end=101, err=101)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    't(i,k,jrow,n) = t(',i,',',k,',',jrow,',',n, ') = '
     &,   tnew,' : local min was ',tsml, ' : local max was ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 101
          endif
        enddo
101     continue
        call relunit (iostab)

!       show T and S outside allowable bounds used for density coeffs

        call getunit (iobadt, 'iobadt', 'formatted sequential rewind')
        rewind iobadt
        write (stdout,'(/,10x,a/)')
     &  'Temperatures (if any) outside allowable ranges follow:'
        do num=1,1000
          read (iobadt,'(i4, i4, i4, i2, 3g14.7)', end=201, err=201)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    't(i,k,jrow,n) = t(',i,',',k,',',jrow,',',n, ') = '
     &,   tnew,' : tmin is ',tsml, ' : tmax is ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 201
          endif
        enddo
201     continue
        call relunit (iobadt)

        call getunit (iobads, 'iobads', 'formatted sequential rewind')
        rewind iobads
        write (stdout,'(/,10x,a/)')
     &  'Salinities (if any) outside allowable ranges follow:'
        do num=1,1000
          read (iobads,'(i4, i4, i4, i2, 3g14.7)', end=301, err=301)
     &           i, k, jrow, n, tnew, tsml, tbig
          write (stdout,'(1x,a,i4,a,i4,a,i4,a,i2,3(a,g14.7))')
     &    'Error condition: t(i,k,jrow,n) = t(',i,',',k,',',jrow
     &,   ',',n, ') = ', tnew,' : smin is ',tsml, ' : smax is ', tbig
          if (num .eq. 100) then
            write (stdout,'(/a/)') 'Bailing out after 100 locations...'
            go to 301
          endif
        enddo
301     continue
        call relunit (iobads)
        write (stdout,'(/60x,a/)') ' End Stability Analysis'
      endif
10300 format (1x,'Maximim zonal Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10310 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10400 format (1x,'Maximim meridional Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10410 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10500 format (1x,'Maximim vertical Reynolds number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10510 format (1x,' local Wu =',1pe9.2, ' and  mixing =',e9.2)
10600 format (1x,'Maximim zonal Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10610 format (1x,' local U =',1pe9.2, ' and  mixing =',e9.2)
10700 format (1x,'Maximim meridional Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10710 format (1x,' local V =',1pe9.2, ' and  mixing =',e9.2)
10800 format (1x,'Maximim vertical Peclet number is ',1pe9.2
     &,       ' at location: (i,j,k) = (',i4,',',i4,',',i4,'),'
     &,       ' (lon,lat,dpt) = (',e9.2,',',e9.2,',',e9.2,')')
10810 format (1x,' local Wt =',1pe9.2, ' and  mixing =',e9.2)
# endif
# if defined O_energy_analysis

!-----------------------------------------------------------------------
!     add external mode component of total work done
!-----------------------------------------------------------------------

      if (glents) then
        call ge3 (c2dtuv)
      endif

!-----------------------------------------------------------------------
!     integrate previously computed energy components vertically
!-----------------------------------------------------------------------

      if (glents) then
        jte  = 1
        jue  = 1
        jwte = 1
        jwue = 1
        do k=1,km
          wtlev(k,0) = c0
          wulev(k,0) = c0
        enddo
        do jrow=2,jmt-1
          do k=km,1,-1
            buoy(0,1) = buoy(0,1) + buoy(k,jrow)
            wtlev(k,0) = wtlev(k,0) + wtlev(k,jrow)
            wulev(k,0) = wulev(k,0) + wulev(k,jrow)
          enddo
          do ll=1,8
            do k=km,1,-1
              engint(0,ll,1) = engint(0,ll,1) + engint(k,ll,jrow)
            enddo
            engext(ll,1) = engext(ll,1) + engext(ll,jrow)
          enddo
          if (abs(tcerr(jrow)) .gt. abs(tcerr(jte))) jte = jrow
          if (abs(ucerr(jrow)) .gt. abs(ucerr(jte))) jue = jrow
          if (abs(wtbot(jrow)) .gt. abs(wtbot(jwte))) jwte = jrow
          if (abs(wubot(jrow)) .gt. abs(wubot(jwue))) jwue = jrow
        enddo
        buoy(0,1) = buoy(0,1)/ucellv
        do ll=1,8
          engint(0,ll,1) = engint(0,ll,1)/ucellv
          engext(ll,1)   = engext(ll,1)/ucellv
        enddo

        plicin = engint(0,1,1) - engint(0,2,1) - engint(0,3,1)
     &           - engint(0,4,1) - engint(0,5,1) - engint(0,6,1)
        plicex = engext(1,1) - engext(2,1) - engext(3,1)
     &            - engext(4,1)  - engext(5,1) - engext(6,1)
        buoerr = buoy(0,1) - engint(0,6,1) - engext(6,1)
        enleak = engint(0,2,1) + engint(0,3,1) + engext(2,1)
     &         + engext(3,1)

        if (ioglen .eq. stdout .or. ioglen .lt. 0) then
          write (stdout,'(///40x,a,/)') 'E N E R G Y    A N A L Y S I S'
          write (stdout,9100)
     &              'Globally averaged work done'
     &,                     ' for ts =', itt, stamp, ucellv, ucella(1)
          write (stdout,9101) ' time rate of change ',engint(0,1,1)
     &,                      engint(0,1,1), engext(1,1), engext(1,1)
          write (stdout,9101) ' horizontal advection', engint(0,2,1)
     &,                      engint(0,2,1), engext(2,1), engext(2,1)
          write (stdout,9101) ' vertical advection  ',engint(0,3,1)
     &,                      engint(0,3,1), engext(3,1), engext(3,1)
          write (stdout,9101) ' horizontal friction ',engint(0,4,1)
     &,                      engint(0,4,1), engext(4,1), engext(4,1)
          write (stdout,9101) ' vertical friction   ',engint(0,5,1)
     &,                      engint(0,5,1), engext(5,1), engext(5,1)
          write (stdout,9101) ' pressure forces     ',engint(0,6,1)
     &,                      engint(0,6,1), engext(6,1), engext(6,1)
          write (stdout,9101) ' ficticious sources  ',plicin
     &,                       plicin, plicex, plicex
          write (stdout,9101) ' work by wind        ',engint(0,7,1)
     &,                      engint(0,7,1), engext(7,1), engext(7,1)
          write (stdout,9101) ' bottom drag         ',engint(0,8,1)
     &,                      engint(0,8,1), engext(8,1), engext(8,1)
          write (stdout,9110) buoy(0,1), buoy(0,1), buoerr, buoerr
     &,                      enleak, enleak
          write (stdout,9111) tcerr(jte), itcerr(jte), jtcerr(jte)
     &,                       ktcerr(jte)
     &,                       ucerr(jue), iucerr(jue), jucerr(jue)
     &,                       kucerr(jue)
          write (stdout,9112) wtbot(jwte), iwtbot(jwte), jwtbot(jwte)
     &,                       kwtbot(jwte)
     &,                       wubot(jwue), iwubot(jwue), jwubot(jwue)
     &,                       kwubot(jwue)

          write (stdout,'(/a,a,//a,a,a,a,a,/1x,a,a)')
     &    'Average vertical velocity through bottom of "T" and "u" '
     &,   'cells at each level','level', '   adv_vbt err '
     &, '  "T" cell area ', '    adv_vbu err ','  "U" cell area'
     &, '(Note: adv_vbu err only goes to zero when non-zero values on'
     &, ' boundary cells are taken into account)'
          do k=1,km
            if (tcella(k) .ne. c0) then
              wtlev(k,0) = wtlev(k,0)/tcella(k)
            else
              wtlev(k,0) = c0
            endif
            if (ucella(k) .ne. c0) then
              wulev(k,0) = wulev(k,0)/ucella(k)
            else
              wulev(k,0) = c0
            endif
            write (stdout,'(i4,4(2x,e14.7))')
     &      k, wtlev(k,0), tcella(k), wulev(k,0), ucella(k)
          enddo
        endif

        if (ioglen .ne. stdout .or. ioglen .lt. 0) then
          reltim = prelyr

          call getunit (io, 'energy_int.dta'
     &,                'unformatted sequential append ieee')

          write (stdout, *)
     &      ' ==> Global energy integrals written unformatted'
     &,     ' to file energy_int.dta on ts =', itt, stamp
          iotext = 'read(ioglen) reltim,(engint(i),engext(i),i=1,8)'
          iotext(46:) = ',plicin,plicex,buoy,buoerr,enleak'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, (engint(0,i,1),engext(i,1),i=1,8)
     &,                 plicin, plicex, buoy(0,1), buoerr, enleak
          call relunit (io)
        endif
      endif
9100  format(///,1x,
     &/1x,a,a,i10,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', ocean surface area =',e16.9,' cm**2'//' work by:',14x
     &,     'internal mode                         external mode'/)
9101  format(a21,2(1pe15.6, ' (',z16,' hex)'))
9110  format(/' work by buoyancy forces   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' energy conversion error   =',1pe14.6, ' (',z16,' hex)'/
     &,       ' nonlinear error           =',1pe14.6, ' (',z16,' hex)'/)
9111  format(/' max "t" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max "u" cell continuity error =',1pe14.6, ' at location'
     &,       ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')')
9112  format(/' max bottom "adv_vbt" (error)       =',1pe14.6
     &, ' at location', ' (i,jrow,k) = ','(',i4,',',i4,',',i4,')'
     &,      /' max bottom "adv_vbu" (slope vel)   =',1pe14.6
     &, ' at location',' (i,jrow,k) = ','(',i4,',',i4,',',i4,')')
# endif
# if defined O_term_balances

!-----------------------------------------------------------------------
!     add the external mode part of d/dt into the momentum balance,
!     the external mode part of the implicit coriolis term, and the
!     surface pressure gradients into the specified volumes
!-----------------------------------------------------------------------

      if (trmbts) then

        call utb3

!-----------------------------------------------------------------------
!     integrate previously computed term balance components vertically
!-----------------------------------------------------------------------

        do n=0,numreg
          if (n .gt. 0) then
            nv = (n-1)/nhreg + 1
            ks = llvreg(nv,1)
            ke = llvreg(nv,2)
            do ll=1,17
              do k=ke,ks,-1
                termbm(0,ll,1,n) = termbm(0,ll,1,n) + termbm(k,ll,1,n)
                termbm(0,ll,2,n) = termbm(0,ll,2,n) + termbm(k,ll,2,n)
              enddo
            enddo
          else
            ks = 1
            ke = km
          endif
          do m=1,nt
            do ll=1,15
              do k=ke,ks,-1
                termbt(0,ll,m,n) = termbt(0,ll,m,n) + termbt(k,ll,m,n)
              enddo
            enddo

!           construct change due to convection and filtering

            dtconv = termbt(0,10,m,n) - termbt(0,9,m,n)
            dtfilt = termbt(0,1,m,n) - termbt(0,10,m,n)
            termbt(0,9,m,n) = dtconv
            termbt(0,10,m,n) = dtfilt
          enddo
        enddo

!       normalize integrals by appropriate volume (or area)

        do n=0,numreg
          do m=1,nt
            if (n .le. nhreg)  then
              stflx(m,n) = stflx(m,n)*rareat(n)
              asst(m,n)  = asst(m,n)*rareat(n)
            endif
            do ll=1,15
              termbt(0,ll,m,n) = termbt(0,ll,m,n)*rvolt(n)
            enddo
          enddo
        enddo

        do n=0,numreg
          if (n .le. nhreg) then
            smflx(1,0) = smflx(1,0) + smflx(1,n)
            smflx(2,0) = smflx(1,0) + smflx(2,n)
            smflx(1,n) = smflx(1,n)*rareau(n)
            smflx(2,n) = smflx(2,n)*rareau(n)
          endif
          if (n .gt. 0) then
            avgw(n) = avgw(n)*rvolu(n)
            do ll=1,17
              termbm(0,ll,1,n) = termbm(0,ll,1,n)*rvolu(n)
              termbm(0,ll,2,n) = termbm(0,ll,2,n)*rvolu(n)
            enddo
          endif
        enddo
        smflx(1,0) = smflx(1,0)*rareau(0)
        smflx(2,0) = smflx(2,0)*rareau(0)

        if (iotrmb .eq. stdout .or. iotrmb .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'T E R M    B A L A N C E S'

          n = 0
          taux = smflx(1,n)
          tauy = smflx(2,n)
          write (stdout,10106)
     &                   'All regions added together for ts ='
     &,                     itt, stamp, volt(n), areat(n)
     &,                     volu(n), areau(n)
          write (stdout,10104)
          write (stdout,10098) n, ' smf(1)   = ', taux,' dynes/cm**2   '
          write (stdout,10098) n, ' smf(2)   = ', tauy,' dynes/cm**2   '
          do m=1,nt
            write (stdout,10098) n, ustf(m,1), stflx(m,n), ustf(m,2)
          enddo
          write (stdout,10098) n, ' tot heat = '
     &,                           (termbt(0,15,1,n)*volt(n))
     &,                                        ' deg C * cm**3 '
          write (stdout,10098) n, ' sst      = ',asst(1,n)
     &,                                       ' deg C         '

          do n=0,numreg
            if (n .eq. 1) then
              write (stdout,10050)
     &  'Regional averaged Momentum & Tracer Term Balances for ts =   '
     &,       itt, stamp
            endif
            iv = 0
            if (n .gt. 0) then
              iv = (n-1)/nhreg + 1
              ih = n - (iv-1)*nhreg
              write (stdout,10100)
     &        'Momentum terms averaged over region #'
     &,          n, ': ', hregnm(ih), vregnm(iv), volu(n), areau(n)
              write (stdout,10104)
              write (stdout,10101) n,' dU/dt   = ', termbm(0,1,1,n)
     &,                              ' dV/dt   = ', termbm(0,1,2,n)
              write (stdout,10101) n,' -Px     = ', termbm(0,2,1,n)
     &,                              ' -Py     = ', termbm(0,2,2,n)
              write (stdout,10101) n,' -surf Px= ', termbm(0,12,1,n)
     &,                              ' -surf Py= ', termbm(0,12,2,n)
              write (stdout,10101) n,' -(UU)x  = ', termbm(0,3,1,n)
     &,                              ' -(UV)x  = ', termbm(0,3,2,n)
              write (stdout,10101) n,' -(VU)y  = ', termbm(0,4,1,n)
     &,                              ' -(VV)y  = ', termbm(0,4,2,n)
              write (stdout,10101) n,' -(WU)z  = ', termbm(0,5,1,n)
     &,                              ' -(WV)z  = ', termbm(0,5,2,n)
              write (stdout,10101) n,'ADV_Umet = ', termbm(0,13,1,n)
     &,                              '-ADV_Vmet= ', termbm(0,13,2,n)
              write (stdout,10101) n,'  DIFF_Ux= ', termbm(0,6,1,n)
     &,                              '  DIFF_Vx= ', termbm(0,6,2,n)
              write (stdout,10101) n,'  DIFF_Uy= ', termbm(0,7,1,n)
     &,                              '  DIFF_Vy= ', termbm(0,7,2,n)
              write (stdout,10101) n,'  DIFF_Uz= ', termbm(0,8,1,n)
     &,                              '  DIFF_Vz= ', termbm(0,8,2,n)
              write (stdout,10101) n,'DIFF_Umet= ', termbm(0,9,1,n)
     &,                              'DIFF_Vmet= ', termbm(0,9,2,n)
              write (stdout,10101) n,'  fV     = ', termbm(0,10,1,n)
     &,                              ' -fU     = ', termbm(0,10,2,n)
              write (stdout,10101) n,'  source = ', termbm(0,11,1,n)
     &,                              '  source = ', termbm(0,11,2,n)
              erru = c0
              errv = c0
              do lll=2,13
                erru = erru + termbm(0,lll,1,n)
                errv = errv + termbm(0,lll,2,n)
              enddo
              write (stdout,10101) n,'  error  = ', termbm(0,1,1,n)-erru
     &,                              '  error  = ', termbm(0,1,2,n)-errv

              write (stdout,*) ' '
              write (stdout,10101) n,' -U(U)x  = ', termbm(0,14,1,n)
     &,                              ' -U(V)x  = ', termbm(0,14,2,n)
              write (stdout,10101) n,' -V(U)y  = ', termbm(0,15,1,n)
     &,                              ' -V(V)y  = ', termbm(0,15,2,n)
              write (stdout,10101) n,' -W(U)z  = ', termbm(0,16,1,n)
     &,                              ' -W(V)z  = ', termbm(0,16,2,n)

!             mass conservation within volume: Ux + Vy + Wz

              contu = (termbm(0,3,1,n) + termbm(0,4,1,n) +
     &                termbm(0,5,1,n)) - (termbm(0,14,1,n) +
     &                termbm(0,15,1,n) + termbm(0,16,1,n))
              write (stdout,10101) n,' mass err= ', contu

              write (stdout,10101) n,'  ubar   = ', termbm(0,17,1,n)
     &,                              '  vbar   = ', termbm(0,17,2,n)
     &,                              '  wbar   = ', avgw(n)
            endif
            if (iv .eq. 1) then
              write (stdout,10101) n,'  surf Uz= ', smflx(1,n)
     &,                              '  surf Vz= ', smflx(2,n)
            endif

            if (n .eq. 0) then
              write (stdout,10051)
     &   'Global averaged (all basins) Tracer Term Balances for ts = '
     &,       itt, stamp
            else
              write (stdout,10100)
     &        'Tracer terms averaged over region   #'
     &,        n, ': ', hregnm(ih), vregnm(iv), volt(n), areat(n)
            endif

            do m=1,nt
              dchg            = termbt(0,2,m,n) + termbt(0,3,m,n) +
     &                          termbt(0,4,m,n) + termbt(0,5,m,n) +
     &                          termbt(0,6,m,n) + termbt(0,7,m,n) +
     &                          termbt(0,8,m,n) + termbt(0,9,m,n) +
     &                          termbt(0,10,m,n)
              terr(m)        = termbt(0,1,m,n) - dchg
            enddo
            maxm = (nt-1)/7 + 1
            do mloop=1,maxm
              ms = (mloop-1)*7 + 1
              me = min(ms + 6,nt)
              write (stdout,10103) (trname(m),m=ms,me)
              write (stdout,10102) n,' dT/dt   = ', (termbt(0,1,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(UT)x  = ', (termbt(0,2,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(VT)y  = ', (termbt(0,3,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -(WT)z  = ', (termbt(0,4,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tx= ', (termbt(0,5,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Ty= ', (termbt(0,6,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  DIFF_Tz= ', (termbt(0,7,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  source = ', (termbt(0,8,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  convct = ', (termbt(0,9,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  filter = ', (termbt(0,10,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  error  = ', (terr(m)
     &,                                           m=ms,me)

              write (stdout,*) ' '
              write (stdout,10102) n,' -U(T)x  = ', (termbt(0,11,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -V(T)y  = ', (termbt(0,12,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,' -W(T)z  = ', (termbt(0,13,m,n)
     &,                                           m=ms,me)

!             mass conservation within volume: Ux + Vy + Wz

              cont = (termbt(0,2,1,n) + termbt(0,3,1,n) +
     &               termbt(0,4,1,n)) - (termbt(0,11,1,n) +
     &               termbt(0,12,1,n) + termbt(0,13,1,n))
              write (stdout,10102) n,' mass err= ', cont

              write (stdout,10102) n,'  chg var= ', (termbt(0,14,m,n)
     &,                                           m=ms,me)
              write (stdout,10102) n,'  tbar   = ', (termbt(0,15,m,n)
     &,                                           m=ms,me)
              if (iv .eq. 1) then
                write (stdout,10102) n,'  surflx = ', (stflx(m,n)
     &,                                           m=ms,me)
                taux = smflx(1,n)
                tauy = smflx(2,n)
                write (stdout,10105) ' Regionally averaged quantities:'
                write (stdout,'(1x)')
                write (stdout,10098) n
     &,            ' smf(1)   = ', taux,' dynes/cm**2   '
                write (stdout,10098) n
     &,            ' smf(2)   = ', tauy,' dynes/cm**2   '
                do m=1,nt
                  write (stdout,10098) n
     &,            ustf(m,1), stflx(m,n), ustf(m,2)
                enddo
                if (ms .eq. 1) then
                  write (stdout,10098) n, ' tot heat = '
     &,             (termbt(0,15,1,n)*volt(n)),' deg C * cm**3 '
                  write (stdout,10098) n
     &,            ' sst      = ', asst(1,n),' deg C         '
                endif
              endif
            enddo
          enddo
        endif

        if (iotrmb .ne. stdout .or. iotrmb .lt. 0) then

          reltim  = prelyr
          tperiod = 0.0
          call getunit (io, 'term_bal.dta'
     &,                'unformatted sequential append ieee')

          write (stdout, *)
     &     ' ==> Term balances written unformatted to file term_bal.dta'
     &,    ' for ts =', itt, ' ',stamp

          iotext =' read (iotrmb) reltim, nt, numreg, nhreg, km'
          iotext(45:) =', ((ustf*15(n,i),n=1,nt),i=1,2)'
          write (io) pstamp, iotext, expnam
          write (io)  reltim, nt, numreg, nhreg, km, ustf

          iotext ='read (iotrmb) (hregnm*40(n),n=1,nhreg)'
          iotext(39:)=', (vregnm*20(n),n=1,nvreg)'
          write (io) pstamp, iotext, expnam
          write (io)  hregnm, vregnm

          iotext ='read (iotrmb) (trname*12(n),n=1,nt)'
          write (io) pstamp, iotext, expnam
          write (io)  trname

          iotext ='read (iotrmb) (volu(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volu, numreg+1)

          iotext ='read (iotrmb) (volt(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, volt, numreg+1)

          iotext ='read (iotrmb) (areau(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areau, numreg+1)

          iotext ='read (iotrmb) (areat(l),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, areat, numreg+1)

          iotext ='read (iotrmb) ((smflx(i,l),i=1,2),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, smflx, 2*(nhreg+1))

          iotext ='read (iotrmb) ((stflx(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, stflx, nt*(nhreg+1))

          iotext ='read (iotrmb) ((asst(n,l),n=1,nt),l=0,nhreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, asst, nt*(nhreg+1))

          iotext ='read (iotrmb) (avgw(l),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, avgw, numreg)

          iotext = 'read (iotrmb) (((termbm'
          iotext(24:) = '(i,n,l),i=1,17),n=1,2),l=1,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbm(0,i,n,l),i=1,17),n=1,2),l=1,numreg)

          iotext = 'read (iotrmb) (((termbt'
          iotext(24:) = '(i,n,l),i=1,15),n=1,nt),l=0,numreg)'
          write (io) pstamp, iotext, expnam
          write (io) (((termbt(0,i,n,l),i=1,15),n=1,nt),l=0,numreg)

          call relunit (io)
        endif
      endif
10050 format(///,1x,/1x,a,i10,a/)
10051 format(///,1x,a,i10,a/)
10098 format(1x,'(',i4,')',a12,(1pe16.8,a15))
10100 format(/1x,a,i5,a,a,a/1x,'ocean volume =',e16.9,' cm**3'
     &, ', horizontal ocean area =',e16.9,' cm**2'/)
10101 format(1x,'(',i4,')',a11,1pe15.7, 2x, a11,1pe15.7, 2x, a11
     &,      1pe15.7)
10102 format(1x,'(',i4,')',a11,7(1pe15.7))
10103 format(' Region',11x,8a15)
10104 format(' Region')
10105 format(/8x,a32)
10106 format(///,1x,
     &/1x,a,i10,a/1x,'ocean "t" volume =',e16.9,' cm**3 '
     &, ', ocean "t" surface area =',e16.9,' cm**2 '
     &/ 1x,'ocean "u" volume =',e16.9,' cm**3, '
     &, ' ocean "u" surface area =',e16.9,' cm**2'/)
# endif
# if defined O_gyre_components

!-----------------------------------------------------------------------
!     write out the gyre_components
!-----------------------------------------------------------------------

!     convert heat transport to petawatts,
!     salt transport to 10**10 cm**3/sec

      if (gyrets) then
        pwatts = 4.186e-15
        csalt  = 1.e-10
        do jrow=1,jmt
          do ll=1,8
            ttn(ll,jrow,1)=ttn(ll,jrow,1)*pwatts
            ttn(ll,jrow,2)=ttn(ll,jrow,2)*csalt
          enddo
        enddo
        if (iogyre .eq. stdout .or. iogyre .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'G Y R E   C O M P O N E N T S'

          write (stdout,8195)
          do jrow=2,jmtm2
            l = jmt - jrow
            write (stdout,8196) l, (ttn(i,l,1),i=1,8)
     &,                            (ttn(i,l,2),i=1,8)
          enddo

          do m=0,nhreg
            do jrow=1,jmt
              do ll=6,8
                ttn2(ll,jrow,1,m) = ttn2(ll,jrow,1,m)*4.186e-15
                ttn2(ll,jrow,2,m) = ttn2(ll,jrow,2,m)*1.e-10
              enddo
            enddo
          enddo

          write (stdout,82001)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
              write (stdout,8204) jrow,(ttn2(i,jrow,1,m),i=6,8)
            enddo
          enddo

          write (stdout,8205)
          do m=0,nhreg
            if (m .ne. 0) then
              write (stdout,8201) hregnm(m)
            else
              write (stdout,8202)
            endif
            write (stdout,8203)
            do jrow=1,jmt
              write (stdout,8204) jrow,(ttn2(i,jrow,2,m),i=6,8)
            enddo
          enddo
        endif
        if (iogyre .ne. stdout .or. iogyre .lt. 0) then
          reltim = prelyr

          call getunit (io, 'gyre_comp.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Gyre components written'
     &,      ' unformatted to file gyre_comp.dta on ts=',itt, ' ',stamp

          iotext = 'read (iogyre) jmt, ntmin2, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, ntmin2, reltim

          iotext =
     &      'read (iogyre) (((ttn(l,j,n),l=1,8),j=1,jmt),n=1,ntmin2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, ttn, 8*jmt*ntmin2)

          iotext =
     &    'read (iogyre) (((ttn2(6..8,1..jmt,1..2,0..nhreg)'
          write (io) pstamp, iotext, expnam
          len = 3*jmt*nt*(1+nhreg)
          call wrufio (io, ttn2, len)

          call relunit (io)
        endif
      endif
8195  format(//,30x,'Gyre Components'/
     &,6x,'northward transport of heat (x10**15 watts)'
     &       ,22x,'northward transport of salt (x10**10 cm**3/sec)',/,
     &       3x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     &       'diffus   total'))
8196  format(i4,8f8.3,1x,8f8.3)
82001 format (/,6x,'northward transport of heat (x10**15 watts)',/)
8201  format (/,22x,a40,/)
8202  format (/,22x,' Global ',/)
8203  format (8x,'total advection',2x,'total diffusion',2x,'total')
8204  format (2x,i3,3x,3(e12.5,5x))
8205  format (/,6x,'northward transport of salt (x10**10 cm**3/sec)',/)
# endif
# if defined O_meridional_overturning

!-----------------------------------------------------------------------
!     write out the meridional overturning (mass)
!-----------------------------------------------------------------------

      if (vmsfts) then
        if (iovmsf .eq. stdout .or. iovmsf .lt. 0) then

          write (stdout,'(///,40x,a,/)')
     &    'M E R I D I O N A L    O V E R T U R N I N G'

          scl = 1.e12
          write (stdout,8194)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
          call matrix (vmsf, jmt, js, je, 1, km, scl)
        endif
        if (iovmsf .ne. stdout .or. iovmsf .lt. 0) then
          reltim = prelyr

          call getunit (io, 'overturn.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Meridional overturning of mass'
     &,    ' written unformatted to file overturn.dta on ts=',itt
     &,    ' ',stamp

          iotext = 'read (iovmsf) jmt, km, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, km, reltim

          iotext = 'read (iovmsf)  (zw(k),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zw, km)

          iotext = 'read (iovmsf)  (yu(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yu, jmt)

          iotext = 'read (iovmsf)  ((vmsf(j,k),j=1,jmt),k=1,km)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, vmsf, jmt*km)

          call relunit (io)
        endif
      endif
8194  format(/,' meridional overturning stream function (sverdrups)')
# endif
# if defined O_show_external_mode
#  if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface

!-----------------------------------------------------------------------
!     write out the external mode velocities and surface pressure
!-----------------------------------------------------------------------

      if (extts .and. eots) then
        if (ioext .eq. stdout .or. ioext .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'E X T E R N A L    M O D E'

          is = indp (slonxy, xt, imt)
          ie = indp (elonxy, xt, imt)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)
          scl=1.0

          write (stdout,8000) ' ubar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ubar(1,1,1), imt, is, ie, -js, -je, scl)
          write (stdout,8000) ' vbar '
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ubar(1,1,2), imt, is, ie, -js, -je, scl)
          scl = 0.0
          write (stdout,8000) ' div ps'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (divf(1,1), imt, is, ie, -js, -je, scl)
          scl = grav
!          write (stdout,8000) ' surface press(gm/cm/sec**2)'
          write (stdout,8000)
     &   'sea surface height (surface_press/grav) in cm'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (ps(1,1,1), imt, is, ie, -js, -je, scl)
        endif
        if (ioext .ne. stdout .or. ioext .lt. 0) then
          reltim = relyr

#   if defined O_rigid_lid_surface_pressure
          call getunit (io, 'surf_press.dta'
     &,                'unformatted sequential append ieee')
#   endif
#   if defined O_implicit_free_surface
          call getunit (io, 'ifree_surf.dta'
     &,                'unformatted sequential append ieee')
#   endif

          iotext = 'read (ioext) imt, jmt, reltim'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, reltim

          iotext = 'read (ioext)  (xt(i),i=1,imt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, xt, imt)

          iotext = 'read (ioext)  (yt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, yt, jmt)

          iotext = 'read (ioext)  ((ps(i,j),i=1,imt),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, ps(1,1,1), imt*jmt)
          write (stdout,*) ' => Surface pressure written unformatted'
     &,              ' to unit ',io, ' on ts=',itt, stamp
#   if defined O_rigid_lid_surface_pressure
     &,              ' to file surf_press.dta on ts=',itt, stamp
#   endif
#   if defined O_implicit_free_surface
     &,              ' to file ifree_surf.dta on ts=',itt, stamp
#   endif

          call relunit (io)
        endif
      endif
#  endif
#  if defined O_stream_function

!-----------------------------------------------------------------------
!     write out the stream function
!-----------------------------------------------------------------------

      if (extts .and. eots) then
        if (ioext .eq. stdout .or. ioext .lt. 0) then

          write (stdout,'(///,40x,a,/)') 'E X T E R N A L    M O D E'

          is = indp (slonxy, xt, imt)
          ie = indp (elonxy, xt, imt)
          js = indp (slatxy, yt, jmt)
          je = indp (elatxy, yt, jmt)

          scl=1.e12
          write (stdout,8000) ' stream function (sverdrups)'
     &,   itt, xt(is), xt(ie), yt(js), yt(je), scl
          call matrix (psi(1,1,1), imt, is, ie, -js, -je, scl)
        endif
        if (ioext .ne. stdout .or. ioext .lt. 0) then
          reltim = relyr

          call getunit (io, 'psi.dta'
     &,                'unformatted sequential append ieee')

          iotext = 'read (ioext) imt, jmt, reltim'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, reltim

          iotext = 'read (ioext)  (xt(i),i=1,imt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, xt, imt)

          iotext = 'read (ioext)  (yt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, yt, jmt)

          iotext = 'read (ioext)  ((psi(i,j),i=1,imt),j=1,jmt)'
          write (io) stamp, iotext, expnam
          call wrufio (io, psi(1,1,1), imt*jmt)
          write (stdout,*) ' => Stream function written unformatted'
     &,              ' to file psi.dta on ts=',itt, stamp

          call relunit (io)
        endif
      endif
#  endif
8000  format(1x,a,1x,'ts=',i7
     &,', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     &,', scaling=',1pg10.3)
# endif
# if defined O_meridional_tracer_budget

!-----------------------------------------------------------------------
!     construct and write out the averaged components for the
!     meridional tracer balance
!-----------------------------------------------------------------------

      if (tmbts) then

!       quantities have been weighted by volume at each cell and
!       summed over all ocean longitudes and depths for each latitude.
!       average quantities over "numtmb" days and do units conversion
!       to petawatts for temperature and giga cm**3/sec for salt

        tmbavg = numtmb * dtts
        rnum = c1 / float(numtmb)

        do mask=1,ntmbb
          do jrow=1,jmt
            smdvol(jrow,mask) = rnum*smdvol(jrow,mask)
          enddo
        enddo

        do mask=1,ntmbb
          do n=1,nt
            if (n .eq. 1) then
              convrt = rnum * 4.186e-15
            elseif (n .eq. 2) then
              convrt = rnum * 1.e-9
            else
              convrt = c1
            endif
            do jrow=1,jmt
              tstor(jrow,n,mask) = convrt*tstor(jrow,n,mask)
              tdiv(jrow,n,mask)  = convrt*tdiv(jrow,n,mask)
              tflux(jrow,n,mask) = convrt*tflux(jrow,n,mask)
              tdif(jrow,n,mask)  = convrt*tdif(jrow,n,mask)
              tsorc(jrow,n,mask) = convrt*tsorc(jrow,n,mask)
            enddo
          enddo
        enddo

!       accumulate terms for all basins and store in basin # 0

        do jrow=1,jmt
          smdvol(jrow,0) = c0
          do mask=1,ntmbb
            smdvol(jrow,0) = smdvol(jrow,0) + smdvol(jrow,mask)
          enddo
        enddo

        do n=1,nt
          do jrow=1,jmt
            tstor(jrow,n,0) = c0
            tdiv(jrow,n,0)  = c0
            tflux(jrow,n,0) = c0
            tdif(jrow,n,0)  = c0
            tsorc(jrow,n,0) = c0
          enddo
        enddo

        do mask=1,ntmbb
          do n=1,nt
            do jrow=1,jmt
              tstor(jrow,n,0) = tstor(jrow,n,0) + tstor(jrow,n,mask)
              tdiv(jrow,n,0)  = tdiv(jrow,n,0)  + tdiv(jrow,n,mask)
              tflux(jrow,n,0) = tflux(jrow,n,0) + tflux(jrow,n,mask)
              tdif(jrow,n,0)  = tdif(jrow,n,0)  + tdif(jrow,n,mask)
              tsorc(jrow,n,0) = tsorc(jrow,n,0) + tsorc(jrow,n,mask)
            enddo
          enddo
        enddo

!       write out the results

        if (iotmb .eq. stdout .or. iotmb .lt. 0) then

          write (stdout,'(///,20x,a,/)')
     &     'M E R I D I O N A L    T R A C E R    B U D G E T'

          do n=1,nt

!           construct sum in latitude for all basins

            sumsto = c0
            sumdiv = c0
            sumflx = c0
            sumdif = c0
            sumsor = c0
            sumvol = c0
            do jrow=1,jmt
              sumsto = sumsto + tstor(jrow,n,0)
              sumdiv = sumdiv + tdiv(jrow,n,0)
              sumflx = sumflx + tflux(jrow,n,0)
              sumdif = sumdif + tdif(jrow,n,0)
              sumsor = sumsor + tsorc(jrow,n,0)
              sumvol = sumvol + smdvol(jrow,0)
            enddo

            write (stdout,*) ' '
            if (n .eq. 1) then
              write (stdout,*) '(Tracer # 1 results are in petawatts)'
            elseif (n .eq. 2) then
              write (stdout,*)
     &       '(Tracer # 2 results are in giga cm**3/sec)'
            else
              write (stdout,*) ' '
            endif
            write (stdout,8050) n, tmbavg*secday, stamp

            do jrow=jmt,1,-1
              terror = tstor(jrow,n,0) - (tdiv(jrow,n,0) +
     &               tflux(jrow,n,0) + tdif(jrow,n,0) + tsorc(jrow,n,0))
              write (stdout,8100) jrow, yt(jrow), tstor(jrow,n,0)
     &,       tdiv(jrow,n,0), tflux(jrow,n,0), tdif(jrow,n,0)
     &,       tsorc(jrow,n,0), terror, smdvol(jrow,0)
            enddo
            serror = sumsto - (sumdiv + sumflx + sumdif + sumsor)
            write (stdout,8200)  sumsto, sumdiv, sumflx, sumdif
     &,                          sumsor, serror, sumvol
          enddo
        endif
        if (iotmb .ne. stdout .or. iotmb .lt. 0) then
          reltim = relyr

          call getunit (io, 'tracer_bud.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Meridional Tracer Budget written '
     &,   ' unformatted to file tracer_bud.dta on ts=',itt, ' ',stamp

          iotext = 'read (iotmb) imt, jmt, nt, ntmbb, reltim, avgper'
          write (io) stamp, iotext, expnam
          write (io) imt, jmt, nt, ntmbb, reltim, tmbavg*secday

          iotext = 'read (iotmb) (yt(j),j=1,jmt), (dyt(j),j=1,jmt)'
          write (io) stamp, iotext, expnam
          write (io) (yt(j),j=1,jmt), (dyt(j),j=1,jmt)

          iotext =
     &   'read (iotmb) (((tstor(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tstor, jmt*nt*(ntmbb+1))

          iotext =
     &    'read (iotmb) (((tdiv(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tdiv, jmt*nt*(ntmbb+1))

          iotext =
     &    'read (iotmb) (((tflux(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tflux, jmt*nt*(ntmbb+1))

          iotext =
     &    'read (iotmb) (((tdif(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tdif, jmt*nt*(ntmbb+1))

          iotext =
     &    'read (iotmb) (((tsorc(j,n,l),j=1,jmt),n=1,nt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, tsorc, jmt*nt*(ntmbb+1))

          iotext =
     &    'read (iotmb) ((smdvol(j,l),j=1,jmt),l=0,ntmbb)'
          write (io) stamp, iotext, expnam
          call wrufio (io, smdvol, jmt*(ntmbb+1))

          call relunit (io)
        endif

!       zero quantities in preparation for next average

        numtmb = 0
        do mask=0,ntmbb
          do jrow=1,jmt
            smdvol(jrow,mask)  = c0
          enddo
        enddo
        do mask=0,ntmbb
          do m=1,nt
            do jrow=1,jmt
              tstor(jrow,m,mask) = c0
              tdiv(jrow,m,mask)  = c0
              tflux(jrow,m,mask) = c0
              tdif(jrow,m,mask)  = c0
              tsorc(jrow,m,mask) = c0
            enddo
          enddo
        enddo
      endif
8050  format (/1x,'Meridional Tracer Balance for tracer #',i2
     &,' averaged over depth and longitude',/1x
     &,' and averaged over a period of ', f10.2
     &,' days ending on ', a32//
     &,' jrow    lat      (T)t      (VT)y   DIFF_Tz     DIFF_Ty'
     &,'      source       error    volume'/)
8100  format (1x,i4,2x,f6.2,1x,7(1pe10.3,1x))
8200  format (1x,'global sum =',1x,7(1pe10.3,1x))
# endif
# if defined O_show_zonal_mean_of_sbc

!-----------------------------------------------------------------------
!     print the zonal mean boundary conditions and related items
!-----------------------------------------------------------------------

      if (zmbcts) then
        cmmday = 10.0*daylen
        cwatts = 4.186e4
        do jrow=1,jmt
          if (zmau(jrow) .ne. c0) then
            zmsmf(jrow,1) = zmsmf(jrow,1) / zmau(jrow)
            zmsmf(jrow,2) = zmsmf(jrow,2) / zmau(jrow)
            zmsm(jrow,1)  = zmsm(jrow,1)  / zmau(jrow)
            zmsm(jrow,2)  = zmsm(jrow,2)  / zmau(jrow)
          endif
          if (zmat(jrow) .ne. c0) then
            zmstf(jrow,1) = cwatts*zmstf(jrow,1)
            zmstf(jrow,2) = cmmday*zmstf(jrow,2)
            zmst(jrow,2)  = 1.0e3*zmst(jrow,2)
            do m=1,nt
              zmstf(jrow,m) = zmstf(jrow,m) / zmat(jrow)
              zmst(jrow,m)  = zmst(jrow,m)  / zmat(jrow)
            enddo
          endif
        enddo

        if (iozmbc .eq. stdout .or. iozmbc .lt. 0) then

          write (stdout,'(///,27x,a,/)')
     &  'Z O N A L   M E A N   O F    S U R F A C E   B.  C.'

          write (stdout,10200)
          do jj=1,jmt
            jrow = jmt + 1 - jj
            write (stdout,10201) jrow, yt(jrow), zmsmf(jrow,1)
     &,     zmsmf(jrow,2), zmstf(jrow,1), zmstf(jrow,2)
     &,     zmsm(jrow,1), zmsm(jrow,2), zmst(jrow,1), zmst(jrow,2)
          enddo
        endif
        if (iozmbc .ne. stdout .or. iozmbc .lt. 0) then
          reltim = prelyr
          call getunit (io, 'zmean_sbc.dta'
     &,                'unformatted sequential append ieee')

          write (stdout,*) ' => Zonal mean surface bc written'
     &,   ' unformatted to file zmean_sbc.dta on ts=',itt, ' ',stamp

          iotext = 'read (iozmbc) jmt, nt, reltim'
          write (io) pstamp, iotext, expnam
          write (io) jmt, nt, reltim

          iotext = 'read (iozmbc) (yt(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yt, jmt)

          iotext = 'read (iozmbc) (yu(j),j=1,jmt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, yu, jmt)

          iotext = 'read (iozmbc) ((zmsmf(j,n),j=1,jmt),n=1,2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmsmf, jmt*2)

          iotext = 'read (iozmbc) ((zmstf(j,n),j=1,jmt),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmstf, jmt*nt)

          iotext = 'read (iozmbc) ((zmsm(j,n),j=1,jmt),n=1,2)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmsm, jmt*2)

          iotext = 'read (iozmbc) ((zmst(j,n),j=1,jmt),n=1,nt)'
          write (io) pstamp, iotext, expnam
          call wrufio (io, zmst, jmt*nt)

          call relunit (io)
        endif

      endif
10200 format(//
     &,26x,'Zonal Mean of Surface Boundary Condition diagnostic'
     &, //,1x,' row',1x,'  lat ', '   Taux    ','   Tauy    '
     &,' Heat Flux ',' Prec-Evap ',' Surface u ',' Surface v '
     &,'   SST     ','   SSS    ',/,12x
     &,' dyn/cm**2 ',' dyn/cm**2 ',' watts/m**2','   mm/day  '
     &,'   cm/sec  ','   cm/sec  ','   deg C   ','   ppt-35  ',/)
10201 format (1x, i4, 1x, f6.2, 8(1pe11.3))
# endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!     save the time mean averages on the "averaging" grid
!-----------------------------------------------------------------------

      if (timavgts .and. timavgint .gt. c0) then
        call avgout
      endif
# endif
# if defined O_xbts

!-----------------------------------------------------------------------
!     add in surface pressure gradients and missing external mode
!     parts of d/dt and implicit coriolis term for all XBTs
!-----------------------------------------------------------------------

      if (xbtperts) then
        call uxbt3
      endif

!-----------------------------------------------------------------------
!     save all time averaged XBT data at end of averaging period
!-----------------------------------------------------------------------

      if ((xbtts .or. eorun) .and. xbtint .gt. c0) then
        call xbto
      endif
# endif
# if defined O_mom_tbt

      if ((tbtts .or. eorun) .and. tbtint > c0 .and. ntbtts > 0) then
!-----------------------------------------------------------------------
!       construct change due to convection and filtering
!-----------------------------------------------------------------------

        do j=1,jmt
          do i=2,imtm1
            do k=1,kmt(i,j)
              do n=1,nt
                tconv = tbt(i,j,k,n,10) - tbt(i,j,k,n,9)
                tfilt = tbt(i,j,k,n,1) - tbt(i,j,k,n,10)
                tbt(i,j,k,n,9)  = tconv
                tbt(i,j,k,n,10) = tfilt
              enddo
            enddo
            enddo
        enddo

!-----------------------------------------------------------------------
!       average the data, write it out, then zero the accumulators
!-----------------------------------------------------------------------

        rnavg = c1/ntbtts
        tm(:,:,:) = 0.
        errt(:) = 0.
        do j=1,jmt
          do i=2,imtm1
            do k=1,kmt(i,j)
              tm(i,j,k) = 1.
              do n=1,nt
                do m=1,11
                  tbt(i,j,k,n,m) = rnavg*tbt(i,j,k,n,m)
                enddo
                do m=2,10
                  errt(n) = errt(n) + tbt(i,j,k,n,m)
                enddo
                errt(n) = errt(n) - tbt(i,j,k,n,1)
              enddo
            enddo
            do n=1,nt
              tbtsf(i,j,n) = rnavg*tbtsf(i,j,n)
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       save all time averaged data at end of the averaging period
!-----------------------------------------------------------------------

        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        if (ftbt .eq. ' ') then
          call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
          nyear = time
          call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
          ftbt = file_stamp ('mom_tbt',nstamp,'.nc')
          ftbt = new_file_name (ftbt)
          call inqdefined (ftbt, defined)
          if (defined) then
            call opennext (ftbt, relyr, ntrec, iou)
          else
            call opennew (ftbt, iou)
          endif

          call mom_tbt_def (ftbt, imt, jmt, km, nt, xt, yt
     &,                     calendar, expnam, runstamp, mapt)
        endif
# if !defined O_save_time_relyear0
!       make output time relative to year 1
        time = time - 1.
# endif
        avgper = tbtper*accel
        if (avgper .le. 1e-6) avgper = 0.
# if defined O_save_time_endper
        tmp = 0.
# elif defined O_save_time_startper
        tmp = 1.
# else
        tmp = 0.5
# endif
# if defined O_units_time_years
#  if defined O_calendar_360_day
        time = time - tmp*avgper/360.
#  elif defined O_calendar_gregorian
        time = time - tmp*avgper/365.25
#  else
        time = time - tmp*avgper/365.
#  endif
# else
#  if defined O_calendar_360_day
        time = time*360. - tmp*avgper
#  elif defined O_calendar_gregorian
        time = time*365.25 - tmp*avgper
#  else
        time = time*365. - tmp*avgper
#  endif
# endif
        is = 1
        ie = imt
        js = 1
        je = jmt
        call mom_tbt_out (ftbt, is, ie, js, je, imt, jmt, km
     &,                   nt, xt, yt, zt, xu, yu, zw, dxt, dyt, dzt
     &,                   dxu, dyu, dzw, avgper, time, nstamp, mapt
     &,                   tbt, tbtsf, tlat, tlon, ulat, ulon, kmt
     &,                   mskhr, tm)
        write (stdout,'(a,i4,a,a,a,i10,a,a)')
     &    '=> Ocn tracer term balances #', ntrec, ' written to '
     &,   trim(ftbt),' on ts = ',itt, ', ', stamp

!-----------------------------------------------------------------------
!       zero the accumulators
!-----------------------------------------------------------------------

        do j=1,jmt
          do i=1,imt
            do n=1,nt
              do k=1,kmt(i,j)
                do m=1,11
                  tbt(i,j,k,n,m) = 0.
                enddo
              enddo
              tbtsf(i,j,n) = 0.
            enddo
          enddo
        enddo
        ntbtts = 0

      endif
# endif
# if defined O_save_carbon_carbonate_chem

!-----------------------------------------------------------------------
!     accumulate sea surface time averages
!-----------------------------------------------------------------------

      if (timavgperts .and. eots) then
        nta_sscar = nta_sscar + 1
        do j=2,jmtm1
          do i=2,imtm1
            ta_sspH(i,j) = ta_sspH(i,j) + sspH(i,j)
            ta_ssCO3(i,j) = ta_ssCO3(i,j) + ssCO3(i,j)
            ta_ssOc(i,j) = ta_ssOc(i,j) + ssOc(i,j)
            ta_ssOa(i,j) = ta_ssOa(i,j) + ssOa(i,j)
            ta_sspCO2(i,j) = ta_sspCO2(i,j) + sspCO2(i,j)
          enddo
        enddo
      endif
# endif
#endif

      return
      end

#if defined O_mom && defined O_time_step_monitor
      subroutine ta_mom_tsi (m)

!=======================================================================
!     ocean data time integral averaging

!     input:
!       m = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer i, j, k, ib(10), ic(10), m, n

      real dens, tq, sq, drodt, drods, drhodt, drhods, ddensdtdt
      real ddensdtds, ddensdsds, rntatio, rnv_otsf, rnt_slh, tarea
      real cstdyt, area, s, tins, tinsit, dins, sum, d, tmp, tmp1, tmp2

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "coord.h"
      include "cregin.h"
      include "grdvar.h"
      include "diag.h"
      include "diaga.h"
      include "emode.h"
      include "levind.h"
      include "mw.h"
      include "iounit.h"
# if defined O_mobi
      include "mobi.h"
# endif
# if defined O_sed && !defined O_sed_uncoupled
      include "sed.h"
# endif
      include "state.h"
      include "dens.h"

# if defined O_tai_otsf
      real otsf(km)
# endif
      real dmsk(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

      if (m .eq. 0.) then

!       zero
        ntatio = 0
# if defined O_tai_otsf_from_averages
        nv_otsf = 0
        v_otsf(:,:) = 0.
# endif
# if defined O_tai_slh_from_averages
        nt_slh = 0
        t_slh(:,:,:,:) = 0.
# endif
        tai_otmax = 0.
        tai_otmin = 0.
        tai_slh = 0.
        tai_ek = 0.
        tai_t = 0.
        tai_s = 0.
        tai_tvar = 0.
        tai_svar = 0.
        tai_dt = 0.
        tai_ds = 0.
        tai_scan = 0.
        tai_hflx = 0.
        tai_sflx = 0.
        tai_dic = 0.
        tai_dicflx = 0.
        tai_dic13 = 0.
        tai_dic13flx = 0.
        tai_dicwflx = 0.
        tai_dic13wflx = 0.
        tai_alk = 0.
        tai_o2 = 0.
        tai_o2flx = 0.
        tai_po4 = 0.
        tai_dop = 0.
        tai_phyt = 0.
        tai_zoop = 0.
        tai_detr = 0.
        tai_no3 = 0.
        tai_dfe = 0.
        tai_ddfe = 0.
# if defined O_mobi_silicon
        tai_sil = 0.
        tai_silflx = 0.
        tai_opl = 0.
# endif        
        tai_don = 0.
        tai_diaz = 0.
        tai_din15 = 0.
        tai_don15 = 0.
        tai_phytn15 = 0.
        tai_diatn15 = 0.        
        tai_zoopn15 = 0.
        tai_detrn15 = 0.
        tai_diazn15 = 0.
        tai_doc13 = 0.
        tai_phytc13 = 0.
        tai_diatc13 = 0.
        tai_zoopc13 = 0.
        tai_detrc13 = 0.
        tai_diazc13 = 0.
        tai_c14 = 0.
        tai_dc14 = 0.
        tai_c14flx = 0.
        tai_sspH = 0.
        tai_ssCO3 = 0.
        tai_ssOc = 0.
        tai_ssOa = 0.
        tai_sspCO2 = 0.
        tai_c = 0.
        tai_caco3 = 0.
        tai_caco3c13 = 0.
        tai_d_B = 0.
# if defined O_mobi_silicon
        tai_diat = 0.
# endif                
      elseif (m .eq. 1) then

!       accumulate
        ntatio = ntatio + 1

# if defined O_tai_otsf && !defined O_tai_otsf_from_averages

        tmp1 = 0.
        tmp2 = 0.
        if (nv_otsf .gt. 0) then
          tmp1 = -1.e20
          tmp2 = 1.e20
          rnv_otsf = 1./float(nv_otsf)
          do j=jsot,jeot
            do k=1,km
              v_otsf(j,k) = v_otsf(j,k)*rnv_otsf*dzt(k)*csu(j)
            enddo
!           integrate up from the bottom
            otsf(km) = 0.
            do k=km-1,1,-1
              otsf(k) = otsf(k+1) - v_otsf(j,k)
            enddo
            do k=ksot,keot
              if (otsf(k) .gt. tmp1) tmp1 = otsf(k)
              if (otsf(k) .lt. tmp2) tmp2 = otsf(k)
            enddo
          enddo
        endif
        tai_otmax = tai_otmax + tmp1
        tai_otmin = tai_otmin + tmp2
        nv_otsf = 0
        v_otsf(:,:) = 0.

# endif
# if defined O_tai_slh && !defined O_tai_slh_from_averages

        tmp = 0.
        if (nt_slh .gt. 0) then
          t_slh(:,:,:,:) = t_slh(:,:,:,:)/float(nt_slh)
          tarea = 0.
          do j=2,jmtm1
            cstdyt = cst(j)*dyt(j)
            do i=2,imtm1
              if (kmt(i,j) .ge. 1) then
                area = dxt(i)*cstdyt
                tarea = tarea + area
                do k=1,kmt(i,j)
                  s = (t_slh(i,j,k,2)*1000.) + 35.
                  d = zt(k)*0.01
                  tins = tinsit(t_slh(i,j,k,1), s, d)
                  dins = dens(tins-to(k), t_slh(i,j,k,2)-so(k), k)
!                 volume integrate (rho - rho0)/rho0
!                 assume rho0 = 1g/cm3 so difference is just delta rho
!                 dslh = sum[((rho1-rho0)/rho0)*dvol]/area
!                      - sum[((rho2-rho0)/rho0)*dvol}/area
!                      = sum[((rho1-rho2)/rho0)*dvol]/area
                  tmp = tmp - (dins - d_slh(i,j,k))*area*dzt(k)
                enddo
              endif
            enddo
          enddo
          if (tarea .gt. 0) tmp = tmp/tarea
        endif
        tai_slh = tai_slh + tmp
        nt_slh = 0
        t_slh(:,:,:,:) = 0.

# endif
        do j=2,jmtm1
          do k=km,1,-1
            ektot(0,1) = ektot(0,1) + ektot(k,j)
          enddo
          do n=1,nt
            do k=km,1,-1
              tbar(0,n,1)   = tbar(0,n,1) + tbar(k,n,j)
              travar(0,n,1) = travar(0,n,1) + travar(k,n,j)
              dtabs(0,n,1)  = dtabs(0,n,1) + dtabs(k,n,j)
            enddo
          enddo
        enddo
        ektot(0,1) = ektot(0,1)/ucellv

        do n=1,nt
          tbar(0,n,1) = tbar(0,n,1)/tcellv
          travar(0,n,1) = travar(0,n,1)/tcellv - tbar(0,n,1)**2
          dtabs(0,n,1) = dtabs(0,n,1)/tcellv
        enddo

        tai_ek = tai_ek + ektot(0,1)
        tai_t = tai_t + tbar(0,itemp,1)
        tai_s = tai_s + tbar(0,isalt,1)
        tai_tvar = tai_tvar + travar(0,itemp,1)
        tai_svar = tai_svar + travar(0,isalt,1)
        tai_dt = tai_dt + dtabs(0,itemp,1)
        tai_ds = tai_ds + dtabs(0,isalt,1)
        tai_scan = tai_scan + mscan

        dmsk(:,:) = 1.
        where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.

        call areaavg (sbc(1,1,ihflx), dmsk, tmp)
        tai_hflx = tai_hflx + tmp
        call areaavg (sbc(1,1,isflx), dmsk, tmp)
        tai_sflx = tai_sflx + tmp
# if defined O_carbon
        tai_dic = tai_dic + tbar(0,idic,1)
        call areaavg (sbc(1,1,idicflx), dmsk, tmp)
        tai_dicflx = tai_dicflx + tmp
#  if defined O_carbon_13
        tai_dic13 = tai_dic13 + tbar(0,idic13,1)
        call areaavg (sbc(1,1,idic13flx), dmsk, tmp)
        tai_dic13flx = tai_dic13flx + tmp
#  endif
#  if defined O_sed && !defined O_sed_uncoupled
        tai_dicwflx = tai_dicwflx + dicwflx
#   if defined O_carbon_13
        tai_dic13wflx = tai_dic13wflx ! + dic13wflx
#   endif
#  endif
#  if defined O_carbon_14
        tai_c14 = tai_c14 + tbar(0,ic14,1)
        tai_dc14 = tai_dc14 + dc14bar/tcellv
        call areaavg (sbc(1,1,ic14flx), dmsk, tmp)
        tai_c14flx = tai_c14flx + tmp
#  endif
# endif
# if defined O_mobi_silicon
        tai_sil = tai_sil + tbar(0,isil,1)
        call areaavg (sbc(1,1,isilflx), dmsk, tmp)
        tai_silflx = tai_silflx + tmp
        tai_opl = tai_opl + tbar(0,iopl,1)
# endif
# if defined O_mobi_alk
        tai_alk = tai_alk + tbar(0,ialk,1)
# endif
# if defined O_mobi_o2
        tai_o2 = tai_o2 + tbar(0,io2,1)
        call areaavg (sbc(1,1,io2flx), dmsk, tmp)
        tai_o2flx = tai_o2flx + tmp
# endif
# if defined O_mobi
        tai_po4 = tai_po4 + tbar(0,ipo4,1)
        tai_phyt = tai_phyt + tbar(0,iphyt,1)
        tai_zoop = tai_zoop + tbar(0,izoop,1)
        tai_detr = tai_detr + tbar(0,idetr,1)
#  if defined O_kk_ballast
        tai_d_B = tai_d_B + tbar(0,idetr_B,1)
#  endif
#  if defined O_mobi_caco3
        tai_caco3 = tai_caco3 + tbar(0,icaco3,1)
#  endif
#  if defined O_mobi_silicon
        tai_diat = tai_diat + tbar(0,idiat,1)
#  endif        
#  if defined O_mobi_nitrogen
        tai_dop = tai_dop + tbar(0,idop,1)
        tai_no3 = tai_no3 + tbar(0,ino3,1)
        tai_don = tai_don + tbar(0,idon,1)
        tai_diaz = tai_diaz + tbar(0,idiaz,1)
#   if defined O_mobi_nitrogen_15
        tai_din15 = tai_din15 + tbar(0,idin15,1)
        tai_don15 = tai_don15 + tbar(0,idon15,1)
        tai_phytn15 = tai_phytn15 + tbar(0,iphytn15,1)
#    if defined O_mobi_silicon
        tai_diatn15 = tai_diatn15 + tbar(0,idiatn15,1)
#    endif        
        tai_zoopn15 = tai_zoopn15 + tbar(0,izoopn15,1)
        tai_detrn15 = tai_detrn15 + tbar(0,idetrn15,1)
        tai_diazn15 = tai_diazn15 + tbar(0,idiazn15,1)
#   endif
#  endif
#  if defined O_mobi_iron
        tai_dfe = tai_dfe + tbar(0,idfe,1)
        tai_ddfe = tai_ddfe + tbar(0,idetrfe,1)
#  endif
#  if defined O_carbon_13
        tai_doc13 = tai_doc13 + tbar(0,idoc13,1)
        tai_phytc13 = tai_phytc13 + tbar(0,iphytc13,1)
#   if defined O_mobi_silicon
        tai_diatc13 = tai_diatc13 + tbar(0,idiatc13,1)
#   endif        
#   if defined O_mobi_caco3
        tai_caco3c13 = tai_caco3c13 + tbar(0,icaco3c13,1)
#   endif        
        tai_zoopc13 = tai_zoopc13 + tbar(0,izoopc13,1)
        tai_detrc13 = tai_detrc13 + tbar(0,idetrc13,1)
#   if defined O_mobi_nitrogen
        tai_diazc13 = tai_diazc13 + tbar(0,idiazc13,1)
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
        tai_cfc11 = tai_cfc11 + tbar(0,icfc11,1)
        call areaavg (sbc(1,1,icfc11flx), dmsk, tmp)
        tai_cfc11flx = tai_cfc11flx + tmp
        tai_cfc12 = tai_cfc12 + tbar(0,icfc12,1)
        call areaavg (sbc(1,1,icfc12flx), dmsk, tmp)
        tai_cfc12flx = tai_cfc12flx + tmp
# endif
# if defined O_save_carbon_carbonate_chem
        dmsk(:,:) = 1.
        where (sspH(:,:).gt.11.5 .or. sspH(:,:).le.0.5) dmsk(:,:) = 0.
        call areaavg (sspH, dmsk, tmp)
        tai_sspH = tai_sspH + tmp
        call areaavg (ssOc, dmsk, tmp)
        tai_ssOc = tai_ssOc + tmp
        call areaavg (ssOa, dmsk, tmp)
        tai_ssOa = tai_ssOa + tmp
        call areaavg (ssCO3, dmsk, tmp)
        tai_ssCO3 = tai_ssCO3 + tmp
        call areaavg (sspCO2, dmsk, tmp)
        tai_sspCO2 = tai_sspCO2 + tmp
# endif

      elseif (m .eq. 2 .and. ntatio .ne. 0) then

!       average
        rntatio = 1./float(ntatio)
# if defined O_tai_otsf_from_averages

        if (nv_otsf .gt. 0) then
          rnv_otsf = 1./float(nv_otsf)
          tai_otmax = -1.e20
          tai_otmin = 1.e20
          do j=jsot,jeot
            do k=1,km
              v_otsf(j,k) = v_otsf(j,k)*rnv_otsf*dzt(k)*csu(j)
            enddo
!           integrate up from the bottom
            otsf(km) = 0.
            do k=km-1,1,-1
              otsf(k) = otsf(k+1) - v_otsf(j,k)
            enddo
            do k=ksot,keot
              if (otsf(k) .gt. tai_otmax) tai_otmax = otsf(k)
              if (otsf(k) .lt. tai_otmin) tai_otmin = otsf(k)
            enddo
          enddo
        endif

# else
        tai_otmax = tai_otmax*rntatio
        tai_otmin = tai_otmin*rntatio
# endif
# if defined O_tai_slh_from_averages

        if (nt_slh .gt. 0) then
          rnt_slh = 1./float(nt_slh)
          t_slh(:,:,:,:) = t_slh(:,:,:,:)*rnt_slh
          tai_slh = 0.0
          tarea = 0.0
          do j=2,jmtm1
            cstdyt = cst(j)*dyt(j)
            do i=2,imtm1
              if (kmt(i,j) .ge. 1) then
                area = dxt(i)*cstdyt
                tarea = tarea + area
                do k=1,kmt(i,j)
                  s = (t_slh(i,j,k,2)*1000.) + 35.
                  d = zt(k)*0.01
                  tins = tinsit(t_slh(i,j,k,1), s, d)
                  dins = dens(tins-to(k), t_slh(i,j,k,2)-so(k), k)
!                 volume integrate (rho - rho0)/rho0
!                 assume rho0 = 1g/cm3 so difference is just delta rho
!                 dslh = sum[((rho1-rho0)/rho0)*dvol]/area
!                      - sum[((rho2-rho0)/rho0)*dvol}/area
!                      = sum[((rho1-rho2)/rho0)*dvol]/area
                  tai_slh = tai_slh - (dins - d_slh(i,j,k))*area*dzt(k)
                enddo
              endif
            enddo
          enddo
          if (tarea .gt. 0) tai_slh = tai_slh/tarea
        endif

# else
        tai_slh = tai_slh*rntatio
# endif
        tai_ek = tai_ek*rntatio
        tai_t = tai_t*rntatio
        tai_s = tai_s*rntatio
        tai_tvar = tai_tvar*rntatio
        tai_svar = tai_svar*rntatio
        tai_dt = tai_dt*rntatio
        tai_ds = tai_ds*rntatio
        tai_scan = tai_scan*rntatio
        tai_hflx = tai_hflx*rntatio
        tai_sflx = tai_sflx*rntatio
        tai_dic = tai_dic*rntatio
        tai_dicflx = tai_dicflx*rntatio
        tai_dic13 = tai_dic13*rntatio
        tai_dic13flx = tai_dic13flx*rntatio
        tai_dicwflx = tai_dicwflx*rntatio
        tai_dic13wflx = tai_dic13wflx*rntatio
        tai_alk = tai_alk*rntatio
        tai_o2 = tai_o2*rntatio
        tai_o2flx = tai_o2flx*rntatio
        tai_po4 = tai_po4*rntatio
# if defined O_mobi_silicon
        tai_sil = tai_sil*rntatio
        tai_silflx = tai_silflx*rntatio
        tai_opl = tai_opl*rntatio
# endif        
        tai_dop = tai_dop*rntatio
        tai_phyt = tai_phyt*rntatio
        tai_zoop = tai_zoop*rntatio
        tai_detr = tai_detr*rntatio
        tai_no3 = tai_no3*rntatio
        tai_don = tai_don*rntatio
        tai_diaz = tai_diaz*rntatio
        tai_din15 = tai_din15*rntatio
        tai_don15 = tai_don15*rntatio
        tai_phytn15 = tai_phytn15*rntatio
        tai_diatn15 = tai_diatn15*rntatio
        tai_zoopn15 = tai_zoopn15*rntatio
        tai_detrn15 = tai_detrn15*rntatio
        tai_diazn15 = tai_diazn15*rntatio
        tai_doc13 = tai_doc13*rntatio
        tai_phytc13 = tai_phytc13*rntatio
        tai_zoopc13 = tai_zoopc13*rntatio
        tai_detrc13 = tai_detrc13*rntatio
        tai_diazc13 = tai_diazc13*rntatio
        tai_c14 = tai_c14*rntatio
        tai_dc14 = tai_dc14*rntatio
        tai_c14flx = tai_c14flx*rntatio
        tai_cfc11 = tai_cfc11*rntatio
        tai_cfc11flx = tai_cfc11flx*rntatio
        tai_cfc12 = tai_cfc12*rntatio
        tai_cfc12flx = tai_cfc12flx*rntatio
        tai_sspH = tai_sspH*rntatio
        tai_ssCO3 = tai_ssCO3*rntatio
        tai_ssOc = tai_ssOc*rntatio
        tai_ssOa = tai_ssOa*rntatio
        tai_sspCO2 = tai_sspCO2*rntatio
        tai_dfe = tai_dfe*rntatio
        tai_ddfe = tai_ddfe*rntatio
# if defined O_mobi_silicon
        tai_diatc13 = tai_diatc13*rntatio
# endif
# if defined O_mobi_caco3
        tai_caco3 = tai_caco3*rntatio
        tai_caco3c13 = tai_caco3c13*rntatio
# endif
# if defined O_kk_ballast
        tai_d_B = tai_d_B*rntatio
# endif
# if defined O_mobi_silicon
        tai_diat = tai_diat*rntatio
# endif        
      endif

      return
      end
#endif
#if defined O_mom && defined O_tai_slh
      real function tinsit (ti, si, di)
!=======================================================================
!     function to calculate insitu temperature using Newtons method

!     input:
!       ti = potential temperature (C)
!       si = salinity (psu)
!       di = depth (m)

!     output:
!       tinsit = institu temperataure
!=======================================================================

      implicit none

      integer iter

      real, intent(in) :: di, si, ti
      real(kind=8) d, dfdt, error, s, t, ta, tb, tol

      error = 1.
      iter = 0
      tol = 1.e-5
      t = ti
      d = di
      s= si
      do while (abs(error) .gt. tol .and. iter .lt. 25)
        iter = iter + 1
        call potem(t, s, d, ta)
        call potem(t+1.e-3 ,s, d, tb)
        dfdt = (tb - ta)*1.e3
        error = (ta - ti)/dfdt
        t = t - error
      enddo
      if (abs(error) .gt. tol) then
        print*, "=> Error: insitu temperature did not converge"
        print*, "          ti = ", ti
        print*, "          si = ", si
        print*, "          di = ", di
        print*, "          error = ", abs(error)
        print*, "          tinsit = ", t
        print*, "   Stopped in function tinsit in diago.F"
        stop
      endif
      tinsit = t

      return
      end
#endif
