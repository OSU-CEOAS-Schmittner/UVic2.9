      subroutine globalsum (ind)

#if defined O_global_sums || defined O_co2emit_diag
!=======================================================================
!     calculate the global sums of heat and fresh water

!     input:
!       ind = index for calculation (1=initial, 2=current, 3=final)
!=======================================================================

      implicit none

      character(120) :: fname
      character(32) :: nstamp

      integer i, iou, id_time, ind, j, jrow, k, l, n, nc, ntrec
      integer nyear, nmonth, nday, nhour, nmin, nsec, it(10)
      integer layer

      logical defined

      real area, cstdyt, dtoih1, dtoih2, dtoic1, dtoic2
      real dtah1, dtah2, dtaf1, dtaf2, dtac1, dtac2
      real dtsh1, dtsh2, dtsf1, dtsf2, dtsc1, dtsc2
      real dtih1, dtih2, dtif1, dtif2, dtic1, dtic2
      real dtlh1, dtlh2, dtlf1, dtlf2, dtlc1, dtlc2
      real dtoh1, dtof1, dtoh2, dtof2, dtoc1, dtoc2
      real dth1,  dth2,  dtf1,  dtf2,  dtc1,  dtc2
      real tah(3), taf(3), tac(3), tsh(3), tsf(3), tsc(3), tih(3)
      real tif(3), tic(3), tlh(3), tlf(3), tlc(3), toh(3), tof(3)
      real toc(3), th, tf, tc, RFTIME, time, tmp, c1e20, vol, avgper
# if defined O_mobi_silicon
      real dtos1, dtos2, dts1, dts2, dtois1, dtois2, tos(3)
      save tos, dtois1, ts
# endif
      save tah, taf, tac, tsh, tsf, tsc, tih, tif, tic, tlh, tlf, tlc
      save toh, tof, toc, dtoih1, dtoic1

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
# if defined O_embm
      include "atm.h"
# endif
# if defined O_ice && defined O_embm
#  if defined O_ice_cpts
      include "cpts.h"
#  endif
      include "ice.h"
# endif
      include "calendar.h"
      include "coord.h"
      include "grdvar.h"
      include "levind.h"
      include "csbc.h"
      include "iounit.h"
      include "switch.h"
# if defined O_mom
      include "mw.h"
#  if defined O_mobi
      include "mobi.h"
#  endif
# endif
# if defined O_mtlm
      include "mtlm.h"
# endif
      include "tmngr.h"

!-----------------------------------------------------------------------
!     sum heat, water, carbon in atmosphere, snow, ice, land, ocean
!-----------------------------------------------------------------------

      c1e20 = 1.e20
      if (ind .eq. 1) then
        dtoih = 0.0
        dtoih1 = 0.0
        dtoic = 0.
        dtoic1 = 0.0
      endif
# if defined O_mobi_silicon
      if (ind .eq. 1) then
        dtois = 0.
        dtois1 = 0.0
      endif
# endif
!     zero accumulators

      do n=min(3,ind*2-1),3
        tah(n) = 0.
        taf(n) = 0.
        tac(n) = 0.
        tsh(n) = 0.
        tsf(n) = 0.
        tsc(n) = 0.
        tih(n) = 0.
        tif(n) = 0.
        tic(n) = 0.
        tlh(n) = 0.
        tlf(n) = 0.
        tlc(n) = 0.
        toh(n) = 0.
        tof(n) = 0.
        toc(n) = 0.
# if defined O_mobi_silicon
        tos(n) = 0.
# endif
      enddo

# if defined O_mtlm
#  if defined O_mtlm_segday
      RFTIME = DAY_TRIF*segtim/DAY_YEAR
#  else
      RFTIME = DAY_TRIF/DAY_YEAR
#  endif
# endif

!     calculate totals

      L = 0
      do jrow=2,jmtm1
# if defined O_mom
        if (wide_open_mw) then
          j = jrow
        else
          j = jmw
          call getrow (latdisk(taup1disk), nslab, jrow
     &,                u(1,1,j,1,taup1), t(1,1,j,1,taup1))
        endif
# endif
        cstdyt = cst(jrow)*dyt(jrow)
        do i=2,imtm1
          area = cstdyt*dxt(i)
          if (kmt(i,jrow) .eq. 0) L = L + 1

!         atmosphere
# if defined O_embm
          tah(3) = tah(3) + at(i,jrow,2,isat)*area
          taf(3) = taf(3) + at(i,jrow,2,ishum)*area
#  if defined O_carbon_co2_2d
          tac(3) = tac(3) + at(i,jrow,2,ico2)*area
#  else
          tac(3) = tac(3) + co2ccn*area
#  endif
# endif

!         snow and ice
# if defined O_ice && defined O_embm
          tsf(3) = tsf(3) + hsno(i,jrow,2)*area
          tif(3) = tif(3) + hice(i,jrow,2)*area
#  if defined O_ice_cpts
          do nc=1,ncat
            tsf(3) = tsf(3) + hseff(i,jrow,2,nc)*area
            tif(3) = tif(3) + heff(i,jrow,2,nc)*area
          enddo
          do layer=1,ntilay
            tih(3) = tih(3) + E(i,jrow,2,layer)*area
     &             + hice(i,jrow,2)*area
          enddo
#  endif
# endif
# if defined O_mtlm && defined O_embm
          if (kmt(i,jrow) .eq. 0)
     &      tsf(3) = tsf(3) + LYING_SNOW(L)*area*0.1/rhosno
# endif

# if defined O_embm
!         land
          tlf(3) = tlf(3) + soilm(i,jrow,2)*area*1.0e-3
#  if defined O_mtlm
          if (kmt(i,jrow) .eq. 0) then
            tlh(3) = tlh(3) + (TS1(L)*HCAP_SOIL*ROOTDEP
!           add in most recent fluxes not felt by atmosphere yet
     &             + (sbc(i,jrow,ievap)*vlocn + sbc(i,jrow,isens)
     &             + sbc(i,jrow,ilwr))*segtim*SEC_DAY*1.e-3)*area*1.e-4
            tlf(3) = tlf(3) + (M(L) + MNEG(L)
!           add back in most recent fluxes not felt by atmosphere yet
     &             + sbc(i,jrow,ievap)*segtim*SEC_DAY*10.)*area*1.e-4
            tlc(3) = tlc(3) + (CV(L) + CS(L)
#   if defined O_carbon
!           add in most recent fluxes not felt by atmosphere yet
#    if defined O_mtlm_segday
     &             + (sbc(i,jrow,isr) - sbc(i,jrow,inpp)
     &             + sbc(i,jrow,iburn)*(DAY_TRIF-ntlbc))*segtim*SEC_DAY
#    else
     &             + ((sbc(i,jrow,isr) - sbc(i,jrow,inpp))*segtim
     &             + sbc(i,jrow,iburn)*(DAY_TRIF-ntlbc*segtim))*SEC_DAY
#    endif
#   endif
!           add in driving fluxes not felt by triffid yet
     &             - RESP_S_DR(L)*RFTIME
     &               )*area*1.e-4
            do N=1,NPFT
              tlc(3) = tlc(3) + NPP_DR(L,N)*FRAC(L,N)*RFTIME*area*1.0e-4
            enddo
          endif
#  endif
# endif

!         ocean
# if defined O_mom
          do k=1,km
            vol = dzt(k)*area
            toh(3) = toh(3) + t(i,k,j,itemp,taup1)*vol
            tof(3) = tof(3) + t(i,k,j,isalt,taup1)*vol
#  if defined O_carbon
            toc(3) = toc(3) + t(i,k,j,idic,taup1)*vol
#   if defined O_mobi
            toc(3) = toc(3) + t(i,k,j,iphyt,taup1)*vol*redctn
            toc(3) = toc(3) + t(i,k,j,izoop,taup1)*vol*redctn
            toc(3) = toc(3) + t(i,k,j,idetr,taup1)*vol*redctn
#    if defined O_mobi_caco3
            toc(3) = toc(3) + t(i,k,j,icaco3,taup1)*vol*1.e-3
! attached calcite
            toc(3) = toc(3) + t(i,k,j,izoop,taup1)*vol*redctn*capr
#     if defined O_kk_ballast
            toc(3) = toc(3) + t(i,k,j,idetr_B,taup1)*vol*redctn
#     endif
#    endif
#    if defined O_mobi_silicon
            toc(3) = toc(3) + t(i,k,j,idiat,taup1)*vol*redctn
#    endif
#    if defined O_mobi_nitrogen
            toc(3) = toc(3) + t(i,k,j,idiaz,taup1)*vol*redctn
#    endif
#   endif
#  endif
#  if defined O_mobi_silicon
            tos(3) = tos(3) + t(i,k,j,isil,taup1)*vol
#  endif
          enddo
#  if defined O_sealev_data && defined O_sealev_salinity
          if (kmt(i,jrow) .ne. 0)
     &      tof(3) = tof(3) + sealev*area*socn
#  endif
# endif
        enddo
      enddo

!     convert units to Joules and kilograms

# if defined O_embm
!     atmosphere
      tah(3) = (taf(3)*rhoatm*shq*vlocn+tah(3)*cpatm*rhoatm*sht)*1.0e-7
      taf(3) = taf(3)*rhoatm*shq*1.0e-3
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      tac(3) = tac(3)*4.138e-7*rhoatm*shc*1.0e-3

# endif
# if defined O_ice && defined O_embm
!     snow and ice
      tsh(3) = -tsf(3)*rhosno*flice*1.0e-7
      tsf(3) = tsf(3)*rhosno*1.0e-3
#  if defined O_ice_cpts
      tih(3) = tih(3)*1.0e-7   ! cpts ice energy of melting is negative
#  else
      tih(3) = -tif(3)*rhoice*flice*1.0e-7
#  endif
      tif(3) = tif(3)*rhoice*1.0e-3

# endif
!     land
      tlh(3) = tlh(3)
      tlf(3) = tlf(3)
      tlc(3) = tlc(3)

# if defined O_mom
!     ocean
      toh(3) = toh(3)/0.2389
      tof(3) = -tof(3)/(socn*1.e3)
!     12e-9 kg/umol
      toc(3) = toc(3)*12.e-9
#   if defined O_mobi_silicon
      tos(3) = tos(3)*28.0855e-9
#   endif
# endif
!     dtoic includes any sediment burial, weathering or emissions
!     12e-9 kg/umol
      dtoic = dtoic*12e-9

!     calculate differences from start of run

      dtah1 = tah(3) - tah(1)
      dtaf1 = taf(3) - taf(1)
      dtac1 = tac(3) - tac(1)
      dtsh1 = tsh(3) - tsh(1)
      dtsf1 = tsf(3) - tsf(1)
      dtsc1 = tsc(3) - tsc(1)
      dtih1 = tih(3) - tih(1)
      dtif1 = tif(3) - tif(1)
      dtic1 = tic(3) - tic(1)
      dtlh1 = tlh(3) - tlh(1)
      dtlf1 = tlf(3) - tlf(1)
      dtlc1 = tlc(3) - tlc(1)
      dtoh1 = toh(3) - toh(1)
      dtof1 = tof(3) - tof(1)
      dtoc1 = toc(3) - toc(1)
      dtoih1 = dtoih1 + dtoih
      dtoic1 = dtoic1 + dtoic
# if defined O_mobi_silicon
      dtois = dtois*28.0855e-9 !28 kg/umol
      dtos1 = tos(3) - tos(1)
      dtois1 = dtois1 + dtois
# endif
!     calculate differences from last time step

      dtah2 = tah(3) - tah(2)
      dtaf2 = taf(3) - taf(2)
      dtac2 = tac(3) - tac(2)
      dtsh2 = tsh(3) - tsh(2)
      dtsf2 = tsf(3) - tsf(2)
      dtsc2 = tsc(3) - tsc(2)
      dtih2 = tih(3) - tih(2)
      dtif2 = tif(3) - tif(2)
      dtic2 = tic(3) - tic(2)
      dtlh2 = tlh(3) - tlh(2)
      dtlf2 = tlf(3) - tlf(2)
      dtlc2 = tlc(3) - tlc(2)
      dtoh2 = toh(3) - toh(2)
      dtof2 = tof(3) - tof(2)
      dtoc2 = toc(3) - toc(2)
      dtoih2 = dtoih
      dtoic2 = dtoic
# if defined O_mobi_silicon
      dtos2 = tos(3) - tos(2)
      dtois2 = dtois
      dtois = 0.0
# endif

      dtoih = 0.0
      dtoic = 0.0
      do n=ind,2
        tah(n) = tah(3)
        taf(n) = taf(3)
        tac(n) = tac(3)
        tsh(n) = tsh(3)
        tsf(n) = tsf(3)
        tsc(n) = tsc(3)
        tih(n) = tih(3)
        tif(n) = tif(3)
        tic(n) = tic(3)
        tlh(n) = tlh(3)
        tlf(n) = tlf(3)
        tlc(n) = tlc(3)
        toh(n) = toh(3)
        tof(n) = tof(3)
        toc(n) = toc(3)
# if defined O_mobi_silicon
        tos(n) = tos(3)
# endif
      enddo

!     write differences or totals

      if (ind .eq. 2) then

        dth1 = dtah1 + dtsh1 + dtih1 + dtlh1 + dtoh1 + dtoih1
        dtf1 = dtaf1 + dtsf1 + dtif1 + dtlf1 + dtof1
        dtc1 = dtac1 + dtsc1 + dtic1 + dtlc1 + dtoc1 + dtoic1
        dth2 = dtah2 + dtsh2 + dtih2 + dtlh2 + dtoh2 + dtoih2
        dtf2 = dtaf2 + dtsf2 + dtif2 + dtlf2 + dtof2
        dtc2 = dtac2 + dtsc2 + dtic2 + dtlc2 + dtoc2 + dtoic2
# if defined O_mobi_silicon
        dts1 = dtos1 + dtois1
        dts2 = dtos2 + dtois2
# endif
        if (tsits) then

          time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
          call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
          nyear = time
          call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

          call def_tsi
          call def_tsi_gsums (fname)

#  if !defined O_save_time_relyear0
!         make output time relative to year 1
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

          call opennext (fname, time, ntrec, iou)
          if (ntrec .le. 0) ntrec = 1

          call putvars ('time', iou, ntrec, time, c1, c0)
          call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
          call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
# if defined O_co2emit_diag
          call putvars ('F_co2diag', iou, ntrec
#  if defined O_save_carbon_totals
     &,     dtc2, tsiint*1.e12/yrlen, c0)
#  else
     &,     dtc2, tsiint*daylen, c0)
#  endif
          carbemit = carbemit + dtc2*1e-12
          dtoic1 = dtoic1 - dtc1
          dtoic2 = dtoic2 - dtc2
          dtc1 = 0.
          dtc2 = 0.
# endif
# if defined O_global_sums
          call putvars ('G_dtheat', iou, ntrec, dth1, c1, c0)
          call putvars ('G_dtwater', iou, ntrec, dtf1, c1, c0)
          call putvars ('G_dtcarb', iou, ntrec, dtc1, c1, c0)
#   if defined O_mobi_silicon
          call putvars ('G_dtsilica', iou, ntrec, dts1, c1, c0)
#   endif
# endif

# if defined O_global_sums
          if (iotsi .eq. stdout .or. iotsi .lt. 0) then
            write (*,'(a,a)') 'Changes in heat, fresh water and carbon '
     &,       'from start of run and last calculation'
            write (*,'(a,e22.14,a,e22.14,a,e22.14,a
     &,                  e22.14,a,e22.14,a,e22.14,a)')
     &        ' d atm    ',dtah1, ' J ',dtaf1,' kg ',dtac1, ' kg'
     &,                    dtah2, ' J ',dtaf2,' kg ',dtac2, ' kg'
     &,       ' d snow   ',dtsh1, ' J ',dtsf1,' kg ',dtsc1, ' kg'
     &,                    dtsh2, ' J ',dtsf2,' kg ',dtsc2, ' kg'
     &,       ' d ice    ',dtih1, ' J ',dtif1,' kg ',dtic1, ' kg'
     &,                    dtih2, ' J ',dtif2,' kg ',dtic2, ' kg'
     &,       ' d lnd    ',dtlh1, ' J ',dtlf1,' kg ',dtlc1, ' kg'
     &,                    dtlh2, ' J ',dtlf2,' kg ',dtlc2, ' kg'
     &,       ' d ocn    ',dtoh1, ' J ',dtof1,' kg ',dtoc1, ' kg'
     &,                    dtoh2, ' J ',dtof2,' kg ',dtoc2, ' kg'
     &,       ' d out-in ',dtoih1,' J ',0.0,  ' kg ',dtoic1,' kg'
     &,                    dtoih2,' J ',0.0,  ' kg ',dtoic2,' kg'
     &,       ' d total  ',dth1,  ' J ',dtf1, ' kg ',dtc1,  ' kg'
     &,                    dth2,  ' J ',dtf2, ' kg ',dtc2,  ' kg'
#   if defined O_save_silica_totals
!KK This doesn't make any sense as long as Si is only in the ocean
            write (*,'(a,a)') 'Changes in silica '
     &,       'from start of run and last calculation'
            write (*,'(a,e22.14,a,e22.14,a)')
     &        ' d ocn    ',dtos1, ' kg'
     &,                    dtos2, ' kg'
     &,       ' d out-in ',dtois1,' kg'
     &,                    dtois2,' kg'
     &,       ' d total  ',dts1,  ' kg'
     &,                    dts2,  ' kg'
#   endif
          endif
        endif

      else

        th = tah(ind) + tsh(ind) + tih(ind) + tlh(ind) + toh(ind)
        tf = taf(ind) + tsf(ind) + tif(ind) + tlf(ind) + tof(ind)
        tc = tac(ind) + tsc(ind) + tic(ind) + tlc(ind) + toc(ind)
#  if defined O_mobi_silicon
        ts = tos(ind)
#  endif
        write (*,'(/,a,a)') 'Total heat (in Joules referenced to 0 C'
     &,   ' and no ice or snow) and fresh water (in kg)'
#  if defined O_embm
        write (*,'(a,a)') 'Total ocean fresh water is the equivalent'
     &,   ' difference from the ocean volume referenced to socn'
#  else
        write (*,'(a,a)') 'Total ocean fresh water is the equivalent'
     &,   ' difference from the ocean volume referenced to 34.9 ppt'
#  endif
        write (*,'(a,e22.14,a,e22.14,a,e22.14,a)')
     &    ' t atm    ', tah(ind),' J ',taf(ind),' kg ',tac(ind),' kg'
     &,   ' t snow   ', tsh(ind),' J ',tsf(ind),' kg ',tsc(ind),' kg'
     &,   ' t ice    ', tih(ind),' J ',tif(ind),' kg ',tic(ind),' kg'
     &,   ' t lnd    ', tlh(ind),' J ',tlf(ind),' kg ',tlc(ind),' kg'
     &,   ' t ocn    ', toh(ind),' J ',tof(ind),' kg ',toc(ind),' kg'
     &,   ' t total  ', th,      ' J ',tf,      ' kg ',tc,      ' kg'
#   if defined O_mobi_silicon
!     &,   ' ocn sil  ', tos(ind),' kg tot sil', ts, ' kg'
#   endif
# else
        endif
# endif
      endif

      return
      end

      subroutine gsums_tsi_def (fname, calendar, expnam, runstamp)

!=======================================================================
!     output routine for atmospheric global sums step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou

      real c0, c1, c1e3, c1e20

      c0 = 0.
      c1 = 1.
      c1e3 = 1.e3
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      id(1) = id_time

!-----------------------------------------------------------------------
!     define data
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'
# if defined O_units_time_years
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'years since 1-1-1')
#  else
     &, 'time', 'time', 'years since 0-1-1')
#  endif
# else
#  if !defined O_save_time_relyear0
     &, 'time', 'time', 'days since 1-1-1')
#  else
     &, 'time', 'time', 'days since 0-1-1')
#  endif
# endif
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_avgper', iou, 1, id, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','day')
# if defined O_co2emit_diag
      call defvar ('F_co2diag', iou, 1, id, -c1e20, c1e20, ' '
#  if defined O_save_carbon_totals
     &, 'F', 'diagnosed carbon emissions', ' ', 'Pg yr-1')
#  else
     &, 'F', 'diagnosed carbon emissions', ' ', 'kg s-1')
#  endif
# endif
# if defined O_global_sums
      call defvar ('G_dtheat', iou, 1, id, -c1e20
     &, c1e20, ' ', 'F', 'heat conservation error', ' ', 'J')
      call defvar ('G_dtwater', iou, 1, id, -c1e20
     &, c1e20, ' ', 'F', 'water conservation error', ' ', 'kg')
      call defvar ('G_dtcarb', iou, 1, id, -c1e20
     &, c1e20, ' ', 'F', 'carbon conservation error', ' ', 'kg')
#   if defined O_mobi_silicon
      call defvar ('G_dtsilica', iou, 1, id, -c1e20
     &, c1e20, ' ', 'F', 'silica conservation error', ' ', 'kg')
#   endif
# endif

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

#endif
      return
      end
