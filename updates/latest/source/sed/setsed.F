      subroutine setsed (is, ie, js, je)

#if defined O_sed
!-----------------------------------------------------------------------
!     set up everything which must be done only once per run
!-----------------------------------------------------------------------

      implicit none

      integer is, ie, js, je

      character(120) :: fname, new_file_name

      integer i, ioun, j, k

      logical exists, inqvardef

      real expb, db, difo2, difc(3)

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "csbc.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#include "sed.h"
#include "switch.h"
#include "tmngr.h"

      dtsed = dtsedyr*yrlen*daylen
      nsedacc = max(nsedacc, 1)
# if !defined O_sed_uncoupled
      if (nsedacc .gt. 1) then
        print*, '==> Warning: sediment model accelerated and coupled'
      endif
# endif
!-----------------------------------------------------------------------
!     initialise sediment model
!-----------------------------------------------------------------------
      weathflx = 0.  ! weathering flux in umol s-1
      sed_year = 0.
      carb(:,1,:) = 20.E-6
      carb(:,2,:) = 2000.E-6
      carb(:,3,:) = 80.E-6
      dcpls(:,:,:) = 0.
      dcmin(:,:,:) = 0.
      pore(:,:) = 0.
      form(:,:) = 0.
      o2(:,:) = 150.e-6
      orgml(:,:) = 0.
      calml(:,:) = 0.
      dopls(:,:) = 0.
      domin(:,:) = 0.
      dbpls(:,:) = 0.
      dbmin(:,:) = 0.
      co3_p(:) = 0.
      atsed = 0
      ipsed = 0
      map_sed(:,:) = 0
      calgg(:,:) = 0.1
      orggg(:,:) = 0.0002
      buried_mass(:,:) = 0.
      buried_calfrac(:,:) = 0.1
      buried_mass(1,:) = 500.
      kmin = max(1,kmin)
      ttrcal(:) = 0.
      rain_cal_p(:) = 0.
      sedsa = 0.
      carblith = 0.
#ifndef O_TMM
      do j=2,jmtm1
        do i=2,imtm1
#else
      do j=js,je
        do i=is,ie
#endif        
          if (kmt(i,j) .ge. kmin) then
            ipsed = ipsed + 1
            map_sed(i,j) = ipsed
            imap(ipsed) = i
            jmap(ipsed) = j
            water_z_p(ipsed) = zw(kmt(i,j))*0.01
 	    sedsa = sedsa + dxt(i)*dyt(j)*cst(j)
          endif
        enddo
      enddo

      delz(:) = 3.
      delz(1:5) = (/0.,.5,.5,1.,2./)
      kmax = 7
      zsed(1) = 0.
      do k=2,kmax
        zsed(k) = zsed(k-1) + delz(k)
      enddo
      zrct(:) = zsed(kmax)
      dissc = 1.1574e-5
      dissn = 4.5

!     depth_age is only used with the option "sed_profile".
!     make sure you have enough levels for the profile you want to simulate.
!     sediment is stored inverse to the depth_age. the latest sediment
!     is stored incrementing toward ibmax. once ibmax is reached all new
!     sediment is stored at ibmax. for example if sed_year is 120000 and
!     ibmax is 20 with a time increment of 4000, then the sediment in
!     level 1 is older than 120000 (usually the initial condition), in
!     level 2 it is between 120000 and 116000, level 3 between 116000 and
!     112000, level ibmax-1 between 76000 and 72000 and level ibmax between
!     76000 and 0. if sed_year is 20000 then the sediment in level 1 is
!     older than 20000, in level 2 it is between 20000 and 16000, in
!     level 5 between 8000 and 4000, in level 6 between 4000 and 0 and
!     levels higher than 6 would be available for future accumulation.
      do k=1,ibmax
        depth_age(k) = 4000.*(k-1)
      enddo

      expb = 3.
      difo2 = 12.1e-6
      difc(1) = 10.5e-6
      difc(2) = 6.4e-6
      difc(3) = 5.2e-6
      db = 0.15

      call set_pore (calgg, zsed, pore, kmax, ipsed, ipmax, nzmax)
      call pore_2_form (pore, form, kmax, expb, ipsed, ipmax, nzmax)
      call sldcon (orgml, orggg, 12., pore, kmax, ipsed, ipmax, nzmax)
      call sldcon (calml,calgg, 100., pore, kmax, ipsed, ipmax, nzmax)
      call calc_do2 (difo2, form, pore, delz, kmax, dopls, domin, ipsed
     &,              ipmax, nzmax)
      call calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin, ipsed
     &,             ipmax, nzmax)
      call calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin, ipsed
     &,             ipmax, nzmax)
!     calculate mass of mixed layer (g/cm2)
      call get_sed_ml_mass (delz, pore, sed_ml_mass, nzmax, ipsed
     &,                     ipmax, kmax)

!     zero boundary condition accumulators
      sbc(:,:,ibtemp) = 0.
      sbc(:,:,ibsalt) = 0.
      sbc(:,:,ircal) = 0.
      sbc(:,:,irorg) = 0.
# if defined O_carbon
      sbc(:,:,ibdic) = 0.
# endif
# if defined O_npzd_alk
      sbc(:,:,ibalk) = 0.
# endif
# if defined O_npzd_o2
      sbc(:,:,ibo2) = 0.
# endif

#ifndef O_TMM
      if (.not. init) then
        fname = new_file_name ("restart_sed.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) call sed_rest_in (fname, is, ie, js, je)
# if defined O_restart_2
        fname = new_file_name ("restart_2_sed.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) call sed_rest_in (fname, is, ie, js, je)
# endif
      endif
#endif ! !O_TMM
      
# if !defined O_sed_weath_diag

!-----------------------------------------------------------------------
!     if defined, use namelist value for constant weathering flux
!-----------------------------------------------------------------------

!     convert from kg/s to umol/s
      if (weath .lt. 1.e20) weathflx = weath*1.e9/12.

# endif
# if defined O_time_averages
!-----------------------------------------------------------------------
!     zero time average accumulators
!-----------------------------------------------------------------------

      call ta_sed_tavg (is, ie, js, je, 0)

# endif
# if defined O_time_step_monitor
!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_sed_tsi (is, ie, js, je, 0)

# endif
#endif

      return
      end
