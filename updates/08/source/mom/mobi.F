      subroutine mobi_init

#if defined O_mom && defined O_mobi 
!     initialize MOBI parameters
      implicit none
#include "size.h"
#include "mobi.h"
#include "calendar.h"
#include "coord.h"
#include "stdunits.h"
#include "scalar.h"

      integer k, ioun, imobi
      real alpha, par
      real alpha_C, sumzpref

# if defined O_mobi_iron || defined O_mobi_silicon
!     IO for fe_hydr
#include "param.h"
#include "pconst.h"
      integer ib(10), ic(10), iou, j, i
      character (120) :: fname, new_file_name
      logical exists
      real c1e3
      real tmpijkm(1:imtm2,1:jmtm2,1:km)
! end IO for fe_hydr
# endif
      namelist /mobi/   alpha, kw, kc, abio_P, bbio, cbio, nup
     &,                 gamma1, gbio, epsbio, nuz, nud0, dbct_D
     &,                 wd0, par, dtnpzd, redctn, ki, redptn, caprmax
     &,                 kcapr, dcaco3, nupt0, redotn, jdiar, kzoo, geZ
     &,                 zprefP, zprefZ, zprefDet, zprefDiaz, kfe, kfe_D
     &,                 knmax, knmin
# if defined O_mobi_iron
     &,                 kfeleq, kfecol, alphamax, alphamin
     &,                 thetamaxhi, rfeton, fetopsed, o2min
     &,                 lig, thetamaxlo, mc, kfeorg
     &,                 kfemax, kfemin, pmax
#  if defined O_mobi_silicon
     &,                 knmax_Diat, knmin_Diat
     &,                 kfemin_Diat, kfemax_Diat
#  endif
# endif
     &,                 nupt0_D, sgbdfac, diazntp, nup_D, dfr, dfrt
     &,                 nudon0, nudop0, eps_assim, eps_excr, eps_nfix
     &,                 eps_wcdeni, eps_bdeni0, eps_recy, hdop
     &,                 mw, mwz
# if defined O_kk_ballast
     &,                 bapr
# endif
# if defined O_mobi_silicon
     &,                 abiodiat, nudt0, nu_diat, alpha_Diat
     &,                 zprefDiat, k1si, sildustfluxfac
     &,                 sipr0, dopal, si_sol, si_h_sol, opl_disk0, wo0
# endif
     &,                 wc0, kcal, dissk0, kc_c
!     set defaults for namelist npzd

      alpha = 0.16  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      kw = 0.04  ! Light attenuation due to water [1/m]
      kc = 0.047  ! Light atten. by phytoplankton [1/(m*mmol/m^3)]
      ki = 5.0  ! Light attenuation through sea ice & snow
      abio_P = 0.6  ! a; Maximum growth rate parameter [1/day]
      bbio = 1.066  ! b
      cbio = 1.0  ! c [1/deg_C]
      k1n = 0.7  ! Half saturation constant for N uptake [mmol/m^3]
      nup = 0.03  ! Specific mortality rate (Phytoplankton) [1/day]
      nup_D = 0.0001 !  Specific mortaility rate (Diazotrophs [1/day]
      nupt0 = 0.015  ! Fast-recycling  mortality rate (Phytoplankton) [1/day]
      nupt0_D = 0.001 ! Fast-recyling mortality rate (Diazotrophs) [1/day]
      gamma1 = 0.70  ! gama1; Assimilation efficiency (zpk)
      gbio = 0.38  ! Maximum grazing rate [1/day]
      epsbio = 1.6  ! Prey capture rate [(mmol/m^3)^(-2)day^(-1)]
      nuz = 0.06  ! Quadratic mortality (zpk) [(mmol/m^3)^(-2)day^(-1)]
      nud0 = 0.07  ! detritus remineralization rate [1/day]
      nudon0 = 2.33e-5 ! DON remineralization rate [1/day]
      nudop0 = 7.e-5  ! DOP remineralization rate [1/day]
      wd0 = 16.0  ! Sinking speed of detritus at surface [m/day]
      mwz = 100000. !Depth where sinking below remains constant (cm)
      mw = 0.02 !Sinking speed increase with depth (s-1)
      mw_c = 0.06 !Calcite sinking speed increase with depth (s-1)
      par = 0.43 ! fraction of photosythetically active radiation
      dtnpzd = dtts/4.  ! time step of biology [s]
      redctn = 7.1  ! C/N Redfield ratio (Schneider et al., 2003, GBC)
      redptn = 1./16.  ! P/N Redfield ratio
      caprmax = 0.022  ! carbonate to carbon production ratio
      kcapr = 0.4 ! half saturation for PIC/POC ratio (Gehlen et al. 2007, BGS)
      dcaco3 = 650000.0  ! remineralisation depth of calcite [cm]
      redotn = 10.6  ! O2/N Redfield ratio (Anderson & Sarmiento, 1994, GBC)
      jdiar = 0.08  ! factor reducing the growth rate of diazotrophs
      dbct_D = 2.6 ! subtracted from bct for diazotrophs
      kzoo = 0.15  ! half sat. constant for Z grazing [mmol N m^3]
      geZ = 0.6   ! Zooplankton growth efficiency
      sgbdfac = 1.0 ! sub-grid benthic denitrification rate factor 
      diazntp = 28. !diazotroph N:P ratio
      dfr = 0.08   ! phyt mortality refractory/semi-labile DOM fraction
      dfrt = 0.01  ! phyt fast-recy refracotyr/semi-labile DOM fraction
      hdop = 0.4   ! DOP growth rate handicap
# if defined O_mobi_silicon
      alpha_Diat = 0.1  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      abiodiat = 3.45  ! a; Maximum growth rate parameter diatoms [1/day]
      k1si = 2.0!0.3  ! Half sat const for Si uptake (diatoms) [mmol/m^3]
      k1n_Diat = 0.003  ! Half sat const for N uptake (diatoms) [mmol/m^3]
      nu_diat = 0.03  ! Specific mortality rate (diatoms) [1/day]
      nudt0 = 0.015  ! Specific mortality rate (diatoms) [1/day]
!      redstn = 0.9375  ! S/N Redfield ratio
      wo0 = 50. !Opal sinking rate in m/d
c      si_dis0 = 5.07e-7 !!S&G 2006 dissl rate per s
c     the following was estimated from an assumed e-folding depth
c     z` = 6 km = wo0/si_dis0
      opl_disk0 = 8.3e-3 !Andreas dissl rate per day = 9.7e-8 1/s
      si_h_sol = 1.0 !global average solubility ratio of hydrothermal input
      si_sol = 0.03 !global average solubility ratio of dust input 
      dopal = 60000.0  ! remineralisation depth of opal [cm]      
      sipr0 = 0.13  ! opal to carbon production ratio from Aumont et al GBC 2003
      sildustfluxfac = 1. ! scaling factor for silicon dust flux
# endif      
# if defined O_kk_ballast
      bapr = 0.05 ! detritus to carbonate ratio [mg POC/mg PIC]
# endif
      kc_c = 0.047  ! Light atten. by calcite [1/(m*mmol/m^3)]
      wc0 = 35. ! constant calcite sinking speed [m/day]
!      Gehlen et al 2007 use a surface sink speed of 50 m/day inc w/ depth
      kcal = 100.  !half sat for PIC dissolution from Planktom 10
      dissk0 = 0.013 !initial dissolution rate parameter [1/day] 

      zprefP = 0.18 ! Zooplankton preference for P
      zprefDiat = 0.18 ! Zooplankton preference for Diatoms
      zprefDiaz = 0.1 ! Zooplankton preference for diazotrophs
      zprefZ = 0.18 ! Zooplankton preference for other Z
      zprefDet = 0.18 ! Zooplankton preference for Detritus
# if !defined O_mobi_nitrogen
      zprefDiaz = 0.
# endif
# if !defined O_mobi_silicon
      zprefDiat = 0.
# endif

# if defined O_mobi_nitrogen
#  if !defined O_mobi_o2
      print*,'error in mobi_init:'
      print*,'O_mobi_o2 must be defined if O_mobi_nitrogen is used'
      stop 
#  endif
# endif

# if defined O_mobi_nitrogen_15
#  if !defined O_mobi_nitrogen
      print*,'error in mobi_init:'
      print*,'O_mobi_nitrogen must be on if O_mobi_nitrogen_15 is used'
      stop
#  endif      
      eps_assim = 6.
      eps_excr = 4.
      eps_nfix = 1.
      eps_wcdeni = 25.
      eps_bdeni0 = 6.
      eps_recy = 1.
# endif
# if defined O_mobi_iron
#  if !defined O_mobi_o2
      print*,'error in mobi_init:'
      print*,'O_mobi_o2 must be defined if O_mobi_iron is used'
      stop 
#  endif
      kfemin = 0.04e-3 ! Min. half sat. const. for Fe lim. [mmol Fe / m-3]
      kfemax = 0.2e-3  ! Max. half sat. const. for Fe lim. [mmol Fe / m-3]
      knmin = 0.15 ! Min. half sat. const. for N lim. [mmol N / m-3]
      knmax = 1.5  ! Max. half sat. const. for N lim. [mmol N / m-3]
      pmax = 0.15 ! Phyt biomass above which kfe increases [mmol N / m-3]
      kfe_D = 0.1e-3 ! Half sat. const. for Diaz Fe limit. [mmol Fe / m-3]
c     juan: Levin's value was 0.1e-3
#  if defined O_mobi_silicon
      kfemin_Diat = 0.04e-3
      kfemax_Diat = 0.8e-3
      knmin_Diat = 0.3
      knmax_Diat = 3.0
      pmax_Diat = 0.15
      alpha_Diat = 4.0  ! Chlorophyll specific initial slope P-I curve [(gC/gChl)^(-1) (W/m^2)^(-1) /day]     
      kfe_Diat = 1.87e-5 ! Half saturation constant for Diat Fe limitation      
#  endif
      kfeleq = 10.**5.5   !  Fe-ligand stability constant [m^3/(mmol ligand)]
      lig = 1.0e-3      !  Ligand concentration  [mmol/m^3]
      thetamaxhi = 0.04    !  Maximum Chl:C ratio, abundant iron [gChl/(gC)]
      thetamaxlo = 0.01    !  Maximum Chl:C ratio, extreme iron limitation [gChl/(gC)]
      alphamax = 73.6e-6*86400  !   Maximum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
      alphamin = 18.4e-6*86400   !   Minimum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
      mc = 12.011           ! Molar mass of carbon [g/mol]
      fetopsed = 0.004  !  Fe:P for sedimentary iron source [molFe/molP]
      o2min = 5.      !  Minimum O2 concentration for aerobic respiration [mmolO_2/m^3]
      kfeorg = 0.45*(1./86400.)     !  Organic-matter dependent scavenging rate [(m^3/(gC s))^0.58]
      rfeton = 10.e-6*6.625        ! Uptake ratio of iron to nitrogen [mol Fe/mol N] = 10 micromol Fe /mol C
      kfecol = 0.005/86400.     ! Colloidal production and precipitation rate [s^-1]
# endif
!     read namelist
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, mobi, end=108)
108   continue
      write (stdout, mobi)
      call relunit (ioun)

!     convert units of NPZD parameters to MOM units
!     carbon, alk, o2 are in mol, nitrogen, phosporous in mmol  
      redctn = redctn*1.e-3 ! mol/mmol
      redotn = redotn*1.e-3 ! mol/mmol
      redotp = redotn/redptn ! mol/mmol
      redctp = redctn/redptn ! mol/mmol
      redotc = redotn/redctn ! mol/mol
      redntp = 1./redptn ! mol/mol
      redntc = 1./redctn ! mmol/mol
      diazptn = 1./diazntp ! mol/mol
c      k1p_P   = k1n*redptn
# if defined O_mobi_caco3
#  if !defined O_carbon
      print*,'error in mobi_init:'
      print*,'O_carbon must be defined if O_mobi_caco3 is used'
      stop
#  endif      
      wc0 = wc0*1.e2
      kc_c = kc_c*1.e-2
      dissk0 = dissk0/daylen
# endif
#  if defined O_mobi_silicon
c      k1p_Diat  = k1n_Diat*redptn
      alpha_Diat = alpha_Diat/daylen
      abiodiat = abiodiat/daylen
      nu_diat = nu_diat/daylen
      nudt0 = nudt0/daylen
#  endif
#  if defined O_mobi_silicon
      wo0 = wo0*1.e2
      opl_disk0 = opl_disk0/daylen
#  endif
      kw = kw*1.e-2
      kc = kc*1.e-2
      ki = ki*1.e-2
      wd0 = wd0*1.e2
      alpha = alpha/daylen
      abio_P = abio_P/daylen
      nup = nup/daylen
      nup_D = nup_D/daylen
      nupt0 = nupt0/daylen
      nupt0_D = nupt0_D/daylen
      gbio = gbio/daylen
      epsbio = epsbio/daylen
      nuz = nuz/daylen
      nud0 = nud0/daylen
      nudop0 = nudop0/daylen
      nudon0 = nudon0/daylen
# if defined O_mobi_iron
      alphamax = alphamax/daylen
      alphamin = alphamin/daylen
      tap = 2.*par   ! alpha is iron dependent and calculated in npzd_src
# else
      tap = 2.*alpha*par
# endif


!     calculate sinking speed of detritus divided by layer thickness
      do k=1,km
!     linear increase wd0-200m with depth
         if (zt(k) .lt. mwz) then
            wd(k) = (wd0+mw*zt(k))/daylen/dzt(k) ! [s-1]
# if defined O_mobi_caco3
            wc(k) = (wc0+mw_c*zt(k))/daylen/dzt(k)    ! [s-1]
# endif            
         else
            wd(k) = (wd0+mw*mwz)/daylen/dzt(k) ! [s-1]
# if defined O_mobi_caco3
            wc(k) = (wc0+mw_c*mwz)/daylen/dzt(k)    ! [s-1]
# endif 
         endif
# if defined O_mobi_silicon
         wo(k) = wo0/daylen/dzt(k) ! [s-1]
# endif 
      enddo
      ztt(1)=0.0
      do k=1,km-1
         ztt(k+1)=(-1)*zw(k)
      enddo

!     check grazing preferences = 1 for N case
      sumzpref = zprefP + zprefDet + zprefZ + zprefDiaz + zprefC 
     &         + zprefDiaz
      if (sumzpref.ne.1.) then
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC, zprefDiat"
     &,          zprefP, zprefDet, zprefZ, zprefDiaz, zprefC, zprefDiat
         zprefP = zprefP/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         zprefDiaz = zprefDiaz/sumzpref
         zprefC = zprefC/sumzpref
         zprefDiat = zprefDiat/sumzpref
         print*,"to"
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC, zprefDiat"
     &,          zprefP, zprefDet, zprefZ, zprefDiaz, zprefC, zprefDiat
      endif

# ifndef O_TMM
#  if defined O_mobi_iron
! read the hydrothermal Fe input from O_fe_hydr.nc
      fe_hydr(:,:,:) = 0.
      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = km

      fname = new_file_name ("O_fe_hydr.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         c1e3 = 1000
         call openfile (trim(fname), iou)
         call getvara ('O_fe_hydr', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1e3, c0)
         fe_hydr(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)
            do k=1,km
               do j=1,jmt
                  fe_hydr(1,j,k) = fe_hydr(imtm1,j,k)
                  fe_hydr(imt,j,k) = fe_hydr(2,j,k)
               enddo
               do i=1,imt
                  fe_hydr(i,1,k) = fe_hydr(i,2,k)
                  fe_hydr(i,jmt,k) = fe_hydr(2,j,k)
               enddo
            enddo
      else
         print*,"Warning => Cannot find", trim(fname)
      endif
#  endif
#  if defined O_mobi_silicon
! read in Silicate hydrothermal input file as annual average mask
! read in as mol Si per m^3 per sec 
      si_hydr(:,:,:) = 0.
      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = km

      fname = new_file_name ("O_si_hydr.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (trim(fname), iou)
        call getvara ('O_SIL_HYDR', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1, c0)
        si_hydr(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2,1:jmtm2,:)
         where (si_hydr(:,:,:) .lt. -1.e30)
           si_hydr(:,:,:) = 0.0
         endwhere
        do k=1,km
          do j=1,jmt
            si_hydr(1,j,k) = si_hydr(imtm1,j,k)
            si_hydr(imt,j,k) = si_hydr(2,j,k)
          enddo
          do i=1,imt
            si_hydr(i,1,k) = si_hydr(i,2,k)
            si_hydr(i,jmt,k) = si_hydr(i,jmtm1,k)
          enddo
        enddo
      else
        print*,"Warning => Cannot find" , trim(fname)
      endif      
#  endif
# endif ! not O_TMM
# if defined O_mobi_silicon
!---------------------------------------------------------------------
!     calculate variables used in Opal remineralization
!---------------------------------------------------------------------      
!     Test a few options. All ropk(k) unitless. From
!     HAMOCC5 (Aumont et al. 2003 from Maier-Reimer et al. 1993)
      ropk(1) = -(exp(-zw(1)/dopal)-1.0)
      do k=2,km
        ropk(k) = -(exp(-zw(k)/dopal)) + (exp(-zw(k-1)/dopal))
      enddo        
# endif   
# if defined O_carbon || defined O_mobi_alk

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = 1./dzt(1)
      do k=2,km
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
# endif
      imobi=1
      call setimobi('imobipo4',imobipo4,imobi)
      call setimobi('imobiphyt',imobiphyt,imobi)
      call setimobi('imobizoop',imobizoop,imobi)
      call setimobi('imobidetr',imobidetr,imobi)
# if defined O_carbon
      call setimobi('imobidic',imobidic,imobi)
#  if defined O_carbon_13
      call setimobi('imobidic13',imobidic13,imobi)
      call setimobi('imobiphytc13',imobiphytc13,imobi)
      call setimobi('imobizoopc13',imobizoopc13,imobi)
      call setimobi('imobidetrc13',imobidetrc13,imobi)
#   if defined O_mobi_nitrogen
      call setimobi('imobidoc13',imobidoc13,imobi)
      call setimobi('imobidiazc13',imobidiazc13,imobi)
#   endif
#   if defined O_mobi_silicon
      call setimobi('imobidiatc13',imobidiatc13,imobi)
#   endif
#   if defined O_mobi_caco3
      call setimobi('imobicaco3c13',imobicaco3c13,imobi)
#   endif
#  endif
# endif
# if defined O_mobi_nitrogen
      call setimobi('imobidop',imobidop,imobi)
      call setimobi('imobino3',imobino3,imobi)
      call setimobi('imobidon',imobidon,imobi)
      call setimobi('imobidiaz',imobidiaz,imobi)
#  if defined O_mobi_nitrogen_15
      call setimobi('imobidin15',imobidin15,imobi)
      call setimobi('imobidon15',imobidon15,imobi)
      call setimobi('imobiphytn15',imobiphytn15,imobi)
      call setimobi('imobizoopn15',imobizoopn15,imobi)
      call setimobi('imobidetrn15',imobidetrn15,imobi)
      call setimobi('imobidiazn15',imobidiazn15,imobi)
#   if defined O_mobi_silicon
      call setimobi('imobidiatn15',imobidiatn15,imobi)
#   endif
#  endif
# endif
# if defined O_mobi_caco3
      call setimobi('imobicaco3',imobicaco3,imobi)
# endif
# if defined O_kk_ballast
      call setimobi('imobidetr_B',imobidetr_B,imobi)
# endif
# if defined O_mobi_silicon
      call setimobi('imobidiat',imobidiat,imobi)
      call setimobi('imobisil',imobisil,imobi)
      call setimobi('imobiopl',imobiopl,imobi)
# endif
# if defined O_mobi_iron
      call setimobi('imobidfe',imobidfe,imobi)
      call setimobi('imobidetrfe',imobidetrfe,imobi)
# endif
      if (imobi-1 .eq. ntnpzd) then
         print*,'total number of MOBI tracers =', imobi-1
      else
         print*,'error in mobi_init:'
         print*,'imobi-1=',imobi-1,' not equal to ntnpzd', ntnpzd 
         stop
      endif
      return
      end
!     END mobi_init

      subroutine setimobi (name, index, inc)
      implicit none
      character(*) :: name
      integer index, inc
      index = inc
      inc = index + 1
      print*,name,' = ',index
      return
      end

      subroutine mobi_driver(
     &                 kmx, twodt, rctheta, dayfrac, swr, tnpzd
     &,                t_in
# if defined O_mobi_o2
     &,                o2_in
# endif
# if defined O_carbon_13 || defined O_mobi_caco3
     &,                s_in, dic_in, alk_in, co2_in
# endif
# if defined O_carbon_13      
c     &,                dic13_in
# endif
# if defined O_mobi_silicon
c     &,                sil_in
# endif
     &,                sgb_in
# if defined O_mobi_nitrogen
c     &,                no3_in
#  if defined O_mobi_nitrogen_15
c     &,                din15_in
#  endif
# endif
# if defined O_sed
     &,                sedcorgflx, sedcalflx
# endif
     &,                src)

!     main driver for Model of Ocean Biogeochemistry and Isotopes (MOBI)
!     input: kmx, twodt, rctheta, dayfrac, swr, tnpzd, t_in, 
!            po4_in, (o2_in, s_in, dic_in, alk_in, co2_in, dic13_in, no3_in,
!            sgb_in, din15_in)
!     output: src, (sedcorgflx, sedcalflx)
!             if O_save_mobi_fluxes is switched on additional output (rnpp etc.)
!                is in mobi.h
!     copied from tracer.F (Andreas Schmittner, Oct 2, 2013)
!     modified to a vertical column model 
      implicit none
      integer k,kmx
#include "size.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "mobi.h"

      real snpzd(ntnpzd), tnpzd(km,ntnpzd), src(km,nsrc)
      real dic_npzd_sms(km), twodt, rctheta, gl, impo, expo
      real npp, remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real graz_Det, graz_Z, phin, dz, prca, dprca, nud, bct, bctz
      real t_in(km), po4_in(km), sedrr, nudop, nudon, npp_dop
      real impofe, expofe, remife, nfix(km)
      real fesed, bfe, o2_in(km)
# if defined O_save_mobi_diagnostics
      real avej, avej_D, gmax, no3P, po4P, po4_D
# endif
# if defined O_mobi_o2
      real fo2, so2, o2flag
# endif
# if defined O_mobi_iron
#  if defined O_save_mobi_diagnostics
      real feorgads, thetamax, deffe, feprime, fecol
      real deffe_C, thetamax_C
#  endif      
# endif
# if defined O_kk_ballast
      real expo_B, impo_B, remi_B
      real graz_Det_B, g_Det_B, biod_B
      real graz_Z_B, graz_C_B, morp_C_B, morz_B
# endif
      real calpro
# if defined O_mobi_caco3
      real dissl, caco3in, atmpres1, pHlo1, pHhi1, pH1, p_in
      real co2star1, dco2star1, pCO21, dpco21, CO31, omegaca, omegaar
      real c_in, wwc, expocaco3, impocaco3, dissk1, calatt
# endif
# if defined O_carbon_13 || defined O_mobi_caco3
      real s_in(km), dic_in(km), alk_in(km), co2_in
      real rc13impo, rc13expo, ac13_DIC_aq
      real ac13_aq_POC, ac13b, prca13, rdic13, rtdic13(km)
      real pH, co2star, dco2star, pCO2, dpco2, CO3
      real Omega_c, Omega_a, pCO2_old, pCO2_new
      real atmpres, depth
#  if defined O_mobi_caco3
      real rcaco3c13impo, rcaco3c13expo, rtcaco3c13(km)
#  endif
# endif
# if defined O_mobi_silicon
      real npp_Diat_dop
      real npp_Diat, graz_Diat, morp_Diat, morpt_Diat,felimit_Diat
#   if defined O_mobi_iron
      real deffe_Diat, thetamax_Diat
#   endif
# endif
# if defined O_mobi_silicon
      real wwo, impoopl, expoopl
      real settle, constant, temp_2(km)
      real bctd(km), opl_disk1, opl_dis, opl_pro
# endif
      real sgb_in(km)
# if defined O_mobi_nitrogen
      real lno3, sg_bdeni
      real morp_D, lntp, npp_D_dop
      real npp_D, graz_D, morpt_D, no3flag, wcdeni, bdeni(km)
      real din15_in(km), din15flag, eps_bdeni
      real rn15impo, rn15expo, uno3, rno3, bwcdeni, bbdeni
# endif   
# if defined O_sed
      real sedcorgflx, sedcalflx
      sedcalflx = 0.0
      sedcorgflx = 0.0
# endif

      expo = 0.0
      impo = 0.0
      phin = 0.0                ! integrated phytoplankton
      prca = 0.0                ! integrated production of calcite

      sedrr = 0.0
# if defined O_save_mobi_fluxes
      rsedrr = 0.0
      rprca = 0.0
# endif
# if defined O_mobi_nitrogen_15
      rn15impo = 0.0
      rn15expo = 0.0
# endif
# if defined O_carbon_13
      rc13impo = 0.0
      rc13expo = 0.0
      prca13 = 0.0
#  if defined O_mobi_caco3
      rcaco3c13impo = 0.0
      rcaco3c13expo = 0.0
#  endif
# endif

# if defined O_mobi_iron
      expofe = 0.0
      impofe = 0.0
# endif
# if defined O_kk_ballast
      expo_B = 0.0              ! detritus from coccs and zoop
      impo_B = 0.0              ! detritus from coccs and zoop
# endif
      calpro = 0.0
# if defined O_mobi_caco3
      caco3in = 0.0
      calatt = 0.0
      dissl = 0.0
      impocaco3 = 0.0
      expocaco3 = 0.0
      dissk1 = 0.0
# endif
# if defined O_mobi_silicon
      expoopl = 0.0             ! export of opal
      impoopl = 0.0             ! import of opal 
      opl_dis = 0.0             ! opal dissolution rate 
      opl_pro = 0.0             ! opal production rate 
      bctd(:) = 0.0             ! temperature scaling of microbial dissolution
# endif  

      src(:,:) = 0.0
      snpzd(:) = 0.0
      rcalpro(:) = 0.0
# if defined O_mobi_caco3
      rcalatt(:) = 0.0
      rdissl(:) = 0.0
      rexpocaco3(:) = 0.0
      romega_c(:) = 0.0
# endif
# if defined O_save_mobi_fluxes
      rexpo(:) = 0.0
      rgraz(:) = 0.0
      rgraz_Det(:) = 0.0
      rgraz_Z(:) = 0.0
      rmorp(:) = 0.0
      rmorz(:) = 0.0
      rnpp(:) = 0.0
      rmorpt(:) = 0.0
      rremi(:) = 0.0
      rexcr(:) = 0.0
#  if defined O_mobi_nitrogen
      rnpp_dop(:) = 0.0
      rnpp_D(:) = 0.0
      rnpp_D_dop(:) = 0.0
      rgraz_D(:) = 0.0
      rmorp_D(:) = 0.0
      rmorpt_D(:) = 0.0
      rnfix(:) = 0.0
#   if defined O_mobi_silicon
      rnpp_Diat_dop(:) = 0.0
      rdisopl(:) = 0.0
      rexpoopl(:) = 0.0
#   endif
#  endif
#  if defined O_kk_ballast
      rexpo_B(:) = 0.0          !ballast detritus
      rgraz_Det_B(:) = 0.0
      rremi_B(:) = 0.0
#  endif
#  if defined O_mobi_silicon
      rnpp_Diat(:) = 0.0
      rgraz_Diat(:) = 0.0
      rmorp_Diat(:) = 0.0
      rmorpt_Diat(:) = 0.0
#  endif
#  if defined O_mobi_iron
      rexpofe(:) = 0.0
      rremife(:) = 0.0
#   if defined O_save_mobi_diagnostics
      rfeorgads(:) = 0.0
      rdeffe(:) = 0.0
      rfeprime(:) = 0.0
      rfecol(:) = 0.0
#    if defined O_mobi_silicon
      rdeffe_Diat(:) = 0.0
#    endif
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
      ravej(:) = 0.0
      ravej_D(:) = 0.0
      rgmax(:) = 0.0
      rno3P(:) = 0.0
      rpo4P(:) = 0.0
      rpo4_D(:) = 0.0
#  endif
# endif

!1111 start of main k-loop
      do k=1,kmx
# if defined O_mobi_nitrogen_15
         rn15impo = rn15expo
# endif
# if defined O_carbon_13 || defined O_mobi_caco3
         atmpres = 1.0          !atm
         depth = zt(k)/100.     !m
CSPKKK still need to implement passing in sil_in when O_mobi_silicon is defined
         call co2calc_SWS (t_in(k), s_in(k), dic_in(k), alk_in(k)
     &,                    co2_in, atmpres, depth, pH, co2star, dco2star
     &,                    pCO2, dpco2
     &,                    CO3, Omega_c, Omega_a)
#  if defined O_carbon_13
!        c13 biological fractionation 
         ac13_DIC_aq = -1.0512994e-4*t_in(k)+1.011765
!        Popp et al. (1989) Am. J. Sci.
         ac13_aq_POC = -0.017*log10(min(max(co2star*1000.,2.),74.))
     &                 +1.0034
         ac13b = ac13_aq_POC/ac13_DIC_aq
         rc13impo = rc13expo*dztr(k)
#   if defined O_mobi_caco3         
         rcaco3c13impo = rcaco3c13expo*dztr(k)
#   endif
#  endif
#  if defined O_mobi_caco3
c         dissk0 =(1-(CO3-CO3/Omega_c))/(kcal+(CO3-CO3/Omega_C))
c         dissk1 = min(1.,max(0.,dissk0)) ! Andreas
         dissk1 = dissk0*max(0.,(1.-Omega_c)) ! Andreas 11 d^-1 from Gehlen 2006 
!AAA add omega dependence of calcite production
         capr = caprmax*max(0.,(Omega_c - 1.)/(kcapr + Omega_c - 1.)) ! Gehlen et al. (2007, BGS) eq. (3)
#  endif
#  if defined O_mobi_silicon
         opl_disk1 = opl_disk0 ! Andreas constant dissolution rate
#  endif
# endif
         swr = swr*exp(-kc*phin
# if defined O_mobi_caco3
     &             - kc_c*caco3in
# endif
cAAA add opal
     &              )         

         phin = max(tnpzd(k,imobiphyt), trcmin)*dzt(k)
# if defined O_mobi_nitrogen
     &        + max(tnpzd(k,imobidiaz), trcmin)*dzt(k)
# endif
# if defined O_mobi_silicon
     &        + max(tnpzd(k,imobidiat), trcmin)*dzt(k)
# endif
# if defined O_mobi_caco3     
         caco3in = caco3in + tnpzd(k,imobicaco3)*dzt(k)         
         impocaco3 = expocaco3*dztr(k)
# endif
# if defined O_kk_ballast
         impo_B = expo_B*dztr(k) !ballasted detritus
# endif
         gl = swr*exp(ztt(k)*rctheta)
         impo = expo*dztr(k)
# if defined O_mobi_iron
         impofe = expofe*dztr(k)
# endif

         bct = bbio**(cbio*t_in(k))
#  if defined O_mobi_silicon
         bctd(k) = bct          !store it for using below
         impoopl = expoopl*dztr(k)
#  endif
# if defined O_mobi_o2
         bctz = (0.5*(tanh(o2_in(k) - 8.)+1))
c     &        *bbio**(cbio*min(t_in(k),20.))
cAAA remove limit
     &        *bbio**(cbio*t_in(k))
!        decrease remineralisation rate in oxygen minimum zone
         nud = nud0*(0.65+0.35*tanh(o2_in(k) - 3.0))         
# else
         bctz = bbio**(cbio*min(t_in(k),20.))
         nud = nud0         
# endif
# if defined O_mobi_nitrogen         
         nudon = nudon0
         nudop = nudop0
# endif
!-----------------------------------------------------------------------
!        call the ecosystem dynamics model
!-----------------------------------------------------------------------

         call mobi_src (tnpzd(k,:), gl, bct, impo, dzt(k)
     &,                 dayfrac, wd(k), nud
# if defined O_mobi_caco3
     &,                 impocaco3, wc(k), dissk1, expocaco3
# endif        
# if defined O_kk_ballast
     &,                 impo_B, expo_B, remi_B, graz_Det_B
# endif
# if defined O_mobi_silicon
     &,                 impoopl, wo(k), opl_disk1, expoopl
# endif
# if defined O_mobi_nitrogen
     &,                 nudop, nudon, nfix(k)
# endif
     &,                 snpzd, expo
# if defined O_carbon
     &,                 calpro
#  if defined O_mobi_caco3
     &,                 dissl, calatt
#  endif
# endif
# if defined O_save_mobi_fluxes
     &,                 npp, morpt, remi, excr, graz, morp, morz
     &,                 graz_Det, graz_Z
#  if defined O_mobi_silicon
     &,                 npp_Diat, morpt_Diat, graz_Diat, morp_Diat
     &,                 opl_dis, opl_pro
#  endif
#  if defined O_mobi_nitrogen
     &,                 npp_D, npp_D_dop, graz_D, morp_D, morpt_D
     &,                 npp_dop
#   if defined O_mobi_silicon
     &,                 npp_Diat_dop
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
     &,                 avej, avej_D, gmax, no3P
     &,                 po4P, po4_D
#  endif
# endif
     &,                 bctz
# if defined O_mobi_nitrogen_15
     &,                 rn15impo, rn15expo
# endif
# if defined O_carbon_13
     &,                 rc13impo, rc13expo, ac13b
#  if defined O_mobi_caco3
     &,                 rcaco3c13impo, rcaco3c13expo
#  endif
# endif
# if defined O_mobi_iron
     &,                 expofe, impofe, remife
#  if defined O_save_mobi_diagnostics
     &,                 feorgads, thetamax, deffe
     &,                 feprime, fecol
#   if defined O_mobi_silicon
     &,                 deffe_Diat, thetamax_Diat
#   endif         
#  endif
# endif         
# if defined O_mobi_o2
     &,                 o2_in(k)            
# endif
     &                      )
! These are source/sink terms 
         snpzd = snpzd*rdtts
# if defined O_mobi_iron
         expofe = expofe*rnbio
# endif
# if defined O_mobi_caco3
         expocaco3 = expocaco3*rnbio
# endif
# if defined O_mobi_silicon
         expoopl = expoopl*rnbio
# endif
# if defined O_kk_ballast
         expo_B = expo_B*rnbio !ballast detritus
# endif
         expo = expo*rnbio
# if defined O_mobi_nitrogen_15
         rn15expo = rn15expo*rnbio
# endif
# if defined O_carbon_13
         rc13expo = rc13expo*rnbio
#  if defined O_mobi_caco3
         rcaco3c13expo = rcaco3c13expo*rnbio
#  endif
# endif
         rcalpro(k) = calpro*rnbio
# if defined O_mobi_caco3
         rcalatt(k) = calatt*rnbio
         rdissl(k) = dissl*rnbio
         rexpocaco3(k) = expocaco3
         romega_c(k) = Omega_c
# endif
# if defined O_save_mobi_fluxes
         rexpo(k) = expo
         rgraz(k) = graz*rnbio 
         rgraz_Det(k) = graz_Det*rnbio
         rgraz_Z(k) = graz_Z*rnbio
         rmorp(k) = morp*rnbio
         rmorz(k) = morz*rnbio
         rnpp(k) = npp*rnbio
         rmorpt(k) = morpt*rnbio
         rremi(k) = remi*rnbio
         rexcr(k) = excr*rnbio
#  if defined O_mobi_nitrogen
         rnpp_dop(k) = npp_dop*rnbio         
         rnpp_D(k) = npp_D*rnbio
         rnpp_D_dop(k) = npp_D_dop*rnbio
         rgraz_D(k) = graz_D*rnbio
         rmorp_D(k) = morp_D*rnbio
         rmorpt_D(k) = morpt_D*rnbio
         rnfix(k) = nfix(k)*rnbio
#   if defined O_mobi_silicon
         rdisopl(k) = opl_dis*rnbio
         rproopl(k) = opl_pro*rnbio
         rexpoopl(k) = expoopl
         rnpp_Diat_dop(k) = npp_Diat_dop*rnbio
#   endif
#  endif
#  if defined O_kk_ballast
         rexpo_B(k) = expo_B  !ballast detritus
         rgraz_Det_B(k) = graz_Det_B*rnbio
         rremi_B(k) = remi_B*rnbio
#  endif
#  if defined O_mobi_silicon
         rnpp_Diat(k) = npp_Diat*rnbio
         rgraz_Diat(k) = graz_Diat*rnbio
         rmorp_Diat(k) = morp_Diat*rnbio
         rmorpt_Diat(k) = morpt_Diat*rnbio
#  endif
#  if defined O_mobi_iron
         rexpofe(k) = expofe
         rremife(k) = remife*rnbio
#   if defined O_save_mobi_diagnostics
         rfeorgads(k) = feorgads*rnbio
         rdeffe(k) = deffe*rnbio
         rfeprime(k) = feprime*rnbio
         rfecol(k) = fecol*rnbio
#    if defined O_mobi_silicon
         rdeffe_Diat(k) = deffe_Diat*rnbio
#    endif         
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
         ravej(k) = avej*rnbio
         ravej_D(k) = avej_D*rnbio
         rgmax(k) = gmax*rnbio
         rno3P(k) = no3P*rnbio
         rpo4P(k) = po4P*rnbio
         rpo4_D(k) = po4_D*rnbio
#  endif
# endif
# if defined O_mobi_caco3
!----------------------------------------------------------------------
!             calculate calcite at the bottom and dissolve
!---------------------------------------------------------------------
#  if defined O_sed
c                if (addflxo .and. eots)
c     &           sbc(i,jrow,ircal) = sbc(i,jrow,ircal)
c     &           + expocaco3*twodt(1)*dzt(k)*1.e-3 !mol C
#  endif
# endif

!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

# if defined O_kk_ballast
         sedrr = sgb_in(k)*(expo+expo_B)*dzt(k)
# else               
         sedrr = sgb_in(k)*expo*dzt(k)
# endif
# if defined O_mobi_nitrogen
!-----------------------------------------------------------------------
!        benthic denitrification model of Bohlen et al., 2012, GBC (eq. 10)
!        NO3 is removed out of bottom water nitrate.
!        See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!        limit denitrification as nitrate approaches 0 uM
c         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
         no3flag = 0.5+sign(0.5,tnpzd(k,imobino3)-trcmin)
#  if defined O_mobi_nitrogen_15
c         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
         din15flag = 0.5+sign(0.5,tnpzd(k,imobidin15)-trcmin)
#  else
         din15flag = 1.
#  endif
c         lno3 = 0.5*tanh(no3_in(k)*10 - 5.0)
         lno3 = 0.5*tanh(tnpzd(k,imobino3)*10 - 5.0)

         sg_bdeni = (0.06 + 0.19*0.99**(max(o2_in(k),trcmin)
     &                 - max(tnpzd(k,imobino3),trcmin)))
#  if defined O_kk_ballast
     &              *max((expo+expo_B)*sgb_in(k),trcmin)*redctn*1.e3
         sg_bdeni = min(sg_bdeni, sgb_in(k)*(expo+expo_B))
#  else               
     &              *max(expo*sgb_in(k),trcmin)*redctn*1.e3
         sg_bdeni = min(sg_bdeni, sgb_in(k)*expo)
#  endif
         sg_bdeni = max(sg_bdeni, 0.)

c Andreas is this line correct? Why 2000 m?
c         if (zt(k) .le. 200000.) sg_bdeni = sgbdfac*sg_bdeni

         sg_bdeni = sg_bdeni*(0.5 + lno3)*no3flag*din15flag
         bdeni(k) = sg_bdeni
         snpzd(imobino3) = snpzd(imobino3) + sgb_in(k)*expo
     &                                           - sg_bdeni
#  if defined O_save_mobi_fluxes
         rbdeni(k) = bdeni(k)
#  endif
#  if defined O_mobi_nitrogen_15
         rno3 = max(tnpzd(k,imobidin15), trcmin*rn15std/(1+rn15std))
     &            /max(tnpzd(k,imobino3)-tnpzd(k,imobidin15),
     &                 trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         eps_bdeni = eps_bdeni0*exp(-2.5e-6*(zt(k)))
         bbdeni = rno3 - eps_bdeni*rno3/1000.

         snpzd(imobidin15) = snpzd(imobidin15) 
     &                     + rn15expo*sgb_in(k)*expo 
     &                     - bbdeni/(1+bbdeni)*sg_bdeni
#  endif              
# endif
# if defined O_mobi_iron
         fesed=fetopsed*bct*redptn*expo*sgb_in(k)
         bfe=fesed
         snpzd(imobidfe) = snpzd(imobidfe) + fesed                
         o2flag = 0.5 - sign(0.5, o2_in(k)-o2min)
         snpzd(imobidfe) = snpzd(imobidfe) + expofe*sgb_in(k)*o2flag
         bfe = bfe + expofe*sgb_in(k)*o2flag
         expofe = expofe - sgb_in(k)*expofe*o2flag
#  if defined O_save_mobi_fluxes
         rremife(k) = rremife(k) + bfe
#  endif
#  if defined O_save_mobi_diagnostics
         rfesed=fesed
         rbfe=bfe 
#  endif
# endif
         snpzd(imobipo4) = snpzd(imobipo4) + sgb_in(k)*expo*redptn
# if defined O_carbon
         snpzd(imobidic) = snpzd(imobidic) + sgb_in(k)*expo*redctn
# endif
# if defined O_carbon_13
         snpzd(imobidic13) = snpzd(imobidic13) 
     &                     + rc13expo*sgb_in(k)*redctn
         rc13expo = rc13expo - sgb_in(k)*rc13expo
# endif
         expo = expo - sgb_in(k)*expo
# if defined O_save_mobi_fluxes
#  if defined O_kk_ballast
         expo_B = expo_B - sgb_in(k)*expo_B               
         rremi(k) = rremi(k) + sgb_in(k)*(expo+expo_B)
#  else
         rremi(k) = rremi(k) + sgb_in(k)*expo
#  endif               
         rsedrr = rsedrr + sedrr
# endif

!-----------------------------------------------------------------------
!             set source/sink terms
!-----------------------------------------------------------------------

         src(k,ispo4) = snpzd(imobipo4)
         src(k,isphyt) = snpzd(imobiphyt)
         src(k,iszoop) = snpzd(imobizoop)
         src(k,isdetr) = snpzd(imobidetr)
# if defined O_carbon         
         src(k,isdic) = snpzd(imobidic)
# endif         
# if defined O_mobi_nitrogen
         src(k,isdop) = snpzd(imobidop)   
         src(k,isno3) = snpzd(imobino3)
         src(k,isdon) = snpzd(imobidon)
         src(k,isdiaz) = snpzd(imobidiaz)
#  if defined O_mobi_nitrogen_15
         src(k,isdin15) = snpzd(imobidin15)
         src(k,isdon15) = snpzd(imobidon15)
         src(k,isphytn15) = snpzd(imobiphytn15)
         src(k,iszoopn15) = snpzd(imobizoopn15)
         src(k,isdetrn15) = snpzd(imobidetrn15)
         src(k,isdiazn15) = snpzd(imobidiazn15)
#   if defined O_mobi_silicon
         src(k,isdiatn15) = snpzd(imobidiatn15)
#   endif         
#  endif
# endif
# if defined O_carbon_13
         src(k,isdic13) = snpzd(imobidic13)
         src(k,isphytc13) = snpzd(imobiphytc13)
         src(k,iszoopc13) = snpzd(imobizoopc13)
         src(k,isdetrc13) = snpzd(imobidetrc13)
#  if defined O_mobi_nitrogen
         src(k,isdoc13) = snpzd(imobidoc13)
         src(k,isdiazc13) = snpzd(imobidiazc13)
#  endif         
#  if defined O_mobi_silicon
         src(k,isdiatc13) = snpzd(imobidiatc13)
#  endif
#  if defined O_mobi_caco3
         src(k,iscaco3c13) = snpzd(imobicaco3c13)
#  endif
# endif
# if defined O_mobi_caco3
         src(k,iscaco3) = snpzd(imobicaco3)
# endif
# if defined O_kk_ballast
         src(k,isdetr_B) = snpzd(imobidetr_B)
# endif
#  if defined O_mobi_silicon
         src(k,isdiat) = snpzd(imobidiat)
         src(k,issil) = snpzd(imobisil)
         src(k,isopl) = snpzd(imobiopl)
#  endif
# if defined O_mobi_iron
         src(k,isdfe) = snpzd(imobidfe)
         src(k,isdetrfe) = snpzd(imobidetrfe)
# endif

!-----------------------------------------------------------------------
!             production, export, and dissolution of opal
!-----------------------------------------------------------------------
# if defined O_mobi_silicon
c         sipr = 0.
#  if defined O_kk_variable_sipr && defined O_mobi_iron
!     From HAMOCC5 (Aumont et al, GBC, 2003) values 0-0.5, mostly ~0.1
c         sipr = sipr0 * min(1.,sil_in(k)/k1si)
c     &                  * (4.-3.*min(1.,tnpzd(k,imobidfe)/kfe_Diat))
#  else
c         sipr = sipr0
#  endif
!        production of opal. 
c         dprop = roplpro(k)*1e-3
c         dprop = (morp_Diat + graz_Diat*(1.-gamma1))
c     &              *sipr*redctn*rnbio !umol Si cm-3 s-1
c         prop = prop + dprop*dzt(k) !umol Si cm-2 s-1
CSPKKK I've moved the hydrothermal and silwflux stuff in kk's version to tracer.F
c         src(k,issil) = -dprop !umol Si cm-3 s-1
# endif

!        These are organic sources and sinks of DIC (i.e. remin - pp)
!        all are based on po4 uptake and remineralization
# if !defined O_carbon
         dic_npzd_sms(k) = snpzd(imobipo4)*redctp
# else
         dic_npzd_sms(k) = snpzd(imobidic)
!        production of calcite
         dprca = rcalpro(k)*1e-3
         prca = prca + dprca*dzt(k)

#  if !defined O_mobi_caco3         
         src(k,isdic) = snpzd(imobidic) - dprca             
#  endif
#  if defined O_carbon_13
c         rtdic13(k) = max(dic13_in(k),trcmin*rc13std/(1+rc13std))
         rtdic13(k) = max(tnpzd(k,imobidic13)
     &            ,trcmin*rc13std/(1+rc13std)) / max(dic_in(k),trcmin)
         rtdic13(k) = min(rtdic13(k), 2.*rc13std/(1+rc13std))
         rtdic13(k) = max(rtdic13(k), 0.5*rc13std/(1+rc13std))

         prca13 = prca13 + dprca*dzt(k)*rtdic13(k)
#   if defined O_mobi_caco3          
         rtcaco3c13(k) = max(tnpzd(k,imobicaco3c13)
     &                     ,trcmin*rc13std/(1+rc13std))
     &             / max(tnpzd(k,imobicaco3),trcmin)
         rtcaco3c13(k) = min(rtcaco3c13(k), 2.*rc13std/(1+rc13std))
         rtcaco3c13(k) = max(rtcaco3c13(k), 0.5*rc13std/(1+rc13std))
#   endif
#   if !defined O_mobi_caco3          
         src(k,isdic13) = src(k,isdic13) - rtdic13(k)*dprca
#   endif
#  endif
# endif       
# if defined O_mobi_alk
         src(k,isalk) = -snpzd(imobidic)*redntc*1.e-3
#  if !defined O_mobi_caco3
     &        - 2.*dprca
#  endif
# endif
!        calculate total export to get total import for next layer
         expo = expo*dzt(k)
# if defined O_carbon_13
         rc13expo = rc13expo*dzt(k)
#  if defined O_mobi_caco3
         rcaco3c13expo = rcaco3c13expo*dzt(k)
#  endif         
# endif         
# if defined O_kk_ballast
         expo_B = expo_B*dzt(k)
# endif         
# if defined O_mobi_iron
         expofe = expofe*dzt(k)
# endif
# if defined O_mobi_caco3
         expocaco3 = expocaco3*dzt(k)
# endif
# if defined O_mobi_silicon
         expoopl = expoopl*dzt(k)
# endif
      enddo
!1111 end of main k-loop
# if defined O_save_mobi_fluxes
      rprca = prca
# endif
# if defined O_sed
      sedcorgflx = sedrr*redctn*twodt
#  if defined O_mobi_caco3
      sedcalflx = expocaco3*1.e-3*twodt
#  endif
# endif

# if defined O_mobi_o2
!2222 start of second k-loop
      do k=1,kmx
!        limit oxygen consumption below concentrations of
!        5umol/kg as recommended in OCMIP
         fo2 = 0.5*tanh(o2_in(k) - 2.5)
!        sink of oxygen
         so2 = dic_npzd_sms(k)*redotc
#  if defined O_mobi_nitrogen
! O2 is needed to generate the equivalent of NO3 from N2 during N2 fixation
! 0.5 H2O + 0.5 N2+1.25O2 -> HNO3
! note that so2 is -dO2/dt
     &          + nfix(k)*rnbio*1.25e-3
!        add denitrification as source term for NO3
         no3flag = 0.5+sign(0.5,tnpzd(k,imobino3)-trcmin)
#   if defined O_mobi_nitrogen_15
         din15flag = 0.5+sign(0.5,tnpzd(k,imobidin15)-trcmin)
#   else
         din15flag = 1.
#   endif
         lno3 = 0.5*tanh(tnpzd(k,imobino3) - 2.5)
         lntp = 0.5*tanh(tnpzd(k,imobino3)/
     &                   (redntp*tnpzd(k,imobipo4))*100. - 60.)
!        800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
         wcdeni = 800.*no3flag*so2*(0.5 - fo2)*(0.5 + lno3)
#   if defined O_mobi_nitrogen_15
     &           *din15flag
#   endif
         wcdeni = max(wcdeni, 0.)
         src(k,isno3) = src(k,isno3) - wcdeni
#   if defined O_mobi_nitrogen_15
!        calculate isotope effect of water column denitrification
         uno3 = wcdeni*twodt/tnpzd(k,imobino3)
         uno3 = min(uno3, 0.999)
         uno3 = max(uno3, trcmin)
         rno3 = max(tnpzd(k,imobidin15),trcmin*rn15std/(1+rn15std))
     &             / max(tnpzd(k,imobino3)-tnpzd(k,imobidin15),
     &                   trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         bwcdeni = rno3 + eps_wcdeni*(1-uno3)/uno3*log(1-uno3)*rno3
     &                    /1000.

         src(k,isdin15) = src(k,isdin15) - (bwcdeni/(1+bwcdeni))*wcdeni
#   endif
#   if defined O_mobi_alk
! Correct the ALK stoichiometry to account for N2 fixation
         src(k,isalk) = src(k,isalk) + wcdeni*1.e-3
!bdeni
         src(k,isalk) = src(k,isalk) + bdeni(k)*1.e-3
! Now add account for N2 fixation (ALK production is tied to PO4 change and
! thus in the case of N2 fixation is not correct).
         src(k,isalk) = src(k,isalk) - nfix(k)*rnbio*1.e-3
#   endif
#   if defined O_save_mobi_fluxes
         rwcdeni(k) = wcdeni
#   endif
#  endif
         src(k,iso2) = -so2*(0.5 + fo2)
      enddo
!2222 end of second k-loop
# endif

!-----------------------------------------------------------------------
!           remineralize calcite and opal
!-----------------------------------------------------------------------
!3333 start of third k-loop
      do k=1,kmx-1
# if defined O_carbon
#  if defined O_mobi_caco3
         src(k,isdic) = src(k,isdic)
     &        + rdissl(k)*1.e-3
     &        - rcalpro(k)*1.e-3
#  else
         src(k,isdic) = src(k,isdic) + prca*rcak(k)
#  endif
#  if defined O_carbon_13
#   if defined O_mobi_caco3
         src(k,isdic13) = src(k,isdic13)
     &        + rdissl(k)*1.e-3*rtcaco3c13(k) ! AS
     &        - rcalpro(k)*1.e-3*rtdic13(k)
#   else              
         src(k,isdic13) = src(k,isdic13) + prca13*rcak(k)
#   endif
#  endif
# endif
# if defined O_mobi_alk
#  if defined O_mobi_caco3
         src(k,isalk) = src(k,isalk)
     &        + 2.*rdissl(k)*1.e-3
     &        - 2.*rcalpro(k)*1.e-3
#  else
         src(k,isalk) = src(k,isalk) + 2.*prca*rcak(k)
#  endif
# endif
      enddo
!3333 end of third k-loop
# if defined O_carbon
#  if defined O_mobi_caco3
      src(kmx,isdic) = src(kmx,isdic)
     &     + rdissl(kmx)*1.e-3
     &     - rcalpro(kmx)*1.e-3
     &     + rexpocaco3(kmx)*1.e-3
#  else
#   if defined O_sed
      sedcalflx = (prca*rcab(kmx)
     &     - prca*rcak(kmx))*twodt*dzt(kmx)
#   endif
      src(kmx,isdic) = src(kmx,isdic) + prca*rcab(kmx)
#  endif
#  if defined O_carbon_13
#   if defined O_mobi_caco3
      src(kmx,isdic13) = src(kmx,isdic13)
     &                 + rdissl(kmx)*1.e-3*rtcaco3c13(kmx)
     &                 - rcalpro(kmx)*1.e-3*rtdic13(kmx)
     &                 + rexpocaco3(kmx)*1.e-3*rtcaco3c13(kmx)       
#   else  
      src(kmx,isdic13) = src(kmx,isdic13) + prca13*rcab(kmx)
#   endif
#  endif
# endif
# if defined O_mobi_alk
#  if defined O_mobi_caco3
      src(kmx,isalk) = src(kmx,isalk)
     &               + 2.*rdissl(kmx)*1.e-3
     &               - 2.*rcalpro(kmx)*1.e-3
     &               + 2.*rexpocaco3(kmx)*1.e-3
#  else
      src(kmx,isalk) = src(kmx,isalk) + 2.*prca*rcab(kmx)
#  endif
# endif

CSPKKK one more k loop
#  if defined O_mobi_silicon
!  Test some dissolution parameters, a default exponential from Aumont et al. (2003)
!  can be found in setmom.F    
c      expo_op = prop  !flux

cAAA opal dissolution now calculated in npzd_src
c      do k=1,kmx-1
c       temp_2(k) = t_in(k)+273.15 !convert to deg K 
!    (1) From Gnanadesikan 1999 !lower rates than HAMOCC
!     might over-estimate dissl in deep ocean (discussed in Ridgwell et al 2002)       
!       ropk(k) = dzt(k)*si_dis/ws0*exp(-11481./temp_2(k)) !unitless theoretical dissl
!     (2) From PlankTOM10, includes oxygen. Similar to HAMOCC except it always takes the minimum=1 and you end up with high dissolution.Plus there is no evidence of an oxygen dependency! 
!       ropk(k) = min(1.,si_dis*exp(-11200./temp_2(k)))
 !    &           *avail_oxy(k)/ws0*dzt(k)      
!     (3) Let's include the microbial loop instead. I made this up.
c       ropk(k) = dzt(k)*si_dis0/wo0*bctd(k) !exp(-11481./temp_2(k))
c       si_dis = si_dis0
!      if(rmorpt_Diat(i,k,j).gt.0)then
!       ropk(k) = ropk(k)+ropk(k)*bctd(i,k,j)!*0.5!0.1     
!       ropk(k)=ropk(k)*bctd(i,k,j)
!      endif
c       if(ropk(k).gt.1) ropk(k) = 1.
c       op_diss = expo_op*ropk(k) !umol Si cm-2 s-1
c       expo_op = expo_op - op_diss
c       src(k,issil) = src(k,issil) + op_diss*dztr(k) !umol Si cm-3 s-1
CSPKKK I deleted kk's not O_mobi_subgridbathy block
c       if (sgb_in(k) .gt. 0) then
! Use a limit (2 mmol Si /m^2/day) and ratios from Sarmiento & Gruber (2006)
! Buried Si is lost from the system!!!
c          if (expo_op .le. 0.2/24/60/60) then
c             bsi(k) = 0.05*expo_op*sgb_in(k) !buried sil umol/cm^2/s
c          else
c             bsi(k) = 0.3*expo_op*sgb_in(k)
c          endif
c          expo_op = expo_op - bsi(k)
CSPKKK I've moved the silwflux stuff in kk's version to tracer.F		 
c       endif
c      enddo ! K loop
c      src(kmx,issil) = src(kmx,issil) + expo_op*dztr(kmx) !put the leftovers back into ocean
      src(kmx,issil) = src(kmx,issil) + rexpoopl(kmx) !put the leftovers back into ocean
#  endif

      return
      end
!     END mobi_driver

      subroutine mobi_src (bioin, gl, bct, impo, dzt
     &,                    dayfrac, wwd, nud
# if defined O_mobi_caco3
     &,                    impocaco3, wwc, dissk1, expocaco3out
# endif
# if defined O_kk_ballast
     &,                    impo_B, expoout_B, remiout_B, graz_Det_out_B
# endif
# if defined O_mobi_silicon
     &,                    impoopl, wwo, opl_disk1, expooplout
# endif
# if defined O_mobi_nitrogen
     &,                    nudop, nudon, nfixout
# endif      
     &,                    bioout, expoout
# if defined O_carbon
     &,                    calproout
#  if defined O_mobi_caco3
     &,                    disslout, calattout
#  endif      
# endif      
# if defined O_save_mobi_fluxes
     &,                    nppout, morptout, remiout, excrout, grazout
     &,                    morpout, morzout, graz_Det_out, graz_Zout
#  if defined O_mobi_silicon
     &,                    npp_Diatout, morpt_Diatout, graz_Diatout
     &,                    morp_Diatout, opldisout, oplproout
#  endif      
#  if defined O_mobi_nitrogen
     &,                    npp_Dout, npp_D_dopout, graz_Dout, morp_Dout
     &,                    morpt_Dout, npp_dopout
#   if defined O_mobi_silicon
     &,                    npp_Diat_dopout
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
     &,                    avej_out, avej_D_out, gmax_out, no3P_out
     &,                    po4P_out, po4_D_out
#  endif
# endif
     &,                    bctz
# if defined O_mobi_nitrogen_15
     &,                    rn15impo, rn15expoout
# endif
# if defined O_carbon_13
     &,                    rc13impo, rc13expoout, ac13b
#  if defined O_mobi_caco3
     &,                    rcaco3c13impo, rcaco3c13expoout
#  endif
# endif
# if defined O_mobi_iron
     &,                    expofeout, impofe, remifeout
#  if defined O_save_mobi_diagnostics
     &,                    feorgadsout, thetamaxout, deffeout
     &,                    feprimeout, fecolout
#   if defined O_mobi_silicon
     &,                    deffeout_Diat, thetamaxout_Diat
#   endif  
#  endif
# endif       
# if defined O_mobi_o2
     &,                    o2
# endif
     &                     ) 

!=======================================================================
!     computes source terms of the NPZD model
!     initial version of code adapted from Xavier Giraud:
!     Giraud et al. 2000, J Mar Res, 58, 609-630
!     original model reference:
!     Oeschlies and Garcon 1999, Global Biogeochem. Cycles 13, 135-160
!     Schmittner et al. 2005,  Global Biogeochem. Cycles 19, GB3004,
!     doi:10.1029/2004GB002283.
!     Schmittner et al. 2008, Global Biogeochem. Cycles 22, GB1013
!
!     This version was modified by David Keller and corrects the zooplankton
!     grazing formulation.  Note that zooplankton are now allowed to graze
!     on themselves and detritus, in addition to phyt. and diazotrophs.
!     The calculation of light has also been corrected.
!     Keller et al. 2012, Geosci. Model Dev., 5, 1195-1220
!
!     The nitrogen isotope model (nitrogen_15) has been developed by 
!     Chris Somes and is described in
!     Somes et al. 2010, Global Biogeochem. Cycles 24, GB4019
!     Somes et al. 2010, Geophys. Res. Lett. 37, L23605 and
!     Somes et al. 2013, Biogeosc. 10, 5889-5910
!
!     The carbon isotope model (carbon_13) has been written by Andreas 
!     Schmittner and is documented in
!     Schmittner et al. 2013 Biogeosc. 10, 5793-5816
!     Chris Somes modified the code converting from alpha to beta formulation
!
!     Karen Kvale has included prognostic equations for CaCO3, coccolithphores, and
!     ballast described in
!     Kvale et al. 2015 Atmos.-Ocean 53, doi:10.1080/07055900.2015.1049112
!     Andreas has removed coccs
!
!     Levin Nickelsen has written the iron model:
!     Nickelsen et al. 2015, Geosc. Model Dev. 8, 1357-1381, doi:10.5194/gmd-8-1357-2015
!
!     Juan Muglia has modified the subgrid-scale scheme used for benthic denitrification
!     and iron release from sediments. He has also included geothermal sources of iron:
!     Muglia et al. 2017, doi:10.1002/2016PA003077.
!
!     Note that nutrient (N) represents phosphate
!
!     input variables:
!
!       bioin(1:ntnpzd) [mmol m-3]
!     
!       gl         = 2.*light at top of grid box
!       nbio       = number of time steps
!       dtbio        = time step [s]
!       bct        = bbio**(cbio*temperature)
!       impo       = import of detritus from above [mmol m-3]
!       dzt        = depth of grid box [cm]
!       dayfrac    = day length (fraction: 0 < dayfrac < 1)
!       wwd        = sinking speed of detritus/dzt
!       nud        = remineralisation rate of detritus [s-1]
!
!     output variables:
!
!       bioout     = change from bioin [mmol m-3]
!       nppout     = net primary production [mmol m-3]
!       grazout    = grazing [mmol m-3]
!       morpout    = quadratic mortality of phytoplankton [mmol m-3]
!       morptout   = specific mortality of phytoplankton [mmol m-3]
!       morzout    = mortality of zooplankton [mmol m-3]
!       remiout    = remineralisation [mmol m-3]
!       excrout    = excretion [mmol m-3]
!       expoout    = detrital export [mmol m-3]
!       npp_Dout   = NPP of diazotrophs
!       graz_Dout  = grazing of diazotrophs
!       morp_Dout  = mortality of diazotrophs
!       nfixout    = rate of N2 fixation
!       graz_Det_out = grazing of detritus
!       graz_Zout   = grazing on othe zooplankton
!       avej_out    = light-depend phyt. growth rate
!       avej_D_out  = light-depend Diaz growth rate
!       gmax_out    = temp-depend. zoo growth rate
!       no3P_out    = no3 depend. phyt growth rate
!       po4P_out    = po4 depend. phyt growth rate
!       po4_D_out   = po4 depend. Diaz growth rate
!
!       New grazing formulation variables and parameters
!
!       The following terms determine ingestion according to a
!       a Holling II curve (i.e. Michaelis Menten):
!
!       Ingestion = max_graz_rate * (Ft/(Ft + kzoo))
!
!       where Ft is the weighted measure of the total food available
!       and equals the sum of the different prey types times the
!       preference of Z for that type of prey 
!
!       zprefP   = Z preference for P
!       zprefDiaz   = Z preference for Diaz
!       zprefDet = Z preference for detritus
!       zprefZ   = Z preference for other Z
!       kzoo = half saturation coefficienct for Z ingestion mmol N m-3
!       ing_P    = zooplankton ingestion of phytoplankon
!       ing_D    = zooplankton ingestion of diazotrophs
!       ing_Det  = zooplankton ingestion of detritus
!       ing_Z    = zooplankton ingestion of other zooplankton
!       thetaZ   = Michaelis-Menten denominator
!
!=======================================================================

      implicit none

      integer n

      real gl, f1, biopo4, biophyt, biozoop, biodetr, jmax, u_P, g_P
      real morp, morpt, morz, remi, excr, expo, impo, nppout, grazout
      real morpout, morptout, morzout, remiout, excrout, expoout
      real avej_out, avej_D_out, gmax_out, no3P_out, po4P_out, po4_D_out
      real dzt, po4flag, phytflag, zoopflag, detrflag, wwd, gd
      real nupt, nud, biono3, u_D,npp_D, npp_Dout, no3flag, biodiaz, bct
      real diazflag, g_D,graz_D, morpt_D, jmax_D, gd_D, avej_D, no3upt_D
      real morpt_Dout, graz_Dout, nfixout, biop2, u1, u2, phi1, phi2
      real avej, graz_Det_out, graz_Zout, thetaZ, ing_P, ing_D, dayfrac
      real ing_Det, ing_Z, g_Z, g_Det, graz_Z, graz_Det, gmax
      real no3P, po4P, po4_D, bctz
      real excrdiaz, biodop, biodon, dopflag, donflag, recy_don
      real recy_dop, npp, graz, biodin15, biodon15, biophytn15
      real biodetrn15, biodiazn15, bassim, fcassim, bexcr, fcexcr
      real bnfix, fcnfix, rn15impo, rn15expoout, dig, dig_P, dig_Z
      real dig_Det, dig_D, excr_P, excr_Z, excr_Det, excr_D, sf, sf_P
      real sf_Z, sf_Det, sf_D, nr_excr_D, uno3, rno3, rzoop
      real rtdin15, rtphytn15, rtzoopn15, rtdetrn15, rtdiazn15, rtdon15
      real din15flag, biozoopn15, morp_Dout
      real don15flag, phytn15flag, zoopn15flag, detrn15flag, diazn15flag
      real limP, limP_D, dopupt_D_flag, dopupt_D, morp_D, nupt_D
      real limP_dop, limP_po4
      real udon, rdon, brecy, fcrecy, biodic, dopupt_flag, dopupt
      real biodic13, biophytc13, biozoopc13, biodetrc13, biodiazc13
      real dic13flag, phytc13flag, zoopc13flag, detrc13flag, diazc13flag
      real biodoc13, doc13flag, biodoc, rcaco3c13impo, rcaco3c13expoout
      real rdic13, rtphytc13, rtzoopc13, rtdetrc13, rtdiazc13, ac13b
      real rc13impo, rc13expoout, dicflag, bc13npp, rtdic13, fcnpp
      real rtdoc13, nudop, nudon, npp_dopout, npp_D_dopout

      real expo_B, impo_B, remi_B, remiout_B, expoout_B, dflag_B
      real graz_Det_out_B, graz_Det_B, g_Det_B, biod_B, dig_Det_B
      real graz_Z_B, graz_C_B, morp_C_B, morz_B, excr_Det_B, sf_Det_B

      real diatn15flag, biodiatn15, rtdiatn15
      real biodiatc13, diatc13flag, rtdiatc13
      real biocaco3c13, caco3c13flag, rtcaco3c13

      real caco3flag, wwc, expocaco3
      real expocaco3out,calpro,calproout
      real biocaco3, impocaco3, dissk1, dissl, disslout
      real calatt,calattout
      
      real biodfe, biodetrfe, dfeflag, detrfeflag, o2flag
      real expofe, impofe, feorgads, remife, thetamax, deffe, fepa
      real thetamax_D, deffe_D, thetachl, thetachl_D, chl, feprime
      real fecol, kirr, fesed, gl_O, gl_D, alpha_O
      real alpha_D, deffe_C, thetamax_C, thetamaxout_C, deffeout_C
      real alpha_C

      real expofeout, impofeout, feorgadsout, fecolout
      real remifeout, thetamaxout, deffeout, thetachlout
      real chlout, feprimeout, o2, chl_D_out

      real p1, p2, kfevar, kfevar_C

# if defined O_mobi_silicon
      real biodiat,jmax_Diat,u_Diat,avej_Diat,gd_Diat,npp_Diat
      real npp_Diatout,morp_Diat,morpt_Diat, biosil, bioopl, gl_Diat
      real nudt,g_Diat,diatflag,graz_Diat,graz_Diatout,morp_Diatout
      real morpt_Diatout, ing_Diat, felimit_Diat
      real oplpro, opldis, opldisout, oplproout, opl_disk1
      real limSi, limP_Diat, dopupt_Diat_flag, dopupt_Diat
      real dig_Diat, excr_Diat, npp_Diat_dopout, kfevar_Diat
      real silflag, oplflag, impoopl, expoopl, expooplout, wwo
#  if defined O_mobi_iron
      real thetamaxout_Diat
      real thetamax_Diat, deffe_Diat, sf_Diat
      real deffeout_Diat
#  endif
# endif

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "mobi.h"

      real bioin(ntnpzd), bioout(ntnpzd)

      bioout(:) = 0.0
      biopo4 = bioin(imobipo4)
      biophyt = bioin(imobiphyt)
      biozoop = bioin(imobizoop)
      biodetr = bioin(imobidetr)
# if defined O_carbon      
      biodic = bioin(imobidic)
# endif
# if defined O_mobi_nitrogen
      biodop = bioin(imobidop)
      biono3 = bioin(imobino3)
      biodon = bioin(imobidon)
      biodiaz = bioin(imobidiaz)
#  if defined O_mobi_nitrogen_15
      biodin15 = bioin(imobidin15)
      biodon15 = bioin(imobidon15)
      biophytn15 = bioin(imobiphytn15)
      biozoopn15 = bioin(imobizoopn15)
      biodetrn15 = bioin(imobidetrn15)
      biodiazn15 = bioin(imobidiazn15)
#   if defined O_mobi_silicon
      biodiatn15 = bioin(imobidiatn15)
#   endif      
#  endif
# endif
# if defined O_carbon_13
      biodic13 = bioin(imobidic13)
      biophytc13 = bioin(imobiphytc13)
      biozoopc13 = bioin(imobizoopc13)
      biodetrc13 = bioin(imobidetrc13)
#  if defined O_mobi_nitrogen
      biodoc13 = bioin(imobidoc13)
      biodiazc13 = bioin(imobidiazc13)
#  endif
#  if defined O_mobi_silicon
      biodiatc13 = bioin(imobidiatc13)
#  endif
#  if defined O_mobi_caco3
      biocaco3c13 = bioin(imobicaco3c13)
#  endif
# endif
# if defined O_mobi_caco3
      biocaco3 = bioin(imobicaco3)
# endif
# if defined O_kk_ballast
      biod_B = bioin(imobidetr_B)
# endif
# if defined O_mobi_silicon
      biodiat = bioin(imobidiat)
      biosil = bioin(imobisil)
      bioopl = bioin(imobiopl)
# endif
# if defined O_mobi_iron
      biodfe = bioin(imobidfe)
      biodetrfe = bioin(imobidetrfe)
# endif

!       flags prevent negative values by setting outgoing fluxes to
!     zero if tracers are lower than trcmin     
      po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
      phytflag = 0.5 + sign(0.5,biophyt - trcmin)
      zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
      detrflag = 0.5 + sign(0.5,biodetr - trcmin)
!     set defaults
      no3flag = 1.
      dopflag = 1.
      donflag = 1.
      diazflag = 1.
      din15flag = 1.
      don15flag = 1.
      phytn15flag = 1.
      zoopn15flag = 1.
      detrn15flag = 1.
      diazn15flag = 1.
      dic13flag = 1.
      doc13flag = 1.
      phytc13flag = 1.
      zoopc13flag = 1.
      detrc13flag = 1.
      diazc13flag = 1.
      dfeflag = 1.
      detrfeflag = 1.
      diatn15flag = 1.
      diatc13flag = 1.      
      caco3c13flag = 1.      
# if defined O_mobi_nitrogen
      dopflag = 0.5 + sign(0.5,biodop - trcmin)      
      no3flag = 0.5 + sign(0.5,biono3 - trcmin)
      donflag = 0.5 + sign(0.5,biodon - trcmin)
      diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_mobi_nitrogen_15
      din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
      don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
      phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
#   if defined O_mobi_silicon
      diatn15flag = 0.5 + sign(0.5, biodiatn15 - trcmin)
#   endif
      zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
      detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
      diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_carbon_13
      dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
      phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
#  if defined O_mobi_silicon
      diatc13flag = 0.5 + sign(0.5,biodiatc13 - trcmin)
#  endif      
#  if defined O_mobi_caco3
      caco3c13flag = 0.5 + sign(0.5,biocaco3c13 - trcmin)
#  endif      
      zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
      detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
#  if defined O_mobi_nitrogen
      doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
      diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
#  endif
# endif
# if defined O_mobi_iron
      dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
      detrfeflag = 0.5 + sign(0.5,biodetrfe - trcmin)
# endif
# if defined O_kk_ballast
      dflag_B = 0.5 + sign(0.5,biod_B - trcmin)
# endif
# if defined O_mobi_caco3
      caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)
# endif
# if defined O_mobi_silicon
      diatflag = 0.5 + sign(0.5,biodiat - trcmin)
      silflag = 0.5 + sign(0.5,biosil - trcmin)
      oplflag = 0.5 + sign(0.5,bioopl - trcmin)
# endif

!     limit tracers to positive values
      bioin(:) = max(bioin(:), trcmin)
      biopo4 = max(biopo4, trcmin)
      biophyt = max(biophyt, trcmin)
      biozoop = max(biozoop, trcmin)
      biodetr = max(biodetr, trcmin)
# if defined O_carbon
      biodic = max(biodic, trcmin)
# endif
# if defined O_mobi_nitrogen
      biono3 = max(biono3, trcmin)
      biodop = max(biodop, trcmin)      
      biodon = max(biodon, trcmin)
      biodiaz = max(biodiaz, trcmin)
#  if defined O_mobi_nitrogen_15
      biodin15 = max(biodin15, trcmin)
      biodon15 = max(biodon15, trcmin)
      biophytn15 = max(biophytn15, trcmin)
#   if defined O_mobi_silicon
      biodiatn15 = max(biodiatn15, trcmin)
#   endif      
      biozoopn15 = max(biozoopn15, trcmin)
      biodetrn15 = max(biodetrn15, trcmin)
      biodiazn15 = max(biodiazn15, trcmin)
#  endif
# endif
# if defined O_carbon_13
      biodic13 = max(biodic13, trcmin)
      biophytc13 = max(biophytc13, trcmin)
#   if defined O_mobi_silicon
      biodiatc13 = max(biodiatc13, trcmin)
#   endif 
#   if defined O_mobi_caco3
      biocaco3c13 = max(biocaco3c13, trcmin)
#   endif 
      biozoopc13 = max(biozoopc13, trcmin)
      biodetrc13 = max(biodetrc13, trcmin)    
#   if defined O_mobi_nitrogen
      biodoc13 = max(biodoc13, trcmin)
      biodiazc13 = max(biodiazc13, trcmin)
#   endif 
# endif
# if defined O_kk_ballast
      biod_B = max(biod_B, trcmin)
# endif
# if defined O_mobi_caco3
      biocaco3 = max(biocaco3, trcmin)      
# endif
# if defined O_mobi_silicon
      biodiat = max(biodiat, trcmin)
      biosil = max(biosil, trcmin)
      bioopl = max(bioopl, trcmin)
# endif
# if defined O_mobi_iron
      biodfe = max(biodfe, trcmin)
      biodetrfe = max(biodetrfe, trcmin)
# endif

CAAA I've moved this block from above and replaced bioin(itrc) with biotrc
# if defined O_mobi_iron
      p1 = min(biophyt, pmax) !phytoplankton
      p2 = max(0.0, biophyt-pmax)
      kfevar = (kfemin*p1 + kfemax*p2)/(p1+p2) 
      deffe = biodfe/(kfevar + biodfe)
      thetamax = thetamaxlo + (thetamaxhi - thetamaxlo)*deffe
      alpha_O = alphamin + (alphamax - alphamin)*deffe
      gl_O = gl*thetamax*alpha_O
#  if defined O_mobi_silicon
CAAA I added variable iron half sat
      p1 = min(biodiat, pmax_Diat) !diatoms
      p2 = max(0.0, biodiat-pmax_Diat)
      kfevar_Diat = (kfemin_Diat*p1 + kfemax_Diat*p2)
     &            /(p1 + p2)
      deffe_Diat = biodfe/(kfevar_Diat + biodfe)
c      deffe_Diat = biodfe/(kfe_Diat + biodfe)
      thetamax_Diat = thetamaxlo + (thetamaxhi-thetamaxlo)*deffe_Diat
      alpha_Diat = alphamin + (alphamax-alphamin)*deffe_Diat
      gl_Diat = gl*thetamax_Diat*alpha_Diat
#  else
      kfevar=(kfemin*p1+kfemax*p2)/(p1+p2) 
#  endif
#  if defined O_mobi_nitrogen
      deffe_D = biodfe/(kfe_D + biodfe)
      thetamax_D = thetamaxlo + (thetamaxhi-thetamaxlo)*deffe_D
      alpha_D = alphamin + (alphamax-alphamin)*deffe_D
      gl_D = gl*thetamax_D*alpha_D
#  endif
# endif

!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      kirr=-kw - kc*(biophyt
#  if defined O_mobi_nitrogen        
     &     + biodiaz
#  endif
#  if defined O_mobi_silicon
     &     + biodiat
#  endif
     &     )
#  if defined O_mobi_caco3
     &     - kc_c*biocaco3
#  endif
#  if defined O_mobi_silicon
c     &     - kc_o*bioopl
#  endif
    
      f1 = exp(kirr*dzt)

# if defined O_mobi_iron
      jmax = abio_P*bct*deffe
# else
      jmax = abio_P*bct
# endif
      gd = jmax*dayfrac
# if defined O_mobi_iron
      u1 = max(gl_O/gd,1.e-6)
# else
      u1 = max(gl/gd,1.e-6)
# endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2

      avej = gd*(phi1 - phi2)/(-kirr*dzt)
!     Make the max grazing rate a function of temperature
!     bctz sets an upper limit on the effects of temp on grazing
!     in contrast to phytoplankton growth rates bct, which are unlimited
      gmax = gbio*bctz

# if defined O_mobi_nitrogen
#  if defined O_mobi_iron
      jmax_D = max(0.,abio_P*(bct - dbct_D)*deffe_D)*jdiar
#  else
      jmax_D = max(0.,abio_P*(bct - dbct_D))*jdiar
#  endif
!
      gd_D = max(1.e-14,jmax_D*dayfrac)
#  if defined O_mobi_iron
      u1 = max(gl_D/gd_D,1.e-6)
#  else
      u1 = max(gl/gd_D,1.e-6)
#  endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*(phi1 - phi2)/(-kirr*dzt)
# endif

# if defined O_mobi_silicon
#  if defined O_mobi_iron
      jmax_Diat = abiodiat*bct*deffe_Diat
#  else
      jmax_Diat = abiodiat*bct
#  endif
      gd_Diat = jmax_Diat*dayfrac
#  if defined O_mobi_iron
      u1 = max(gl_Diat/gd_Diat,1.e-6)
#  else
      u1 = max(gl*alpha_Diat/gd_Diat,1.e-6)
#  endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_Diat = gd_Diat*(phi1 - phi2)/(-kirr*dzt)
# endif

      nupt = nupt0*bct
# if defined O_mobi_nitrogen
      nupt_D = nupt0_D*bct
      nfixout = 0.0
# endif
      expoout = 0.0
      grazout = 0.0
      morpout = 0.0
      morzout = 0.0
      graz_Det_out = 0.0
      graz_Zout = 0.0
# if defined O_mobi_nitrogen_15
      rn15expoout = 0.0
# endif
# if defined O_carbon_13
      rc13expoout = 0.0
#  if defined O_mobi_caco3
      rcaco3c13expoout = 0.0
#  endif
# endif
      calproout = 0.0
# if defined O_mobi_caco3
      disslout = 0.0
      expocaco3out = 0.0
# endif
# if defined O_kk_ballast
      expoout_B = 0.0
      graz_Det_out_B = 0.0
      remiout_B = 0.0
# endif
# if defined O_mobi_silicon
      nudt = nudt0*bct
      npp_Diatout = 0.0
      graz_Diatout = 0.0
      morp_Diatout = 0.0
      morpt_Diatout = 0.0
      thetamaxout_Diat = 0.0
      expooplout = 0.0
      opldisout = 0.0
      oplproout = 0.0
# endif
# if defined O_mobi_iron
      expofeout = 0.0
      remifeout = 0.0
# endif
# if defined O_save_mobi_fluxes
      nppout = 0.0
      npp_dopout = 0.0
      morptout = 0.0
      remiout = 0.0
      excrout = 0.0
#  if defined O_mobi_nitrogen
      npp_Dout = 0.0
      npp_D_dopout = 0.0
      graz_Dout = 0.0
      morpt_Dout = 0.0
      morp_Dout = 0.0
#    if defined O_mobi_silicon
      npp_Diat_dopout = 0.0
#    endif
#  endif
#  if defined O_mobi_iron
#   if defined O_save_mobi_diagnostics
      feorgadsout = 0.0
      deffeout = 0.0
      feprimeout = 0.0
      fecolout = 0.0
#    if defined O_mobi_silicon
      deffeout_Diat = 0.0
#    endif      
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
      avej_out = 0.0
      avej_D_out = 0.0
      gmax_out = 0.0
      no3P_out = 0.0
      po4P_out = 0.0
      po4_D_out = 0.0
#  endif
# endif

      do n=1,nbio

        p1 = min(biophyt,pmax)
        p2 = max(0.0,biophyt - pmax)
        k1n = (knmin*p1 + knmax*p2)/(p1 + p2)
        k1p_P = k1n*redptn
# if defined O_mobi_iron
        kfevar = (kfemin*p1 + kfemax*p2)/(p1 + p2)
        deffe = biodfe/(kfevar + biodfe)
        jmax = abio_P*bct*deffe
#  if defined O_mobi_silicon
        p1 = min(biodiat, pmax_Diat) !diatoms
        p2 = max(0.0, biodiat-pmax_Diat)
        kfevar_Diat = (kfemin_Diat*p1 + kfemax_Diat*p2)
     &               /(p1 + p2)
        k1n_Diat = (knmin_Diat*p1 + knmax_Diat*p2)/(p1 + p2)
        k1p_Diat = k1n_Diat*redptn
        deffe_Diat = biodfe/(kfevar_Diat + biodfe)
        jmax_Diat = abiodiat*bct*deffe_Diat 
#  endif
#  if defined O_mobi_nitrogen
        deffe_D = biodfe/(kfe_D + biodfe)
        jmax_D = max(0.,abio_P*(bct - dbct_D)*deffe_D)*jdiar
#  endif
# endif 

# if defined O_mobi_nitrogen
!       growth rate of phytoplankton
!       consume DOP when it is more efficient
        limP_dop = hdop*biodop/(k1p_P + biodop)
        limP_po4 = biopo4/(k1p_P + biopo4)
        dopupt_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP = limP_dop*dopupt_flag + limP_po4*(1.-dopupt_flag)
        u_P = min(avej, jmax*limP)
#  if defined O_mobi_silicon
!       variable half saturation for Si (Aumont et al, 2003) in mmol m-3
c        k1si = 0.8 + 7.2*(biosil/(30.+biosil)) 
! the above leads to small values for k1si and thus almost no Si limitation
! anywhere
        k1si = 5.e-3 ! mol m-3
        limSi = biosil/(k1si + biosil)
        limP_dop = hdop*biodop/(k1p_Diat + biodop)
        limP_po4 = biopo4/(k1p_Diat + biopo4)
        dopupt_Diat_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP_Diat = limP_dop*dopupt_Diat_flag 
     &              + limP_po4*(1.-dopupt_Diat_flag)
        u_Diat = min(avej_Diat, jmax_Diat*limSi)
        u_Diat = min(u_Diat, jmax_Diat*limP_Diat)
#  endif
# else
        limP = biopo4/(k1p_P + biopo4)
        u_P = min(avej, jmax*limP)        
#  if defined O_mobi_silicon
        limSi = biosil/(k1si + biosil)
        limP_Diat = biopo4/(k1p_Diat + biopo4)
        u_Diat = min(avej_Diat, jmax_Diat*limSi)
        u_Diat = min(u_Diat, jmax_Diat*limP_Diat)
#  endif
# endif        
        po4P = jmax*biopo4/(k1p_P + biopo4) 

# if defined O_mobi_nitrogen
!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
#  if defined O_mobi_silicon
        u_Diat = min(u_Diat, jmax_Diat*biono3/(k1n_Diat + biono3))
#  endif
        no3P = jmax*biono3/(k1n + biono3)
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        u_D = min(avej_D, jmax_D*limP)
c        u_D = u_P
        dopupt_D_flag = dopupt_flag
        po4_D = jmax_D*biopo4/(k1p_P + biopo4)
!       Set the grazing coefficients for the N case
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop
     &         + zprefDiaz*biodiaz + kzoo
#  if defined O_kk_ballast
     &         + zprefDet*biod_B
#  endif
#  if defined O_mobi_silicon
     &         + zprefDiat*biodiat
#  endif
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
        ing_D = zprefDiaz/thetaZ
#  if defined O_mobi_silicon
        ing_Diat = zprefDiat/thetaZ
#  endif
# else
!       If "else" then set the grazing coefficients for the no N case
!       note kzoo is in terms of N so convert to P
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop 
     &         + kzoo
#  if defined O_kk_ballast
     &         + zprefDet*biod_B
#  endif
#  if defined O_mobi_silicon
     &         + zprefDiat*biodiat
#  endif
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
#  if defined O_mobi_silicon
        ing_Diat = zprefDiat/thetaZ
#  endif        
# endif
        npp = u_P*biophyt
# if defined O_mobi_silicon
        npp_Diat = u_Diat*biodiat
# endif
# if defined O_mobi_nitrogen
        dopupt = npp*dopupt_flag        
#  if defined O_mobi_silicon
        dopupt_Diat = npp_Diat*dopupt_Diat_flag
#  endif
        npp_D = max(0.,u_D*biodiaz)
!       grazing on diazotrophs
        g_D = gmax*ing_D*biodiaz
        graz_D = g_D*biozoop
        morpt_D = nupt_D*biodiaz ! linear mortality
        morp_D = nup_D*biodiaz*biodiaz
c        no3upt_D = biono3/(k1n + biono3)*npp_D ! nitrate uptake
        no3upt_D = (0.5+0.5*tanh(biono3-5.))*npp_D ! nitrate uptake
        dopupt_D = npp_D*dopupt_D_flag
# endif
!       grazing on P
        g_P = gmax*ing_P*biophyt
        graz = g_P*biozoop
!       grazing on Z
        g_Z = gmax*ing_Z*biozoop
        graz_Z = g_Z*biozoop
!       grazing on Detritus
        g_Det = gmax*ing_Det*biodetr
        graz_Det = g_Det*biozoop
!
        morp = nup*biophyt
        morpt = nupt*biophyt
# if defined O_mobi_nitrogen        
        recy_don = nudon*bct*biodon
        recy_dop = nudop*bct*biodop
# endif        
        morz = nuz*biozoop*biozoop
        remi = nud*bct*biodetr
        expo = wwd*biodetr
# if defined O_mobi_caco3
        dissl = biocaco3*dissk1
        expocaco3 = wwc*biocaco3
# endif
# if defined O_kk_ballast
        g_Det_B = gmax*ing_Det*biod_B
        graz_Det_B = g_Det_B*biozoop
        remi_B = max(bapr*dissl/(capr*redctn*1.e3), 0.) !mmol N/m3
        expo_B = wwc*biod_B !ballasted by CaCO3
# endif
# if defined O_mobi_silicon
        g_Diat =  gmax*ing_Diat*biodiat
        graz_Diat = g_Diat*biozoop 
        morp_Diat = nu_diat*biodiat
        morpt_Diat = nudt*biodiat
        opldis = bioopl*opl_disk1
        expoopl = wwo*bioopl
# endif
# if defined O_mobi_iron
!       remineralization of iron from organic matter
        remife=nud*bct*biodetrfe

!       Scavenging of dissolved iron is based on Honeymoon et al. (1988) 
!       and Parekh et al. (2004). 
!       o2flag is zero for o2 < o2min and one otherwise
        o2flag = 0.5 + sign(0.5, o2-o2min)
        fepa = (1.0 + kfeleq * (lig - biodfe))*o2flag
     
        feprime = ((-fepa +(fepa
     &         * fepa + 4.0 * kfeleq 
     &         * biodfe)**(0.5)) /(2.0 * kfeleq))*o2flag

        feorgads = (kfeorg*(((biodetr * detrflag
#  if defined O_kk_ballast
     &         + biod_B
#  endif          
     &         )*mc*redctn)**0.58)*feprime)*o2flag
        fecol = kfecol*feprime*o2flag

        expofe = wwd * biodetrfe
# endif
        graz = graz*phytflag
# if defined O_mobi_nitrogen_15
     &         *phytn15flag
# endif
        graz_Z = graz_Z*zoopflag
# if defined O_mobi_nitrogen_15
     &          *zoopn15flag
# endif
        graz_Det = graz_Det*detrflag
# if defined O_mobi_nitrogen_15
     &            *detrn15flag
# endif
        morp = morp*phytflag
# if defined O_mobi_nitrogen_15
     &        *phytn15flag
# endif
        morpt = morpt*phytflag
# if defined O_mobi_nitrogen_15
     &        *phytn15flag
# endif
        morz = morz*zoopflag
# if defined O_mobi_nitrogen_15
     &        *zoopn15flag
# endif
        remi = remi*detrflag
# if defined O_mobi_nitrogen_15
     &        *detrn15flag
# endif
        expo = expo*detrflag
# if defined O_mobi_nitrogen_15
     &        *detrn15flag
# endif
# if defined O_mobi_nitrogen
        recy_dop = recy_dop*dopflag               
        npp = npp*no3flag*(dopupt_flag*dopflag
     &         + (1.-dopupt_flag)*po4flag)
#  if defined O_mobi_nitrogen_15
     &           *din15flag
#  endif
#  if defined O_mobi_silicon
        npp_Diat = npp_Diat*no3flag*(dopupt_Diat_flag*dopflag
     &           + (1.-dopupt_Diat_flag)*po4flag)
#   if defined O_mobi_nitrogen_15
     &             *din15flag
#   endif
#  endif
        npp_D = npp_D*(dopupt_D_flag*dopflag
     &           + (1.-dopupt_D_flag)*po4flag)
#   if defined O_mobi_nitrogen_15
     &             *din15flag
#   endif
        graz_D = graz_D*diazflag
#  if defined O_mobi_nitrogen_15
     &          *diazn15flag
#  endif
        morpt_D = morpt_D*diazflag
#  if defined O_mobi_nitrogen_15
     &           *diazn15flag
#  endif
        morp_D = morp_D*diazflag
#  if defined O_mobi_nitrogen_15
     &           *diazn15flag
#  endif
        no3upt_D = no3upt_D*no3flag
#  if defined O_mobi_nitrogen_15
     &            *din15flag
#  endif
        recy_don = recy_don*donflag
#  if defined O_mobi_nitrogen_15
     &           *don15flag
#  endif
# else
        npp = npp*po4flag
#  if defined O_mobi_silicon
        npp_Diat = npp_Diat*po4flag
#  endif
# endif
# if defined O_kk_ballast
        graz_Det_B = graz_Det_B*dflag_B
#  if defined O_mobi_nitrogen_15
     &            *zoopn15flag
#  endif
        remi_B = remi_B*dflag_B
        expo_B = expo_B*dflag_B
# endif
# if defined O_mobi_caco3
        dissl = dissl*caco3flag
        expocaco3 = expocaco3*caco3flag
# endif
# if defined O_mobi_silicon
        graz_Diat = graz_Diat*diatflag
        morp_Diat = morp_Diat*diatflag
        morpt_Diat = morpt_Diat*diatflag
# endif
# if defined O_mobi_iron
        remife=remife*detrfeflag
        feorgads=feorgads*dfeflag
        expofe=expofe*detrfeflag
        fecol=fecol*dfeflag
# endif
!     digestion of grazed material
        dig_P = gamma1*graz
        dig_Z = gamma1*graz_Z
        dig_Det = gamma1*graz_Det
# if defined O_kk_ballast
        dig_Det_B = gamma1*graz_Det_B
# endif        
# if defined O_mobi_silicon
        dig_Diat = gamma1*graz_Diat
# endif        
        dig = dig_P + dig_Z + dig_Det
# if defined O_kk_ballast
     &      + dig_Det_B
# endif        
# if defined O_mobi_silicon
     &      + dig_Diat
# endif
!     excretion based on growth efficiency
        excr_P = gamma1*(1-geZ)*graz
        excr_Z = gamma1*(1-geZ)*graz_Z
        excr_Det = gamma1*(1-geZ)*graz_Det
# if defined O_kk_ballast
       excr_Det_B = gamma1*(1-geZ)*graz_Det_B        
# endif
CSPKKK the following block is based on the above block
# if defined O_mobi_silicon
       excr_Diat = gamma1*(1-geZ)*graz_Diat
# endif
        excr = excr_P + excr_Z + excr_Det
# if defined O_kk_ballast
     &       + excr_Det_B
# endif  
# if defined O_mobi_silicon
     &       + excr_Diat
# endif
!     sloppy feeding
        sf_P = (1.-gamma1)*graz
        sf_Z = (1.-gamma1)*graz_Z
        sf_Det = (1.-gamma1)*graz_Det
# if defined O_kk_ballast
        sf_Det_B = (1.-gamma1)*graz_Det_B
# endif
CSPKKK the following block is based on the above block
# if defined O_mobi_silicon
        sf_Diat = (1.-gamma1)*graz_Diat
# endif
        sf = sf_P + sf_Z + sf_Det   
# if defined O_kk_ballast
     &       + sf_Det_B
# endif  
CSPKKK the following block is based on the above block
# if defined O_mobi_silicon
     &       + sf_Diat
# endif
# if defined O_mobi_nitrogen
!     digestion of grazed material
        dig_D = gamma1*graz_D*(redntp/diazntp)
        dig = dig + dig_D

!     excretion based on growth efficiency
        excr_D = gamma1*(1-geZ)*graz_D*(redntp/diazntp)
        excr = excr + excr_D
!       excrete "extra" non-Redfield diazotroph N
        nr_excr_D = gamma1*graz_D*(1-(redntp/diazntp)) 
     &       + (1-gamma1)*graz_D*(1-(redntp/diazntp))

!     sloppy feeding
        sf_D = (1-gamma1)*graz_D*(redntp/diazntp)
        sf = sf + sf_D
# endif
# if defined O_mobi_nitrogen_15
!       calculate isotope parameters
!       See Somes et al., 2010, GBC for details/results
        uno3 = npp*dtbio/biono3
        uno3 = min(uno3, 0.999)
        uno3 = max(uno3, trcmin)
        rno3 = biodin15/(biono3-biodin15)
        rno3 = min(rno3, 2*rn15std)
        rno3 = max(rno3, rn15std/2.)
        bassim = rno3 + eps_assim*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
        fcassim = bassim/(1+bassim)

        udon = recy_don*dtbio/biodon
        udon = min(udon, 0.999)
        udon = max(udon, trcmin)
        rdon = biodon15/(biodon-biodon15)
        rdon = min(rdon, 2*rn15std)
        rdon = max(rdon, rn15std/2.)
        brecy = rdon + eps_recy*(1-udon)/udon*log(1-udon)*rdon/1000.
        fcrecy = brecy/(1+brecy)

        rzoop = biozoopn15/(biozoop-biozoopn15)
        rzoop = min(rzoop, 2.*rn15std)
        rzoop = max(rzoop, rn15std/2.)
        bexcr = rzoop - eps_excr*rzoop/1000.
        fcexcr = bexcr/(1+bexcr)

        bnfix = rn15std - eps_nfix*rn15std/1000.
        fcnfix = bnfix/(1+bnfix)

        rtdin15 = biodin15/biono3
        rtdin15 = min(rtdin15, 2*rn15std/(1+rn15std))
        rtdin15 = max(rtdin15, rn15std/(1+rn15std)/2.)

        rtdon15 = biodon15/biodon
        rtdon15 = min(rtdon15, 2*rn15std/(1+rn15std))
        rtdon15 = max(rtdon15, rn15std/(1+rn15std)/2.)

        rtphytn15 = biophytn15/biophyt
        rtphytn15 = min(rtphytn15, 2.*rn15std/(1+rn15std))
        rtphytn15 = max(rtphytn15, rn15std/(1+rn15std)/2.)

#  if defined O_mobi_silicon
        rtdiatn15 = biodiatn15/biodiat
        rtdiatn15 = min(rtdiatn15, 2.*rn15std/(1+rn15std))
        rtdiatn15 = max(rtdiatn15, rn15std/(1+rn15std)/2.)
#  endif
        
        rtzoopn15 = biozoopn15/biozoop
        rtzoopn15 = min(rtzoopn15, 2.*rn15std/(1+rn15std))
        rtzoopn15 = max(rtzoopn15, rn15std/(1+rn15std)/2.)

        rtdetrn15 = biodetrn15/biodetr
        rtdetrn15 = min(rtdetrn15, 2.*rn15std/(1+rn15std))
        rtdetrn15 = max(rtdetrn15, rn15std/(1+rn15std)/2.)

        rtdiazn15 = biodiazn15/biodiaz
        rtdiazn15 = min(rtdiazn15, 2.*rn15std/(1+rn15std))
        rtdiazn15 = max(rtdiazn15, rn15std/(1+rn15std)/2.)
# endif
# if defined O_carbon_13

        rdic13 = biodic13/(biodic-biodic13)
        rdic13 = min(rdic13, 2.*rc13std)
        rdic13 = max(rdic13, 0.5*rc13std)
        bc13npp = ac13b*rdic13
        fcnpp = bc13npp/(1+bc13npp)
        
        rtdic13 = biodic13/biodic
        rtdic13 = min(rtdic13, 2*rc13std/(1+rc13std))
        rtdic13 = max(rtdic13, 0.5*rc13std/(1+rc13std))

        rtphytc13 = biophytc13/(biophyt*redctn)
        rtphytc13 = min(rtphytc13, 2.*rc13std/(1+rc13std))
        rtphytc13 = max(rtphytc13, 0.5*rc13std/(1+rc13std))

#  if defined O_mobi_silicon
        rtdiatc13 = biodiatc13/(biodiat*redctn)
        rtdiatc13 = min(rtdiatc13, 2.*rc13std/(1+rc13std))
        rtdiatc13 = max(rtdiatc13, 0.5*rc13std/(1+rc13std))
#  endif        

#  if defined O_mobi_caco3
        rtcaco3c13 = biocaco3c13/biocaco3
        rtcaco3c13 = min(rtcaco3c13, 2.*rc13std/(1+rc13std))
        rtcaco3c13 = max(rtcaco3c13, 0.5*rc13std/(1+rc13std))
#  endif        

        rtzoopc13 = biozoopc13/(biozoop*redctn)
        rtzoopc13 = min(rtzoopc13, 2.*rc13std/(1+rc13std))
        rtzoopc13 = max(rtzoopc13, 0.5*rc13std/(1+rc13std))

        rtdetrc13 = biodetrc13/(biodetr*redctn)
        rtdetrc13 = min(rtdetrc13, 2.*rc13std/(1+rc13std))
        rtdetrc13 = max(rtdetrc13, 0.5*rc13std/(1+rc13std))

#  if defined O_mobi_nitrogen
        rtdoc13 = biodoc13/(biodon*redctn)
        rtdoc13 = min(rtdoc13, 2*rc13std/(1+rc13std))
        rtdoc13 = max(rtdoc13, 0.5*rc13std/(1+rc13std))

        rtdiazc13 = biodiazc13/(biodiaz*redctn)
        rtdiazc13 = min(rtdiazc13, 2.*rc13std/(1+rc13std))
        rtdiazc13 = max(rtdiazc13, 0.5*rc13std/(1+rc13std))

#  endif        
# endif
# if defined O_mobi_caco3
!       CaCO3 production
        calpro = ((sf_Z + morz)*capr+(sf_P+morp)*capr)*redctn*1.e3 !stay in mmol
c        calpro = (sf_C + sf_Z + morp_C + morz)*capr*redctn*1.e3 !stay in mmol
# else
        calpro = (morp+morz+(graz+graz_Z)*(1.-gamma1))*capr*redctn*1.e3
# endif
# if defined O_mobi_silicon
!       opal production
        oplpro = (morp_Diat + sf_Diat)*sipr0*redctn*silflag ! mol Si m-3 s-1
        opldis = opldis*oplflag
        expoopl = expoopl*oplflag
# endif
!       Update prognostic equations
# if defined O_mobi_nitrogen
!       nutrients equation
        biopo4 = biopo4 + dtbio*(redptn*(excr + remi 
     &           + (1.-dfrt)*morpt - (npp-dopupt)
#  if defined O_mobi_silicon
     &           + (1.-dfrt)*morpt_Diat - (npp_Diat-dopupt_Diat)
#  endif
     &            )
     &         + diazptn*(morpt_D - (npp_D-dopupt_D)) + recy_dop
     &          )

!       DOP equation
        biodop = biodop + dtbio*(redptn*(
     &             dfr*morp + dfrt*morpt - dopupt
#  if defined O_mobi_silicon
     &           + dfr*morp_Diat + dfrt*morpt_Diat - dopupt_Diat
#  endif
     &            )
     &         - diazptn*dopupt_D - recy_dop)

!       phytoplankton equation
        biophyt = biophyt + dtbio*(npp - morp - graz - morpt)
!       zooplankton equation
        biozoop = biozoop + dtbio*(dig - morz - graz_Z - excr)

!       detritus equation
#  if defined O_kk_ballast
        biod_B = biod_B + dtbio*(bapr*((1.-dfr)*morp + sf + morz)
     &         - remi_B - graz_Det_B - expo_B + impo_B
     &         + bapr*morp_D*(redntp/diazntp)
#   if defined O_mobi_silicon
     &          + bapr*(1.-dfr)*morp_Diat
#   endif
     &           )
        biodetr = biodetr + dtbio*((1.-bapr)*((1.-dfr)*morp + sf + morz)
     &          - remi - graz_Det - expo + impo
     &          + (1.-bapr)*morp_D*(redntp/diazntp)
#   if defined O_mobi_silicon
     &          + (1.-bapr)*(1.-dfr)*morp_Diat
#   endif
     &          + remi_B)
#  else        
        biodetr = biodetr + dtbio*((1.-dfr)*morp + sf + morz - remi 
     &          - graz_Det - expo + impo + morp_D*(redntp/diazntp)
#   if defined O_mobi_silicon
     &          + (1.-dfr)*morp_Diat
#   endif
     &           )
#  endif

#  if defined O_carbon
!       DIC equation
        biodic = biodic + dtbio*redctn*(excr + remi 
     &         + (1.-dfrt)*morpt - npp 
#   if defined O_mobi_silicon
     &         + (1.-dfrt)*morpt_Diat - npp_Diat
#   endif
     &         + morpt_D - npp_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )        
#  endif
!       nitrate (NO3) equation
        biono3 = biono3 + dtbio*(excr + remi   
     &         + (1.-dfrt)*morpt - npp 
#  if defined O_mobi_silicon
     &         + (1.-dfrt)*morpt_Diat - npp_Diat
#  endif
     &         + morpt_D - no3upt_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )
!       DON equation
        biodon = biodon + dtbio*(dfr*morp + dfrt*morpt - recy_don
#  if defined O_mobi_silicon
     &         + dfr*morp_Diat + dfrt*morpt_Diat
#  endif
     &       )
!       diazotroph equation
        biodiaz = biodiaz + dtbio*(npp_D - morp_D - morpt_D - graz_D)
# else
!       nutrients equation
        biopo4 = biopo4 + dtbio*redptn*(remi + excr - npp + morpt
#  if defined O_mobi_silicon
     &         - npp_Diat + morpt_Diat
#  endif
     &       )
!       phytoplankton equation
        biophyt = biophyt + dtbio*(npp - morp - graz - morpt)
!     zooplankton equation
        biozoop = biozoop + dtbio*(dig - morz - graz_Z - excr)
!     detritus equation
#  if defined O_kk_ballast
        biod_B = biod_B + dtbio*(bapr*(morp + sf + morz)
     &         - remi_B - graz_Det_B - expo_B + impo_B
     &          + bapr*morp_C
#   if defined O_mobi_silicon
     &          + bapr*morp_Diat
#   endif
     &           )
        biodetr = biodetr + dtbio*((1.-bapr)*(morp + sf + morz)
     &          - remi - graz_Det - expo + impo
     &          + (1.-bapr)*morp_C
#   if defined O_mobi_silicon
     &          + (1.-bapr)*morp_Diat
#   endif
     &          + remi_B) 
#  else        
cAAA this block can be removed if bapr=0 in previous block
        biodetr = biodetr + dtbio*(morp + sf + morz - remi 
     &          - graz_Det - expo + impo
#   if defined O_mobi_silicon
     &          + morp_Diat
#   endif
     &           )
#  endif        
#  if defined O_carbon
!     DIC equation
        biodic = biodic + dtbio*redctn*(morpt + excr + remi - npp        
#   if defined O_mobi_silicon
     &          + morpt_Diat - npp_Diat
#   endif
     &           )  
#  endif
# endif
# if defined O_mobi_caco3
        biocaco3 = biocaco3 + dtbio*(calpro-dissl-expocaco3+impocaco3)
# endif        
# if defined O_mobi_silicon
!       diatoms equation
        biodiat = biodiat + dtbio*(npp_Diat - morp_Diat - graz_Diat
     &       - morpt_Diat)
!       silicate equation in mol
        biosil = biosil + dtbio*(opldis - oplpro)
!       opal equation in mol
        bioopl = bioopl + dtbio*(oplpro - opldis - expoopl + impoopl)
# endif
# if defined O_mobi_iron
!       dissolved iron equation
#  if defined O_mobi_nitrogen
        biodfe = biodfe + dtbio*(rfeton*(excr + (1.-dfrt)*morpt 
     &            - npp + morpt_D - npp_D + recy_don + nr_excr_D 
     &            + morp_D*(1-(redntp/diazntp))) - feorgads 
     &            + remife - fecol
#   if defined O_mobi_silicon
     &            + rfeton*((1.-dfrt)*morpt_Diat - npp_Diat)
#   endif
     &            )   
!     particulate iron equation
        biodetrfe = biodetrfe + dtbio*(rfeton*(sf + (1.-dfr)*morp
     &            + morp_D*(redntp/diazntp) + morz - graz_Det) 
     &            + feorgads + fecol - remife - expofe + impofe        
#   if defined O_mobi_silicon
     &            + rfeton*(1.-dfr)*morp_Diat
#   endif        
     &       )
#  else
        biodfe = biodfe + dtbio*(rfeton*(excr + morpt 
     &            - npp) - feorgads + remife - fecol        
#   if defined O_mobi_silicon
     &            + rfeton*(morpt_Diat - npp_Diat)
#   endif
     &       )
!     particulate iron equation
        biodetrfe = biodetrfe + dtbio*(rfeton*(sf + morp
     &             + morz - graz_Det) 
     &            + feorgads + fecol - remife - expofe + impofe
#   if defined O_mobi_silicon
     &            + rfeton*morp_Diat
#   endif
#   if defined O_kk_ballast
     &            - rfeton*graz_Det_B
#   endif        
     &       )
#  endif
# endif
# if defined O_mobi_nitrogen_15
!       isotope equations
        biodin15 = biodin15 + dtbio*(rtphytn15*(1.-dfrt)*morpt 
#  if defined O_mobi_silicon
     &       + rtdiatn15*(1.-dfrt)*morpt_Diat - fcassim*npp_Diat
#  endif
     &       + fcexcr*excr + rtdiazn15*morpt_D + rtdiazn15*nr_excr_D 
     &       + rtdiazn15*morp_D*(1.-(redntp/diazntp)) + rtdetrn15*remi 
     &       + fcrecy*recy_don - fcassim*npp - fcassim*no3upt_D)

        biodon15 = biodon15 + dtbio*(dfr*rtphytn15*morp 
#  if defined O_mobi_silicon
     &       + dfr*rtdiatn15*morp_Diat + dfrt*rtdiatn15*morpt_Diat
#  endif
     &       + dfrt*rtphytn15*morpt - fcrecy*recy_don)

        biophytn15 = biophytn15 + dtbio*(fcassim*npp - rtphytn15*morp 
     &       - rtphytn15*graz - rtphytn15*morpt)

#  if defined O_mobi_silicon
        biodiatn15 = biodiatn15 + dtbio*(fcassim*npp_Diat
     &       - rtdiatn15*morp_Diat - rtdiatn15*graz_Diat 
     &       - rtdiatn15*morpt_Diat)
#  endif
        biozoopn15 = biozoopn15 + dtbio*(rtphytn15*dig_P 
#  if defined O_mobi_silicon
     &       + rtdiatn15*dig_Diat
#  endif        
     &       + rtzoopn15*dig_Z + rtdetrn15*dig_Det + rtdiazn15*dig_D 
     &       - rtzoopn15*morz - rtzoopn15*graz_Z - fcexcr*excr)

        biodetrn15 = biodetrn15 + dtbio*(rtphytn15*(1.-dfr)*morp 
#  if defined O_mobi_silicon
     &       + rtdiatn15*(1.-dfr)*morp_Diat + rtdiatn15*sf_Diat
#  endif
     &       + rtphytn15*sf_P + rtzoopn15*sf_Z + rtdetrn15*sf_Det 
     &       + rtdiazn15*sf_D + rtzoopn15*morz - rtdetrn15*remi 
     &       - rtdetrn15*graz_Det - rtdetrn15*expo + rn15impo*impo 
     &       + rtdiazn15*morp_D*(redntp/diazntp))

        biodiazn15 = biodiazn15 + dtbio*(fcnfix*(npp_D-no3upt_D)
     &       + fcassim*no3upt_D - rtdiazn15*morp_D - rtdiazn15*graz_D 
     &       - rtdiazn15*morpt_D)
# endif
# if defined O_carbon_13
!     !!!!!!!!!!!!!!!!!! isotope equations       
#  if defined O_mobi_nitrogen
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*(1.-dfrt)*morpt 
     &       + rtzoopc13*excr + rtdiazc13*morpt_D + rtdiazc13*nr_excr_D 
     &       + rtdiazc13*morp_D*(1-(redntp/diazntp)) + rtdetrc13*remi 
#   if defined O_mobi_silicon
     &       + rtdiatc13*(1.-dfrt)*morpt_Diat - fcnpp*npp_Diat
#   endif        
     &       + rtdoc13*recy_don - fcnpp*npp - fcnpp*npp_D)

        biodoc13 = biodoc13 + dtbio*redctn*(dfr*rtphytc13*morp 
#   if defined O_mobi_silicon
     &       + rtdiatc13*(dfr*morp_Diat + dfrt*morpt_Diat)
#   endif       
     &       + rtphytc13*dfrt*morpt - rtdoc13*recy_don)        

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp 
     &       - rtphytc13*morp - rtphytc13*graz - rtphytc13*morpt)
        
        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P 
#   if defined O_mobi_silicon
     &      + rtdiatc13*dig_Diat
#   endif
     &      + rtzoopc13*dig_Z + rtdetrc13*dig_Det + rtdiazc13*dig_D 
     &       - rtzoopc13*morz - rtzoopc13*graz_Z - rtzoopc13*excr)

        biodetrc13 = biodetrc13 + dtbio*redctn*(rtphytc13*(1.-dfr)*morp 
#   if defined O_mobi_silicon
     &       + rtdiatc13*(1.-dfr)*morp_Diat + rtdiatc13*sf_Diat
#   endif
     &       + rtphytc13*sf_P + rtzoopc13*sf_Z + rtdetrc13*sf_Det 
     &       + rtdiazc13*sf_D
     &       + rtzoopc13*morz - rtdetrc13*remi 
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo 
     &       + rtdiazc13*morp_D*(redntp/diazntp))
       
        biodiazc13 = biodiazc13 + dtbio*redctn*(fcnpp*npp_D 
     &       - rtdiazc13*(morp_D + graz_D + morpt_D))
#  else ! no nitrogen
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*morpt 
     &       + rtzoopc13*excr + rtdetrc13*remi - fcnpp*npp
#   if defined O_mobi_silicon
     &       + rtdiatc13*morpt_Diat - fcnpp*npp_Diat
#   endif
     &        )

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp 
     &       - rtphytc13*(morp + graz + morpt))
        
        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P 
#   if defined O_mobi_silicon
     &      + rtdiatc13*dig_Diat
#   endif        
     &      + rtdetrc13*dig_Det 
     &      + rtzoopc13*(dig_Z - morz - graz_Z - excr))
      
        biodetrc13 = biodetrc13 + dtbio*redctn*(
#   if defined O_kk_ballast
     &         rtphytc13*(1.-bapr)*morp 
#   else
     &         rtphytc13*morp 
#   endif        
     &       + rtphytc13*sf_P + rtzoopc13*sf_Z + rtdetrc13*sf_Det 
     &       + rtdiazc13*sf_D
     &       + rtzoopc13*morz - rtdetrc13*remi 
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo 
#   if defined O_mobi_silicon
     &       + rtdiatc13*morp_Diat + rtdiatc13*sf_Diat
#   endif
     &       )

#  endif
#  if defined O_mobi_caco3
        biocaco3c13 = biocaco3c13 + dtbio*(rtdic13*calpro
     &       - rtcaco3c13*dissl - rtcaco3c13*expocaco3 
     &                + rcaco3c13impo)
#  endif
#  if defined O_mobi_silicon
        biodiatc13 = biodiatc13 + dtbio*redctn*(fcnpp*npp_Diat
     &       - rtdiatc13*(morp_Diat + graz_Diat + morpt_Diat))
#  endif
# endif
!       Accumulate output fields
        expoout = expoout + expo
# if defined O_mobi_nitrogen_15
        rn15expoout = rn15expoout + rtdetrn15
# endif
# if defined O_carbon_13
        rc13expoout = rc13expoout + rtdetrc13*expo
#  if defined O_mobi_caco3
        rcaco3c13expoout = rcaco3c13expoout + rtcaco3c13*expocaco3
#  endif
# endif
        calproout = calproout + calpro
# if defined O_mobi_caco3
        calattout = calattout + calatt
        disslout = disslout + dissl
        expocaco3out = expocaco3out + expocaco3
# endif
# if defined O_mobi_silicon
        opldisout = opldisout + opldis
        oplproout = oplproout + oplpro
        expooplout = expooplout + expoopl
# endif
# if defined O_mobi_nitrogen
        nfixout = nfixout + npp_D - no3upt_D
# endif
# if defined O_mobi_iron
        expofeout = expofeout + expofe
        remifeout = remifeout + remife
# endif
# if defined O_save_mobi_fluxes
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
        graz_Det_out = graz_Det_out + graz_Det
        graz_Zout = graz_Zout + graz_Z
        nppout = nppout + npp
        morptout = morptout + morpt
        remiout = remiout + remi
        excrout = excrout + excr
#  if defined O_mobi_nitrogen
        npp_dopout = npp_dopout + npp*dopupt_flag        
#   if defined O_mobi_silicon
        npp_Diat_dopout = npp_Diat_dopout + npp_Diat*dopupt_Diat_flag  
#   endif        
        npp_Dout = npp_Dout + npp_D
        npp_D_dopout = npp_D_dopout + npp_D*dopupt_D_flag
        graz_Dout = graz_Dout + graz_D
        morpt_Dout = morpt_Dout + morpt_D
        morp_Dout = morp_Dout + morp_D
#  endif
#  if defined O_kk_ballast
        expoout_B = expoout_B + expo_B
        remiout_B = remiout_B + remi_B
        graz_Det_out_B = graz_Det_out_B + graz_Det_B
#  endif
#  if defined O_mobi_silicon
        npp_Diatout = npp_Diatout + npp_Diat
        graz_Diatout = graz_Diatout + graz_Diat
        morp_Diatout = morp_Diatout + morp_Diat
        morpt_Diatout = morpt_Diatout + morpt_Diat
#  endif
#  if defined O_mobi_iron
#   if defined O_save_mobi_diagnostics
        feorgadsout = feorgadsout + feorgads
        deffeout = deffeout + deffe
        feprimeout = feprimeout + feprime
        fecolout = fecolout + fecol
#    if defined O_mobi_silicon
        deffeout_Diat = deffeout_Diat + deffe_Diat
        thetamaxout_Diat = thetamax_Diat
#    endif        
#   endif
#  endif
#  if defined O_save_mobi_diagnostics
        avej_out = avej_out + avej
        gmax_out = gmax_out + gmax
        po4P_out = po4P_out + po4P
#   if defined O_mobi_nitrogen
        avej_D_out = avej_D_out + avej_D
        no3P_out = no3P_out + no3P
        po4_D_out = po4_D_out + po4_D
#   endif
#  endif
# endif
!     flags prevent negative values by setting outgoing fluxes to
!     zero if tracers are lower than trcmin
        if (po4flag .eq. 1) po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
        if (phytflag .eq. 1) phytflag = 0.5 + sign(0.5,biophyt - trcmin)
        if (zoopflag .eq. 1) zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
        if (detrflag .eq. 1) detrflag = 0.5 + sign(0.5,biodetr - trcmin)
# if defined O_mobi_nitrogen
        if (no3flag .eq. 1) no3flag = 0.5 + sign(0.5,biono3 - trcmin)
        if (dopflag .eq. 1) dopflag = 0.5 + sign(0.5,biodop - trcmin)        
        if (donflag .eq. 1) donflag = 0.5 + sign(0.5,biodon - trcmin)
        if (diazflag .eq. 1) diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_mobi_nitrogen_15
        if (din15flag .eq. 1) 
     &     din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
        if (don15flag .eq. 1)
     &     don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
        if (phytn15flag .eq. 1) 
     &     phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
#   if defined O_mobi_silicon
        if (diatn15flag .eq. 1) 
     &     diatn15flag = 0.5 + sign(0.5, biodiatn15 - trcmin)
#   endif
        if (zoopn15flag .eq. 1)
     &     zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
        if (detrn15flag .eq. 1)
     &     detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
        if (diazn15flag .eq. 1)
     &     diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_mobi_caco3
        if (caco3flag .eq. 1) 
     &     caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)
#  if defined O_kk_ballast
        if (dflag_B .eq. 1) 
     &     dflag_B = 0.5 + sign(0.5,biod_B - trcmin)
#  endif
#  if defined O_mobi_silicon
        if (diatflag .eq. 1) 
     &     diatflag = 0.5 + sign(0.5,biodiat - trcmin)
        if (silflag .eq. 1) 
     &     silflag = 0.5 + sign(0.5,biosil - trcmin)
        if (oplflag .eq. 1) 
     &     oplflag = 0.5 + sign(0.5,bioopl - trcmin)
#  endif
# endif        
# if defined O_mobi_iron
        if (dfeflag .eq. 1) dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
        if (detrfeflag .eq. 1) detrfeflag = 0.5 
     &                                 + sign(0.5,biodetrfe - trcmin)
# endif
# if defined O_carbon_13
        if (dic13flag .eq. 1)
     &       dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
        if (phytc13flag .eq. 1)
     &       phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
#  if defined O_mobi_silicon
        if (diatc13flag .eq. 1)
     &       diatc13flag = 0.5 + sign(0.5,biodiatc13 - trcmin)
#  endif
#  if defined O_mobi_caco3
        if (caco3c13flag .eq. 1)
     &       caco3c13flag = 0.5 + sign(0.5,biocaco3c13 - trcmin)
#  endif        
        if (zoopc13flag .eq. 1)
     &       zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
        if (detrc13flag .eq. 1)
     &       detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
#  if defined O_mobi_nitrogen
        if (doc13flag .eq. 1)
     &       doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
        if (diazc13flag .eq. 1)
     &       diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
#  endif        
# endif
      enddo
      
      bioout(imobipo4) = biopo4 - bioin(imobipo4)
      bioout(imobiphyt) = biophyt - bioin(imobiphyt)
      bioout(imobizoop) = biozoop - bioin(imobizoop)
      bioout(imobidetr) = biodetr - bioin(imobidetr)
# if defined O_carbon
      bioout(imobidic) = biodic - bioin(imobidic)
# endif
# if defined O_mobi_nitrogen
      bioout(imobidop) = biodop - bioin(imobidop)
      bioout(imobino3) = biono3 - bioin(imobino3)
      bioout(imobidon) = biodon - bioin(imobidon)
      bioout(imobidiaz) = biodiaz - bioin(imobidiaz)
#  if defined O_mobi_nitrogen_15
      bioout(imobidin15) = biodin15 - bioin(imobidin15)
      bioout(imobidon15) = biodon15 - bioin(imobidon15)
      bioout(imobiphytn15) = biophytn15 - bioin(imobiphytn15)
      bioout(imobizoopn15) = biozoopn15 - bioin(imobizoopn15)
      bioout(imobidetrn15) = biodetrn15 - bioin(imobidetrn15)
      bioout(imobidiazn15) = biodiazn15 - bioin(imobidiazn15)
#   if defined O_mobi_silicon
      bioout(imobidiatn15) = biodiatn15 - bioin(imobidiatn15)
#   endif
#  endif
# endif
# if defined O_mobi_caco3
      bioout(imobicaco3) = biocaco3 - bioin(imobicaco3)
#  if defined O_kk_ballast
      bioout(imobidetr_B) = biod_B - bioin(imobidetr_B)
#  endif
# endif
#  if defined O_mobi_silicon
      bioout(imobidiat) = biodiat - bioin(imobidiat)
      bioout(imobisil) = biosil - bioin(imobisil)
      bioout(imobiopl) = bioopl - bioin(imobiopl)
#  endif
# if defined O_mobi_iron
      bioout(imobidfe) = biodfe - bioin(imobidfe)
      bioout(imobidetrfe) = biodetrfe - bioin(imobidetrfe)
# endif
# if defined O_carbon_13
      bioout(imobidic13) = biodic13 - bioin(imobidic13)
      bioout(imobiphytc13) = biophytc13 - bioin(imobiphytc13)
      bioout(imobizoopc13) = biozoopc13 - bioin(imobizoopc13)
      bioout(imobidetrc13) = biodetrc13 - bioin(imobidetrc13)
#  if defined O_mobi_nitrogen
      bioout(imobidoc13) = biodoc13 - bioin(imobidoc13)
      bioout(imobidiazc13) = biodiazc13 - bioin(imobidiazc13)
#  endif      
#  if defined O_mobi_silicon
      bioout(imobidiatc13) = biodiatc13 - bioin(imobidiatc13)
#  endif
#  if defined O_mobi_caco3
      bioout(imobicaco3c13) = biocaco3c13 - bioin(imobicaco3c13)
#  endif
# endif
#endif
      return
      end
