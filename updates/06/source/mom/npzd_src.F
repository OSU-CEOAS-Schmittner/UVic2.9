
      subroutine mobi_init

#if defined O_mom && defined O_npzd 
!     initialize MOBI parameters
      implicit none
#include "size.h"
#include "npzd.h"
#include "calendar.h"
#include "coord.h"
#include "stdunits.h"
#include "scalar.h"

      integer k, ioun, imobi
      real alpha, par
      real alpha_C, sumzpref

# if defined O_npzd_iron
!     IO for fe_hydr
#include "param.h"
#include "pconst.h"
      integer ib(10), ic(10), iou, j, i
      character (120) :: fname, new_file_name
      logical exists
      real c1e3
      real tmpijkm(1:100,1:100,1:19)
! end IO for fe_hydr
# endif
      namelist /npzd/   alpha, kw, kc, abio_P, bbio, cbio, k1n, nup
     &,                 gamma1, gbio, epsbio, nuz, nud0
     &,                 wd0, par, dtnpzd, redctn, ki, redptn, capr
     &,                 dcaco3, nupt0, redotn, jdiar, kzoo, geZ
     &,                 zprefP, zprefZ, zprefDet, zprefDiaz, kfe, kfe_D
# if defined O_npzd_iron
     &,                 kfemin, kfemax
#  if defined O_npzd_caco3
     &,                 kfemin_C, kfemax_C
#  endif
     &,                 kfeorg
# endif
     &,                 nupt0_D, sgbdfac, diazntp, nup_D, dfr, dfrt
     &,                 nudon0, nudop0, eps_assim, eps_excr, eps_nfix
     &,                 eps_wcdeni, eps_bdeni0, eps_recy, hdop
     &,                 mw, mwz, mw_c
     &,                 abio_C, k1n_C, nuct0,nuc, alpha_C
     &,                 zprefC
# if defined O_kk_ballast
     &,                 bapr
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
      capr = 0.022  ! carbonate to carbon production ratio
      dcaco3 = 650000.0  ! remineralisation depth of calcite [cm]
      redotn = 10.6  ! O2/N Redfield ratio (Anderson & Sarmiento, 1994, GBC)
      jdiar = 0.08  ! factor reducing the growth rate of diazotrophs
      kzoo = 0.15  ! half sat. constant for Z grazing [mmol N m^3]
      geZ = 0.6   ! Zooplankton growth efficiency
      sgbdfac = 1.0 ! sub-grid benthic denitrification rate factor 
      diazntp = 28. !diazotroph N:P ratio
      dfr = 0.08   ! phyt mortality refractory/semi-labile DOM fraction
      dfrt = 0.01  ! phyt fast-recy refracotyr/semi-labile DOM fraction
      hdop = 0.4   ! DOP growth rate handicap
      alpha_C = 0.06  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      abio_C = 0.52  ! a; Maximum growth rate parameter coccolithophore [1/day]
      k1n_C = 0.4  ! Half sat const for N uptake (coccolithophores) [mmol/m^3]
      nuc = 0.03  ! Specific mortality rate (coccolithophores) [1/day]
      nuct0 = 0.015  ! Specific mortality rate (coccolithophores) [1/day]
# if defined O_kk_ballast
      bapr = 0.05 ! detritus to carbonate ratio [mg POC/mg PIC]
# endif
      kc_c = 0.047  ! Light atten. by calcite [1/(m*mmol/m^3)]
      wc0 = 35. ! constant calcite sinking speed [m/day]
!      Gehlen et al 2007 use a surface sink speed of 50 m/day inc w/ depth
      kcal = 100.  !half sat for PIC dissolution from Planktom 10
      dissk0 = 0.013 !initial dissolution rate parameter [1/day] 

# if defined O_npzd_nitrogen
#  if !defined O_npzd_o2
      print*,'error in mobi_init:'
      print*,'O_npzd_o2 must be defined if O_npzd_nitrogen is used'
      stop 
#  endif      
#  if defined O_npzd_caco3
      zprefP = 0.225 ! Zooplankton preference for P
      zprefC = 0.225 ! Zooplankton preference for C
      zprefZ = 0.225 ! Zooplankton preference for other Z
      zprefDet = 0.225 ! Zooplankton preference for Detritus
      zprefDiaz = 0.1 ! Zooplankton preference for diazotrophs
#  else      
      zprefP = 0.3       ! Zooplankton preference for P
      zprefZ = 0.3       ! Zooplankton preference for other Z
      zprefDet = 0.3     ! Zooplankton preference for Detritus
      zprefDiaz = 0.1    ! Zooplankton preference for diazotrophs
#  endif      
# else
#  if defined O_npzd_caco3
      zprefP = 0.30 ! Zooplankton preference for P
      zprefC = 0.30 ! Zooplankton preference for C
      zprefZ = 0.30 ! Zooplankton preference for other Z
      zprefDet = 0.10 ! Zooplankton preference for Detritus
#  else      
      zprefP = 0.35       ! Zooplankton preference for P
      zprefZ = 0.35       ! Zooplankton preference for other Z
      zprefDet = 0.30           ! Zooplankton preference for Detritus
#  endif      
# endif
# if defined O_npzd_nitrogen_15
#  if !defined O_npzd_nitrogen
      print*,'error in mobi_init:'
      print*,'O_npzd_nitrogen must be on if O_npzd_nitrogen_15 is used'
      stop
#  endif      
      eps_assim = 6.
      eps_excr = 4.
      eps_nfix = 1.
      eps_wcdeni = 25.
      eps_bdeni0 = 6.
      eps_recy = 1.
# endif
# if defined O_npzd_iron
#  if !defined O_npzd_o2
      print*,'error in mobi_init:'
      print*,'O_npzd_o2 must be defined if O_npzd_iron is used'
      stop 
#  endif      
      kfemin = 0.04e-3    ! Minimum half saturation constant for Fe limitation [mmol Fe / m-3]
      kfemax = 0.2e-3    ! Maximum half saturation constant for Fe limitation [mmol Fe / m-3]
      pmax = 0.15    ! Phytoplankton biomass above which kfe increases [mmol N / m-3]
      kfe_D = 0.1e-3 ! Half saturation constant for Diaz Fe limitation [mmol Fe / m-3]
c     juan: Levin's value was 0.1e-3
#  if defined O_npzd_caco3
      kfemin_C = 0.04e-3    ! Minimum half saturation constant for Fe limitation [mmol Fe / m-3]
      kfemax_C = 0.4e-3    ! Maximum half saturation constant for Fe limitation [mmol Fe / m-3]
      pmax_C = 0.15    ! Phytoplankton biomass above which kfe increases [mmol N / m-3]
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
      read  (ioun, npzd, end=108)
108   continue
      write (stdout, npzd)
      call relunit (ioun)

!     convert units of NPZD parameters to MOM units
      redctn = redctn*1.e-3
      redotn = redotn*1.e-3
      redotp = redotn/redptn
      redctp = redctn/redptn
      redotc = redotn/redctn
      redntp = 1./redptn
      redntc = 1./redctn
      diazptn = 1./diazntp
      k1p_P   = k1n*redptn
# if defined O_npzd_caco3
#  if !defined O_carbon
      print*,'error in mobi_init:'
      print*,'O_carbon must be defined if O_npzd_caco3 is used'
      stop
#  endif      
      k1p_C  = k1n_C*redptn
      alpha_C = alpha_C/daylen
      abio_C = abio_C/daylen
      nuc = nuc/daylen
      nuct0 = nuct0/daylen
      wc0 = wc0*1.e2
      kc_c = kc_c*1.e-2
      dissk0 = dissk0/daylen
# endif
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
# if defined O_npzd_iron
      alphamax = alphamax/daylen
      alphamin = alphamin/daylen
      tap = 2.*par   ! alpha is iron dependent and calculated in npzd_src
# else
      tap = 2.*alpha*par
# endif


!     calculate sinking speed of detritus divided by grid width
      do k=1,km
!     linear increase wd0-200m with depth
         if (zt(k) .lt. mwz) then
            wd(k) = (wd0+mw*zt(k))/daylen/dzt(k) ! [s-1]
# if defined O_npzd_caco3
            wc(k) = (wc0+mw_c*zt(k))/daylen/dzt(k)    ! [s-1]
# endif            
         else
            wd(k) = (wd0+mw*mwz)/daylen/dzt(k) ! [s-1]
# if defined O_npzd_caco3
            wc(k) = (wc0+mw_c*mwz)/daylen/dzt(k)    ! [s-1]
# endif 
         endif
         rkwz(k) = 1./(kw*dzt(k))
      enddo
      ztt(1)=0.0
      do k=1,km-1
         ztt(k+1)=(-1)*zw(k)
      enddo

# if defined O_npzd_nitrogen
!     check grazing preferences = 1 for N case
#  if !defined O_npzd_caco3
      sumzpref = zprefP + zprefDet + zprefZ + zprefDiaz
      IF (sumzpref.ne.1.) THEN
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz",
     &           zprefP, zprefDet, zprefZ, zprefDiaz
         zprefP = zprefP/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         zprefDiaz = zprefDiaz/sumzpref
         print*,"to"
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz",
     &           zprefP, zprefDet, zprefZ, zprefDiaz
      END IF
#  else
      sumzpref = zprefP + zprefDet + zprefZ + zprefDiaz + zprefC
      IF (sumzpref.ne.1.) THEN
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefDiaz, zprefC
         zprefP = zprefP/sumzpref
         zprefC = zprefC/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         zprefDiaz = zprefDiaz/sumzpref
         print*,"to"
         print*,"zprefP, zprefDet, zprefZ, zprefDiaz, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefDiaz, zprefC
      END IF
#  endif
# else
!     check that grazing preferences = 1 for no N case
#  if !defined O_npzd_caco3
      sumzpref = zprefP + zprefDet + zprefZ
      IF (sumzpref.ne.1.) THEN
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ",
     &           zprefP, zprefDet, zprefZ
         zprefP = zprefP/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         print*,"to"
         print*,"zprefP, zprefDet, zprefZ",
     &           zprefP, zprefDet, zprefZ
      END IF
#  else
      sumzpref = zprefP + zprefDet + zprefZ + zprefC
      IF (sumzpref.ne.1.) THEN
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefC
         zprefP = zprefP/sumzpref
         zprefC = zprefC/sumzpref
         zprefZ = zprefZ/sumzpref
         zprefDet = zprefDet/sumzpref
         print*,"to"
         print*,"mobi_init: adjust zprefs from" 
         print*,"zprefP, zprefDet, zprefZ, zprefC",
     &           zprefP, zprefDet, zprefZ, zprefC
      END IF
#  endif
# endif
# if defined O_npzd_iron
#ifndef O_TMM
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
#endif ! not O_TMM
# endif
# if defined O_carbon || defined O_npzd_alk

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
#   if defined O_npzd_nitrogen
      call setimobi('imobidoc13',imobidoc13,imobi)
      call setimobi('imobidiazc13',imobidiazc13,imobi)
#   endif
#   if defined O_npzd_caco3
      call setimobi('imobicoccc13',imobicoccc13,imobi)
      call setimobi('imobicaco3c13',imobicaco3c13,imobi)
#   endif
#  endif
# endif
# if defined O_npzd_nitrogen
      call setimobi('imobidop',imobidop,imobi)
      call setimobi('imobino3',imobino3,imobi)
      call setimobi('imobidon',imobidon,imobi)
      call setimobi('imobidiaz',imobidiaz,imobi)
#  if defined O_npzd_nitrogen_15
      call setimobi('imobidin15',imobidin15,imobi)
      call setimobi('imobidon15',imobidon15,imobi)
      call setimobi('imobiphytn15',imobiphytn15,imobi)
      call setimobi('imobizoopn15',imobizoopn15,imobi)
      call setimobi('imobidetrn15',imobidetrn15,imobi)
      call setimobi('imobidiazn15',imobidiazn15,imobi)
#   if defined O_npzd_caco3
      call setimobi('imobicoccn15',imobicoccn15,imobi)
#   endif
#  endif
# endif
# if defined O_npzd_caco3
      call setimobi('imobicocc',imobicocc,imobi)
      call setimobi('imobicaco3',imobicaco3,imobi)
# endif
# if defined O_kk_ballast
      call setimobi('imobidetr_B',imobidetr_B,imobi)
# endif
# if defined O_npzd_iron
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
     &,                t_in, po4_in
# if defined O_npzd_o2
     &,                o2_in
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
     &,                s_in, dic_in, alk_in, co2_in
# endif
# if defined O_carbon_13      
     &,                dic13_in
# endif
     &,                sgb_in
# if defined O_npzd_nitrogen
     &,                no3_in
#  if defined O_npzd_nitrogen_15
     &,                din15_in
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
!             if O_save_npzd is switch on additional output (rnpp etc.)
!                is in npzd.h 
!     copied from tracer.F (Andreas Schmittner, Oct 2, 2013)
!     modified to a vertical column model 
      implicit none
      integer k,kmx
#include "size.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "npzd.h"

      real snpzd(ntnpzd), tnpzd(km,ntnpzd), src(km,nsrc)
      real dic_npzd_sms(km), twodt, rctheta, gl, impo, expo
      real npp, remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real graz_Det, graz_Z, phin, dz, prca, dprca, nud, bct, bctz
      real t_in(km), po4_in(km), sedrr, nudop, nudon, npp_dop
      real impofe, expofe, remife, nfix(km)
      real fesed, bfe, sgb, o2_in(km)
# if defined O_npzd_extra_diagnostics
      real avej, avej_D, gmax, no3P, po4P, po4_D
# endif
# if defined O_npzd_o2
      real fo2, so2, o2flag
# endif
# if defined O_npzd_iron
#  if defined O_npzd_iron_diagnostics
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
# if defined O_npzd_caco3
      real biococc,jmax_C,u_C,avej_C,gd_C,npp_C, morp_C,morpt_C
      real nuct, g_C, graz_C, ing_C, felimit_C, gl_C, npp_C_dop
      real dissl, caco3in, sil_in, atmpres1, pHlo1, pHhi1, pH1, p_in
      real co2star1, dco2star1, pCO21, dpco21, CO31, omegaca, omegaar
      real c_in, wwc,expocaco3, impocaco3, dissk1, calatt
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
      real rc13impo, rc13expo, ac13_DIC_aq
      real ac13_aq_POC, ac13b, prca13, rdic13, rtdic13(km)
      real pH, co2star, dco2star, pCO2, dpco2, CO3
      real Omega_c, Omega_a, pCO2_old, pCO2_new
      real atmpres, depth
# if defined O_npzd_caco3
      real rcaco3c13impo, rcaco3c13expo, rtcaco3c13(km)
# endif
# endif
      real sgb_in(km)
# if defined O_npzd_nitrogen
      real no3_in(km), lno3, sg_bdeni
      real morp_D, lntp, npp_D_dop
      real npp_D, graz_D, morpt_D, no3flag, wcdeni, bdeni(km)
      real din15_in(km), din15flag, eps_bdeni
      real rn15impo, rn15expo, uno3, rno3, bwcdeni, bbdeni
# endif   
# if defined O_sed
      real sedcorgflx, sedcalflx
# endif

      expo = 0.0
      impo = 0.0
      phin = 0.0                ! integrated phytoplankton
      prca = 0.0                ! integrated production of calcite

      sedrr = 0.0
# if defined O_save_npzd
      rsedrr = 0.0
      rprca = 0.0
# endif
# if defined O_npzd_nitrogen_15
      rn15impo = 0.0
      rn15expo = 0.0
# endif
# if defined O_carbon_13
      rc13impo = 0.0
      rc13expo = 0.0
      prca13 = 0.0
#  if defined O_npzd_caco3
      rcaco3c13impo = 0.0
      rcaco3c13expo = 0.0
#  endif
# endif

# if defined O_npzd_iron
      expofe = 0.0
      impofe = 0.0
# endif
# if defined O_kk_ballast
      expo_B = 0.0              ! detritus from coccs and zoop
      impo_B = 0.0              ! detritus from coccs and zoop
# endif
      calpro = 0.0
# if defined O_npzd_caco3
      caco3in = 0.0
      calatt = 0.0
      dissl = 0.0
      impocaco3 = 0.0
      expocaco3 = 0.0
      dissk1 = 0.0
# endif
      src(:,:) = 0.0
      snpzd(:) = 0.0
      rcalpro(:) = 0.0
# if defined O_npzd_caco3
      rcalatt(:) = 0.0
      rdissl(:) = 0.0
      rexpocaco3(:) = 0.0
# endif
# if defined O_save_npzd
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
#  if defined O_npzd_nitrogen
      rnpp_dop(:) = 0.0
      rnpp_D(:) = 0.0
      rnpp_D_dop(:) = 0.0
      rgraz_D(:) = 0.0
      rmorp_D(:) = 0.0
      rmorpt_D(:) = 0.0
      rnfix(:) = 0.0
#   if defined O_npzd_caco3
      rnpp_C_dop(:) = 0.0
#   endif
#  endif
#  if defined O_kk_ballast
      rexpo_B(:) = 0.0          !ballast detritus
      rgraz_Det_B(:) = 0.0
      rremi_B(:) = 0.0
#  endif
#  if defined O_npzd_caco3
      rnpp_C(:) = 0.0
      rgraz_C(:) = 0.0
      rmorp_C(:) = 0.0
      rmorpt_C(:) = 0.0
#  endif
#  if defined O_npzd_iron
      rexpofe(:) = 0.0
      rremife(:) = 0.0
#   if defined O_npzd_iron_diagnostics
      rfeorgads(:) = 0.0
      rdeffe(:) = 0.0
      rfeprime(:) = 0.0
      rfecol(:) = 0.0
#    if defined O_npzd_caco3
      rdeffe_C(:) = 0.0
#    endif         
#   endif
#  endif
#  if defined O_npzd_extra_diagnostics
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
# if defined O_npzd_nitrogen_15
         rn15impo = rn15expo
# endif
# if defined O_carbon_13 || defined O_npzd_caco3
         atmpres = 1.0          !atm
         depth = zt(k)/100.     !m

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
#   if defined O_npzd_caco3         
         rcaco3c13impo = rcaco3c13expo*dztr(k)
#   endif
#  endif
#  if defined O_npzd_caco3         
c         dissk0 =(1-(CO3-CO3/Omega_c))/(kcal+(CO3-CO3/Omega_C))
c         dissk1 = min(1.,max(0.,dissk0)) ! Andreas
         dissk1 = dissk0*(1-Omega_c) ! Andreas 11 d^-1 from Gehlen 2006 
#  endif         
# endif
         swr = swr*exp(-kc*phin
# if defined O_npzd_caco3
     &             - kc_c*caco3in
# endif
     &              )         

         phin = max(tnpzd(k,imobiphyt), trcmin)*dzt(k)
# if defined O_npzd_nitrogen
     &        + max(tnpzd(k,imobidiaz), trcmin)*dzt(k)
# endif
# if defined O_npzd_caco3
     &        + max(tnpzd(k,imobicocc), trcmin)*dzt(k)         
         caco3in = caco3in + tnpzd(k,imobicaco3)*dzt(k)         
         impocaco3 = expocaco3*dztr(k)
# endif
# if defined O_kk_ballast
         impo_B = expo_B*dztr(k) !ballasted detritus
# endif
         gl = swr*exp(ztt(k)*rctheta)
         impo = expo*dztr(k)
# if defined O_npzd_iron
         impofe = expofe*dztr(k)
# endif

         bct = bbio**(cbio*t_in(k))
# if defined O_npzd_o2
         bctz = (0.5*(tanh(o2_in(k) - 8.)+1))
     &        *bbio**(cbio*min(t_in(k),20.))
!        decrease remineralisation rate in oxygen minimum zone
         nud = nud0*(0.65+0.35*tanh(o2_in(k) - 3.0))         
# else
         bctz = bbio**(cbio*min(t_in(k),20.))
         nud = nud0         
# endif
# if defined O_npzd_nitrogen         
         nudon = nudon0
         nudop = nudop0
# endif
!-----------------------------------------------------------------------
!        call the npzd ecosystem dynamics model
!-----------------------------------------------------------------------

         call mobi_src (tnpzd(k,:), gl, bct, impo, dzt(k)
     &,                 dayfrac, wd(k), rkwz(k), nud
# if defined O_npzd_caco3
     &,                 impocaco3,wc(k),dissk1, expocaco3
# endif        
# if defined O_kk_ballast
     &,                 impo_B, expo_B, remi_B, graz_Det_B
# endif
# if defined O_npzd_nitrogen
     &,                 nudop, nudon, nfix(k)
# endif
     &,                 snpzd, expo
# if defined O_carbon
     &,                 calpro
#  if defined O_npzd_caco3
     &,                 dissl, calatt
#  endif
# endif
# if defined O_save_npzd
     &,                 npp, morpt, remi, excr, graz, morp, morz
     &,                 graz_Det, graz_Z
#  if defined O_npzd_nitrogen
     &,                 npp_D, npp_D_dop, graz_D, morp_D, morpt_D
     &,                 npp_dop
#   if defined O_npzd_caco3
     &,                 npp_C_dop
#   endif
#  endif
#  if defined O_npzd_caco3
     &,                 npp_C, morpt_C, graz_C
     &,                 morp_C
#  endif        
#  if defined O_npzd_extra_diagnostics
     &,                 avej, avej_D, gmax, no3P
     &,                 po4P, po4_D
#  endif
# endif
     &,                 bctz
# if defined O_npzd_nitrogen_15
     &,                 rn15impo, rn15expo
# endif
# if defined O_carbon_13
     &,                 rc13impo, rc13expo, ac13b
#  if defined O_npzd_caco3
     &,                 rcaco3c13impo, rcaco3c13expo
#  endif
# endif
# if defined O_npzd_iron
     &,                    expofe, impofe, remife
#  if defined O_npzd_iron_diagnostics
     &,                    feorgads, thetamax, deffe, feprime
     &,                    fecol
#   if defined O_npzd_caco3
     &,                    deffe_C, thetamax_C
#   endif         
#  endif
c     juan: use o2_in(k) instead of t(i,k,j,io2,taum1)      
# endif         
# if defined O_npzd_o2
     &,                    o2_in(k)            
# endif
     &                      )
! These are source/sink terms 
         snpzd = snpzd*rdtts
# if defined O_npzd_iron
         expofe = expofe*rnbio
# endif
# if defined O_npzd_caco3
         expocaco3 = expocaco3*rnbio
# endif
# if defined O_kk_ballast
         expo_B = expo_B*rnbio !ballast detritus
# endif
         expo = expo*rnbio
# if defined O_npzd_nitrogen_15
         rn15expo = rn15expo*rnbio
# endif
# if defined O_carbon_13
         rc13expo = rc13expo*rnbio
#  if defined O_npzd_caco3
         rcaco3c13expo = rcaco3c13expo*rnbio
#  endif
# endif
         rcalpro(k) = calpro*rnbio
# if defined O_npzd_caco3
         rcalatt(k) = calatt*rnbio
         rdissl(k) = dissl*rnbio ! AS changed Aug 17, 2017
         rexpocaco3(k) = expocaco3
# endif
# if defined O_save_npzd
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
#  if defined O_npzd_nitrogen
         rnpp_dop(k) = npp_dop*rnbio         
         rnpp_D(k) = npp_D*rnbio
         rnpp_D_dop(k) = npp_D_dop*rnbio
         rgraz_D(k) = graz_D*rnbio
         rmorp_D(k) = morp_D*rnbio
         rmorpt_D(k) = morpt_D*rnbio
         rnfix(k) = nfix(k)*rnbio
#   if defined O_npzd_caco3
         rnpp_C_dop(k) = npp_C_dop*rnbio
#   endif
#  endif
#  if defined O_kk_ballast
         rexpo_B(k) = expo_B  !ballast detritus
         rgraz_Det_B(k) = graz_Det_B*rnbio
         rremi_B(k) = remi_B*rnbio
#  endif
#  if defined O_npzd_caco3
         rnpp_C(k) = npp_C*rnbio
         rgraz_C(k) = graz_C*rnbio
         rmorp_C(k) = morp_C*rnbio
         rmorpt_C(k) = morpt_C*rnbio
#  endif
#  if defined O_npzd_iron
         rexpofe(k) = expofe
         rremife(k) = remife*rnbio
#   if defined O_npzd_iron_diagnostics
         rfeorgads(k) = feorgads*rnbio
         rdeffe(k) = deffe*rnbio
         rfeprime(k) = feprime*rnbio
         rfecol(k) = fecol*rnbio
#    if defined O_npzd_caco3
         rdeffe_C(k) = deffe_C*rnbio
#    endif         
#   endif
#  endif
#  if defined O_npzd_extra_diagnostics
         ravej(k) = avej*rnbio
         ravej_D(k) = avej_D*rnbio
         rgmax(k) = gmax*rnbio
         rno3P(k) = no3P*rnbio
         rpo4P(k) = po4P*rnbio
         rpo4_D(k) = po4_D*rnbio
#  endif
# endif
# if defined O_npzd_caco3
!----------------------------------------------------------------------
!             calculate calcite at the bottom and dissolve
!---------------------------------------------------------------------
#  if defined O_sed
                if (addflxo .and. eots)
     &           sbc(i,jrow,ircal) = sbc(i,jrow,ircal)
     &           + expocaco3*twodt(1)*dzt(k)*1.e-3 !mol C
#  endif
# endif

!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

# if defined O_sed
         if (k .eq. kmx) then
#  if defined O_kk_ballast                   
            sedcorgflx = sedcorgflx + (expo+expo_B)*redctn*twodt*dzt(k)
#  else
            sedcorgflx = sedcorgflx + expo*redctn*twodt*dzt(k)
#  endif
         endif
# endif
# if defined O_kk_ballast
         sedrr = sgb_in(k)*(expo+expo_B)*dzt(k)
# else               
         sedrr = sgb_in(k)*expo*dzt(k)
# endif
# if defined O_npzd_nitrogen
!-----------------------------------------------------------------------
!        benthic denitrification model of Bohlen et al., 2012, GBC (eq. 10)
!        NO3 is removed out of bottom water nitrate.
!        See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!        limit denitrification as nitrate approaches 0 uM
         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
#  if defined O_npzd_nitrogen_15
         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
#  else
         din15flag = 1.
#  endif
         lno3 = 0.5*tanh(no3_in(k)*10 - 5.0)

         sg_bdeni = (0.06 + 0.19*0.99**(max(o2_in(k),trcmin)             
     &                 - max(no3_in(k),trcmin)))
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
#  if defined O_save_npzd
         rbdeni(k) = bdeni(k)
#  endif
#  if defined O_npzd_nitrogen_15
         rno3 = max(din15_in(k), trcmin*rn15std/(1+rn15std))
     &            /max(no3_in(k)-din15_in(k),trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         eps_bdeni = eps_bdeni0*exp(-2.5e-6*(zt(k)))
         bbdeni = rno3 - eps_bdeni*rno3/1000.

         snpzd(imobidin15) = snpzd(imobidin15) 
     &                     + rn15expo*sgb_in(k)*expo 
     &                     - bbdeni/(1+bbdeni)*sg_bdeni
#  endif              
# endif
# if defined O_npzd_iron
         fesed=fetopsed*bct*redptn*expo*sgb_in(k)
         bfe=fesed
         snpzd(imobidfe) = snpzd(imobidfe) + fesed                
         o2flag = 0.5 - sign(0.5, o2_in(k)-o2min)
         snpzd(imobidfe) = snpzd(imobidfe) + expofe*sgb_in(k)*o2flag
         bfe = bfe + expofe*sgb_in(k)*o2flag
         expofe = expofe - sgb_in(k)*expofe*o2flag
#  if defined O_save_npzd
         rremife(k) = rremife(k) + bfe
#  endif
#  if defined O_npzd_iron_diagnostics
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
# if defined O_save_npzd
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
# if defined O_npzd_nitrogen
         src(k,isdop) = snpzd(imobidop)   
         src(k,isno3) = snpzd(imobino3)
         src(k,isdon) = snpzd(imobidon)
         src(k,isdiaz) = snpzd(imobidiaz)
#  if defined O_npzd_nitrogen_15
         src(k,isdin15) = snpzd(imobidin15)
         src(k,isdon15) = snpzd(imobidon15)
         src(k,isphytn15) = snpzd(imobiphytn15)
         src(k,iszoopn15) = snpzd(imobizoopn15)
         src(k,isdetrn15) = snpzd(imobidetrn15)
         src(k,isdiazn15) = snpzd(imobidiazn15)
#   if defined O_npzd_caco3
         src(k,iscoccn15) = snpzd(imobicoccn15)
#   endif         
#  endif
# endif
# if defined O_carbon_13
         src(k,isdic13) = snpzd(imobidic13)
         src(k,isphytc13) = snpzd(imobiphytc13)
         src(k,iszoopc13) = snpzd(imobizoopc13)
         src(k,isdetrc13) = snpzd(imobidetrc13)
#  if defined O_npzd_nitrogen
         src(k,isdoc13) = snpzd(imobidoc13)
         src(k,isdiazc13) = snpzd(imobidiazc13)
#  endif         
#  if defined O_npzd_caco3
         src(k,iscoccc13) = snpzd(imobicoccc13)
         src(k,iscaco3c13) = snpzd(imobicaco3c13)
#  endif
# endif
# if defined O_npzd_caco3
         src(k,iscocc) = snpzd(imobicocc)
         src(k,iscaco3) = snpzd(imobicaco3)
# endif
# if defined O_kk_ballast
         src(k,isdetr_B) = snpzd(imobidetr_B)
# endif
# if defined O_npzd_iron
         src(k,isdfe) = snpzd(imobidfe)
         src(k,isdetrfe) = snpzd(imobidetrfe)
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

#  if !defined O_npzd_caco3         
         src(k,isdic) = snpzd(imobidic) - dprca             
#  endif         
#  if defined O_carbon_13
         rtdic13(k) = max(dic13_in(k),trcmin*rc13std/(1+rc13std))
     &             / max(dic_in(k),trcmin)
         rtdic13(k) = min(rtdic13(k), 2.*rc13std/(1+rc13std))
         rtdic13(k) = max(rtdic13(k), 0.5*rc13std/(1+rc13std))

         prca13 = prca13 + dprca*dzt(k)*rtdic13(k)
#   if defined O_npzd_caco3          
         rtcaco3c13(k) = max(tnpzd(k,imobicaco3c13)
     &                     ,trcmin*rc13std/(1+rc13std))
     &             / max(tnpzd(k,imobicaco3),trcmin)
         rtcaco3c13(k) = min(rtcaco3c13(k), 2.*rc13std/(1+rc13std))
         rtcaco3c13(k) = max(rtcaco3c13(k), 0.5*rc13std/(1+rc13std))
#   endif
#   if !defined O_npzd_caco3          
         src(k,isdic13) = src(k,isdic13) - rtdic13(k)*dprca
#   endif
#  endif
# endif       
# if defined O_npzd_alk
         src(k,isalk) = -snpzd(imobidic)*redntc*1.e-3
#  if !defined O_npzd_caco3
     &        - 2.*dprca
#  endif
# endif
!        calculate total export to get total import for next layer
         expo = expo*dzt(k)
# if defined O_carbon_13
         rc13expo = rc13expo*dzt(k)
#  if defined O_npzd_caco3
         rcaco3c13expo = rcaco3c13expo*dzt(k)
#  endif         
# endif         
# if defined O_kk_ballast
         expo_B = expo_B*dzt(k)
# endif         
# if defined O_npzd_iron
         expofe = expofe*dzt(k)
# endif
# if defined O_npzd_caco3
         expocaco3 = expocaco3*dzt(k)
# endif
      enddo
!1111 end of main k-loop
# if defined O_save_npzd
      rprca = prca
# endif

# if defined O_npzd_o2
!2222 start of second k-loop
      do k=1,kmx
!        limit oxygen consumption below concentrations of
!        5umol/kg as recommended in OCMIP
         fo2 = 0.5*tanh(o2_in(k) - 2.5)
!        sink of oxygen
         so2 = dic_npzd_sms(k)*redotc
#  if defined O_npzd_nitrogen
! O2 is needed to generate the equivalent of NO3 from N2 during N2 fixation
! 0.5 H2O + 0.5 N2+1.25O2 -> HNO3
! note that so2 is -dO2/dt
     &          + nfix(k)*rnbio*1.25e-3
!        add denitrification as source term for NO3
         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
#   if defined O_npzd_nitrogen_15
         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
#   else
         din15flag = 1.
#   endif
         lno3 = 0.5*tanh(no3_in(k) - 2.5)
         lntp = 0.5*tanh(no3_in(k)/(redntp*po4_in(k))*100 - 60.)
!        800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
         wcdeni = 800.*no3flag*so2*(0.5 - fo2)*(0.5 + lno3)
#   if defined O_npzd_nitrogen_15
     &           *din15flag
#   endif
         wcdeni = max(wcdeni, 0.)
         src(k,isno3) = src(k,isno3) - wcdeni
#   if defined O_npzd_nitrogen_15
!        calculate isotope effect of water column denitrification
         uno3 = wcdeni*twodt/no3_in(k)
         uno3 = min(uno3, 0.999)
         uno3 = max(uno3, trcmin)
         rno3 = max(din15_in(k),trcmin*rn15std/(1+rn15std))
     &             / max(no3_in(k)-din15_in(k),
     &                   trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         bwcdeni = rno3 + eps_wcdeni*(1-uno3)/uno3*log(1-uno3)*rno3/1000.

         src(k,isdin15) = src(k,isdin15) - (bwcdeni/(1+bwcdeni))*wcdeni
#   endif
#   if defined O_npzd_alk
! Correct the ALK stoichiometry to account for N2 fixation
         src(k,isalk) = src(k,isalk) + wcdeni*1.e-3
!bdeni
         src(k,isalk) = src(k,isalk) + bdeni(k)*1.e-3
! Now add account for N2 fixation (ALK production is tied to PO4 change and
! thus in the case of N2 fixation is not correct).
         src(k,isalk) = src(k,isalk) - nfix(k)*rnbio*1.e-3
#   endif
#   if defined O_save_npzd
         rwcdeni(k) = wcdeni
#   endif
#  endif
         src(k,iso2) = -so2*(0.5 + fo2)
      enddo
!2222 end of second k-loop
# endif

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
!3333 start of third k-loop
      do k=1,kmx-1
# if defined O_carbon
#  if defined O_npzd_caco3
         src(k,isdic) = src(k,isdic)
     &        + rdissl(k)*1.e-3
     &        - rcalpro(k)*1.e-3
     &        - rcalatt(k)*1.e-3
#  else
         src(k,isdic) = src(k,isdic) + prca*rcak(k)
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
         src(k,isdic13) = src(k,isdic13)
     &        + rdissl(k)*1.e-3*rtcaco3c13(k) ! AS
     &        - rcalpro(k)*1.e-3*rtdic13(k)
     &        - rcalatt(k)*1.e-3*rtdic13(k)
#   else              
         src(k,isdic13) = src(k,isdic13) + prca13*rcak(k)
#   endif
#  endif
# endif
# if defined O_npzd_alk
#   if defined O_npzd_caco3
         src(k,isalk) = src(k,isalk)
     &        + 2.*rdissl(k)*1.e-3
     &        - 2.*rcalpro(k)*1.e-3
     &        - 2.*rcalatt(k)*1.e-3
#   else
         src(k,isalk) = src(k,isalk) + 2.*prca*rcak(k)
#   endif
# endif
      enddo
!3333 end of third k-loop
# if defined O_sed
      sedcalflx = sedcalflx + (prca*rcab(kmx)
     &     - prca*rcak(kmx))*twodt*dzt(kmx)
# endif
# if defined O_carbon
#  if defined O_npzd_caco3
      src(kmx,isdic) = src(kmx,isdic)
     &     + rdissl(kmx)*1.e-3
     &     - rcalpro(kmx)*1.e-3
     &     - rcalatt(kmx)*1.e-3
     &     + rexpocaco3(kmx)*1.e-3
#  else
      src(kmx,isdic) = src(kmx,isdic) + prca*rcab(kmx)
#  endif
#  if defined O_carbon_13
#   if defined O_npzd_caco3
      src(kmx,isdic13) = src(kmx,isdic13)
     &                 + rdissl(kmx)*1.e-3*rtcaco3c13(kmx)
     &                 - rcalpro(kmx)*1.e-3*rtdic13(kmx)
     &                 - rcalatt(kmx)*1.e-3*rtdic13(kmx)
     &                 + rexpocaco3(kmx)*1.e-3*rtcaco3c13(kmx)       
#   else  
      src(kmx,isdic13) = src(kmx,isdic13) + prca13*rcab(kmx)
#   endif
#  endif
# endif
# if defined O_npzd_alk
#  if defined O_npzd_caco3
      src(kmx,isalk) = src(kmx,isalk)
     &               + 2.*rdissl(kmx)*1.e-3
     &               - 2.*rcalpro(kmx)*1.e-3
     &               - 2.*rcalatt(kmx)*1.e-3
     &               + 2.*rexpocaco3(kmx)*1.e-3
#  else
      src(kmx,isalk) = src(kmx,isalk) + 2.*prca*rcab(kmx)
#  endif
# endif
      return
      end
!     END mobi_driver

      subroutine mobi_src (bioin, gl, bct, impo, dzt
     &,                    dayfrac, wwd, rkw, nud
# if defined O_npzd_caco3
     &,                    impocaco3, wwc, dissk1, expocaco3out
# endif
# if defined O_kk_ballast
     &,                    impo_B, expoout_B, remiout_B, graz_Det_out_B
# endif
# if defined O_npzd_nitrogen
     &,                    nudop, nudon, nfixout
# endif      
     &,                    bioout, expoout
# if defined O_carbon
     &,                    calproout
#  if defined O_npzd_caco3
     &,                    disslout, calattout
#  endif      
# endif      
# if defined O_save_npzd
     &,                    nppout, morptout, remiout, excrout, grazout
     &,                    morpout, morzout, graz_Det_out, graz_Zout
#  if defined O_npzd_nitrogen
     &,                    npp_Dout, npp_D_dopout, graz_Dout, morp_Dout
     &,                    morpt_Dout, npp_dopout
#   if defined O_npzd_caco3
     &,                    npp_C_dopout
#   endif
#  endif
#  if defined O_npzd_caco3
     &,                    npp_Cout, morpt_Cout,graz_Cout
     &,                    morp_Cout
#  endif
#  if defined O_npzd_extra_diagnostics
     &,                    avej_out, avej_D_out, gmax_out, no3P_out
     &,                    po4P_out, po4_D_out
#  endif
# endif
     &,                    bctz
# if defined O_npzd_nitrogen_15
     &,                    rn15impo, rn15expoout
# endif
# if defined O_carbon_13
     &,                    rc13impo, rc13expoout, ac13b
#  if defined O_npzd_caco3
     &,                    rcaco3c13impo, rcaco3c13expoout
#  endif
# endif
# if defined O_npzd_iron
     &,                    expofeout, impofe, remifeout
#  if defined O_npzd_iron_diagnostics
     &,                    feorgadsout, thetamaxout, deffeout
     &,                    feprimeout, fecolout
#   if defined O_npzd_caco3
     &,                    deffeout_C, thetamaxout_C
#   endif  
#  endif
# endif       
# if defined O_npzd_o2
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
!
!     Levin Nickelsen has written the iron model:
!     Nickelsen et al. 2015, Geosc. Model Dev. 8, 1357-1381, doi:10.5194/gmd-8-1357-2015
!
!     Juan Muglia has modified the subgrid-scale scheme used for benthic denitrification
!     and iron release from sediments. He has also included geothermal sources of iron.
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
!       rkw        = reciprical of kw*dzt(k)
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
      real dzt, po4flag, phytflag, zoopflag, detrflag, wwd, rkw, gd
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

      real biococc,jmax_C,u_C,avej_C,gd_C,npp_C,npp_Cout,morp_C,morpt_C
      real nuct,g_C,coccflag,graz_C,graz_Cout,morp_Cout,coccn15flag
      real morpt_Cout, ing_C, gl_C, limP_C, dig_C, excr_C, sf_C
      real dopupt_C_flag, dopupt_C, npp_C_dopout, biococcn15,rtcoccn15
      real biococcc13, coccc13flag, rtcoccc13, biocaco3c13, caco3c13flag
      real rtcaco3c13

      real caco3flag, wwc,expocaco3
      real expocaco3out,calpro,calproout
      real biocaco3, impocaco3, dissk1, dissl, disslout
      real calatt,calattout
      
      real biodfe, biodetrfe, dfeflag, detrfeflag, o2flag
      real expofe, impofe, feorgads, remife, thetamax, deffe, fepa
      real thetamax_D, deffe_D, thetachl, thetachl_D, chl, feprime
      real fecol, irrtop, kirr, aveirr, fesed, gl_O, gl_D, alpha_O
      real alpha_D, par, deffe_C, thetamax_C, thetamaxout_C, deffeout_C
      real alpha_C

      real expofeout, impofeout, feorgadsout, fecolout
      real remifeout, thetamaxout, deffeout, thetachlout
      real chlout, feprimeout, o2, chl_D_out

      real p1,p2,kfevar, p1_C, p2_C, kfevar_C

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "npzd.h"

      real bioin(ntnpzd), bioout(ntnpzd)

# if defined O_npzd_iron
      p1=min(bioin(imobiphyt),pmax) !phytoplankton
      p2=max(0.0,bioin(imobiphyt)-pmax)
      kfevar=(kfemin*p1+kfemax*p2)/(p1+p2) 
      deffe=bioin(imobidfe)/(kfevar+bioin(imobidfe))
      thetamax=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe
      alpha_O=alphamin+(alphamax-alphamin)*deffe
      gl_O=gl*thetamax*alpha_O
#  if defined O_npzd_caco3
      p1_C=min(bioin(imobicocc),pmax_C) !calcifiers
      p2_C=max(0.0,bioin(imobicocc)-pmax_C)
      kfevar_C=(kfemin_C*p1_C+kfemax_C*p2_C)/(p1_C+p2_C) 
      deffe_C=bioin(imobidfe)/(kfevar_C+bioin(imobidfe))
      thetamax_C=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe_C
      alpha_C=alphamin+(alphamax-alphamin)*deffe_C
      gl_C=gl*thetamax_C*alpha_O
#  endif
#  if defined O_npzd_nitrogen
      deffe_D=bioin(imobidfe)/(kfe_D+bioin(imobidfe))
      thetamax_D=thetamaxlo+(thetamaxhi-thetamaxlo)*deffe_D
      alpha_D=alphamin+(alphamax-alphamin)*deffe_D
      gl_D=gl*thetamax_D*alpha_D
#  endif
# endif

      bioout(:) = 0.0
      biopo4 = bioin(imobipo4)
      biophyt = bioin(imobiphyt)
      biozoop = bioin(imobizoop)
      biodetr = bioin(imobidetr)
# if defined O_carbon      
      biodic = bioin(imobidic)
# endif
# if defined O_npzd_nitrogen
      biodop = bioin(imobidop)
      biono3 = bioin(imobino3)
      biodon = bioin(imobidon)
      biodiaz = bioin(imobidiaz)
#  if defined O_npzd_nitrogen_15
      biodin15 = bioin(imobidin15)
      biodon15 = bioin(imobidon15)
      biophytn15 = bioin(imobiphytn15)
      biozoopn15 = bioin(imobizoopn15)
      biodetrn15 = bioin(imobidetrn15)
      biodiazn15 = bioin(imobidiazn15)
#   if defined O_npzd_caco3
      biococcn15 = bioin(imobicoccn15)
#   endif      
#  endif
# endif
# if defined O_carbon_13
      biodic13 = bioin(imobidic13)
      biophytc13 = bioin(imobiphytc13)
      biozoopc13 = bioin(imobizoopc13)
      biodetrc13 = bioin(imobidetrc13)
#  if defined O_npzd_nitrogen
      biodoc13 = bioin(imobidoc13)
      biodiazc13 = bioin(imobidiazc13)
#  endif
#  if defined O_npzd_caco3
      biococcc13 = bioin(imobicoccc13)
      biocaco3c13 = bioin(imobicaco3c13)
#  endif
# endif
# if defined O_npzd_caco3
      biococc = bioin(imobicocc)
      biocaco3 = bioin(imobicaco3)
# endif
# if defined O_kk_ballast
      biod_B = bioin(imobidetr_B)
# endif
# if defined O_npzd_iron
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
      coccflag = 1.
      detrfeflag = 1.
      coccn15flag = 1.
      coccc13flag = 1.      
      caco3c13flag = 1.      
# if defined O_npzd_nitrogen
      dopflag = 0.5 + sign(0.5,biodop - trcmin)      
      no3flag = 0.5 + sign(0.5,biono3 - trcmin)
      donflag = 0.5 + sign(0.5,biodon - trcmin)
      diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_npzd_nitrogen_15
      din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
      don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
      phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
#   if defined O_npzd_caco3
      coccn15flag = 0.5 + sign(0.5, biococcn15 - trcmin)
#   endif      
      zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
      detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
      diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_carbon_13
      dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
      phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
#  if defined O_npzd_caco3
      coccc13flag = 0.5 + sign(0.5,biococcc13 - trcmin)
      caco3c13flag = 0.5 + sign(0.5,biocaco3c13 - trcmin)
#  endif      
      zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
      detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
#  if defined O_npzd_nitrogen
      doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
      diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
#  endif
# endif
# if defined O_npzd_iron
      dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
      detrfeflag = 0.5 + sign(0.5,biodetrfe - trcmin)
# endif
# if defined O_kk_ballast
      dflag_B = 0.5 + sign(0.5,biod_B - trcmin)
# endif
# if defined O_npzd_caco3
      coccflag = 0.5 + sign(0.5,biococc - trcmin)
      caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)
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
# if defined O_npzd_nitrogen
      biono3 = max(biono3, trcmin)
      biodop = max(biodop, trcmin)      
      biodon = max(biodon, trcmin)
      biodiaz = max(biodiaz, trcmin)
#  if defined O_npzd_nitrogen_15
      biodin15 = max(biodin15, trcmin)
      biodon15 = max(biodon15, trcmin)
      biophytn15 = max(biophytn15, trcmin)
#   if defined O_npzd_caco3
      biococcn15 = max(biococcn15, trcmin)
#   endif      
      biozoopn15 = max(biozoopn15, trcmin)
      biodetrn15 = max(biodetrn15, trcmin)
      biodiazn15 = max(biodiazn15, trcmin)
#  endif
# endif
# if defined O_carbon_13
      biodic13 = max(biodic13, trcmin)
      biophytc13 = max(biophytc13, trcmin)
#   if defined O_npzd_caco3
      biococcc13 = max(biococcc13, trcmin)
      biocaco3c13 = max(biocaco3c13, trcmin)
#   endif 
      biozoopc13 = max(biozoopc13, trcmin)
      biodetrc13 = max(biodetrc13, trcmin)    
#   if defined O_npzd_nitrogen
      biodoc13 = max(biodoc13, trcmin)
      biodiazc13 = max(biodiazc13, trcmin)
#   endif 
# endif
# if defined O_kk_ballast
      biod_B = max(biod_B, trcmin)
# endif
# if defined O_npzd_caco3
      biococc = max(biococc, trcmin)
      biocaco3 = max(biocaco3, trcmin)      
# endif
# if defined O_npzd_iron
      biodfe = max(biodfe, trcmin)
      biodetrfe = max(biodetrfe, trcmin)
# endif


!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      f1 = exp((-kw - kc*(biophyt
# if defined O_npzd_nitrogen
     &     + biodiaz
# endif
# if defined O_npzd_caco3
     &     + biococc
# endif
     &     )      
# if defined O_npzd_caco3
     &     - kc_c*biocaco3
# endif
     &     )*dzt)

# if defined O_npzd_iron
      jmax = abio_P*bct*deffe
# else
      jmax = abio_P*bct
# endif
      gd = jmax*dayfrac
# if defined O_npzd_iron
      u1 = max(gl_O/gd,1.e-6)
# else
      u1 = max(gl/gd,1.e-6)
# endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2

      avej = gd*(phi1 - phi2)/((kw+kc*(biophyt
# if defined O_npzd_nitrogen
     &     + biodiaz  
# endif      
# if defined O_npzd_caco3
     &     + biococc
# endif
     &     )          
# if defined O_npzd_caco3
     &     + kc_c*biocaco3
# endif
     &       )*dzt)
!     Make the max grazing rate a function of temperature
!     bctz sets an upper limit on the effects of temp on grazing
!     in contrast to phytoplankton growth rates bct, which are unlimited
      gmax = gbio*bctz
# if defined O_npzd_nitrogen
#  if defined O_npzd_iron
      jmax_D = max(0.,abio_P*(bct - 2.6)*deffe_D)*jdiar
#  else
      jmax_D = max(0.,abio_P*(bct - 2.6))*jdiar
#  endif
!
      gd_D = max(1.e-14,jmax_D*dayfrac)
#  if defined O_npzd_iron
      u1 = max(gl_D/gd_D,1.e-6)
#  else
      u1 = max(gl/gd_D,1.e-6)
#  endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*(phi1 - phi2)/((kw+kc*(biophyt+biodiaz
#  if defined O_npzd_caco3
     &       + biococc
#  endif
     &        )          
#  if defined O_npzd_caco3
     &       + kc_c*biocaco3
#  endif
     &         )*dzt)
# endif
# if defined O_npzd_caco3
#  if defined O_npzd_iron
      jmax_C = abio_C*bct*deffe_C
#  else
      jmax_C = abio_C*bct
#  endif
      gd_C = jmax_C*dayfrac
#  if defined O_npzd_iron
      u1 = max(gl_C/gd_C,1.e-6)
#  else
      u1 = max(gl/gd_C,1.e-6)
#  endif
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_C = gd_C*(phi1 - phi2)/((kw+kc*(biophyt
#  if defined O_npzd_nitrogen
     &       +  biodiaz
#  endif
     &       +  biococc)
#  if defined O_npzd_caco3
     &       + kc_c*biocaco3
#  endif
     &       )*dzt)
# endif
      nupt = nupt0*bct
# if defined O_npzd_nitrogen
      nupt_D = nupt0_D*bct
      nfixout = 0.0
# endif
      expoout = 0.0
      grazout = 0.0
      morpout = 0.0
      morzout = 0.0
      graz_Det_out = 0.0
      graz_Zout = 0.0
# if defined O_npzd_nitrogen_15
      rn15expoout = 0.0
# endif
# if defined O_carbon_13
      rc13expoout = 0.0
#  if defined O_npzd_caco3
      rcaco3c13expoout = 0.0
#  endif
# endif
      calproout = 0.0
# if defined O_npzd_caco3
      nuct = nuct0*bct
      npp_Cout = 0.0
      graz_Cout = 0.0
      morp_Cout = 0.0
      morpt_Cout = 0.0
      calattout = 0.0
      disslout = 0.0
      expocaco3out = 0.0
# endif
# if defined O_kk_ballast
      expoout_B = 0.0
      graz_Det_out_B = 0.0
      remiout_B = 0.0
# endif
# if defined O_npzd_iron
      expofeout = 0.0
      remifeout = 0.0
# endif
# if defined O_save_npzd
      nppout = 0.0
      npp_dopout = 0.0
      morptout = 0.0
      remiout = 0.0
      excrout = 0.0
#  if defined O_npzd_nitrogen
      npp_Dout = 0.0
      npp_D_dopout = 0.0
      graz_Dout = 0.0
      morpt_Dout = 0.0
      morp_Dout = 0.0
#    if defined O_npzd_caco3
      npp_C_dopout = 0.0
#    endif
#  endif
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics
      feorgadsout = 0.0
      deffeout = 0.0
      feprimeout = 0.0
      fecolout = 0.0
#    if defined O_npzd_caco3
      deffeout_C = 0.0
#    endif      
#   endif
#  endif
#  if defined O_npzd_extra_diagnostics
      avej_out = 0.0
      avej_D_out = 0.0
      gmax_out = 0.0
      no3P_out = 0.0
      po4P_out = 0.0
      po4_D_out = 0.0
#  endif
# endif

      do n=1,nbio

# if defined O_npzd_iron
        p1 = min(biophyt,pmax)
        p2 = max(0.0,biophyt - pmax)
        kfevar = (kfemin*p1+ kfemax*p2)/(p1 + p2)
        deffe = biodfe/(kfevar + biodfe)
        jmax = abio_P*bct*deffe
#  if defined O_npzd_caco3
        p1_C = min(biococc,pmax_C)
        p2_C = max(0.0,biococc - pmax_C)
        kfevar_C = (kfemin_C*p1_C + kfemax_C*p2_C)/(p1_C+p2_C) 
        deffe_C = biodfe/(kfevar_C + biodfe)
        jmax_C = abio_C*bct*deffe_C 
#  endif
#  if defined O_npzd_nitrogen
        deffe_D = biodfe/(kfe_D + biodfe)
        jmax_D = max(0.,abio_P*(bct - 2.6)*deffe_D)*jdiar
#  endif
# endif 
# if defined O_npzd_nitrogen
!       growth rate of phytoplankton
!       consume DOP when it is more efficient
        limP_dop = hdop*biodop/(k1p_P + biodop)
        limP_po4 = biopo4/(k1p_P + biopo4)
        dopupt_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP = limP_dop*dopupt_flag + limP_po4*(1.-dopupt_flag)
        u_P = min(avej, jmax*limP)
#  if defined O_npzd_caco3
        limP_dop = hdop*biodop/(k1p_C + biodop)
        limP_po4 = biopo4/(k1p_C + biopo4)
        dopupt_C_flag = 0.5 + sign(0.5, limP_dop - limP_po4)
        limP_C = limP_dop*dopupt_C_flag + limP_po4*(1.-dopupt_C_flag)
        u_C = min(avej, jmax_C*limP_C)
#  endif        
# else
        limP = biopo4/(k1p_P + biopo4)
        u_P = min(avej, jmax*limP)        
#  if defined O_npzd_caco3
        limP_C = biopo4/(k1p_C + biopo4)
        u_C = min(avej_C, jmax_C*limP_C)
#  endif
# endif        
        po4P = jmax*biopo4/(k1p_P + biopo4) 

# if defined O_npzd_nitrogen
!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
#  if defined O_npzd_caco3
        u_C = min(u_C, jmax_C*biono3/(k1n_C + biono3))
#  endif
        no3P = jmax*biono3/(k1n + biono3)
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        u_D = min(avej_D, jmax_D*limP)
        dopupt_D_flag = dopupt_flag
        po4_D = jmax_D*biopo4/(k1p_P + biopo4)
!       Set the grazing coefficients for the N case
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop
     &         + zprefDiaz*biodiaz + kzoo
#  if defined O_npzd_caco3
     &         + zprefC*biococc
#  endif
#  if defined O_kk_ballast
     &         + zprefDet*biod_B
#  endif
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
        ing_D = zprefDiaz/thetaZ
#  if defined O_npzd_caco3
        ing_C = zprefC/thetaZ
#  endif
# else
!       If "else" then set the grazing coefficients for the no N case
!       note kzoo is in terms of N so convert to P
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop 
     &         + kzoo*redptn
#  if defined O_npzd_caco3
     &         + zprefC*biococc
#  endif
#  if defined O_kk_ballast
     &         + zprefDet*biod_B
#  endif
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
#  if defined O_npzd_caco3
        ing_C = zprefC/thetaZ
#  endif        
# endif
        npp = u_P*biophyt
# if defined O_npzd_caco3
        npp_C = u_C*biococc
# endif
# if defined O_npzd_nitrogen
        dopupt = npp*dopupt_flag        
#  if defined O_npzd_caco3
        dopupt_C = npp_C*dopupt_C_flag         
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
# if defined O_npzd_nitrogen        
        recy_don = nudon*bct*biodon
        recy_dop = nudop*bct*biodop
# endif        
        morz = nuz*biozoop*biozoop
        remi = nud*bct*biodetr
        expo = wwd*biodetr
# if defined O_npzd_caco3
        g_C =  gmax*ing_C*biococc
        graz_C = g_C*biozoop
        morp_C = nuc*biococc
        morpt_C = nuct*biococc
        dissl = biocaco3*dissk1
        expocaco3 = wwc*biocaco3
# endif
# if defined O_kk_ballast
        g_Det_B = gmax*ing_Det*biod_B
        graz_Det_B = g_Det_B*biozoop
        remi_B = max(bapr*dissl/(capr*redctn*1.e3), 0.) !mmol N/m3
        expo_B = wwc*biod_B !ballasted by CaCO3
# endif
# if defined O_npzd_iron
!   Calculation of average light in mixed layer for calculation of Chl diagnostic
        par = 0.43  ! fraction of photosythetically active radiation
        irrtop=gl/2/par
        kirr=-kw - kc*(bioin(imobiphyt)
#  if defined O_npzd_nitrogen        
     &      + bioin(imobidiaz)
#  endif
#  if defined O_npzd_caco3
     &      + bioin(imobicocc)
#  endif
     &       )  
#  if defined O_npzd_caco3
     &      - kc_c*bioin(imobicaco3)
#  endif        
    
        aveirr=-1/dzt/kirr*(irrtop-irrtop*exp(kirr*dzt))
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
!

        graz = graz*phytflag*zoopflag
# if defined O_npzd_nitrogen_15
     &         *phytn15flag*zoopn15flag       
# endif
        graz_Z = graz_Z*zoopflag
# if defined O_npzd_nitrogen_15
     &          *zoopn15flag
# endif
        graz_Det = graz_Det*detrflag*zoopflag
# if defined O_npzd_nitrogen_15
     &            *detrn15flag*zoopn15flag
# endif
        morp = morp*phytflag
# if defined O_npzd_nitrogen_15
     &        *phytn15flag
# endif
        morpt = morpt*phytflag
# if defined O_npzd_nitrogen_15
     &        *phytn15flag
# endif
        morz = morz*zoopflag
# if defined O_npzd_nitrogen_15
     &        *zoopn15flag
# endif
        remi = remi*detrflag
# if defined O_npzd_nitrogen_15
     &        *detrn15flag
# endif
        expo = expo*detrflag
# if defined O_npzd_nitrogen_15
     &        *detrn15flag
# endif
# if defined O_npzd_nitrogen
        recy_dop = recy_dop*dopflag               
           npp = npp*no3flag*(dopupt_flag*dopflag
     &         + (1.-dopupt_flag)*po4flag)
#  if defined O_npzd_nitrogen_15
     &           *din15flag
#  endif
#  if defined O_npzd_caco3
           npp_C = npp_C*no3flag*(dopupt_C_flag*dopflag
     &           + (1.-dopupt_C_flag)*po4flag)
#   if defined O_npzd_nitrogen_15
     &             *din15flag
#   endif
#  endif
           npp_D = npp_D*(dopupt_D_flag*dopflag
     &           + (1.-dopupt_D_flag)*po4flag)
#   if defined O_npzd_nitrogen_15
     &             *din15flag
#   endif
        graz_D = graz_D*diazflag*zoopflag
#  if defined O_npzd_nitrogen_15
     &          *diazn15flag*zoopn15flag
#  endif
        morpt_D = morpt_D*diazflag
#  if defined O_npzd_nitrogen_15
     &           *diazn15flag
#  endif
        morp_D = morp_D*diazflag
#  if defined O_npzd_nitrogen_15
     &           *diazn15flag
#  endif
        no3upt_D = no3upt_D*no3flag
#  if defined O_npzd_nitrogen_15
     &            *din15flag
#  endif
        recy_don = recy_don*donflag
#  if defined O_npzd_nitrogen_15
     &           *don15flag
#  endif
# else
        npp = npp*po4flag
#  if defined O_npzd_caco3
        npp_C = npp_C*po4flag
#  endif
# endif
# if defined O_kk_ballast
        graz_Det_B = graz_Det_B*dflag_B*zoopflag
#  if defined O_npzd_nitrogen_15
     &            *zoopn15flag
#  endif
        remi_B = remi_B*dflag_B
        expo_B = expo_B*dflag_B
# endif
# if defined O_npzd_caco3
c Andreas commented out (don't think we need these)
c        calpro = calpro*caco3flag*coccflag*zoopflag
        dissl = dissl*caco3flag ! AS add caco3c13flag
        expocaco3 = expocaco3*caco3flag
c        calatt = calatt*caco3flag*coccflag*zoopflag
        graz_C = graz_C*coccflag*zoopflag
#  if defined O_npzd_nitrogen_15
     &        *coccn15flag*zoopn15flag 
#  endif        
        morp_C = morp_C*coccflag
#  if defined O_npzd_nitrogen_15
     &        *coccn15flag
#  endif
        morpt_C = morpt_C*coccflag
#  if defined O_npzd_nitrogen_15
     &        *coccn15flag
#  endif
# endif
# if defined O_npzd_iron
        remife=remife*detrfeflag
        feorgads=feorgads*dfeflag
        expofe=expofe*detrfeflag
        fecol=fecol*dfeflag
# endif

!     digestion of grazed material
        dig_P = gamma1*graz
        dig_Z = gamma1*graz_Z
        dig_Det = gamma1*graz_Det
# if defined O_npzd_caco3
        dig_C = gamma1*graz_C
# endif
# if defined O_kk_ballast
        dig_Det_B = gamma1*graz_Det_B
# endif        
        dig = dig_P + dig_Z + dig_Det
# if defined O_npzd_caco3
     &      + dig_C
# endif
# if defined O_kk_ballast
     &      + dig_Det_B
# endif        
!     excretion based on growth efficiency
        excr_P = gamma1*(1-geZ)*graz
        excr_Z = gamma1*(1-geZ)*graz_Z
        excr_Det = gamma1*(1-geZ)*graz_Det
# if defined O_npzd_caco3
        excr_C = gamma1*(1-geZ)*graz_C        
# endif
# if defined O_kk_ballast
       excr_Det_B = gamma1*(1-geZ)*graz_Det_B        
# endif
        excr = excr_P + excr_Z + excr_Det
# if defined O_npzd_caco3
     &       + excr_C
# endif
# if defined O_kk_ballast
     &       + excr_Det_B
# endif  
!     sloppy feeding
        sf_P = (1-gamma1)*graz
        sf_Z = (1-gamma1)*graz_Z
        sf_Det = (1-gamma1)*graz_Det
# if defined O_npzd_caco3
        sf_C = (1-gamma1)*graz_C        
# endif
# if defined O_kk_ballast
        sf_Det_B = (1-gamma1)*graz_Det_B        
# endif
        sf = sf_P + sf_Z + sf_Det   
# if defined O_npzd_caco3
     &       + sf_C
# endif
# if defined O_kk_ballast
     &       + sf_Det_B
# endif  
# if defined O_npzd_nitrogen
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
# if defined O_npzd_nitrogen_15
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

#  if defined O_npzd_caco3
        rtcoccn15 = biococcn15/biococc
        rtcoccn15 = min(rtcoccn15, 2.*rn15std/(1+rn15std))
        rtcoccn15 = max(rtcoccn15, rn15std/(1+rn15std)/2.)
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

#  if defined O_npzd_caco3
        rtcoccc13 = biococcc13/(biococc*redctn)
        rtcoccc13 = min(rtcoccc13, 2.*rc13std/(1+rc13std))
        rtcoccc13 = max(rtcoccc13, 0.5*rc13std/(1+rc13std))

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

#  if defined O_npzd_nitrogen
        rtdoc13 = biodoc13/(biodon*redctn)
        rtdoc13 = min(rtdoc13, 2*rc13std/(1+rc13std))
        rtdoc13 = max(rtdoc13, 0.5*rc13std/(1+rc13std))

        rtdiazc13 = biodiazc13/(biodiaz*redctn)
        rtdiazc13 = min(rtdiazc13, 2.*rc13std/(1+rc13std))
        rtdiazc13 = max(rtdiazc13, 0.5*rc13std/(1+rc13std))

#  endif        
# endif
# if defined O_npzd_caco3
!  net formation of attached (living) tests
!     total primary production by coccs
        calatt = (npp_C - morpt_C - morp_C - graz_C
! total growth of zooplankton
     &         + dig - morz - graz_Z - excr     
#  if defined O_kk_ballast
!     &         + graz_Det_B
#  endif
! convert to carbon units mmol C
     &          )*capr*redctn*1.e3

!   formation of detached (dead) tests, or PIC
        calpro = (sf_C + sf_Z + morp_C + morz)*capr*redctn*1.e3 !stay in mmol
# else
        calpro = (morp+morz+(graz+graz_Z)*(1.-gamma1))*capr*redctn*1.e3
# endif        
# if defined O_npzd_nitrogen
!       nutrients equation
        biopo4 = biopo4 + dtbio*(redptn*(excr + remi 
     &           + (1.-dfrt)*morpt - (npp-dopupt)
#  if defined O_npzd_caco3
     &           + (1.-dfrt)*morpt_C - (npp_C-dopupt_C)
#  endif
     &            )
     &         + diazptn*(morpt_D - (npp_D-dopupt_D)) + recy_dop
     &          )

!       DOP equation
        biodop = biodop + dtbio*(redptn*(
     &             dfr*morp + dfrt*morpt - dopupt
#  if defined O_npzd_caco3
     &           + dfr*morp_C + dfrt*morpt_C - dopupt_C
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
#   if defined O_npzd_caco3
     &          + bapr*(1.-dfr)*morp_C
#   endif
     &           )
        biodetr = biodetr + dtbio*((1.-bapr)*((1.-dfr)*morp + sf + morz)
     &          - remi - graz_Det - expo + impo
     &          + (1.-bapr)*morp_D*(redntp/diazntp)
#   if defined O_npzd_caco3
     &          + (1.-bapr)*(1.-dfr)*morp_C
#   endif
     &          + remi_B)
#  else        
        biodetr = biodetr + dtbio*((1.-dfr)*morp + sf + morz - remi 
     &          - graz_Det - expo + impo + morp_D*(redntp/diazntp)
#   if defined O_npzd_caco3
     &          + (1.-dfr)*morp_C
#   endif
     &           )
#  endif

#  if defined O_carbon
!       DIC equation
        biodic = biodic + dtbio*redctn*(excr + remi 
     &         + (1.-dfrt)*morpt - npp 
#   if defined O_npzd_caco3
     &         + (1.-dfrt)*morpt_C - npp_C
#   endif
     &         + morpt_D - npp_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )        
#  endif
!       nitrate (NO3) equation
        biono3 = biono3 + dtbio*(excr + remi   
     &         + (1.-dfrt)*morpt - npp 
#  if defined O_npzd_caco3
     &         + (1.-dfrt)*morpt_C - npp_C
#  endif
     &         + morpt_D - no3upt_D + recy_don + nr_excr_D
     &         + morp_D*(1.-(redntp/diazntp))
     &           )
!       DON equation
        biodon = biodon + dtbio*(dfr*morp + dfrt*morpt - recy_don
#  if defined O_npzd_caco3
     &         + dfr*morp_C + dfrt*morpt_C
#  endif
     &       )
!       diazotroph equation
        biodiaz = biodiaz + dtbio*(npp_D - morp_D - morpt_D - graz_D)
# else
!       nutrients equation
        biopo4 = biopo4 + dtbio*redptn*(remi + excr - npp + morpt
#  if defined O_npzd_caco3
     &         - npp_C + morpt_C
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
     &           )
        biodetr = biodetr + dtbio*((1.-bapr)*(morp + sf + morz)
     &          - remi - graz_Det - expo + impo
     &          + (1.-bapr)*morp_C
     &          + remi_B) 
#  else        
        biodetr = biodetr + dtbio*(morp + sf + morz - remi 
     &          - graz_Det - expo + impo
#   if defined O_npzd_caco3
     &          + morp_C
#   endif
     &           )
#  endif        
#  if defined O_carbon
!     DIC equation
        biodic = biodic + dtbio*redctn*(morpt + excr + remi - npp        
#   if defined O_npzd_caco3
     &          + morpt_C - npp_C
#   endif
     &           )  
#  endif
# endif
# if defined O_npzd_caco3
!       coccolithophores equation
        biococc = biococc + dtbio*(npp_C - morp_C - graz_C
     &       - morpt_C)
!     calcite equation
c        biocaco3 = biocaco3 - dissl + dtbio*(calpro-expocaco3+impocaco3)
c     Andreas put dissl in bracket
        biocaco3 = biocaco3 + dtbio*(calpro-dissl-expocaco3+impocaco3)
# endif        
# if defined O_npzd_iron
!       dissolved iron equation
#  if defined O_npzd_nitrogen
        biodfe = biodfe + dtbio*(rfeton*(excr + (1.-dfrt)*morpt 
     &            - npp + morpt_D - npp_D + recy_don + nr_excr_D 
     &            + morp_D*(1-(redntp/diazntp))) - feorgads 
     &            + remife - fecol
#   if defined O_npzd_caco3
     &            + rfeton*((1.-dfrt)*morpt_C - npp_C)
#   endif
     &            )   
!     particulate iron equation
        biodetrfe = biodetrfe + dtbio*(rfeton*(sf + (1.-dfr)*morp
     &            + morp_D*(redntp/diazntp) + morz - graz_Det) 
     &            + feorgads + fecol - remife - expofe + impofe        
#   if defined O_npzd_caco3
     &            + rfeton*(1.-dfr)*morp_C         
#   endif        
     &       )
#  else
        biodfe = biodfe + dtbio*(rfeton*(excr + morpt 
     &            - npp) - feorgads + remife - fecol        
#   if defined O_npzd_caco3
     &            + rfeton*(morpt_C - npp_C)
#   endif
     &       )
!     particulate iron equation
        biodetrfe = biodetrfe + dtbio*(rfeton*(sf + morp
     &             + morz - graz_Det) 
     &            + feorgads + fecol - remife - expofe + impofe
#   if defined O_npzd_caco3
     &            + rfeton*morp_C         
#   endif
#   if defined O_kk_ballast
     &            - rfeton*graz_Det_B
#   endif        
     &       )
#  endif
# endif
# if defined O_npzd_nitrogen_15
!       isotope equations
        biodin15 = biodin15 + dtbio*(rtphytn15*(1.-dfrt)*morpt 
#  if defined O_npzd_caco3
     &       + rtcoccn15*(1.-dfrt)*morpt_C - fcassim*npp_C
#  endif        
     &       + fcexcr*excr + rtdiazn15*morpt_D + rtdiazn15*nr_excr_D 
     &       + rtdiazn15*morp_D*(1-(redntp/diazntp)) + rtdetrn15*remi 
     &       + fcrecy*recy_don - fcassim*npp - fcassim*no3upt_D)

        biodon15 = biodon15 + dtbio*(dfr*rtphytn15*morp 
#  if defined O_npzd_caco3
     &       + dfr*rtcoccn15*morp_C + dfrt*rtcoccn15*morpt_C
#  endif
     &       + dfrt*rtphytn15*morpt - fcrecy*recy_don)

        biophytn15 = biophytn15 + dtbio*(fcassim*npp - rtphytn15*morp 
     &       - rtphytn15*graz - rtphytn15*morpt)
#  if defined O_npzd_caco3
        biococcn15 = biococcn15 + dtbio*(fcassim*npp_C-rtcoccn15*morp_C 
     &       - rtcoccn15*graz_C - rtcoccn15*morpt_C)
#  endif        
        biozoopn15 = biozoopn15 + dtbio*(rtphytn15*dig_P 
#  if defined O_npzd_caco3
     &       + rtcoccn15*dig_C
#  endif        
     &       + rtzoopn15*dig_Z + rtdetrn15*dig_Det + rtdiazn15*dig_D 
     &       - rtzoopn15*morz - rtzoopn15*graz_Z - fcexcr*excr)

        biodetrn15 = biodetrn15 + dtbio*(rtphytn15*(1.-dfr)*morp 
#  if defined O_npzd_caco3
     &       + rtcoccn15*(1.-dfr)*morp_C + rtcoccn15*sf_C
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
#  if defined O_npzd_nitrogen
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*(1.-dfrt)*morpt 
     &       + rtzoopc13*excr + rtdiazc13*morpt_D + rtdiazc13*nr_excr_D 
     &       + rtdiazc13*morp_D*(1-(redntp/diazntp)) + rtdetrc13*remi 
#   if defined O_npzd_caco3
     &       + rtcoccc13*(1.-dfrt)*morpt_C - fcnpp*npp_C
#   endif        
     &       + rtdoc13*recy_don - fcnpp*npp - fcnpp*npp_D)

        biodoc13 = biodoc13 + dtbio*redctn*(dfr*rtphytc13*morp 
#   if defined O_npzd_caco3
     &       + rtcoccc13*(dfr*morp_C + dfrt*morpt_C)
#   endif       
     &       + rtphytc13*dfrt*morpt - rtdoc13*recy_don)        

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp 
     &       - rtphytc13*morp - rtphytc13*graz - rtphytc13*morpt)
        
        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P 
#   if defined O_npzd_caco3
     &      + rtcoccc13*dig_C
#   endif        
     &      + rtzoopc13*dig_Z + rtdetrc13*dig_Det + rtdiazc13*dig_D 
     &       - rtzoopc13*morz - rtzoopc13*graz_Z - rtzoopc13*excr)
      
        biodetrc13 = biodetrc13 + dtbio*redctn*(rtphytc13*(1.-dfr)*morp 
#   if defined O_npzd_caco3
     &       + rtcoccc13*(1.-dfr)*morp_C + rtcoccc13*sf_C
#   endif        
     &       + rtphytc13*sf_P + rtzoopc13*sf_Z + rtdetrc13*sf_Det 
     &       + rtzoopc13*morz - rtdetrc13*remi 
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo 
     &       + rtdiazc13*morp_D*(redntp/diazntp))
       
        biodiazc13 = biodiazc13 + dtbio*redctn*(fcnpp*npp_D 
     &       - rtdiazc13*(morp_D + graz_D + morpt_D))
#  else ! no nitrogen
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*morpt 
     &       + rtzoopc13*excr + rtdetrc13*remi - fcnpp*npp
#   if defined O_npzd_caco3
     &       + rtcoccc13*morpt_C - fcnpp*npp_C
#   endif        
     &        )

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp 
     &       - rtphytc13*(morp + graz + morpt))
        
        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P 
#   if defined O_npzd_caco3
     &      + rtcoccc13*dig_C
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
     &       + rtzoopc13*morz - rtdetrc13*remi 
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo 
#   if defined O_npzd_caco3
     &       + rtcoccc13*morp_C + rtcoccc13*sf_C
#   endif        
     &       )

#  endif
#  if defined O_npzd_caco3
        biococcc13 = biococcc13 + dtbio*redctn*(fcnpp*npp_C 
     &       - rtcoccc13*(morp_C + graz_C + morpt_C))
        biocaco3c13 = biocaco3c13 + dtbio*(rtdic13*calpro
     &       - rtcaco3c13*dissl - rtcaco3c13*expocaco3 
     &                + rcaco3c13impo)
#  endif
# endif
        expoout = expoout + expo
# if defined O_npzd_nitrogen_15
        rn15expoout = rn15expoout + rtdetrn15
# endif
# if defined O_carbon_13
        rc13expoout = rc13expoout + rtdetrc13*expo
#  if defined O_npzd_caco3
        rcaco3c13expoout = rcaco3c13expoout + rtcaco3c13*expocaco3
#  endif
# endif
        calproout = calproout + calpro
# if defined O_npzd_caco3
        calattout = calattout + calatt
        disslout = disslout + dissl
        expocaco3out = expocaco3out + expocaco3
# endif
# if defined O_npzd_nitrogen
        nfixout = nfixout + npp_D - no3upt_D
# endif
# if defined O_npzd_iron
        expofeout = expofeout + expofe
        remifeout = remifeout + remife
# endif
# if defined O_save_npzd
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
        graz_Det_out = graz_Det_out + graz_Det
        graz_Zout = graz_Zout + graz_Z
        nppout = nppout + npp
        morptout = morptout + morpt
        remiout = remiout + remi
        excrout = excrout + excr
#  if defined O_npzd_nitrogen
        npp_dopout = npp_dopout + npp*dopupt_flag        
#   if defined O_npzd_caco3
        npp_C_dopout = npp_C_dopout + npp_C*dopupt_C_flag  
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
#  if defined O_npzd_caco3
        npp_Cout = npp_Cout + npp_C
        graz_Cout = graz_Cout + graz_C
        morp_Cout = morp_Cout + morp_C
        morpt_Cout = morpt_Cout + morpt_C
#  endif        
#  if defined O_npzd_iron
#   if defined O_npzd_iron_diagnostics
        feorgadsout = feorgadsout + feorgads
        deffeout = deffeout + deffe
        feprimeout = feprimeout + feprime
        fecolout = fecolout + fecol
#    if defined O_npzd_caco3
        deffeout_C = deffeout_C + deffe_C
#    endif        
#   endif
#  endif
#  if defined O_npzd_extra_diagnostics
        avej_out = avej_out + avej
        gmax_out = gmax_out + gmax
        po4P_out = po4P_out + po4P
#   if defined O_npzd_nitrogen
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
# if defined O_npzd_nitrogen
        if (no3flag .eq. 1) no3flag = 0.5 + sign(0.5,biono3 - trcmin)
        if (dopflag .eq. 1) dopflag = 0.5 + sign(0.5,biodop - trcmin)        
        if (donflag .eq. 1) donflag = 0.5 + sign(0.5,biodon - trcmin)
        if (diazflag .eq. 1) diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_npzd_nitrogen_15
        if (din15flag .eq. 1) 
     &     din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
        if (don15flag .eq. 1)
     &     don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
        if (phytn15flag .eq. 1) 
     &     phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
#   if defined O_npzd_caco3
        if (coccn15flag .eq. 1) 
     &     coccn15flag = 0.5 + sign(0.5, biococcn15 - trcmin)
#   endif        
        if (zoopn15flag .eq. 1)
     &     zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
        if (detrn15flag .eq. 1)
     &     detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
        if (diazn15flag .eq. 1)
     &     diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_npzd_caco3
        if (caco3flag .eq. 1) 
     &     caco3flag = 0.5 + sign(0.5,biocaco3 - trcmin)
        if (coccflag .eq. 1) 
     &     coccflag = 0.5 + sign(0.5,biococc - trcmin)
#  if defined O_kk_ballast
        if (dflag_B .eq. 1) 
     &     dflag_B = 0.5 + sign(0.5,biod_B - trcmin)
#  endif
# endif        
# if defined O_npzd_iron
        if (dfeflag .eq. 1) dfeflag = 0.5 + sign(0.5,biodfe - trcmin)
        if (detrfeflag .eq. 1) detrfeflag = 0.5 
     &                                 + sign(0.5,biodetrfe - trcmin)
# endif
# if defined O_carbon_13
        if (dic13flag .eq. 1)
     &       dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
        if (phytc13flag .eq. 1)
     &       phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
#  if defined O_npzd_caco3
        if (coccc13flag .eq. 1)
     &       coccc13flag = 0.5 + sign(0.5,biococcc13 - trcmin)
        if (caco3c13flag .eq. 1)
     &       caco3c13flag = 0.5 + sign(0.5,biocaco3c13 - trcmin)
#  endif        
        if (zoopc13flag .eq. 1)
     &       zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
        if (detrc13flag .eq. 1)
     &       detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
#  if defined O_npzd_nitrogen
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
# if defined O_npzd_nitrogen
      bioout(imobidop) = biodop - bioin(imobidop)
      bioout(imobino3) = biono3 - bioin(imobino3)
      bioout(imobidon) = biodon - bioin(imobidon)
      bioout(imobidiaz) = biodiaz - bioin(imobidiaz)
#  if defined O_npzd_nitrogen_15
      bioout(imobidin15) = biodin15 - bioin(imobidin15)
      bioout(imobidon15) = biodon15 - bioin(imobidon15)
      bioout(imobiphytn15) = biophytn15 - bioin(imobiphytn15)
      bioout(imobizoopn15) = biozoopn15 - bioin(imobizoopn15)
      bioout(imobidetrn15) = biodetrn15 - bioin(imobidetrn15)
      bioout(imobidiazn15) = biodiazn15 - bioin(imobidiazn15)
#   if defined O_npzd_caco3
      bioout(imobicoccn15) = biococcn15 - bioin(imobicoccn15)      
#   endif
#  endif
# endif
# if defined O_npzd_caco3
      bioout(imobicocc) = biococc - bioin(imobicocc)
      bioout(imobicaco3) = biocaco3 - bioin(imobicaco3)
#  if defined O_kk_ballast
      bioout(imobidetr_B) = biod_B - bioin(imobidetr_B)
#  endif
# endif
# if defined O_npzd_iron
      bioout(imobidfe) = biodfe - bioin(imobidfe)
      bioout(imobidetrfe) = biodetrfe - bioin(imobidetrfe)
# endif
# if defined O_carbon_13
      bioout(imobidic13) = biodic13 - bioin(imobidic13)
      bioout(imobiphytc13) = biophytc13 - bioin(imobiphytc13)
      bioout(imobizoopc13) = biozoopc13 - bioin(imobizoopc13)
      bioout(imobidetrc13) = biodetrc13 - bioin(imobidetrc13)
#  if defined O_npzd_nitrogen
      bioout(imobidoc13) = biodoc13 - bioin(imobidoc13)
      bioout(imobidiazc13) = biodiazc13 - bioin(imobidiazc13)
#  endif      
#  if defined O_npzd_caco3
      bioout(imobicoccc13) = biococcc13 - bioin(imobicoccc13)
      bioout(imobicaco3c13) = biocaco3c13 - bioin(imobicaco3c13)
#  endif
# endif
#endif
      return
      end
