!====================== include file "mobi.h" =========================

!   variables for Model of Ocean Biogeochemistry and Isotopes (MOBI)

!   ntnpzd   = number of npzd tracers
!   nbio     = number of npzd timesteps per ocean timestep
!   trcmin   = minimum tracer for flux calculations
!   tap      = 2*alpha*par with
!   alpha    = initial slope P-I curve [(W/m^2)^(-1)/day] and
!   par      = fraction of photosythetically active radiation
!   kw       = light attenuation due to water [1/m]
!   kc       = light attenuation by phytoplankton [1/(m*mmol m-3)]
!   ki       = light attenuation through sea ice & snow
!   abio     = maximum growth rate parameter [1/day]
!   bbio     = b
!   cbio     = [1/deg_C]
!   k1n      = half saturation constant for N uptake [mmol m-3]
!   nup      = specific mortality rate (Phytoplankton) [day-1]
!   gamma1   = assimilation efficiency (zpk)
!   gbio     = maximum grazing rate at 0 deg C [day-1]
!   nuz      = quadratic mortality (zpk)
!   nud0     = remineralization rate [day-1]
!   nudop0   = DON remineralization rate [day-1]
!   nudon0   = DOP remineralization rate [day-1]
!   wd       = sinking speed of detritus [m day-1]
!   ztt      = depth to top of grid cell [cm]
!   dtnpzd   = time step of biology
!   caprmax  = maximum carbonate to carbon production ratio
!   capr     = carbonate to carbon production ratio
!   kcapr    = half saturation for capr
!   dcaco3   = remineralisation depth of calcite [cm]
!   rcak     = array used in calculating calcite remineralization
!   rcab     = array used in calculating bottom calcite remineralization
!   nupt0    = specific mortality rate (Phytoplankton) [1/day]
!   wd0      = sinking speed of detritus at surface [m/day]
!   mw       = sinking speed increase with depth [1/day]
!   mw_c     = calcite sinking speed increase with depth [1/day]
!   mwz      = sinking speed increase depth cut-off (cm)
!   k1p_P    = half saturation constant for P uptake phytoplankton
!   jdiar    = factor reducing the growth rate of diazotrophs
!   dbct_D   = subtracted from bct for diazotrophs
!   redctn   = C/N Redfield ratio (includes mol to mmol conversion)
!   redctp   = C/P Redfield ratio (includes mol to mmol conversion)
!   redptn   = P/N Redfield ratio
!   redntp   = N/P Redfield ratio
!   redotn   = O/N Redfield ratio (includes mol to mmol conversion)
!   redotp   = O/P Redfield ratio (includes mol to mmol conversion)
!   rnbio    = reciprical of nbio
!   rdtts    = reciprical of dtts [s-1]
!   dtbio    = npzd time step [s]
!   rnpp     = rate of net primary production [nmol cm-3 s-1]
!   rgraz    = rate of grazing [nmol cm-3 s-1]
!   rmorp    = rate of mortality of phytoplankton [nmol cm-3 s-1]
!   rmorz    = rate of mortality of zooplankton [nmol cm-3 s-1]
!   rremi    = rate of remineralization [nmol cm-3 s-1]
!   rexcr    = rate of excretion [nmol cm-3 s-1]
!   rexpo    = rate of export through the bottom [nmol cm-3 s-1]
!   rnpp_D   = npp for diazotraphs [nmol cm-3 s-1]
!   rgraz_D  = rgraz for diazotraphs [nmol cm-3 s-1]
!   rmorpt_D = rmorp for diazotraphs [nmol cm-3 s-1]
!   rnfix    = rate of nitrogen fixation [nmol cm-3 s-1]
!   rdeni    = rate of denitrification [nmol cm-3 s-1]
!   kzoo     = half saturation constant for Z grazing
!   zprefP   = Z preference for grazing on P
!   zprefDiaz   = Z preference for grazing on Diaz
!   zprefZ   = Z preference for grazing on other Z
!   zprefDet = Z preference for grazing on Detritus
!   rgraz_Det = rate of grazing on Detritus [nmol cm-3 s-1]
!   rgraz_Z   = rate of grazing on other Zooplankton [nmol cm-3 s-1]
!   geZ      = growth efficiency of zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! New diagnostic output
!   ravej      = light-dependant growth rate of P
!   ravej_D    = light-dependant growth rate of Diaz
!   rgmax      = temp-dependant growth rate of zoo
!   rno3P      = nitrate-dependant growth rate of P
!   rpo4P       = phosphate-dependant growth rate of P
!   rpo4_D     = phosphate-dependant growth rate of D
!
!   fe_dissolved = dissolved iron concentration
!   kfe = Fe limitation half saturation parameter
!   kfe_D = Fe limitation half sat. param. for diaz.
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Dynamic iron cycle
!   kfeleq = Fe-ligand stability constant [m^3/(mmol ligand)]
!   lig = Ligand concentration  [mmol/m^3]
!   thetamaxhi = Maximum Chl:C ratio, abundant iron [gChl/(gC)]
!   thetamaxlo = Maximum Chl:C ratio, extreme iron limitation [gChl/(gC)]
!   alphamax = Maximum initial slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   alphamin = Minimum intital slope in PI-curve [day^-1 (W/m^2)^-1 (gChl/(gC))^-1]
!   mc = Molar mass of carbon [g/mol]
!   fetopsed = Fe:P ratio for sedimentary iron source [molFe/molP]
!   o2min = Minimum O2 concentration for aerobic respiration [mmolO_2/m^3]
!   kfeorg = Organic-matter dependent scavenging rate [(m^3/(gC s))^0.58]
!   kfecol = Colloidal production and precipitation rate [s^-1]
			 
      integer ntnpzd, nbio
      parameter (ntnpzd = 4 ! po4, phyt, zoop, detr
#if defined O_carbon
     &                  + 1 ! dic
# if defined O_carbon_13
     &                  + 4 ! dic13, phytc13, zoopc13, detrc13
#  if defined O_mobi_nitrogen
     &                  + 2 ! doc13, diazc13
#  endif		 
#  if defined O_mobi_caco3
     &                  + 1 ! caco3c13
#  endif		 
#  if defined O_mobi_silicon
     &                  + 1 ! diatc13
#  endif		 
# endif
#endif		 
#if defined O_mobi_nitrogen
     &                  + 4 ! no3, diaz, don, dop
# if defined O_mobi_nitrogen_15
     &                  + 6 ! din15, don15, phytn15, zoopn15, detrn15, diazn15
#  if defined O_mobi_silicon
     &                  + 1 ! diatn15
#  endif		 
# endif
#endif
#if defined O_mobi_caco3
     &                   +1 ! caco3
#endif
#if defined O_kk_ballast
     &                   +1 ! detr_B
#endif
#if defined O_mobi_silicon
     &                   +3 ! diat, sil, opl
#endif
#if defined O_mobi_iron
     &                  + 2 ! dfe, detrfe
#endif
     &                     )
      common /npzd_i/ nbio

      integer imobipo4, imobiphyt, imobizoop, imobidetr
      common /npzd_i/ imobipo4, imobiphyt, imobizoop, imobidetr
#if defined O_carbon
      integer imobidic
      common /npzd_i/ imobidic
# if defined O_carbon_13
      integer imobidic13, imobidoc13, imobiphytc13
      integer imobizoopc13, imobidetrc13, imobidiazc13
      common /npzd_i/ imobidic13, imobidoc13, imobiphytc13
      common /npzd_i/ imobizoopc13, imobidetrc13, imobidiazc13
#  if defined O_mobi_silicon
      integer imobidiatc13
      common /npzd_i/ imobidiatc13
#  endif
#  if defined O_mobi_caco3
      integer imobicaco3c13
      common /npzd_i/ imobicaco3c13
#  endif
# endif
#endif
#if defined O_mobi_nitrogen
      integer imobidop, imobino3, imobidon, imobidiaz
      common /npzd_i/ imobidop, imobino3, imobidon, imobidiaz
# if defined O_mobi_nitrogen_15
      integer imobidin15, imobidon15, imobiphytn15
      integer imobizoopn15, imobidetrn15, imobidiazn15
      common /npzd_i/ imobidin15, imobidon15, imobiphytn15
      common /npzd_i/ imobizoopn15, imobidetrn15, imobidiazn15
#  if defined O_mobi_silicon
      integer imobidiatn15
      common /npzd_i/ imobidiatn15
#  endif
# endif
#endif
#if defined O_mobi_caco3
      integer imobicaco3
      common /npzd_i/ imobicaco3
#endif
#if defined O_kk_ballast
      integer imobidetr_B
      common /npzd_i/ imobidetr_B
#endif
# if defined O_mobi_silicon
      integer imobidiat, imobisil, imobiopl
      common /npzd_i/ imobidiat, imobisil, imobiopl
# endif
#if defined O_mobi_iron
      integer imobidfe, imobidetrfe
      common /npzd_i/ imobidfe, imobidetrfe
#endif

      real trcmin
      parameter (trcmin=5e-12)

#if defined O_mobi_nitrogen_15
      real rn15std
      parameter (rn15std=0.0036765)
#endif
#if defined O_carbon_13
      real rc13std
      parameter (rc13std=0.0112372)
#endif
#if defined O_carbon_14
!     rc14std   = standard c14/c12 ratio
      real rc14std
      parameter (rc14std=1.176e-12)
#endif

#if defined O_mobi
      real            tap, kw, kc, ki, abio_P, bbio, cbio, k1n, nup
      common /npzd_r/ tap, kw, kc, ki, abio_P, bbio, cbio, k1n, nup
      real            gamma1, gbio, epsbio, nuz, nud0, LFe, dfr
      common /npzd_r/ gamma1, gbio, epsbio, nuz, nud0, LFe, dfr
      real            wd,     ztt,     dtnpzd, capr, caprmax, kcapr
      common /npzd_r/ wd(km), ztt(km), dtnpzd, capr, caprmax, kcapr
      real            dcaco3, rcak,     rcab,     nupt0, wd0, k1p_P
      common /npzd_r/ dcaco3, rcak(km), rcab(km), nupt0, wd0, k1p_P
      real            redptc, redctn, redctp, redptn, redntp, redotn
      common /npzd_r/ redptc, redctn, redctp, redptn, redntp, redotn
      real            redotp, redotc, redntc, rnbio, rdtts, dtbio, geZ
      common /npzd_r/ redotp, redotc, redntc, rnbio, rdtts, dtbio, geZ
      real            kzoo, zprefP, zprefDiaz, zprefZ, zprefDet
      common /npzd_r/ kzoo, zprefP, zprefDiaz, zprefZ, zprefDet
      real            zprefC, zprefDiat
      common /npzd_r/ zprefC, zprefDiat
      real            kfe, kfe_D, sgbdfac, nupt0_D, diazntp, diazptn
      common /npzd_r/ kfe, kfe_D, sgbdfac, nupt0_D, diazntp, diazptn
      real            jdiar, dbct_D, nup_D, dfrt, nudop0, nudon0
      common /npzd_r/ jdiar, dbct_D, nup_D, dfrt, nudop0, nudon0
      real            eps_bdeni0, eps_recy, hdop, mw, mwz, mw_c
      common /npzd_r/ eps_bdeni0, eps_recy, hdop, mw, mwz, mw_c
      real            eps_assim, eps_excr, eps_nfix, eps_wcdeni
      common /npzd_r/ eps_assim, eps_excr, eps_nfix, eps_wcdeni
      real            abio_C, k1n_C, k1p_C, nuc, nuct0, kfe_C, tap_C
      common /npzd_r/ abio_C, k1n_C, k1p_C, nuc, nuct0, kfe_C, tap_C
      real            kcal, wc0, dissk0, rdissl,     wc
      common /npzd_r/ kcal, wc0, dissk0, rdissl(km), wc(km)
      real            rcalatt,       rexpocaco3,     rimpocaco3
      common /npzd_r/ rcalatt(kpzd), rexpocaco3(km), rimpocaco3(km)
      real            kc_c, rcalpro, romega_c
      common /npzd_r/ kc_c, rcalpro(kpzd), romega_c(kpzd)

# if defined O_mobi_iron
      real            kfeleq, alphamax, alphamin, lig, kfeorg, rfeton
      common /npzd_r/ kfeleq, alphamax, alphamin, lig, kfeorg, rfeton
      real            thetamaxhi, thetamaxlo, mc, fetopsed, o2min
      common /npzd_r/ thetamaxhi, thetamaxlo, mc, fetopsed, o2min
      real            kfecol, kfemax, kfemin, pmax
      common /npzd_r/ kfecol, kfemax, kfemin, pmax
      real            knmax, knmin
      common /npzd_r/ knmax, knmin
#  if defined O_mobi_caco3
      real            kfemax_C, kfemin_C, pmax_C
      common /npzd_r/ kfemin_C, kfemax_C, pmax_C
      real            knmax_C, knmin_C
      common /npzd_r/ knmax_C, knmin_C
#  endif
#  if defined O_mobi_silicon
      real            kfemax_Diat, kfemin_Diat, pmax_Diat
      common /npzd_r/ kfemin_Diat, kfemax_Diat, pmax_Diat
      real            knmax_Diat, knmin_Diat
      common /npzd_r/ knmax_Diat, knmin_Diat
#  endif
      real fe_hydr
      common /fe_hydr/ fe_hydr(imt,jmt,km)
# endif
# if defined O_kk_ballast
      real bapr
      common /npzd_r/ bapr
# endif
# if defined O_mobi_silicon
      real abiodiat, k1n_Diat, k1p_Diat, nu_diat, nudt0
      real alpha_Diat, k1si, si_msk
      real kfe_Diat
      common /npzd_r/ abiodiat, k1n_Diat, k1p_Diat, nu_diat, nudt0
      common /npzd_r/ alpha_Diat, k1si
      common /npzd_r/ kfe_Diat, si_msk(imt,jmt,km)      
      real sipr0, ropk, dopal, si_sol, si_h_sol, opl_disk0
      real globalsilwflx, sildustflux, sildustfluxfac, wo0, rivsil
      real prop, bsi, wo
      real si_hydr
      common /npzd_r/ sipr0, ropk(km), dopal, si_sol, si_h_sol
      common /npzd_r/ globalsilwflx, sildustflux, sildustfluxfac
      common /npzd_r/ prop, bsi(km), wo(km), wo0, rivsil, opl_disk0
      common /npzd_r/ si_hydr(imt,jmt,km)
      real            rexpoopl,     rdisopl,     rproopl
      common /npzd_r/ rexpoopl(km), rdisopl(km), rproopl(km)
# endif
# if defined O_save_mobi_fluxes
      real rnpp, rgraz, rmorp, rmorpt, rmorz, rexcr, rremi, rexpo
      common /npzd_r/ rnpp(kpzd), rgraz(kpzd), rmorp(kpzd), rmorpt(kpzd)
      common /npzd_r/ rmorz(kpzd), rexcr(kpzd), rremi(km), rexpo(km)
      real rgraz_Det, rgraz_Z, rsedrr, rprca, rnpp_dop
      common /npzd_r/ rgraz_Det(kpzd), rgraz_Z(kpzd), rsedrr, rprca
      common /npzd_r/ rnpp_dop(kpzd)
#  if defined O_mobi_caco3
      real            rnpp_C,       rgraz_C,       rmorp_C
      common /npzd_r/ rnpp_C(kpzd), rgraz_C(kpzd), rmorp_C(kpzd)
      real            rmorpt_C,       rnpp_C_dop
      common /npzd_r/ rmorpt_C(kpzd), rnpp_C_dop(kpzd)
#  endif
#  if defined O_kk_ballast
      real            rgraz_Det_B,       rremi_B,     rexpo_B
      common /npzd_r/ rgraz_Det_B(kpzd), rremi_B(km), rexpo_B(km)
#  endif
#  if defined O_mobi_silicon
      real            rnpp_Diat,       rgraz_Diat
      common /npzd_r/ rnpp_Diat(kpzd), rgraz_Diat(kpzd)
      real            rmorp_Diat,       rmorpt_Diat
      common /npzd_r/ rmorp_Diat(kpzd), rmorpt_Diat(kpzd)
      real            rnpp_Diat_dop
      common /npzd_r/ rnpp_Diat_dop(kpzd)
#  endif
#  if defined O_save_mobi_diagnostics
      real ravej, ravej_D, rgmax, rno3P, rpo4P, rpo4_D
      common /npzd_r/ ravej(kpzd), ravej_D(kpzd), rgmax(kpzd)
      common /npzd_r/ rno3P(kpzd), rpo4P(kpzd), rpo4_D(kpzd)
#  endif
#  if defined O_mobi_iron
      real            rremife,     rexpofe
      common /npzd_r/ rremife(km), rexpofe(km)
#   if defined O_save_mobi_diagnostics
      real            rfeorgads,     rdeffe,     rfeprime
      common /npzd_r/ rfeorgads(km), rdeffe(km), rfeprime(km)
      real            rfesed,     rbfe,     rfecol
      common /npzd_r/ rfesed(km), rbfe(km), rfecol(km)
#    if defined O_mobi_caco3
      real            rdeffe_C
      common /npzd_r/ rdeffe_C(km)
#    endif     
#    if defined O_mobi_silicon
      real            rdeffe_Diat
      common /npzd_r/ rdeffe_Diat(km)
#    endif     
#   endif
#  endif
#  if defined O_mobi_nitrogen
      real rnpp_D, rgraz_D, rmorp_D, rmorpt_D, rnfix, rwcdeni, rbdeni
      real rnpp_D_dop
      common /npzd_r/ rnpp_D(kpzd), rgraz_D(kpzd), rmorp_D(kpzd)	 
      common /npzd_r/ rmorpt_D(kpzd), rnfix(kpzd), rwcdeni(km)
      common /npzd_r/ rbdeni(km), rnpp_D_dop(kpzd)

#  endif
# endif
#endif
