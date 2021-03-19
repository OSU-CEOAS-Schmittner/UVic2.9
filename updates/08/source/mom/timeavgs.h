! source file: /Users/dkeller/Desktop/UVic_ESCM/2.9/source/mom/timeavgs.h
!====================== include file "timeavgs.h" ======================

!     imtav =  # of longitudes used for the time averages grid
!     jmtav =  # of latitudes used for the time averages grid
!     kmav  =  # of levels used for the time averages grid

      integer imtav, jmtav, kmav
      parameter (imtav=imt, jmtav=jmt-2, kmav=km)
#if defined O_isopycmix && defined O_gent_mcwilliams && defined O_time_averages
      real ta_vetiso, ta_vntiso, ta_vbtiso
      common /ta_gm_r/ ta_vetiso(imt,km,jmt), ta_vntiso(imt,km,jmt)
      common /ta_gm_r/ ta_vbtiso(imt,km,jmt)
#endif
#if defined O_carbon_14
      real ta_dc14
      common /ta_dc14/ ta_dc14(imt,km,jmt)
#endif
#if defined O_save_kv
      real ta_diff_cbt
      common /ta_kv/ ta_diff_cbt(imt,km,jmt)
#endif
#if defined O_KGMvar
      real ta_kgm, ta_isop, ta_Lr, ta_LRhi
      common /ta_gm_r/ ta_kgm(imt,km,jmt,1)
      common /ta_gm_r/ ta_isop(imt,km,jmt,1)
      common /ta_gm_r/ ta_Lr(imt,jmt)
      common /ta_gm_r/ ta_LRhi(imt,km,jmt)
#endif
#if defined O_save_mobi_fluxes
      real ta_rnpp, ta_rgraz, ta_rmorp, ta_rmorpt, ta_rmorz, ta_rexcr
      real ta_rremi, ta_rexpo, ta_rexpocal, ta_rsedrr
      real ta_rgraz_Det, ta_rgraz_Z
      common /ta_npzd_r/ ta_rnpp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexcr(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremi(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpocal(imt,km,jmt)
      common /ta_npzd_r/ ta_rsedrr(imt,jmt)
      common /ta_npzd_r/ ta_rgraz_Det(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_Z(imt,kpzd,jmt)
# if defined O_save_mobi_diagnostics
      real ta_ravej, ta_ravej_D, ta_rgmax, ta_rno3P
      real ta_rpo4P, ta_rpo4_D
      common /ta_npzd_r/ ta_ravej(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_ravej_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgmax(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rno3P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4P(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rpo4_D(imt,kpzd,jmt)
# endif
# if defined O_pipe_co2
      real ta_rpo4P
      common /ta_npzd_r/ ta_kpipe(imt,jmt)
# endif
# if defined O_mobi_nitrogen
      real ta_rnpp_dop, ta_rnpp_D, ta_rnpp_D_dop
      real ta_rgraz_D, ta_rmorp_D, ta_rmorpt_D
      real ta_rnfix, ta_rwcdeni, ta_rbdeni
      common /ta_npzd_r/ ta_rnpp_dop(imt,kpzd,jmt) 
      common /ta_npzd_r/ ta_rnpp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp_D_dop(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnfix(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rwcdeni(imt,km,jmt)
      common /ta_npzd_r/ ta_rbdeni(imt,km,jmt)
#  if defined O_mobi_caco3
      real ta_rnpp_C_dop
      common /ta_npzd_r/ ta_rnpp_C_dop(imt,kpzd,jmt) 
#  endif
#  if defined O_mobi_silicon
      real ta_rnpp_Diat_dop
      common /ta_npzd_r/ ta_rnpp_Diat_dop(imt,kpzd,jmt) 
#  endif
# endif
# if defined O_kk_ballast
      real ta_rremi_B, ta_rexpo_B, ta_rgraz_Det_B
# endif
      real ta_rprocal
# if defined O_mobi_caco3
      real ta_rnpp_C, ta_rgraz_C, ta_rmorp_C, ta_rmorpt_C
      real ta_rexpocaco3, ta_rimpocaco3, ta_rdissl, ta_rcalatt
      real ta_romega_c
      common /ta_npzd_r/ ta_rnpp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_C(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpocaco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rimpocaco3(imt,km,jmt)
      common /ta_npzd_r/ ta_rdissl(imt,km,jmt)
      common /ta_npzd_r/ ta_romega_c(imt,km,jmt)
      common /ta_npzd_r/ ta_rcalatt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rprocal(imt,kpzd,jmt)
# else
      common /ta_npzd_r/ ta_rprocal(imt,jmt)
# endif
# if defined O_kk_ballast
      common /ta_npzd_r/ ta_rremi_B(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo_B(imt,km,jmt)
      common /ta_npzd_r/ ta_rgraz_Det_B(imt,kpzd,jmt)
# endif
# if defined O_mobi_silicon
      real ta_rmorpt_Diat
      real ta_rnpp_Diat, ta_rgraz_Diat, ta_rmorp_Diat
      common /ta_npzd_r/ ta_rmorp_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnpp_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_Diat(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt_Diat(imt,kpzd,jmt)
      real ta_rproopl, ta_rrivsil, ta_rexpoopl, ta_rdisopl
c      real ta_rbur_sil
      common /ta_npzd_r/ ta_rproopl(imt,km,jmt)
      common /ta_npzd_r/ ta_rrivsil(imt,jmt)
      common /ta_npzd_r/ ta_rexpoopl(imt,km,jmt)
      common /ta_npzd_r/ ta_rdisopl(imt,km,jmt)
c      common /ta_npzd_r/ ta_rbur_sil(imt,km,jmt)
# endif
#endif
#if defined O_save_convection
      integer nta_conv
      common /ta_conv_i/ nta_conv

      real ta_totalk, ta_vdepth, ta_pe
      common /ta_conv_r/ ta_totalk(imt,jmt), ta_vdepth(imt,jmt)
      common /ta_conv_r/ ta_pe(imt,jmt)
#endif
      integer nta_sscar
      common /ta_car_i/ nta_sscar
#if defined O_save_carbon_carbonate_chem
      real ta_sspH, ta_ssCO3, ta_ssOc, ta_ssOa, ta_sspCO2
      common /ta_car_r/ ta_sspH(imt,jmt), ta_ssCO3(imt,jmt)
      common /ta_car_r/ ta_ssOc(imt,jmt), ta_ssOa(imt,jmt)
      common /ta_car_r/ ta_sspCO2(imt,jmt)
#endif

#if defined O_mobi_iron
      real ta_rremife, ta_rexpofe
      common /ta_npzd_r/ ta_rremife(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexpofe(imt,kpzd,jmt)
# if defined O_save_mobi_diagnostics
      real ta_rfeorgads, ta_rdeffe, ta_rfecol
      real ta_rfeprime
      real ta_rfesed, ta_rbfe
      common /ta_npzd_r/ ta_rfeorgads(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfecol(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeffe(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfeprime(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rfesed(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rbfe(imt,kpzd,jmt)
#  if defined O_mobi_caco3
      real ta_rdeffe_C
      common /ta_npzd_r/ ta_rdeffe_C(imt,kpzd,jmt)
#  endif
#  if defined O_mobi_silicon
      real ta_rdeffe_Diat
      common /ta_npzd_r/ ta_rdeffe_Diat(imt,kpzd,jmt)
#  endif
# endif
#endif
