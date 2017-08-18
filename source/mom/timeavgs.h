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
#if defined O_save_npzd
      real ta_rnpp, ta_rgraz, ta_rmorp, ta_rmorpt, ta_rmorz, ta_rexcr
      real ta_rremi, ta_rexpo, ta_rnpp_D, ta_rgraz_D, ta_rmorp_D
      real ta_rnfix, ta_rdeni, ta_rexpocal, ta_rprocal
      common /ta_npzd_r/ ta_rnpp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorpt(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorz(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rexcr(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rremi(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpo(imt,km,jmt)
      common /ta_npzd_r/ ta_rexpocal(imt,km,jmt)
      common /ta_npzd_r/ ta_rprocal(imt,jmt)
# if defined O_npzd_nitrogen
      common /ta_npzd_r/ ta_rnpp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rgraz_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rmorp_D(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rnfix(imt,kpzd,jmt)
      common /ta_npzd_r/ ta_rdeni(imt,km,jmt)
# endif
#endif
#if defined O_save_convection
      integer nta_conv
      common /ta_conv_i/ nta_conv

      real ta_totalk, ta_vdepth, ta_pe
      common /ta_conv_r/ ta_totalk(imt,jmt), ta_vdepth(imt,jmt)
      common /ta_conv_r/ ta_pe(imt,jmt)
#endif
#if defined O_save_carbon_carbonate_chem
      integer nta_sscar
      common /ta_car_i/ nta_sscar

      real ta_sspH, ta_ssCO3, ta_ssOc, ta_ssOa, ta_sspCO2
      common /ta_car_r/ ta_sspH(imt,jmt), ta_ssCO3(imt,jmt)
      common /ta_car_r/ ta_ssOc(imt,jmt), ta_ssOa(imt,jmt)
      common /ta_car_r/ ta_sspCO2(imt,jmt)
#endif
