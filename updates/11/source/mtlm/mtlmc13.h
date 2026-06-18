!======================== include file "mtlmc13.h" ========================
! rc13std = standard c13/c12 ratio
! ac13npp = fractionation factor "alpha" for NPP
! rc13a   = c13/c13 ratio of atmospheric CO2

      real rc13std
      parameter (rc13std=0.0112372)

      real ac13npp, rc13a

      common /land_c13_r/ ac13npp(NPFT), rc13a

!-----------------------------------------------------------------------
! IVM Prognostics
!-----------------------------------------------------------------------
! CS13   = Soil carbon 13 (kg C/m2).

      real CS13

      common /land_c13_r/ CS13(POINTS)

!-----------------------------------------------------------------------
! IVM Fluxes and diagnostics
!-----------------------------------------------------------------------
! C_VEG13     = Vegetation carbon (kg C/m2).
! CV13        = Gridbox mean vegetation carbon (kg C/m2).
! NPP13       = Net Primary Productivity (kg C/m2/s).
! RESP_S13    = Soil respiration rate (kg C/m2/s).

      real C_VEG13, CV13, NPP13, RESP_S13

      common /land_c13_r/ C_VEG13(POINTS,NPFT), CV13(POINTS)
      common /land_c13_r/ NPP13(POINTS,NPFT), RESP_S13(POINTS)

!-----------------------------------------------------------------------
! Driving variables for TRIFFID
!-----------------------------------------------------------------------
! NPP_DR13     = Accumulated Net Primary Productivity (kg C/m2/yr).
! RESP_S_DR13  = Accumulated soil respiration rate (kg C/m2/yr).

      real NPP_DR13, RESP_S_DR13

      common /land_c13_r/ NPP_DR13(POINTS,NPFT), RESP_S_DR13(POINTS)

# if defined O_time_averages
!-----------------------------------------------------------------------
! time averaged arrays
!-----------------------------------------------------------------------
! ta_CS13         = time averaged CS13
! ta_C_VEG13      = time averaged C_VEG13

      real ta_CS13, ta_C_VEG13
      real ta_RESP_S13, ta_BURN13, ta_NPP13

      common /land_c13_r/ ta_CS13(POINTS), ta_C_VEG13(POINTS,NPFT)
      common /land_c13_r/ ta_RESP_S13(POINTS), ta_BURN13(POINTS)
      common /land_c13_r/ ta_NPP13(POINTS,NPFT)

# endif
# if defined O_time_step_monitor
!-----------------------------------------------------------------------
! time average integrals
!-----------------------------------------------------------------------
! tai_CS13            = time averaged integral CS13
! tai_CV13            = time averaged integral CV13

      real tai_CS13, tai_CV13
      real tai_RESP_S13, tai_BURN13, tai_NPP13

      common /land_c13_r/ tai_CS13, tai_CV13
      common /land_c13_r/ tai_RESP_S13, tai_BURN13, tai_NPP13
# endif
