!======================== include file "mtlmc14.h" ========================
! rc14std = standard c14/c12 ratio
! ac14npp = fractionation factor "alpha" for NPP
! rc14a   = c14/c14 ratio of atmospheric CO2

!     rc14std   = standard c14/c12 ratio
      real rc14std
      parameter (rc14std=1.176e-12)

      real ac14npp, rc14a

      common /land_c14_r/ ac14npp(NPFT), rc14a

!-----------------------------------------------------------------------
! IVM Prognostics
!-----------------------------------------------------------------------
! CS14   = Soil carbon 14 (kg C/m2).

      real CS14

      common /land_c14_r/ CS14(POINTS)

!-----------------------------------------------------------------------
! IVM Fluxes and diagnostics
!-----------------------------------------------------------------------
! C_VEG14     = Vegetation carbon (kg C/m2).
! CV14        = Gridbox mean vegetation carbon (kg C/m2).
! NPP14       = Net Primary Productivity (kg C/m2/s).
! RESP_S14    = Soil respiration rate (kg C/m2/s).

      real C_VEG14, CV14, NPP14, RESP_S14

      common /land_c14_r/ C_VEG14(POINTS,NPFT), CV14(POINTS)
      common /land_c14_r/ NPP14(POINTS,NPFT), RESP_S14(POINTS)

!-----------------------------------------------------------------------
! Driving variables for TRIFFID
!-----------------------------------------------------------------------
! NPP_DR14     = Accumulated Net Primary Productivity (kg C/m2/yr).
! RESP_S_DR14  = Accumulated soil respiration rate (kg C/m2/yr).

      real NPP_DR14, RESP_S_DR14

      common /land_c14_r/ NPP_DR14(POINTS,NPFT), RESP_S_DR14(POINTS)

# if defined O_time_averages
!-----------------------------------------------------------------------
! time averaged arrays
!-----------------------------------------------------------------------
! ta_CS14         = time averaged CS14
! ta_C_VEG14      = time averaged C_VEG14

      real ta_CS14, ta_C_VEG14
      real ta_RESP_S14, ta_BURN14, ta_NPP14

      common /land_c14_r/ ta_CS14(POINTS), ta_C_VEG14(POINTS,NPFT)
      common /land_c14_r/ ta_RESP_S14(POINTS), ta_BURN14(POINTS)
      common /land_c14_r/ ta_NPP14(POINTS,NPFT)

# endif
# if defined O_time_step_monitor
!-----------------------------------------------------------------------
! time average integrals
!-----------------------------------------------------------------------
! tai_CS14            = time averaged integral CS14
! tai_CV14            = time averaged integral CV14
! tai_clnd14          = average total carbon 14 in land
! tai_cfa2l14         = average total c14 flux atmosphere to land

      real tai_CS14, tai_CV14
      real tai_RESP_S14, tai_BURN14, tai_NPP14

      common /land_c14_r/ tai_CS14, tai_CV14
      common /land_c14_r/ tai_RESP_S14, tai_BURN14, tai_NPP14
# endif
