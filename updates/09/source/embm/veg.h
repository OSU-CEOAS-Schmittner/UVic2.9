!========================= include file "veg.h" ========================

!   variables for the vegetation data

!   Vegetation types
!     1 = tropical forest
!     2 = temperate/boreal forest
!     3 = grass
!     4 = shrub
!     5 = tundra
!     6 = desert
!     7 = ice

!     nveg      = number of vegetation classes
!     iveg      = vegetation class
!     veg_rl    = roughness length
!     veg_alb   = albedo
!     veg_rs    = stomatal resistance
!     veg_smd   = snow masking depth (m)
!     veg_dalt  = dalton number over land
!     idesert   = index for desert
!     iagric    = index for agricultural land
!     icrops    = this variable is obsolete and its value ignored
!   land use ("crop_data" "pasture_data" "agric_data")
!     agric     = agricultural extent as percentage of grid cell

      integer nveg, iveg, idesert, iagric, icrops, iice
      real veg_rl, veg_alb, veg_rs, veg_smd, veg_dalt, agric

      parameter (nveg=7)
      common /veg_i/ iveg(imt,jmt), idesert, iagric, icrops, iice

      common /veg_r/ veg_rl(nveg), veg_alb(nveg), veg_rs(nveg)
      common /veg_r/ veg_smd(nveg)
#if defined O_crop_data || defined O_crop_data_transient || defined O_pasture_data || defined O_agric_data || defined O_agric_data_transient
      common /veg_r/ veg_dalt(imt,jmt)
      common /veg_r/ agric(imt,jmt,1:3)
#else
      common /veg_r/ veg_dalt(nveg)
#endif
