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
!     icrops    = index for cropland
!   land use ("crop_data")
!     crops     = cropland extent as percentage of grid cell
!     crops_yr  = year for crops (1700 to 1992)

      integer nveg, iveg, idesert, icrops, iice
      real veg_rl, veg_alb, veg_rs, veg_smd, veg_dalt, crops, crops_yr

      parameter (nveg=7)
      common /veg_r/ veg_rl(nveg), veg_alb(nveg), veg_rs(nveg)
      common /veg_r/ veg_smd(nveg)
      common /veg_i/ iveg(imt,jmt), idesert, icrops, iice
#if defined O_crop_data
      common /veg_r/ crops(imt,jmt,1:3), crops_yr
      common /veg_r/ veg_dalt(imt,jmt)
#else
      common /veg_r/ veg_dalt(nveg)
#endif
