!======================== include file "riv.h" =========================

!     parameters for use in the river model

!     maxnb   = maximum number of basins
!     nb      = number of basins
!     nbp     = number of basin points per basin
!     ndp     = number of discharge points per basin
!     ndis    = map of discharge points
!     nrfill  = map of filled river numbers (extrapolated over water)
!     nriv    = map of river (basin) numbers
!     psum    = total discharge for a basin
!     wdar    = discharge weights over discharge area
!     ta_psum = time average total discharge for a basin

      integer maxnb
      parameter (maxnb=200)

      integer nb, nbp, ndp, ndis, nrfill, nriv
      common /river_i/ nb, nbp(maxnb), ndp(maxnb), ndis(imt,jmt)
      common /river_i/ nrfill(imt,jmt), nriv(imt,jmt)

      real psum, wdar, ta_psum
      common /river_r/ psum(maxnb), wdar(imt,jmt)
#if defined O_time_averages
      common /river_r/ ta_psum(maxnb)
#endif
