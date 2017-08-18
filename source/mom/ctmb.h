!====================== include file "ctmb.h" ==========================

!                      meridional tracer balance

!     to calculate the meridional tracer balance, each term in the
!     tracer equation is integrated in depth,longitude and time for each
!     latitude. the time integration is over "tmbint"(see switch.h) days
!     for the calculation to make sense, the longitudinal
!     integration must go around the world (in the case of an open
!     latitude {eg: 60 deg S}) or between land masses (in the case of
!     a closed latitude {eg: 30 deg S from south america to africa})
!     the domain may be divided into arbitrary basins (atlantic, indian
!     pacific) but for certain latitudes, results from basins may have
!     to be combined to satisfy the above conditions.

!     ntmbb    = number of basins over which the tracer meridional
!                balance will be calculated. the test case assumes 1
!                basin as defined in "setocn.F"
!     tstor    = average latitudinal storage of tracer
!     tdiv     = average latitudinal divergence of tracer
!     tflux    = average latitudinal surface flux of tracer
!     tdif     = average latitudinal diffusion of tracer
!     tsorc    = average latitudinal additional source of tracer
!     smdvol   = ocean volume of basin at latitude
!     msktmb   = mask for basin numbers (1..ntmbb. 0 is reserved for
!                sum of all basins)
!     numtmb   = number of time steps over which the terms have been
!                accumulated
      integer ntmbb, msktmb, numtmb
      parameter (ntmbb=1)

      real tstor, tdiv, tflux, tdif, tsorc, smdvol

      common /tmb_r/ tstor(jmt,nt,0:ntmbb), tdiv(jmt,nt,0:ntmbb)
      common /tmb_r/ tflux(jmt,nt,0:ntmbb), tdif(jmt,nt,0:ntmbb)
      common /tmb_r/ tsorc(jmt,nt,0:ntmbb), smdvol(jmt,0:ntmbb)
      common /tmb_i/ msktmb(imt,jmt), numtmb
