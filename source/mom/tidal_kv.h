!===================== include file "tidal_kv.h" =======================

!     quantities for parameterization of tidal mixing

!     Simmons, H.L., S.R. Jayne, L.C. St. Laurent and A.J. Weaver, 2004:
!     Tidally driven mixing in a numerical model of the ocean general
!     circulation. Ocean Modelling, 6, 245-263.

# if defined O_tidal_kv

      real zetar, ogamma, gravrho0r
      common /tdr/ zetar, ogamma, gravrho0r

!     edr = energy dissipation rate due to tides

      real edr
      common /edr/ edr(imt,jmt)
#endif
