!===================== include file "cnep.h" ===========================

!     the option "neptune" provides a subgridscale parameterization
!     for the interaction of eddies and topography

!     reference:
!       Holloway, G., 1992: representing topographic stress for large
!     scale ocean models, J. Phys. Oceanogr., 22, 1033-1046

!     neptune is calculated as an equilibrium streamfunction given by
!     pnep=-f*snep*snep*hnep and is applied through eddy viscosity

!     hnep = model streamfunction depth
!     snep = spnep + (senep-spnep)*(0.5 + 0.5*cos(2.0*latitude))
!       the neptune length scale snep has a value of senep at the
!       equator and smoothly changes to spnep at the poles

!     variables used in applying neptune

!     spnep = neptune length scale at the pole
!     senep = neptune length scale at the equator
!     unep  = neptune velocity
!     pnep  = neptune streamfunction

      real spnep, senep, pnep, unep
      common /cnep_r/ spnep, senep, pnep(imt,jmt), unep(imt,jmt,2)
