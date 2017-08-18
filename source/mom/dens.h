!====================== include file "dens.h" ==========================

!-----------------------------------------------------------------------
!     statement function
!-----------------------------------------------------------------------

#if defined O_linearized_density
!     approximating rho = 1.035*(1-alpha*tq). The one is removed by
!     gradients and the 1.035 is absorbed into alpha.

      dens(tq,sq,k) = -2.e-4*tq
#else
      dens (tq, sq, k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     &                   (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     &                   (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq
#endif
#if defined O_isopycmix
      drodt (tq, sq, k) = c(k,1) + (c(k,4) + c(k,7)*sq)*sq + (2.0*c(k,3)
     &                  + 2.0*c(k,8)*sq + 3.0*c(k,6)*tq)*tq

      drods (tq, sq, k) = (c(k,4) + 2.0*c(k,7)*sq + c(k,8)*tq)*tq
     &                  + c(k,2) + (2.0*c(k,5) + 3.0*c(k,9)*sq)*sq
#endif
