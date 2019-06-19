!====================== include file "fdifm.h" =========================

!     finite difference numerics for momentum
!=======================================================================

!-----------------------------------------------------------------------
!     advective terms
!-----------------------------------------------------------------------

#if defined O_linearized_advection
      ADV_Ux(i,k,j) = 0.0
      ADV_Uy(i,k,j,jrow,n) = 0.0
      ADV_Uz(i,k,j) = 0.0
      ADV_metric(i,k,j,jrow,n) = 0.0
#else
      ADV_Ux(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*csudxu2r(i,j)
      ADV_Uy(i,k,j,jrow,n) = (adv_vnu(i,k,j)*(u(i,k,j,n,tau)
     &  + u(i,k,j+1,n,tau)) - adv_vnu(i,k,j-1)*(u(i,k,j-1,n,tau)
     &  + u(i,k,j,n,tau)))*csudyu2r(jrow)
      ADV_Uz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
      ADV_metric(i,k,j,jrow,n) = advmet(jrow,n)*u(i,k,j,1,tau)
     &  *u(i,k,j,3-n,tau)
#endif

!-----------------------------------------------------------------------
!     diffusive terms
!-----------------------------------------------------------------------

      DIFF_Ux(i,k,j) = (diff_fe(i,k,j) - diff_fe(i-1,k,j))
     &  *csudxur(i,j)
      DIFF_Uz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
#if defined O_implicitvmix
     &  *(c1-aidif)
#endif
#if defined O_consthmix && !defined O_biharmonic
# if defined O_neptune
      DIFF_Uy(i,k,j,jrow,n) = amc_north(jrow)*(u(i,k,j+1,n,taum1)
     &  - u(i,k,j,n,taum1) - unep(i,jrow+1,n)*umask(i,k,j+1)
     $  + unep(i,jrow,n)*umask(i,k,j)) - amc_south(jrow)
     &  *(u(i,k,j,n,taum1) - u(i,k,j-1,n,taum1) - unep(i,jrow,n)
     $  *umask(i,k,j) + unep(i,jrow-1,n)*umask(i,k,j-1))
# else
#  if defined O_anisotropic_viscosity
      DIFF_Uy(i,k,j,jrow,n) = amc_north(i,k,jrow)*(u(i,k,j+1,n,taum1)
     &  - u(i,k,j,n,taum1)) - amc_south(i,k,jrow)*(u(i,k,j,n,taum1)
     &  - u(i,k,j-1,n,taum1))
#  else
      DIFF_Uy(i,k,j,jrow,n) = amc_north(jrow)*(u(i,k,j+1,n,taum1)
     &  - u(i,k,j,n,taum1)) - amc_south(jrow)*(u(i,k,j,n,taum1)
     &  - u(i,k,j-1,n,taum1))
#  endif
# endif
#else
      DIFF_Uy(i,k,j,jrow,n) = (diff_fn(i,k,j) - diff_fn(i,k,j-1))
     &  *csudyur(jrow)
#endif

!-----------------------------------------------------------------------
!     metric term
!-----------------------------------------------------------------------

#if defined O_consthmix
# if defined O_biharmonic
      DIFF_metric(i,k,j,jrow,n) = am3(jrow)*del2(i,k,j,n) + am4(jrow,n)
     &  *(del2(i+1,k,j,3-n) - del2(i-1,k,j,3-n))*dxmetr(i)
# else
#  if defined O_neptune
      DIFF_metric(i,k,j,jrow,n) = am3(jrow)*(u(i,k,j,n,taum1)
     $  - unep(i,jrow,n)*umask(i,k,j)) + am4(jrow,n)*dxmetr(i)
     &  *(u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1)
     $  - unep(i+1,jrow,3-n)*umask(i+1,k,j) + unep(i-1,jrow,3-n)
     $  *umask(i-1,k,j))
#  else
      DIFF_metric(i,k,j,jrow,n) = am3(jrow)*u(i,k,j,n,taum1)
     &  + am4(jrow,n)*dxmetr(i)*(u(i+1,k,j,3-n,taum1)
     &  - u(i-1,k,j,3-n,taum1))
#  endif
# endif
#else
# if defined O_smagnlmix
      DIFF_metric(i,k,j,jrow,n) = smag_metric(i,k,j)
# endif
#endif

!-----------------------------------------------------------------------
!     coriolis term
!-----------------------------------------------------------------------

#if defined O_stream_function
# if defined O_damp_inertial_oscillation
      CORIOLIS(i,k,j,jrow,n) = (c1-acor)*cori(i,jrow,n)
     &  *u(i,k,j,3-n,taum1)
# else
      CORIOLIS(i,k,j,jrow,n) = cori(i,jrow,n)*u(i,k,j,3-n,tau)
# endif
#else
      CORIOLIS(i,k,j,jrow,n) = cori(i,jrow,n)*(gcor*u(i,k,j,3-n,tau)
     &  + (c1-gcor)*u(i,k,j,3-n,taum1))
#endif
