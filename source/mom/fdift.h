!====================== include file "fdift.h" =========================

!     finite difference numerics for tracers
!=======================================================================

      T_i(i,k,j,n,ip) = t(i+ip,k,j,n,taum1)
      T_j(i,k,j,n,jp) = t(i,k,j+jp,n,taum1)
      dz_t2r(i,k,j) = dzt2r(k)
      dz_tr(i,k,j)  = dztr(k)
      dz_wtr(i,k,j) = dzwr(k)
      dx_t2r(i,k,j) = cstdxt2r(i,j)
      dx_tr(i,k,j)  = cstdxtr(i,j)
      dy_t2r(i,k,j) = cstdyt2r(jrow)
      dy_tr(i,k,j)  = cstdytr(jrow)

!-----------------------------------------------------------------------
!     advective terms
!-----------------------------------------------------------------------

#if defined O_linearized_advection
      ADV_Tx(i,k,j) = 0.0
      ADV_Ty(i,k,j,jrow,n) = 0.0
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
#else
      ADV_Tx(i,k,j) = (adv_fe(i,k,j) - adv_fe(i-1,k,j))*cstdxt2r(i,j)
# if defined O_fourth_order_tracer_advection || defined O_quicker
      ADV_Ty(i,k,j,jrow,n) = (adv_f4n(i,k,j,n) - adv_f4n(i,k,j-1,n))
     &  *cstdyt2r(jrow)
# else
#  if defined O_fct
      ADV_Ty(i,k,j,jrow,n) = (adv_fn(i,k,j) - adv_fn(i,k,j-1))
     &  *cstdyt2r(jrow)
#  else
      ADV_Ty(i,k,j,jrow,n) = (adv_vnt(i,k,j)*(t(i,k,j,n,tau)
     &  + t(i,k,j+1,n,tau)) - adv_vnt(i,k,j-1)*(t(i,k,j-1,n,tau)
     &  + t(i,k,j,n,tau)))*cstdyt2r(jrow)
#  endif
# endif
      ADV_Tz(i,k,j) = (adv_fb(i,k-1,j) - adv_fb(i,k,j))*dzt2r(k)
#endif
#if defined O_gent_mcwilliams && defined O_isopycmix

!     gent_mcwilliams isopycnal advective terms simulating the effect
!     of eddies on the isopycnals

      ADV_Txiso(i,k,j,n) = cstdxt2r(i,j)*(adv_vetiso(i,k,j)
     &  *(t(i+1,k,j,n,taum1) + t(i,k,j,n,taum1)) - adv_vetiso(i-1,k,j)
     &  *(t(i,k,j,n,taum1) + t(i-1,k,j,n,taum1)))
      ADV_Tyiso(i,k,j,jrow,n) = cstdyt2r(jrow)*(adv_vntiso(i,k,j)
     &  *(t(i,k,j+1,n,taum1) + t(i,k,j,n,taum1)) - adv_vntiso(i,k,j-1)
     &  *(t(i,k,j,n,taum1) + t(i,k,j-1,n,taum1)))
      ADV_Tziso(i,k,j) = dzt2r(k)*(adv_fbiso(i,k-1,j)-adv_fbiso(i,k,j))
#endif

!-----------------------------------------------------------------------
!     diffusive terms
!-----------------------------------------------------------------------

!     zonal component

      DIFF_Tx(i,k,j) = (diff_fe(i,  k,j)*tmask(i+1,k,j)
     &  - diff_fe(i-1,k,j)*tmask(i-1,k,j))*cstdxtr(i,j)

!     meridional component

#if defined O_consthmix && !defined O_biharmonic && !defined O_isopycmix
# if defined O_bryan_lewis_horizontal
      DIFF_Ty(i,k,j,jrow,n) = ahc_north(jrow,k)*tmask(i,k,j+1)
     &  *(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1)) - ahc_south(jrow,k)
     &  *tmask(i,k,j-1)*(t(i,k,j,n,taum1) - t(i,k,j-1,n,taum1))
# else
      DIFF_Ty(i,k,j,jrow,n) = ahc_north(jrow)*tmask(i,k,j+1)
     &  *(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1)) - ahc_south(jrow)
     &  *tmask(i,k,j-1)*(t(i,k,j,n,taum1) - t(i,k,j-1,n,taum1))
# endif
#else
      DIFF_Ty(i,k,j,jrow,n) = (diff_fn(i,k,j  )*tmask(i,k,j+1)
     &  - diff_fn(i,k,j-1)*tmask(i,k,j-1))*cstdytr(jrow)
#endif

!     vertical component

      DIFF_Tz(i,k,j) = (diff_fb(i,k-1,j) - diff_fb(i,k,j))*dztr(k)
#if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
     &  *(c1-aidif)
#endif
#if defined O_isopycmix
     &  + (diff_fbiso(i,k-1,j) - diff_fbiso(i,k,j))*dztr(k)
#endif
