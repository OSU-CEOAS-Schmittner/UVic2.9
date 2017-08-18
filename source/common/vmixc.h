!====================== include file "vmixc.h" =========================

!         vertical mixing coefficients and related variables

!     kappa_h = constant vertical diffusion coefficient (cm**2/sec)
!     kappa_m = constant vertical viscosity coefficient (cm**2/sec)

!     visc_cbu  = viscosity coeff at bottom of "u" cell
!     diff_cbt  = diffusion coeff at bottom of "T" cell
!     visc_cbu_limit = largest allowable "visc_cbu"
!     diff_cbt_limit = largest allowable "diff_cbt"
!     aidif = coefficient for implicit time differencing for
!             vertical diffusion. aidif=1 gives the fully implicit
!             case. aidif=0 gives the fully explicit case
!             note: not used unless "implicitvmix" or "isopycmix"
!                   is enabled
!=======================================================================

      real visc_cbu_limit, diff_cbt_limit, aidif, kappa_h, kappa_m
      real wndmix, fricmx, diff_cbt_back, visc_cbu_back, rhom1z
      real uzsq, diff_cbt, visc_cbt, visc_cbu, Ahv

      common /vmix_r/ visc_cbu_limit, diff_cbt_limit, aidif
      common /vmix_r/ kappa_h, kappa_m
#if !defined O_ppvmix
      common /vmix_r/ visc_cbu(imt,km,jsmw:jemw)
      common /vmix_r/ diff_cbt(imt,km,jsmw:jemw)
#endif

#if defined O_bryan_lewis_vertical
      common /vmix_r/ Ahv(km)
#endif

#if defined O_ppvmix

!     variables for pacanowski-philander vertical diffusion

!     rhom1z = rho(k)-rho(k+1)
!     uzsq   = (u(k)-u(k+1))**2 + (v(k)-v(k+1))**2
!     visc_cbt  = viscosity coeff at bottom of "T" cell
!     fricmx = max vertical mixing coefficient
!     wndmix = min vertical mixing in level 1 to simulate wind mixing
!     diff_cbt_back = background "diff_cbt"
!     visc_cbu_back = background "visc_cbu"

      common /vmix_r/ wndmix, fricmx, diff_cbt_back, visc_cbu_back
      common /vmix_r/ rhom1z(imt,km,jmw), uzsq(imt,km,jmw)
      common /vmix_r/ diff_cbt(imt,km,jsmw:jmw)
      common /vmix_r/ visc_cbt(imt,km,jsmw:jmw)
      common /vmix_r/ visc_cbu(imt,km,jsmw:jemw)
#endif
