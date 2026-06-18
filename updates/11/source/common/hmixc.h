!======================= include file "hmixc.h" ========================

!                    horizontal mixing coefficients

!     visc_cnu = viscosity coeff for northern face of "u" cell
!     visc_ceu = viscosity coeff for eastern face of "u" cell
!     diff_cnt = diffusion coeff for northern face of "T" cell
!     diff_cet = diffusion coeff for eastern face of "T" cell

!     am     = constant lateral viscosity coeff for momentum
!     ah     = constant lateral diffusion coeff for tracers
!     am3    = viscosity coeff for metric term on "u" cell
!     am4    = another viscosity coeff for metric term on "u" cell
!     ambi   = constant lateral biharmonic viscosity coeff for momentum
!     ahbi   = constant lateral biharmonic diffusion coeff for tracers
!=======================================================================

      real am, ambi, am3, am4, ah, ahbi, visc_ceu, visc_cnu, amc_north
      real amc_south, Ahh(km), diff_cnt, diff_cet, ahc_north, ahc_south
      real strain, am_lambda, am_phi, smag_metric, diff_c_back
      real hl_depth, hl_back, hl_max, hl_u, hl_n, hl_e, hl_b
      real droz, rich_inv

#if defined O_consthmix
      common /diffus_r/ am, ambi, am3(jmt), am4(jmt,2)
      common /diffus_r/ ah, ahbi
# if defined O_anisotropic_viscosity
      common /diffus_r/ visc_ceu(imt,km,jmt)
      common /diffus_r/ visc_cnu(imt,km,jmt)
      common /diffus_r/ amc_north(imt,km,jmt)
      common /diffus_r/ amc_south(imt,km,jmt)
# else
      common /diffus_r/ visc_ceu, visc_cnu
      common /diffus_r/ amc_north(jmt), amc_south(jmt)
# endif
# if defined O_bryan_lewis_horizontal

!     bryan_lewis mixing case

      common /diffus_r/ Ahh
      common /diffus_r/ diff_cnt(km), diff_cet(km)
      common /diffus_r/ ahc_north(jmt,km), ahc_south(jmt,km)
# else
      common /diffus_r/ diff_cnt, diff_cet
      common /diffus_r/ ahc_north(jmt), ahc_south(jmt)
# endif
#else
# if defined O_smagnlmix

!     non-linear horizontal viscosity after Smagorinsky 1963,
!     as described in Rosati & Miyakoda (jpo,vol 18,#11,1988)
!     see Smagorinsky 1963, Mon Wea Rev, 91, 99-164.
!     Also see Deardorff 1973 J. Fluid Eng. Sep., 429-438.

!     strain = tension(1) and shearing(2) rates of strain
!     smag_metric  = metric term
!     diff_c_back = background diffusion coeff for "t" cell (cm**2/sec)

      common /diffus_r/ strain(imt,km,1:jemw,2)
      common /diffus_r/ am_lambda(imt,km,1:jemw), am_phi(imt,km,1:jemw)
      common /diffus_r/ smag_metric(imt,km,jsmw:jemw)
      common /diffus_r/ diff_c_back
      common /diffus_r/ visc_ceu(imt,km,jsmw:jemw)
      common /diffus_r/ visc_cnu(imt,km,1:jemw)
      common /diffus_r/ diff_cet(imt,km,jsmw:jemw)
      common /diffus_r/ diff_cnt(imt,km,1:jemw)
# endif
#endif
