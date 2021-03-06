!======================== include file "isopyc.h" ======================
#if defined O_mom && defined O_isopycmix

!     isopycnal diffusion variables:

!     ahisop = isopycnal tracer mixing coefficient (cm**2/sec)
!     beta   = df/dy where f is related to coriolis force
!     drodx  = d(rho)/dx local to east face of T cell
!     drody  = d(rho)/dy local to north face of T cell
!     drodz  = d(rho)/dz local to bottom face of T cell
!     Ai_e   = diffusion coefficient on eastern face of T cell
!     Ai_n   = diffusion coefficient on northern face of T cell
!     Ai_bx  = diffusion coefficient on bottom face of T cell
!     Ai_by  = diffusion coefficient on bottom face of T cell

!     fisop  = structure function for isopycnal diffusion coefficient.
!     slmxr  = reciprocal of maximum allowable slope of isopycnals for
!              small angle approximation

      real alphai, betai, beta
      common /cisop_r/ alphai(imt,km,jmw), betai(imt,km,jmw), beta

      real addisop
      real ddxt, ddyt, ddzt, Ai_ez, Ai_nz, Ai_bx, Ai_by, K11, K22, K33
      real ahisop, fisop, slmxr, del_dm, s_dmr
# if defined O_full_tensor
      common /cisop_r/ ddxt(imt,km,jmw,2)
# else
      common /cisop_r/ ddxt(imt,km,jsmw:jemw,2)
# endif
      common /cisop_r/ ddyt(imt,km,1:jemw,2)
      common /cisop_r/ ddzt(imt,0:km,jmw,2)

      common /cisop_r/ Ai_ez(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ Ai_nz(imt,km,1:jemw,0:1,0:1)
      common /cisop_r/ Ai_bx(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ Ai_by(imt,km,jsmw:jemw,0:1,0:1)
      common /cisop_r/ K11(imt,km,jsmw:jemw)
      common /cisop_r/ K22(imt,km,1:jemw)
      common /cisop_r/ K33(imt,km,jsmw:jemw)
!     Changing fisop here... hopefully not a mistakE!
      common /cisop_r/ ahisop, fisop(imt,jmt,km), slmxr
      common /cisop_r/ addisop(imt,km,jsmw:jemw)
# if defined O_dm_taper
      common /cisop_r/ del_dm, s_dmr
# endif
      real delta_iso, s_minus, s_plus
      common /cisop_r/ delta_iso, s_minus, s_plus
# if defined O_gent_mcwilliams
!     adv_vetiso = zonal isopycnal mixing velocity computed at the
!                  center of the eastern face of the "t" cells
!     adv_vntiso = meridional isopycnal mixing velocity computed at
!                  the center of the northern face of the "t" cells
!                  (Note: this includes the cosine as in "adv_vnt")
!     adv_vbtiso = vertical isopycnal mixing velocity computed at the
!                  center of the top face of the "t" cells
!     adv_fbiso  = "adv_vbtiso" * (tracer) evaluated at the center of
!                  the bottom face of the "t" cells
!     athkdf = isopycnal thickness diffusivity (cm**2/sec)

      real athkdf, adv_vetiso, adv_vntiso, adv_vbtiso, adv_fbiso
      common /cisop_r/ athkdf
      common /cisop_r/ adv_vetiso(imt,km,jsmw:jemw)
      common /cisop_r/ adv_vntiso(imt,km,1:jemw)
      common /cisop_r/ adv_vbtiso(imt,0:km,jsmw:jemw)
      common /cisop_r/ adv_fbiso(imt,0:km,jsmw:jemw)
# endif

!     Define variables related to calculating K_gm mesoscale eddy 
!     diffiusivity as outlined in Gent an McWilliams Paper (1989). 

      real drodxte, drodxbe
      real drodytn, drodybn
      real drodzte, drodzbe
      real drodztn, drodzbn
      
      integer countx, county
      real abs_grd_rho2, grd_rho_x, grd_rho_y, abs_drho_dz
      
      common /cisop_r/ drodxte(imt,km,jmt), drodxbe(imt,km,jmt)
      common /cisop_r/ drodytn(imt,km,jmt), drodybn(imt,km,jmt)
      common /cisop_r/ drodzte(imt,km,jmt), drodzbe(imt,km,jmt)
      common /cisop_r/ drodztn(imt,km,jmt), drodzbn(imt,km,jmt)

# if defined O_KGMvar
      !     Oleg and Geoff
      !     niso = number of indices in kgm. =2 for anisotropic GM coeff

      integer niso
      parameter (niso = 1)

!     Further refinement from Eden 2009. 
!     LRhi       = Rhines scale. Defined as sigma/beta. Where sigma is 
!                  the Eady Growth rate of baroclinic instability
!     Lr         = 1st baroclinic Rossby Radius
!     Lr1        = 1st baroclinic Rossby Radius - placeholder
!     Lr2        = 1st baroclinic Rossby Radius - placeholder
!     c_eden     = Tuning parameter to ensure some O_KGM bounds from Eden 2009
!     gamma_eden = Tuning parameter used in sigma calculation from Eden 2009
!     kgm        = Thickness diffisivity constant
!     ahisop_var = Isopycnal diffusivity constant which derives from Kgm 
!     sigma_ave  = The average of Eady growth rate
      real Lm, Lr, Lr1, Lr2, LRhi, kgm
      real ahisop_var, gridsum_area
      real c_eden, coef, gamma_eden, pii
      real ahisop_ave, kgm_ave, sigma_ave
      real stratif_int, sum_zz
      real eddy_min, eddy_max

      common /kgm3d_r/ kgm(imt,km,jmt,niso), ahisop_var(imt,km,jmt,niso)
      common /kgm3d_r/ Lr(imt,jmt), LRhi(imt,km,jmt)
      common /cisop_r/ c_eden, gamma_eden

# endif  

      real drodxe, drodze, drodyn, drodzn, drodxb, drodyb, drodzb
      real drodye, drodxn

!     statement functions

      drodxe(i,k,j,ip) =    alphai(i+ip,k,j)*ddxt(i,k,j,1) +
     &                      betai(i+ip,k,j)*ddxt(i,k,j,2)
      drodze(i,k,j,ip,kr) = alphai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,1) +
     &                      betai(i+ip,k,j)*ddzt(i+ip,k-1+kr,j,2)

      drodyn(i,k,j,jq) =    alphai(i,k,j+jq)*ddyt(i,k,j,1) +
     &                      betai(i,k,j+jq)*ddyt(i,k,j,2)
      drodzn(i,k,j,jq,kr) = alphai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,1) +
     &                      betai(i,k,j+jq)*ddzt(i,k-1+kr,j+jq,2)

      drodxb(i,k,j,ip,kr) = alphai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,1) +
     &                      betai(i,k+kr,j)*ddxt(i-1+ip,k+kr,j,2)
      drodyb(i,k,j,jq,kr) = alphai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,1) +
     &                      betai(i,k+kr,j)*ddyt(i,k+kr,j-1+jq,2)
      drodzb(i,k,j,kr) =    alphai(i,k+kr,j)*ddzt(i,k,j,1) +
     &                      betai(i,k+kr,j)*ddzt(i,k,j,2)

# if defined O_full_tensor
      drodye(i,k,j,ip,jq) = alphai(i+ip,k,j)*ddyt(i+ip,k,j-1+jq,1) +
     &                      betai(i+ip,k,j)*ddyt(i+ip,k,j-1+jq,2)
      drodxn(i,k,j,ip,jq) = alphai(i,k,j+jq)*ddxt(i-1+ip,k,j+jq,1) +
     &                      betai(i,k,j+jq)*ddxt(i-1+ip,k,j+jq,2)
# endif

#endif
