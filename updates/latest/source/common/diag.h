!====================== include file "diag.h" ==========================
#if defined O_mom

!     variables used for computing diagnostics:

# if defined O_time_step_monitor

!     ektot     = "total" kinetic energy per unit volume at "tau". units
!                 ergs/cm**3 = dyn/cm**2 = g/cm/sec**2 = 10**-7 J/cm**3.
!                 ektot is the "total" ke in the sense that it considers
!                 both the internal and external modes summed over the
!                 entire ocean volume. The contributions of
!                 vertical motions are neglected on the basis of scaling
!                 arguments (i.e., w**2 << (u**2 + v**2).
!     dtabs     = absolute value of rate of change of tracer per unit
!                 volume centered at "tau"
!     tbar      = first moment of tracer at "tau"
!     travar    = variance = second moment of tracer about mean at "tau"
!     isot1     = starting i index of section 1 for max/min overturning
!     ieot1     = ending i index of section 1 for max/min overturning
!     isot2     = starting i index of section 2 for max/min overturning
!     ieot2     = ending i index of section 2 for max/min overturning
!     jsot      = starting j index for max/min overturning
!     jeot      = ending j index for max/min overturning
!     ksot      = starting k index for max/min overturning
!     keot      = ending k index for max/min overturning
!     mrot      = regional mask region for max/min overturning
!     v_otsf    = velocity field for calculating max/min overturning
!     t_slh     = tracer fields for calculating sea level height
!     d_slh     = density field for calculating sea level height

      real ektot, dtabs, travar, tbar
      common /cdiag_r/ ektot(0:km,jmt), dtabs(0:km,nt,jmt)
      common /cdiag_r/ travar(0:km,nt,jmt), tbar(0:km,nt,jmt)

      integer isot1, ieot1, isot2, ieot2, jsot, jeot, ksot, keot, mrot
      common /cdiag_i/ isot1, ieot1, isot2, ieot2, jsot, jeot
      common /cdiag_i/ ksot, keot, mrot
#  if defined O_tai_otsf
      integer nv_otsf
      common /cdiag_i/ nv_otsf

      real v_otsf
      common /cdiag_r/ v_otsf(jmt,km)
#  endif
#  if defined O_tai_slh
      integer nt_slh
      common /cdiag_i/ nt_slh

      real t_slh, d_slh
      common /cdiag_r/ t_slh(imt,jmt,km,2), d_slh(imt,jmt,km)
#  endif

!     ntatio        = number of time averaged time step integrals
!     tai_ek        = average integrated kinetic energy
!     tai_t         = average integrated temperature
!     tai_s         = average integrated salinity
!     tai_tvar      = average integrated second moment of temperature
!     tai_svar      = average integrated second moment of salinity
!     tai_dt        = average integrated rate of change of temperature
!     tai_ds        = average integrated rate of change of salinity
!     tai_scan      = average scans for ocean solver
!     tai_otmax     = average maximum overturning
!     tai_otmin     = average minimum overturning
!     tai_slh       = average integrated sea level height
!     tai_hflx      = average heat flux
!     tai_sflx      = average salt flux
!     tai_dic       = average carbon
!     tai_dicflx    = average carbon flux
!     tai_alk       = average alkalinity
!     tai_o2        = average oxygen
!     tai_o2flx     = average oxygen flux
!     tai_po4       = average phosphate
!     tai_dop       = average dop
!     tai_phyt      = average phytoplankton
!     tai_phyt_phos      = average phytoplankton phosphorous
!     tai_zoop      = average zooplankton
!     tai_detr      = average detritus
!     tai_detr_phos      = average detritus phosphorous
!     tai_no3       = average nitrate
!     tai_don       = average don
!     tai_diaz      = average diazotrophs
!     tai_din15     = average nitrate 15
!     tai_don15     = average don15
!     tai_phytn15   = average phytoplankton n15
!     tai_diatn15   = average diatom n15			     
!     tai_zoopn15   = average zooplankton n15
!     tai_detrn15   = average detritus n15
!     tai_diazn15   = average diazotrophs n15
!     tai_dic13     = average DIC13
!     tai_doc13     = average DOC13
!     tai_phytc13   = average phytoplankton c13
!     tai_diatc13   = average diatom c13			     
!     tai_zoopc13   = average zooplankton c13
!     tai_detrc13   = average detritus c13
!     tai_diazc13   = average diazotrophs c13
!     tai_c14       = average carbon 14
!     tai_dc14      = average delta carbon 14
!     tai_c14flx    = average carbon 14 flux
!     tai_cfc11     = average CFC11
!     tai_cfc11flx  = average CFC11 flux
!     tai_cfc12     = average CFC12
!     tai_cfc12flx  = average CFC12 flux
!     tai_sspH      = average sea surface pH
!     tai_ssCO3     = average sea surface CO3
!     tai_ssOc      = average sea surface calcite
!     tai_ssOa      = average sea surface omega aragonite
!     tai_sspCO2    = average sea surface pCO2
!     tai_cocn      = average total carbon in ocean
!     tai_cfa2o     = average total flux atmosphere to ocean
!     tai_dicwflx   = average carbon flux from weathering
!     tai_dfe       = average iron
!     tai_ddfe      = average particulate fe
!     tai_caco3     = average detached calcite
!     tai_caco3c13  = average detached calcite c13
!     tai_diat      = average diatom
!     tai_sil       = average silica
!     tai_silflx    = average silica flux
!     tai_opl       = average opal
			     
      integer ntatio
      common /cdiagi/ ntatio

      real tai_ek, tai_t, tai_s, tai_tvar, tai_svar, tai_dt
      real tai_ds, tai_scan, tai_otmax, tai_otmin, tai_slh
      real tai_hflx, tai_sflx, tai_dic, tai_dicflx, tai_alk
      real tai_o2, tai_o2flx, tai_po4, tai_phyt, tai_zoop, tai_detr
      real tai_no3, tai_diaz, tai_c14, tai_dc14, tai_c14flx
      real tai_phyt_phos, tai_detr_phos	
      real tai_cfc11, tai_cfc11flx, tai_cfc12, tai_cfc12flx
      real tai_sspH, tai_ssCO3, tai_ssOc, tai_ssOa, tai_sspCO2
      real tai_cocn, tai_cfa2o, tai_dicwflx
      real tai_dop, tai_don
      real tai_din15, tai_don15, tai_phytn15, tai_zoopn15, tai_detrn15
      real tai_diazn15, tai_diatn15, tai_diatc13
      real tai_dic13, tai_doc13, tai_phytc13, tai_zoopc13, tai_detrc13
      real tai_diazc13, tai_dic13flx, tai_dic13wflx, tai_dfe, tai_ddfe
      real tai_d_B, tai_c, tai_caco3, tai_caco3c13, tai_calatt, tcalatt
			     
      common /cdiag_r/ tai_ek, tai_t, tai_s, tai_tvar, tai_svar, tai_dt
      common /cdiag_r/ tai_ds, tai_scan, tai_otmax, tai_otmin, tai_slh
      common /cdiag_r/ tai_hflx, tai_sflx, tai_dic, tai_dicflx, tai_alk
      common /cdiag_r/ tai_o2, tai_o2flx, tai_po4, tai_phyt, tai_zoop
      common /cdiag_r/ tai_no3, tai_diaz, tai_c14, tai_dc14, tai_c14flx
      common /cdiag_r/ tai_cfc11, tai_cfc11flx, tai_cfc12, tai_cfc12flx
      common /cdiag_r/ tai_sspH, tai_ssCO3, tai_ssOc, tai_ssOa
      common /cdiag_r/ tai_sspCO2, tai_cocn, tai_cfa2o, tai_dicwflx
      common /cdiag_r/ tai_dop, tai_don, tai_diatn15
      common /cdiag_r/ tai_din15, tai_don15, tai_phytn15, tai_zoopn15
      common /cdiag_r/ tai_detrn15, tai_diazn15, tai_detr
      common /cdiag_r/ tai_phyt_phos, tai_detr_phos
      common /cdiag_r/ tai_dic13, tai_doc13, tai_phytc13, tai_zoopc13
      common /cdiag_r/ tai_detrc13, tai_diazc13, tai_dic13wflx, tai_dfe
      common /cdiag_r/ tai_dic13flx, tai_ddfe, tai_diatc13, tai_caco3c13
      common /cdiag_r/ tai_d_B, tai_c, tai_caco3, tai_calatt
      common /cdiag_r/ tcalatt(0:km,jmt)	
      real tai_diat
      common /cdiag_r/ tai_diat
      real tai_sil, tai_silflx, tai_opl
      common /cdiag_r/ tai_sil, tai_silflx, tai_opl
# endif
# if defined O_energy_analysis

!     engint    = volume averaged internal mode energy integral
!                 components
!     engext    = volume averaged external mode energy integral
!                 components
!     buoy      = volume averaged buoyancy

!     tcerr     = maximum "t" cell continuity error
!     ucerr     = maximum "u" cell continuity error
!     itcerr    = "i" index corresponding to "tcerr"
!     jtcerr    = "jrow" index corresponding to "tcerr"
!     ktcerr    = "k" index corresponding to "tcerr"
!     iucerr    = "i" index corresponding to "ucerr"
!     jucerr    = "jrow" index corresponding to "ucerr"
!     kucerr    = "k" index corresponding to "ucerr"

!     wtbot     = maximum "adv_vbt" error at ocean bottom
!     iwtbot    = "i" index corresponding to "wtbot"
!     jwtbot    = "jrow" index corresponding to "wtbot"
!     kwtbot    = "k" index corresponding to "wtbot"
!     wubot     = maximum "adv_vbu" at ocean bottom
!     iwubot    = "i" index corresponding to "wubot"
!     jwubot    = "jrow" index corresponding to "wubot"
!     kwubot    = "k" index corresponding to "wubot"

!     wtlev     = zonally integrated adv_vbt for each level
!     wulev     = zonally integrated adv_vbu for each level

      integer itcerr, jtcerr, ktcerr, iucerr, jucerr, kucerr
      integer iwtbot, jwtbot, kwtbot, iwubot, jwubot, kwubot
      common /cdiag_i/ itcerr(jmt), jtcerr(jmt), ktcerr(jmt)
      common /cdiag_i/ iucerr(jmt), jucerr(jmt), kucerr(jmt)
      common /cdiag_i/ iwtbot(jmt), jwtbot(jmt), kwtbot(jmt)
      common /cdiag_i/ iwubot(jmt), jwubot(jmt), kwubot(jmt)

      real buoy, engint, engext, tcerr, ucerr, wtbot, wubot, wtlev
      real wulev

      common /cdiag_r/ buoy(0:km,jmt), engint(0:km,8,jmt), engext(8,jmt)
      common /cdiag_r/ tcerr(jmt), ucerr(jmt)
      common /cdiag_r/ wtbot(jmt), wubot(jmt)
      common /cdiag_r/ wtlev(km,0:jmt), wulev(km,0:jmt)
# endif
# if defined O_gyre_components

!     ttn       = northward transport of tracer components

!     ttn2      = northward transport of tracers for ocean basins
!                  (.,.,.,0)       Global
!                  (.,.,.,1:nhreg) Ocean basins
!                also,
!                  (6,.,.,.) total transport due to advection
!                  (7,.,.,.) total transport due to diffusion
!                  (8,.,.,.) total transport
#  if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
!                  (9,.,.,.) total transport due to isopycnal advection
#  endif
      real ttn
      common /gyres_r/ ttn(8,jmt,ntmin2)
#  if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
      real ttn2
      common /gyres_r/ ttn2(6:9,jmt,nt,0:nhreg)
#  else
      real ttn2
      common /gyres_r/ ttn2(6:8,jmt,nt,0:nhreg)
#  endif
# endif
# if defined O_meridional_overturning

!     vmsf      = vertical_meridional stream function

      real vmsf
      common /cdiag_r/ vmsf(jmt,km)
# endif
# if defined O_show_zonal_mean_of_sbc

!     zmsmf     = zonal mean surface momentum flux
!     zmstf     = zonal mean surface tracer flux
!     zmsm      = zonal mean surface momentum
!     zmst      = zonal mean surface tracers
!     zmau      = surface area weighting for "u" latitudes
!     zmat      = surface area weighting for "t" latitudes

      real zmsmf, zmstf, zmau, zmat, zmsm, zmst
      common /cdiag/ zmsmf(jmt,2), zmstf(jmt,nt), zmau(jmt), zmat(jmt)
      common /cdiag/ zmsm(jmt,2), zmst(jmt,nt)
# endif

# if defined O_tracer_yz

!     tyz(,,,1) = zonal mean of tracer T
!     tyz(,,,2) = zonal mean of d(T)/dt
!     tyz(,,,3) = zonal mean of advection of T
!     tyz(,,,4) = zonal mean of diffusion of T
!     tyz(,,,5) = zonal mean of source of T

      real tyz
      common /cdiag/ tyz(jmt,km,nt,5)
# endif
# if defined O_term_balances

!     term balances are instantaneous breakdowns of all terms in the
!     momentum & tracer equations. They are averaged over ocean volumes
!     defined by horizontal and vertical regional masks:

!     termbt   = term balance components for time rate of change of
!                tracers within a  volume. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (10) where each
!                term has the units of "tracer units/sec" using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the tracer
!                (2)  = change due to zonal nonlinear term: (UT)x
!                (3)  = change due to meridional nonlinear term: (VT)y
!                (4)  = change due to vertical nonlinear term: (WT)z
!                (5)  = change due to zonal diffusion: Ah*Txx
!                (6)  = change due to meridional diffusion: Ah*Tyy
!                (7)  = change due to vertical diffusion:  kappa_h*Tzz
!                (8)  = change due to source term
!                (9)  = change due to explicit convection
!                (10) = change due to filtering
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of tracer "A" is -U(A)x = A(U)x - (UA)x
!                (11) = zonal advection U(Ax)
!                (12) = meridional advection V(Ay)
!                (13) = vertical advection W(Az)
!                (14) = change of tracer variance (tracer**2 units)
!                (15) = average tracer within volume (tracer units)
!     terr     = error term = (1) - sum (2) ... (10)
!     asst     = average sea surface tracer for regional surface areas
!     stflx    = average surface tracer flux for regional surface areas
!                tracer (#1,#2) units = (cal/cm**2/sec, gm/cm**2/sec)

!     termbm   = term balance components for time rate of change of
!                momentum within a volume. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (13) where each
!                term has the units of "cm/sec**2" and "Q" is the
!                momentum component {zonal or meridional} using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the momentum
!                (2)  = change due to the pressure gradient: grad_p
!                       without the surface pressure gradients
!                       (i.e., for computing the internal modes)
!                (3)  = change due to zonal nonlinear term: (UQ)x
!                (4)  = change due to meridional nonlinear term: (VQ)y
!                (5)  = change due to vertical nonlinear term: (wQ)z
!                (6)  = change due to zonal viscosity: Am*Qxx
!                (7)  = change due to meridional viscosity: Am*Qyy
!                (8)  = change due to vertical viscosity: kappa_m*Qzz
!                (9)  = change due to metric diffusion terms
!                (10) = change due to coriolis terms: fQ
!                (11) = change due to source terms
!                (12) = change due to surface pressure gradient
!                       this is obtained after solving the external mode
!                       in the stream function technique. It is solved
!                       directly from the elliptic equation for the
!                       prognostic surface pressure technique
!                (13) = change due to metric advection
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of vel component "Q" is -U(Q)x = Q(U)x - (UQ)x
!                (14) = zonal advection U(Qx)
!                (15) = meridional advection V(Qy)
!                (16) = vertical advection W(Qz)
!                (17) = average velocity component
!     smflx    = average surface momentum flux for regional surf areas
!                in dynes/cm**2
!     avgw     = average vertical velocity (cm/sec)

!     ustf     = names & units for surface tracer fluxes

      integer ntterms, nuterms
      parameter (ntterms=15, nuterms=17)

      character(15) :: ustf(nt,2)
      common /termb_c/ ustf

      real asst, avgw, termbt, termbm, smflx, stflx, terr
      common /termb_r/ asst(nt,0:nhreg), avgw(numreg)
      common /termb_r/ termbt(0:km,ntterms,nt,0:numreg)
      common /termb_r/ termbm(0:km,nuterms,2,numreg)
      common /termb_r/ smflx(2,0:nhreg), stflx(nt,0:nhreg), terr(nt)
# endif
#endif
