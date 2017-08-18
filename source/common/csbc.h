!======================= include file "csbc.h" =========================

!     surface boundary conditions (S.B.C.)

!     numsbc    = total number of surface boundary conditions
!     indices and names set in UVic_ESCM.F:
!       itaux     = x component of wind stress (dynes cm-2)
!       itauy     = y component of wind stress (dynes cm-2)
!       iws       = surface wind speed (cm s-1)
!       iaca      = atmospheric coalbedo
!       isca      = surface coalbedo (%/100)
!       ihflx     = heat flux
!       isflx     = salt flux
!       isst      = ocean model SST (C)
!       isss      = ocean model SSS (PSU*0.001-0.035)
!       iwa       = surface wind angle (degrees)
!       iro       = surface runoff (g cm-2 s-1)
!       iwxq      = x component of wind for moisture advection (cm s-1)
!       iwyq      = y component of wind for moisture advection (cm s-1)
!       iwxt      = x component of wind for temperature advection (cm s-1)
!       iwyt      = y component of wind for temperature advection (cm s-1)
!       iwxc      = x component of wind for carbon advection (cm s-1)
!       iwyc      = y component of wind for carbon advection (cm s-1)
!       ipsw      = penetrating shortwave (cal cm-2 s-1)
!       isu       = x component of ocean first layer velocity (cm s-1)
!       isv       = y component of ocean first layer velocity (cm s-1)
!       igu       = x component of ocean second layer velocity (cm s-1)
!       igv       = y component of ocean second layer velocity (cm s-1)
!       issdic    = sea surface concentration of dic (umol cm-3)
!       idicflx   = sea surface flux of carbon (umol cm-2 s-1)
!       issalk    = sea surface concentration of alkalinity (umol cm-3)
!       ialkflx   = sea surface flux of alkalinity (umol cm-2 s-1)
!       isso2     = sea surface concentration of oxygen (umol cm-3)
!       io2flx    = sea surface flux of oxygen (umol cm-2 s-1)
!       isspo4    = sea surface concentration of phosphate (nmol cm-3)
!       ipo4flx   = sea surface flux of phosphate (nmol cm-2 s-1)
!       issphyt   = sea surface concentration of phytoplankton (nmol P cm-3)
!       iphytflx  = sea surface flux of phytoplankton (nmol P cm-2 s-1)
!       isszoop   = sea surface concentration of zooplankton (nmol P cm-3)
!       izoopflx  = sea surface flux of zooplankton (nmol P cm-2 s-1)
!       issdetr   = sea surface concentration of detritus (nmol P cm-3)
!       idetrflx  = sea surface flux of detritus (nmol P cm-2 s-1)
!       issno3    = sea surface concentration of nitrate (nmol P cm-3)
!       ino3flx   = sea surface flux of nitrate (nmol P cm-2 s-1)
!       issdiaz   = sea surface concentration of diazotraphs (nmol P cm-3)
!       idiazflx  = sea surface flux of diazotraphs (nmol P cm-2 s-1)
!       issc14    = sea surface concentration of carbon 14 (umol cm-3)
!       ic14flx   = sea surface flux of carbon 14 (umol cm-2 s-1)
!       isscfc11  = sea surface concentration of cfc11 (umol cm-3)
!       icfc11flx = sea surface flux of cfc11 (umol cm-2 s-1)
!       isscfc12  = sea surface concentration of cfc12 (umol cm-3)
!       icfc12flx = sea surface flux of cfc12 (umol cm-2 s-1)
!       iat       = surface air temperature (C)
!       irh       = surface relative humidity (%/100)
!       ipr       = precipitation as rain (g cm-2 s-1)
!       ips       = precipitation as snow (g cm-2 s-1)
!       iaws      = averaged surface wind speed (cm s-1)
!       iswr      = surface shortwave radiation (g s-3)
!       ilwr      = surface longwave radiation (g s-3)
!       isens     = surface sensible heat flux (g s-3)
!       ievap     = surface evaporation (g cm-2 s-1)
!       idtr      = diurnal temperature range (C)
!       inpp      = carbon flux from net primary production (kg C m-2 s-1)
!       isr       = carbon flux from soil respiration (kg m-2 s-1)
!       iburn     = carbon flux from burning vegetation (kg m-2 s-1)
!       ibtemp    = ocean bottom temperature (C)
!       ibsalt    = ocean bottom temperature (PSU*0.001-0.035)
!       ibdic     = ocean bottom dissolved organic carbon (umol cm-3)
!       ibdicfx   = ocean bottom flux of carbon (umol cm-2 s-1)
!       ibalk     = ocean bottom alkalinity (umol cm-3)
!       ibalkfx   = ocean bottom flux of alkalinity (umol cm-2 s-1)
!       ibo2      = ocean bottom oxygen (umol cm-3)
!       ircal     = rain rate of calcite (umol cm-2 s-1)
!       irorg     = rain rate of organic carbon (umol cm-2 s-1)
!     mapsbc      = surface boundary conditions names
!     sbc         = surface boundary condition data.
!     gaost       = global average ocean surface tracer
!     socn        = average ocean sea surface salinity
!     ntspas      = the number of ocean time steps per atmosphere segment
!     ntspls      = the number of ocean time steps per land segment
!     ntspos      = the number of ocean time steps per ocean segment
!     nats        = number of atm time steps since mixing time step
!     nots        = number of ocn time steps since mixing time step
!     addflxa     = logical flag for adding only even mode fluxes from atm
!     addflxo     = logical flag for adding only even mode fluxes from ocn
!     dtatm       = atmosphere time step (s)
!     dtism       = icesheet time step (s)
!     dtismyr     = icesheet time step (years)
!     nismacc     = icesheet acceleration time (years)
!     dtlnd       = land time step (s)
!     dtocn       = ocean time step (s)
!     dtsed       = sediment time step (s)
!     dtsedyr     = sediment time step (years)
!     nsedacc     = sediment acceleration time (years)
!     atatm       = atmosphere boundary condition averaging time (s)
!     atlnd       = land boundary condition averaging time (s)
!     atocn       = ocean boundary condition averaging time (s)
!     atsed       = sediment boundary condition averaging time (s)
!     atism       = icesheet boundary condition averaging time (s)
!     dampts      = time scale for damping surface tracers (days) to data
!     dampdz      = ocean level thickness for converting Newtonian damping
!                   to a surface flux
!     land_map    = map with indices for coupling to land arrays
!     dtoih       = total system heat lost minus heat gained
!     dtoic       = total system carbon lost minus carbon gained

      integer numsbc
      parameter (numsbc = 14
#if defined O_embm_awind
     &                  + 1
#endif
#if defined O_carbon_co2_2d
     &                  + 2
#endif
#if defined O_shortwave
     &                  + 1
#endif
#if defined O_ice_evp
     &                  + 4
#endif
#if defined O_carbon
     &                  + 2
# if defined O_carbon_14
     &                  + 2
# endif
#endif
#if defined O_npzd_alk
     &                  + 2
#endif
#if defined O_npzd_o2
     &                  + 2
#endif
#if defined O_npzd
     &                  + 2
# if !defined O_npzd_no_vflux
     &                  + 6
# endif
#endif
#if defined O_npzd_nitrogen
     &                  + 2
# if !defined O_npzd_no_vflux
     &                  + 2
# endif
#endif
#if defined O_cfcs_data_transient
     &                  + 4
#endif
#if defined O_mtlm
     &                  + 10
#endif
#if defined O_mtlm && defined O_carbon
     &                  + 3
#endif
#if defined O_sed
     &                  + 7

# if defined O_carbon
     &                  + 1
# endif
# if defined O_npzd_alk
     &                  + 1
# endif
#endif
     &                     )

      integer itaux, itauy, iws, iaca, isca, ihflx, isflx, isst, isss
      integer iwa, iro, iwxq, iwyq, iwxt, iwyt, iwxc, iwyc, ipsw, isu
      integer isv, igu, igv, issdic, idicflx, issalk, ialkflx, isso2
      integer io2flx, isspo4, ipo4flx, issphyt, iphytflx, isszoop
      integer izoopflx, issdetr, idetrflx, issno3, ino3flx, issdiaz
      integer idiazflx, issc14, ic14flx, isscfc11, icfc11flx, isscfc12
      integer icfc12flx, iat, irh, ipr, ips, iaws, iswr, ilwr, isens
      integer ievap, idtr, inpp, isr, iburn, ibtemp, ibsalt, ibdic
      integer ibdicfx, ibalk, ibalkfx, ibo2, ircal, irorg

      common /csbc_i/ itaux, itauy, iws, iaca, isca, ihflx, isflx, isst
      common /csbc_i/ isss, iwa, iro, iwxq, iwyq, iwxt, iwyt, iwxc, iwyc
      common /csbc_i/ ipsw, isu, isv, igu, igv, issdic, idicflx, issalk
      common /csbc_i/ ialkflx, isso2, io2flx, isspo4, ipo4flx, issphyt
      common /csbc_i/ iphytflx, isszoop, izoopflx, issdetr, idetrflx
      common /csbc_i/ issno3, ino3flx, issdiaz, idiazflx, issc14
      common /csbc_i/ ic14flx, isscfc11, icfc11flx, isscfc12, icfc12flx
      common /csbc_i/ iat, irh, ipr, ips, iaws, iswr, ilwr, isens, ievap
      common /csbc_i/ idtr, inpp, isr, iburn, ibtemp, ibsalt, ibdic
      common /csbc_i/ ibdicfx, ibalk, ibalkfx, ibo2, ircal, irorg

      character(20) :: mapsbc
      common /csbc_c/ mapsbc(numsbc)

      real sbc
      common /csbc_r/ sbc(imt,jmt,numsbc)

      real gaost, socn

      common /csbc_r/ gaost(nt), socn

      integer ntspas, ntspls, ntspos, nats, nots, nismacc, nsedacc
      common /csbc_i/ ntspas, ntspls, ntspos, nats, nots, nismacc
      common /csbc_i/ nsedacc

      logical addflxa, addflxo
      common /csbc_l/ addflxa, addflxo

      real dtatm, dtism, dtismyr, dtlnd, dtocn, dtsed, dtsedyr, atatm
      real atism, atlnd, atocn, atsed
      common /csbc_r/ dtatm, dtism, dtismyr, dtlnd, dtocn, dtsed
      common /csbc_r/ dtsedyr, atatm, atism, atlnd, atocn, atsed

      real dampts, dampdz, subflux, subz
      common /csbc_r/ dampts(nt), dampdz(nt)
#if defined O_plume
      common /csbc_r/ subflux(imt,jmt,nt), subz(imt,jmt)
#endif
#if defined O_mtlm
      integer land_map

      common /csbc_i/ land_map(imt,jmt)
#endif
#if defined O_global_sums || defined O_co2emit_diag
      real dtoih, dtoic

      common /csbc_r/ dtoih, dtoic
#endif
