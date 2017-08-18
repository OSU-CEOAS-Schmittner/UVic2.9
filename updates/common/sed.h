!     sed include file

!     ipmax  = maximum number of sediment points
!     nzmax  = maximum number of mixed layers
!     ibmax  = maximum number of buried layers

      integer ipmax, nzmax, ibmax, imap, jmap, map_sed, ipsed, kmax
      integer kmin

      parameter (ipmax=imt*jmt, nzmax=8, ibmax=20)

!     imap    = 1d maping index for 2d map (i index)
!     jmap    = 1d maping index for 2d map (j index)
!     map_sed = 2d map of sediment indices
!     ipsed   = number of sediment points (must be >= to ipmax)
!     kmax    = number of mixed layers (must be less than nzmax)
!     kmin    = minimum kmt level for sediments

      common /sed_i/ imap(ipmax), jmap(ipmax), map_sed(imt,jmt)
      common /sed_i/ ipsed, kmax, kmin

!   mixed layer fields
!     carb           = forms of carbon, 1=co2,2=hco3,3=co3 (mol l-1)
!     dcpls          = diffusion coefficient above, 1=co2,2=hco3,3=co3
!     dcmin          = diffusion coefficient below, 1=co2,2=hco3,3=co3
!     pore           = porosity
!     form           = pore**3.
!     o2             = bottom water oxygen (mol l-1)
!     orggg          = organic carbon fraction
!     orgml          = organic carbon mass (g cm-2)
!     calgg          = calcite fraction
!     calml          = calcite mass (g CaCO3 cm-2)
!     dopls          = diffusion coefficient for oxygen above
!     domin          = diffusion coefficient for oxygen below
!     dbpls          = diffusion coefficient for organic carbon above
!     dbmin          = diffusion coefficient for organic carbon below
!   buried layer fields
!     buried_mass    = buried mass (g CaCO3 cm-2)
!     buried_calfrac = buried calcite fraction
!     depth_age      = age of buried level (years)
!   sediment surface or integrated fields
!     zrct           = maximum respiration depth (cm)
!     water_z_p      = ocean depth (m)
!     k1             = first dissociation constant for carbonic acid
!     k2             = second dissociation constant for carbonic acid
!     k3             = dissociation constant for boric acid
!     csat           = carbonate saturation at sediment surface (mol l-1)
!     rc             = respiration coefficient
!     sed_ml_mass    = mixed layer mass (g cm-2)
!     ttrorg         = dissolution rate organic carbon (mol cm-2 dtsed-1)
!     ttrcal         = dissolution rate of calcite (mol cm-2 dtsed-1)
!     c_advect       = burial flux [mol cm-2 dtsed-1]
!     zsed           = mixed layer grid depth (cm)
!     delz           = mixed layer grid layer thickness (cm)
!     rain_org_p     = rain rate of organic carbon (mol cm-2 dtsed-1)
!     rain_cal_p     = rain rate of calcite (mol cm-2 dtsed-1)
!     co3_p          = carbonate at sediment surface (mol l-1)
!   scalars
!     dissc          = calcite dissolution constant
!     dissn          = calcite dissolution constant
!     weath          = weathering flux (kg s-1)
!     weathflx       = weathering flux (umol s-1)
!     sed_year       = year for sediment profile
!     sedsa          = surface area of potential deep sediments (cm2)
!     carblith       = change in lithosphere carbon (umol)
!                      set to zero if time is initialized

      real carb, dcpls, dcmin, pore, form, o2, orggg, orgml, calgg
      real calml, dopls, domin, dbpls, dbmin, buried_mass
      real buried_calfrac, depth_age, zrct, water_z_p, k1, k2, k3, csat
      real rc, sed_ml_mass, ttrorg, ttrcal, c_advect, zsed, delz
      real rain_org_p, rain_cal_p, co3_p, dissc, dissn, weath, weathflx
      real dicwflx, sed_year, sedsa, carblith

      common /sed_r/ carb(nzmax,3,ipmax), dcpls(nzmax,3,ipmax)
      common /sed_r/ dcmin(nzmax,3,ipmax)
      common /sed_r/ pore(nzmax,ipmax), form(nzmax,ipmax)
      common /sed_r/ o2(nzmax,ipmax), orggg(nzmax,ipmax)
      common /sed_r/ orgml(nzmax,ipmax), calgg(nzmax,ipmax)
      common /sed_r/ calml(nzmax,ipmax), dopls(nzmax,ipmax)
      common /sed_r/ domin(nzmax,ipmax), dbpls(nzmax,ipmax)
      common /sed_r/ dbmin(nzmax,ipmax), buried_mass(ibmax,ipmax)
      common /sed_r/ buried_calfrac(ibmax,ipmax), depth_age(ibmax)
      common /sed_r/ zrct(ipmax), water_z_p(ipmax), k1(ipmax), k2(ipmax)
      common /sed_r/ k3(ipmax), csat(ipmax), rc(ipmax), ttrorg(ipmax)
      common /sed_r/ ttrcal(ipmax), sed_ml_mass(ipmax), c_advect(ipmax)
      common /sed_r/ zsed(nzmax), delz(nzmax), rain_org_p(ipmax)
      common /sed_r/ rain_cal_p(ipmax), co3_p(ipmax)
      common /sed_r/ dissc, dissn, weath, weathflx, dicwflx, sed_year
      common /sed_r/ sedsa, carblith

!     ntatss         = time step counter for time averaging
!     ta_ttrcal      = time average dissolution rate of calcite
!     ta_rain_cal    = time average rain rate of calcite
!     ta_cal         = time average mixed layer calcite fraction
!     ta_calmass     = time average mixed layer calcite mass
!     ta_calmass_bur = time average buried calcite mass
!     ta_co3         = time average carbonate at sediment surface
!     ta_co3sat      = time average carbonate saturation at surface
!     ta_rainr       = time average rain ratio

      integer ntatss

      common /sed_i/ ntatss

      real ta_ttrcal, ta_rain_cal, ta_cal, ta_calmass, ta_calmass_bur
      real ta_co3, ta_co3sat, ta_rainr

      common /sed_r/ ta_ttrcal(imt,jmt), ta_rain_cal(imt,jmt)
      common /sed_r/ ta_cal(imt,jmt), ta_calmass(imt,jmt)
      common /sed_r/ ta_calmass_bur(imt,jmt), ta_co3(imt,jmt)
      common /sed_r/ ta_co3sat(imt,jmt), ta_rainr(imt,jmt)

!     ntatis          = number of time averaged time step integrals
!     tai_ttrcal      = time step integral of ttrcal
!     tai_rain_cal    = time step integral of rain_cal
!     tai_cal         = time step integral of cal
!     tai_weathflx    = time step integral of weathflx
!     tai_calmass     = time step integral of mixed layer calcite mass
!     tai_calmass_bur = time step integral of buried calcite mass
!     tai_co3         = time step integral of co3
!     tai_co3sat      = time step integral of co3sat
!     tai_rainr       = time step integral of rain ratio
!     tai_carblith    = time step integral of carblith
!     tai_cfo2s       = average total flux ocean to sediments
!     tai_cfl2o       = average total flux land to ocean (weathering)

      integer ntatis

      common /sed_i/ ntatis

      real tai_ttrcal, tai_rain_cal, tai_cal, tai_weathflx, tai_calmass
      real tai_calmass_bur, tai_co3, tai_co3sat, tai_rainr, tai_carblith
      real tai_cfo2s, tai_cfl2o

      common /sed_r/ tai_ttrcal, tai_rain_cal, tai_cal, tai_weathflx
      common /sed_r/ tai_calmass, tai_calmass_bur, tai_co3, tai_co3sat
      common /sed_r/ tai_rainr, tai_carblith, tai_cfo2s, tai_cfl2o
