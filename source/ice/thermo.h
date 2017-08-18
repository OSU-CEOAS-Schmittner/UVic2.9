!====================== include file "thermo.h" ========================

! thermodynamic parameters for model in cgs units
! temperatures are in Celsius

! error tolerances
      real maxedif ! error tolerance for energy difference in thermo
      real maxwdif ! error tolerance for water difference in thermo
      real errmax  ! error tolerance for temp difference in tstm
      parameter (maxedif=5.e2*kilo)
      parameter (maxwdif=1.e-7*centi)
# if defined O_global_sums
      parameter (errmax=5.0e-8)
# else
      parameter (errmax=5.0e-6)
# endif

      real cpice,cpsno
      parameter (cpice=2054.e+04)  ! heat capacity of fresh ice
      parameter (cpsno=0.0)        ! heat capacity of snow

      common /handy/ rflice, rflsno, rslice, rslsno, rvlice, rvlsno
      common /handy/ rcpice, rcpsno, rcpatm, rvlatm, rslatm

      real rflice, rflsno, rslice, rslsno, rvlice, rvlsno
      real rcpice, rcpsno, rcpatm, rvlatm, rslatm

      real gamma ! param for heat capacity (J deg/kg)
      common /gam/ gamma

! for computing melting temp as a function of salinity, -alpha*salinity (deg)
      real alpha
      parameter(alpha=0.054)

      real tffresh ! freezing temperature of freshwater (K)
      real tmelt   ! melting temperature of ice top surface (C)
      real tsmelt  ! melting temperature of snow top surface (C)
      parameter (tffresh=273.16)
      parameter (tmelt=0.0)
      parameter (tsmelt=0.0)

! coefficients for computing saturation vapour pressure
      real aoc,boc,ai,bi
      parameter (aoc=17.2694)          ! over ocean
      parameter (boc=35.86-tffresh)    ! over ocean
      parameter (ai=21.8746)           ! over ice
      parameter (bi=7.66-tffresh)      ! over ice
      real qs1
!(mol weight of water:dry air)/(surface pressure in mb)*6.11
      parameter (qs1=0.622*6.11/1013.)
      real kappa ! solar extinction coef in ice (1/m)
      parameter (kappa=1.5/centi)

! conductivity of sea ice, k=kappai+beta*salinity/T
      real kappai ! thermal conductivity of fresh ice (W/m/deg)
      real kappas ! thermal conductivity of snow (W/m/deg)
      real kimin  ! minimum conductivity in ice (W/m/deg)
      real beta   ! param for conductivity (W/m)
      parameter (kappai = 2.0340*kilo*centi)
      parameter (kappas = 0.3100*kilo*centi)
      parameter (kimin  = 0.1000*kilo*centi)
      parameter (beta   = 0.1172*kilo*centi)

!  parameters for computing bottom melting rate lead from Steele, 199X
!  give a relaxation time constant of about 11 days
      real ut, ch
      parameter (ut=.01*centi, ch=0.0058)

!  parameters for computing lateral melting rate lead
!  from Maykut and Perovich, 199X
!  and floe perimeter parameters from Steele, 199X
      real  m1, m2, wclead, pi_eta
      parameter (m1=3.0e-4,m2=1.36, wclead=1.0e+3*centi)
      parameter (pi_eta=3.14159265358973/0.66)
