      subroutine co2calc_SWS (t, s, dic_in, ta_in, xco2_in, pt_in
     &,                       sit_in, atmpres, depth, phlo, phhi, ph
     &,                       co2star, dco2star, pco2, hSWS
     &,                       CO3, Omega_c, Omega_a)
!-------------------------------------------------------------------------
! Modified from co2calc.f (RCS version 1.8, OCMIP-2)
! Constants are given on "seawater" H scale (hSWS) except for the "free" 
! H scale dissociation constants for S and F (necessary to work with hSWS).

! PURPOSE: Calculate delta co2* (dco2star)

! INPUT:
!        t       = temperature (degrees C)
!        s       = salinity (PSU)
!        dic_in  = total inorganic carbon (mol/m^3)
!        ta_in   = total alkalinity (eq/m^3)
!        xco2_in = atmospheric CO2 concentration (ppmv)
!        pt_in   = inorganic phosphate (mol/m^3)
!        sit_in  = inorganic silicate (mol/m^3)
!        atmpres = atmospheric pressure in atmospheres (1 atm==1013.25mbar)
!        depth   = ocean depth (m)
!        phlo    = lower limit of pH range
!        phhi    = upper limit of pH range

! OUTPUT:
!        co2star  = ocean CO2* (mol/m^3)
!        dco2star = delta (atm-ocn) CO2* (mol/m^3)
!        pco2     = oceanic pCO2 (uatm)
!        hsws     = H+ on seawater H scale
!--------------------------------------------------------------------------

      implicit none

      real permil, pt, pt_in, sit, sit_in, ta, ta_in, dic, dic_in
      real permeg, xco2, xco2_in, tk, t, tk100, tk1002, dlogtk, s
      real sqrtis, s2, sqrts, s15, scl, bt, st, ft, ff, x1, phhi
      real x2, phlo, xacc, hSWS, drtsafe, hSWS2, co2star, co2starair
      real atmpres, dco2star, ph, pco2, invtk, is, is2, k0, k1
      real k12, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, C2K
      real Kspc, Kspa, ca, co3, omega_c, omega_a, depth, pres
      real t2, dvc, dva, dk, pitkR, p2itkR, Rgas, fugf

!--------------------------------------------------------------------------
! Models carry tracers in mol/m^3 but this routine uses umol/kg
! Change units from the input of mol/m^3 -> mol/kg:
! (1 mol/m^3) x (1 m^3/1024.5 kg)
! where the ocean's mean surface density is 1024.5 kg/m^3

      permil = 1./1024.5    ! to convert mol/m^3 to mol/kg
      permeg = 1.e-6        ! to convert uatm to atm.
      Rgas = 83.1451        ! ideal gas constant R in cm3 bar k-1
      C2K = 273.15          ! to cinvert C to K

!     Convert input
      dic = dic_in*permil   ! change mol/m^3 to mol/kg
      ta = ta_in*permil     ! change mol/m^3 to mol/kg
      xco2 = xco2_in*permeg ! change units from ppmv to parts by volume.
      pt = pt_in*permil     ! change mol/m^3 to mol/kg
      sit = sit_in*permil   ! change mol/m^3 to mol/kg
      pres = depth*0.1      ! change from meters to decibars

!--------------------------------------------------------------------------
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Version 2 of "Handbook of Methods for the Analysis of the Various
! Parameters of the Carbon Dioxide System in Seawater", DOE, 1994
! (SOP No. 3, p25-26).

! Derive simple terms used more than once
      tk = C2K + t
      tk100 = tk/100.
      tk1002 = tk100*tk100
      invtk = 1./tk
      dlogtk = log(tk)
      is = 19.924*s/(1000. - 1.005*s)
      is2 = is*is
      sqrtis = sqrt(is)
      s2 = s*s
      t2 = t*t
      sqrts = sqrt(s)
      s15 = s**1.5
      scl = s/1.80655
      pitkR = pres/tk/Rgas
      p2itkR = pres*pitkR

!------------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride
! Uppstrom (1974)
      bt = 0.000232*scl/10.811
! Morris & Riley (1966)
      st = 0.14*scl/96.062
! Riley (1965)
      ft = 0.000067*scl/18.9984

!------------------------------------------------------------------------
! f = k0(1-pH2O)*correction term for non-ideality
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
      ff = exp(-162.8301 + 218.2968/tk100  + 90.9241*log(tk100)
     &   - 1.47696*tk1002 + s*(.025695 - .025225*tk100
     &   + 0.0049867*tk1002))

!------------------------------------------------------------------------
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3] on hSWS
! Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale
! (Original reference: Dickson and Millero, DSR, 1987)
      k1 = 10**(-1.*(3670.7*invtk - 62.008 + 9.7944*dlogtk
     &   - 0.0118*s + 0.000116*s2))
      k2 = 10**(-1*(1394.7*invtk + 4.777 - 0.0184*s + 0.000118*s2))

!------------------------------------------------------------------------
! k1p = [H][H2PO4]/[H3PO4] on hSWS
! Millero p.670 (1995)
      k1p = exp(-4576.752*invtk + 115.540 - 18.453*dlogtk
     &    + (-106.736*invtk + 0.69171)*sqrts
     &    + (-0.65643*invtk - 0.01844)*s)

!------------------------------------------------------------------------
! k2p = [H][HPO4]/[H2PO4] on hSWS
! Millero p.670 (1995)
      k2p = exp(-8814.715*invtk + 172.1033 - 27.927*dlogtk
     &    + (-160.340*invtk + 1.3566)*sqrts
     &    + (0.37335*invtk - 0.05778)*s)

!------------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4] on hSWS
! Millero p.670 (1995)
      k3p = exp(-3070.75*invtk - 18.126
     &    + (17.27039*invtk + 2.81197)*sqrts
     &    + (-44.99486*invtk - 0.09984)*s)

!------------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
! Millero p.671 (1995) using data from Yao and Millero (1995)
! change to (mol/ kg soln)
! depth dependancy assumed to be the same as boric acid
! typo in Millero 1994 corrected in sign of 0.1622

       ksi = exp(-8904.2*invtk + 117.400 - 19.334*dlogtk
     &     + (-458.79*invtk + 3.5913)*sqrtis
     &     + (188.74*invtk - 1.5998)*is
     &     + (-12.1652*invtk + 0.07871)*is2
     &     + log(1. - 0.001005*s))

!------------------------------------------------------------------------
! kw = [H][OH] on hSWS
! Millero p.670 (1995) using composite data
! pressure dependancy in Millero 1994 corrected for sea water from
! Millero 1983
       kw = exp(-13847.26*invtk + 148.9802 - 23.6521*dlogtk
     &    + (118.67*invtk - 5.977 + 1.0495*dlogtk)*sqrts - 0.01615*s)

!------------------------------------------------------------------------
! ks = [H][SO4]/[HSO4] on free H scale
! Dickson (1990, J. chem. Thermodynamics 22, 113)
! change to (mol/kg soln)
      ks = exp(-4276.1*invtk + 141.328 - 23.093*dlogtk
     &   + (-13856*invtk + 324.57 - 47.986*dlogtk)*sqrtis
     &   + (35474*invtk - 771.54 + 114.723*dlogtk)*is
     &   - 2698*invtk*is**1.5 + 1776*invtk*is2
     &   + log(1. - 0.001005*s))

!------------------------------------------------------------------------
! kf = [H][F]/[HF] on free H scale
! Dickson and Riley (1979)
! change to (mol/ kg soln)
      kf = exp(1590.2*invtk - 12.641 + 1.525*sqrtis
     &   + log(1. - 0.001005*s))

!------------------------------------------------------------------------
! kb = [H][BO2]/[HBO2] on hSWS
! Dickson p.673 (1990)
! change from htotal to hSWS
! typo in Millero 1994 corrected in sign of 0.1622
      kb = exp((-8966.90 - 2890.53*sqrts - 77.942*s
     &   + 1.728*s15 - 0.0996*s2)*invtk
     &   + (148.0248 + 137.1942*sqrts + 1.62142*s)
     &   + (-24.4344 - 25.085*sqrts - 0.2474*s)*dlogtk
     &   + 0.053105*sqrts*tk
     &   + log((1+(st/ks)+(ft/kf))/(1+(st/ks))))

# if defined O_save_carbon_carbonate_depth
!------------------------------------------------------------------------
! apply pressure correction terms for calculations at depth

      k1 = k1*exp((25.5 - 0.1271*t)*pitkR
     &   + 0.5*(-3.08e-3 + 8.77e-5*t)*p2itkR)

      k2 = k2*exp((15.82 + 0.0219*t)*pitkR
     &   + 0.5*(1.13e-3 - 1.475e-4*t)*p2itkR)

      k1p = k1p*exp((14.51 - 0.1211*t + 3.21e-4*t2)*pitkR
     &    + 0.5*(-2.67e-3 + 4.27e-5*t)*p2itkR)

      k2p = k2p*exp((23.12 - 0.1758*t + 2.647e-3*t2)*pitkR
     &    + 0.5*(-5.15e-3 + 9.e-5*t)*p2itkR)

      k3p = k3p*exp((26.57 - 0.202*t + 3.042e-3*t2)*pitkR
     &    + 0.5*(-4.08e-3 + 7.14e-5*t)*p2itkR)

      ksi = ksi*exp((29.48 - 0.1622*t + 2.608e-3*t2)*pitkR
     &    + 0.5*(-2.84e-3)*p2itkR)

      kw = kw*exp((20.02 - 0.1119*t + 1.409e-3*t2)*pitkR
     &   + 0.5*(-5.13e-3 + 7.94e-5*t)*p2itkR)

      ks = ks*exp((18.03 - .0466*t - 3.16e-4*t2)*pitkR
     &   + 0.5*(-4.53e-3 + 9.e-5*t)*p2itkR)

      kf = kf*exp((9.78 + 9.e-3*t + 9.42e-4*t2)*pitkR
     &   + 0.5*(-3.91e-3 + 5.4e-5*t)*p2itkR)

      kb = kb*exp((29.48 - 0.1622*t + 2.608e-3*t2)*pitkR
     &   + 0.5*(-2.84e-3)*p2itkR)
# endif

!-------------------------------------------------------------------------
! Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
# if defined O_carbon_carbonate_approx
! Use a non-iterative solver

      call calc_hSWS_approx (dic, ta, pt, sit, bt, k1, k2
     &,                      k1p, k2p, k3p, kb, kw, ksi, hSWS)

# else
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.

! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate hSWS. In this case we use the Newton-Raphson
! "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
! error trapping removed).

! As currently set, this procedure iterates about 12 times. The x1 and x2
! values set below will accomodate ANY oceanographic values. If an initial
! guess of the pH is known, then the number of iterations can be reduced to
! about 5 by narrowing the gap between x1 and x2. It is recommended that
! the first few time steps be run with x1 and x2 set as below. After that,
! set x1 and x2 to the previous value of the pH +/- ~0.5. The current
! setting of xacc will result in co2star accurate to 3 significant figures
! (xx.y). Making xacc bigger will result in faster convergence also, but this
! is not recommended (xacc of 10**-9 drops precision to 2 significant figures).

      x1 = 10.**(-phhi)
      x2 = 10.**(-phlo)
      xacc = 1.e-10
      hSWS = drtsafe (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st
     &,               ft, sit, pt, dic, ta, x1, x2, xacc)

#endif
!-------------------------------------------------------------------------
! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
      hSWS2 = hSWS*hSWS
      co2star = (dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2))
      co2starair = xco2*ff*atmpres
      dco2star = co2starair - co2star
# if defined O_save_carbon_carbonate_chem
      ph = -log10(hSWS)
! Calculate [CO3] in mol/kg from co2star and hSWS
      CO3 = k1*k2*co2star/hSWS2
! K0 and fugf from Weiss 1974 Marine Chemistry
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585*log(tk100)
     &   + s*(.023517 - 0.023656*tk100 + 0.0047036*tk1002))
! Fugacity factor (assumes pressure of 1 atm or 1.01325 bars)
      fugf = exp((-1636.75 + 12.0408*tk - 0.0327957*tk**2 
     &     + 3.16528e-5*tk**3 + 2.*(57.7 - 0.118*tk))*1.01325/(Rgas*tk))
! Calculate the partial pressures
      pCO2 = co2star/(k0*fugf)/permeg ! convert atm to uatm

!-------------------------------------------------------------------------
! Solubility products of calcite and aragonite at sea level
! CalciteSolubility:
! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
      Kspc = 10**(-171.9065 - 0.077993*tk + 2839.319/tk
     &       + 71.595*dlogtk/log(10.) + (-0.77712 + 0.0028426*tk 
     &       + 178.34/tk)*sqrts - 0.07711*s + 0.0041249*sqrts*s)

! Aragonite Solubility
! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983
      Kspa = 10**(-171.945 - 0.077993*tk + 2903.293/tk 
     &       + 71.595*dlogtk/log(10.) 
     &       + (-0.068393 + 0.0017276*tk + 88.135/tk)*sqrts
     &       - 0.10018*s + 0.0059415*sqrts*s)
#  if defined O_save_carbon_carbonate_depth

! Pressure Correction For Calcite:
! Ingle, Marine Chemistry 3:301-319, 1975
! Same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
! has typos (-.5304, -.3692, and 10^3 for Kappa factor)
      DVc = -48.76 + 0.5304*t
      DK  = (-11.76 + 0.3692*t)/1000

! Pressure Correction For Aragonite:
! Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
! same as Millero, GCA 1995 except for typos (-.5304, -.3692,
! and 10^3 for Kappa factor)
      DVa = DVc + 2.8

! KspP for calcite, pressure (pres) in bars
      Kspc = Kspc*exp(-DVc*pitkR + 0.5*DK*p2itkR)

! KspP for aragonite, pressure (pres) in bars
      Kspa = Kspa*exp(-DVa*pitkR + 0.5*DK*p2itkR)
#  endif

! Calculate saturation state (Omega)
! Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967
      Ca = (0.02128/40.087)*(s/1.80655)

! calcite (Omgc)
      Omega_c = Ca*CO3/Kspc

! aragonite (Omga)
      Omega_a = Ca*CO3/Kspa

      CO3 = CO3/permil ! convert from mol/kg -> mol/m^3
# else
! not used unless option save_carbon_carbonate_chem is set
      ph = 0.
      pCO2 = 0.
      CO3 = 0.
      Omega_c = 0.
      Omega_a = 0.
# endif

!-------------------------------------------------------------------------
      co2star = co2star/permil ! convert from mol/kg -> mol/m^3
      dco2star = dco2star/permil ! convert from mol/kg -> mol/m^3
 
      return
      end


      real function drtsafe (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                      bt, st, ft, sit, pt, dic, ta, x1, x2, xacc)
!-------------------------------------------------------------------------
!     File taken from Numerical Recipes. Modified  R.M.Key 4/94

      implicit none

      integer maxit, j

      real k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st, ft, sit
      real pt, dic, ta
      real x1, fl, df, x2, fh, xl, xh, swap, dxold, dx, f, temp, xacc

      maxit = 100
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                   bt, st, ft, sit, pt, dic, ta, x1, fl, df)
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                   bt, st, ft, sit, pt, dic, ta, x2, fh, df)
      if (fl .lt. 0.) then
        xl = x1
        xh = x2
      else
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap
      endif
      drtsafe = 0.5*(x1 + x2)
      dxold = abs(x2 - x1)
      dx = dxold
      call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                 bt, st, ft, sit, pt, dic, ta, drtsafe, f, df)
      do j=1,maxit
        if (((drtsafe - xh)*df - f)*((drtsafe - xl)*df - f) .ge. 0.
     &    .or. abs(2.*f) .gt. abs(dxold*df)) then
          dxold = dx
          dx = 0.5*(xh - xl)
          drtsafe = xl + dx
          if (xl .eq. drtsafe) return
        else
          dxold = dx
          dx = f/df
          temp = drtsafe
          drtsafe = drtsafe - dx
          if (temp .eq. drtsafe) return
        endif
        if (abs(dx) .lt. xacc) return
        call ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                   bt, st, ft, sit, pt, dic, ta, drtsafe, f, df)
        if (f .lt. 0.) then
          xl = drtsafe
          fl = f
        else
          xh = drtsafe
          fh = f
        endif
      enddo

      return
      end


      subroutine ta_iter_SWS (k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi
     &,                       bt, st, ft, sit, pt, dic, ta, x, fn, df)
!-------------------------------------------------------------------------
! This routine expresses TA as a function of DIC, hSWS and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS
! Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)

      implicit none

      real k1, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, bt, st, ft, sit
      real pt, dic, ta, x, fn, df
      real rx, x2, x3, k12, k12p, k123p, c, rc, a, ra, da, b, rb, db
      real rkb, rksi, kwrx, ptra, xrc, crx, crx2, t1, t2, t3, t4, t5, t6

      rx = 1./x
      x2 = x*x
      x3 = x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1. + st/ks + ft/kf
      rc = 1./c
      a = x3 + k1p*x2 + k12p*x + k123p
      ra = 1./a
      da = 3.*x2 + 2.*k1p*x + k12p
      b = x2 + k1*x + k12
      rb = 1./b
      db = 2.*x + k1
      rkb = 1./kb
      rksi = 1./ksi
      kwrx = kw*rx
      ptra = pt*ra
      xrc = x*rc
      crx = c*rx
      crx2 = crx*rx
      t1 = k1*dic*rb
      t2 = 2.*dic*k12*rb
      t3 = 1. + x*rkb
      t4 = 1. + x*rksi 
      t5 = 1. + ks*crx
      t6 = 1. + kf*crx

!-------------------------------------------------------------------------
! fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
! df = dfn/dx

      fn = t1*x + t2 + bt/t3 + kwrx + ptra*k12p*x + 2.*ptra*k123p 
     &   + sit/t4 - xrc - st/t5 - ft/t6 - ptra*x3 - ta

      df = t1 - t1*x*db*rb - t2*db*rb - bt*rkb/t3**2 - kwrx*rx 
     &   + (ptra*k12p*(a - x*da))*ra - 2.*ptra*k123p*da*ra
     &   - sit*rksi/t4**2 - rc - st*t5**(-2)*(ks*crx2) 
     &   - ft*t6**(-2)*(kf*crx2) - ptra*x2*(3.*a - x*da)*ra

      return
      end


      subroutine calc_hSWS_approx (dic, ta, pt, sit, bt, k1, k2, 
     &                             k1p, k2p, k3p, kb, kw, ksi, H)
!-------------------------------------------------------------------------
! solve carbonate system for H+
! based on code from M. Follows, T. Ito, S. Dutkiewicz
! see: Ocean Modelling 12 (2006), 290-301
!-------------------------------------------------------------------------

      implicit none

      real dic, ta, pt, sit, bt, k1, k2, k1p, k2p, k3p, kb, kw, ksi
      real H, gamm, hg, cag, bohg, h3po4g, h2po4g, hpo4g, po4g
      real siooh3g, recip, dummy

! dic           = dissolved inorganic carbon
! ta            = total alkalinity
! pt            = dissolved inorganic phosphorus 
! sit           = dissolved inorganic silica
! bt            = dissolved inorganic boron
! ta            = total alkalinity
! k1, k2        = carbonate equilibrium coefficients
! kw            = dissociation of water
! klp, k2p, k3p = phosphate equilibrium coefficients
! ksi, kb       = silicate and borate equilibrium coefficients
! H             = [H+]

! First guess of [H+]: from last timestep
      hg = H
! estimate contributions to total alk from borate, silicate, phosphate
      bohg = bt*kb/(hg + kb)
      siooh3g = sit*ksi/(ksi + hg)
      recip = 1./(hg*hg*hg + k1p*hg*hg + k1p*k2p*hg + k1p*k2p*k3p)
      h3po4g = (pt*hg*hg*hg)*recip
      h2po4g = (pt*k1p*hg*hg)*recip
      hpo4g = (pt*k1p*k2p*hg)*recip
      po4g = (pt*k1p*k2p*k3p)*recip
! estimate carbonate alkalinity
      cag = ta -bohg -(kw/hg) +hg -hpo4g -(2.*po4g) +h3po4g -siooh3g
! improved estimate of hydrogen ion concentration
      gamm = dic/cag
      dummy = (1. - gamm)*(1. - gamm)*k1*k1 - 4.*k1*k2*(1. - 2.*gamm)
      H = 0.5*((gamm - 1.)*k1 + sqrt(dummy))

      return
      end



