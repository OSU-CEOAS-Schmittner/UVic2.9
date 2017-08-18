      subroutine npzd_src (bioin, ntsb, tsb, gl, bct, impo, dzt
     &,                    dayfrac, wwd, rkw, nud, bioout, expoout
     &,                    grazout, morpout, morzout
# if defined O_save_npzd
     &,                    nppout, morptout, remiout, excrout
#  if defined O_npzd_nitrogen
     &,                    npp_Dout, graz_Dout, morp_Dout, nfixout
#  endif
# endif
     &                     )

#if defined O_mom && defined O_npzd
!=======================================================================
!     computes source terms of the NPZD model
!     initial version of code adapted from Xavier Giraud:
!     Giraud et al. 2000, J Mar Res, 58, 609-630
!     original model reference:
!     Oeschlies and Garcon 1999, Global Biogeochem. Cycles 13, 135-160
!     Schmittner et al. 2005,  Global Biogeochem. Cycles 19, GB3004,
!     doi:10.1029/2004GB002283.
!     This is the optimized model version of Schartau & Oschlies, 2003,
!     J. Mar. Res. 62, 765-793 and Oschlies & Schartau, in prep.
!     as described in Schmittner and Oschlies, in prep.
!     note that nutrient now represents phosphate

!     input variables:

!       bioin(1:4) = N,P,Z,D [mmol m-3]
!       bioin(5)   = nitrate [mmol m-3]
!       bioin(6)   = diazotrophs [mmol m-3]

!       gl         = 2.*light at top of grid box
!       ntsb       = number of time steps
!       tsb        = time step [s]
!       bct        = bbio**(cbio*temperature)
!       impo       = import of detritus from above [mmol m-3]
!       dzt        = depth of grid box [cm]
!       dayfrac    = day length (fraction: 0 < dayfrac < 1)
!       wwd        = sinking speed of detritus/dzt
!       rkw        = reciprical of kw*dzt(k)
!       nud        = remineralisation rate of detritus [s-1]

!     output variables:

!       bioout     = change from bioin [mmol m-3]
!       nppout     = net primary production [mmol m-3]
!       grazout    = grazing [mmol m-3]
!       morpout    = quadratic mortality of phytoplankton [mmol m-3]
!       morptout   = specific mortality of phytoplankton [mmol m-3]
!       morzout    = mortality of zooplankton [mmol m-3]
!       remiout    = remineralisation [mmol m-3]
!       excrout    = excretion [mmol m-3]
!       expoout    = detrital export [mmol m-3]
!       npp_Dout   = NPP of diazotrophs
!       graz_Dout  = grazing of diazotrophs
!       morp_Dout  = mortality of diazotrophs
!       nfixout    = rate of N2 fixation
!=======================================================================

      implicit none

      integer n, ntsb

      real gl, f1, bion, biop, bioz, biod, jmax, u_P, g_P, npp, graz
      real morp, morpt, morz, remi, excr, expo, impo, nppout, grazout
      real morpout, morptout, morzout, remiout, excrout, expoout, tsb
      real dzt, nflag, pflag, zflag, dflag, wwd, rkw, gd, dayfrac, bct
      real nupt, nud, biono3, u_D,npp_D, npp_Dout, no3flag, biodiaz
      real diazflag, g_D,graz_D, morp_D, jmax_D, gd_D, avej_D, no3upt_D
      real morp_Dout, graz_Dout, nfixout, biop2, u1, u2, phi1, phi2
      real avej

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "npzd.h"

      real bioin(ntnpzd), bioout(ntnpzd)

!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      f1 = exp((-kw - kc*bioin(2))*dzt)
      jmax = abio*bct
      gd = jmax*dayfrac
      u1 = max(gl/gd,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej = gd*rkw*(phi1 - phi2)
# if defined O_npzd_nitrogen
      jmax_D = max(0.,abio*(bct - 2.6))*jdiar
      gd_D = max(1.e-14,jmax_D*dayfrac)
      u1 = max(gl/gd_D,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*rkw*(phi1 - phi2)
# endif
      nupt = nupt0*bct

      bioout(:) = 0.0
      bion = bioin(1)
      biop = bioin(2)
      bioz = bioin(3)
      biod = bioin(4)
# if defined O_npzd_nitrogen
      biono3 = bioin(5)
      biodiaz = bioin(6)
# endif
      expoout = 0.0
      grazout = 0.0
      morpout = 0.0
      morzout = 0.0
# if defined O_save_npzd
      nppout = 0.0
      morptout = 0.0
      remiout = 0.0
      excrout = 0.0
#  if defined O_npzd_nitrogen
      npp_Dout = 0.0
      graz_Dout = 0.0
      morp_Dout = 0.0
      nfixout = 0.0
#  endif
# endif

      do n=1,ntsb

!       growth rate of phytoplankton
        u_P = min(avej, jmax*bion/(k1p + bion))
# if defined O_npzd_nitrogen
!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        u_D = min(avej_D, jmax_D*bion/(k1p + bion))
# endif
        npp = u_P*biop
# if defined O_npzd_nitrogen
        npp_D = max(0.,u_D*biodiaz)
!       grazing function for diazotrophs
        g_D = gbio*epsbio*biodiaz*biodiaz/(gbio+epsbio*biodiaz*biodiaz)
        graz_D = g_D*bioz
        morp_D = nupt*biodiaz ! linear mortality
        no3upt_D = biono3/(k1n + biono3)*npp_D ! nitrate uptake
# endif
!       grazing function
        biop2 = biop*biop
        g_P = gbio*epsbio*biop2/(gbio+epsbio*biop2)
        graz = g_P*bioz
        morp = nup*biop2
        morpt = nupt*biop
        morz = nuz*bioz*bioz
        remi = nud*bct*biod
        excr = gamma2*bct*bioz
        expo = wwd*biod
!       flags prevent negative values by setting outgoing fluxes to
!       zero if tracers are lower than trcmin
        nflag = 0.5 + sign(0.5,bion - trcmin)
        pflag = 0.5 + sign(0.5,biop - trcmin)
        zflag = 0.5 + sign(0.5,bioz - trcmin)
        dflag = 0.5 + sign(0.5,biod - trcmin)
# if defined O_npzd_nitrogen
        no3flag = 0.5 + sign(0.5,biono3 - trcmin)
        diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
# endif
        graz = graz*pflag
        morp = morp*pflag
        morpt = morpt*pflag
        morz = morz*zflag
        remi = remi*dflag
        excr = excr*zflag
        expo = expo*dflag
# if defined O_npzd_nitrogen
        npp = npp*nflag*no3flag
        npp_D = npp_D*nflag
        graz_D = graz_D*diazflag
        morp_D = morp_D*diazflag
        no3upt_D = no3upt_D*no3flag
# else
        npp = npp*nflag
# endif

# if defined O_npzd_nitrogen
!       nutrients equation
        bion = bion + tsb*redptn*(remi + excr - (npp + npp_D) + morpt)
!       phytoplankton equation
        biop = biop + tsb*(npp - morp - graz - morpt)
!       zooplankton equation
        bioz = bioz + tsb*(gamma1*(graz + graz_D) - excr - morz)
!       detritus equation
        biod = biod + tsb*((1.-gamma1)*(graz + graz_D) + morp + morp_D
     &       + morz - remi - expo + impo)
!       nitrate (NO3) equation
        biono3 = biono3 + tsb*(remi + excr - npp + morpt - no3upt_D)
!       diazotroph equation
        biodiaz = biodiaz + tsb*(npp_D - morp_D - graz_D)
# else
!       nutrients equation
        bion = bion + tsb*redptn*(remi + excr - npp + morpt)
!       phytoplankton equation
        biop = biop + tsb*(npp - morp - graz - morpt)
!       zooplankton equation
        bioz = bioz + tsb*(gamma1*graz - excr - morz)
!       detritus equation
        biod = biod + tsb*((1.-gamma1)*graz + morp
     &       + morz - remi - expo + impo)
# endif
        expoout = expoout + expo
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
# if defined O_save_npzd
        nppout = nppout + npp
        morptout = morptout + morpt
        remiout = remiout + remi
        excrout = excrout + excr
#  if defined O_npzd_nitrogen
        npp_Dout = npp_Dout + npp_D
        graz_Dout = graz_Dout + graz_D
        morp_Dout = morp_Dout + morp_D
        nfixout = nfixout + npp_D - no3upt_D
#  endif
# endif

      enddo

      bioout(1) = bion - bioin(1)
      bioout(2) = biop - bioin(2)
      bioout(3) = bioz - bioin(3)
      bioout(4) = biod - bioin(4)
# if defined O_npzd_nitrogen
      bioout(5) = biono3 - bioin(5)
      bioout(6) = biodiaz - bioin(6)
# endif
#endif
      return
      end
