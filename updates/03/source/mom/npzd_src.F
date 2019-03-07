      subroutine npzd_src (bioin, ntsb, tsb, gl, bct, impo, dzt
     &,                    dayfrac, wwd, rkw, nud, bioout, expoout
     &,                    grazout, morpout, morzout, graz_Det_out
     &,                    graz_Zout
# if defined O_save_npzd
     &,                    nppout, morptout, remiout, excrout
#  if defined O_npzd_nitrogen
     &,                    npp_Dout, graz_Dout, morp_Dout, nfixout
#  endif
#  if defined O_npzd_extra_diagnostics
     &,                    avej_out, avej_D_out, gmax_out, no3P_out
     &,                    po4P_out, po4_D_out
#  endif
# endif
#  if defined O_npzd_fe_limitation
     &,                    felimit, felimit_D
#  endif
     &,                    bctz
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
!     Schmittner et al. 2008, Global Biogeochem. Cycles 22, GB1013
!
!     This version was modified by David Keller and corrects the zooplankton
!     grazing formulation.  Note that zooplankton are now allowed to graze
!     on themselves and detritus, in addition to phyt. and diazotrophs.
!     The calculation of light has also been corrected.
!
!     Note that nutrient now represents phosphate

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
!       graz_Det_out = grazing of detritus
!       graz_Zout   = grazing on othe zooplankton
!       avej_out    = light-depend phyt. growth rate
!       avej_D_out  = light-depend Diaz growth rate
!       gmax_out    = temp-depend. zoo growth rate
!       no3P_out    = no3 depend. phyt growth rate
!       po4P_out    = po4 depend. phyt growth rate
!       po4_D_out   = po4 depend. Diaz growth rate

!      New grazing formulation variables and parameters

!       The following terms determine ingestion according to a
!       a Holling II curve (i.e. Michaelis Menten):
!
!       Ingestion = max_graz_rate * (Ft/(Ft + kzoo))
!
!       where Ft is the weighted measure of the total food available
!       and equals the sum of the different prey types times the
!       preference of Z for that type of prey 
!
!       zprefP   = Z preference for P
!       zprefD   = Z preference for Diaz
!       zprefDet = Z preference for detritus
!       zprefZ   = Z preference for other Z
!       kzoo = half saturation coefficienct for Z ingestion mmol N m-3
!       ing_P    = zooplankton ingestion of phytoplankon
!       ing_D    = zooplankton ingestion of diazotrophs
!       ing_Det  = zooplankton ingestion of detritus
!       ing_Z    = zooplankton ingestion of other zooplankton
!       thetaZ   = Michaelis-Menten denominator
!
!       felimit = Fe limitation parameter
!       felmit_D = Fe limitation parameter for diazotrophs
!
!=======================================================================

      implicit none

      integer n, ntsb

      real gl, f1, bion, biop, bioz, biod, jmax, u_P, g_P, npp, graz
      real morp, morpt, morz, remi, excr, expo, impo, nppout, grazout
      real morpout, morptout, morzout, remiout, excrout, expoout, tsb
      real avej_out, avej_D_out, gmax_out, no3P_out, po4P_out, po4_D_out
      real dzt, nflag, pflag, zflag, dflag, wwd, rkw, gd, dayfrac, bct
      real nupt, nud, biono3, u_D,npp_D, npp_Dout, no3flag, biodiaz
      real diazflag, g_D,graz_D, morp_D, jmax_D, gd_D, avej_D, no3upt_D
      real morp_Dout, graz_Dout, nfixout, biop2, u1, u2, phi1, phi2
      real avej, graz_Det_out, graz_Zout, thetaZ, ing_P, ing_D
      real ing_Det, ing_Z, g_Z, g_Det, graz_Z, graz_Det, gmax
      real no3P, po4P, po4_D, felimit, bctz, felimit_D

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "npzd.h"

      real bioin(ntnpzd), bioout(ntnpzd)

!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      f1 = exp((-kw - kc*(bioin(2)+bioin(6)))*dzt)
# if defined O_npzd_fe_limitation
! In the following "felimit" is determined by an iron mask and
! is used to limit phytoplankton growth in HNLC regions
      jmax = abio*bct*felimit
# else
      jmax = abio*bct
# endif
      gd = jmax*dayfrac
      u1 = max(gl/gd,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
# if defined O_npzd_cdom_attenuation
      avej = gd*(phi1 - phi2)/(1.2*(kw+kc*(bioin(2)+bioin(6)))*dzt)
# else
      avej = gd*(phi1 - phi2)/((kw+kc*(bioin(2)+bioin(6)))*dzt)
# endif
! Make the max grazing rate a function of temperature
      gmax = gbio*bctz
! Note bctz, which sets an upper limit on the effects of temp on the
! grazing rate, is set in tracers.F
# if defined O_npzd_nitrogen
# if defined O_npzd_fe_limitation
      jmax_D = max(0.,abio*(bct - 2.6)*felimit_D)*jdiar
# else
      jmax_D = max(0.,abio*(bct - 2.6))*jdiar
# endif
      gd_D = max(1.e-14,jmax_D*dayfrac)
      u1 = max(gl/gd_D,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*(phi1 - phi2)/((kw+kc*(bioin(2)+bioin(6)))*dzt)
# endif
# if defined O_npzd_nitrogen
!     check grazing preferences = 1 for N case
        IF ((zprefP + zprefDet + zprefZ + zprefD).ne.1) THEN
           zprefP = 0.30
           zprefZ = 0.30
           zprefDet = 0.30
           zprefD = 0.10
        END IF
# else
!      check that grazing preferences = 1 for no N case
        IF ((zprefP + zprefDet + zprefZ).ne.1) THEN
           zprefP = 0.35
           zprefZ = 0.35
           zprefDet = 0.30
        END IF
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
      graz_Det_out = 0.0
      graz_Zout = 0.0
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
# if defined O_npzd_extra_diagnostics
      avej_out = 0.0
      avej_D_out = 0.0
      gmax_out = 0.0
      no3P_out = 0.0
      po4P_out = 0.0
      po4_D_out = 0.0
# endif
# endif

      do n=1,ntsb

!       growth rate of phytoplankton
        u_P = min(avej, jmax*bion/(k1p + bion))
        po4P = jmax*bion/(k1p + bion) 

# if defined O_npzd_nitrogen
!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
        no3P = jmax*biono3/(k1n + biono3)
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        u_D = min(avej_D, jmax_D*bion/(k1p + bion))
        po4_D = jmax_D*bion/(k1p + bion)
!       Set the grazing coefficients for the N case
        thetaZ = zprefP*biop+zprefDet*biod+zprefZ*bioz+zprefD*biodiaz
     &         + kzoo
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
        ing_D = zprefD/thetaZ
# else
!       If "else" then set the grazing coefficients for the no N case
!       note kzoo is in terms of N so convert to P
        thetaZ = zprefP*biop+zprefDet*biod+zprefZ*bioz+kzoo*redptn
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
# endif
        npp = u_P*biop
# if defined O_npzd_nitrogen
        npp_D = max(0.,u_D*biodiaz)
!       grazing on diazotrophs
        g_D = gmax*ing_D*biodiaz
        graz_D = g_D*bioz
        morp_D = nupt*biodiaz ! linear mortality
        no3upt_D = biono3/(k1n + biono3)*npp_D ! nitrate uptake 
# endif
!       grazing on P
        g_P = gmax*ing_P*biop
        graz = g_P*bioz
!       grazing on Z
        g_Z = gmax*ing_Z*bioz
        graz_Z = g_Z*bioz
!       grazing on Detritus
        g_Det = gmax*ing_Det*biod
        graz_Det = g_Det*bioz
!
        morp = nup*biop
        morpt = nupt*biop
        morz = nuz*bioz*bioz
        remi = nud*bct*biod
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
        graz = graz*pflag*zflag
        graz_Z = graz_Z*zflag
        graz_Det = graz_Det*dflag*zflag
        morp = morp*pflag
        morpt = morpt*pflag
        morz = morz*zflag
        remi = remi*dflag
        expo = expo*dflag
# if defined O_npzd_nitrogen
        npp = npp*nflag*no3flag
        npp_D = npp_D*nflag
        graz_D = graz_D*diazflag*zflag
        morp_D = morp_D*diazflag
        no3upt_D = no3upt_D*no3flag
# else
        npp = npp*nflag
# endif

!   Excretion is the difference between
!   the zooplankton assimilation and growth efficiencies
# if defined O_npzd_nitrogen
        excr = (gamma1-geZ)*(graz+graz_Z+graz_Det+graz_D)
# else
        excr = (gamma1-geZ)*(graz+graz_Z+graz_Det)
#endif

# if defined O_npzd_nitrogen
!       nutrients equation
        bion = bion + tsb*redptn*(remi + excr - (npp + npp_D) + morpt)
!       phytoplankton equation
        biop = biop + tsb*(npp - morp - graz - morpt)
!       zooplankton equation
        bioz = bioz + tsb*(geZ*(graz + graz_D + graz_Det + graz_Z) 
     &       - morz - graz_Z)
!       detritus equation
        biod = biod + tsb*((1.-gamma1)*(graz + graz_D + graz_Det 
     &       + graz_Z) + morp + morp_D + morz - remi - graz_Det 
     &       - expo + impo)
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
        bioz = bioz + tsb*(geZ*(graz + graz_Det + graz_Z) - morz
     &       - graz_Z)
!       detritus equation
        biod = biod + tsb*((1.-gamma1)*(graz + graz_Det + graz_Z) + morp
     &       + morz - remi - graz_Det -expo + impo)
# endif
        expoout = expoout + expo
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
        graz_Det_out = graz_Det_out + graz_Det
        graz_Zout = graz_Zout + graz_Z
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
# if defined O_npzd_extra_diagnostics
        avej_out = avej_out + avej
        avej_D_out = avej_D_out + avej_D
        gmax_out = gmax_out + gmax
        no3P_out = no3P_out + no3P
        po4P_out = po4P_out + po4P
        po4_D_out = po4_D_out + po4_D
# endif
# endif
      enddo

      bioout(1) = bion - max(bioin(1),trcmin)
      bioout(2) = biop - max(bioin(2),trcmin)
      bioout(3) = bioz - max(bioin(3),trcmin)
      bioout(4) = biod - max(bioin(4),trcmin)
# if defined O_npzd_nitrogen
      bioout(5) = biono3 - max(bioin(5),trcmin)
      bioout(6) = biodiaz - max(bioin(6),trcmin)
# endif
#endif
      return
      end
