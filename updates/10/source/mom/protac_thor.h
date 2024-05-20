!====================== include file "protac_thor.h" =========================

!   variables for Pa-Th module

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			 
#if defined O_PaTh

      real pathtrcmin
      parameter (pathtrcmin=5e-12)

      real rhosw ! nominal density of seawater kg/m^3
      parameter (rhosw=1024.5)

      real avogradroNumber
      parameter (avogradroNumber=6.02214076e23)
      
      real MwPa, MwTh ! atomic mass in g/mol (https://ciaaw.org/)
      parameter (MwPa=231.03588, MwTh=230.033132)

!     MwOpal from https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GB005186
      real MwC, MwCaCO3, MwOpal ! molar mass of C, CaCO3 and biogenic silica in kg
      parameter (MwC=12.01e-3, MwCaCO3=100.1e-3, MwOpal=67.3e-3)

      real betaPa, betaTh ! production rate in ocean
      real lambdaDecayPa, lambdaDecayTh ! radioactive decay constant

      real OrgMatterToOrgCarbonMassRatio, detrToPOM
      real fmoltodpmPa, fmoltodpmTh

      real PaKref, ThKref, PaSPMexponent, ThSPMexponent
      real KPaPOMFac,KPaCaCO3Fac,KPaOpalFac,KPaDustFac,KPaLithFac
      real KThPOMFac,KThCaCO3Fac,KThOpalFac,KThDustFac,KThLithFac      
      real KPaPOM(km),KPaCaCO3(km),KPaOpal(km),KPaDust(km),KPaLith(km)
      real KThPOM(km),KThCaCO3(km),KThOpal(km),KThDust(km),KThLith(km)
      real wPOM(km), wCaCO3(km), wOpal(km), wDust, wLith(km)
      real dzmr(km)

      real PaTh_lith(imt,jmt,km)
#if !defined O_mobi
      real PaTh_pom(imt,jmt,km), PaTh_caco3(imt,jmt,km), 
     &     PaTh_opal(imt,jmt,km)
#endif      

      common /path_r/ betaPa, betaTh
      common /path_r/ lambdaDecayPa, lambdaDecayTh
      common /path_r/ OrgMatterToOrgCarbonMassRatio, detrToPOM
      common /path_r/ fmoltodpmPa, fmoltodpmTh 
      common /path_r/ PaKref, ThKref, PaSPMexponent, ThSPMexponent
      common /path_r/ KPaPOMFac,KPaCaCO3Fac,KPaOpalFac,KPaDustFac
      common /path_r/ KPaLithFac
      common /path_r/ KThPOMFac,KThCaCO3Fac,KThOpalFac,KThDustFac
      common /path_r/ KThLithFac
      common /path_r/ KPaPOM,KPaCaCO3,KPaOpal,KPaDust,KPaLith
      common /path_r/ KThPOM,KThCaCO3,KThOpal,KThDust,KThLith
      common /path_r/ wPOM,wCaCO3,wOpal,wDust,wLith
      common /path_r/ dzmr
      common /path_r/ PaTh_lith
#if !defined O_mobi
      common /path_r/ PaTh_pom, PaTh_caco3, PaTh_opal
#endif
      
#endif
