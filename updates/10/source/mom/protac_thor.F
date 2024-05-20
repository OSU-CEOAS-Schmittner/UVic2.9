#include "CPP_OPTIONS_CHECK.h"

      subroutine protac_thor_init

#if defined O_PaTh

!     initialize Pa-Th parameters
      implicit none
#include "size.h"
#if defined O_mobi
#include "mobi.h"
#endif
#include "protac_thor.h"
#include "calendar.h"
#include "coord.h"
#include "grdvar.h"
#include "stdunits.h"
#include "scalar.h"
!     IO
#include "param.h"
#include "pconst.h"
      integer ib(10), ic(10), iou, j, i
      character (120) :: fname, new_file_name
      logical exists
      real tmpijkm(1:imtm2,1:jmtm2,1:km)
! end IO

      integer k, ioun

      namelist /path/ PaKref, ThKref, PaSPMexponent, ThSPMexponent
     &,        KPaPOMFac,KPaCaCO3Fac,KPaOpalFac,KPaDustFac,KPaLithFac
     &,        KThPOMFac,KThCaCO3Fac,KThOpalFac,KThDustFac,KThLithFac
     &,        wDust 
     &,        OrgMatterToOrgCarbonMassRatio

      write(*,*)'Initializing Pa-Th module ...'

#if defined O_mobi && !defined O_mobi_caco3
	  print*,"Error in S/R path_init: O_mobi_caco3 must be switched on"
	  stop 'path_init'
#endif
      
      betaPa=2.33e-3 ! production rate in ocean in dpm/m^3/y
      betaTh=2.52e-2 ! production rate in ocean in dpm/m^3/y    
      lambdaDecayPa=2.13e-5 ! y^-1
      lambdaDecayTh=9.22e-6 ! y^-1

!     Convert units
      lambdaDecayPa=lambdaDecayPa/(86400.0*365.0) ! y^-1 -> s^-1
      lambdaDecayTh=lambdaDecayTh/(86400.0*365.0) ! y^-1 -> s^-1

      fmoltodpmPa=((1e-15)*avogradroNumber*lambdaDecayPa*60.0) ! fmol Pa to dpm
      fmoltodpmTh=((1e-15)*avogradroNumber*lambdaDecayTh*60.0) ! fmol Th to dpm

      betaPa=betaPa/(86400.0*365.0)/fmoltodpmPa ! production rate: dpm/m^3/y -> fmol/m^3/sec
      betaTh=betaTh/(86400.0*365.0)/fmoltodpmTh ! production rate: dpm/m^3/y -> fmol/m^3/sec

      write(*,*)'betaPa=',betaPa,' fmol/m^3/sec'
      write(*,*)'betaTh=',betaTh,' fmol/m^3/sec'
      write(*,*)'lambdaDecayPa=',lambdaDecayPa,' 1/sec'
      write(*,*)'lambdaDecayTh=',lambdaDecayTh,' 1/sec'
      write(*,*)'fmoltodpmPa=',fmoltodpmPa,' fmol/dpm'
      write(*,*)'fmoltodpmTh=',fmoltodpmTh,' fmol/dpm'
      
      do k=1,km
        dzmr(k)=dztr(k)*100. ! reciprocal of cell thickness in m^-1
      enddo  
      
!     set defaults for namelist path
      PaKref=1.e7
      ThKref=1.e7
C     ThCref=1.e-4 ! kg/m^3
C     Note: these are specified as positive values; we multiply by -1 below
      PaSPMexponent=0.66 ! Hayes et al., Marine Chem., 2015
      ThSPMexponent=0.66 ! Hayes et al., Marine Chem., 2015
      KPaPOMFac=1.0
      KPaCaCO3Fac=1.0/40.0
      KPaOpalFac=1.0/6.0
      KPaDustFac=1.0/20.0      
      KPaLithFac=1.0/20.0
      KThPOMFac=1.0
      KThCaCO3Fac=1.0
      KThOpalFac=1.0/20.0
      KThDustFac=1.0/20.0
      KThLithFac=1.0/20.0
      
!
      wDust=1000.0 ! m/y (converted to m/s below)
!     mass ratio of org C to org matter :      
      OrgMatterToOrgCarbonMassRatio=2.199 ! Klaas and Archer (2002)
!     An alternate value from SLAMS3D (Anna Rufas): 1.0/0.54
!     Or: 2.7 ! Panel 6.2.1 of Sarmiento & Gruber

!     read namelist
      call getunit(ioun, 'control.in', 'f s r')
      read(ioun, path, end=108)
108   continue
      write(stdout, path)
      call relunit(ioun)

C     write(*,*)'Constant reference scavenging coefficient for Pa used'
C     write(*,*)' PaKref=',PaKref
C     write(*,*)'  PaKref will be multipled by:'
C     write(*,*)'    KPaPOMFac=',KPaPOMFac
C     write(*,*)'    KPaCaCO3Fac=',KPaCaCO3Fac
C     write(*,*)'    KPaOpalFac=',KPaOpalFac
C     write(*,*)'    KPaDustFac=',KPaDustFac
C     write(*,*)'    KPaLithFac=',KPaLithFac
C     KPaPOM(:)=PaKref*KPaPOMFac
C     KPaCaCO3(:)=PaKref*KPaCaCO3Fac
C     KPaOpal(:)=PaKref*KPaOpalFac
C     KPaDust(:)=PaKref*KPaDustFac
C     KPaLith(:)=PaKref*KPaLithFac
C     write(*,*)'  Pa scavenging coefficients: '
C     write(*,*)'    KPaPOM=',KPaPOM(1)
C     write(*,*)'    KPaCaCO3=',KPaCaCO3(1)
C     write(*,*)'    KPaOpal=',KPaOpal(1)
C     write(*,*)'    KPaDust=',KPaDust(1)
C     write(*,*)'    KPaLith=',KPaLith(1)

#if defined O_Pa_variable_scavenging_coeff
      write(*,*)'Variable reference scavenging coefficient for Pa used'
C     write(*,*)' PaCref [kg/m^3]=',PaCref
      PaSPMexponent=-PaSPMexponent;
      write(*,*)' PaSPMexponent=',PaSPMexponent      
#else
      write(*,*)'Constant reference scavenging coefficient for Pa used'
      write(*,*)' PaKref=',PaKref
#endif
      write(*,*)'  Reference lithogenic scavenging coefficient for Pa'
      write(*,*)'    PaKref=',PaKref
      write(*,*)'  PaKref will be multipled by:'
      write(*,*)'    KPaPOMFac=',KPaPOMFac
      write(*,*)'    KPaCaCO3Fac=',KPaCaCO3Fac
      write(*,*)'    KPaOpalFac=',KPaOpalFac
      write(*,*)'    KPaDustFac=',KPaDustFac
      write(*,*)'    KPaLithFac=',KPaLithFac
#if !defined O_Pa_variable_scavenging_coeff
      KPaPOM(:)=PaKref*KPaPOMFac
      KPaCaCO3(:)=PaKref*KPaCaCO3Fac
      KPaOpal(:)=PaKref*KPaOpalFac
      KPaDust(:)=PaKref*KPaDustFac
      write(*,*)'  Pa scavenging coefficients: '      
      write(*,*)'    KPaPOM=',KPaPOM(1)
      write(*,*)'    KPaCaCO3=',KPaCaCO3(1)
      write(*,*)'    KPaOpal=',KPaOpal(1)
      write(*,*)'    KPaDust=',KPaDust(1)
#endif
      KPaLith(:)=PaKref*KPaLithFac
      write(*,*)'  Pa lithogenic scavenging coefficient: '
      write(*,*)'    KPaLith=',KPaLith(1)

#if defined O_Th_variable_scavenging_coeff
      write(*,*)'Variable reference scavenging coefficient for Th used'
C     write(*,*)' ThCref [kg/m^3]=',ThCref
      ThSPMexponent=-ThSPMexponent
      write(*,*)' ThSPMexponent=',ThSPMexponent      
#else
      write(*,*)'Constant reference scavenging coefficient for Th used'
      write(*,*)' ThKref=',ThKref
#endif
      write(*,*)'  Reference lithogenic scavenging coefficient for Th'
      write(*,*)'    ThKref=',ThKref
      write(*,*)'  ThKref will be multipled by:'
      write(*,*)'    KThPOMFac=',KThPOMFac
      write(*,*)'    KThCaCO3Fac=',KThCaCO3Fac
      write(*,*)'    KThOpalFac=',KThOpalFac
      write(*,*)'    KThDustFac=',KThDustFac
      write(*,*)'    KThLithFac=',KThLithFac
#if !defined O_Th_variable_scavenging_coeff
      KThPOM(:)=ThKref*KThPOMFac
      KThCaCO3(:)=ThKref*KThCaCO3Fac
      KThOpal(:)=ThKref*KThOpalFac
      KThDust(:)=ThKref*KThDustFac
      write(*,*)'  Th scavenging coefficients: '      
      write(*,*)'    KThPOM=',KThPOM(1)
      write(*,*)'    KThCaCO3=',KThCaCO3(1)
      write(*,*)'    KThOpal=',KThOpal(1)
      write(*,*)'    KThDust=',KThDust(1)
#endif
      KThLith(:)=ThKref*KThLithFac
      write(*,*)'  Th lithogenic scavenging coefficient: '
      write(*,*)'    KThLith=',KThLith(1)

      wDust=wDust/(86400.0*365.0) ! dust sinking speed in m/s

#if defined O_mobi
!     detritus in MOBI is in mmolN/m^3. We multiply it by redctn (molC/mmolN) 
!     to convert to molC/m^3.
      detrToPOM=redctn*MwC*OrgMatterToOrgCarbonMassRatio ! convert mmolN -> molC -> kg C -> kg POM
      write(*,*)'Using a DETR to POM conversion factor of ',detrToPOM

      do k=1,km
C       wPOM(k)=1000.0/(86400.0*365.0) !(wd(k)*1.e-2)*dzt(k) ! sinking speed in m/s
C       wCaCO3(k)=1000.0/(86400.0*365.0) !(wc(k)*1.e-2)*dzt(k) ! sinking speed in m/s
C       wOpal(k)=1000.0/(86400.0*365.0) !(wo(k)*1.e-2)*dzt(k) ! sinking speed in m/s
        wPOM(k)=(wd(k)*1.e-2)*dzt(k) ! sinking speed in m/s
        wCaCO3(k)=(wc(k)*1.e-2)*dzt(k) ! sinking speed in m/s
        wOpal(k)=(wo(k)*1.e-2)*dzt(k) ! sinking speed in m/s
      enddo  
# if !defined O_TMM
C       write out w profiles to netcdf?
# endif      
#else
! If using TMM, the sinking speeds and particle fields are set in S/R MOBI_INI.
# if !defined O_TMM 
      PaTh_pom(:,:,:) = 0.0
      PaTh_caco3(:,:,:) = 0.0
      PaTh_opal(:,:,:) = 0.0
      wPOM(:) = 1000.0/(86400.0*365.0) ! sinking speed in m/s
      wCaCO3(:) = 1000.0/(86400.0*365.0) ! sinking speed in m/s
      wOpal(:) = 1000.0/(86400.0*365.0) ! sinking speed in m/s

      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = km

      fname = new_file_name ('O_particle_fields.nc')
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         call openfile (trim(fname), iou)
         call getvara ('O_pom', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1, c0)
         PaTh_pom(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)
		 do k=1,km
			do j=1,jmt
			   PaTh_pom(1,j,k) = PaTh_pom(imtm1,j,k)
			   PaTh_pom(imt,j,k) = PaTh_pom(2,j,k)
			enddo
			do i=1,imt
			   PaTh_pom(i,1,k) = PaTh_pom(i,2,k)
			   PaTh_pom(i,jmt,k) = PaTh_pom(2,j,k)
			enddo
		 enddo

         call openfile (trim(fname), iou)
         call getvara ('O_caco3', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1, c0)
         PaTh_caco3(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)
		 do k=1,km
			do j=1,jmt
			   PaTh_caco3(1,j,k) = PaTh_caco3(imtm1,j,k)
			   PaTh_caco3(imt,j,k) = PaTh_caco3(2,j,k)
			enddo
			do i=1,imt
			   PaTh_caco3(i,1,k) = PaTh_caco3(i,2,k)
			   PaTh_caco3(i,jmt,k) = PaTh_caco3(2,j,k)
			enddo
		 enddo

         call openfile (trim(fname), iou)
         call getvara ('O_opal', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1, c0)
         PaTh_opal(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)
		 do k=1,km
			do j=1,jmt
			   PaTh_opal(1,j,k) = PaTh_opal(imtm1,j,k)
			   PaTh_opal(imt,j,k) = PaTh_opal(2,j,k)
			enddo
			do i=1,imt
			   PaTh_opal(i,1,k) = PaTh_opal(i,2,k)
			   PaTh_opal(i,jmt,k) = PaTh_opal(2,j,k)
			enddo
		 enddo

		 call closefile(iou)

      else
         print*,"Error in S/R path_init: Cannot find ", trim(fname)
         stop 'path_init'
      endif

      fname = new_file_name ('O_particle_fields.nc')
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         call openfile (trim(fname), iou)
		 ib(:) = 1
		 ic(:) = km
		 call getvara ('wPOM', iou, km, ib, ic, wPOM, c1, c0)
		 call getvara ('wCaCO3', iou, km, ib, ic, wCaCO3, c1, c0)
		 call getvara ('wOpal', iou, km, ib, ic, wOpal, c1, c0)

		 call closefile(iou)

      else
         print*,"Error in S/R path_init: Cannot find ", trim(fname)
         stop 'path_init'
      endif
# endif ! O_TMM
#endif        

! We always read the lith field
#if !defined O_TMM
      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = km

      PaTh_lith(:,:,:) = 0.0
      wLith(:) = 1000.0/(86400.0*365.0) ! sinking speed in m/s
      
      fname = new_file_name ('O_particle_fields.nc')
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         call openfile (trim(fname), iou)
         call getvara ('O_lith', iou, ic(1)*ic(2)*ic(3)
     &,                 ib, ic, tmpijkm, c1, c0)
         PaTh_lith(2:imtm1,2:jmtm1,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,:)
		 do k=1,km
			do j=1,jmt
			   PaTh_lith(1,j,k) = PaTh_lith(imtm1,j,k)
			   PaTh_lith(imt,j,k) = PaTh_lith(2,j,k)
			enddo
			do i=1,imt
			   PaTh_lith(i,1,k) = PaTh_lith(i,2,k)
			   PaTh_lith(i,jmt,k) = PaTh_lith(2,j,k)
			enddo
		 enddo

		 call closefile(iou)
            
      else
         print*,"Error in S/R path_init: Cannot find ", trim(fname)
         stop 'path_init'
      endif

      fname = new_file_name ('O_particle_fields.nc')
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         call openfile (trim(fname), iou)

		 ib(:) = 1
		 ic(:) = km
		 call getvara ('wLith', iou, km, ib, ic, wLith, c1, c0)

		 call closefile(iou)
            
      else
         print*,"Error in S/R path_init: Cannot find ", trim(fname)
         stop 'path_init'
      endif      
#endif /* O_TMM */

	  write(*,*)'Sinking speed profile for Pa-Th:'
	  write(*,*)'k  wPOM   wCaCO3 wOpal wDust wLith'
      do k=1,km
        write(*,*)'k=',k,wPOM(k),wCaCO3(k),wOpal(k),wDust,wLith(k)        
      enddo  

#endif
      return
      end
!     END path_init

      subroutine protac_thor_driver(
     &                 kmx, twodt
     &,                Pa, Th
     &,                POM, CaCO3, Opal, Lith, DustFlux
     &,                srcPa, srcTh
     &,                Pad, Pab, Thd, Thb
     &,                FPa_bot, FTh_bot
     &                 )

#if defined O_PaTh

      implicit none
#include "size.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#if defined O_mobi
#include "mobi.h"
#endif
#include "protac_thor.h"

      integer kmx
      real twodt
      real Pa(kmx), Th(kmx)
      real POM(kmx), CaCO3(kmx), Opal(kmx), Lith(kmx), DustFlux
      real srcPa(kmx), srcTh(kmx)
      real Pad(kmx), Pab(kmx), Thd(kmx), Thb(kmx)
      real FPa_bot, FTh_bot

      integer k, it
      real Krefloc, Ctotal, SPM
      real Dust(km)
      real Pam, Thm
      real SPa_pom, SPa_caco3, SPa_opal, SPa_lith, SPa_dust, SPa
      real STh_pom, STh_caco3, STh_opal, STh_lith, STh_dust, STh
      real Pa_pom, Pa_caco3, Pa_opal, Pa_lith, Pa_dust
      real Th_pom, Th_caco3, Th_opal, Th_lith, Th_dust
      real FPa_in, FPa_out, FTh_in, FTh_out
      real Paold(km), Thold(km), dtpath
      integer ntpath
      
!     Input POM, CaCO3, Opal and Lith are in kg/m^3
!     DustFlux is the surface flux in kg/m^2/s
      
      Dust(1:kmx)=DustFlux/wDust ! kg/m^3
      
      ntpath = 2
      dtpath = twodt/ntpath
      
      Paold(1:kmx)=Pa(1:kmx)
      Thold(1:kmx)=Th(1:kmx)
      
      FPa_bot=0.0
      FTh_bot=0.0

#if defined O_Th_variable_scavenging_coeff || defined O_Pa_variable_scavenging_coeff
      do k=1,kmx
        Ctotal = POM(k) + CaCO3(k) + Opal(k) + Lith(k) + Dust(k)
        SPM = 1.e9*Ctotal/rhosw ! normalized particle concentration in ug/kg (kg/m3 -> kg/kg -> ug/kg)
C       Krefloc = ((Ctotal/ThCref)**ThSPMexponent)*1.e7
C	KPaPOM(k)=Krefloc*KPaPOMFac
C	KPaCaCO3(k)=Krefloc*KPaCaCO3Fac
C	KPaOpal(k)=Krefloc*KPaOpalFac
C	KPaDust(k)=Krefloc*KPaDustFac
C	KPaLith(k)=Krefloc*KPaLithFac
#if defined O_Th_variable_scavenging_coeff
        Krefloc = (SPM**ThSPMexponent)*1.e7
		KThPOM(k)=Krefloc*KThPOMFac
		KThCaCO3(k)=Krefloc*KThCaCO3Fac
		KThOpal(k)=Krefloc*KThOpalFac
		KThDust(k)=Krefloc*KThDustFac
C	KThPOM(k)=Krefloc*KThPOMFac
C	KThCaCO3(k)=Krefloc*KThCaCO3Fac
C	KThOpal(k)=Krefloc*KThOpalFac
C	KThDust(k)=Krefloc*KThDustFac
#endif
#if defined O_Pa_variable_scavenging_coeff
        Krefloc = (SPM**PaSPMexponent)*1.e7
		KPaPOM(k)=Krefloc*KPaPOMFac
		KPaCaCO3(k)=Krefloc*KPaCaCO3Fac
		KPaOpal(k)=Krefloc*KPaOpalFac
		KPaDust(k)=Krefloc*KPaDustFac
#endif
C	KThLith(k)=Krefloc*KThLithFac
      enddo      
#endif

      do it=1,ntpath
C	Pam=max(Pa(1), pathtrcmin) ! total concentration in fmol/m^3
C	SPa_pom=KPaPOM(1)*max(POM(1),pathtrcmin)/rhosw
C	SPa_caco3=KPaCaCO3(1)*max(CaCO3(1),pathtrcmin)/rhosw
C	SPa_opal=KPaOpal(1)*max(Opal(1),pathtrcmin)/rhosw
C	SPa_lith=KPaLith(1)*max(Lith(1),pathtrcmin)/rhosw
C	SPa_dust=KPaDust(1)*max(Dust(1),pathtrcmin*1.e-3)/rhosw
C	SPa=SPa_pom+SPa_caco3+SPa_opal+SPa_lith+SPa_dust
C	Pad(1)=(1.0/(1.0+SPa))*Pam ! dissolved concentration in fmol/m^3
C	Pab(1)=Pam-Pad(1) ! particle-associated concentration in fmol/m^3
C	Pa_pom=SPa_pom*Pad(1) ! particle-associated concentration in fmol/m^3
C	Pa_caco3=SPa_caco3*Pad(1) ! particle-associated concentration in fmol/m^3
C	Pa_opal=SPa_opal*Pad(1) ! particle-associated concentration in fmol/m^3
C	Pa_lith=SPa_lith*Pad(1) ! particle-associated concentration in fmol/m^3
C	Pa_dust=SPa_dust*Pad(1) ! particle-associated concentration in fmol/m^3      
C	FPa_in=0.0
C	FPa_out = wPOM(1)*Pa_pom + wCaCO3(1)*Pa_caco3 + 
C    &            wOpal(1)*Pa_opal + wLith(1)*Pa_lith +
C    &            wDust*Pa_dust
C	srcPa(1) = betaPa - lambdaDecayPa*Pam + 
C    &             (FPa_in-FPa_out)*dzmr(1) ! source term in fmol/m^3/s
C
C	Thm=max(Th(1), pathtrcmin) ! total concentration in fmol/m^3
C	STh_pom=KThPOM(1)*max(POM(1),pathtrcmin)/rhosw
C	STh_caco3=KThCaCO3(1)*max(CaCO3(1),pathtrcmin)/rhosw
C	STh_opal=KThOpal(1)*max(Opal(1),pathtrcmin)/rhosw
C	STh_lith=KThLith(1)*max(Lith(1),pathtrcmin)/rhosw
C	STh_dust=KThDust(1)*max(Dust(1),pathtrcmin*1.e-3)/rhosw
C	STh=STh_pom+STh_caco3+STh_opal+STh_lith+STh_dust
C	Thd(1)=(1.0/(1.0+STh))*Thm ! dissolved concentration in fmol/m^3
C	Thb(1)=Thm-Thd(1) ! particle-associated concentration in fmol/m^3
C	Th_pom=STh_pom*Thd(1) ! particle-associated concentration in fmol/m^3
C	Th_caco3=STh_caco3*Thd(1) ! particle-associated concentration in fmol/m^3
C	Th_opal=STh_opal*Thd(1) ! particle-associated concentration in fmol/m^3
C	Th_lith=STh_lith*Thd(1) ! particle-associated concentration in fmol/m^3
C	Th_dust=STh_dust*Thd(1) ! particle-associated concentration in fmol/m^3      
C	FTh_in=0.0
C	FTh_out = wPOM(1)*Th_pom + wCaCO3(1)*Th_caco3 + 
C    &            wOpal(1)*Th_opal + wLith(1)*Th_lith +
C    &            wDust*Th_dust
C	srcTh(1) = betaTh - lambdaDecayTh*Thm + 
C    &             (FTh_in-FTh_out)*dzmr(1)

		FPa_in=0.0
		FTh_in=0.0
		do k=1,kmx
		  Pam=max(Pa(k), pathtrcmin) ! total concentration in fmol/m^3		
		  SPa_pom=KPaPOM(k)*max(POM(k),pathtrcmin)/rhosw
		  SPa_caco3=KPaCaCO3(k)*max(CaCO3(k),pathtrcmin)/rhosw
		  SPa_opal=KPaOpal(k)*max(Opal(k),pathtrcmin)/rhosw
		  SPa_lith=KPaLith(k)*max(Lith(k),pathtrcmin)/rhosw
		  SPa_dust=KPaDust(k)*max(Dust(k),pathtrcmin*1.e-3)/rhosw
		  SPa=SPa_pom+SPa_caco3+SPa_opal+SPa_lith+SPa_dust
		  Pad(k)=(1.0/(1.0+SPa))*Pam ! dissolved concentration in fmol/m^3
		  Pab(k)=Pam-Pad(k) ! particle-associated concentration in fmol/m^3
		  Pa_pom=SPa_pom*Pad(k) ! particle-associated concentration in fmol/m^3
		  Pa_caco3=SPa_caco3*Pad(k) ! particle-associated concentration in fmol/m^3
		  Pa_opal=SPa_opal*Pad(k) ! particle-associated concentration in fmol/m^3
		  Pa_lith=SPa_lith*Pad(k) ! particle-associated concentration in fmol/m^3
		  Pa_dust=SPa_dust*Pad(k) ! particle-associated concentration in fmol/m^3
		  if (k.gt.1) FPa_in = FPa_out
		  FPa_out = wPOM(k)*Pa_pom + wCaCO3(k)*Pa_caco3 + 
     &              wOpal(k)*Pa_opal + wLith(k)*Pa_lith +
     &              wDust*Pa_dust
C       if (k.eq.kmx) FPa_out = bottomFlux*FPa_out
		  srcPa(k) = betaPa - lambdaDecayPa*Pam + 
     &               (FPa_in-FPa_out)*dzmr(k) ! source term in fmol/m^3/s

		  Thm=max(Th(k), pathtrcmin) ! total concentration in fmol/m^3		
		  STh_pom=KThPOM(k)*max(POM(k),pathtrcmin)/rhosw
		  STh_caco3=KThCaCO3(k)*max(CaCO3(k),pathtrcmin)/rhosw
		  STh_opal=KThOpal(k)*max(Opal(k),pathtrcmin)/rhosw
		  STh_lith=KThLith(k)*max(Lith(k),pathtrcmin)/rhosw
		  STh_dust=KThDust(k)*max(Dust(k),pathtrcmin*1.e-3)/rhosw
		  STh=STh_pom+STh_caco3+STh_opal+STh_lith+STh_dust
		  Thd(k)=(1.0/(1.0+STh))*Thm ! dissolved concentration in fmol/m^3
		  Thb(k)=Thm-Thd(k) ! particle-associated concentration in fmol/m^3
		  Th_pom=STh_pom*Thd(k) ! particle-associated concentration in fmol/m^3
		  Th_caco3=STh_caco3*Thd(k) ! particle-associated concentration in fmol/m^3
		  Th_opal=STh_opal*Thd(k) ! particle-associated concentration in fmol/m^3
		  Th_lith=STh_lith*Thd(k) ! particle-associated concentration in fmol/m^3
		  Th_dust=STh_dust*Thd(k) ! particle-associated concentration in fmol/m^3
		  if (k.gt.1) FTh_in = FTh_out
		  FTh_out = wPOM(k)*Th_pom + wCaCO3(k)*Th_caco3 + 
     &              wOpal(k)*Th_opal + wLith(k)*Th_lith +
     &              wDust*Th_dust
C       if (k.eq.kmx) FTh_out = bottomFlux*FTh_out
		  srcTh(k) = betaTh - lambdaDecayTh*Thm + 
     &               (FTh_in-FTh_out)*dzmr(k)

		enddo

		Pa(:) = Pa(:) + dtpath*srcPa(:)
		Th(:) = Th(:) + dtpath*srcTh(:)

		FPa_bot = FPa_bot + dtpath*FPa_out
		FTh_bot = FTh_bot + dtpath*FTh_out
      
      enddo
 
      srcPa(1:kmx) = (Pa(1:kmx)-Paold(1:kmx))/twodt ! fmol/m^3/s
      srcTh(1:kmx) = (Th(1:kmx)-Thold(1:kmx))/twodt ! fmol/m^3/s
      FPa_bot = FPa_bot/twodt ! fmol/m^2/s
      FTh_bot = FTh_bot/twodt ! fmol/m^2/s
      
#endif                  
      return
      end
