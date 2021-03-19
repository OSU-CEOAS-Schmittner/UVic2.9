
#if defined O_mom

      subroutine tracer (joff, js, je, is, ie)

!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer istrt, iend, i, k, j, ip, kr, jq, n, jp, jrow, iou, js
      integer je, limit, joff, is, ie, kmx, m, kb, idiag, index
      integer it(10), iu(10), ib(10), ic(10), nfnpzd, mfnpzd, mxfnpzd
      integer id_time, id_xt, id_yt, id_zt,fe_jlo,fe_m,fe_k,fe_n

      logical inqvardef, exists

      parameter (fe_n = 14)

      real time
      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, zmax, cont, drho
      real drhom1, wt, ahbi_cstr, ahbi_csu_dyur, gamma, fy, fyz
      real fe_dy, fe_conc, fe_x(fe_n), fe_y(fe_n), bctz

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "timeavgs.h"
#include "tmngr.h"
#include "vmixc.h"
# if defined O_save_convection || defined O_carbon_14
#include "diaga.h"
      real rrc14std
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
# if defined O_mobi || defined O_carbon_14
#include "mobi.h"
# endif
# if defined O_mobi
      real t_in(km), po4_in(km)
#  if defined O_mobi_o2
      real o2_in(km)
#  endif
#  if defined O_carbon_13 || defined O_mobi_caco3
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
#  endif
# if defined O_mobi_silicon
      real sil_in(km)
      real opalburial, silhydrflux
# endif
      real sgb_in(km)
#  if defined O_mobi_nitrogen
      real no3_in(km)
#   if defined O_mobi_nitrogen_15
      real din15_in(km)
#   endif
#  endif     
      real expo, tnpzd(km,ntnpzd)
#  if defined O_embm
#include "atm.h"
#   if defined O_carbon_13 || defined O_mobi_caco3
#include "cembm.h"
#   endif
#  endif
# endif

# if defined O_carbon_fnpzd
#include "calendar.h"
# endif

# if defined O_matrix
#include "matrix.h"
      integer ntile
      logical isTile
# endif

#ifndef O_TMM
      parameter (istrt=2, iend=imt-1)
#else
      parameter (istrt=1, iend=imt)
#endif

      real twodt(km)
# if defined O_plume
      real  work(imt,km,jsmw:jemw)
# endif
# if defined O_mobi || defined O_carbon_14
      real ai, hi, hs, rctheta, dayfrac, swr, declin
      real src(imt,km,jsmw:jemw,nsrc)
#  ifdef O_TMM
      common/mobicomm/src
#  endif 
#  if defined O_sed
      real sedcorgflx, sedcalflx
#  endif
# endif
# if defined O_carbon_fnpzd
      real tmpijk(imtm2,jmtm2,km), tmpi(imtm2), tmpj(jmtm2)
      real, allocatable :: fnpzd(:,:,:,:)
      save fnpzd, nfnpzd, mfnpzd, mxfnpzd
# endif
# if defined O_isopycmix
#include "isopyc.h"
# endif
#include "fdift.h"

# if defined O_carbon_fnpzd
!-----------------------------------------------------------------------
!     create file for writing npzd 3d fluxes of carbon and alkalinity
!-----------------------------------------------------------------------
      if (.not. allocated (fnpzd)) then
        allocate ( fnpzd(imt,jmt,km,2) )
        fnpzd(:,:,:,:) = 0.
!       initialize time record counter
        mfnpzd = 0
!       initialize records written counter
        nfnpzd = 0
!       set the repeating time interval to one year
        mxfnpzd = nint(daylen*yrlen/dtts)
        fname = new_file_name ("O_fnpzd.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) then
!         limit time interval to what exists
          call openfile (fname, iou)
          call getdimlen ('time', iou, nfnpzd)
          mxfnpzd = min(mxfnpzd, nfnpzd)
#  if defined O_carbon
          if (inqvardef('O_fnpzd1', iou)) nfnpzd = mxfnpzd
#  endif
#  if defined O_mobi_alk
          if (inqvardef('O_fnpzd2', iou)) nfnpzd = mxfnpzd
#  endif
        else
!         define file if not there
          call openfile (fname, iou)
          call redef (iou)
          call defdim ('time', iou, 0, id_time)
          call defdim ('longitude', iou, imtm2, id_xt)
          call defdim ('latitude', iou, jmtm2, id_yt)
          call defdim ('depth', iou, km, id_zt)
          it(1) = id_time
          call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,     'time', 'time', 'years since 0-1-1')
          call putatttext (iou, 'time', 'calendar', calendar)
          it(1) = id_xt
          call defvar ('longitude', iou, 1, it, c0, c0, 'X', 'D'
     &,     'longitude', 'longitude', 'degrees')
          it(1) = id_yt
          call defvar ('latitude', iou, 1, it, c0, c0, 'Y', 'D'
     &,     'latitude', 'latitude', 'degrees')
          it(1) = id_zt
          call defvar ('depth', iou, 1, it, c0, c0, 'Z', 'D'
     &,     'depth', 'depth', 'm')
          it(1) = id_xt
          it(2) = id_yt
          it(3) = id_zt
          it(4) = id_time
#  if defined O_carbon
          call defvar ('O_fnpzd1', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
#  if defined O_mobi_alk
          call defvar ('O_fnpzd2', iou, 4, it, c0, c0, ' ', 'D'
     &,       ' ', ' ', ' ')
#  endif
          call enddef (iou)
          ib(:) = 1
          ic(:) = imtm2
          tmpi(1:imtm2) = xt(2:imtm1)
          call putvara ('longitude', iou, imtm2, ib, ic, xt, c1, c0)
          ib(:) = 1
          ic(:) = jmtm2
          tmpj(1:jmtm2) = yt(2:jmtm1)
          call putvara ('latitude', iou, jmtm2, ib, ic, yt, c1, c0)
          ib(:) = 1
          ic(:) = km
          call putvara ('depth', iou, km, ib, ic, zt, c1, c0)
        endif
      endif
# endif

# ifndef O_TMM
!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

#  if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic || defined O_isopycmix
      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
#  else
      do j=js,je
#  endif
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
#  if defined O_consthmix && !defined O_bryan_lewis_horizontal && !defined O_biharmonic
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
#  endif
        enddo
      enddo
#  if defined O_plume

!-----------------------------------------------------------------------
!     find depth of penetration for the llnl simple plume model
!-----------------------------------------------------------------------

!     set parameters
      zmax = 160.0e2 ! use constant depth penetration (if cont=0)
!      cont = 0.0    ! use rho contrast for penetration (if cont>0)
      cont = 0.4e-3  ! 0.4 kg/m3 = 0.4e-3g/cm3

      do j=js,je
        jrow = j + joff
        do i=is,ie
          subz(i,jrow) = zmax
        enddo
        if (cont .gt. c0) then
          call state_ref (t(1,1,1,1,tau), t(1,1,1,2,tau)
     &,                   work(1,1,jsmw), j, j, 1, imt, 1)
          do i=is,ie
            kmx = kmt(i,jrow)
            subz(i,jrow) = zw(1)
            do k=2,kmx
              drho = work(i,k,j) - work(i,1,j)
              drhom1 = work(i,k-1,j) - work(i,1,j)
              if (drho .ge. cont .and. drhom1 .lt. cont)
     &          subz(i,jrow) = zt(k-1) + (zt(k) - zt(k-1))
     &                      *(cont - drhom1)/(drho - drhom1)
            enddo
            if (drho .le. cont) subz(i,jrow) = zw(kmx)
          enddo
        endif

!-----------------------------------------------------------------------
!     calculate "3d stf" multiplier for the llnl simple plume model
!-----------------------------------------------------------------------

        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            subz(i,jrow) = max(min(subz(i,jrow),zw(kmt(i,jrow))),zw(1))
            work(i,1,j) = c1/subz(i,jrow)
            do k=2,km
              wt =  min(max(c0,(subz(i,jrow) - zw(k-1))),dzt(k))*dztr(k)
              work(i,k,j) = wt*work(i,1,j)
            enddo
          else
            do k=1,km
              work(i,k,j) = c0
            enddo
          endif
        enddo
      enddo
#  endif

# endif ! SPK not O_TMM

# if defined O_mobi
!-----------------------------------------------------------------------
!     ocean biogeochemistry and isotopes
!-----------------------------------------------------------------------
      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      nbio = c2dtts/dtnpzd 
      dtbio = c2dtts/nbio
      rdtts = 1./c2dtts
      rnbio = 1./nbio

#  if defined O_mobi_silicon
      rivsil = 0.0        !river Si deposition mol si/s??
CSPKKK I've removed the variable silwflx because it was confusing. In it's place 
C      I've created two new ones opalburial and silhydrflux. 
C      kk's silwflx = opalburial - silhydrflux
C      The units really need to be checked!
      opalburial = 0.0    ! Opal flux through bottom umol si/s
      silhydrflux = 0.0       !hydrothermal input umol si/s
#  endif  

      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
#  if defined O_ice
#   if defined O_ice_cpts
            ai = 0.
            hi = 0.
            hs = 0.
            do n=1,ncat
               ai =  ai + A(i,jrow,2,n)
               hi =  hi + heff(i,jrow,2,n)
               hs =  hs + hseff(i,jrow,2,n)
            enddo
#   else
            ai = aice(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsno(i,jrow,2)
#   endif
#  else
            ai = 0.
            hi = 0.
            hs = 0.
#  endif
!     calculate day fraction and incoming solar
!     angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
#  if defined O_embm
            swr = tap*dnswr(i,jrow)*1e-3*(1. + ai*(exp(-ki*(hi + hs))
     &                                          - 1.))
#  else
            swr = tap*200.            
#  endif

!        set prognostic variables
            tnpzd(:,imobipo4) = t(i,:,j,ipo4,taum1)
            tnpzd(:,imobiphyt) = t(i,:,j,iphyt,taum1)
            tnpzd(:,imobizoop) = t(i,:,j,izoop,taum1)
            tnpzd(:,imobidetr) = t(i,:,j,idetr,taum1)
#  if defined O_carbon            
            tnpzd(:,imobidic) = t(i,:,j,idic,taum1)
#  endif            
#  if defined O_mobi_nitrogen
            tnpzd(:,imobidop) = t(i,:,j,idop,taum1)            
            tnpzd(:,imobino3) = t(i,:,j,ino3,taum1)
            tnpzd(:,imobidon) = t(i,:,j,idon,taum1)
            tnpzd(:,imobidiaz) = t(i,:,j,idiaz,taum1)
#   if defined O_mobi_nitrogen_15
            tnpzd(:,imobidin15) = t(i,:,j,idin15,taum1)
            tnpzd(:,imobidon15) = t(i,:,j,idon15,taum1)
            tnpzd(:,imobiphytn15) = t(i,:,j,iphytn15,taum1)
            tnpzd(:,imobizoopn15) = t(i,:,j,izoopn15,taum1)
            tnpzd(:,imobidetrn15) = t(i,:,j,idetrn15,taum1)
            tnpzd(:,imobidiazn15) = t(i,:,j,idiazn15,taum1)
#    if defined O_mobi_silicon
            tnpzd(:,imobidiatn15) = t(i,:,j,idiatn15,taum1)
#    endif
#   endif
#  endif
#  if defined O_mobi_caco3
            tnpzd(:,imobicaco3) = t(i,:,j,icaco3,taum1)
#  endif
#  if defined O_kk_ballast
            tnpzd(:,imobidetr_B) = t(i,:,j,idetr_B,taum1) 
#  endif
#  if defined O_mobi_silicon
            tnpzd(:,imobidiat) = t(i,:,j,idiat,taum1)
            tnpzd(:,imobisil) = t(i,:,j,isil,taum1)
            tnpzd(:,imobiopl) = t(i,:,j,iopl,taum1)
#  endif
#  if defined O_mobi_iron
            tnpzd(:,imobidfe) = t(i,:,j,idfe,taum1)
            tnpzd(:,imobidetrfe) = t(i,:,j,idetrfe,taum1)
#  endif
#  if defined O_carbon_13
            tnpzd(:,imobidic13) = t(i,:,j,idic13,taum1)
            tnpzd(:,imobidoc13) = t(i,:,j,idoc13,taum1)
            tnpzd(:,imobiphytc13) = t(i,:,j,iphytc13,taum1)
            tnpzd(:,imobizoopc13) = t(i,:,j,izoopc13,taum1)
            tnpzd(:,imobidetrc13) = t(i,:,j,idetrc13,taum1)
            tnpzd(:,imobidiazc13) = t(i,:,j,idiazc13,taum1)
#   if defined O_mobi_silicon
            tnpzd(:,imobidiatc13) = t(i,:,j,idiatc13,taum1)
#   endif            
#   if defined O_mobi_caco3
            tnpzd(:,imobicaco3c13) = t(i,:,j,icaco3c13,taum1) 
#   endif            
#  endif
!        set other input variables
            t_in(:) = t(i,:,j,itemp,taum1) !degree Celsius
            po4_in(:) = t(i,:,j,ipo4,taum1)
#  if defined O_mobi_o2
            o2_in(:) = t(i,:,j,io2,taum1)*1000. ! uM
#  endif
#  if defined O_mobi_silicon
cAAA now in tnpzd
c            sil_in(:) =  t(i,:,j,isil,taum1)*1.e3 ! convert to mmol m-3 for npzd_src!
#  endif
#  if defined O_carbon_13 || defined O_mobi_caco3
            s_in(:) = 1.e3*t(i,:,j,isalt,taum1) + 35.0
            dic_in(:) = t(i,:,j,idic,taum1) !mol/m^3
            alk_in(:)  = t(i,:,j,ialk,taum1) !eq/m^3
#   if defined O_carbon_co2_2d
            co2_in = at(i,j,2,ico2)
#   else
            co2_in = co2ccn
#   endif
#   if defined O_carbon_13
            dic13_in(:) = tnpzd(:,imobidic13) !mol/m^3
#   endif
#  endif
#  if defined O_mobi_nitrogen
cAAA is this necessary (already in tnpzd)?
            no3_in(:) = tnpzd(:,imobino3)
#   if defined O_mobi_nitrogen_15
            din15_in(:) = tnpzd(:,imobidin15)
#   endif
#  endif
            sgb_in(:) = sg_bathy(i,j,:)
!  end input
            call mobi_driver(
     &                 kmt(i,jrow), c2dtts, rctheta, dayfrac, swr, tnpzd
     &,                t_in
#  if defined O_mobi_o2
     &,                o2_in
#  endif
#  if defined O_carbon_13 || defined O_mobi_caco3
     &,                s_in, dic_in, alk_in, co2_in
#  endif
#  if defined O_carbon_13
c     &,                dic13_in
#  endif
#  if defined O_mobi_silicon
c     &,                sil_in
#  endif
     &,                sgb_in
#  if defined O_mobi_nitrogen
c     &,                no3_in
#   if defined O_mobi_nitrogen_15
c     &,                din15_in
#   endif
#  endif
#  if defined O_sed
     &,                sedcorgflx, sedcalflx
#  endif
     &,                src(i,:,j,:))

!juan added dust iron here instead of in the set source/sinks part:
#  if defined O_mobi_iron
! comes in mol/m2/s gets converted to mmol/m3/s
                src(i,1,j,isdfe) = src(i,1,j,isdfe) +    
     &          sbc(i,j,idfeadep)*1000./(dzt(1)/100.)  !dzt is in cm
!Add hydrothermal Fe:
               do k=1,kmt(i,j)
                src(i,k,j,isdfe) = src(i,k,j,isdfe) +
     &          fe_hydr(i,j,k)
               enddo
#  endif
# if defined O_mobi_silicon
CSPKKK      sink due to opal production term has already been applied to src 
C           in S/R mobi_driver; we add hydrothermal fluxes here
cAAA commented out for now
c			src(i,k,j,issil) = src(i,k,j,issil) 
c     &                         + si_hydr(i,j,k)*si_h_sol !umol Si cm-3 s-1
CSPKKK     accumulate hydrothermal flux and opal burial terms
 	        silhydrflux = silhydrflux 
     &                          + si_hydr(i,j,k)*si_h_sol
     &                            *dzt(k)*dxt(i)*dyt(j)*cst(j) !umol Si/s
 	        opalburial = opalburial + bsi(k)*dxt(i)*dyt(j)*cst(j)
# endif

#  if defined O_sed
            if (addflxo .and. eots) then
               sbc(i,jrow,irorg) = sbc(i,jrow,irorg) + sedcorgflx
               sbc(i,jrow,ircal) = sbc(i,jrow,ircal) + sedcalflx
            endif
#  endif
#  if (defined O_time_averages || defined O_TMM) && defined O_save_mobi_fluxes

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------
            if (timavgperts .and. .not. euler2) then
             expo = rprca
#   if !defined O_mobi_caco3
             ta_rprocal(i,jrow) =  ta_rprocal(i,jrow) + rprca
#   endif
             do k = 1, kmt(i,jrow)
              ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(k)
              ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(k)
              ta_rgraz_Z(i,k,jrow) = ta_rgraz_Z(i,k,jrow) + rgraz_Z(k)
              ta_rgraz_Det(i,k,jrow) = ta_rgraz_Det(i,k,jrow)
     &                                 + rgraz_Det(k)
              ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(k)
              ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(k)
              ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(k)
              ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(k)
#   if defined O_mobi_nitrogen
              ta_rnpp_dop(i,k,jrow) = ta_rnpp_dop(i,k,jrow) 
     &             + rnpp_dop(k)              
              ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow) + rnpp_D(k)
              ta_rnpp_D_dop(i,k,jrow) = ta_rnpp_D_dop(i,k,jrow) 
     &             + rnpp_D_dop(k)
              ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow) + rgraz_D(k)
              ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow) + rmorp_D(k)
              ta_rmorpt_D(i,k,jrow) = ta_rmorpt_D(i,k,jrow)
     &                               + rmorpt_D(k)
              ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(k)
              ta_rbdeni(i,k,jrow) = ta_rbdeni(i,k,jrow) + rbdeni(k)
#    if defined O_mobi_caco3
              ta_rnpp_C_dop(i,k,jrow) = ta_rnpp_C_dop(i,k,jrow) 
     &             + rnpp_C_dop(k)   
#    endif              
#    if defined O_mobi_silicon
              ta_rnpp_Diat_dop(i,k,jrow) = ta_rnpp_Diat_dop(i,k,jrow) 
     &             + rnpp_Diat_dop(k)
#    endif              
#   endif
#   if defined O_kk_ballast
              ta_rgraz_Det_B(i,k,jrow) = ta_rgraz_Det_B(i,k,jrow)
     &                                 + rgraz_Det_B(k)
#   endif
#   if defined O_mobi_caco3
!       print*,'KK tracer', rcalatt(k), rcalpro(k)
              ta_rcalatt(i,k,jrow) = ta_rcalatt(i,k,jrow)
     &                               + rcalatt(k)
              ta_rprocal(i,k,jrow) = ta_rprocal(i,k,jrow)
     &                               + rcalpro(k)
              ta_rnpp_C(i,k,jrow) = ta_rnpp_C(i,k,jrow)
     &                              + rnpp_C(k)
              ta_rgraz_C(i,k,jrow) = ta_rgraz_C(i,k,jrow)
     &                               + rgraz_C(k)
              ta_rmorp_C(i,k,jrow) = ta_rmorp_C(i,k,jrow)
     &                               + rmorp_C(k)
              ta_rmorpt_C(i,k,jrow)= ta_rmorpt_C(i,k,jrow)
     &                               + rmorpt_C(k)
#   endif
#   if defined O_mobi_silicon
              ta_rnpp_Diat(i,k,jrow) = ta_rnpp_Diat(i,k,jrow)
     &                                 + rnpp_Diat(k)
              ta_rgraz_Diat(i,k,jrow) = ta_rgraz_Diat(i,k,jrow)
     &                                  + rgraz_Diat(k)
              ta_rmorp_Diat(i,k,jrow) = ta_rmorp_Diat(i,k,jrow)
     &                                  + rmorp_Diat(k)
              ta_rmorpt_Diat(i,k,jrow)= ta_rmorpt_Diat(i,k,jrow)
     &                                  + rmorpt_Diat(k)
#   endif
#   if defined O_mobi_iron
#    if defined O_save_mobi_diagnostics
              ta_rfeorgads(i,k,jrow) = ta_rfeorgads(i,k,jrow)
     &                              + rfeorgads(k)
              ta_rdeffe(i,k,jrow) = ta_rdeffe(i,k,jrow)
     &                              + rdeffe(k)
              ta_rfeprime(i,k,jrow) = ta_rfeprime(i,k,jrow)
     &                              + rfeprime(k)
              ta_rfesed(i,k,jrow) = ta_rfesed(i,k,jrow)
     &                              + rfesed(k)
              ta_rbfe(i,k,jrow) = ta_rbfe(i,k,jrow)
     &                              + rbfe(k)
              ta_rfecol(i,k,jrow) = ta_rfecol(i,k,jrow)
     &                              + rfecol(k)
#     if defined O_mobi_caco3
              ta_rdeffe_C(i,k,jrow) = ta_rdeffe_C(i,k,jrow)
     &                              + rdeffe_C(k)
#     endif                
#     if defined O_mobi_silicon
              ta_rdeffe_Diat(i,k,jrow) = ta_rdeffe_Diat(i,k,jrow)
     &                              + rdeffe_Diat(k)
#     endif                
#    endif
#   endif
#   if defined O_save_mobi_diagnostics
              ta_ravej(i,k,jrow) = ta_ravej(i,k,jrow) + ravej(k)
              ta_ravej_D(i,k,jrow) = ta_ravej_D(i,k,jrow) + ravej_D(k)
              ta_rgmax(i,k,jrow) = ta_rgmax(i,k,jrow) + rgmax(k)
              ta_rno3P(i,k,jrow) = ta_rno3P(i,k,jrow) + rno3P(k)
              ta_rpo4P(i,k,jrow) = ta_rpo4P(i,k,jrow) + rpo4P(k)
              ta_rpo4_D(i,k,jrow) = ta_rpo4_D(i,k,jrow) + rpo4_D(k)
#   endif
#   if defined O_mobi_caco3
              ta_rdissl(i,k,jrow) = ta_rdissl(i,k,jrow) + rdissl(k)
              ta_romega_c(i,k,jrow) = ta_romega_c(i,k,jrow)
     &                              + romega_c(k)
              ta_rexpocaco3(i,k,jrow) =
     &          ta_rexpocaco3(i,k,jrow) + rexpocaco3(k)
#   endif
#   if defined O_kk_ballast
              ta_rremi_B(i,k,jrow) = ta_rremi_B(i,k,jrow)
     &                                 + rremi_B(k)
              ta_rexpo_B(i,k,jrow) = ta_rexpo_B(i,k,jrow)
     &                                 + rexpo_B(k) !ballast detritus
#   endif
#   if defined O_mobi_silicon
              ta_rexpoopl(i,k,jrow) = ta_rexpoopl(i,k,jrow) 
     &                                 + rexpoopl(k) 
              ta_rdisopl(i,k,jrow) = ta_rdisopl(i,k,jrow)
     &                                 + rdisopl(k)
              ta_rproopl(i,k,jrow) = ta_rproopl(i,k,jrow)
     &                                 + rproopl(k)
#   endif
              expo = expo*dztr(k)
              ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(k)
              ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(k)
              expo = expo - rprca*rcak(k)
              ta_rexpocal(i,k,jrow) = ta_rexpocal(i,k,jrow) + expo
              expo = expo*dzt(k)
#   if defined O_mobi_nitrogen && defined O_mobi_o2
              ta_rwcdeni(i,k,jrow) = ta_rwcdeni(i,k,jrow) 
     &           + rwcdeni(k)
#   endif
#   if defined O_mobi_iron
                ta_rremife(i,k,jrow) = ta_rremife(i,k,jrow) 
     &               + rremife(k)
                ta_rexpofe(i,k,jrow) = ta_rexpofe(i,k,jrow)
     &               + rexpofe(k)
#   endif 
!
             enddo
            endif
#  endif
          endif
        enddo
      enddo

#  if defined O_embm
#   if defined O_mobi_silicon
#    if defined O_kk_si_compensating_sources
CSPKKK This will not work with TMM. opalburial and silhydrflux as computed above 
C      are the spatially integrated local (for the current CPU) terms. But we need 
C      the global ones below. How to get them here from the TMM? We could use values 
C      from the previous time step? Or call MPI_sum via a C wrapper? Also, opalburial 
C      and silhydrflux need to be passed back to the TMM to compute the global sum for 
C      the next time step. Ugh.
	  rivsil = (opalburial - silhydrflux - sildustflux)/globaldisch !scale Si over total discharge
#    else
	  rivsil = globalsilwflx/globaldisch
#    endif
#    if defined O_save_mobi_fluxes
	 do j=js,je
	   jrow = j + joff
	   do i=is,ie
		 if (kmt(i,jrow) .gt. 0) then      
cAAA comment out for now
c		 src(i,1,j,issil) = src(i,1,j,issil)
c     &                        + rivsil*disch(i,j)*dztr(1)
		  if (timavgperts .and. .not. euler2) then
			 ta_rrivsil(i,jrow) = ta_rrivsil(i,jrow)
     &                           + rivsil*disch(i,j)
		  endif
		 endif
	  enddo
	 enddo 
#    endif
#   endif
#  endif ! O_embm
# endif
# if defined O_carbon_fnpzd

!-----------------------------------------------------------------------
!    option to write and read fixed biological fluxes for carbon and
!    alkalinity. code assumes that restarts line up with the first
!    record. the repeating time period and restarts should probably
!    line up with the start of the year.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff

!       adjust counters and read fluxes, if first row of slab
        if (jrow .eq. 2) then
          mfnpzd = mfnpzd + 1
          if (mfnpzd .gt. mxfnpzd) mfnpzd = 1
          if (nfnpzd .le. mxfnpzd) nfnpzd = nfnpzd + 1
!         read fluxes, if enough records written
          if (nfnpzd .gt. mxfnpzd) then
            fname = new_file_name ("O_fnpzd.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              ib(:) = 1
              ic(:) = 1
              ic(1) = imtm2
              ic(2) = jmtm2
              ic(3) = km
              ib(4) = mfnpzd
              print*, "Reading fnpzd from: ", trim(fname),
     &          " record:", mfnpzd
              call openfile (fname, iou)
#  if defined O_carbon
              call getvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,1)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
#  if defined O_mobi_alk
              call getvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,2)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
            endif
          endif
        endif

        if (nfnpzd .le. mxfnpzd) then

!         store fluxes, if enough records written
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              fnpzd(i,jrow,k,1) = src(i,k,j,isdic)
#  endif
#  if defined O_mobi_alk
              fnpzd(i,jrow,k,2) = src(i,k,j,isalk)
#  endif
            enddo
          enddo
!         write fluxes, if last row of slab
          if (jrow .eq. jmt-1) then
            fname = new_file_name ("O_fnpzd.nc")
            call openfile (fname, iou)
            print*, "Writing fnpzd to: ", trim(fname),
     &          " record:", mfnpzd
            ib(:) = 1
            ic(:) = 1
            ic(1) = imtm2
            ic(2) = jmtm2
            ic(3) = km
            ib(4) = mfnpzd
            time = relyr
            call putvars ('time', iou, mfnpzd, time, c1, c0)
#  if defined O_carbon
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,1)
            call putvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
#  if defined O_mobi_alk
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,2)
            call putvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
          endif

        else

!         set fluxes to source terms
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              src(i,k,j,isdic) = fnpzd(i,jrow,k,1)
#  endif
#  if defined O_mobi_alk
              src(i,k,j,isalk) = fnpzd(i,jrow,k,2)
#  endif
            enddo
          enddo

        endif

      enddo
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     set source for c14
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
#  if defined O_mobi
              src(i,k,j,isc14) = src(i,k,j,isdic)*rc14std
     &                         - 3.836e-12*t(i,k,j,ic14,taum1)
#  else
              src(i,k,j,isc14) = - 3.836e-12*t(i,k,j,ic14,taum1)
#  endif
            enddo
          endif
        enddo
      enddo
# endif
# if defined O_sed && !defined O_sed_uncoupled

!-----------------------------------------------------------------------
!     set source terms from sediment model
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            k = kmt(i,jrow)
#  if defined O_carbon && defined O_mobi
            src(i,k,j,isdic) = src(i,k,j,isdic) + sbc(i,j,ibdicfx)
#   if defined O_global_sums
            if (addflxo .and. eots) dtoic = dtoic - sbc(i,j,ibdicfx)
     &        *c2dtts*dtxcel(k)*dxt(i)*dyt(jrow)*cst(jrow)*dzt(k)
#   endif
#  endif
#  if defined O_mobi_alk
            src(i,k,j,isalk) = src(i,k,j,isalk) + sbc(i,j,ibalkfx)
#  endif
          endif
        enddo
      enddo
# endif

#ifndef O_TMM
!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt
# if defined O_matrix
        isTile=.false.
        if ((n.ge.ntilestart).and.(n.le.ntileend)) then
          isTile=.true.
		  ntile=n-ntilestart+1  ! index for tiles
!         reset tracer to tiles/initial condition
		  t(:,:,:,n,tau)=tile(:,:,:,ntile)
		  t(:,:,:,n,taum1)=tile(:,:,:,ntile)         
        endif
# endif

!-----------------------------------------------------------------------
!       calculate advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

# if defined O_consthmix
#  if !defined O_biharmonic || defined O_bryan_lewis_horizontal

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
#   if defined O_bryan_lewis_horizontal
     &                         diff_cet(k)*cstdxur(i,j)*
#   else
     &                         ah_cstdxur(i,j)*
#   endif
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#   if defined O_isopycmix

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
#    if defined O_bryan_lewis_horizontal
     &                         diff_cnt(k)*
#    else
     &                         diff_cnt*
#    endif
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#   else

!       diffusive flux on northern face of "T" cells is not calculated.
!       instead, it is incorporated into DIFF_Ty for performance issues.

#   endif
#  else

!-----------------------------------------------------------------------
!       diffusive flux across eastern face of "T" cell
!       diffusive flux across northern face of "T" cell
!-----------------------------------------------------------------------

        m = 2 + n
        do j=js,je
          jrow = j + joff
          ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (del2(i+1,k,j,m)-del2(i,k,j,m))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                        (del2(i,k,j+1,m) - del2(i,k,j,m))
            enddo
          enddo
        enddo
#  endif
# else
#  if defined O_smagnlmix

!       diffusive flux on eastern and northern faces of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) = diff_cet(i,k,j)*cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) = diff_cnt(i,k,j)*csu_dyur(jrow)*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
#  endif
# endif
!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

# if defined O_isopycmix

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)
# endif

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
# if defined O_replacst
            diff_fb(i,0,j)  = c0
# else
            diff_fb(i,0,j)  = stf(i,j,n)
# endif
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

# if defined O_source_term || defined O_mobi || defined O_carbon_14

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

#  if defined O_mobi || defined O_carbon_14
        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif
c        print*,'src(alk)'
c,src(50,5,50,itrc(ipo4))
#  endif
#  if defined O_shortwave

!-----------------------------------------------------------------------
!       incorporate short wave penetration into source
!-----------------------------------------------------------------------

        if (n .eq. 1) then
          call swflux0 (joff, js, je, istrt, iend, source)
        endif
#  endif
# endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
# if defined O_isopycmix && defined O_gent_mcwilliams && !defined O_fct && !defined O_quicker
     &          - ADV_Txiso(i,k,j,n) - ADV_Tyiso(i,k,j,jrow,n)
     &          - ADV_Tziso(i,k,j)
# endif
# if defined O_source_term || defined O_mobi || defined O_carbon_14
     &          + source(i,k,j)
# endif
# if defined O_plume
     &          + work(i,k,j)*subflux(i,jrow,n)
# endif
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo

# if defined O_matrix
!     accumulate explicit matrix and reset to tiles/initial condition
	  if (isTile) then
	    call MATRIX_STORE_EXPLICIT(ntile,t(1,1,1,n,taup1),twodt)
      endif      
# endif        
# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion

!       2nd: add in portion of vertical diffusion handled implicitly
#  if defined O_matrix
		if (isTile) then
          call ivdift (joff, js, je, istrt, iend, n, dtmatrix)
        else
          call ivdift (joff, js, je, istrt, iend, n, twodt)        
		endif          
#  else
        call ivdift (joff, js, je, istrt, iend, n, twodt)
#  endif

# endif

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

# if defined O_replacst

      gamma = zw(1)/twodt(1)
      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,1) = gamma*(sbc(i,jrow,isst) - t(i,1,j,1,taup1))
          stf(i,j,2) = gamma*(sbc(i,jrow,isss) - t(i,1,j,2,taup1))
          t(i,1,j,1,taup1) = sbc(i,jrow,isst)
          t(i,1,j,2,taup1) = sbc(i,jrow,isss)
        enddo
      enddo
# endif
# if defined O_convect_brine

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!     convect brine rejected under all ice categories
!-----------------------------------------------------------------------

      call convect_brine (joff, js, je, is, ie)
# else
#  if !defined O_implicitvmix || defined O_isopycmix

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

#   if defined O_fullconvect
      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
#   else
      call convct (t(1,1,1,1,taup1), ncon, joff, js, je, istrt, iend
     &,            kmt)
#   endif
#  endif
# endif
# if defined O_save_convection
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif
# endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

# if defined O_matrix
!     accumulate implicit matrix
      do n=ntilestart,ntileend
        ntile=n-ntilestart+1
        call MATRIX_STORE_IMPLICIT(ntile,t(1,1,1,n,taup1))
      enddo
# endif

# if defined O_fourfil || defined O_firfil

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
# endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

# if !defined O_replacst

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------
	  do n=1,2
		if (trsbcindex(n) .ne. 0) then
		  call set_sbc(joff, js, je, istrt, iend, trsbcindex(n), n, 
     &                 taup1, .true.)
		endif
	  enddo
# endif

	  do n=3,nt
		if (trsbcindex(n) .ne. 0) then
		  call set_sbc(joff, js, je, istrt, iend, trsbcindex(n), n, 
     &                 taup1, .true.)
		endif
	  enddo

#endif ! not O_TMM
# if defined O_sed
!-----------------------------------------------------------------------
!     accumulate bottom tracer values for the sediment model
!-----------------------------------------------------------------------

      if (addflxo .and. eots) then
C SPK: this used to be twodt(1) but in offline mode this is not set anywhere
        if (joff .eq. 0) atsed = atsed + c2dtts
        do j=js,je
          jrow = j + joff
          do i=is,ie
            if (kmt(i,jrow) .gt. 0) then
              k = kmt(i,jrow)
              sbc(i,jrow,ibtemp) = sbc(i,jrow,ibtemp)
     &                           + t(i,k,j,itemp,taup1)*c2dtts
              sbc(i,jrow,ibsalt) = sbc(i,jrow,ibsalt)
     &                           + t(i,k,j,isalt,taup1)*c2dtts
#  if defined O_carbon
              sbc(i,jrow,ibdic) = sbc(i,jrow,ibdic)
     &                          + t(i,k,j,idic,taup1)*c2dtts
#  endif
#  if defined O_mobi_alk
              sbc(i,jrow,ibalk) = sbc(i,jrow,ibalk)
     &                          + t(i,k,j,ialk,taup1)*c2dtts
#  endif
#  if defined O_mobi_o2
              sbc(i,jrow,ibo2) = sbc(i,jrow,ibo2)
     &                         + t(i,k,j,io2,taup1)*c2dtts
#  endif
            endif
          enddo
        enddo
      endif
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     calculate diagnostic delta carbon 14
!-----------------------------------------------------------------------

      if (tsiperts .or. timavgperts) then
        rrc14std = 1000./rc14std
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              dc14(i,k,j) = (rrc14std*t(i,k,j,ic14,taup1)
     &                      /(t(i,k,j,idic,taup1) + epsln) - 1000.)
     &                      *tmask(i,k,j)
            enddo
          enddo 
        enddo
      endif
      if (tsiperts .and. eots) then
        if (js+joff .eq. 2) dc14bar = 0.
        do j=js,je
          jrow = j + joff
          fy = cst(jrow)*dyt(jrow)
          do k=1,km
            fyz = fy*dzt(k)
            do i=istrt,iend
              dc14bar = dc14bar + dc14(i,k,j)*dxt(i)*fyz*tmask(i,k,j)
            enddo
          enddo
        enddo
      endif
      if (timavgperts .and. .not. euler2) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              ta_dc14(i,k,jrow) = ta_dc14(i,k,jrow) + dc14(i,k,j)
            enddo
          enddo
        enddo
      endif
# endif

      return
      end

#ifndef O_TMM

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, is, ie
      integer mask, m

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, dtdx, dtdy, dtdz
      real r2dt, cosdyt, fx, darea, boxar, rtwodt, sumdx, delx
      real sumdxr, dxdy, dxdydz

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "coord.h"
#include "cregin.h"
#include "csbc.h"
# if defined O_tracer_averages
#include "ctavg.h"
# endif
# if defined O_mobi_caco3
#include "mobi.h"
# endif      
#include "diag.h"
#include "diaga.h"
#include "emode.h"
#include "grdvar.h"
#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "vmixc.h"

# if defined O_meridional_tracer_budget
#include "ctmb.h"
      real stor(imt,km), div(imt,km), sorc(imt,km), dif(imt,km)
      real t1(imt), t2(imt), t3(imt), t4(imt)
# endif
# if defined O_time_step_monitor
      real temp1(imt,km), temp2(imt,km), temp3(imt,km)
# endif
      real twodt(km)

# if defined O_isopycmix
#include "isopyc.h"
# endif
#include "fdift.h"

# if defined O_save_mixing_coeff

!-----------------------------------------------------------------------
!     diagnostic: estimate mixing coefficients on east, north, and
!                 bottom face of T cells from the flux
!-----------------------------------------------------------------------

      if (cmixts .and. n .eq. 1 .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=2,imt-1
              dtdx = (t(i+1,k,j,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i+1,k,j)*tmask(i,k,j)
     &                *cstr(jrow)*dxur(i) + epsln
              ce(i,k,j,2) = diff_fe(i,k,j)/dtdx
#  if !defined O_consthmix || defined O_biharmonic || defined O_isopycmix
              dtdy = (t(i,k,j+1,1,taum1)-t(i,k,j,1,taum1))
     &                *tmask(i,k,j+1)*tmask(i,k,j)
     &                *dyur(jrow) + epsln
              cn(i,k,j,2) = csur(jrow)*diff_fn(i,k,j)/dtdy
#  else
              cn(i,k,j,2) = csur(jrow)*ah
#  endif
            enddo
          enddo
        enddo
        do j=js,je
          jrow = j + joff
          do k=1,km-1
            do i=2,imt-1
              dtdz = (t(i,k,j,1,taum1)-t(i,k+1,j,1,taum1))
     &                *tmask(i,k,j)*tmask(i,k+1,j)
     &                *dzwr(k) + epsln
              cb(i,k,j,2) = diff_fb(i,k,j)/dtdz
#  if defined O_isopycmix
     &                      + diff_fbiso(i,k,j)/dtdz
#  endif
            enddo
          enddo
          do i=2,imt-1
            cb(i,km,j,2) = 0.0
          enddo
        enddo

        do j=js,je
          call setbcx (ce(1,1,j,2), imt, km)
          call setbcx (cn(1,1,j,2), imt, km)
          call setbcx (cb(1,1,j,2), imt, km)
        enddo
      endif
# endif

# if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: initialize temperature before convection
!-----------------------------------------------------------------------

      if (exconvts .and. n .eq. 1 .and. eots) then
        do j=js,je
          do k=1,km
            do i=1,imt
              excnv0(i,k,j) = t(i,k,j,1,taup1)
            enddo
          enddo
        enddo
      endif

# endif

# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif
# endif

# if defined O_tracer_averages

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif
# endif

# if defined O_tracer_yz

!----------------------------------------------------------------------
!     diagnostic: integrate tracer equations in longitude
!----------------------------------------------------------------------

      if (tyzts .and. eots) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            rtwodt = c1/twodt(k)
            sumdx = c0
            do m=1,5
              tyz(jrow,k,n,m) = c0
            enddo
            do i=is,ie
              delx  = dxt(i)*tmask(i,k,j)
              sumdx = sumdx + delx
              tyz(jrow,k,n,1) = tyz(jrow,k,n,1)+ delx*t(i,k,j,n,tau)
              tyz(jrow,k,n,2) = tyz(jrow,k,n,2) + rtwodt*delx*
     &          (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              tyz(jrow,k,n,3) =  tyz(jrow,k,n,3)
     &          - delx*(ADV_Ty(i,k,j,jrow,n) + ADV_Tz(i,k,j))
              tyz(jrow,k,n,4) =  tyz(jrow,k,n,4)
     &          + delx*(DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j))
#  if defined O_source_term || defined O_mobi || defined O_carbon_14
              tyz(jrow,k,n,5) = tyz(jrow,k,n,5) + delx*source(i,k,j)
#  endif
            enddo
            sumdxr = c1/(sumdx + epsln)
            do m=1,5
              tyz(jrow,k,n,m) = tyz(jrow,k,n,m)*sumdxr
            enddo
          enddo
        enddo
        if (js+joff .eq. 2) then
          do m=1,5
            do k=1,km
              tyz(1,k,n,m)   = c0
              tyz(jmt,k,n,m) = c0
            enddo
          enddo
        endif
      endif
# endif

# if defined O_meridional_tracer_budget

!----------------------------------------------------------------------
!     diagnostic: integrate equations in depth, lon, and time for
!                 each basin and jrow. basin #  0 is used to catch
!                 values in land areas
!----------------------------------------------------------------------

      if (tmbperts .and. eots) then
        do j=js,je
          jrow = j + joff
          if (jrow .eq. 2 .and. n .eq. 1) numtmb = numtmb + 1
          cosdyt = cst(jrow)*dyt(jrow)
          do k=1,km
            do i=is,ie
              stor(i,k) = 0.0
              div(i,k)  = 0.0
              sorc(i,k) = 0.0
              dif(i,k)  = 0.0
            enddo
          enddo
          do k=1,km
            rtwodt = c1/twodt(k)
            do i=is,ie
              dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
              dxdydz          = dxdy*dzt(k)
              stor(i,k) = stor(i,k) + rtwodt*dxdydz*
     &                            (t(i,k,j,n,taup1) - t(i,k,j,n,taum1))
              div(i,k)  = div(i,k)  - dxdydz*ADV_Ty(i,k,j,jrow,n)
#  if defined O_source_term || defined O_mobi || defined O_carbon_14
              sorc(i,k) = sorc(i,k) + dxdydz*source(i,k,j)
#  endif
              dif(i,k)  = dif(i,k)  + dxdydz*DIFF_Ty(i,k,j,jrow,n)
            enddo
          enddo

!         the accumulation is done this way for issues of speed

          do i=is,ie
            t1(i) = stor(i,1)
            t2(i) = div(i,1)
            t3(i) = sorc(i,1)
            t4(i) = dif(i,1)
          enddo

          do k=2,km
            do i=is,ie
              t1(i) = t1(i) + stor(i,k)
              t2(i) = t2(i)  + div(i,k)
              t3(i) = t3(i) + sorc(i,k)
              t4(i) = t4(i)  + dif(i,k)
            enddo
          enddo

          do i=is,ie
            m               = msktmb(i,jrow)
            tstor(jrow,n,m) = tstor(jrow,n,m) + t1(i)
            tdiv(jrow,n,m)  = tdiv(jrow,n,m)  + t2(i)
            tsorc(jrow,n,m) = tsorc(jrow,n,m) + t3(i)
            tdif(jrow,n,m)  = tdif(jrow,n,m)  + t4(i)
          enddo

          k = 1
          do i=is,ie
            m               = msktmb(i,jrow)
            dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
            tflux(jrow,n,m) = tflux(jrow,n,m) + dxdy*stf(i,j,n)
          enddo
          if (n .eq. 1) then
            do k=1,km
              do i=is,ie
                m               = msktmb(i,jrow)
                dxdy            = cosdyt*dxt(i)*tmask(i,k,j)
                dxdydz          = dxdy*dzt(k)
                smdvol(jrow,m)  = smdvol(jrow,m)  + dxdydz
              enddo
            enddo
          endif
        enddo
      endif
# endif

# if defined O_gyre_components

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)
# endif

# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt1 (joff, js, je, n)
# endif
# if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate r.h.s. terms in the tracer equations
!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt1 (joff, js, je, n)
# endif

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

      implicit none

      integer idiag, j, js, je, k, i, joff, iocv, jrow, is, ie

      real rdt, reltim, period

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "diaga.h"
#include "iounit.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "timeavgs.h"

# if defined O_save_convection_full

!-----------------------------------------------------------------------
!     diagnostic: save tracer time change due to explicit convection
!     idiag = 10 signifies diagnostics immediately after convection
!-----------------------------------------------------------------------

      if (exconvts .and. idiag .eq. 10 .and. eots) then
        rdt = c1/c2dtts

!       excnv1 = epsln over land cells and d(convect)/dt over ocean

        do j=js,je
          do k=1,km
            do i=1,imt
              excnv1(i,k,j) = tmask(i,k,j)*
     &         (t(i,k,j,1,taup1)-excnv0(i,k,j))*rdt
     &                       + (c1-tmask(i,k,j))*epsln
            enddo
          enddo
        enddo

        reltim = relyr
        if (joff + js .eq. 2) then
          write (stdout,*) ' =>Writing explicit convection at ts=',itt
     &  , ' ',stamp
          call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

          period = 0.0
          iotext = 'read(iocv) reltim, period, imt, jmt, km, flag'
          write (iocv) stamp, iotext, expnam
          write (iocv) reltim, period, imt, jmt, km, epsln

          iotext = 'read(iocv) (xt(i),i=1,imt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, xt, imt)

          iotext = 'read(iocv) (yt(j),j=1,jmt)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, yt, jmt)

          iotext = 'read(iocv) (zt(k),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, zt, km)

          call relunit (iocv)
        endif

        call getunit (iocv, 'cvct.dta'
     &,                'unformatted sequential append ieee')

        do j=js,je
          jrow = j+joff

          write(iotext,'(a10,i4)') ' for jrow=',jrow
          iotext(15:) = ': read (iocv) (convect(i,k),i=1,imt),k=1,km)'
          write (iocv) stamp, iotext, expnam
          call wrufio (iocv, excnv1(1,1,j), imt*km)

        enddo
        call relunit(iocv)
      endif
# endif

# if defined O_term_balances

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)
# endif

# if defined O_xbts

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering
!-----------------------------------------------------------------------

      if (xbtperts .and. eots) call txbt2 (joff, js, je, idiag)
# endif
# if defined O_mom_tbt

!-----------------------------------------------------------------------
!     diagnostic: accumulate d/dt(tracer) after convection
!                 and filtering

!-----------------------------------------------------------------------

      if (tbtperts .and. eots) call tbt2 (joff, js, je, idiag)
# endif

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie

      real rts

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

      subroutine ivdift (joff, js, je, is, ie, n, twodt)

# if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real rc2dt

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "levind.h"
#include "mw.h"
#include "switch.h"
#include "vmixc.h"

      real twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

#  if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#  else
#   if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif
#   endif
#  endif

      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)

!     compute residual implicit vertical mixing

#  if defined O_xbts || defined O_mom_tbt
      if (eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#  else
#   if defined O_term_balances
      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif
#   endif
#  endif
# endif

      return
      end

      subroutine swflux0 (joff, js, je, is, ie, source)

# if defined O_mom && defined O_shortwave
!=======================================================================
!       incorporate short wave penetration via the "source"
!       term. note that the divergence of shortwave for the first
!       level "divpen(1)" is compensating for the effect of having
!       the shortwave component already included in the total
!       surface tracer flux "stf(i,j,1)"

!       incorporating the penetrative shortwave radiative flux into
!       the vertical diffusive flux term directly is not correct when
!       vertical diffusion is solved implicitly.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

!     output:
!       source = source term penetrative short wave
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, k, i, is, ie

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "cshort.h"

      real source(imt,km,jsmw:jemw)

      if (iswr .ne. 0) then
        do j=js,je
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              source(i,k,j) = source(i,k,j)
     &                       + sbc(i,jrow,ipsw)*divpen(k)
            enddo
          enddo
        enddo
      endif
# endif

      return
      end
!

#endif ! SPK not O_TMM

#endif
