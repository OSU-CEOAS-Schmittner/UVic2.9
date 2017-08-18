      subroutine dh (sal1d, Fx, Fmix, Fneti, ultnt, condb, Tw, Ts, Ti
     &,              Tbot, Tf, qi, hi, hs, dhib, subi, subs, dhit, dhs
     &,              dhif, dhsf, ni)

#if defined O_ice && defined O_ice_cpts && defined O_embm

! Energy-conserving sea ice model
! Routines to grow/melt ice and adjust temperature profile

! See Bitz, C.M., and W.H. Lipscomb, 1999:
! A new energy-conserving sea ice model for climate study,
! J. Geophys. Res., in press.

! The author grants permission to the public to copy and use this
! software without charge, provided that this Notice and any statement
! of authorship are reproduced on all copies.

!     Computes the thickness changes at the top and bottom
!     and adjusts layer specific enthalpy
!     This algorithm does not allow h<0

!     Fx= actual flux of heat from the mixed layer under sea ice
!     (equal to Fmix unless all the ice melts away)
!     parameterized according to the
!     Tw as described in papers by McPhee/Steel
!     compensates for rare case of melting entire slab through

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"

! input variables
      real sal1d(NMAX) ! ice salinity (ppt)
      real Ts          ! top srf temp (C)
      real Ti(0:NMAX)  ! snow/ice internal temp (C)
      real Tbot        ! ice bottom in (C)
      real Tw          ! water temperature below ice (C)
      real Tf          ! freezing temperature of sea water
      real hi          ! initial ice thickness (m)
      real hs          ! initial snow thickness (m)
      real Fmix        ! Flux from mixed layer from relaxation (W/m**2)
      real Fneti       ! net flux at upper surface (W/m**2)
      real condb       ! conductive flux at bottom (w/m**2)
      real ultnt       ! latent heat flux, for sublimation (w/m**2)
      integer ni       ! number of layers

! output variables
!     thickness changes from grow/melt (default) or sublimate/flooding
      real dhib         ! ice bot, dhib<0 if melt (m)
      real dhit         ! ice top, dhit<=0 (m)
      real dhs          ! snow top, dhit<=0 (m)
      real subi         ! ice top, subi<0 if sublimating (m)
      real subs         ! snow, subs<0 if sublimating (m)
      real dhif         ! ice top from flooding, dhif>0 (m)
      real dhsf         ! snow from flooding, dhsf>0 (m)

      real qi(NMAX)     ! energy of melt of ice per unit vol. (J/m**3)
      real Fx           ! actual flx of heat used from ocn (w/m**2)

! local variables
      real delti(NMAX)  ! evolving ice layer thickness
      real delts        ! evolving snow thickness
      real sumr         ! dummy for adjusting delti
      real rmvi         ! another dummy for adjusting delti
      real dtop         ! dhit+subi for adjust

      real qigrow       ! energy of melt of ice that grows
      real Ebot         ! heat available to grow/melt at bottom
      real Etop         ! heat available to melt at top

      real qiflood      ! energy of melt of flooded ice (W/m**2)
      real zintfc       ! height of snow/ice interf. wrt ocn (m)
      real qs           ! energy of melting of snow per unit volume
      real Enet         ! sum of the energy of melting relative to melt
      real Evnet        ! sum of the energy of melting relative to vapour
      real energ        ! function
      integer layer
      logical alarm2,alarm

      alarm2=.false.
      alarm=.false.

      dhib  = 0.
      dhit  = 0.
      dhs   = 0.
      dhif  = 0.
      dhsf  = 0.
      subi  = 0.
      subs  = 0.

! thickness of snow and each ice layer
      delti(1) = hi/ni
      do layer = 2,ni
        delti(layer) = delti(1)
      enddo
      delts = hs

! energy of melt per unit vol for snow and ice for each layer and sum
      qs = -RFLSNO
      enet = 0.
      do layer = 1,ni
        qi(layer) = energ(Ti(layer),sal1d(layer))
        enet = enet + qi(layer)
      enddo
      enet = enet*delti(1) + qs*hs

      if (ultnt.gt.0.0) then
! etop is positive is sublimating and negative if condensing
         Etop  = ultnt*dtau
         Evnet = -RVLSNO*hs-RVLICE*hi+Enet+Etop
         if ( Evnet .ge. 0.0 ) then
            subi=-hi
            subs=-hs
            Fx=condb-Evnet/dtau  ! maybe wrong, irrelevant though
# if defined O_ice_cpts_warnings
            print*,'sublimated though all snow and ice',Fx,hi,hs
            stop 'sublimate too much is not allowed, model will stop'
!            return
# endif
         endif
         call srfsub(    qi,   qs, delti, delts,    ni,
     $                 subi, subs,  etop,  enet,  alarm2 )

! adjust the layer thickness to reflect subl/cond changes
         delts = hs + subs
         rmvi = subi
         do layer = 1,ni
           sumr = max( -delti(layer), rmvi )
           rmvi = rmvi-sumr
           delti(layer) = delti(layer) + sumr
         enddo
      endif

      if ( Fneti .gt. 0.0 ) then
!     it is possible to melt at srf when Ts is below melting
!     but it never melts more than a negligible amount
         Etop=Fneti*dtau
         Enet=Enet+Etop
         if (Enet.ge.0.0) then
            dhit=-(hi+subi)
            dhs=-(hs+subs)
            Fx=condb-Enet/dtau
            return
         endif
         call srfmelt (qi, qs, delti, delts, ni,
     $                 dhit, dhs, etop, alarm2)
!       adjust the layer thickness to reflect melt changes
        delts = delts + dhs
        rmvi = dhit
        do layer = 1,ni
          sumr = max( -delti(layer), rmvi )
          rmvi = rmvi - sumr
          delti(layer) = delti(layer) + sumr
        enddo
      endif

      dtop=dhit+subi

      Fx=Fmix
      Ebot=dtau*(Fmix-condb)

      if (Ebot.le.0.0) then
         qigrow=energ(Tbot,SALNEW)
         dhib=Ebot/qigrow
      else
!     melt at bottom, melting temperature = Tbot
         qigrow=0.              ! on purpose
         if ((Enet+Ebot).ge.0.0) then
            dhib=-(hi+dtop)
            dhs=-(hs+subs)
            Fx=condb-Enet/dtau
            return
         endif
         call botmelt (qi, qs, delti, delts, ni,
     $                 dhib, dhs, ebot, alarm2)
      endif

! print out helpful variables if something went wrong
      if (alarm2) then
         print*,'Ti ',(Ti(layer),layer=1,ni)
         print*,'Tw ',Tw,condb,Fx,Fneti
         print*,'Enet, Ebot',Enet,Ebot,Fneti*dtau,dtau*(Fmix-condb)
         Enet=0.
         do layer=1,ni
           Enet=Enet+qi(layer)
         enddo
         Enet=Enet*hi/ni+qs*hs
         print*,'orig Enet',Enet
         print*,'delti',delts, (delti(layer),layer=1,ni)
         print*,dhib,dhit,subs,subi,dhs
         alarm=.true.
      endif

! flooding
      zintfc = hi - (rhosno*hs +  rhoice*hi)/rhoocn
      qiflood = 0.
      if ((zintfc .lt. 0.) .and. (hs+dhs.gt.0.)) then
         dhsf   = -rhoice/rhosno*zintfc
         dhsf   = min(dhsf,hs+dhs)
         dhif   = rhosno/rhoice*dhsf
         if (dhif .gt. 0) qiflood = qs*dhsf/dhif
      endif

      call adjust(hi,dhib,dtop,dhif,dhsf,qiflood,qigrow,ni,qi)

      return
      end

      subroutine srfsub( qi   ,qs   ,delti ,delts ,ni    ,
     $                   subi ,subs ,etop  ,enet  ,alarm )

!     compute the sea ice and snow thickness changes from
!     sublimation/condensation
! Etop>0     qi<0  qs<=0 and   Enet<0

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"

! input variables
      real qi (1:nmax), qs    ! energy of melt per vol  (J/m**3)
      real delti(NMAX), delts ! thickness of ice/snow layer (m)
      integer ni              ! number of layers

! output variables
      real subi, subs ! subl/cond. amount for ice/snow (m)

! input/output variables
      real etop ! energy avail to sub/cond ice/snow   (J/m**2)
      real enet ! energy needed to melt all ice/snow  (J/m**2)

! local variables
      real u     ! energy of melt + vapour for ice/snow (J/m**3)
      real subil ! subl from this layer
      integer layer
      logical alarm

      if ( delts .gt. 0. ) then
!       convert etop into snow thickness
!       positive if cond/negative if subl
        u = qs - rvlsno
        subs = etop/u
        if ( (delts + subs) .ge. 0. ) then
!         either all condensation becomes snow
!         or sublimate some snow but no ice
          etop = 0.0
          enet = enet + subs * qs
          return
        else
!         sublimate all of the snow and some ice too
          subs = -delts
          etop = etop + delts * u
          enet = enet + subs  * qs
        endif
      endif

      do layer = 1,ni
!       convert etop into ice thickness
!       positive if cond/negative if subl
        u = qi(layer) - rvlice
        subil = etop/u
        if ( (delti(layer)+subil) .ge. 0. ) then
!         either all condensation becomes ice
!         or sublimate some ice from this layer
          subi = subi + subil
          etop = 0.0
          enet = enet + subil * qi(layer)
          return
        else
!         sublimate all of the layer
          subi = subi - delti(layer)
          etop = etop + delti(layer) * u
          enet = enet - delti(layer) * qi(layer)
        endif
      enddo

! model should never get here, probably a only minor error if it does
      print*, 'ERROR in srfsub',Etop
      alarm=.true.

      return
      end

      subroutine srfmelt(   qi,   qs, delti, delts, ni,
     $                    dhit,  dhs,  etop, alarm  )

! melt ice/snow from the top srf
! Etop>0     qi<0  and  qs<=0

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"

! input variables
      real qi (1:nmax), qs    ! energy of melt per vol (J/m**3)
      real delti(NMAX), delts ! thickness of ice/snow layer (m)
      integer ni              ! number of layers

! output variables
      real dhit, dhs ! ice/snow thickness change(m)

! input/output variables
      real etop ! energy avail to melt ice and snow (J/m**2)

! local variables
      real dhitl ! melt from this layer (m)
      integer layer
      logical alarm

      if ( delts .gt. 0. ) then
!       convert etop into snow thickness
        dhs  =  etop/qs
        if ( (delts + dhs) .ge. 0. ) then
!         melt only some of the snow
          etop = 0.0
          return
        else
!         melt all of the snow and some ice too
          dhs  = -delts
          etop = etop+qs*delts
        endif
      endif

      do layer = 1,ni
!       convert etop into ice thickness
        dhitl = etop/qi(layer)
        if ( (delti(layer)+dhitl) .ge. 0. ) then
!         melt some ice from this layer
          dhit = dhit + dhitl
          etop = 0.0
          return
        else
!         melt all of the ice in this layer
          dhit = dhit - delti(layer)
          etop = etop + delti(layer) * qi(layer)
        endif
      enddo

! model should never get here, probably a only minor error if it does
      print*, 'ERROR in srfmelt',Etop
      alarm=.true.
      return
      end

      subroutine botmelt(   qi,   qs, delti, delts, ni,
     $                    dhib,  dhs,  ebot, alarm  )

! melt from bottom
! Ebot>0     qi<0  and  qs<=0

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"

! input variables
      real qi (1:nmax), qs    ! energy of melt per vol (J/m**3)
      real delti(NMAX), delts ! thickness of ice/snow layer (m)
      integer ni              ! number of layers

! output variables
      real dhib ! ice thickness change (m)

! input/output variables
      real ebot ! energy avail to melt ice and snow (J/m**2)
      real dhs  ! snow thickness change (m)

! local variables
      real dhibl ! melt from this ice layer (m)
      real dhsl  ! melt from this snow layer (m)
      integer layer
      logical alarm

      do layer = ni,1,-1
         dhibl = ebot/qi(layer)
         if ( (delti(layer)+dhibl) .ge. 0. ) then
            dhib = dhib + dhibl
            ebot = 0.0
            return
         else
            dhib = dhib - delti(layer)
            ebot = ebot + delti(layer) * qi(layer)
         endif
      enddo

!     finally melt snow if necessary
      dhsl = ebot/qs
      if ( (delts + dhsl) .ge. 0. ) then
         dhs = dhs + dhsl
         ebot = 0.0
         return
      endif

! model should never get here, probably a only minor error if it does
      print*,'melted completely through all ice and snow layer',delts
      print*,'some diagnostics',(delti(layer),layer=1,ni),(qi(layer)
     $      ,layer=1,ni)
      print*,'ERROR in botmelt',Ebot
      alarm=.true.

      return
      end

      subroutine adjust(hi,dhib,dhit,dhif,dhsf,qiflood,qigrow,ni,qi_tw)

! Adjusts temperature profile after melting/growing.

! E_tw = E at start
! E is the energy density after updating Ti from
! the heat equation without regard dhit and dhib

! h is the thickness from previous time step
! h_tw is the NEW thickness h(i,j)
! dhib is negative if there is melt at the bottom
! dhit is negative if there is melt at the top

! generally _tw is a suffix to label the new layer spacing variables

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"

! input variables
      real hi   ! initial ice thickness (m)
      real dhib ! ice bot, dhib<0 if melt (m)
      real dhit ! ice top, dhit<=0 (m)
      real dhif ! ice top from flooding, dhif>0 (m)
      real dhsf ! snow from flooding, dhsf>0 (m)

      real qiflood ! qi for flooded ice (J/m**3)
      real qigrow  ! qi for ice growing on bot (J/m**3)
      integer ni

! input and output variables
      real qi_tw(nmax) ! qi for adjusted layers (J/m**3)

! local variables
!     the following pairs are before & after adjusting
      real     h, h_tw      ! ice thickness (m)
      real     z(nmax+1)    ! vertical layer position (m)
      real     z_tw(nmax+1) ! vertical layer position (m)
      real     D, D_tw      ! layer thickness (m)
      integer  k, k_tw      ! index of layer

      real qi(nmax)         ! qi of initial layers (J/m**3)
      real fract(nmax,nmax) ! fract of layer k that makes up layer k_tw
      integer layer

! first check if there is any ice left after top/bottom melt/growth
      h = hi
      h_tw = hi+dhib+dhit

      if ( h_tw .le. 0. ) then
         do layer = 1,ni
            qi_tw(layer) = 0.
         enddo
         dhsf = 0.
         dhif = 0.
         return
      endif

! there is ice left so allow snow-ice conversion from flooding
      h_tw = h_tw+dhif

      if ( (abs(dhib) .gt. tiny) .or. (abs(dhit) .gt. tiny)
     $     .or. (dhif.gt.0)) then
         do layer = 1,ni
            qi(layer) = qi_tw(layer)
         enddo

! layer thickness
         D = h/ni
         D_tw = h_tw/ni

! z is positive down and zero is relative to the top
! of the ice from the old time step
         z(1) = -dhit          ! necessary
         z_tw(1) = -dhit-dhif
         do layer = 2,ni
           z(layer) = D*(layer-1)
           z_tw(layer) = z_tw(1)+D_tw*(layer-1)
         enddo
         z(ni+1) = h + min(dhib,0.)
         z_tw(ni+1) = z_tw(1)+h_tw

         do k_tw = 1,ni
           qi_tw(k_tw) = 0.
           do k = 1,ni
             fract(k,k_tw) = ( min( z_tw(k_tw+1), z(k+1)) -
     $                         max( z_tw(k_tw)  , z(k)  )   )
             if (fract(k,k_tw).gt.0.)
     $            qi_tw(k_tw) = qi_tw(k_tw)+fract(k,k_tw)*qi(k)
           enddo
         enddo

         if (dhif.gt.D_tw) then
           do k = 1,ni
             qi_tw(k) = qi_tw(k) + qiflood*max(0.,
     $                   min(D_tw,dhif-D_tw*(k-1)))
           enddo
         else
           qi_tw(1) = qi_tw(1) + dhif*qiflood
         endif

         if (dhib.gt.D_tw) then
           z(ni+2) = h + dhib
           do k = 1,ni
             qi_tw(k) = qi_tw(k) +  qigrow*max(0.,
     $                   (min(z(ni+2),Z_tw(k+1))-max(h,Z_tw(k))))
           enddo
         else
           qi_tw(ni) = qi_tw(ni) + max(dhib,0.)*qigrow
         endif

         do k_tw = 1,ni
           qi_tw(k_tw) = qi_tw(k_tw)/D_tw
         enddo

      endif

      return
      end

      subroutine init_ice_cpts (is, ie, js, je)

! if the cpts ice option is used without the embm, figure out where to
! call this routine

      implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "csbc.h"
      include "cpts.h"
      include "ice.h"

      integer i, ie, is, j, je, js, nc, layers

! use double/single Euler steps
      dtau=dts
      idx=1
      if (mod(nats,2).eq.0) idx=2

      if (nats.eq.1) then   ! mixing step, mode 2 is kept
        do nc=1,ncat
          do j=js,je
            do i=is,ie
              A(i,j,1,nc)=A(i,j,2,nc)
              heff(i,j,1,nc)=heff(i,j,2,nc)
              hseff(i,j,1,nc)=hseff(i,j,2,nc)
              Ts(i,j,1,nc)=Ts(i,j,2,nc)
            enddo
          enddo
        enddo
        do layers=1,ntilay
          do j=js,je
            do i=is,ie
              E(i,j,1,layers)=E(i,j,2,layers)
            enddo
          enddo
        enddo
      endif

      return
      end

      subroutine adv_ridge_cpts (is, ie, js, je)

      implicit none

# if defined O_ice_evp
      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "cpts.h"
      include "ice.h"

      integer i, ie, is, j, je, js, layer, nc

      real A0(is:ie,js:je), ah, A_Ts(is:ie,js:je,ncat)

      do nc=1,ncat
         call embmbc (heff(1,1,idx,nc))
         call embmbc (A(1,1,idx,nc))
         call embmbc (Ts(1,1,idx,nc))
         call embmbc (hseff(1,1,idx,nc))
      enddo
      do layer=1,ntilay
         call embmbc (E(1,1,idx,layer))
      enddo

!     compute open water area so it can be advected
      do j=js,je
        do i=is,ie
          A0(i,j)=1.0
          do nc=1,ncat
            A0(i,j)=A0(i,j)-A(i,j,idx,nc)
          enddo
          A0(i,j)=max(0.,A0(i,j))
        enddo
      enddo

!     proper surface temperature advection is done by advecting
!     the area-temperature product
      do nc=1,ncat
        do j=js,je
          do i=is,ie
            A_Ts(i,j,nc) = min(A(i,j,idx,nc)*Ts(i,j,idx,nc),0.)
          enddo
        enddo
      enddo

! At this time mechanical redistribution for this model assumes
! you are using upstream advection.
      call advupb(uice, vice ,A0, is, ie, js, je)
      do nc=1,ncat
        call advupb(uice, vice, A(1,1,idx,nc), is, ie, js, je)
      enddo
      do nc=1,ncat
        call advupb(uice, vice, A_Ts(1,1,nc), is, ie, js, je)
      enddo
      do nc=1,ncat
        call advupb(uice, vice, heff(1,1,idx,nc), is, ie, js, je)
      enddo
      do nc=1,ncat
        call advupb(uice, vice, hseff(1,1,idx,nc), is, ie, js, je)
      enddo
      do layer=1,ntilay
        call advupb(uice, vice, E(1,1,idx,layer), is, ie, js, je)
      enddo

      do nc=1,ncat
        do j=js,je
          do i=is,ie
            Ts(i,j,idx,nc)=frzpt(i,j)
            if (A(i,j,idx,nc).gt.tiny) then
              Ts(i,j,idx,nc)=A_Ts(i,j,nc)/A(i,j,idx,nc)
              Ts(i,j,idx,nc)=min(Ts(i,j,idx,nc),0.)
            endif
          enddo
        enddo
      enddo

      call mechred(A0, is, ie, js, je)

# endif

      return
      end

! Energy-conserving sea ice model
! Routines to do basic energy/temperature conversion

! See Bitz, C.M., and W.H. Lipscomb, 1999:
! A new energy-conserving sea ice model for climate study,
! J. Geophys. Res., in press.

! The author grants permission to the public to copy and use this
! software without charge, provided that this Notice and any statement
! of authorship are reproduced on all copies.

      function energ(Tmp,sal)
!     compute the energy of melting for non-new ice
!     relative to melting (negative quantity)
!     assuming Tmp is in Celsius
      implicit none
      include "cpts.h"
      include "thermo.h"

      real Tmp,energ,sal

      energ=-rflice
     $        -rcpice*(-alpha*sal-Tmp)-gamma*sal/Tmp

      return
      end

      subroutine getTmp(Tmp,E,sal1d,ni)

!     compute the midpoint temperature from the
!     energy of melting for each layer

      implicit none

      include "cpts.h"
      include "thermo.h"

      real Tmp(nmax),E(nmax),sal1d(nmax),quad
      integer ni,i

      real(kind=8) q, B, C

      do i=1,ni
        q=E(i)+rflice-rcpice*alpha*sal1d(i)
        B=-q/rcpice
        C=-gamma*sal1d(i)/rcpice
        Tmp(i)= quad(B,C)
      enddo

      return
      end

      function quad(B,C)
      implicit none
      include "cpts.h"
      include "thermo.h"

      real quad,root,T1,T2
      real(kind=8) B, C, B_2

      B_2=B/2.0
      root=sqrt(B_2*B_2-C)
!      T1=-B_2+root
      T2=-B_2-root
      quad=T2
!      print*,'quad T1=',T1,T2
!      if (T1.le.(TMELT-0.1)) print*,'quad T1=',T1,T2
      return
      end

      subroutine grownew(fsh_new,Io_new,flo,focean,Fbot,
     $            hnew,tnew,tair,qair,ug,dtp,fice,frzpt,evp,
     $            dswr,ultnt,usens,ulwr)

!     Uses a Newton_Raphson (non-linear)
!     method to compute new ice surface temp
!     and thickness.
!     This algorithm tends to fail when the correct solution is
!     to grow a very thin layer of ice. In this case I just set
!     the new ice surface temperature to 0.5*(Tair+frzpt) and
!     grow accordingly. It still conserves energy.

      implicit none

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"
      include "calendar.h"

      real  coltnt,cosens
      real  fsh_new,flo,Io_new
      real  Fbot,focean,hnew,tnew,frzpt
      real  dswr,ultnt,usens,ulwr
      real  tair,qair,evp,qnew,rday
      real  ug,dt
      real  dtp
      real  b,fice,e1,e2,e3,e4,e5,terror,herror
      real  u
      real  J11,J12,J21,J22,F1,F2,detJ
      real  rand,energ,kthin
      real  ult_new,use_new,ulw_new,evp_new
      integer iter

      kthin=kappai+beta*salnew/(frzpt-5.)
      Tnew=0.5*(frzpt+Tair) ! initial guess
      Tnew=min(frzpt,Tnew)

# if defined O_ice_cpts_simple_growth
       fice=focean
       u=0.5*(frzpt+Tnew)
       u=energ(u,salnew)
       hnew=(focean+Fbot)*dtp/u
# else
      dt=Tnew-Tair
      cosens=0.94*rcpatm*dalt_i*ug

      ult_new=ultnt
      e1=fsh_new-ult_new+flo+cosens*tair

      e2=-cosens*Tnew-ESICE*(Tnew+TFFRESH)**4
      fice=e1+e2

      if ((fice+focean+2*Fbot).lt.0.0) then
         rday=daylen/dtau
         hnew=centi*0.0133/rday*(frzpt-Tnew)

         iter=1

 100     continue

            e2=-cosens*Tnew-esice*Tnew**4
            e3=-cosens-4.*ESICE*(Tnew+TFFRESH)**3

            if (Tnew.gt.frzpt) then
               u=frzpt
            else
               u=0.5*(frzpt+Tnew)
            endif
            e4=-0.5*(rcpice-gamma*salnew/(u**2))
            u=-energ(u,salnew)

            J11=hnew*e3-kthin
            J12=e1+e2
            J21=0.5*e3+hnew/dtp*e4
            J22=u/dtp
            F1=(J12-Io_new)*hnew+kthin*(frzpt-Tnew)
            F2=0.5*(focean+J12)+Fbot+hnew*J22
            detJ=J11*J22-J12*J21

            terror=(J22*F1-J12*F2)/detJ
            herror=(-J21*F1+J11*F2)/detJ
            Tnew= Tnew-terror
            hnew=hnew-herror
            iter=iter+1
            if ((iter.gt.10).or.(hnew.lt.0).or.(Tnew.gt.0.)) then
!     when cannot converge, assume a temperature and continue
               Tnew=0.5*(frzpt+Tair)
               Tnew=min(frzpt,Tnew)
               go to 200
            endif
            if (abs(terror).gt.0.005) go to 100

 200        continue

            u=0.5*(frzpt+Tnew)
            u=energ(u,salnew)

            ulw_new=esice*(Tnew+TFFRESH)**4-Flo
            use_new=cosens*(Tnew-tair)

            fice=fsh_new-ulw_new-use_new-ult_new
            hnew=(0.5*(fice+focean)+Fbot)*dtp/u

            dswr = 0.5*(dswr+fsh_new)
            ulwr = 0.5*(ulwr +ulw_new)
            usens = 0.5*(usens+use_new)

      else
!     it is not clear what to do with the negative
!     ocean energy balance if in the presence of
!     open water the energy balance is negative and
!     would allow ice to grow, yet once it grows a little
!     ice the energy balance swings to positive

!     I am going to set h=0 in this case and balance the fluxes
!     by messing with the sensible heat

         usens = dswr-ultnt-ulwr+Fbot
         fice=-focean-2*Fbot
         Tnew=frzpt
         hnew=0
      endif

# endif
      return
      end

      subroutine mechred (A0, is, ie, js, je)

!  it is perfectly fine to call this routine when NCAT=1
!  in fact, it will provide a source of open water from shear
!  deformation. See Stern et al, 1995 for information about
!  how and why this is done

      implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "evp.h"
      include "atm.h"
      include "cpts.h"
      include "ice.h"
      include "grdvar.h"

!   input variables
      real A0(is:ie,js:je)

!    local variables

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1, n, nc
      integer layer, l

      real En(NMAX,ncat)
      real Esn(ncat)
      real heffn(ncat)
      real hseffn(ncat)
      real An(ncat)
      real Tsn(ncat)
      real closing       ! lead closing rate
      real M(ncat,ncat)  ! area fraction of ice from i into j
      real HN(ncat,ncat) ! H*volume fraction of ice from i into j
      real rpress

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

# if defined O_ice_cpts_roth_press
      do j=js,je
        do i=isp1,ie
          strength(i,j,idx)=0.
        enddo
      enddo
# endif

      do n=1,nseg
        do j=jsi(n),jei(n)+1
          do i=isp1,iem1
            if (tmsk(i,j) .ge. 0.5) then

!    ice velocity divergence and lead closing rate, following FH95
              closing=0.25*(del(i,j)-abs(eI(i,j)))-min(eI(i,j),0.)
              closing=max(closing,-eI(i,j))

              l=1
              do nc=1,ncat
                do layer=1,nilay(nc)
                   En(layer,nc)=E(i,j,idx,l)
                   l=l+1
                enddo
                Tsn(nc)=Ts(i,j,idx,nc)
                An(nc)=A(i,j,idx,nc)
                heffn(nc)=heff(i,j,idx,nc)
                hseffn(nc)=hseff(i,j,idx,nc)
              enddo

              call ridge(heffn,hseffn,A0(i,j),An,En,Esn,Tsn,
     $            closing,eI(i,j),M,HN,frzpt(i,j))

# if defined O_ice_cpts_roth_press
              strength(i,j,idx)=rpress(heffn,A0(i,j),An,HN,M)
# endif
              l=1
              do nc=1,ncat
                do layer=1,nilay(nc)
                   E(i,j,idx,l)=En(layer,nc)
                   l=l+1
                 enddo
                 Ts(i,j,idx,nc)=Tsn(nc)
                 A(i,j,idx,nc)=An(nc)
                 heff(i,j,idx,nc)=heffn(nc)
                 hseff(i,j,idx,nc)=hseffn(nc)
              enddo

            endif
          enddo
        enddo
      enddo

      do nc=1,ncat
         call embmbc (heff(1,1,idx,nc))
         call embmbc (A(1,1,idx,nc))
         call embmbc (Ts(1,1,idx,nc))
         call embmbc (hseff(1,1,idx,nc))
      enddo
      do layer=1,ntilay
         call embmbc (E(1,1,idx,layer))
      enddo
      return
      end

      subroutine ridge(heffn,hseffn,A0,An,En,Esn,Tsn,closed,divu,
     $     M,HN,frzpt)

!     note the time available for ridging between the cats is
!     dtau, however have to subcycle if a category runs
!     out of ice then dtn is less than dtau

!     M(nc,k) is the fraction of ice that participates in ridging

      implicit none

      include "cpts.h"
      include "thermo.h"

      real heffn(ncat)
      real hseffn(ncat)
      real An(ncat)
      real Tsn(ncat),dtf(ncat),Atsn(ncat)
      real En(nmax,ncat)
      real Esn(ncat)
      real closed        ! the closing as defined by Stern et al
      real divu, A0
      real frzpt
      real strg          ! ice strength

      real Asum
      real W(0:ncat)     ! ridging for each category
      real Wa(0:ncat)    ! the area that participates
      real M(ncat,ncat)  ! area fraction of ice from i into j
      real N(ncat,ncat)  ! volume fraction of ice from i into j
      real HN(ncat,ncat) ! dummy
      real Hmean(ncat)   ! Hi+sqrt(cK*Hi)

      real dtn,dtnn,dtp
      integer k,layer,index,isub,nc,layk,lay,i,j

      real Qi(nmax,ncat) ! energy of melting ice per unit volume
      real Qs(ncat)      ! energy of melting of snow per unit area
      real delt(ncat),delts(ncat),deltEs(ncat)
      real deltE(nmax,ncat)

      real Psi, closdtn

      real dsum,ssum,esum,fsum
      real Hi(ncat),Hs(ncat), eta

      dtn=dtau

      dtp=0

      A0=A0
      Asum=A0
      do nc=1,ncat
         Asum = Asum + An(nc)
      enddo

! unfortunately the second order accurate advection
! schemes do not have Asum=1-divu*dt after advection
! owing to diffusion.
! Adjust opening rate here to compensate.
      divu=(1-Asum)/dtn

      closed=-min(divu,0.) ! from Bill Hibler suggestions
      closed=max(closed,-divu)

      A0=A0+dtn*(divu+closed)
      Asum=Asum+dtn*(divu+closed)

      if (Asum.le.1.0) then
        A0=A0+1.-Asum
# if defined O_ice_cpts_roth_press
        call ridge_matrices(heffn,An,Hi,Hmean,M,N,HN)
# endif
        return
      endif

 50   continue

      call ridge_matrices(heffn,An,Hi,Hmean,M,N,HN)
      call ridging_mode(A0,An,W,Wa,M)

# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
      do nc=1,ncat
         Hs(nc)=0.0
         Qs(nc)=0.0
         do layer=1,nilay(nc)
            Qi(layer,nc)=0.0
         enddo

         if ((An(nc).gt.TINY) .and. (heffn(nc).gt.0.))  then
            Hs(nc)=hseffn(nc)/An(nc)
            Qs(nc)=Esn(nc)/An(nc)
            do layer=1,nilay(nc)
               Qi(layer,nc)=En(layer,nc)*nilay(nc)/heffn(nc)
            enddo
         endif
      enddo

      dsum=0.
      ssum=0.
      esum=0.
      do nc=1,ncat-1
         delt(nc)  =-Hi(nc)*Wa(nc)
         delts(nc) =-Hs(nc)*Wa(nc)
         deltEs(nc)=-Qs(nc)*Wa(nc)
         do k=1,nc
            delt(nc)  =delt(nc)+Wa(k)*N(k,nc)
            delts(nc) =delts(nc)+Wa(k)*M(k,nc)*Hs(k)*Hmean(k)/hi(k)
            deltEs(nc)=deltEs(nc)+Wa(k)*M(k,nc)*Qs(k)*Hmean(k)/hi(k)
         enddo
         dsum=dsum+delt(nc)
         ssum=ssum+delts(nc)
         esum=esum+deltEs(nc)
      enddo

      delt(ncat)  =-dsum
      delts(ncat) =-ssum
      deltEs(ncat)=-esum

      fsum=0.
      do nc=1,ncat-1
         do layer=1,nilay(nc)
            deltE(layer,nc)=-Wa(nc)*Hi(nc)*Qi(layer,nc) /nilay(nc)
            do k=1,nc
               layk=(layer+ncrel(k,nc)-1)/ncrel(k,nc)
               deltE(layer,nc)=deltE(layer,nc)+Qi(layk,k)*Wa(k)*N(k,nc)
     $                        /nilay(nc)
            enddo
            fsum=fsum+deltE(layer,nc)
         enddo
      enddo
      do layer=1,nilay(ncat)
         deltE(layer,ncat)=0.
         do k=1,ncat-1
            layk=(layer+ncrel(k,ncat)-1)/ncrel(k,ncat)
            deltE(layer,ncat)=deltE(layer,ncat)+Qi(layk,k)*Wa(k)
     $                       *N(k,ncat)/nilay(ncat)
         enddo
         fsum=fsum+deltE(layer,ncat)
      enddo

!     geometrical factor for treating ice/snow surface temp when ridging
      do nc = 1,ncat
         dtf(nc) = - tsn(nc)*Wa(nc)
         do k=1,nc
           dtf(nc) = dtf(nc)+ tsn(k)*Wa(k)*M(k,nc)
         enddo
      enddo

# endif

!     make sure that only ridge as much ice as is available, otherwise
!     reduce timestep
      Psi=closed*Wa(0)
      if ((Psi*dtn.gt.A0) .and. (Psi.gt.0.)) then
         dtnn=A0/(Psi)
         dtn=min(dtn,dtnn)
      endif
      do nc=1,ncat
         Psi=closed*Wa(nc)
         if ((Psi*dtn.gt.An(nc)) .and. (Psi.gt.0.)) then
            dtnn=An(nc)/(Psi)
            dtn=min(dtn,dtnn)
         endif
      enddo
      dtp=dtp+dtn

      closdtn=closed*dtn

      A0=A0+closdtn*W(0)
      Asum=A0
      do nc=1,ncat
         Atsn(nc) = An(nc)*tsn(nc) + closdtn*dtf(nc)
         An(nc)=An(nc)+closdtn*W(nc)
         Asum=Asum+An(nc)
# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
         heffn(nc)=heffn(nc)+closdtn*delt(nc)
         hseffn(nc)=hseffn(nc)+closdtn*delts(nc)
         Esn(nc)=Esn(nc)+closdtn*deltEs(nc)
         do layer=1,nilay(nc)
            En(layer,nc)=En(layer,nc)+closdtn*deltE(layer,nc)
         enddo
# endif
      enddo
      do nc = 1,ncat
        if (An(nc).gt.tiny) then
          tsn(nc) = Atsn(nc)/An(nc)
        else
          tsn(nc) = frzpt
        endif
      enddo

      dtn=dtau-dtp

      if (dtn.gt.TINY) go to 50

      return
      end

      subroutine ridging_mode(A0,An,W,Wa,M)

!     compute W and Wa
!     M(nc,k) is the fraction of ice that participates in ridging
!     from category nc that goes into category k

      implicit none

      include "cpts.h"

      real An(ncat),A0

      real W(0:ncat)       ! ridging for each category
      real Wa(0:ncat)      ! the area that participates
      real M(ncat,ncat)    ! area fraction of ice from i into j

      real V(1:ncat+1)     ! the area that participates
      real Asum,eta,Wsum
      integer k,layer,index,isub,nc,layk,lay,i,j

      Asum=0.
      Wa(0)=min(gstar-Asum,A0)
      Asum=Asum+Wa(0)
!      print*,0,Wa(0),Asum,(1.0 + 0.5*Wa(0)/gstar - Asum/gstar)
      Wa(0)=Wa(0)*(1.0 + 0.5*Wa(0)/gstar - Asum/gstar)
      Wa(0)=max(Wa(0),0.)
      do nc=1,ncat
         if (An(nc).gt.0.) then
            Wa(nc)=min(gstar-Asum,An(nc))
            Asum=Asum+Wa(nc)
            Wa(nc)=Wa(nc)*(1.0 + 0.5*Wa(nc)/gstar - Asum/gstar)
            Wa(nc)=max(Wa(nc),0.)
         else
            Wa(nc)=0.
         endif
      enddo

      Wsum=-Wa(0)
      do nc=1,ncat
         W(nc)=-Wa(nc)
         do k=1,nc
            W(nc)=W(nc)+Wa(k)*M(k,nc)
         enddo
         Wsum=Wsum+W(nc)
      enddo

      eta=-Wsum
      do nc=1,ncat
         W(nc)=W(nc)/eta
         Wa(nc)=Wa(nc)/eta
      enddo
      Wa(0)=Wa(0)/eta
      W(0)=-Wa(0)

      return
      end

      subroutine ridge_matrices(heffn,An,Hi,Hmean,M,N,HN)

!     compute M,N,HN for ridge and press for true thicknesses

      implicit none

      include "cpts.h"

      real heffn(ncat)
      real An(ncat)

      real Hi(ncat)
      real M(ncat,ncat)    ! area fraction of ice from i into j
!     N is used in ridge and HN is used in rpress
      real N(ncat,ncat)    ! volume fraction of ice from i into j
      real HN(ncat,ncat)   ! H*volume fraction of ice from i into j

      real rcKH            ! sqrt(cK*Hi)
      real Hmean(ncat)     ! Hi+sqrt(cK*Hi)
      real bet             ! 2*(K-Hi)
      integer n1,n2        ! see above
      integer k,nc
      logical row_flag(ncat)

      do nc=1,ncat
         if ((An(nc).gt.TINY) .and. (heffn(nc).gt.0.))  then
            Hi(nc)=heffn(nc)/An(nc)
            row_flag(nc)=.true.
            do k=1,ncat
              M(nc,k)=0.
              N(nc,k)=0.
              HN(nc,k)=0.
            enddo
         else
            Hi(nc)=hstar(nc-1)
            row_flag(nc)=.false.
            do k=1,ncat
              M(nc,k)=M_def(nc,k)
              N(nc,k)=N_def(nc,k)
              HN(nc,k)=HN_def(nc,k)
            enddo
         endif
      enddo

      call comp_matrices(row_flag,Hi,Hmean,M,N,HN)

      return
      end

      subroutine comp_matrices(row_flag,Hi,Hmean,M,N,HN)

!     compute M,N,HN for ridge and rpress for true thicknesses
!     N is used in ridge and HN is used in rpress

      implicit none

      include "cpts.h"

      logical row_flag(ncat)
      real Hi(ncat)
      real M(ncat,ncat)    ! area fraction of ice from i into j
      real N(ncat,ncat)    ! volume fraction of ice from i into j
      real HN(ncat,ncat)   ! H*volume fraction of ice from i into j

      real rcKH            ! sqrt(cK*Hi)
      real Hmean(ncat)     ! Hi+sqrt(cK*Hi)
      real bet             ! 2*(K-Hi)
      integer n1,n2        ! see above
      integer k,layer,index,isub,nc,layk,lay,i,j

      do nc=1,ncat
         rcKH=sqrt(cK*Hi(nc))
         Hmean(nc)=Hi(nc)+rcKH
         bet=2*(cK-Hi(nc))

         if (row_flag(nc)) then
! ice from category nc is ridging into a distribution from n1 to n2
            n1=nc  ! first guess
            n2=n1  ! first guess
            do k=nc+1,ncat
               if (2*Hi(nc).ge.hstar(k-1)) n1=k
            enddo
            if (n1.eq.ncat) then
!              case when all ice ridges into thickest cat
               M(nc,n1)=Hi(nc)/Hmean(nc)
               N(nc,n1)=Hi(nc)
               HN(nc,n1)=Hi(nc)*Hmean(nc)
            elseif (2*rcKH.le.hstar(n1)) then
!              case when all ice ridges into just one category
               M(nc,n1)=Hi(nc)/Hmean(nc)
               N(nc,n1)=Hi(nc)
               HN(nc,n1)=Hi(nc)*Hmean(nc)
            else
               n2=ncat
               do k=ncat-1,n1+1,-1
                  if (2*rcKH.lt.hstar(k)) n2=k
               enddo
               M(nc,n1)=(hstar(n1)-2*Hi(nc))/bet
               N(nc,n1)=M(nc,n1)*(Hi(nc)+hstar(n1)/2.)
               HN(nc,n1)=N(nc,n1)*(Hi(nc)+hstar(n1)/2.)

               M(nc,n2)=(2.*rcKH-hstar(n2-1))/bet
               N(nc,n2)=M(nc,n2)*(hstar(n2-1)/2.+rcKH)
               HN(nc,n2)=N(nc,n2)*(hstar(n2-1)/2.+rcKH)

               do k=n1+1,n2-1
                  M(nc,k)=(hstar(k)-hstar(k-1))/bet
                  N(nc,k)=M(nc,k)*(hstar(k)+hstar(k-1))/2.
                  HN(nc,k)=N(nc,k)*(hstar(k)+hstar(k-1))/2.
               enddo
            endif
         endif
      enddo

      return
      end

      function rpress(heffn,A0,An,HN,M)
!     Compute the ice strength based on Rothrock, 1975 and also see FH95
!     does not make sense to do this unless the ice that participates in
!     ridging is well resolved -- must have about 5 categories, 10 would
!     be better

!     M(nc,k) is the fraction of ice that participates in ridging
!     from category nc that goes into category k
!     n1 is the smallest k with nonzero fraction
!     n2 is the largest  k with nonzero fraction
!     the ice that participates has thickness Hi(nc)
!     and it ridges up to linear distribution between
!     2*Hi(nc) and 2*sqrt(cK*Hi(nc))

      implicit none
      include "cpts.h"

      real rpress
      real heffn(ncat)
      real An(ncat)

      real W(0:ncat)         ! ridging for each category
      real Wa(0:ncat)        ! the area that participates
      real M(ncat,ncat)    ! area fraction of ice from i into j
      real HN(ncat,ncat)   ! H*volume fraction of ice from i into j

      real A0

      integer n1,n2        ! see above
      integer k,nc,i,j

      real Hi(ncat),Hs(ncat)

      real peterm, frterm

!     ideally cpe should be a function of density (a very sensitive one)
!     see FH95, I set it equal to a number so you are more likely to get
!     a reasonable ice strength.
!     It has units (kg/m**3)*(m/s**2) = (kg/m**2/s)

      rpress=0.

# if defined O_ice_cpts_roth_press
      if (A0.ge.gstar) return

      do nc=1,ncat
         if ((An(nc).gt.TINY) .and. (heffn(nc).gt.0.))  then
            Hi(nc)=heffn(nc)/An(nc)
         else
            Hi(nc)=hstar(nc-1)
         endif
      enddo

!     strictly speaking, should recompute the M and HN matrices here,
!     but they would be only very slightly different from what they were
!     before ridging and it does not affect the mass conservation
      call ridging_mode(A0,An,W,Wa,M)

      peterm=0.
      do nc=1,ncat
         peterm=peterm-Wa(nc)*Hi(nc)*Hi(nc)
         do k=1,nc
            peterm=peterm+Wa(k)*HN(k,nc)
         enddo
      enddo
      peterm=max(0.,cpe*peterm)

      rpress=peterm*17.

# endif
      return
      end

      subroutine movedown(frzpt,heff1,hseff1,A1,hsnow1,hice1,h1min,Es1,
     $     ts1,nn1,heff2,hseff2,A2,hsnow2,hice2,h2min,Es2,ts2,nn2,E1,E2)

!     hice and hsnow are outputs (not needed for input)

      implicit none

      include "cpts.h"
      real heff1,hseff1,A1,hsnow1,hice1,h1min,Es1,E1(NMAX),ts1
      real heff2,hseff2,A2,hsnow2,hice2,h2min,Es2,E2(NMAX),ts2
      real frzpt
      integer nn1,nn2,M,layer,index,isub

      real x(NMAX)
      real Atemp

      Atemp=A1+A2
      if (Atemp.le.0.0) return

      heff1=heff1+heff2
      hseff1=hseff1+hseff2
      ts1=ts1*A1+ts2*A2

      A1=Atemp

      hsnow1=hseff1/A1
      hice1=heff1/A1
      ts1=ts1/A1

      Es1=Es1+Es2

      M=nn2/nn1
!     combine the N2 layers of energy in category 2
!     into N1 portions
      do 1000 layer=1,nn1
         x(layer)=0.
         do 1010 isub=1,M
            index=(layer-1)*M+isub
            x(layer)=x(layer)+E2(index)
 1010    continue
 1000 continue

!     sum the energies
      do layer=1,nn1
         E1(layer)=E1(layer)+x(layer)
      enddo

      call zerocat(heff2,hseff2,A2,hsnow2,hice2,h2min,
     $     Es2,E2,ts2,frzpt,nn2)

      return
      end

      subroutine moveup(frzpt,heff1,hseff1,A1,hsnow1,hice1,h1min,Es1,
     $   ts1,nn1,heff2,hseff2,A2,hsnow2,hice2,h2min,Es2,ts2,nn2,E1,E2)

!     hice and hsnow are outputs (not needed for input)

      implicit none

      include "cpts.h"
      real heff1,hseff1,A1,hsnow1,hice1,h1min,Es1,E1(NMAX),ts1
      real heff2,hseff2,A2,hsnow2,hice2,h2min,Es2,E2(NMAX),ts2
      real frzpt
      integer nn1,nn2,M,layer,index,isub

      real x(NMAX)
      real Atemp

      Atemp=A1+A2
      if (Atemp.le.0.0) return

      heff2=heff1+heff2
      hseff2=hseff1+hseff2
      ts2=ts1*A1+ts2*A2

      A2=Atemp

      hsnow2=hseff2/A2
      hice2=heff2/A2
      ts2=ts2/A2

      Es2=Es1+Es2

      M=nn2/nn1
!     move energy of nn1 layers from 1 into the nn2 of 2
      do index=1,nn1
         do isub =1,M
            layer=M*(index-1)+isub
            x(layer)=E1(index)/M
         enddo
      enddo

!     sum the energies
      do layer=1,nn2
         E2(layer)=E2(layer)+x(layer)
      enddo

      call zerocat(heff1,hseff1,A1,hsnow1,hice1,h1min,Es1,E1,ts1,frzpt,
     $             nn1)

      return
      end

      subroutine zerocat(heff,hseff,A,hsnow,hice,hmin,Es,E,ts,frzpt,nng)

      implicit none

      include "cpts.h"
      real heff,hseff,A,hsnow,hice,hmin,Es,E(NMAX)
      real ts, frzpt
      integer nng,layer

      heff=0.0
      hseff=0.0
      A=0.0
      hsnow=0.0
      hice=hmin
      Es=0.0
      do layer=1,nng
         E(layer)=0.0
      enddo
      ts=frzpt
      return
      end

      subroutine thermo (is, ie, js, je)

!     heat budget

!     compute fluxes, temp and thickness  over open water and ice

!     flux is positive toward the atm/ice or atm/ocean surface even if
!     the flux is below the surface

      implicit none

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
      include "cpts.h"
      include "ice.h"
      include "thermo.h"
      include "coord.h"
      include "grdvar.h"

      real  tocn           ! surface ocean temperature in Celsius
      real  tair           ! surface air temperature in celcius
      real  qair           ! surface air humidity
      real  ug             ! geostrophic wind speed

      real  evp(0:ncat)    ! evapourative water flux for each category
      real  wff            ! ice-ocean water flux for each category

      real  dswr(0:ncat)   ! downward shortwave for each category
      real  ultnt(0:ncat)  ! upward latent heat for each category
      real  usens(0:ncat)  ! upward sensible heat for each category
      real  ulwr(0:ncat)   ! upward longwave for each category

      real  Io             ! solar transmitted through top surface
      real  Ib             ! solar transmitted through bottom surface
      real  absorb         ! sum of Ib*A for each category

      real  Focean         ! net flux into ocean from atmosphere
      real  Flwatm         ! downward longwave from atmosphere (not net)
      real  Fbot           ! flux from ocean due to thermal relaxation
      real  Flat           ! flux from ocean due to lateral melt

      real  Fice           ! net flux into ice from atmosphere
      real  Fneti          ! Fice + conductive flux
      real  Fx             ! adjusted Fbot for when ice melts away
      real  condb          ! conductive flux in ice at bottom surface

      real  A0             ! open water fraction
      real  Fsh_new        ! solar incident on new ice (over open water)
      real  Io_new         ! solar transmitted through top of new ice
      real  hnew           ! thickness of new ice
      real  Tnew           ! surface temperature of new ice

      real  hi(ncat)       ! ice thickness for each category
      real  hs(ncat)       ! snow thickness for each category
      real  En(1:NMAX,ncat)! energy of melt of each layer and category
      real  qi(1:NMAX,ncat)! like Ei but per unit volume
      real  Ti(0:NMAX,ncat)! Temperature of each layer and category (C)
      real  Tbot(ncat)     ! Temperature of ice bottom and category (C)

      real  delhib         ! ice thickness change at bottom
      real  delhit         ! ice thickness change (no sublimation)
      real  delhs          ! snow thickness change (no sublimation)
      real  subi           ! ice thickness change from sublimation
      real  subs           ! snow thickness change from sublimation
      real  delf           ! ice thickness change from flooding
      real  delfs          ! snow thickness change from flooding

      integer i,j,l,layer,nc,nr  ! counters

      real sblmt            ! heat to ice-snow from sublimation (J/m**2)
      real htatm            ! heat to ice-snow from atm (J/m**2)
      real wtatm            ! water to ice-snow from atm (m/s)
      real htocn            ! heat to ice-snow from ocn (J/m**2)
      real wtocn            ! waterto ice-snow from ocn (m/s)

      real einit            ! initial energy in ice/snow (J/m**2)
      real winit            ! initial water equiv. in ice/snow (m)
      real efinl            ! final energy in ice/snow (J/m**2)
      real wfinl            ! final water equiv. in ice/snow (m)
      real ediff            ! energy difference in ice/snow (J/m**2)
      real wdiff            ! water difference  in ice/snow (m)

      logical alarm

! vars for small stuff check and computing Ti, hi, hs and albedo
      real Aleft, tio, Enet

! lateral melt vars
      real Esum   ! sum of energy in ice and snow (J/m**2)
      real mrate  ! lateral melt rate (m/s)
      real arate  ! rate ice concentration is melting laterally (s)
      real DAice  ! change in concentration from lateral melt
      real Vfrac  ! ratio of volume after to before melting
      real Dhlatm, Dslatm, sla, fracsnow

! redistribute the ice categories vars
      real h0prime,A0prime,h1prime,A1prime,E0,Ad0
      real energ,Tmp,epsilon

      real ca, as, aicet, coaloff, C2K

      integer ie, iem1, is, isp1, je, jem1, js, jsp1

!     dummy energy of melting in snow for calls to movedown/up
      real esdum1, esdum2
!-----------------------------------------------------------------------

      sla = zw(1)*secday/dampice/2.389e-8
      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1
      C2K = 273.15

      do j=jsp1,jem1
        do i=isp1,iem1
          if (tmsk(i,j) .ge. 0.5) then

            Io_new=0.0
!     adjust albedo for new ice growth
# if defined O_sulphate_data_transient
            fsh_new = solins(i,j)*sbc(i,j,iaca)*pass
     &                *max(c0, sbc(i,j,isca) - sulph(i,j,2))
# else
            fsh_new = solins(i,j)*sbc(i,j,iaca)*pass*sbc(i,j,isca)
# endif
            tocn = sbc(i,j,isst)
            tair = at(i,j,2,isat)
# if defined O_sealev || defined O_sealev_data
     &           - elev_sealev(i,j)*rlapse
# endif
            qair = at(i,j,2,ishum)
            Fbot = -sla*(frzpt(i,j) - tocn)
            ug = sbc(i,j,iws)
            Flwatm = esatm*(tair + C2K)**4
            Focean = dnswr(i,j) - uplwr(i,j) - upsens(i,j) - upltnt(i,j)
            dswr(0) = dnswr(i,j)
            ultnt(0) = upltnt(i,j)
            usens(0) = upsens(i,j)
            ulwr(0) = uplwr(i,j)
            evp(0) = evap(i,j)

!-----------------------------------------------------------------------
!     make simple energy of melt matrix for speed
!-----------------------------------------------------------------------
            do nc=1,ncat
              do layer=1,nilay(nc)
                En(layer,nc)=E(i,j,idx,layer1(nc)+layer-1)
              enddo
            enddo

!-----------------------------------------------------------------------
!     initial energy and water diagnostics
!     could be eliminated to save time
!-----------------------------------------------------------------------
            einit = 0.
            winit = 0.
            do nc=1,ncat
              do layer=1,nilay(nc)
                einit=einit+En(layer,nc)
              enddo
              einit = einit-hseff(i,j,idx,nc)*rflsno
              winit = winit+(heff(i,j,idx,nc)*rhoice
     $                        +hseff(i,j,idx,nc)*rhosno)
            enddo

# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
!-----------------------------------------------------------------------
!     get rid of small area by moving to other category (if possible)
!     does not fix negative ice area/volumes (that comes later)

!     push small bit of ice to thicker category
!     hi and hs are outputs of moveup (not inputs)
!-----------------------------------------------------------------------
            do nc=1,(ncat-1)
              if ((A(i,j,idx,nc) .le. ASMALL(nc))
     $             .and.(A(i,j,idx,nc) .gt. 0.0) ) then

                A(i,j,idx,nc)=heff(i,j,idx,nc)/hstar(nc)
                call moveup(frzpt(i,j),heff(i,j,idx,nc),
     $               hseff(i,j,idx,nc),
     $               A(i,j,idx,nc),hs(nc),hi(nc),hstar(nc-1),
     $               esdum1,Ts(i,j,idx,nc),nilay(nc),
     $               heff(i,j,idx,nc+1),hseff(i,j,idx,nc+1),
     $               A(i,j,idx,nc+1),hs(nc+1),hi(nc+1),hstar(nc),
     $               esdum2,Ts(i,j,idx,nc+1),nilay(nc+1),
     $               En(1,nc),En(1,nc+1))
              endif
            enddo

!     push small bit of ice to thinner category
            do nc=ncat,2,-1
              if ((A(i,j,idx,nc) .le. ASMALL(nc))
     $             .and.(A(i,j,idx,nc) .gt. 0.0) ) then
!     reshape the thicker ice to thickness of thinner category
!     while making certain not to make area > 1
                Aleft=1.0
                do  nr=1,(nc-1)
                  Aleft=Aleft-A(i,j,idx,nr)
                enddo
                A(i,j,idx,nc)=min(heff(i,j,idx,nc)/hstar(nc-1),Aleft)
                call movedown(frzpt(i,j),heff(i,j,idx,nc-1),
     $               hseff(i,j,idx,nc-1),
     $               A(i,j,idx,nc-1),hs(nc-1),hi(nc-1),hstar(nc-2),
     $               esdum2,Ts(i,j,idx,nc-1),nilay(nc-1),
     $               heff(i,j,idx,nc),hseff(i,j,idx,nc),A(i,j,idx,nc),
     $               hs(nc),hi(nc),hstar(nc-1),
     $               esdum1,Ts(i,j,idx,nc),nilay(nc),
     $               En(1,nc-1),En(1,nc))
              endif
            enddo
# endif

!-----------------------------------------------------------------------
!     get rid of small open water area
!-----------------------------------------------------------------------

            A0=1.0
            do nc=1,ncat
              A0=A0-max(A(i,j,idx,nc),0.0)
            enddo
            if (A0.lt.A0SMALL) then
!     put tiny open water areas out of misery by reshaping
!     the ice starting with the thin category
              if (A0.ge.0.) then
                do nc=1,ncat
                  if (A(i,j,idx,nc).gt.0.) then
                    A(i,j,idx,nc)=A(i,j,idx,nc)+A0
                    A0=0.0
                  endif
                enddo
              else
!     if the open water fraction is negative then reduce the area
!     of the category with the thin ice first
                nr=1
                do nc=2,ncat
                  if ((A(i,j,idx,nr)+A0).lt.ASMALL(nr)) nr=nc
                enddo
                A(i,j,idx,nr)=A(i,j,idx,nr)+A0
              endif
            endif

!-----------------------------------------------------------------------
!     fix negative ice/snow area/volumes and compute hi, hs, ti
!     if ice volume is too small
!     energy and water in ice-snow will be transferred into open water
!     negative ice or snow volume will increase the temperature
!     of the open water
!-----------------------------------------------------------------------

            htatm = 0.
            wtatm = 0.
            htocn = 0.
            wtocn = 0.

            Aicet=0.0

            do nc=1,ncat
              if ((A(i,j,idx,nc) .le. ASMALL(nc))
     $             .or.(heff(i,j,idx,nc) .le. 0.0)) then

                Enet=0.0
                do layer=1,nilay(nc)
                  Enet=Enet-En(layer,nc)
                enddo

                htocn=htocn+(Enet+hseff(i,j,idx,nc)*rflsno)
                wtocn=wtocn+(heff(i,j,idx,nc)*rhoice
     $                      +hseff(i,j,idx,nc)*rhosno)

                groice(i,j,nc)=-heff(i,j,idx,nc)
                A(i,j,idx,nc)=0.0
                heff(i,j,idx,nc)=0.0
                hseff(i,j,idx,nc)=0.0
                Ts(i,j,idx,nc)=frzpt(i,j)
                do layer=0,nilay(nc)
                  Ti(layer,nc)=frzpt(i,j)
                enddo
                hi(nc)=HSTAR(nc-1)
                hs(nc)=0.0

              else

                groice(i,j,nc)=0.
                hi(nc)=heff(i,j,idx,nc)/A(i,j,idx,nc)
                hi(nc)=max(HSTAR(nc-1),hi(nc))
                A(i,j,idx,nc)=heff(i,j,idx,nc)/hi(nc)

                do layer=1,nilay(nc)
                  qi(layer,nc)=En(layer,nc)*nilay(nc)/heff(i,j,idx,nc)
                enddo
                call getTmp(Ti(1,nc),qi(1,nc),saltz(1,nc),nilay(nc))
                Ti(0,nc)=Ti(1,nc) ! need a starting value for tstm

                if (hseff(i,j,idx,nc) .le. HSSTAR ) then
                  htocn=htocn+hseff(i,j,idx,nc)*rflsno
                  wtocn=wtocn+hseff(i,j,idx,nc)*rhosno
                  hseff(i,j,idx,nc)=0.0
                  hs(nc)=0.0
                else
                  hs(nc)=hseff(i,j,idx,nc)/A(i,j,idx,nc)
                endif
              endif
              Aicet=Aicet+A(i,j,idx,nc)
            enddo

!-----------------------------------------------------------------------
!     finally see if aicet is greater than 1, it never should
!     be after all the trouble above
!-----------------------------------------------------------------------
            if (Aicet.gt.(1.0+TINY)) then
# if defined O_ice_cpts_warnings
              print*,'warning too much ice area by',Aicet
              print*,nats,i,j,(A(i,j,idx,nc),nc=1,ncat)
# endif
              A(i,j,idx,ncat)=min(1.0,A(i,j,idx,ncat))
              Aleft=1.0-A(i,j,idx,ncat)
              do nc=(ncat-1),1,-1
                A(i,j,idx,nc)=min(A(i,j,idx,nc),Aleft)
                Aleft=Aleft-A(i,j,idx,nc)
              enddo
            endif

            A0=1.0
            do nc=1,ncat
              A0=A0-A(i,j,idx,nc)
            enddo

!-----------------------------------------------------------------------
!     heat conduction in ice and change thickness from growth/melt
!     this is only section that needs Ti
!-----------------------------------------------------------------------

            sblmt = 0.
            do nc=1,ncat
              if (A(i,j,idx,nc).ge.TINY) then
!               if snow is less than 25 cm linearly reduce to ice albedo
                as = min(hseff(i,j,idx,nc)*0.04, c1)
                ca = ice_calb*(c1 - as) + sno_calb*as
# if defined O_sulphate_data_transient
                dswr = solins(i,j)*sbc(i,j,iaca)*pass
     &                 *max(c0, ca - sulph(i,j,2))
# else
                dswr = solins(i,j)*sbc(i,j,iaca)*pass*ca
# endif
                delhit=0.
                delhib=0.
                delhs=0.
                subi=0.
                subs=0.
                delf=0.
                delfs=0.

                fracsnow = hs(nc)/(hs(nc) + 0.1*centi) ! ~ snow param
                Io=dswr(nc)*0.3*(1-fracsnow) ! solar penetrating surf

                call tstm(tmelz(1,nc),saltz(1,nc),
     $               ts(i,j,idx,nc),ti(0,nc),tocn,tbot(nc),frzpt(i,j),
     $               hi(nc),hs(nc),tair,Qair,Flwatm,
     $               Io,Ib,ug,A(i,j,idx,nc),
     $               dswr(nc),ultnt(nc),usens(nc),ulwr(nc),
     $               nilay(nc),Fneti,condb,i,j)

                call dh(saltz(1,nc),Fx,Fbot,Fneti,ultnt(nc),condb,
     $               tocn,ts(i,j,idx,nc),ti(0,nc),tbot(nc),frzpt(i,j),
     $               qi(1,nc),hi(nc),hs(nc),delhib,
     $               subi,subs,
     $               delhit,delhs,delf,delfs,nilay(nc))

!     note that Ti is no longer current because of layer adjustment

                hi(nc)=hi(nc)+delhit+delhib+subi+delf
                hs(nc)=hs(nc)+delhs+subs-delfs

                groice(i,j,nc)=groice(i,j,nc)+(delhit+delhib)
     $                        *A(i,j,idx,nc)

!     atm gains lwe from sublimation
                evp(nc) = -(rhosno*subs+rhoice*subi)
!     ocn gains lwe from melt/growth
# if defined O_convect_brine
                if (delhib .le. 0.) then
                  wff = -(rhosno*delhs+rhoice*(delhit+delhib))
                else
                  wff = -(rhosno*delhs + rhoice*delhit)
                  if (addflxa) then
                    cbf(i,j,nc) = cbf(i,j,nc)
     $                          - A(i,j,idx,nc)*rhoice*delhib
                    cba(i,j,nc) = cba(i,j,nc) + A(i,j,idx,nc)*dtau
                   endif
                endif
# else
                wff = -(rhosno*delhs+rhoice*(delhit+delhib))
# endif
!     heat flux from atmosphere to ice
                Fice = dswr(nc)-ultnt(nc)-usens(nc)-ulwr(nc)

                htatm = htatm + Fice * dtau * A(i,j,idx,nc)
                sblmt = sblmt + (ultnt(nc) * dtau +
     $                  (rvlsno*subs+rvlice*subi) ) * A(i,j,idx,nc)
                wtatm = wtatm + evp(nc) * A(i,j,idx,nc)
                htocn = htocn + (Fx-Ib)* dtau * A(i,j,idx,nc)
                wtocn = wtocn + wff * A(i,j,idx,nc)

              else
!     give a reasonable value for ts, even if there is no ice
                ts(i,j,idx,nc)=frzpt(i,j)
                dswr(nc) = 0.
                ultnt(nc) = 0.
                usens(nc) = 0.
                ulwr(nc) = 0.
                evp(nc) = 0.
              endif

            enddo               ! loop over bins
            dnswr(i,j)  = 0.
            upltnt(i,j) = 0.
            upsens(i,j) = 0.
            uplwr(i,j)  = 0.
            evap(i,j)   = 0.
            do nc=1,ncat
              dnswr(i,j)  = dnswr(i,j)  + A(i,j,idx,nc)*dswr(nc)
              upltnt(i,j) = upltnt(i,j) + A(i,j,idx,nc)*ultnt(nc)
              upsens(i,j) = upsens(i,j) + A(i,j,idx,nc)*usens(nc)
              uplwr(i,j)  = uplwr(i,j)  + A(i,j,idx,nc)*ulwr(nc)
              evap(i,j)   = evap(i,j)   + A(i,j,idx,nc)*evp(nc)/dtau
            enddo

!-----------------------------------------------------------------------
!     recalc heff hseff and energy of melting for each layer
!-----------------------------------------------------------------------

            do nc=1,ncat
              heff(I,J,idx,nc)=hi(nc)*A(I,J,idx,nc)
              hseff(I,J,idx,nc)=hs(nc)*A(I,J,idx,nc)
              do layer=1,nilay(nc)
                En(layer,nc)=qi(layer,nc)*heff(i,j,idx,nc)/nilay(nc)
              enddo
            enddo

!-----------------------------------------------------------------------
!     ice growth in lead and upper layer water temp
!-----------------------------------------------------------------------
            hnew=0.0
            Tnew=frzpt(i,j)

            if (A0.le.TINY) then                !   do nothing
            elseif ((Focean+Fbot).lt.0.0) then  !   grow new ice
              call grownew(fsh_new,Io_new,Flwatm,Focean,Fbot,
     $             hnew,Tnew,tair,Qair,ug,
     $             dtau,fice,frzpt(i,j),evp(0),
     $             dswr(0),ultnt(0),usens(0),ulwr(0))

              Focean=0.5*(Focean+Fice) ! avg. flx w/ and w/o new ice

              htatm = htatm + Focean*A0*dtau
              htocn = htocn + A0*Fbot*dtau
# if defined O_convect_brine
              if (hnew .le. 0.) then
                wtocn = wtocn - A0*hnew*rhoice
              else
                if (addflxa) then
                  cbf(i,j,0) = cbf(i,j,0) - A0*hnew*rhoice
                  cba(i,j,0) = cba(i,j,0) + A0*dtau
                endif
              endif
# else
              wtocn = wtocn - A0*hnew*rhoice
# endif
              groice(i,j,1)=groice(i,j,1)+hnew*A0

            elseif ((A0.lt.1.0)) then           ! lateral ice melt
!     melt ice lateral
              mrate=m1*(max(tocn-frzpt(i,j),0.0))**m2
              arate=pi_eta*mrate/wclead
              arate=min(arate,1./dtau)
!     concentration that melts laterally per ice conc
              DAice=dtau*arate
              Vfrac=1.0-DAice
              Flat=0.
              Dhlatm=0.
              Dslatm=0.

              do nc=1,ncat
                if (A(i,j,idx,nc).ge.TINY) then

!     diagnostic for heat flux, negative when lateral melt
                  Esum=0.
                  do layer=1,nilay(nc)
                    Esum=Esum+En(layer,nc)
                  enddo
                  Esum=Esum-rflsno*hseff(i,j,idx,nc)
                  Flat=Flat+Esum*arate

!     diagnostic for water flux, positive when lateral melt
                  Dhlatm=Dhlatm+heff(i,j,idx,nc)
                  Dslatm=Dslatm+hseff(i,j,idx,nc)

!     update the ice state
                  heff(i,j,idx,nc)=heff(i,j,idx,nc)*Vfrac
                  hseff(i,j,idx,nc)=hseff(i,j,idx,nc)*Vfrac
                  A(i,j,idx,nc)=A(i,j,idx,nc)*Vfrac
                  do layer=1,nilay(nc)
                    En(layer,nc)=En(layer,nc)*Vfrac
                  enddo
                  groice(i,j,nc)=groice(i,j,nc)-heff(i,j,idx,nc)*DAice
                endif
              enddo
              Dhlatm=Dhlatm*DAice
              Dslatm=Dslatm*DAice

              htocn = htocn - Flat*dtau
              wtocn = wtocn + (Dhlatm*rhoice+Dslatm*rhosno)
            endif

!-----------------------------------------------------------------------
!     budgets:   ice-ocn fluxes of heat/water for coupled model
!     subtract htatm and add wtatm to compensate
!     for adding them back later in embm
!-----------------------------------------------------------------------
! do not move this block
            if (addflxa) then
              flux(i,j,isat) = flux(i,j,isat) - (htocn + htatm)
              flux(i,j,ishum) = flux(i,j,ishum) + (wtatm + wtocn)
            endif

            dnswr(i,j)  = dnswr(i,j)  + dswr(0)*A0
            upltnt(i,j) = upltnt(i,j) + ultnt(0)*A0
            upsens(i,j) = upsens(i,j) + usens(0)*A0
            uplwr(i,j)  = uplwr(i,j)  + ulwr(0)*A0
            evap(i,j)   = evap(i,j)   + evp(0)*A0

!-----------------------------------------------------------------------
!     add lead ice to first category
!-----------------------------------------------------------------------
!     compute the energy in the new ice,
!     divide it into nilay(1) equal portions
            A0=1                       ! must recompute after meltlat
            do nc=1,ncat
              A0=A0-A(i,j,idx,nc)
            enddo
            H0prime=max(hstar(0),hnew) ! new ice h, made >=HOSTAR
            A0prime=A0*hnew/H0prime    ! area covered by H0prime ice

            A1prime=A0prime+A(i,j,idx,1)
            h1prime=0.0
            if (A1prime.gt.0.0) then
              h1prime=(hnew*A0+heff(i,j,idx,1))/A1prime
              Ts(i,j,idx,1)=(Ts(i,j,idx,1)*A(i,j,idx,1)+Tnew*A0prime)
     $                     /A1prime
              A(i,j,idx,1)=A1prime
              hi(1)=h1prime
              heff(i,j,idx,1)=hi(1)*A(i,j,idx,1)
              Ad0=A0*hnew/nilay(1)
              Tmp=0.5*(Tnew+frzpt(i,j))
              E0=Ad0*energ(Tmp,salnew)

              epsilon=Ad0*(energ(Tnew,salnew)-energ(frzpt(i,j),salnew))
              epsilon=0.0
              E0=E0+epsilon/2.

!     sum the energies
              do layer=1,nilay(1)
                En(layer,1)=E0+En(layer,1)-epsilon*(layer-0.5)/nilay(1)
              enddo
            endif               ! finished adding new ice to category 1

!     make sure thickness of category 1 is at least hstar(0)
            heff(i,j,idx,1)=A(i,j,idx,1)*hi(1)
            if (hi(1).gt.hstar(0)) then
              A(i,j,idx,1)=heff(i,j,idx,1)/hi(1)
            else
              A(i,j,idx,1)=heff(i,j,idx,1)/hstar(0)
              hi(1)=hstar(0)
            endif

# if defined O_ice_cpts3 || defined O_ice_cpts5 || defined O_ice_cpts10
!-----------------------------------------------------------------------
!     redistribute the ice categories
!-----------------------------------------------------------------------
!     see if category 1 has outgrown its boundary
!     if so, add it to category 2, and so forth
            do nc=1,(ncat-1)    ! only does loop if ncat>1
              if (hi(nc).ge.hstar(nc)) then
                call moveup(frzpt(i,j),heff(i,j,idx,nc),
     $               hseff(i,j,idx,nc),
     $               A(i,j,idx,nc),hs(nc),hi(nc),hstar(nc-1),
     $               esdum1,Ts(i,j,idx,nc),nilay(nc),
     $               heff(i,j,idx,nc+1),hseff(i,j,idx,nc+1),
     $               A(i,j,idx,nc+1),hs(nc+1),hi(nc+1),hstar(nc),
     $               esdum2,Ts(i,j,idx,nc+1),nilay(nc+1),En(1,nc),
     $               En(1,nc+1))
              endif
            enddo

!     see if the thickest category has become too thin (and so forth)
            do nc=ncat,2,-1
              if (hi(nc).lt.hstar(nc-1)) then
                call movedown(frzpt(i,j),heff(i,j,idx,nc-1),
     $               hseff(i,j,idx,nc-1),
     $               A(i,j,idx,nc-1),hs(nc-1),hi(nc-1),hstar(nc-2),
     $               esdum2,Ts(i,j,idx,nc-1),nilay(nc-1),
     $               heff(i,j,idx,nc),hseff(i,j,idx,nc),A(i,j,idx,nc),
     $               hs(nc),hi(nc),hstar(nc-1),
     $               esdum1,Ts(i,j,idx,nc),nilay(nc),En(1,nc-1),
     $               En(1,nc))
              endif
            enddo

!     cannot hurt to do this again
            heff(i,j,idx,1)=A(i,j,idx,1)*hi(1)
            hi(1)=max(hstar(0),hi(1))
            A(i,j,idx,1)=heff(i,j,idx,1)/hi(1)
# endif

!-----------------------------------------------------------------------
!     reconstruct 3d energy of melting matrix
!-----------------------------------------------------------------------
            l=1
            do nc=1,ncat
              do layer=1,nilay(nc)
                E(i,j,idx,l)=En(layer,nc)
                l=l+1
              enddo
            enddo

!-----------------------------------------------------------------------
!     final energy and water diagnostic
!     could be eliminated to save time
!-----------------------------------------------------------------------
            efinl = 0.
            wfinl = 0.
            do nc = 1,ncat
              do layer = 1,nilay(nc)
                efinl = efinl+En(layer,nc)
              enddo
              efinl = efinl-hseff(i,j,idx,nc)*rflsno
              wfinl = wfinl+(heff(i,j,idx,nc)*rhoice
     $                     +hseff(i,j,idx,nc)*rhosno)
            enddo

            ediff = (efinl-einit)-(htatm+htocn+sblmt)
            wdiff = (wfinl-winit)+(wtatm+wtocn)
# if defined O_ice_cpts_warnings
            if (abs(ediff).gt.maxedif)
     $           print*,'big energy difference in ice model at',i,j,
     $                  ediff*1.e-9,
     $           (htatm+htocn+sblmt)*1.e-9, (efinl-einit)*1.e-9
#  if !defined O_convect_brine
            if (abs(wdiff).gt.maxwdif)
     $           print*,'big water difference in ice model at',i,j,
     $                  wdiff,(wtatm+wtocn), (wfinl-winit)
#  endif
# endif

          endif
        enddo
      enddo

      do nc=1,ncat
        call embmbc (heff(1,1,idx,nc))
        call embmbc (A(1,1,idx,nc))
        call embmbc (Ts(1,1,idx,nc))
        call embmbc (hseff(1,1,idx,nc))
      enddo
      do layer=1,ntilay
        call embmbc (E(1,1,idx,layer))
      enddo

      return
      end

      subroutine tstm(tmz1d, sal1d, ts,ti,tw,tbot,tfrz,
     $     h,hs,Tair,qair,Flo,
     $     Io,Ib,ug,area,
     $     dswr,ultnt,usens,ulwr,
     $     ni,F,condb,i,j)

      implicit none

!     This subroutine calculates the evolution of the ice interior and
!     surface temperature from the heat equation and surface energy
!     balance

!     the albedo is fixed for this calculation
!     if you think you want it to be computed implicitly, then Iabs
!     will have to depend on ts too and that could get nasty
!     I think it would be alot of trouble because the albedo affects
!     boundary conditions too

!     solves the heat equation which is non-linear for saline ice
!     by linearizing and then iterating a set of equations.
!     Scheme is backward Euler giving a tridiagonal
!     set of equations implicit in ti
!     (tried crank-nicholson but it behaves poorly when
!     thin ice has a weird initial temperature profile. This can happen
!     when I redistribute between cats so I always use backward

!     Must have a sufficient amount of snow to solve heat equation in
!     snow hsmin is the minimum depth of snow in order to solve for
!     ti(0) if snow thickness < hsmin then do not change ti(0)

!     The number of equations that must be solved by the tridiagonal
!     solver depends on whether the surface is melting and whether
!     there is snow. Four cases are possible (see variable case below):
!     1 = freezing w/ snow, 2 = freezing w/ no snow,
!     3 = melting w/ snow, and 4 = melting w/ no snow

      include "cpts.h"
      include "thermo.h"
      include "cembm.h"
      include "tmngr.h"

! input variables
      integer ni       ! number of vertical layers in the ice
      integer i,j      ! grid points
      real tmz1d(NMAX) ! melting temp of each layer (C)
      real sal1d(NMAX) ! salinity of each layer (ppt)
      real tfrz        ! freezing temperature of ocean (C)
      real h,hs        ! ice and snow thickness (m)
      real area        ! area of the ice/snow
      real Tw          ! water temp below ice (C)
      real ug          ! geostrophic wind speed (m)
      real qice,qair   ! vapour at ice/snow srfc and in air above
      real tair        ! air temperature (C)
      real Flo         ! long wave down at srfc
      real dswr        ! above srfc net dnwd shortwave, + down (W/m**2)
      real Io          ! solar penetrating top srfc, + down (W/m**2)

! output variables
      real Tbot        ! bottom ice temp (C)
      real condb       ! conductive flux at bottom srfc, + up (W/m**2)
      real Ib          ! solar passing through bottom of ice (W/m**2)
      real F           ! net flux at top srfc including conductive flux
                       ! in ice/snow, + down (W/m**2)
      real usens       ! upward sensible flux, + up (W/m**2)
      real ultnt       ! upward latent flux, + up (W/m**2)
      real ulwr        ! upward longwave flux, + up (W/m**2)

! input/output variables
      real ts          ! srfc temp of snow or ice (C)
      real ti(0:NMAX)  ! snow(0) and ice(1:NMAX) interior temp (C)

! local variables
      real ts_old      ! initial temperature of the ice/snow srfc (C)
      real ti_old(0:NMAX) ! initial temperature in the ice and snow (C)
      real tinext(-1:NMAX) ! incremented ts and ti (C)
      real dti(-1:NMAX) ! incremental changes to ts and ti (C)
      real Ts_kelv     ! srfc temperature in kelvin (K)
      real z           ! vertical coordinate (m)

      real coltnt,cosens ! bulk coef for latent and sensible heat flux

      real condt       ! conductive flux at top srfc, + down (W/m**2)
      real Fo          ! net flux at top srfc excluding conductive flux
                       ! in ice/snow (F-condt), + down (W/m**2)
      real dFo_dt      ! derivative of Fo wrt temperature (W/m**2/K)
      real iru         ! net upward longwave flux, + up (W/m**2)
      real dultnt      ! derivative of ultnt wrt temperature (W/m**2/K)
      real ultnti      ! latent heat from initial temperature (W/m**2)
      real dultnti     ! derivative of ultnti wrt temperature (W/m**2/K)
      logical zult     ! flag that is true if latent heat is zero
      real Iabs(NMAX+1) ! solar absorbed in each layer (W/m**2)
      real absorb      ! sum of Iabs (W/m**2)

      real melts       ! the srfc melting temp (either TSMELT or TMELT)
      real alph, bet   ! parameters for maintaining 2nd order accurate
                       ! diff at boundary

      real tinterf     ! snow and ice interf temp (C)

! a,b,c are vectors that describe the diagonal and off-diagonal elements
! of the matrix [A], such that [A] ti = r
      real a(-1:NMAX),b(-1:NMAX),c(-1:NMAX),d(-1:1),r(-1:NMAX)
      real ki(0:NMAX+1) ! layer conductivity divided by layer thickness
      real zeta(0:NMAX) ! the terms in heat equation independent of ti
      real eta(0:NMAX) ! time step / ice layer thickness / fresh ice cp
      real cpi(1:NMAX) ! time step / ice layer thickness / fresh ice cp
      real dh          ! ice layer thickness
      real dt_dh,dt_hs ! time step / ice or snow thickness
      integer  N       ! number of equations solved by tridiag solver
      integer layer    ! counter for ice layers
      integer ipc      ! counter for iterations of dti
      real errit       ! the absolute value of the maximum dti

      integer case     ! indicates the case of snow vs. no snow and
                       ! melting vs. freezing

      real cond2i,dheat
      logical convrg   ! flag that is true if temperature converges
      logical complt   ! flag that is true if a layer melts internally

      real hsmin
      parameter ( hsmin = 0.01*centi )

      real EXPMAX
      parameter(EXPMAX=10.)

! -------------------------------------------------------------------------
!     setup helpful parameters
! -------------------------------------------------------------------------
      dh    = h/(ni + 1.e-20)
      dt_dh = dtau/(dh + 1.e-20)
      dt_hs = 0.0
      if (hs.gt.hsmin) dt_hs=dtau/(hs + 1.e-20)

      ts_old = ts
      do layer = 0,ni
         ti_old(layer) = ti(layer)
      enddo
      tbot = tfrz     ! min(tw,TMELT)

      call conductiv(sal1d,ki,ti,tbot,hs,dh,ni)

!     the solar radiation absorbed internally
      z = 0.0
      Iabs(1) = 1.0
      Iabs(1) = 1.0
      do layer = 1,ni
        z = z+dh
        if (z.lt.0.10*centi) then ! no absorption in the first 10 cm
          Iabs(layer+1) = 1.0
        elseif (z.lt.0.20*centi) then ! ramp up the absorption
          Iabs(layer+1) = 1.0-(1.0-exp(-kappa*0.20*centi))*
     $         (z/(0.10*centi)-1.0)**2
        else
          Iabs(layer+1) = exp(-kappa*z)
        endif
        Iabs(layer) = Io*(Iabs(layer)-Iabs(layer+1))
      enddo

      coltnt = rslatm*dalt_i*ug
      cosens = 0.94*rcpatm*dalt_i*ug
      qice = QS1*exp(min(AI*ts/(ts-BI+1.e-20),EXPMAX))

      if (qice.lt.qair) then
        zult     = .true.
        qice     = qair
        dultnti  = 0.
        ultnti   = 0.
      else
        zult    = .false.
        dultnti = -coltnt*qice*ai*bi/(ts-bi+1.e-20)**2
        ultnti  = coltnt*(qice-qair)
      endif
      ultnt  = ultnti
      dultnt = dultnti

      usens   = cosens*(ts-tair)
      ts_kelv = ts+TFFRESH
      iru     = ESICE*(Ts_kelv)**4

      dFo_dt  = -cosens-dultnt-4*ESICE*Ts_kelv**3
      Fo      = dswr-Io+Flo-ultnt-usens-iru

      if (hs.gt.hsmin) then
         alph  = 2.0*(2.*hs + dh)/(hs+dh+1.e-20)
         bet   = -2.0*hs*hs/(2.*hs+dh+1.e-20)/(hs+dh+1.e-20)
         melts = tsmelt
      else
         alph  =  3.
         bet   =  -1./3.
         melts = tmelt
      endif

      ts_old = ts
      ts     = min(ts,melts-0.01)   ! absolutely necessary
! -------------------------------------------------------------------------
!     beginning of iterative procedure
! -------------------------------------------------------------------------

      ipc = 1
 1000 continue

! -------------------------------------------------------------------------
!     setup terms that depend on the iterated temperature
! -------------------------------------------------------------------------
      call get_Cpice(sal1d,tmz1d,Cpi,Ti(1),Ti_old(1),ni)
      eta(0) = 0.0
      if (hs.gt.hsmin) eta(0) = dtau/(hs+1.e-20)
      do layer = 1,ni
         eta(layer) = dt_dh/(cpi(layer)+1.e-20)
      enddo

      zeta(0) = rcpsno*Ti_old(0)
      do layer = 1,ni
        zeta(layer) = Ti_old(layer)+eta(layer)*Iabs(layer)
      enddo

      if (ts.lt.melts-tiny) then
!        solve heat equation for ice/snow using flux BC

         if (hs.gt.hsmin) then
! -------------------------------------------------------------------------
            case = 1      ! case of freezing with snow layer
! -------------------------------------------------------------------------

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,ni,1)

            a(0) = -eta(0)*ki(0)*(alph+bet)
            c(0) = eta(0)*(bet*ki(0)-ki(1))
            b(0) = rcpsno+eta(0)*(alph*ki(0)+ki(1))
            r(0) = zeta(0)

            a(-1) = 0.
            c(-1) = ki(0)*alph
            d(-1) = ki(0)*bet
            b(-1) = dFo_dt-c(-1)-d(-1)
            r(-1) = -Fo+dFo_dt*ts_old

!           row operation to get rid of d(-1)
            b(-1) = c(0)*b(-1)-d(-1)*a(0)
            c(-1) = c(0)*c(-1)-d(-1)*b(0)
            r(-1) = c(0)*r(-1)-d(-1)*r(0)

            N = ni+2
            call tridag(a(-1),b(-1),c(-1),r(-1),tinext(-1),N)
            dti(-1) = tinext(-1)-ts
            ts = tinext(-1)
            do layer = 0,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo

!            write(*,2000) 'I. iterations',ipc*1.0,ts,(ti(layer),layer = 0,4)
!     $           ,(dti(layer),layer = -1,4)

         else
! -------------------------------------------------------------------------
            case = 2      ! case of freezing with no snow layer
! -------------------------------------------------------------------------

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,ni,2)

            a(1) = -eta(1)*ki(1)*(alph+bet)
            c(1) = -eta(1)*(ki(2)-bet*ki(1))
            b(1) = 1+eta(1)*(ki(2)+alph*ki(1))
            r(1) = zeta(1)

            a(0) = 0.
            c(0) = ki(1)*alph
            d(0) = ki(1)*bet
            b(0) = dFo_dt-c(0)-d(0)
            r(0) = -Fo+dFo_dt*ts_old

!           row operation to get rid of d(0)
            b(0) = c(1)*b(0)-d(0)*a(1)
            c(0) = c(1)*c(0)-d(0)*b(1)
            r(0) = c(1)*r(0)-d(0)*r(1)

            N = ni+1
            call tridag(a(0),b(0),c(0),r(0),tinext(0),N)

            dti(-1) = 0.
            dti(0)  = tinext(0)-ts
            ts = tinext(0)
            do layer = 1,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo
            ti(0) = ts

         endif
         ts = min(ts,melts)
      else
         ts = melts
         if (hs.gt.hsmin) then
! -------------------------------------------------------------------------
            case = 3      ! case of melting with snow layer
! -------------------------------------------------------------------------

            a(0) = -eta(0)*ki(0)*(alph+bet)
            c(0) = eta(0)*(bet*ki(0)-ki(1))
            b(0) = rcpsno+eta(0)*(alph*ki(0)+ki(1))
            r(0) = zeta(0)-a(0)*ts
            a(0) = 0.0

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,ni,1)

            N = ni+1
            call tridag(a(0),b(0),c(0),r(0),tinext(0),N)

            dti(-1) = 0.
            do layer = 0,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo

!            write(*,2000)'III. T =  ', ipc*1.0,(ti(layer),layer = 0,4)
!     $           ,(dti(layer),layer = 1,4)

         else
! -------------------------------------------------------------------------
            case = 4      ! case of melting with no snow layer
! -------------------------------------------------------------------------

            call getabc(a,b,c,r,ti,tbot,zeta,ki,eta,ni,2)

            a(1) = -eta(1)*ki(1)*(alph+bet)
            c(1) = -eta(1)*(ki(2)-bet*ki(1))
            b(1) = 1+eta(1)*(ki(2)+alph*ki(1))
            r(1) = zeta(1)-a(1)*Ts
            a(1) = 0.0

            N = ni
            call tridag(a(1),b(1),c(1),r(1),tinext(1),N)

            dti(-1) = 0.
            dti(0) = 0.
            do layer = 1,ni
               dti(layer) = tinext(layer)-ti(layer)
               ti(layer)  = tinext(layer)
            enddo

            ti(0) = ts
         endif
      endif

! -------------------------------------------------------------------------
!     end iterative procedure, see if need to reiterate
! -------------------------------------------------------------------------
! -------------------------------------------------------------------------
!     make certain that latent heat is always positive
!     be sure to iterate one more time if latent heat becomes
!     less than zero during the iteration (and was positive before)
! -------------------------------------------------------------------------
      ultnt  =  ultnti + dultnti*(ts-ts_old)
      if (ultnt.lt.0. .and. .not.zult) then
        zult    = .true.
        ultnt   = 0.
        ultnti  = 0.
        dultnti = 0.
        dFo_dt  = -cosens-4*ESICE*Ts_kelv**3
        Fo      = dswr-Io+Flo-usens-iru
        ipc     = ipc+1
        go to 1000
      endif

      errit = 0.0
      do layer = -1,ni
         errit = max(errit,abs(dti(layer)))
      enddo
      if ((errit.lt.errmax).or.(ipc.gt.20)) go to 6000 ! done iterating

      ipc = ipc+1
      go to 1000    ! to beginning of iterative routine

! -------------------------------------------------------------------------
!     continue from here when done iterative
!     begin with error checking
! -------------------------------------------------------------------------
 6000 continue

      complt = .false.
      do layer = 1,ni
         if (ti(layer).gt.tmz1d(layer)) complt = .true.
      enddo

      convrg = .true.
      if (errit.gt.errmax) then
        if (errit.gt.0.005) then
          print*,'WARNING No CONVERGENCE in icemodel',errit,i,j
          convrg = .false.
        else
          if (errit.gt.5.e-4)
     $      print*,'WARNING POOR CONVERGENCE in icemodel',errit,i,j
        endif
      endif

      if (complt .or. .not.convrg) then
        if (complt .and. convrg)
     $       print*,'WARNING CONVERGES to ti>TMELT in icemodel',i,j
        if (complt .and. .not.convrg)
     $       print*,'& final iteration has ti>TMELT'
        print*,'diagnostics useful for debugging'
        print*,area,h,hs,dswr,Flo
        print*,'tstmnew computed T(z):'
        print*,ts,(ti(layer),layer = 0,ni)
        print*,'tstmnew started with T(z):'
        print*,ts_old,(ti_old(layer),layer = 0,ni),tbot

        ts = min(ts_old,melts)
        tinterf=(h*kappas*Ts+hs*kappai*tbot)/(h*kappas+hs*kappai+1.e-20)
        Ti(0) = 0.5*(Ts+tinterf)
        do layer = 1,ni
          Ti(layer) = tinterf+(layer-0.5)*(tbot-tinterf)/(ni+1.e-20)
        enddo
        print*,'setting the temperature profile to be linear T(z):'
        print*,ts,(ti(layer),layer = 0,ni)
      elseif (ti(0).gt.(TSMELT+TINY)) then
        print*,'WARNING SNOW has ti>TSMELT in icemodel',i,j
        print*,'  likely no problem if hs or area is tiny'
        print*,'  area=',area,' hs=',hs
        print*,'setting snow temp to be TSMELT'
        ti(0) = tsmelt
       endif

! -------------------------------------------------------------------------
!     finish up by updating the surface fluxes, etc.
! -------------------------------------------------------------------------

      if (hs.gt.hsmin) then
         condt = ki(0)*(alph*(ti(0)-ts)+bet*(ti(1)-ts))
      else
         condt = ki(1)*(alph*(ti(1)-ts)+bet*(ti(2)-ts))
      endif
      condb = ki(ni+1)*( 3.*(tbot-ti(ni)) - (tbot-ti(ni-1))/3.0  )

      iru   =  iru + (4*ESICE*(ts_old+tffresh)**3)*(Ts-Ts_old)
      ulwr  =  -Flo + iru
      usens =  usens + cosens*(Ts-Ts_old)

      Fo    =  Fo+dFo_dt*(Ts-Ts_old)
      F     =  Fo + condt

!     in the rare event that F<0 set it equal to 0, adjust sensible heat
      if (F .lt. 0.) then
        Fo = -condt
        F  = 0.
        usens = condt+dswr-Io-ulwr-ultnt
      endif

      absorb = 0.
      do layer = 1,ni
         absorb = absorb+Iabs(layer)
      enddo
      Ib = Io-absorb  ! solar passing through the bottom of the ice

! 2000 format(A15, 50(f8.3))

      return
      end

      subroutine getabc(a,b,c,r,ti,tbot,zeta,k,eta,ni,lfirst)

!     compute elements of tridiagonal matrix

      implicit none
      include "cpts.h"
      include "thermo.h"

!  parameters for maintaining 2nd order accurate diff at bottom boundary
      real       alph, bet
      parameter( alph = 3.0 ,bet = -1.0/3.0 )

! input variables
      integer ni              ! number of layers
     $       ,lfirst          ! start with this layer
      real    ti   (0:nmax)   ! temperature of ice-snow layers
     $       ,tbot            ! temperature of ice bottom surface
     $       ,zeta (0:nmax)
     $       ,k    (0:nmax+1) ! ice-snow conductivity
     $       ,eta  (0:nmax)
! output variables
      real    a    (-1:nmax)  ! sub-diagonal elements
     $       ,b    (-1:nmax)  ! diagonal elements
     $       ,c    (-1:nmax)  ! super-diagonal elements
     $       ,r    (-1:nmax)  ! constants (indep. of ti)

! local variable
      integer layer

! if there is snow lfirst = 1 otherwise it is 2
      do layer = lfirst,(ni-1)
        a(layer) = -eta(layer)*k(layer)
        c(layer) = -eta(layer)*k(layer+1)
        b(layer) = 1-c(layer)-a(layer)
        r(layer) = zeta(layer)
      enddo
      a(ni) = -eta(ni)*(k(ni)-bet*k(ni+1))
      c(ni) = 0.0
      b(ni) = 1+eta(ni)*(k(ni)+alph*k(ni+1))
      r(ni) = zeta(ni)+eta(ni)*(alph+bet)*k(ni+1)*tbot

      return
      end

      subroutine tridag(a,b,c,r,u,n)

      implicit none

      include "cpts.h"
      integer nnmax
      parameter (nnmax = nmax+2) ! don't change this

! input variables
      integer n               ! number of layers
      real    a      (nnmax)  ! sub-diagonal elements
     $       ,b      (nnmax)  ! diagonal elements
     $       ,c      (nnmax)  ! super-diagonal elements
     $       ,r      (nnmax)  ! constants (indep. of ti)

! output variables
      real    u      (nnmax)  ! solution

! local variables
      integer layer
      real    bet             ! work variable
     $       ,gam    (nnmax)  ! work array

      bet = b(1)
      u(1) = r(1)/bet
      do layer = 2,n
        gam(layer) = c(layer-1)/bet
        bet = b(layer)-a(layer)*gam(layer)
        u(layer) = (r(layer)-a(layer)*u(layer-1))/bet
      enddo
      do layer = n-1,1,-1
        u(layer) = u(layer)-gam(layer+1)*u(layer+1)
      enddo

      return
      end

      subroutine conductiv(sal1d,ki,ti,tbot,hs,dh,ni)

      implicit none

      include "cpts.h"
      include "thermo.h"

      real sal1d    (NMAX) ! salinity of each layer
      real ki(0:nmax+1),ti(0:nmax),tbot,hs,dh
      real tmax
      parameter (tmax = -0.1)
      integer ni,layer

      real hsmin
      parameter ( hsmin = 0.01*centi )

      do layer = 2,ni
         ki(layer) = kappai+BETA*(sal1d(layer)+sal1d(layer+1))/2.
     $        /min(tmax,(ti(layer-1)+ti(layer))/2.)
         ki(layer) = max(ki(layer),KIMIN)
         ki(layer) = ki(layer)/dh
      enddo
      ki(ni+1) = kappai+BETA*sal1d(ni+1)/tbot
      ki(ni+1) = max(ki(ni+1),KIMIN)
      ki(ni+1) = ki(ni+1)/dh
      ki(1) = kappai+BETA*sal1d(1)/min(tmax,ti(1))
      ki(1) = max(ki(1),KIMIN)
      ki(1) = ki(1)/dh
!cc      ki(1) = kappai/dh
      if (hs.gt.hsmin) then
         ki(0) = kappas/hs
         ki(1) = 2.0*ki(1)*ki(0)/(ki(1)+ki(0))
      else
         ki(0) = 0.0
      endif

      return
      end

      subroutine get_Cpice(sal1d,tmz1d,cpi,Tnew,Told,ni)

!     compute the heat capacity of each layer in the ICE
!     Tnew is assumed to be the midpoint temperature of each layer

      implicit none

      include "cpts.h"
      include "thermo.h"

      real sal1d    (NMAX) ! salinity of each layer
      real tmz1d    (NMAX) ! melting temp of each layer
      real cpi(nmax),Tnew(nmax),Told(nmax)
      integer ni,layer

      do layer = 1,ni
        cpi(layer) = rcpice+gamma*sal1d(layer)/Told(layer)/
     $       min(Tnew(layer),tmz1d(layer))
      enddo

#endif

      return
      end
