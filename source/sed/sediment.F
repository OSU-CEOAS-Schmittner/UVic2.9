!=======================================================================
      subroutine bury (n_control, zsed, delz, pore, kmax, calgg, orggg
     &,                sed_ml_mass, rain_cal_p, ttrcal, c_advect
     &,                rain_clay_p, nyear, buried_mass, buried_calfrac
     &,                depth_age, ipsed, ipmax, nzmax, ibmax)

#if defined O_sed
      implicit none

      integer n_control, ibmax, ip, ipmax, ipsed, k
      integer kmax, n_coretop, n_depth, nyear, nzmax, n_accum

!     arguments
      real zsed(nzmax), delz(nzmax), pore(nzmax,ipmax)
      real calgg(nzmax,ipmax), orggg(nzmax,ipmax)
      real sed_ml_mass(ipmax), rain_cal_p(ipmax), ttrcal(ipmax)
      real c_advect(ipmax), rain_clay_p(ipmax)

!     results
      real buried_mass(ibmax,ipmax), buried_calfrac(ibmax,ipmax)
      real depth_age(ibmax)

!     internal variables (temporary)
      real buried_mass_step(ipmax), calgg_est(ipmax)
      real sed_ml_mass_new(ipmax), up_cal_frac, cal_total_new

!!     uncomment for debugging
!      real calml(ipmax), calml_new(ipmax), calbur(ipmax)
!      real calbur_new(ipmax), caltot(ipmax), caltot_new(ipmax), tmp

!!     uncomment for debugging, add up total calcite
!      do ip=1,ipsed
!        calml(ip) = sed_ml_mass(ip)*calgg(kmax,ip)
!        calbur(ip) = 0.
!        do k=1,ibmax
!          calbur(ip) = calbur(ip) + buried_mass(k,ip)
!     &                 *buried_calfrac(k,ip)
!        enddo
!        caltot(ip) = calml(ip) + calbur(ip)
!      enddo

!     ==> do burial

!     estimate the change in calgg, then calculate new porosities
      do ip=1,ipsed
        calgg_est(ip) = (calgg(kmax,ip)*sed_ml_mass(ip)
     &                + rain_cal_p(ip)*100. - ttrcal(ip)*100.)
     &                /sed_ml_mass(ip)
      enddo

      call set_est_pore (calgg_est, zsed, pore, kmax, ipsed
     &,                 ipmax, nzmax)

      call get_sed_ml_mass (delz, pore, sed_ml_mass_new, nzmax
     &,                     ipsed, ipmax, kmax)

      n_depth = 2
      do k=ibmax-2,1,-1
        if (float(nyear) - depth_age(k) .gt. -1.e-6) n_depth = k + 1
      enddo

!     calculate mass buried this time step
!     units of g/cm2 per this time step
      do ip=1,ipsed
        buried_mass_step(ip) = rain_cal_p(ip)*100.
     &                       - ttrcal(ip)*100. + rain_clay_p(ip)
     &                       + sed_ml_mass(ip) - sed_ml_mass_new(ip)
      enddo

!     now do some burial
!     this loop won't parallelize very well.  too bad.
      do ip=1,ipsed

        do k=1,n_depth
          if (buried_mass(k,ip) .gt. 0.) n_coretop = k
        enddo

!       find the "latest" box with mass in it, one ip at a time
        if (buried_mass_step(ip) .gt. 0.) then
!         then we're accumulating -- that's easy

!         add buried calcite to the sediment record
          buried_calfrac(n_depth,ip) = (buried_mass(n_depth,ip)
     &                               *buried_calfrac(n_depth,ip)
     &                               + buried_mass_step(ip)
     &                               *calgg(kmax,ip))
     &                               /(buried_mass(n_depth,ip)
     &                               + buried_mass_step(ip))

!         add buried mass to sediment record
          buried_mass(n_depth,ip) = buried_mass(n_depth,ip)
     &                            + buried_mass_step(ip)

!         adjust calgg
          do k=1,kmax
            calgg(k,ip) = (calgg(k,ip)*sed_ml_mass(ip)
!                       [ g calcite (old) in ml ]
     &                  + rain_cal_p(ip)*100. - ttrcal(ip)*100.
     &                  - buried_mass_step(ip)*calgg(kmax,ip))
!                       [ minus calcite buried ]
     &                  /sed_ml_mass_new(ip)
          enddo

!         all business is complete at this point

        else
!         or else there's chemical erosion -- bms < 0

          if (-buried_mass_step(ip) .gt. buried_mass(n_coretop,ip))
     &      then

!           then we're using up the entire coretop box by erosion
!           assume that this can only happen once per year timestep.

!           calculate the fraction of calcite in the upwelling
!           material by combining the fractions of boxes n_coretop
!           (until it's used up) and n_coretop - 1 (below top)

!           write(6,*) "into box depletion in bury"

            up_cal_frac = (buried_mass(n_coretop,ip)
     &                  *buried_calfrac(n_coretop,ip)
     &                  + ((-buried_mass_step(ip)
     &                  - buried_mass(n_coretop,ip))
     &                  *buried_calfrac(n_coretop-1,ip)))
     &                  /(-buried_mass_step(ip))

!           adjust the sediment record to reflect erosion

            buried_mass(n_coretop-1,ip) = buried_mass(n_coretop-1,ip)
     &                                  + (buried_mass_step(ip)
     &                                  + buried_mass(n_coretop, ip))
            buried_mass(n_coretop,ip) = 0.
            n_coretop = n_coretop - 1

          else
!           or else erosion can come entirely out of the box
!           n_coretop, and the upwelling calcite fraction is
!           simply  buried_calfrac(n_coretop,ip)

            up_cal_frac = buried_calfrac(n_coretop,ip)

!           adjust the sediment record

            buried_mass(n_coretop,ip) = buried_mass(n_coretop,ip)
     &                                + buried_mass_step(ip)

          endif

!         now gotta adjust calgg
          cal_total_new = sed_ml_mass(ip)*calgg(kmax,ip)
     &                  + rain_cal_p(ip)*100. - ttrcal(ip)*100.
     &                  - buried_mass_step(ip)*up_cal_frac

          do k=1,kmax
            calgg(k,ip) = cal_total_new/sed_ml_mass_new(ip)
          enddo

        endif
!       end of the "erosion if"

        do k=1,kmax
          if (calgg(k,ip) .lt. 0.) then
            calgg(k,ip) = 1.e-6
            write(6,*) "calgg reached 0 in ip", ip
          endif
        enddo

!!       uncomment for debugging, add up new total calcite
!        calml_new(ip) = sed_ml_mass_new(ip)*calgg(kmax,ip)
!        calbur_new(ip) = 0.
!        do k=1,ibmax
!          calbur_new(ip) = calbur_new(ip) + buried_mass(k,ip)
!     &                   *buried_calfrac(k,ip)
!        enddo
!        caltot_new(ip) = calml_new(ip) + calbur_new(ip)

      enddo

!!     uncomment for debugging
!      do ip=1,ipsed
!        tmp = caltot(ip) + (rain_cal_p(ip) - ttrcal(ip))*100.
!     &      - caltot_new(ip)
!        if (abs(tmp) .gt. 1.e-8) then
!          write(6,*) "conservation error in bury ", ip
!     &,     caltot(ip), rain_cal_p(ip)*100., ttrcal(ip)*100.
!     &,     caltot_new(ip), calml_new(ip), calbur_new(ip)
!     &,     abs(tmp)
!        endif
!      enddo

      do ip=1,ipsed
        sed_ml_mass(ip) = sed_ml_mass_new(ip)
      enddo

      return
      end

!=======================================================================
      subroutine set_pore (calgg, z, pore, kmax, ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipmax, ipsed, k, kmax, nzmax

      real pore(nzmax,ipmax), calgg(nzmax,ipmax), z(nzmax), pore_max
      real exp_pore

      do ip=1,ipsed
!       pore_max = 0.65*calgg(kmax,ip) + 0.85*(1. - calgg(kmax,ip))
        pore_max = 1 - (0.483 + 0.45*calgg(kmax,ip))/2.5
        exp_pore = 0.25*calgg(kmax,ip) + 3.*(1. - calgg(kmax,ip))
        do k=2,kmax
          pore(k,ip) = exp(-z(k)/exp_pore)*(1.-pore_max) + pore_max
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine set_est_pore (calgg, z, pore, kmax, ipsed, ipmax
     &,                        nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipmax, ipsed, k, kmax, nzmax

      real pore (nzmax, ipmax), calgg(ipmax), z(nzmax),  pore_max
      real exp_pore

      do ip=1,ipsed
!       pore_max = 0.65*calgg(ip) + 0.85*(1. - calgg(ip))
        pore_max = 1 - (0.483 + 0.45*calgg(ip))/2.5
        exp_pore = 0.25*calgg(ip) + 3.*(1-calgg(ip))
        do k=2,kmax
          pore(k,ip) = exp(-z(k)/exp_pore)*(1.-pore_max) + pore_max
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine get_sed_ml_mass (delz, pore, sed_ml_mass, nzmax, ipsed
     &,                           ipmax, kmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipmax, ipsed, k, kmax, nzmax

      real delz(nzmax), pore(nzmax,ipmax), sed_ml_mass(ipmax)

      do ip=1,ipsed
        sed_ml_mass(ip) = 0.
        do k=1,kmax
          sed_ml_mass(ip) = sed_ml_mass(ip)+delz(k)*(1.-pore(k,ip))*2.5
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine estimate_rc (rain_org, rc, ipsed)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed

      real rain_org(ipsed), rc(ipsed)

      do ip=1,ipsed
        rc(ip) = 2.e-9
      enddo

      return
      end

!=======================================================================
      subroutine setup_pw (co2, hco3, co3, o2_bw, carb, o2, csat, kmax
     &,                   ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipmax, ipsed, kmax, nzmax

      real co2(ipmax), hco3(ipmax), co3(ipmax), o2_bw(ipmax)
      real carb(nzmax,3,ipmax), o2(nzmax,ipmax), csat(ipmax)

CDIR& NOVECTOR
      do ip=1,ipsed
        o2(1,ip) = o2_bw(ip)
        carb(1,1,ip) = co2(ip)
        carb(1,2,ip) = hco3(ip)
        carb(1,3,ip) = co3(ip)
      enddo

      return
      end

!=======================================================================
      subroutine reset_pw (carb, csat, kmax, ip, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipmax, ipsed, k, kmax, nzmax

      real carb(nzmax,3,ipmax), csat(ipmax)

CDIR& NOVECTOR
      carb(2,3,ip) = (2.*carb(1,3,ip) + (csat(ip) - 15.e-6))/3.
      carb(3,3,ip) = (carb(1,3,ip) + 2.*(csat(ip) - 15.e-6))/3.

      do k=4,kmax
        carb(k,3,ip) = csat(ip) - 10.e-6
      enddo

      do k=2,kmax
        carb(k,1,ip) = carb(1,1,ip)
        carb(k,2,ip) = carb(1,2,ip)
      enddo

      return
      end

!=======================================================================
      subroutine sed_ss (zsed, delz, form, pore, kmax, o2, zrct, carb
     &,                  orgml, orggg, calml, calgg, raincal, rainorg
     &,                  rainclay, r13or, r14or, r13ca, r14ca, rc
     &,                  dissc, dissn, csat, u1, u2, dopls, domin
     &,                  dcpls, dcmin, dbpls, dbmin, resp_c, cal_c
     &,                  ttrorg, ttrcal, ttrtc, ttral, diftc, difal
     &,                  c_advect, buried_mass, buried_calfrac
     &,                  n_buried_depths, ipsed, ipmax, nzmax, ibmax
     &,                  n_control, n_debug)

!-----------------------------------------------------------------------
!  finds the steady state calcite concentration in the sediments
!-----------------------------------------------------------------------

      implicit none

      integer ibmax, ipsed, ipmax, kmax, nzmax, n_control, n_debug

!     arguments
      real zsed(nzmax),delz(nzmax), form(nzmax,ipmax), pore(nzmax,ipmax)
      real o2(nzmax,ipmax), zrct(ipmax), carb(nzmax,3,ipmax)
      real orgml(nzmax,ipmax), calml(nzmax,ipmax), orggg(nzmax,ipmax)
      real calgg(nzmax,ipmax), raincal(ipmax), rainorg(ipmax)
      real rainclay(ipmax), r13or, r14or, r13ca, r14ca, rc(ipmax)
      real dissc, dissn, csat(ipmax), u1(ipmax), u2(ipmax)
      real dopls(nzmax,ipmax), domin(nzmax,ipmax), dcpls(nzmax,3,ipmax)
      real dcmin(nzmax,3,ipmax), dbpls(nzmax,ipmax), dbmin(nzmax,ipmax)
      real resp_c(nzmax,3,ipmax), cal_c(nzmax,ipmax), expb, db, difo2
      real difc(3)

!     results
      real ttrorg(ipmax),ttrcal(ipmax), ttrtc(ipmax),ttral(ipmax)
      real difal(ipmax),diftc(ipmax), c_advect(ipmax)
      real buried_mass(ibmax,ipmax), buried_calfrac(ibmax,ipmax)

      integer n_buried_depths(ipmax)
CDIR& NOVECTOR

      call o2org (rainorg, rc, kmax, zsed, delz, form, pore, dopls
     &,           domin, dbpls, dbmin, o2, zrct, orgml, orggg, resp_c
     &,           ipsed, ipmax, nzmax, n_control)

      call calss (raincal, rainorg, rainclay, resp_c, cal_c, dissc
     &,          dissn, csat, u1, u2, dcpls, dcmin, zsed, delz, form
     &,          pore, kmax, calml, calgg, carb, ttrorg, ttrcal, ttral
     &,          ttrtc, difal, diftc, c_advect, ipsed, ipmax, nzmax
     &,          n_control, n_debug)

      expb = 3.0
      difo2 = 12.1e-6
      difc(1) = 10.5e-6
      difc(2) = 6.4e-6
      difc(3) = 5.2e-6
      db = 0.15

      call set_pore (calgg, zsed, pore, kmax, ipsed, ipmax, nzmax)

      call pore_2_form (pore, form, kmax, expb, ipsed, ipmax, nzmax)

      call calc_do2 (difo2, form, pore, delz, kmax, dopls, domin, ipsed
     &,             ipmax, nzmax)
      call calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin, ipsed
     &,            ipmax, nzmax)
      call calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin, ipsed
     &,            ipmax, nzmax)

      return
      end

!=======================================================================
      subroutine sed_const_cal (zsed, delz, form, pore, kmax, o2, zrct
     &,                         carb, orgml, orggg, calml, calgg
     &,                         raincal, rainorg, rainclay, r13or
     &,                         r14or, r13ca, r14ca, rc, dissc, dissn
     &,                         csat, u1, u2, dopls, domin, dcpls
     &,                         dcmin, dbpls, dbmin, resp_c, cal_c
     &,                         ttrorg, ttrcal, ttrtc, ttral, diftc
     &,                         difal, c_advect, ipsed, ipmax, nzmax
     &,                         n_control, n_debug)

!-----------------------------------------------------------------------
!  finds the steady state pore water chemistry with const calcite
!-----------------------------------------------------------------------

      implicit none

      integer inner_loop_limit, ip, ipsed, ipmax, iter
      integer k, kmax, nzmax, n_control, n_debug

!   arguments
      real zsed(nzmax),delz(nzmax), form(nzmax,ipmax), pore(nzmax,ipmax)
      real o2(nzmax,ipmax), zrct(ipmax), carb(nzmax,3,ipmax)
      real orgml(nzmax,ipmax), calml(nzmax,ipmax), orggg(nzmax,ipmax)
      real calgg(nzmax,ipmax), raincal(ipmax),rainorg(ipmax)
      real rainclay(ipmax), r13or, r14or, r13ca, r14ca, rc(ipmax)
      real dissc, dissn, csat(ipmax), u1(ipmax), u2(ipmax)
      real dopls(nzmax,ipmax), domin(nzmax,ipmax), dcpls(nzmax,3,ipmax)
      real dcmin(nzmax,3,ipmax), dbpls(nzmax,ipmax), dbmin(nzmax,ipmax)
      real resp_c(nzmax,3,ipmax), cal_c(nzmax,ipmax)
      real raincal_cutoff

!   results
      real ttrorg(ipmax), ttrcal(ipmax), ttrtc(ipmax), ttral(ipmax)
      real difal(ipmax),diftc(ipmax), c_advect(ipmax)

!   internal variables
      integer still_wrong(ipmax), sum_wrong
      real expb, db, difo2, difc(3)

CDIR& NOVECTOR

      expb = 3.0
      difo2 = 12.1e-6
      difc(1) = 10.5e-6
      difc(2) = 6.4e-6
      difc(3) = 5.2e-6
      db = 0.15

      call set_pore (calgg, zsed, pore, kmax, ipsed, ipmax, nzmax)

      call pore_2_form (pore, form, kmax, expb, ipsed, ipmax, nzmax)

      call calc_do2 (difo2, form, pore, delz, kmax, dopls, domin, ipsed
     &,             ipmax, nzmax)
      call calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin, ipsed
     &,            ipmax, nzmax)
      call calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin, ipsed
     &,             ipmax, nzmax)
      call o2org (rainorg, rc, kmax, zsed, delz, form, pore, dopls
     &,           domin, dbpls, dbmin, o2, zrct, orgml, orggg, resp_c
     &,           ipsed, ipmax, nzmax, n_control)

      raincal_cutoff = 0.1e-6

      sum_wrong = 0
      do ip=1,ipsed
        if (raincal(ip) .gt. raincal_cutoff) then
          sum_wrong = sum_wrong + 1
          still_wrong(sum_wrong) = ip
        else

!!         uncomment for debugging
!          if (n_debug .ge. 2) write(6,*) "cutting off site ", ip

          ttrcal(ip) = raincal(ip)
          do k=2,kmax
            calgg(k,ip) = 0.
          enddo
        endif
      enddo

!!     uncomment for debugging
!      if (n_debug .gt. 2) write(6,*) "end"

      inner_loop_limit = 200
      call co3ss (resp_c, cal_c, dissc, dissn, csat, u1, u2, zsed, delz
     &,           form, pore, kmax, dcpls, dcmin, calml, calgg, carb
     &,           ttrorg, ttrcal, ttral, ttrtc, difal, diftc, ipsed
     &,           ipmax, nzmax, inner_loop_limit, n_debug, still_wrong
     &,           sum_wrong, ipmax, iter)

      do ip=1,ipsed
        if ((ttrcal(ip).gt.raincal(ip)) .and. (calgg(kmax,ip).lt.0.001))
     &    ttrcal(ip) = raincal(ip)
        c_advect(ip) = raincal(ip) - ttrcal(ip)
      enddo

      return
      end

!=======================================================================
      subroutine calc_k (temp, sal, z, k1, k2, kb, csat, ipsed)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed

      real kprime, delv, rr, dk

      parameter (kprime = 4.75e-7) ! mol2 / kg2
      parameter (delv = -44.)      ! cm3 / mol
      parameter (rr = 83.14)       ! cm3 bar / k mol
      parameter (dk = -.0133)      ! cm3 / bar mol

!     arguments
      real temp(ipsed), sal(ipsed), z(ipsed), k1(ipsed), k2(ipsed)
      real kb(ipsed), csat(ipsed)

!     local variables
      real tk(ipsed), cp(ipsed), prat(ipsed), pres(ipsed), kpres(ipsed)

CDIR& NOVECTOR

      do ip=1,ipsed

        tk(ip) = temp(ip) + 273.15

!       k1 and k2 (apparent), from mehrbach

        k1(ip) = 13.7201 - 0.031334*tk(ip) - 3235.76/tk(ip)
     &         - 1.3e-5*sal(ip)*tk(ip) + 0.1032*sal(ip)**(0.5)

        k1(ip) = 10**(k1(ip))
        cp(ip) = (z(ip)/10.)/83.143/tk(ip)
        prat(ip) = (24.2 - 0.085*temp(ip))*cp(ip)
        prat(ip) = exp(prat(ip))
        k1(ip) = k1(ip)*prat(ip)

        k2(ip) = - 5371.9645 - 1.671221*tk(ip) + 128375.28/tk(ip)
     &         + 2194.3055*log(tk(ip))/2.30259 - 0.22913*sal(ip)
     &         - 18.3802*log(sal(ip))/2.30259 + 8.0944e-4*sal(ip)*tk(ip)
     &         + 5617.11*log(sal(ip))/tk(ip)/2.30259
     &         - 2.136*sal(ip)/tk(ip)

        k2(ip) = 10**(k2(ip))
        prat(ip) = (16.4 - 0.04*temp(ip))*cp(ip)
        prat(ip) = exp(prat(ip))
        k2(ip) = k2(ip)*prat(ip)

!    lyman's kb
        kb(ip) = 2291.9/(temp(ip) + 273) + 0.01756*(temp(ip) + 273)
     &         - 3.385 - .32051*(sal(ip)/1.80655)**(1./3.)
        kb(ip) = 10**(-kb(ip))
        prat(ip) = (27.5 - 0.095*temp(ip))*cp(ip)
        prat(ip) = exp(prat(ip))
        kb(ip) = kb(ip)*prat(ip)

!    calculate sayles calcite saturation state
        pres(ip) = z(ip)/10 ! bar
        kpres(ip) = LOG(4.75e-7) - delv/(rr*(temp(ip)+273.))*(pres(ip))
     &            + 0.5*dk/(rr*(temp(ip) + 273.))*(pres(ip))**2
        kpres(ip) = eXP(kpres(ip))
        csat(ip) = kpres(ip)/0.01

      enddo

      return
      end

!=======================================================================
      subroutine calc_buff (alk, tco2, sal, k1, k2, kb, co2, hco3, co3
     &,                     ipsed)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer icnt, ip, ipsed

      real alk(ipsed), tco2(ipsed), sal(ipsed), k1(ipsed), k2(ipsed)
      real kb(ipsed), co2(ipsed), hco3(ipsed), co3(ipsed)

      real tbor(ipsed), tkt(ipsed), tk(ipsed), c1(ipsed), c2(ipsed)
      real c4(ipsed), a(ipsed), x(ipsed), ah1(ipsed), aht(ipsed)

!     all units of moles/l
CDIR& NOVECTOR

      do ip=1,ipsed
        tbor(ip) = 4.106e-4*sal(ip)/35.
        c1(ip) = k1(ip)/2.0
        c2(ip) = 1.0 - 4.0*k2(ip)/k1(ip)
        c4(ip) = tbor(ip)*kb(ip)
        aht(ip) = 0.74e-8
      enddo

      do icnt=1,100
        do ip=1,ipsed
          a(ip) = alk(ip) - c4(ip)/(kb(ip) + aht(ip))
          x(ip) = a(ip)/tco2(ip)
          ah1(ip) = c1(ip)/x(ip)*(1. - x(ip)
     &            + sqrt(1. + c2(ip)*x(ip)*(-2. + x(ip))))
          aht(ip)=ah1(ip)
        enddo
      enddo

      do ip=1,ipsed
        co3(ip) = (a(ip) - tco2(ip))
     &          /(1.0 - (ah1(ip)*ah1(ip))/(k1(ip)*k2(ip)))
        hco3(ip) = tco2(ip)/(1. + ah1(ip)/k1(ip) + k2(ip)/ah1(ip))
        co2(ip) = tco2(ip)/(1. + k1(ip)/ah1(ip)
     &          + k1(ip)*k2(ip)/(ah1(ip)*ah1(ip)))
      enddo

      return
      end

!=======================================================================
       subroutine o2org (rainorg, rc, kmax, zsed, delz, form, pore
     &,                  dopls, domin, dbpls, dbmin, o2, zrct, orgml
     &,                  orggg, resp_c, ipsed, ipmax, nzmax, n_control)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipmax, ipsed, kmax, loop_limit, m
      integer nzmax, n_control

!     arguments
      real rainorg(ipmax), rc(ipmax), zsed(kmax), delz(kmax)
      real form(nzmax,ipmax), pore(nzmax,ipmax), dopls(nzmax,ipmax)
      real domin(nzmax,ipmax), dbpls(nzmax,ipmax), dbmin(nzmax,ipmax)
      real o2(nzmax,ipmax),zrct(ipmax), orgml(nzmax,ipmax)
      real orggg(nzmax,ipmax)

!     results
      real resp_c(nzmax,3,ipmax)

!     local variables -- for diagnostics
      real rct_c(ipmax), rms_c(ipmax), rct_o2(ipmax), rms_o2(ipmax)

CDIR& NOVECTOR

      if (n_control .eq. 1) then
        loop_limit = 25
      else
        loop_limit = 3
      endif

      do m=1,loop_limit

        call orgc (rainorg, rc, zsed, dbpls, dbmin, pore, zrct, delz
     &,            orggg, orgml, rct_c, rms_c, kmax, ipsed, ipmax
     &,            nzmax)

        call o2ss (zsed, delz, orgml, pore, dopls, domin, rc, kmax
     &,            zrct, o2, rct_o2, rms_o2, ipsed, ipmax, nzmax)

        do ip=1,ipsed
          zrct(ip) = min(zsed(kmax), zrct(ip)*o2(1,ip)/(o2(1,ip)
     &             - o2(kmax,ip) + 1.e-20))

          if (zrct(ip) .lt. 0.1) zrct(ip) = 0.1
        enddo
      enddo

      call get_resp (rc, orgml, kmax, zrct, zsed, delz, resp_c, ipsed
     &,              ipmax, nzmax)

      return
      end

!=======================================================================
      subroutine get_resp (rc, orgml, kmax, zrct, zsed, delz, resp_c
     &,                    ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipmax, ipsed, kmax, nzmax

!     arguments
      real rc(ipmax), orgml(nzmax,ipmax), zrct(ipmax), zsed(kmax)
      real delz(kmax)

!     results
      real resp_c(nzmax,3,ipmax)

CDIR& NOVECTOR

      do i=1,kmax
        do ip=1,ipsed
          if (zsed(i) .le. zrct(ip)) then
            resp_c(i,1,ip) = rc(ip)*orgml(i,ip)
          elseif (zsed(i-1) .le. zrct(ip)) then
            resp_c(i,1,ip) = rc(ip)*orgml(i,ip)
     &                     *(zrct(ip) - zsed(i-1))/delz(i)
          elseif (zsed(i-1) .gt. zrct(ip)) then
            resp_c(i,1,ip) = 0.
          endif
          resp_c(i,2,ip) = 0.
          resp_c(i,3,ip) = 0.
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine orgc (rainorg, rc, zsed, dbpls, dbmin, pore, zrct
     &,                delz, orggg, orgml, smrct, rmserr, kmax, ipsed
     &,                ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipsed, ipmax, kmax, nzmax

!     arguments
      real rainorg(ipmax), rc(ipmax), zrct(ipmax), dbpls(nzmax,ipmax)
      real dbmin(nzmax,ipmax), pore(nzmax,ipmax), zsed(kmax), delz(kmax)

!     results
      real orggg(nzmax,ipmax), orgml(nzmax,ipmax), smrct(ipmax)
      real rmserr(ipmax)

!     local variables
      real res(kmax,ipmax), dres(kmax,3,ipmax), react(kmax,ipmax)
      real dreac(kmax,ipmax)

!     arrays for tridiag
      real a(kmax,ipmax), b(kmax,ipmax), c(kmax,ipmax), r(kmax,ipmax)
      real u(kmax,ipmax)

CDIR& NOVECTOR

      smrct(1:ipsed) = 0.

      do i=2,kmax
        do ip=1,ipsed
          if (zsed(i) .le. zrct(ip) .and. orggg(i,ip) .gt. 0) then
            react(i,ip) = - rc(ip)*orggg(i,ip)* 3.15e7
            dreac(i,ip) = - rc(ip)*3.15e7
          elseif (zsed(i-1) .le. zrct(ip)) then
            react(i,ip) = - rc(ip)*orggg(i,ip)*3.15e7
     &                  *(zrct(ip) - zsed(i-1))/delz(i)
            dreac(i,ip) = - rc(ip)*3.15e7*(zrct(ip) - zsed(i-1))/delz(i)
          elseif (zsed(i-1) .gt. zrct(ip)) then
            react(i,ip) = 0.
            dreac(i,ip) = 0.
          endif
        enddo
      enddo

      do i=3,kmax-1
        do ip=1,ipsed

!         residual(i,ip)
          res(i,ip) = dbpls(i,ip)*(orggg(i+1,ip) - orggg(i,ip))
     &              - dbmin(i,ip)*(orggg(i,ip) - orggg(i-1,ip))
     &              + react(i,ip)
!         dr/dx(i+1,ip)
          dres(i,1,ip) = dbpls(i,ip)
!         dr/dx(i,ip)
          dres(i,2,ip) = -dbpls(i,ip) - dbmin(i,ip) + dreac(i,ip)
!         dr/dx(i-1,ip)
          dres(i,3,ip) = dbmin(i,ip)
        enddo
      enddo

      do ip=1,ipsed
        i=2
!       residual(2,ip)
        res(i,ip) = dbpls(i,ip)*(orggg(i+1,ip) - orggg(i,ip))
     &            + react(i,ip) + rainorg(ip)*12./delz(i)
     &            /(1. - pore(i,ip))/2.5

!       dr/dx(i+1,ip)
        dres(i,1,ip) = dbpls(i,ip)

!       dr/dx(i,ip)
        dres(i,2,ip) = -dbpls(i,ip) + dreac(i,ip)

!       dr/dx(i-1,ip) (not defined anyway for i=2)

        i = kmax
!       residual(kmax,ip)
        res(i,ip) = -dbmin(i,ip)*(orggg(i,ip) - orggg(i-1,ip))
     .            + react(i,ip)

!       dr/dx(i+1) (not defined anyway)

!       dr/dx(i,ip)
        dres(i,2,ip) = -dbmin(i,ip) + dreac(i,ip)

!       dr/dx(i-1,ip)
        dres(i,3,ip) = dbmin(i,ip)

      enddo

!     set up the residual array
      do i=1,kmax-1
        do ip=1,ipsed
!         b(i,1,ip) = - res(i+1,ip)
!         wanted res(2) into b(1), res(3) into b(2) etc.
          r(i,ip) = -res(i+1,ip)
!         want res(2) into r(1), etc
        enddo
      enddo

!     lower off-diagonal, dri/dxi-1
      do i=1,kmax-2
        do ip=1,ipsed
!         a(i+1,i,ip) = dres(i+2,3,ip)
!         wanted dres(3,3) in a(2,1), dres(4,3) in a(3,2) etc
          a(i+1,ip) = dres(i+2,3,ip)
!         now want dres(3,3) in a(2), dres(4,3) in a(3), etc
        enddo
      enddo

!     diagonal, dri/dxi
      do i=1,kmax-1
        do ip=1,ipsed
!         a(i,i,ip) = dres(i+1,2,ip)
!         wanted dres(2,2) in a(1,1), dres(3,2) in a(2,2) etc
          b(i,ip) = dres(i+1,2,ip)
!         now want dres(2,2) in b(1), dres(3,2) in b(2), etc
        enddo
      enddo

!     upper off-diagonal, dri/dxi+1
      do i=1,kmax-2
        do ip=1,ipsed
!         a(i,i+1,ip) = dres(i+1,1,ip)
!         want dres(2,1) in a(1,2), dres(3,1) in a(2,3) etc
          c(i,ip) = dres(i+1,1,ip)
!         now want dres(2,1) in c(1), dres(3,1) in c(2) etc
        enddo
      enddo

      call tridiag (a, b, c, r, u, kmax-1, ipsed, ipmax, kmax)

!     update the concentration array
      do ip=1,ipsed
         rmserr(ip) = 0
      enddo
      do i=1,kmax-1
        do ip=1,ipsed
          orggg(i+1,ip) = orggg(i+1,ip) + u(i,ip)
          rmserr(ip) = rmserr(ip) + res(i+1,ip)**2
        enddo
      enddo

      call sldcon (orgml, orggg, 12., pore, kmax, ipsed, ipmax, nzmax)

      do i=2,kmax
        do ip=1,ipsed

          if (zsed(i) .le. zrct(ip)) then
            react(i,ip) = -rc(ip)*orggg(i,ip)*3.15e7
          elseif (zsed(i-1) .le. zrct(ip)) then
            react(i,ip) = -rc(ip)*orggg(i,ip)*3.15e7/pore(i,ip)
     &                  *(zrct(ip) - zsed(i-1))/delz(i)
          elseif (zsed(i-1) .gt. zrct(ip)) then
            react(i,ip) = 0.
          endif
          smrct(ip) = smrct(ip) + react(i,ip)*2.5*(1. - pore(i,ip))
     &              *delz(i)
        enddo
      enddo

!     into units of moles / cm2 yr
      do ip=1,ipsed
        smrct(ip) = -smrct(ip)/12.
      enddo

      do i=3,kmax-1
        do ip=1,ipsed

!         residual(i)
          res(i,ip) = dbpls(i,ip)*(orggg(i+1,ip) - orggg(i,ip))
     &              - dbmin(i,ip)*(orggg(i,ip) - orggg(i-1,ip))
     &              + react(i,ip)

        enddo
      enddo

      do ip=1,ipsed
        i = 2
!       residual(2)
        res(i,ip) = dbpls(i,ip)*(orggg(i+1,ip) - orggg(i,ip))
     &            + react(i,ip) + rainorg(ip)*12./delz(i)
     &            /(1. - pore(i,ip))/2.5
        i = kmax
!       residual(kmax,ip)
        res(i,ip) = -dbmin(i,ip)*(orggg(i,ip) - orggg(i-1,ip))
     &            + react(i,ip)
      enddo

      do ip=1,ipsed
        rmserr(ip) = 0.
      enddo
      do i=2,kmax
        do ip=1,ipsed
          rmserr(ip) = rmserr(ip) + res(i,ip)**2
          if (orggg(i,ip) .gt. 1) orggg(i,ip) = 1.
          if (orggg(i,ip) .lt. 0) orggg(i,ip) = 0.
        enddo
      enddo
      do ip=1,ipsed
        rmserr(ip) = rmserr(ip)**(0.5)
      enddo

      return
      end

!=======================================================================
      subroutine o2ss (zsed, delz, orgml, pore, dopls, domin, rate
     &,                kmax, zrct, o2, smrct, rmserr, ipsed, ipmax
     &,                nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, j, kmax, nzmax

!     arguments
      real zsed(kmax), delz(kmax), pore(nzmax,ipmax)
      real orgml(nzmax,ipmax), dopls(nzmax,ipmax)
      real domin(nzmax,ipmax), rate(ipmax), zrct(ipmax)

!     results
      real o2(nzmax,ipmax), smrct(ipmax),rmserr(ipmax)

!    internal arrays
      real dfplus(kmax,ipmax), dfzero(kmax,ipmax), dfmins(kmax,ipmax)
      real res(kmax,ipmax)

!    arrays for tridiag
      real a(kmax,ipmax), b(kmax,ipmax), c(kmax,ipmax), r(kmax,ipmax)
      real u(kmax,ipmax)

CDIR& NOVECTOR

      do j=2,kmax-1
        do ip=1,ipsed
          res(j,ip) = (dopls(j,ip)*(o2(j+1,ip) - o2(j,ip))
     &              - domin(j,ip)*(o2(j,ip) - o2(j-1,ip)))
!         units of m / cm2 (total) s
          if (zsed(j) .le. zrct(ip)) then
            res(j,ip) = res(j,ip) - 1.3*rate(ip)*orgml(j,ip)/pore(j,ip)
          elseif (zsed(j-1) .le. zrct(ip)) then
            res(j,ip) = res(j,ip) - 1.3*rate(ip)*orgml(j,ip)/pore(j,ip)
     &                *(zrct(ip) - zsed(j-1))/delz(j)
          endif
!         units of m / cm2 (total) s
          dfplus(j,ip) = dopls(j,ip)
          dfzero(j,ip) = -dopls(j,ip) - domin(j,ip)
          dfmins(j,ip) = domin(j,ip)
!         appropriate to units of m / cm2 s
        enddo
      enddo

      do ip=1,ipsed
        j=kmax
        res(j,ip) = (-domin(j,ip)*(o2(j,ip) - o2(j-1,ip)))
        if (zsed(j) .lt. zrct(ip)) then
          res(j,ip) = res(j,ip) - 1.3*rate(ip)*orgml(j,ip)/pore(j,ip)

        elseif (zsed(j-1) .lt. zrct(ip)) then
          res(j,ip) = res(j,ip) - 1.3*rate(ip)*orgml(j,ip)
     &              *(zrct(ip) - zsed(j-1))/delz(j)
         endif
         dfplus(j,ip) = 0.
         dfzero(j,ip) = -domin(j,ip)
         dfmins(j,ip) = domin(j,ip)
      enddo

      do j=1,kmax-2
        do ip=1,ipsed
          a(j+1,ip) = dfmins(j+2,ip)
!         want dr3/dx2 in a2, dr4/dx3 in a3, etc
          b(j,ip) = dfzero(j+1,ip)
!         want dr2/dx2 in b1, dr3dx3 in b2, etc
          c(j,ip) = dfplus(j+1,ip)
!         want dr2/dx3 in c1, dr3/dx4 in c2, etc
        enddo
      enddo

      do ip=1,ipsed
        b(kmax-1,ip) = dfzero(kmax,ip) + dfplus(kmax,ip)
      enddo

      do j=1,kmax-1
        do ip=1,ipsed
          r(j,ip) = -res(j+1,ip)
!         want res(2) into r(1), etc
        enddo
      enddo

      call tridiag (a, b, c, r, u, kmax-1, ipsed, ipmax, kmax)

      do ip=1,ipsed
        smrct(ip) = 0.
      enddo
      do j=2,kmax
        do ip=1,ipsed
          o2(j,ip) = o2(j,ip) + u(j-1,ip)
          if (zsed(j) .le. zrct(ip)) then
!           convert from [moles O2/l total second] to [moles O2 / cm2 yr]
            smrct(ip) = smrct(ip) - 1.3*rate(ip)*orgml(j,ip)
     &                *3.15e7/1000.*delz(j)
          elseif (zsed(j-1) .le. zrct(ip)) then
            smrct(ip) = smrct(ip) - 1.3*rate(ip)*orgml(j,ip)
     &                *(zrct(ip)-zsed(j-1))/delz(j)*3.15e7/1000.*delz(j)
           endif
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine calc_do2 (difo2, form, pore, delz, kmax, dopls, domin
     &,                    ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipsed, ipmax, kmax, nzmax

      real difo2, form(nzmax,ipmax), pore(nzmax,ipmax)
      real dopls(nzmax,ipmax), domin(nzmax,ipmax), delz(kmax)

CDIR& NOVECTOR

      difo2 = 12.e-6

      do i=3,kmax-1
        do ip=1,ipsed
          dopls(i,ip) = difo2*((form(i+1,ip) + form(i,ip))*0.5)
     &                *1./pore(i,ip)*2./((delz(i+1) + delz(i))*delz(i))
          domin(i,ip) = difo2*((form(i-1,ip) + form(i,ip))*0.5)
     &                *1./pore(i,ip)*2./((delz(i-1) + delz(i))*delz(i))
        enddo
      enddo

      do ip=1,ipsed
        i=kmax
        dopls(i,ip) = 0.
        domin(i,ip) = difo2*((form(i-1,ip) + form(i,ip))*0.5)
     &              *1./pore(i,ip)*2./((delz(i-1) + delz(i))*delz(i))
        i=2
        dopls(i,ip) = difo2*((form(i+1,ip) + form(i,ip))*0.5)
     &              *1./pore(i,ip)*2./((delz(i+1) + delz(i))*delz(i))
        domin(i,ip) = difo2*(form(i,ip) + 1)*0.5*1./pore(i,ip)
     &              *1./delz(i)**2
      enddo

      return
      end

!=======================================================================
      subroutine calss (raincal, rainorg, rainclay, resp_c, cal_c
     &,                 dissc, dissn, csat, u1, u2, dcpls, dcmin, zsed
     &,                 delz, form, pore, kmax, calml, calgg, carb
     &,                 ttrorg, ttrcal, ttral, ttrtc, difal, diftc
     &,                 c_advect, ipsed, ipmax, nzmax, n_control
     &,                 n_debug)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i_loop, i_p, ip, ipsed, ipmax, iter, k, kmax, l
      integer loop_limit, nzmax, n_control, n_debug

!     arguments
      real raincal(ipmax),rainorg(ipmax), rainclay(ipmax)
      real resp_c(nzmax,3,ipmax), cal_c(nzmax,ipmax)
      real csat(ipmax), u1(ipmax), u2(ipmax), dcpls(nzmax,3,ipmax)
      real dcmin(nzmax,3,ipmax), zsed(kmax),delz(kmax)
      real form(nzmax,ipmax), pore(nzmax,ipmax), raincal_cutoff
      real dissc, dissn

!     results
      real calml(nzmax,ipmax), calgg(nzmax,ipmax), carb(nzmax,3,ipmax)
      real ttrorg(ipmax), ttrcal(ipmax), ttral(ipmax), ttrtc(ipmax)
      real difal(ipmax), diftc(ipmax), c_advect(ipmax)

CDIR& NOVECTOR

!     internal variables
      integer still_wrong(ipmax), sum_wrong
      logical done

      loop_limit = 100

      raincal_cutoff = 0.1e-6

      sum_wrong = 0
      do ip=1,ipsed
        if (raincal(ip) .gt. raincal_cutoff) then
          sum_wrong = sum_wrong + 1
          still_wrong(sum_wrong) = ip
        else

!!         uncomment for debugging
!          if (n_debug .ge. 2) write(6,*) "cutting off site ", ip

          ttrcal(ip) = raincal(ip)
          do k=2,kmax
            calgg(k,ip) = 0.
          enddo
        endif
      enddo

!!     uncomment for debugging
!      if (n_debug .ge. 2) write(6,*) "end"

      l = 1
      done = .false.
      do while (l .le. 2000 .and. .not. done)
        l = l + 1

        call co3ss (resp_c, cal_c, dissc, dissn, csat, u1, u2, zsed
     &,             delz, form, pore, kmax, dcpls, dcmin, calml, calgg
     &,             carb, ttrorg, ttrcal, ttral, ttrtc, difal, diftc
     &,             ipsed, ipmax, nzmax, loop_limit, n_debug
     &,             still_wrong, sum_wrong, ipmax, iter)

        do ip=1,ipsed
!         [g / cm2 yr]
          c_advect(ip) = (raincal(ip) + rainclay(ip)*0.01 - ttrcal(ip))
     &                 *calgg(kmax,ip)
!          moles calcite / cm2 yr buried
        enddo

        sum_wrong = 0
        do ip=1,ipsed
          if ((calgg(kmax,ip) .gt. .001) .and. (raincal(ip) .gt.
     &      raincal_cutoff)) then
            if (abs(1. - abs((ttrcal(ip) + c_advect(ip))/raincal(ip)))
     &        .gt. 0.001) then
!             then the system is out of whack
!             outputs = inputs to 0.1%
!             0.05 for quick and dirty runs
              sum_wrong = sum_wrong + 1
              still_wrong(sum_wrong) = ip
            endif
          endif
        enddo

        if (sum_wrong .ne. 0) then

!!         uncomment for debugging, extensive diagnostic output
!          if (n_debug .ge. 1) then
!            write(6,*) 'still', sum_wrong,'  wrong in calss'
!            if (sum_wrong .lt. 10) then
!              do ip=1,sum_wrong
!                write(6,*) still_wrong(ip), ttrcal(still_wrong(ip))
!     $,           c_advect(still_wrong(ip)), raincal(still_wrong(ip))
!              enddo
!            endif
!          endif

          do k=2,kmax
            do i_loop=1,sum_wrong
              i_p = still_wrong(i_loop)
              if (ttrcal(i_p) .lt. 1e-3) then
                 calgg(k,i_p) = calgg(kmax,i_p)*raincal(i_p)
     &                        /(ttrcal(i_p) + c_advect(i_p))
              endif
            enddo
          enddo

!         sometimes get huge initial reaction rates which throw the
!         solution into a wierd mode.
!         calcite = calcite * f rain / (reac+advec)
!         c_advect = (gtotin-greact)*cal(z)

          call sldcon (calml, calgg, 100., pore, kmax, ipsed, ipmax
     &,                nzmax)

        else
          done = .true.
        endif

!       end of do while loop
      enddo

!!     uncomment for debugging, extensive diagnostic output
!      write(6,*) 'Finishing calss after iterations ', l, iter

      do ip=1,ipsed
        if (ttrcal(ip) .gt. raincal(ip)) ttrcal(ip) = raincal(ip)
        c_advect(ip) = raincal(ip) - ttrcal(ip)
      enddo

      return
      end

!=======================================================================
       subroutine co3ss (resp_c, cal_c, dissc, dissn, csat, u1, u2
     &,                  zsed, delz, form, pore, kmax, dcpls, dcmin
     &,                  calml, calgg, carb, ttrorg, ttrcal, ttral
     &,                  ttrtc, difal, diftc, ipsed, ipmax, nzmax
     &,                  loop_limit, n_debug, still_wrong, sum_wrong
     &,                  ipcmax, l)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, i_l, i_loop, ip, ipcmax, ipsed, ipmax
      integer j, k, kmax, l, loop_limit, nzmax, n_debug

!     arguments
      real resp_c(nzmax,3,ipmax), csat(ipmax),u1(ipmax), u2(ipmax)
      real zsed(kmax),delz(kmax), form(nzmax,ipmax)
      real pore(nzmax,ipmax), dcpls(nzmax,3,ipmax)
      real dcmin(nzmax,3,ipmax), dissc, dissn

      integer still_wrong(ipcmax), sum_wrong

!     results
      real cal_c(nzmax,ipmax), calml(nzmax,ipmax), calgg(nzmax,ipmax)
      real carb(nzmax,3,ipmax), ttrorg(ipmax), ttrcal(ipmax)
      real ttral(ipmax),ttrtc(ipmax), difal(ipmax), diftc(ipmax)

      real ttral_o(ipmax), ttrtc_o(ipmax), difal_o(ipmax)
      real diftc_o(ipmax), ttrcal_o(ipmax), ttrorg_o(ipmax)
      real alkbal, tcbal

!     internals
      integer still_wrong_new(ipmax), sum_wrong_new
      logical done

CDIR& NOVECTOR

      l = 1
      done = .false.
      do while (l .le. loop_limit .and. .not. done)
        l = l + 1

!!       uncomment for debugging
!        if (n_debug .ge. 2) write(6,*) 'going into co3'

        call co3 (resp_c, dissc, dissn, csat, u1, u2, dcpls, dcmin
     &,           form, pore, kmax, calml, calgg, carb, cal_c, ipsed
     &,           ipmax, n_debug, nzmax, still_wrong, sum_wrong
     &,           ipcmax)

        call sed_diag (ttral, difal, ttrtc, diftc, ttrcal, ttrorg
     &,                resp_c, cal_c, carb, pore, delz, dcmin, kmax
     &,                ipsed, ipmax, nzmax)

!!       uncomment for debugging
!        if (n_debug .ge. 2) write(6,'(a12,i4)') 'co3ss pass ',l
!        if (n_debug .gt. 3) then
!          write(6,'(a20,2g15.5)') "profile for location", ttral(3)
!     &,     calgg(kmax,3)
!          do k=1,kmax
!            write(6,'(a8,6g15.5)') "profile ", carb(k,1,n_debug)
!     &,       carb(k,2,n_debug), carb(k,3,n_debug), csat(n_debug)
!     &,       resp_c(k,1,n_debug), cal_c(k,n_debug)
!          enddo
!        endif

        sum_wrong_new = 0
        do i_loop=1,sum_wrong
          ip = still_wrong(i_loop)

          tcbal = abs(1. - abs(ttrtc(ip)/(diftc(ip) + 1.e-20)))
          alkbal = abs(1. - abs(ttral(ip)/(difal(ip) + 1.e-20)))

          if (((alkbal .gt. 0.01) .and. (ttral(ip) .gt. 1.e-12))
     &      .or. (tcbal .gt. 0.01)) then
!           then it's out of whack

!!           uncomment for debugging
!            if (n_debug .ge. 3) then
!              write(6,'(a6,3i5,3g15.5)') "rxn", i_loop, ip, l
!     &,         ttral(ip)*1.e6, difal(ip)*1.e6, alkbal
!            endif

            sum_wrong_new = sum_wrong_new + 1
            still_wrong_new(sum_wrong_new) = ip
          endif

!!         uncomment for debugging
!          if (n_debug .ge. 1) then
!            if (ttral(ip) .lt. 0.)
!     &        write(6,*) 'ttral < 0. in pass,ip ',l,ip
!            if (carb(kmax,3,ip) .lt. 0.)
!     &        write(6,*) 'co3-- < 0. in pass,ip ',l,ip
!            if (carb(kmax,3,ip) .gt. 500.e-6)
!     &        write(6,*) 'co3-- too high in pass,ip ',l,ip
!          endif

          if (ttral(ip) .lt. 0.) then
            call reset_pw (carb, csat, kmax, ip, ipmax, nzmax)
          elseif (carb(kmax,3,ip) .lt. 0.) then
            call reset_pw (carb, csat, kmax, ip, ipmax, nzmax)
          elseif (carb(kmax,3,ip) .gt. 500.e-6) then
            call reset_pw (carb, csat, kmax, ip, ipmax, nzmax)
          endif

        enddo

!!       uncomment for debugging
!        if (n_debug .ge. 3) write(6,*) "end"
!        if (n_debug .ge. 2) write(7,*) 'done with loop 30 in co3ss'

        do i_loop=1,sum_wrong_new
          still_wrong(i_loop) = still_wrong_new(i_loop)
        enddo
        sum_wrong = sum_wrong_new

!!       uncomment for debugging
!        if (n_debug .ge. 2) then
!          write(7,*) 'old loop 40 in co3ss'
!          write(6,*) 'still ', sum_wrong, ' out of whack in co3ss'
!        endif

        if (sum_wrong .eq. 0) done = .true.

!     end of do while
      enddo

      if (sum_wrong .gt. 0) then
        do i_l=1,sum_wrong
          ttrcal(still_wrong(i_l)) = 0.

!!         uncomment for debugging
!          if (n_debug .ge. 2)
!     &      write(6,*) 'zeroing fluxes in site ', still_wrong(i_l)

        enddo
      endif

!!     uncomment for debugging, extensive diagnostic output
!      if (n_debug .ge. 2) then
!     &  write (6,'(a20, i4)') 'iterations in co3ss ',l !, still_wrong(1)

      return
      end

!=======================================================================
      subroutine calc_dc (difc, form, pore, delz, kmax, dcpls, dcmin
     &,                   ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipsed, ipmax, j, kmax, nzmax

      real difc(3), form(nzmax,ipmax), pore(nzmax,ipmax), delz(kmax)
      real dcpls(nzmax,3,ipmax), dcmin(nzmax,3,ipmax)

CDIR& NOVECTOR

      difc(1) = 10.5e-6
      difc(2) = 6.4e-6
      difc(3) = 5.2e-6

      do i=3,kmax-1
        do j=1,3
          do ip=1,ipsed
            dcpls(i,j,ip) = difc(j)*(delz(i)*form(i+1,ip) + delz(i+1)
     &                    *form(i,ip))/(delz(i)+delz(i+1))*1./pore(i,ip)
     &                    *(2./((delz(i+1) + delz(i))*delz(i)))
            dcmin(i,j,ip) = difc(j)*(delz(i)*form(i-1,ip) + delz(i-1)
     &                    *form(i,ip))/(delz(i)+delz(i-1))*1./pore(i,ip)
     &                    *(2./((delz(i-1) + delz(i))*delz(i)))
          enddo
        enddo
      enddo

      do j=1,3
        do ip=1,ipsed
          i = kmax
          dcpls(i,j,ip) = 0.
          dcmin(i,j,ip) = difc(j)*(delz(i)*form(i-1,ip) + delz(i-1)
     &                  *form(i,ip))/(delz(i) + delz(i-1))*1./pore(i,ip)
     &                  *(2./((delz(i-1) + delz(i))*delz(i)))
          i = 2
          dcpls(i,j,ip) = difc(j)*(delz(i)*form(i+1,ip) + delz(i+1)
     &                  *form(i,ip))/(delz(i) + delz(i+1))*1./pore(i,ip)
     &                  *(2./((delz(i+1) + delz(i))*delz(i)))
          dcmin(i,j,ip) = difc(j)*(form(i,ip)+1)*0.5*1./pore(i,ip)
     &                  *(1./(delz(i)**2))
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine sed_diag (ttral, difal, ttrtc, diftc, ttrcal, ttrorg
     &,                    resp_c, cal_c, carb, pore, delz, dcmin, kmax
     &,                    ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!     file 'diag.for', which calculates the diffusive fluxes of
!     o2, total co2, and alkalinity at the sediment-water
!     interface, and also the integrated reaction rates of
!     those quantities.  used by co3main to determine when to
!     stop repeating the co3 subroutine.
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, j, k, kmax, nzmax

!     arguments
      real resp_c(nzmax,3,ipmax), cal_c(nzmax,ipmax)
      real carb(nzmax,3,ipmax), pore(nzmax,ipmax), delz(kmax)
      real dcmin(nzmax,3,ipmax)

!     results
      real ttral(ipmax), difal(ipmax), ttrtc(ipmax), diftc(ipmax)
      real ttrcal(ipmax),ttrorg(ipmax)

!     internal arrays
      real ttreac(3,ipmax), difflx(3,ipmax)

CDIR& NOVECTOR

!     zero the diagnostics variables, ttreac and flux
      do j=1,3
        do ip=1,ipsed
          ttreac(j,ip) = 0.
          difflx(j,ip) = 0.
          ttrcal(ip) = 0.
          ttrorg(ip) = 0.
        enddo
      enddo

!     reaction rates are in units of mol species/cm2 (total) y
      do k=1,kmax
        do ip=1,ipsed
!         replaced 3.15e7/1e3 with 3.15e4
          ttreac(1,ip) = ttreac(1,ip) + resp_c(k,1,ip)*delz(k)*3.15e4
          ttreac(2,ip) = ttreac(2,ip) + resp_c(k,2,ip)*delz(k)*3.15e4
          ttreac(3,ip) = ttreac(3,ip) + (resp_c(k,3,ip) + cal_c(k,ip))
     &                 *delz(k)*3.15e4
          ttrcal(ip) = ttrcal(ip) + cal_c(k,ip)*delz(k)*3.15e4
          ttrorg(ip) = ttrorg(ip) + (resp_c(k,1,ip) + resp_c(k,2,ip)
     &               + resp_c(k,3,ip))*delz(k)*3.15e4
        enddo
      enddo

!     the diffusive fluxes

      do ip=1,ipsed
!       replaced 3.15e7/1e3 with 3.15e4
        difflx(1,ip) = dcmin(2,1,ip)*(carb(1,1,ip) - carb(2,1,ip))
     &               *pore(2,ip)*delz(2)*3.15e4
        difflx(2,ip) = dcmin(2,2,ip)*(carb(1,2,ip) - carb(2,2,ip))
     &               *pore(2,ip)*delz(2)*3.15e4
        difflx(3,ip) = dcmin(2,3,ip)*(carb(1,3,ip) - carb(2,3,ip))
     &               *pore(2,ip)*delz(2)*3.15e4

        ttrtc(ip) = ttreac(1,ip) + ttreac(3,ip)
        ttral(ip) = ttreac(3,ip)*2.
        diftc(ip) = difflx(1,ip) + difflx(2,ip) + difflx(3,ip)
        difal(ip) = difflx(2,ip) + difflx(3,ip)*2.
      enddo

      return
      end

!=======================================================================
      subroutine tridiag (a, b, c, r, u, n, ipsed, i_size, n_size)

!-----------------------------------------------------------------------
!     from numerical recipes, page 40, adapted for parallel
!-----------------------------------------------------------------------

      implicit none

      integer ip, i_size, ipsed, j, k, n, nmax, n_size

!     arguments
      real a(n_size,i_size), b(n_size,i_size), c(n_size,i_size)
      real r(n_size,i_size), u(n_size,i_size)

!     local variables
      real gam(n,ipsed), bet(ipsed)

CDIR& NOVECTOR

      do ip=1,ipsed
        bet(ip) = b(1,ip)
        u(1,ip) = r(1,ip)/bet(ip)
      enddo
      do j=2,n
        do ip=1,ipsed
          gam(j,ip) = c(j-1,ip)/bet(ip)
          bet(ip) = b(j,ip) - a(j,ip)*gam(j,ip)
          u(j,ip) =(r(j,ip) - a(j,ip)*u(j-1,ip))/bet(ip)
        enddo
      enddo
      do j=n-1,1,-1
        do ip=1,ipsed
          u(j,ip) = u(j,ip) - gam(j+1,ip)*u(j+1,ip)
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine calc_db (db, pore, zsed, delz, kmax, dbpls, dbmin,
     &                    ipsed, ipmax, nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, k, kmax, nzmax

!     external input variables
      real pore(nzmax,ipmax), zsed(kmax+1), delz(kmax)

!     results; external variables
      real dbpls(nzmax,ipmax), dbmin(nzmax,ipmax)

      real db

CDIR& NOVECTOR

      db = 0.15

      zsed(kmax+1) = zsed(kmax) + 1.

      do k=3,kmax-1
        do ip=1,ipsed
          dbpls(k,ip) = db*2./((delz(k) + delz(k+1))*delz(k))
     &                *(1.-pore(k,ip)+1.-pore(k+1,ip))/(1.-pore(k,ip))
          dbmin(k,ip) = db*2./((delz(k) + delz(k-1))*delz(k) )
     &                *(1.-pore(k,ip)+1.-pore(k-1,ip))/(1.-pore(k,ip))
        enddo
      enddo

      do ip=1,ipsed
        k = 2
        dbpls(k,ip) = db*2./((delz(k) + delz(k+1))*delz(k))
     &              *(2. - pore(k,ip) - pore(k+1,ip))/(1. - pore(k,ip))
            dbmin(k,ip) = 0.
        k = kmax
        dbpls(k,ip) = 0.
        dbmin(k,ip) = db*2./((delz(k) + delz(k-1))*delz(k))
     &              *(2. - pore(k,ip) - pore(k-1,ip))/(1. - pore(k,ip))
      enddo

      return
      end

!=======================================================================
      subroutine pore_2_form (pore, form, kmax, expb, ipsed, ipmax
     &,                       nzmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, k, kmax, nzmax

      real form(nzmax,ipmax), pore(nzmax,ipmax), expb

      do k=1,kmax
        do ip=1,ipsed
          form(k,ip) = pore(k,ip)**expb
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine sldcon (temp_ml, temp_gg, molwt, pore, kmax, ipsed
     &,                  ipmax, nzmax)

!-----------------------------------------------------------------------
!     update the solid concentration accounts
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, k, kmax, nzmax

      real temp_ml(nzmax,ipmax), temp_gg(nzmax,ipmax), pore(nzmax,ipmax)
      real molwt

      do k=2,kmax
        do ip=1,ipsed
          temp_ml(k,ip) = temp_gg(k,ip)*2.5*(1.-pore(k,ip))*1000./molwt
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine sldfrc (temp_ml, temp_gg, molwt, pore, kmax, ipsed
     &,                  ipmax, nzmax)

!-----------------------------------------------------------------------
!     Update the solid wt. pct. accounts
!-----------------------------------------------------------------------

      implicit none

      integer i, ip, ipsed, ipmax, kmax, nzmax

      real temp_gg(nzmax,ipmax), temp_ml(nzmax,ipmax), pore(nzmax,ipmax)
      real molwt

      do i=2,kmax
        do ip=1,ipsed
          temp_gg(i,ip) = temp_ml(i,ip)*molwt
     &                  /(2.5*(1. - pore(i,ip))*1000.)
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine co3 (resp_c, dissc, dissn, csat, disoc1, disoc2, dplus
     &,               dminus, form, pore, kmax, calml, calgg, carb
     &,               cal_c, ipsed, ipmax, n_debug, nzmax, still_wrong
     &,               sum_wrong, ipcmax)

!-----------------------------------------------------------------------
!     routine co3, which calculates a single iteration
!     of the carbonate system chemistry.  must be run
!     several times because of the non-linearity of
!     calcite dissolution kinetics.
!-----------------------------------------------------------------------

      implicit none

      integer i, ib, i_loop, ip, ipsed, ipmax, ipcmax, m, info, i_row
      integer i_col, im, k, kmax, l, lda, nzmax, n_debug

      parameter (lda = 16)

!     arguments
      real resp_c(nzmax,3,ipmax), csat(ipmax), disoc1(ipmax)
      real disoc2(ipmax), dplus(nzmax,3,ipmax), dminus(nzmax,3,ipmax)
      real form(nzmax,ipmax), pore(nzmax,ipmax), calml(nzmax,ipmax)
      real calgg(nzmax,ipmax), dissc, dissn

      integer still_wrong(ipcmax), sum_wrong

!     results
      real cal_c(nzmax,ipmax), carb(nzmax,3,ipmax)

!     local variables
      real r(kmax,3,ipmax), dr(3,3,3,kmax,ipmax)
      real rmstc(ipmax), rmsal(ipmax), rmsph(ipmax), weight(ipmax)
      real trialw

      integer weight_diag(ipmax)

!     linpack variables
      integer ipvt(3*kmax)
      real bbd(ipmax,3*kmax), abd(ipmax,lda,3*kmax)

      do i=1,3
        do i_loop=1,sum_wrong
!         for the bottom boundary condition, no flux
          carb(kmax+1,i,still_wrong(i_loop)) =
     &      carb(kmax,i,still_wrong(i_loop))
        enddo
      enddo

!     the residual terms: array (depth; tc, alk, ph)

      do k=2,kmax
        do i_loop=1,sum_wrong
          ip = still_wrong(i_loop)

!         total co2 equation
          r(k-1,1,ip) = (dplus(k,1,ip)*(carb(k+1,1,ip) - carb(k,1,ip))
     &                - dminus(k,1,ip)*(carb(k,1,ip) - carb(k-1,1,ip)))
     &                + (dplus(k,2,ip)*(carb(k+1,2,ip) - carb(k,2,ip))
     &                - dminus(k,2,ip)*(carb(k,2,ip) - carb(k-1,2,ip)))
     &                + (dplus(k,3,ip)*(carb(k+1,3,ip) - carb(k,3,ip))
     &                - dminus(k,3,ip)*(carb(k,3,ip) - carb(k-1,3,ip)))
!                     units of moles / l *porewater* sec
          r(k-1,1,ip) = r(k-1,1,ip) + resp_c(k,1,ip)/pore(k,ip)
     &                + resp_c(k,2,ip)/pore(k,ip)
     &                + resp_c(k,3,ip)/pore(k,ip)
!                     units of moles / l *porewater* sec
          if (carb(k,3,ip) .lt. csat(ip)) then
            r(k-1,1,ip) = r(k-1,1,ip) + dissc
     &                  *((1. - (carb(k,3,ip)/csat(ip)))**dissn)
     &                  *(1. - pore(k,ip))/pore(k,ip)*calgg(k,ip)
     &                  *(2.5*1000)/(100)
          endif

!       alkalinity equation
          r(k-1,2,ip) = dplus(k,3,ip)*(carb(k+1,3,ip)-carb(k,3,ip))
     &                - dminus(k,3,ip)*(carb(k,3,ip)-carb(k-1,3,ip))
     &                + 0.5*dplus(k,2,ip)*(carb(k+1,2,ip)-carb(k,2,ip))
     &                - 0.5*dminus(k,2,ip)*(carb(k,2,ip)-carb(k-1,2,ip))
          r(k-1,2,ip) = r(k-1,2,ip) + resp_c(k,3,ip)/pore(k,ip)
     &                + 0.5*resp_c(k,2,ip)/pore(k,ip)
          if (carb(k,3,ip) .lt. csat(ip)) then
            r(k-1,2,ip) = r(k-1,2,ip) + dissc
     &                  *((1. - (carb(k,3,ip)/csat(ip)))**dissn)
     &                  *(1. - pore(k,ip))/pore(k,ip)*calgg(k,ip)
     &                  *(2.5*1000)/(100)
         endif

         r(k-1,3,ip) = carb(k,1,ip)*carb(k,3,ip)/carb(k,2,ip)**2
     &               - disoc2(ip)/disoc1(ip)

        enddo
      enddo

!     the derivitive terms: array (function, variable,
!     'k+'= 3 to 'k-' = 1, and depth level k)

      do k=2,kmax-1
        do i_loop=1,sum_wrong
          ip = still_wrong(i_loop)
          dr(1,1,3,k-1,ip) = dplus(k,1,ip)
          dr(1,1,2,k-1,ip) = -dplus(k,1,ip) - dminus(k,1,ip)
          dr(1,1,1,k-1,ip) = dminus(k,1,ip)

          dr(1,2,3,k-1,ip) = dplus(k,2,ip)
          dr(1,2,2,k-1,ip) = -dplus(k,2,ip) - dminus(k,2,ip)
          dr(1,2,1,k-1,ip) = dminus(k,2,ip)

          dr(1,3,3,k-1,ip) = dplus(k,3,ip)
          dr(1,3,2,k-1,ip) = -dplus(k,3,ip) - dminus(k,3,ip)

          if (carb(k,3,ip) .lt. csat(ip))
     &      dr(1,3,2,k-1,ip) = dr(1,3,2,k-1,ip) - dissc*dissn
     &                       *(1-pore(k,ip))/pore(k,ip)*calgg(k,ip)
     &                       *(2.5*1000)/100/csat(ip)
     &                       *((1-(carb(k,3,ip)/csat(ip)))**(dissn-1.))
          dr(1,3,1,k-1,ip) = dminus(k,3,ip)

          dr(2,1,3,k-1,ip) = 0.
          dr(2,1,2,k-1,ip) = 0.
          dr(2,1,1,k-1,ip) = 0.

          dr(2,2,3,k-1,ip) = 0.5*dplus(k,2,ip)
          dr(2,2,2,k-1,ip) = -0.5*dplus(k,2,ip) - 0.5*dminus(k,2,ip)
          dr(2,2,1,k-1,ip) = 0.5*dminus(k,2,ip)

          dr(2,3,3,k-1,ip) = dplus(k,3,ip)
          dr(2,3,2,k-1,ip) = -dplus(k,3,ip) - dminus(k,3,ip)
          if (carb(k,3,ip) .lt. csat(ip)) then
            dr(2,3,2,k-1,ip) = dr(2,3,2,k-1,ip) - dissc*dissn
     &                       *(1.-pore(k,ip))/pore(k,ip)*calgg(k,ip)
     &                       *(2.5*1000)/100/csat(ip)
     &                       *((1.-(carb(k,3,ip)/csat(ip)))**(dissn-1.))
          endif
          dr(2,3,1,k-1,ip) = dminus(k,3,ip)

          dr(3,1,3,k-1,ip) = 0.
          dr(3,1,2,k-1,ip) = carb(k,3,ip)/carb(k,2,ip)**2
          dr(3,1,1,k-1,ip) = 0.

          dr(3,2,3,k-1,ip) = 0.
          dr(3,2,2,k-1,ip) = -0.5*carb(k,1,ip)*carb(k,3,ip)
     &                     /carb(k,2,ip)**3
          dr(3,2,1,k-1,ip) = 0.

          dr(3,3,3,k-1,ip) = 0.
          dr(3,3,2,k-1,ip) = carb(k,1,ip)/carb(k,2,ip)**2
          dr(3,3,1,k-1,ip) = 0.

        enddo
      enddo

!     bottom special conditions

      do l=1,3 ! function
        do i=1,3 ! variable
          do m=1,3 ! above, below
            do i_loop = 1, sum_wrong
              ip = still_wrong(i_loop)
              dr(l,i,m,kmax-1,ip) = dr(l,i,m,kmax-2,ip)
            enddo
          enddo
        enddo
      enddo

!!     uncomment for debugging
!      if (n_debug .gt. 3) then
!        do ip=1,ipsed
!          if (still_wrong(ip) .eq. 3) then
!            write(6,*) "residuals"
!            i_loop = ip
!            do i=1,3
!              do k=2,kmax
!                 write(6,*) r(k-1,1,ip), r(k-1,2,ip), r(k-1,3,ip)
!              enddo
!              write(6,*)
!            enddo
!          endif
!        enddo
!      endif

!     load the big array

      do k=1,3*kmax
        do l=1,lda
          do ip=1,ipmax
            abd(ip,l,k) = 0.
          enddo
        enddo
      enddo

      do k=2,kmax ! depth level
        do l=1,3 ! function
          do m=1,3 ! up, down
            do i=1,3 !variable
              i_row = (k-2)*3 + l ! row number
              i_col = (m+k-4)*3 + i ! column number
              im = 11
              ib = i_row - i_col + im
              do i_loop=1,sum_wrong
                ip = still_wrong(i_loop)
                if ((i_col .gt. 0.) .and. (i_col .le. (kmax-1)*3)) then
                  abd(i_loop,ib,i_col) = dr(l,i,m,k-1,ip)
                  if (k .eq. kmax) then
                    if (m .eq. 2) then
                      abd(i_loop,ib,i_col) = abd(i_loop,ib,i_col)
     &                                     + dr(l,i,m+1,k-1,ip)
                    endif
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo

!     load the residual array

      do k=2,kmax
        do l=1,3
          do i_loop=1,sum_wrong
            ip = still_wrong(i_loop)
            i_row = (k-2)*3 + l
            bbd(i_loop,i_row) = -r(k-1,l,ip)
          enddo
        enddo
      enddo

      call my_sgbfa (abd, lda, (kmax-1)*3, 5, 5, ipvt, info, sum_wrong
     &,              ipmax, 3*kmax)

      call my_sgbsl (abd, lda, (kmax-1)*3, 5, 5, ipvt, bbd, sum_wrong
     &,              ipmax)

      do i_loop=1,sum_wrong
        weight(still_wrong(i_loop)) = 1.0
        weight_diag(still_wrong(i_loop)) = 0
      enddo

      if (sum_wrong .eq. 1) then
        weight(still_wrong(1)) = 0.5
      endif

      do i_loop=1,sum_wrong
        ip = still_wrong(i_loop)
        do k=2,kmax

          i_row = (k-2)*3 + 3
!         CO3= can't go up more than 75%
          trialw = -0.75*carb(k,3,ip)/(bbd(i_loop,i_row) + 1.e-20)
          if ((trialw .gt. 0.) .and. (trialw .lt. weight(ip))) then
            weight(ip) = trialw
            weight_diag(ip) = 1
          endif

          i_row = (k-2)*3 + 1
!         CO3= can't go up more than 75%
          trialw = -0.75*carb(k,1,ip)/(bbd(i_loop,i_row) + 1.e-20)
          if ((trialw .gt. 0.) .and. (trialw .lt. weight(ip))) then
            weight(ip) = trialw
            weight_diag(ip) = 2
          endif

        enddo

!!       uncomment for debugging
!        if (n_debug .ge. 3)
!     &    write(6,*) "weight", ip, weight(ip), weight_diag(ip)

      enddo

!!     uncomment for debugging
!      if (n_debug .ge. 3) write(6,*) "end"
!      if (n_debug .gt. 3) then
!        do ip=1,ipsed
!          if (still_wrong(ip) .eq. n_debug) then
!            write(6,*) "weight criterium"
!            i_loop = ip
!            do i=1,3
!              do k=2,kmax
!                i_row = (k-2)*3 + i
!                write(6,*) carb(k,i,n_debug), carb(k,i,n_debug)
!     &            + bbd(i_loop,i_row)
!              enddo
!              write(6,*)
!            enddo
!          endif
!        enddo
!      endif

      do k=2,kmax
        do i=1,3
          do i_loop=1,sum_wrong
            ip = still_wrong(i_loop)
            i_row = (k-2)*3 + i
            carb(k,i,ip) = carb(k,i,ip) + bbd(i_loop,i_row)*weight(ip)
          enddo
        enddo
      enddo

!     write(6,*) 'adjusted co3 values'

      do k=1,kmax
        do i_loop=1,sum_wrong
          ip = still_wrong(i_loop)
          if (carb(k,3,ip) .lt. csat(ip)) then
            cal_c(k,ip) = dissc*((1. - (carb(k,3,ip)/csat(ip)))**dissn)
     &                  *(1. - pore(k,ip))*calgg(k,ip)*(2.5*1000)/(100)
          else
            cal_c(k,ip) = 0.
          endif
        enddo
      enddo

!     write(6,*) 'calculated cal_c terms'

      return
      end

!=======================================================================
      subroutine calc_cal_c (resp_c, dissc, dissn, csat, disoc1, disoc2
     &,                      dplus, dminus, form, pore, kmax, calml
     &,                      calgg, carb, cal_c, ipsed, ipmax, nzmax
     &,                      still_wrong, sum_wrong, ipcmax)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      implicit none

      integer ip, ipsed, ipmax, ipcmax, k, kmax, nzmax

!     arguments
      real resp_c(nzmax,3,ipmax), csat(ipmax), disoc1(ipmax)
      real disoc2(ipmax), dplus(nzmax,3,ipmax)
      real dminus(nzmax,3,ipmax), form(nzmax,ipmax), pore(nzmax,ipmax)
      real calml(nzmax,ipmax), calgg(nzmax,ipmax), dissc, dissn

      integer still_wrong(ipcmax), sum_wrong

!     results
      real cal_c(nzmax,ipmax), carb(nzmax,3,ipmax)

!     local variables
      real rmstc(ipmax), rmsal(ipmax), rmsph(ipmax), weight(ipmax)

      do k=1,kmax
        do ip=1,ipsed
          if (carb(k,3,ip) .lt. csat(ip)) then
            cal_c(k,ip) = dissc*((1. - (carb(k,3,ip)/csat(ip)))**dissn)
     &                  *(1. - pore(k,ip))*calgg(k,ip)*(2.5*1000)/(100)
          else
            cal_c(k,ip) = 0.
          endif
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine my_sgbfa (abd, lda, n, ml, mu, ipvt, info, n_matrices
     &,                    ipmax, ncolmax)

!-----------------------------------------------------------------------
!     modified by archer, Oct 92, for vectorized co3 solution
!     no pivoting.   operates in vector form on entire array of
!     similar matrices.

!     expects abd to be filled continuously to first index = n_matrices

!     sgbfa factors a real band matrix by elimination.

!     sgbfa is usually called by sgbco, but it can be called
!     directly with a saving in time if  rcond  is not needed.

!     on entry

!        abd     real(lda, n)
!                contains the matrix in band storage.  the columns
!                of the matrix are stored in the columns of  abd  and
!                the diagonals of the matrix are stored in rows
!                ml+1 through 2*ml+mu+1 of  abd .
!                see the comments below for details.

!        lda     integer
!                the leading dimension of the array  abd .
!                lda must be .ge. 2*ml + mu + 1 .

!        n       integer
!                the order of the original matrix.

!        ml      integer
!                number of diagonals below the main diagonal.
!                0 .le. ml .lt. n .

!        mu      integer
!                number of diagonals above the main diagonal.
!                0 .le. mu .lt. n .
!                more efficient if  ml .le. mu .
!     on return

!        abd     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.

!        ipvt    integer(n)
!                an integer vector of pivot indices.

!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgbsl will divide by zero if
!                     called.  use  rcond  in sgbco for a reliable
!                     indication of singularity.
!-----------------------------------------------------------------------

      implicit none

      integer lda, n, ml, mu, ipvt(*), info, nmat, ipmax, n_matrices
      integer ncolmax

      real abd(ipmax,lda,ncolmax), t(ipmax)

      integer i, isamax, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm
      integer nm1

      m = ml + mu + 1
      info = 0

!     zero initial fill-in columns
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .gt. j0) then
        do jz=j0,j1
          i0 = m + 1 - jz
          do i=i0,ml
            do nmat=1,n_matrices
              abd(nmat,i,jz) = 0.
            enddo
          enddo
        enddo
      endif

      jz = j1
      ju = 0

!     gaussian elimination with partial pivoting
      nm1 = n - 1
      if (nm1 .ge. 1) then

        do k=1,nm1
          kp1 = k + 1

!         zero next fill-in column
          jz = jz + 1
          if ((jz .le. n) .and. (ml .ge. 1)) then
            do i=1,ml
              do nmat=1,n_matrices
                abd(nmat,i,jz) = 0.
              enddo
            enddo
          endif

!         hardwire l = pivot index to no pivoting
          lm = min0(ml,n-k)
          l = m
          ipvt(k) = k

!         zero pivot implies this column already triangularized
!         do test on just one element
          if (abd(1,l,k) .ne. 0.) then
!           compute multipliers
            do nmat=1,n_matrices
              if (abd(nmat,m,k) .eq. 0.) write(6,*) 'div by 0'
     &,         nmat, m, k
              t(nmat) = -1./abd(nmat,m,k)
            enddo

            call my_sscal (lm, t(1), abd(1,m+1,k), n_matrices, ipmax)

!           row elimination with column indexing
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .ge. kp1) then
              do j=kp1,ju
                l = l - 1
                mm = mm - 1
                do nmat=1,n_matrices
                  t(nmat) = abd(nmat,l,j)
                enddo
                if (l .ne. mm) then
                  do nmat=1,n_matrices
                    abd(nmat,l,j) = abd(nmat,mm,j)
                    abd(nmat,mm,j) = t(nmat)
                  enddo
                endif
                call my_saxpy (lm, t(1), abd(1,m+1,k), abd(1,mm+1,j)
     &,                        n_matrices, ipmax)
              enddo
            endif
          else
            info = k
          endif

        enddo

      endif

      ipvt(n) = n
      if (abd(1,m,n) .eq. 0.) info = n

      return
      end

!=======================================================================
      subroutine my_sgbsl (abd, lda, n, ml, mu, ipvt, b, n_matrices
     &,                    ipmax)

!-----------------------------------------------------------------------

!     sgbsl solves the real band system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgbco or sgbfa.

!     modified by archer for carbonate solver, Oct 92
!     only job = 0 (solve a*x = b) has been kept

!     on entry

!        abd     real(lda, n)
!                the output from sgbco or sgbfa.

!        lda     integer
!                the leading dimension of the array  abd .

!        n       integer
!                the order of the original matrix.

!        ml      integer
!                number of diagonals below the main diagonal.

!        mu      integer
!                number of diagonals above the main diagonal.

!        ipvt    integer(n)
!                the pivot vector from sgbco or sgbfa.

!        b       real(n)
!                the right hand side vector.

!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b , where
!                            trans(a)  is the transpose.

!     on return

!        b       the solution vector  x .
!-----------------------------------------------------------------------

      implicit none

      integer lda, n, ml, mu, ipvt(*), job, nmat, n_matrices, ipmax

      real abd(ipmax,lda,*), b(ipmax,*), t(ipmax)

      integer k, kb, l, la, lb, lm, m, nm1

      m = mu + ml + 1
      nm1 = n - 1

!     job = 0 , solve  a * x = b
!     first solve l*y = b
      if ((ml .ne. 0) .and. (nm1 .ge. 1)) then
        do k=1,nm1
          lm = min0(ml,n-k)
          l = ipvt(k)
          do nmat=1,n_matrices
            t(nmat) = b(nmat,l)
          enddo
          if (l .ne. k) then
            do nmat=1,n_matrices
               b(nmat,l) = b(nmat,k)
               b(nmat,k) = t(nmat)
            enddo
          endif
          call my_saxpy (lm, t(1), abd(1,m+1,k), b(1,k+1)
     &,                  n_matrices, ipmax)
        enddo
      endif

!     now solve  u*x = y
      do kb=1,n
        k = n + 1 - kb
        do nmat=1,n_matrices
          b(nmat,k) = b(nmat,k)/abd(nmat,m,k)
        enddo
        lm = min0(k,m) - 1
        la = m - lm
        lb = k - lm
        do nmat=1,n_matrices
          t(nmat) = -b(nmat,k)
        enddo
        call my_saxpy (lm, t(1), abd(1,la,k), b(1,lb)
     &,                n_matrices, ipmax)
      enddo

      return
      end

!=======================================================================
      subroutine my_sscal (n, sa, sx, n_elements, n_elem_max)

!-----------------------------------------------------------------------
!     scales a vector by a constant.
!-----------------------------------------------------------------------

      implicit none

      integer i, n, nelem, n_elements, n_elem_max

      real sa(*), sx(n_elem_max,*)

      if (n .le. 0)return
      do i=1,n
        do nelem=1,n_elements
          sx(nelem,i) = sa(nelem)*sx(nelem,i)
        enddo
      enddo

      return
      end

!=======================================================================
      subroutine my_saxpy (n, sa, sx, sy, n_matrices, ipmax)

!-----------------------------------------------------------------------
!     constant times a vector plus a vector.
!-----------------------------------------------------------------------

      implicit none

      integer i, n, n_matrices, nmat, ipmax

      real sx(ipmax,*), sy(ipmax,*), sa(*)

      if (n .le. 0) return
      do i=1,n
        do nmat=1,n_matrices
          sy(nmat,i) = sy(nmat,i) + sa(nmat)*sx(nmat,i)
        enddo
      enddo
#endif

      return
      end
