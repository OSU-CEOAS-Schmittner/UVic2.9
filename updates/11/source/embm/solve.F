      subroutine solve (n)

#if defined O_embm
!=======================================================================
!     solve for tracer distribution after diffusion

!     input:
!       n    = tracer number
!=======================================================================

      implicit none

      integer i, ii, ierr, j, jj, k, n

      logical done

      real afw, afe, afn, atc, ate, atn, ats, atnc, atsc, atw, b, dt
      real dtss, dfw, dfe, dfn, fa, fb, fc, fd, ff, fg, fh, tmp, x

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "solve.h"
      include "atm.h"
      include "cembm.h"
      include "csbc.h"
      include "grdvar.h"
      include "coord.h"
      include "levind.h"
# if defined O_ice_cpts && defined O_ice
      include "cpts.h"
# endif
      include "ice.h"

      real forc(imt,jmt)
# if defined O_embm_solve2x || defined O_embm_solve2y
      real tmp_at(imt,jmt)
# endif

# if defined O_embm_explicit
      dtss = dts/ns

# else
      dtss = dts

# endif
!-----------------------------------------------------------------------
!     set the forcing for each tracer
!-----------------------------------------------------------------------

      if (n .eq. isat) then

!       temperature

        fa = dtss/(cpatm*rhoatm*sht)
        fb = dtss*vlocn/(cpatm*rhoatm*sht)
        fc = dtss*slice/(cpatm*rhoatm*sht) - fb
        fd = scatter*(1. + pass)
        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = fa*(solins(i,j)*sbc(i,j,iaca)*fd
     &                - dnswr(i,j)*scatter - outlwr(i,j)
     &                + uplwr(i,j) + upsens(i,j))
!           latent heat from total precipitation as water
            forc(i,j) = forc(i,j) + precip(i,j)*fb
!           correct for latent heat from snow
            forc(i,j) = forc(i,j) + fc*psno(i,j)
          enddo
        enddo

      elseif (n .eq. ishum) then

!       humidity

        fa = dtss/(rhoatm*shq)
        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = fa*evap(i,j)
          enddo
        enddo

# if defined O_carbon_co2_2d && !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
      elseif (n .eq. ico2) then

        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = -dtss*flux(i,j,ico2)*gtoppm
          enddo
        enddo

# endif
      else

!       other tracers

        do j=2,jmtm1
          do i=2,imtm1
            forc(i,j) = c0
          enddo
        enddo

      endif
      call embmbc (forc)

!-----------------------------------------------------------------------
!     calculate new coefficients if required
!-----------------------------------------------------------------------

      if (newcoef(lf,n)) call coef (n)
# if defined O_embm_explicit
      newcoef(1,n) = .false.
      newcoef(2,n) = .false.
# endif

!-----------------------------------------------------------------------
!     shuffle in time
!-----------------------------------------------------------------------

# if defined O_embm_explicit
      do j=1,jmt
        do i=1,imt
          tmp = at(i,j,2,n)
          at(i,j,2,n) = at(i,j,lf,n)
          at(i,j,1,n) = tmp
        enddo
      enddo
      call embmbc (at(1,1,2,n))

!-----------------------------------------------------------------------
!     forward step through explicit diffusion and upstream advection
!-----------------------------------------------------------------------

      do k=1,ns
        dfs(1:imt) = 0.0
        afs(1:imt) = 0.0
        do j=2,jmtm1
          dfw = dce(1,j,n)*(at(2,j,2,n) - at(1,j,2,n))
          afw = ace(1,j,n)*(at(1,j,2,n) + at(2,j,2,n))
     &        + abs(ace(1,j,n))*(at(1,j,2,n) - at(2,j,2,n))
          do i=2,imtm1
            dfe = dce(i,j,n)*(at(i+1,j,2,n) - at(i,j,2,n))
            dfn = dcn(i,j,n)*(at(i,j+1,2,n) - at(i,j,2,n))
            afe = ace(i,j,n)*(at(i,j,2,n) + at(i+1,j,2,n))
     &          + abs(ace(i,j,n))*(at(i,j,2,n) - at(i+1,j,2,n))
            afn = acn(i,j,n)*(at(i,j,2,n) + at(i,j+1,2,n))
     &          + abs(acn(i,j,n))*(at(i,j,2,n) - at(i,j+1,2,n))
            at(i,j,2,n) = at(i,j,2,n) + forc(i,j) + dtss*(
     &                    (dfe - dfw)*dxtr(i)
     &                  + (dfn - dfs(i))*cstdytr(j)
     &                  - ((afe - afw)*dxt2r(i)
     &                  + (afn - afs(i))*dyt2r(j))*cstr(j))
            dfw = dfe
            dfs(i) = dfn
            afw = afe
            afs(i) = afn
          enddo
        enddo
        call embmbc (at(1,1,2,n))
      enddo

# else
      do j=1,jmt
        do i=1,imt
          tmp = at(i,j,2,n)
          at(i,j,2,n) = at(i,j,lf,n) + forc(i,j)
          at(i,j,1,n) = tmp
        enddo
      enddo

!-----------------------------------------------------------------------
!     load rhs into the solver array
!-----------------------------------------------------------------------

      k = 0
#  if defined O_embm_solve2y
      do jj=1,jjmtm2
        j = jj*2
#  else
      do j=2,jmtm1
#  endif
#  if defined O_embm_solve2x
        do ii=1,iimtm2
          i = ii*2
#  else
        do i=2,imtm1
#  endif
#  if defined O_embm_solve2x || defined O_embm_solve2y
          b = (at(i,j,2,n))*gr(i,j)
          x = at(i,j,1,n)*gr(i,j)
#   if defined O_embm_solve2x
          b = b + (at(i+1,j,2,n))*gr(i+1,j)
          x = x + at(i+1,j,1,n)*gr(i+1,j)
#   endif
#   if defined O_embm_solve2y
          b = b + (at(i,j+1,2,n))*gr(i,j+1)
          x = x + at(i,j+1,2,n)*gr(i,j+1)
#   endif
#   if defined O_embm_solve2x && defined O_embm_solve2y
          b = b + (at(i+1,j+1,2,n))*gr(i+1,j+1)
          x = x + at(i+1,j+1,1,n)*gr(i+1,j+1)
#   endif
          k = k + 1
          bv(k) = b
          xv(k) = x
#  else
          k = k + 1
          bv(k) = at(i,j,2,n)
          xv(k) = at(i,j,1,n)
#  endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     solve for tracer
!-----------------------------------------------------------------------

#  if defined O_embm_adi
      call adi (xv, an(1,lf,n), ans(1,lf,n), as(1,lf,n), ae(1,lf,n)
     &,         aew(1,lf,n), aw(1,lf,n), bv, iimtm2, jjmtm2)
      itout(n) = 1
      epsout(n) = epsin(n)
#  endif
#  if defined O_embm_mgrid
      call mgrid (xv, ap(1,lf,n), an(1,lf,n), as(1,lf,n), ae(1,lf,n)
     &,           aw(1,lf,n), bv, 1, iimtm2, 1, jjmtm2, iimtm2, jjmtm2
     &,           itin(n), levelin, epsin(n), itout(n), levelout
     &,           epsout(n))

#  endif
#  if defined O_embm_slap
      call slap_sslugm (nord, bv, xv, nelm, ia, ja, ar(1,lf,n), 0, 10
     &,                 0, epsin(n), itin(n), itout(n), epsout(n)
     &,                 ierr, 0, raux, nraux, iaux, niaux)
#  endif
#  if defined O_embm_essl
      rparm(1) = epsin(n)
      iparm(1) = itin(n)
      if (newcoef(lf,n)) then
        call dsris ('G', 'I', nord, ar(1,lf,n), ja, ia, bv, xv
     &,             iparm, rparm, aux1(1,lf,n), naux1, aux2, naux2)
      else
        call dsris ('G', 'S', nord, ar(1,lf,n), ja, ia, bv, xv
     &,             iparm, rparm, aux1(1,lf,n), naux1, aux2, naux2)
      endif
      epsout(n) = rparm(1)
      itout = iparm(6)
#  endif
#  if defined O_embm_sparskit
      fpar(1) = epsin(n)
      ipar(6) = itin(n)

      ipar(1) = 0
      done = .false.
      do while (.not. done)
        call bcg (nord, bv, xv, ipar, fpar, work(1,lf,n))
        if (ipar(1) .eq. 1) then
          call amux(n, work(ipar(8),lf,n), work(ipar(9),lf,n)
     &,     ar(1,lf,n), ja, ia)
          done = .false.
        elseif (ipar(1) .eq. 2) then
          call atmux(nord, work(ipar(8),lf,n), work(ipar(9),lf,n)
     &,     ar(1,lf,n), ja, ia)
          done = .false.
        elseif (ipar(1) .eq.3 .or. ipar(1) .eq. 5) then
          call lusol(nord, work(ipar(8),lf,n), work(ipar(9),lf,n)
     &,     aur(1,lf,n), jau, ju)
          done = .false.
        elseif (ipar(1) .eq. 4 .or. ipar(1) .eq. 6) then
          call lutsol(nord, work(ipar(8),lf,n), work(ipar(9),lf,n)
     &,     aur(1,lf,n), jau, ju)
          done = .false.
        elseif (ipar(1) .le. 0) then
          done = .true.
          if (ipar(1) .eq. -1) then
            print *, 'Iterative solver has iterated too many times.'
          elseif (ipar(1) .eq. -2) then
            print *, 'Iterative solver was not given enough work space.'
            print *, 'The work space should at least have ', ipar(4)
     &,       ' elements.'
          elseif (ipar(1) .eq. -3) then
            print *, 'Iterative solver is facing a break-down.'
          endif
        endif
      enddo
      epsout(n) = fpar(6)
      itout = ipar(7)
#  endif
      newcoef(lf,n) = .false.
#  if !defined O_global_sums
      if (epsout(n) .gt. epsin(n)) write(*,*)
     &  '==> Warning:  atmospheric solver not converging in ',
     &  itout(n),' iterations ( eps = ',epsout(n), ' > ',epsin(n),' )'
#  endif

!-----------------------------------------------------------------------
!     copy new solution from left hand side
!-----------------------------------------------------------------------

      k = 0
#  if defined O_embm_solve2y
      do jj=1,jjmtm2
        j = jj*2
#  else
      do j=2,jmtm1
#  endif
#  if defined O_embm_solve2x
        do ii=1,iimtm2
          i = ii*2
#  else
        do i=2,imtm1
#  endif
          k = k + 1
#  if defined O_embm_solve2x || defined O_embm_solve2y
          tmp_at(i,j) = xv(k)
#  else
          at(i,j,2,n) = xv(k)
#  endif
#  if defined O_embm_solve2x
          tmp_at(i+1,j) = xv(k)
#  endif
#  if defined O_embm_solve2y
          tmp_at(i,j+1) = xv(k)
#  endif
#  if defined O_embm_solve2x && defined O_embm_solve2y
          tmp_at(i+1,j+1) = xv(k)
#  endif
        enddo
      enddo

#  if defined O_embm_solve2y || defined O_embm_solve2x
!-----------------------------------------------------------------------
!     interpolate back to the fine atmospheric grid
!-----------------------------------------------------------------------

      call embmbc (tmp_at)

#   if defined O_embm_solve2y
      do jj=1,jjmtm2
        j = jj*2
#   else
      do j=2,jmtm1
#   endif
#   if defined O_embm_solve2x
        do ii=1,iimtm2
          i = ii*2
#   else
        do i=2,imtm1
#   endif

          atc = tmp_at(i,j)
#   if defined O_embm_solve2x
          ff = tmp_at(i+2,j) - tmp_at(i,j)
          fg = tmp_at(i,j) - tmp_at(i-1,j)
          fh = (sign(c1,ff) + sign(c1,fg))*min(abs(ff),abs(fg))
          atw = atc - fh*wti(i)
          ate = atc + fh*wti(i+1)
          dt = (ate - atc)*xgrd(i+1) + (atw - atc)*xgrd(i)
          at(i,j,2,n) = atw - dt
          at(i+1,j,2,n) = ate - dt

#    if defined O_embm_solve2y
          ff = tmp_at(i+2,j+1) - tmp_at(i,j+1)
          fg = tmp_at(i,j+1) - tmp_at(i-1,j+1)
          fh = (sign(c1,ff) + sign(c1,fg))*min(abs(ff),abs(fg))
          atw = atc - fh*wti(i)
          ate = atc + fh*wti(i+1)
          dt = (ate - atc)*xgrd(i+1) + (atw - atc)*xgrd(i)
          at(i,j+1,2,n) = atw - dt
          at(i+1,j+1,2,n) = ate - dt

#    endif
#   endif
#   if defined O_embm_solve2y
#    if defined O_embm_solve2x
          atsc = at(i,j,2,n)
          atnc = at(i,j+1,2,n)
          ff = tmp_at(i,j+2) - tmp_at(i,j)
          fg = tmp_at(i,j) - tmp_at(i,j-1)
          fh = (sign(c1,ff) + sign(c1,fg))*min(abs(ff),abs(fg))
          ats = atsc - fh*wtj(j)
          atn = atnc + fh*wtj(j+1)
          dt = (atn - atnc)*ygrd(j+1) + (ats - atsc)*ygrd(j)
          at(i,j,2,n) = ats - dt
          at(i,j+1,2,n) = atn - dt

          atsc = at(i+1,j,2,n)
          atnc = at(i+1,j+1,2,n)
          ff = tmp_at(i+1,j+2) - tmp_at(i+1,j)
          fg = tmp_at(i+1,j) - tmp_at(i+1,j-1)
          fh = (sign(c1,ff) + sign(c1,fg))*min(abs(ff),abs(fg))
          ats = atsc - fh*wtj(j)
          atn = atnc + fh*wtj(j+1)
          dt = (atn - atnc)*ygrd(j+1) + (ats - atsc)*ygrd(j)
          at(i+1,j,2,n) = ats - dt
          at(i+1,j+1,2,n) = atn - dt

#    else
          ff = tmp_at(i,j+2) - tmp_at(i,j)
          fg = tmp_at(i,j) - tmp_at(i,j-1)
          fh = (sign(c1,ff) + sign(c1,fg))*min(abs(ff),abs(fg))
          ats = atc - fh*wtj(j)
          atn = atc + fh*wtj(j+1)
          dt = (atn - atc)*ygrd(j+1) + (ats - atc)*ygrd(j)
          at(i,j,2,n) = ats - dt
          at(i,j+1,2,n) = atn - dt

#    endif
#   endif
        enddo
      enddo

#  endif
!-----------------------------------------------------------------------
!     set boundary conditions
!-----------------------------------------------------------------------

      call embmbc (at(1,1,2,n))
# endif

      return
      end

      subroutine coef (n)

!=======================================================================
!     compute matrix coefficients

!     input:
!       n    = tracer number
!=======================================================================

      implicit none

      integer i, ide, ii, ielm, iord, j, jdn, jj, jord, n, iwx, iwy

      real acej, acnj, adde, addn, adds, addw, cc, ce, cew, cmax, cn
      real cns, cs, cw, dcej, dcnj, fe, fn, fs, fw, ue, un, uw, ve
      real vn, vs, vw

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "solve.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "csbc.h"
# if defined O_embm_adiff
      real diff_factor

      diff_factor = max(0., (1. + adiff*dtbar))

# endif

      ielm = 0
      iord = 0

      iwx = 0
      iwy = 0
      if (n .eq. isat) then
        iwx = iwxt
        iwy = iwyt
      elseif (n .eq. ishum) then
        iwx = iwxq
        iwy = iwyq
      elseif (n .eq. ico2) then
        iwx = iwxc
        iwy = iwyc
      endif

      if (iwx .gt. 0) call embmbc (sbc(1,1,iwx))
      if (iwy .gt. 0) call embmbc (sbc(1,1,iwy))

      cmax = 3.9e10

# if defined O_embm_explicit
!-----------------------------------------------------------------------
!     set up coefficients for explicit diffusion
!-----------------------------------------------------------------------

      do j=2,jmtm1
        dcej = cstr(j)*cstr(j)*filter(j)
        dcnj = csu(j)*dyur(j)
        do i=2,imtm1
          dce(i,j,n) = de(i,j,n)*dcej*dxur(i)
          dcn(i,j,n) = dn(i,j,n)*dcnj
#  if defined O_embm_adiff
          if (n .eq. isat) then
            dce(i,j,n) = dce(i,j,n)*diff_factor
            dcn(i,j,n) = dcn(i,j,n)*diff_factor
          endif
#  endif
        enddo
      enddo
      dce(1,1:jmt,n) = dce(imtm1,1:jmt,n)
      dcn(1:imt,jmtm1,n) = 0.0
      dcn(1:imt,1,n) = 0.0

!-----------------------------------------------------------------------
!     set up coefficients for explicit advection
!-----------------------------------------------------------------------

      do j=2,jmtm1
        acej = dyt2r(j)*filter(j)
        acnj = csu(j)
        do i=2,imtm1
          ace(i,j,n) = (sbc(i,j-1,iwx)*dyu(j-1)
     &               + sbc(i,j,iwx)*dyu(j))*acej
          acn(i,j,n) = (sbc(i-1,j,iwy)*dxu(i-1)
     &               + sbc(i,j,iwy)*dxu(i))*dxt2r(i)*acnj
        enddo
      enddo
      ace(1,1:jmt,n) = ace(imtm1,1:jmt,n)
      acn(1:imt,jmtm1,n) = 0.0
      acn(1:imt,1,n) = 0.0

# else
#  if defined O_embm_solve2y
      do jj=1,jjmtm2
        j = jj*2
        jdn = j + 1
#  else
      do j=2,jmtm1
        jj = j - 1
        jdn = j
#  endif
#  if defined O_embm_solve2x
        do ii=1,iimtm2
          i = ii*2
          ide = i + 1
#  else
        do i=2,imtm1
          ii = i - 1
          ide = i
#  endif

!-----------------------------------------------------------------------
!         set coefficients for implicit diffusion
!-----------------------------------------------------------------------

#  if defined O_embm_solve2x
          cs = 0.5*(dn(i,j-1,n) + dn(i+1,j-1,n))
          cn = 0.5*(dn(i,jdn,n) + dn(i+1,jdn,n))
#  else
          cs = dn(i,j-1,n)
          cn = dn(i,jdn,n)
#  endif
#  if defined O_embm_solve2y
          cw = 0.5*(de(i-1,j,n) + de(i-1,j+1,n))
          ce = 0.5*(de(ide,j,n) + de(ide,j+1,n))
#  else
          cw = de(i-1,j,n)
          ce = de(ide,j,n)
#  endif
#  if defined O_embm_adiff
          if (n .eq. isat) then
            cn = cn*diff_factor
            cs = cs*diff_factor
            ce = ce*diff_factor
            cw = cw*diff_factor
          endif
#  endif

!-----------------------------------------------------------------------
!         closed north/south boundary conditions for diffusion
!-----------------------------------------------------------------------

          if (j .eq. 2) cs = c0
          if (j .eq. jmtm1) cn = c0

          cs =-dts*cs*dsgrd(j)
          cn =-dts*cn*dngrd(j)
#  if defined O_embm_solve2y
          cw =-dts*cw*csur(j)*csur(j)*dwgrd(i)
          ce =-dts*ce*csur(j)*csur(j)*degrd(i)
#  else
          cw =-dts*cw*cstr(j)*cstr(j)*dwgrd(i)
          ce =-dts*ce*cstr(j)*cstr(j)*degrd(i)
#  endif

#  if defined O_embm_adi
          cns = 1.0 - (cs + cn)
          cew = 1.0 - (ce + cw)
#  else
          cc = 1.0 - cs - cn - cw - ce
#  endif

!-----------------------------------------------------------------------
!         set coefficients for up-stream advection
!-----------------------------------------------------------------------

#  if defined O_embm_solve2y
          vs = 2.0*sbc(i,j-1,iwy)
          vn = 2.0*sbc(i,j+1,iwy)
#  else
          vs = (sbc(i-1,j-1,iwy) + sbc(i,j-1,iwy))
          vn = (sbc(i-1,j,iwy) + sbc(i,j,iwy))
#  endif
#  if defined O_embm_solve2x
          uw = 2.0*sbc(i-1,j,iwx)
          ue = 2.0*sbc(i+1,j,iwx)
#  else
          uw = (sbc(i-1,j-1,iwx) + sbc(i-1,j,iwx))
          ue = (sbc(i,j-1,iwx) + sbc(i,j,iwx))
#  endif

          fs = p5*(c1 + sign(c1,vs))
          fn = p5*(c1 + sign(c1,vn))
          fw = p5*(c1 + sign(c1,uw))
          fe = p5*(c1 + sign(c1,ue))

!-----------------------------------------------------------------------
!         closed north/south boundary conditions for advection
!-----------------------------------------------------------------------

          if (j .eq. 2) vs = c0
          if (j .eq. jmtm1) vn = c0

          cs = cs - dts*fs*vs*asgrd(j)
          cn = cn + dts*(c1-fn)*vn*angrd(j)
#  if defined O_embm_solve2x
          cw = cw - dts*fw*uw*csur(j)*azgrd(i)
          ce = ce + dts*(c1-fe)*ue*csur(j)*azgrd(i)
#  else
          cw = cw - dts*fw*uw*cstr(j)*azgrd(i)
          ce = ce + dts*(c1-fe)*ue*cstr(j)*azgrd(i)
#  endif
#  if defined O_embm_adi
          cns = cns + dts*(fn*vn*angrd(j)-(c1-fs)*vs*asgrd(j))
#   if defined O_embm_solve2x
          cew = cew + dts*(fe*ue -(c1-fw)*uw)*csur(j)*azgrd(i)
#   else
          cew = cew + dts*(fe*ue -(c1-fw)*uw)*cstr(j)*azgrd(i)
#   endif
#  else
          cc = cc + dts*(fn*vn*angrd(j)-(c1-fs)*vs*asgrd(j)
#  if defined O_embm_solve2x
     &       + (fe*ue - (c1-fw)*uw)*csur(j)*azgrd(i))
#   else
     &       + (fe*ue - (c1-fw)*uw)*cstr(j)*azgrd(i))
#   endif
#  endif

          iord = iord + 1
#  if defined O_embm_adi
#   if defined O_embm_solve2x
          jord = jj + (ii-1)*iimtm2
#   else
          jord = jj + (ii-1)*(imtm1-1)
#   endif

!-----------------------------------------------------------------------
!         load the coefficients for the ADI solver
!-----------------------------------------------------------------------

          ans(jord,lf,n) = cns
          an(jord,lf,n)  = -cn
          as(jord,lf,n)  = -cs

          aew(iord,lf,n) = cew
          ae(iord,lf,n)  = -ce
          aw(iord,lf,n)  = -cw

        enddo
      enddo
#  endif
#  if defined O_embm_mgrid

!-----------------------------------------------------------------------
!         load the coefficients for the multigrid solver
!-----------------------------------------------------------------------

          ap(iord,lf,n) = cc
          an(iord,lf,n) = -cn
          as(iord,lf,n) = -cs
          ae(iord,lf,n) = -ce
          aw(iord,lf,n) = -cw

        enddo
      enddo

#  endif
#  if defined O_embm_slap

!-----------------------------------------------------------------------
!         load the coefficients for the slap solver
!-----------------------------------------------------------------------

!         central coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = cc
          ia(ielm) = iord
          ja(ielm) = iord

!         western coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = cw
          ia(ielm) = iord
          if (ii .gt. 1) then
            ja(ielm) = iord - 1
          else
            ja(ielm) = iord + (iimtm2-1)
          endif

!         eastern coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = ce
          ia(ielm) = iord
          if (ii .lt. iimtm2) then
            ja(ielm) = iord + 1
          else
            ja(ielm) = iord - (iimtm2-1)
          endif

!         southern coefficient
          if (jj .gt. 1) then
            ielm = ielm + 1
            ar(ielm,lf,n) = cs
            ia(ielm) = iord
            ja(ielm) = iord - iimtm2
          endif

!         northern coefficient
          if (jj .lt. jjmtm2) then
            ielm = ielm + 1
            ar(ielm,lf,n) = cn
            ia(ielm) = iord
            ja(ielm) = iord + iimtm2
          endif

        enddo
      enddo

#  endif
#  if defined O_embm_essl || defined O_embm_sparskit

!-----------------------------------------------------------------------
!         load the coefficients by rows
!-----------------------------------------------------------------------

          ia(iord) = ielm + 1

!         southern coefficient
          if (jj .gt. 1) then
            ielm = ielm + 1
            ar(ielm,lf,n) = cs
            ja(ielm) = iord - iimtm2
          endif

!         western coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = cw
          if (ii .gt. 1) then
            ja(ielm) = iord - 1
          else
            ja(ielm) = iord + (iimtm2-1)
          endif

!         central coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = cc
          ja(ielm) = iord

!         eastern coefficient
          ielm = ielm + 1
          ar(ielm,lf,n) = ce
          if (ii .lt. iimtm2) then
            ja(ielm) = iord + 1
          else
            ja(ielm) = iord - (iimtm2-1)
          endif

!         northern coefficient
          if (jj .lt. jjmtm2) then
            ielm = ielm + 1
            ar(ielm,lf,n) = cn
            ja(ielm) = iord + iimtm2
          endif

        enddo
      enddo

      ia(iord+1) = ielm + 1

#  endif
# endif
#endif

      return
      end
