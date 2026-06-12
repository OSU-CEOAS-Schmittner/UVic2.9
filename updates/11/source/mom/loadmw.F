      subroutine loadmw (joff, js, je, is, ie, num1, nu, first_mw)

#if defined O_mom
!=======================================================================
!     load variables into the MW for rows "js" through "je"

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       num1 = "tau-1" latitude disk unit
!       nu   = "tau" latitude disk unit
!=======================================================================

      implicit none

      integer istrt, is, iend, ie, j, js, je, jrow, joff, k, i, num1
      integer nu, n

      logical first_mw

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
      include "tmngr.h"

      real taum1old
      save taum1old

!-----------------------------------------------------------------------
!     limit the longitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     for all MW`s after the first, move the MW northward by copying
!     data from the last two rows into the first two rows.
!     (last 3 rows into first 3 rows if using fourth order schemes)
!-----------------------------------------------------------------------

      if (.not. first_mw) then
        call movemw (istrt-1, iend+1)
      endif

!-----------------------------------------------------------------------
!     construct "t" cell and "u" cell land/sea masks
!-----------------------------------------------------------------------

      if (.not. wide_open_mw .or. (wide_open_mw .and. first)) then
        do j=js,je
          jrow = j + joff
          do k=1,km
            do i=istrt-1,iend+1
              if (kmt(i,jrow) .ge. k) then
                tmask(i,k,j) = c1
              else
                tmask(i,k,j) = c0
              endif
              if (kmu(i,jrow) .ge. k) then
                umask(i,k,j) = c1
              else
                umask(i,k,j) = c0
              endif
            enddo
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     read data from "tau-1" and "tau" disk units into the MW
!     or if MW is wide open, copy data to proper time level and add
!     external mode to "tau" only since "tau-1" was done last timestep
!-----------------------------------------------------------------------

      if (wide_open_mw) then
        if (leapfrog) then
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
          if (first) then
            call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
          endif
        elseif (euler2) then
          tau   = taup1
          taup1 = taum1old
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
        elseif (euler1 .or. forward) then
          taum1old = taum1
          taum1    = tau
          call add_ext_mode (joff, js, je, istrt, iend, 'tau')
        endif
      else
        if (leapfrog .or. euler2) then

!         read "tau","tau-1" disk data into "tau","tau-1" MW positions

          call getvar (joff, js, je, istrt-1, iend+1, num1, nu)
        elseif (forward .or. euler1) then

!         read "tau" disk data into "tau" and "tau-1" MW positions

          call getvar (joff, js, je, istrt-1, iend+1, nu, nu)
        endif

!       add external mode to both since only internal modes are on disk

        call add_ext_mode (joff, js, je, istrt, iend, 'tau-1')
        call add_ext_mode (joff, js, je, istrt, iend, 'tau')
      endif
# if defined O_symmetry

!-----------------------------------------------------------------------
!     set symmetry condition at latitude of "u" cell row jmtm1
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        if (jrow .eq. jmt) then
          do k=1,km
            do i=istrt-1,iend+1
              u(i,k,j,1,tau)     = u(i,k,j-2,1,tau)
              u(i,k,j,1,taum1)   = u(i,k,j-2,1,taum1)
              u(i,k,j-1,2,tau)   = c0
              u(i,k,j-1,2,taum1) = c0
              u(i,k,j,2,tau)     = -u(i,k,j-2,2,tau)
              u(i,k,j,2,taum1)   = -u(i,k,j-2,2,taum1)
            enddo
          enddo
          do n=1,nt
            do k=1,km
              do i=istrt-1,iend+1
                t(i,k,j,n,tau)   = t(i,k,j-1,n,tau)
                t(i,k,j,n,taum1) = t(i,k,j-1,n,taum1)
              enddo
            enddo
          enddo
        endif
      enddo
# endif

!-----------------------------------------------------------------------
!     compute density  at "t" cell centers
!-----------------------------------------------------------------------

      call state (t(1,1,1,1,tau), t(1,1,1,2,tau), rho(1,1,jsmw)
     &,           max(jsmw,js), je, istrt-1, iend+1)

!-----------------------------------------------------------------------
!       open boundaries: extrapolate tracer and velocity values onto
!       boundary to prevent diffusion
!-----------------------------------------------------------------------

      return
      end

      subroutine movemw (iss, iee)

!=======================================================================
!     move the MW up (northward) by copying data from the last two rows
!     into the first two rows.
!     (last 3 rows if using fourth order schemes)

!     input:
!      is = starting longitude index in the MW
!      ie = ending longitude index in the MW
!      note: iss,iee bypassed to optimize performance
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, is, ie, nrows, move, jfrom, jto, n
      integer iee, iss

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"
# if defined O_isopycmix || defined O_smagnlmix
      include "hmixc.h"
# endif
      include "vmixc.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif
      parameter (is=1, ie=imt)

      nrows = jmw - ncrows
      do move=1,nrows
        jfrom = jmw - (nrows - move)
        jto   = move

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jmw)
!-----------------------------------------------------------------------

        do k=1,km
          do i=is,ie
            do n=1,nt
              t(i,k,jto,n,taum1) = t(i,k,jfrom,n,taum1)
              t(i,k,jto,n,tau)   = t(i,k,jfrom,n,tau)
            enddo
            do n=1,2
              u(i,k,jto,n,taum1) = u(i,k,jfrom,n,taum1)
              u(i,k,jto,n,tau)   = u(i,k,jfrom,n,tau)
            enddo
            tmask(i,k,jto)   = tmask(i,k,jfrom)
            umask(i,k,jto)   = umask(i,k,jfrom)
            adv_vnt(i,k,jto) = adv_vnt(i,k,jfrom)
          enddo
        enddo

        do i=is,ie
          do n=1,nt
            stf(i,jto,n) = stf(i,jfrom,n)
            btf(i,jto,n) = btf(i,jfrom,n)
          enddo
          smf(i,jto,1) = smf(i,jfrom,1)
          smf(i,jto,2) = smf(i,jfrom,2)
          bmf(i,jto,1) = bmf(i,jfrom,1)
          bmf(i,jto,2) = bmf(i,jfrom,2)
        enddo

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jemw)
!-----------------------------------------------------------------------

        if (jfrom .le. jemw) then

          do k=1,km
            do i=is,ie
              adv_vnu(i,k,jto) = adv_vnu(i,k,jfrom)
            enddo
          enddo
# if defined O_fct

!         row "jmw" is not used when "jmw" < "jmt" so it is not copied
!         treat these arrays as if dimensioned (1:jemw)
          do n=1,nt
            do k=1,km
              do i=is,ie
                anti_fn(i,k,jto,n)  = anti_fn(i,k,jfrom,n)
                R_plusY(i,k,jto,n)  = R_plusY(i,k,jfrom,n)
                R_minusY(i,k,jto,n) = R_minusY(i,k,jfrom,n)
#  if defined O_fct_3d
                R_plus3(i,k,jto,n)  = R_plus3(i,k,jfrom,n)
                R_minus3(i,k,jto,n) = R_minus3(i,k,jfrom,n)
#  endif
              enddo
            enddo
          enddo
# endif
# if defined O_biharmonic

!         row "jmw" is not used when "jmw" < "jmt" so it is not copied
!         treat del2 as if dimensioned (1:jemw)
          do n=1,nvarbh
            do k=1,km
              do i=is,ie
                del2(i,k,jto,n) = del2(i,k,jfrom,n)
              enddo
            enddo
          enddo
# endif
# if defined O_anisotropic_viscosity
          do k=1,km
             do i=is,ie
                visc_cnu(i,k,jto) = visc_cnu(i,k,jfrom)
                visc_ceu(i,k,jto) = visc_ceu(i,k,jfrom)
                amc_north(i,k,jto) = amc_north(i,k,jfrom)
                amc_south(i,k,jto) = amc_south(i,k,jfrom)
             enddo
          enddo
# endif
# if defined O_smagnlmix && !defined O_consthmix

          do k=1,km
            do i=is,ie
              strain(i,k,jto,1)  = strain(i,k,jfrom,1)
              strain(i,k,jto,2)  = strain(i,k,jfrom,2)
              am_lambda(i,k,jto) = am_lambda(i,k,jfrom)
              am_phi(i,k,jto)    = am_phi(i,k,jfrom)
              visc_cnu(i,k,jto)  = visc_cnu(i,k,jfrom)
              diff_cnt(i,k,jto)  = diff_cnt(i,k,jfrom)
            enddo
          enddo
# endif

        endif

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (jsmw:jmw)
!-----------------------------------------------------------------------

        if (jto .ge. jsmw) then
          do k=1,km
            do i=is,ie
              adv_vet(i,k,jto) = adv_vet(i,k,jfrom)
              rho(i,k,jto)  = rho(i,k,jfrom)
            enddo
          enddo
          do k=0,km
            do i=is,ie
              adv_vbt(i,k,jto) = adv_vbt(i,k,jfrom)
            enddo
          enddo
        endif

# if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (jsme:jemw)
!-----------------------------------------------------------------------

        if (jto .ge. jsmw .and. jfrom .le. jemw) then
          do k=1,km
            do i=is,ie
              visc_cbu(i,k,jto)  = visc_cbu(i,k,jfrom)
#  if !defined O_ppvmix
              diff_cbt(i,k,jto)  = diff_cbt(i,k,jfrom)
#  endif
              adv_veu(i,k,jto) = adv_veu(i,k,jfrom)
            enddo
          enddo

#  if defined O_fct
!         row "jmw" is not used when "jmw" < "jmt" so it is not copied
!         treat these arrays as if dimensioned (jsmw:jemw)
          do n=1,nt
            do k=1,km
              do i=is,ie
                anti_fe(i,k,jto,n)  = anti_fe(i,k,jfrom,n)
#   if defined O_fct_dlm2 && !defined O_fct_dlm1
                t_lo(i,k,jto,n) = t_lo(i,k,jfrom,n)
#   endif
              enddo
            enddo
            do k=0,km
              do i=is,ie
                anti_fb(i,k,jto,n) = anti_fb(i,k,jfrom,n)
              enddo
            enddo
          enddo
#  endif

          do k=0,km
            do i=is,ie
              adv_vbu(i,k,jto) = adv_vbu(i,k,jfrom)
            enddo
          enddo
        endif

#  if defined O_fourth_order_tracer_advection || defined O_quicker

!-----------------------------------------------------------------------
!       special case: copy advective flux as if dimensioned (1:jemw-1)
!-----------------------------------------------------------------------

        if (jto .ge. 1 .and. jfrom .le. jemw-1) then
!         row "jemw" is not used when "jmw" < "jmt" so it is not copied
!         treat these arrays as if dimensioned (1:jemw-1)
          do k=1,km
            do i=is,ie
              do n=1,nt
                adv_f4n(i,k,jto,n) = adv_f4n(i,k,jfrom,n)
              enddo
            enddo
          enddo
        endif
#  endif
# if defined O_pressure_gradient_average

!-----------------------------------------------------------------------
!       special case: copy computed tracers from extra row
!-----------------------------------------------------------------------

        if (jto .eq. 2) then
          do k=1,km
            do i=is,ie
              do n=1,nt
                t(i,k,jto,n,taup1) = t(i,k,jfrom,n,taup1)
              enddo
            enddo
          enddo
        endif
#  endif
# endif
# if defined O_isopycmix

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jmw)
!-----------------------------------------------------------------------

        do k=1,km
          do i=is,ie
            alphai(i,k,jto) = alphai(i,k,jfrom)
            betai(i,k,jto)  = betai(i,k,jfrom)
#  if defined O_full_tensor
            ddxt(i,k,jto,1) = ddxt(i,k,jfrom,1)
            ddxt(i,k,jto,2) = ddxt(i,k,jfrom,2)
#  endif
          enddo
        enddo
        do k=0,km
          do i=is,ie
            ddzt(i,k,jto,1) = ddzt(i,k,jfrom,1)
            ddzt(i,k,jto,2) = ddzt(i,k,jfrom,2)
          enddo
        enddo

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jemw)
!-----------------------------------------------------------------------

        if (jfrom .le. jemw) then
          do k=1,km
            do i=is,ie
              K22(i,k,jto)    = K22(i,k,jfrom)
              ddyt(i,k,jto,1) = ddyt(i,k,jfrom,1)
              ddyt(i,k,jto,2) = ddyt(i,k,jfrom,2)
              do kr=0,1
                do jq=0,1
                  Ai_nz(i,k,jto,jq,kr) = Ai_nz(i,k,jfrom,jq,kr)
                enddo
              enddo
#  if defined O_gent_mcwilliams
              adv_vntiso(i,k,jto)  = adv_vntiso(i,k,jfrom)
#  endif
            enddo
          enddo
        endif
#  if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (jsme:jemw)
!-----------------------------------------------------------------------

        if (jto .ge. jsmw .and. jfrom .le. jemw) then
          do k=1,km
            do i=is,ie
              do kr=0,1
                do ip=0,1
                  Ai_ez(i,k,jto,ip,kr)  = Ai_ez(i,k,jfrom,ip,kr)
                enddo
              enddo
              do kr=0,1
                do ip=0,1
                  Ai_bx(i,k,jto,ip,kr) = Ai_bx(i,k,jfrom,ip,kr)
                enddo
              enddo
              do kr=0,1
                do jq=0,1
                  Ai_by(i,k,jto,jq,kr) = Ai_by(i,k,jfrom,jq,kr)
                enddo
              enddo
              K11(i,k,jto)   = K11(i,k,jfrom)
              K33(i,k,jto)   = K33(i,k,jfrom)
#   if !defined O_full_tensor
              ddxt(i,k,jto,1) = ddxt(i,k,jfrom,1)
              ddxt(i,k,jto,2) = ddxt(i,k,jfrom,2)
#   endif
#   if defined O_gent_mcwilliams
              adv_vetiso(i,k,jto) = adv_vetiso(i,k,jfrom)
#   endif
            enddo
          enddo
#   if defined O_gent_mcwilliams
          do k=0,km
            do i=is,ie
              adv_vbtiso(i,k,jto) = adv_vbtiso(i,k,jfrom)
              adv_fbiso(i,k,jto)  = adv_fbiso(i,k,jfrom)
            enddo
          enddo
#   endif
        endif
#  endif
# endif
# if defined O_ppvmix

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jmw)
!-----------------------------------------------------------------------

        do k=1,km
          do i=is,ie
            rhom1z(i,k,jto)  = rhom1z(i,k,jfrom)
            uzsq(i,k,jto)    = uzsq(i,k,jfrom)
          enddo
        enddo

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (1:jemw)
!-----------------------------------------------------------------------

        if (jfrom .le. jemw) then
!         nothing to move
        endif

!-----------------------------------------------------------------------
!       copy quantities with rows dimensioned (jsmw:jmw)
!-----------------------------------------------------------------------

        if (jto .ge. jsmw) then
          do k=1,km
            do i=is,ie
              diff_cbt(i,k,jto)  = diff_cbt(i,k,jfrom)
              visc_cbt(i,k,jto)  = visc_cbt(i,k,jfrom)
            enddo
          enddo
        endif

# endif

      enddo

      return
      end

      subroutine getvar (joff, js, je, is, ie, num1, nu)

!=======================================================================
!     read prognostic quantities from disk units "num1" (tau-1) and
!     "nu" (tau) into the MW for rows "js" through "je"

!     input:
!       joff = offset between "j" in the MW and latitude "jrow"
!       js = starting row in the MW
!       je = ending row in the MW
!       is = starting longitude index in the MW
!       ie = ending longitude index in the MW
!       num1 = "tau-1" latitude disk unit
!       nu   = "tau" latitude disk unit
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, num1, nu, k, i, is, ie, n

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

      do j=js,je
        jrow   = j + joff

!       read "tau-1" data into the "tau-1" portion of the MW

        call getrow (num1, nslab, jrow, u(1,1,j,1,taum1)
     &,                                 t(1,1,j,1,taum1))
        if (nu .ne. num1) then

!         read "tau" data into the "tau" portion of the MW

          call getrow (nu, nslab, jrow, u(1,1,j,1,tau)
     &,                                 t(1,1,j,1,tau))
        else

!         copy "tau" data into "tau-1" portion of the MW

          do k=1,km
            do i=is,ie
              u(i,k,j,1,tau) = u(i,k,j,1,taum1)
              u(i,k,j,2,tau) = u(i,k,j,2,taum1)
            enddo
          enddo
          do n=1,nt
            do k=1,km
              do i=is,ie
                t(i,k,j,n,tau) = t(i,k,j,n,taum1)
              enddo
            enddo
          enddo
        endif
      enddo

      return
      end

      subroutine add_ext_mode (joff, js, je, is, ie, timelev)

!=======================================================================
!     add external mode to velocity for time level "timelev"

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!       timelev = "tau" or "tau-1"
!=======================================================================

      implicit none

      character(*) :: timelev

      integer j, js, je, jrow, joff, i, is, ie, n, k

      real diag1, diag0

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "grdvar.h"
      include "iounit.h"
      include "mw.h"
      include "switch.h"

# if defined O_stream_function
      real ext(imt,1:jmw,2)
# endif

      if (timelev .eq. 'tau') then

!-----------------------------------------------------------------------
!       add external mode to "tau" velocity
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
# if defined O_stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,1) - psi(i  ,jrow,1)
              diag0       = psi(i  ,jrow+1,1) - psi(i+1,jrow,1)
              ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
              ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ext(i,j,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,tau) = (u(i,k,j,n,tau) + ubar(i,jrow,n))
     &                             *umask(i,k,j)
                enddo
              enddo
            enddo
# endif
            do n=1,2
              call setbcx (u(1,1,j,n,tau), imt, km)
            enddo
          endif
        enddo

      elseif (timelev .eq. 'tau-1') then

!-----------------------------------------------------------------------
!       add external mode to "tau-1" velocity
!-----------------------------------------------------------------------

        do j=js,je
          jrow = j + joff
          if (jrow .lt. jmt) then
# if defined O_stream_function
            do i=is,ie
              diag1       = psi(i+1,jrow+1,2) - psi(i  ,jrow,2)
              diag0       = psi(i  ,jrow+1,2) - psi(i+1,jrow,2)
              ext(i,j,1)  = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
              ext(i,j,2)  =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)
     &                       *csur(jrow)
            enddo
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1) + ext(i,j,n))
     &                               *umask(i,k,j)
                enddo
              enddo
            enddo
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
            do n=1,2
              do k=1,km
                do i=is,ie
                  u(i,k,j,n,taum1) = (u(i,k,j,n,taum1)+ubarm1(i,jrow,n))
     &                           *umask(i,k,j)
                enddo
              enddo
            enddo
# endif
            do n=1,2
              call setbcx (u(1,1,j,n,taum1), imt, km)
            enddo
          endif
        enddo
      else
        write (stdout,'(/a,a/)') 'Error: timelev = ',timelev
        stop '=>loadmw'
      endif

      return
      end

      subroutine putmw (joff, js, je, nup1)

!=======================================================================
!     write prognostic quantities from MW to disk unit "nup1" (tau+1)
!     for rows "js" to "je".
!=======================================================================

      implicit none

      integer j, js, je, jrow, joff, nup1

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "mw.h"

!     write all newly computed quantities to disk "tau+1"

      do j=js,je
        jrow   = j + joff
        call putrow (nup1, nslab, jrow, u(1,1,j,1,taup1)
     &,                                 t(1,1,j,1,taup1))
      enddo
#endif

      return
      end
