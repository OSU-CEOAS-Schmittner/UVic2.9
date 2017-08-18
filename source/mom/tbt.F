      subroutine tbt1 (joff, js, je, n)

#if defined O_mom && defined O_mom_tbt
!=======================================================================
!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       n     = tracer component
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, term, r2dt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "grdvar.h"
      include "hmixc.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif
      include "levind.h"
      include "mw.h"
      include "xbt.h"
      include "scalar.h"
      include "tbt.h"
      include "vmixc.h"
      include "fdift.h"

      do j=js,je
        jrow = j + joff

!-----------------------------------------------------------------------
!       accumulate data for one time step
!-----------------------------------------------------------------------

        do i=2,imtm1
          do k=1,kmt(i,jrow)

!-----------------------------------------------------------------------
!           tracer
!-----------------------------------------------------------------------

            term = tmask(i,k,j)*t(i,k,j,n,tau)
!            tbt(i,jrow,k,n,15) = tbt(i,jrow,k,n,15) + term
            tbt(i,jrow,k,n,11) = tbt(i,jrow,k,n,11) + term

!-----------------------------------------------------------------------
!           d(tracer)/dt
!-----------------------------------------------------------------------

            r2dt = c1/(c2dtts*dtxcel(k))
            term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                               t(i,k,j,n,taum1))*r2dt
            tbt(i,jrow,k,n,9) = tbt(i,jrow,k,n,9) + term

!-----------------------------------------------------------------------
!           zonal advection (flux form) of tracer
!-----------------------------------------------------------------------

            term = -tmask(i,k,j)*ADV_Tx(i,k,j)
# if defined O_gent_mcwilliams && !defined O_fct
     &             -tmask(i,k,j)*ADV_Txiso(i,k,j,n)
# endif
            tbt(i,jrow,k,n,2) = tbt(i,jrow,k,n,2) + term

!-----------------------------------------------------------------------
!           pure zonal advection of tracer
!-----------------------------------------------------------------------

!           - U(T)x = T(U)x - (UT)x

!            dudx = (adv_vet(i,k,j)-adv_vet(i-1,k,j))*dxtr(i)*cstr(jrow)
# if defined O_gent_mcwilliams && !defined O_fct
!     &             + (adv_vetiso(i,k,j)-adv_vetiso(i-1,k,j))*dxtr(i)
!     &             *cstr(jrow)
# endif
!            term = tmask(i,k,j)*(t(i,k,j,n,tau)*dudx - ADV_Tx(i,k,j))
# if defined O_gent_mcwilliams && !defined O_fct
!     &             -tmask(i,k,j)*ADV_Txiso(i,k,j,n)
# endif
!            tbt(i,jrow,k,n,11) = tbt(i,jrow,k,n,11) + term

!-----------------------------------------------------------------------
!           meridional advection (flux form) of tracer
!-----------------------------------------------------------------------

            term = -tmask(i,k,j)*ADV_Ty(i,k,j,jrow,n)
# if defined O_gent_mcwilliams && !defined O_fct
     &             -tmask(i,k,j)*ADV_Tyiso(i,k,j,jrow,n)
# endif
            tbt(i,jrow,k,n,3) = tbt(i,jrow,k,n,3) + term

!-----------------------------------------------------------------------
!           pure meridional advection of tracer
!-----------------------------------------------------------------------

!           - V(T)y = T(V)y - (VT)y

!            dvdy = (adv_vnt(i,k,j)-adv_vnt(i,k,j-1))*dytr(jrow)
!     &             *cstr(jrow)
# if defined O_gent_mcwilliams && !defined O_fct
!     &             + (adv_vntiso(i,k,j)-adv_vntiso(i,k,j-1))*dytr(jrow)
!     &             *cstr(jrow)
# endif
!            term = tmask(i,k,j)*(t(i,k,j,n,tau)*dvdy
!     &           - ADV_Ty(i,k,j,jrow,n))
# if defined O_gent_mcwilliams && !defined O_fct
!     &             -tmask(i,k,j)*ADV_Tyiso(i,k,j,jrow,n)
# endif
!            tbt(i,jrow,k,n,12) = tbt(i,jrow,k,n,12) + term

!-----------------------------------------------------------------------
!           vertical advection (flux form) of tracer
!-----------------------------------------------------------------------

            term = -tmask(i,k,j)*ADV_Tz(i,k,j)
#  if defined O_gent_mcwilliams && !defined O_fct
     &             -tmask(i,k,j)*ADV_Tziso(i,k,j)
#  endif
            tbt(i,jrow,k,n,4) = tbt(i,jrow,k,n,4) + term

!-----------------------------------------------------------------------
!           pure vertical advection of tracer
!-----------------------------------------------------------------------

!           - W(T)z = T(W)z - (WT)z

!            dwdz = (adv_vbt(i,k-1,j)-adv_vbt(i,k,j))*dztr(k)
#  if defined O_gent_mcwilliams && !defined O_fct
!     &             + (adv_vbtiso(i,k-1,j)-adv_vbtiso(i,k,j))*dztr(k)
#  endif
!            term = tmask(i,k,j)*(t(i,k,j,n,tau)*dwdz - ADV_Tz(i,k,j))
#  if defined O_gent_mcwilliams && !defined O_fct
!     &             -tmask(i,k,j)*ADV_Tziso(i,k,j)
#  endif
!            tbt(i,jrow,k,n,13) = tbt(i,jrow,k,n,13) + term

!-----------------------------------------------------------------------
!           zonal diffusion of tracer
!-----------------------------------------------------------------------

            term = tmask(i,k,j)*DIFF_Tx(i,k,j)
            tbt(i,jrow,k,n,5) = tbt(i,jrow,k,n,5) + term

!-----------------------------------------------------------------------
!           meridional diffusion of tracer
!-----------------------------------------------------------------------

            term = tmask(i,k,j)*DIFF_Ty(i,k,j,jrow,n)
            tbt(i,jrow,k,n,6) = tbt(i,jrow,k,n,6) + term

!-----------------------------------------------------------------------
!           vertical diffusion of tracer
!-----------------------------------------------------------------------

            term = tmask(i,k,j)*DIFF_Tz(i,k,j)
# if defined O_implicitvmix || defined O_isopycmix
     &             + tmask(i,k,j)*zzi(i,k,j)
# endif
            tbt(i,jrow,k,n,7) = tbt(i,jrow,k,n,7) + term

!-----------------------------------------------------------------------
!           tracer source term
!-----------------------------------------------------------------------

# if defined O_source_term || defined O_npzd || defined O_carbon_14
            term = tmask(i,k,j)*source(i,k,j)
            tbt(i,jrow,k,n,8) = tbt(i,jrow,k,n,8) + term
# else
            tbt(i,jrow,k,n,8) = 0.0
# endif

          enddo

!-----------------------------------------------------------------------
!         surface tracer flux
!-----------------------------------------------------------------------

          tbtsf(i,jrow,n) = tbtsf(i,jrow,n) + tmask(i,1,j)*stf(i,j,n)
        enddo

        if (jrow .eq. 2 .and. n .eq. 1)  ntbtts = ntbtts + 1

      enddo

      return
      end

      subroutine tbt2 (joff, js, je, iterm)
!=======================================================================
!     Accumulate parts of d(tracer)/dt

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!=======================================================================

      implicit none

      integer iterm, j, js, je, jrow, joff, i, n, k

      real r2dt, term

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "tbt.h"

# if !defined O_fourfil && !defined O_firfil
      if (iterm .eq. 10) return
# endif
      if (iterm .ne. 1 .and. iterm .ne. 10) then
        write (stdout,*) '=>Error: iterm=',iterm,' in tbt2'
        stop '=>tbt2'
      endif

      do j=js,je
        jrow  = j + joff
        do i=2,imtm1
          do n=1,nt
            do k=1,kmt(i,jrow)

!-----------------------------------------------------------------------
!             d/dt(tracer)
!-----------------------------------------------------------------------
              r2dt = c1/(c2dtts*dtxcel(k))
              term = tmask(i,k,j)*(t(i,k,j,n,taup1) -
     &                             t(i,k,j,n,taum1))*r2dt
              tbt(i,jrow,k,n,iterm) = tbt(i,jrow,k,n,iterm) + term

!-----------------------------------------------------------------------
!             change in variance of tracer
!-----------------------------------------------------------------------
!              if (iterm .eq. 1) then
!                term = tmask(i,k,j)*(t(i,k,j,n,taup1)**2
!     &                 - t(i,k,j,n,taum1)**2)
!                tbt(i,jrow,k,n,14) = tbt(i,jrow,k,n,14) + term
!              endif

            enddo
          enddo
        enddo
      enddo

#endif

      return
      end
