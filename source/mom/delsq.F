      subroutine delsq (joff, js, je, is, ie)

#if defined O_mom && defined O_biharmonic && !defined O_bryan_lewis_horizontal
!=======================================================================
!     compute del**2 of prognostic variables on MW rows "js" ... "je"
!=======================================================================

      implicit none

      integer i, k, j, n, ip, jp, jrow, js, je, jstrt, jend, m, joff
      integer js, je, jstrt, jend, m, joff, is, ie

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_f4n, adv_txiso
      real adv_vetiso, adv_tyiso, adv_vntiso, adv_tziso, adv_fbiso
      real diff_tx, diff_ty, diff_tz, aidif, adv_ux, adv_uy, adv_uz
      real adv_metric, diff_ux, diff_uz, diff_uy, unep, diff_metric
      real del2, coriolis, ahbi_cstr, ahbi_csu_dyur, ambi_csur
      real ambi_cst_dytr

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "fdift.h"
      include "fdifm.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the MW row indices to what can be calculated
!-----------------------------------------------------------------------

      jstrt = max(2,js-1)
      jend  = je-1

!-----------------------------------------------------------------------
!     compute del**2 of tracers
!-----------------------------------------------------------------------

      do n=1,nt
        m = n+2

!       diffusive flux across eastern face of "T" cell
!       diffusive flux across northern face of "T" cell

        do j=jstrt,jend
          jrow = j + joff
          ahbi_cstr = diff_cet*cstr(jrow)
          do k=1,km
            do i=is-1,ie
              diff_fe(i,k,j) = ahbi_cstr*dxur(i)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=jstrt-1,jend
          jrow = j + joff
          ahbi_csu_dyur = diff_cnt*csu(jrow)*dyur(jrow)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ahbi_csu_dyur*
     &                     (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo

!       compute -ah*del**2 of tracer

        do j=jstrt,jend
          jrow = j + joff
          do i=is-1,ie
            cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          enddo
        enddo
        do j=jstrt,jend
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              del2(i,k,j,m) = -DIFF_Tx(i,k,j) - DIFF_Ty(i,k,j,jrow,n)
            enddo
          enddo
          call setbcx (del2(1,1,j,m), imt, km)

!         set -del**2 = 0 on southern wall

          if (jrow .eq. 2) then
            do k=1,km
              do i=is-1,ie+1
                del2(i,k,j-1,m) = c0
              enddo
            enddo
          endif

!         set -del**2 = 0 northern wall

          if (jrow .eq. jmt-1) then
            do k=1,km
              do i=is-1,ie+1
                del2(i,k,j+1,m) = c0
              enddo
            enddo
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute del**2 of velocities
!-----------------------------------------------------------------------

      do n=1,2

!       diffusive flux across east face of "u" cell
!       diffusive flux across north face of "u" cell

        do j=jstrt,jend
          jrow = j + joff
          ambi_csur = visc_ceu*csur(jrow)
          do k=1,km
            do i=is-1,ie
              diff_fe(i,k,j) = ambi_csur*dxtr(i+1)*
     &                         (u(i+1,k,j,n,taum1) - u(i,k,j,n,taum1))
            enddo
          enddo
        enddo
        do j=jstrt-1,jend
          jrow = j + joff
          ambi_cst_dytr = visc_cnu*cst(jrow+1)*dytr(jrow+1)
          do k=1,km
            do i=is,ie
              diff_fn(i,k,j) = ambi_cst_dytr*
     &                          (u(i,k,j+1,n,taum1) - u(i,k,j,n,taum1))
            enddo
          enddo
        enddo

!       compute -am*del**2

        do j=jstrt,jend
          jrow = j + joff
          do i=is-1,ie
            csudxur(i,j) = csur(jrow)*dxur(i)
          enddo
        enddo
        do j=jstrt,jend
          jrow   = j + joff
          do k=1,km
            do i=is,ie
              del2(i,k,j,n) = (-DIFF_Ux(i,k,j) - DIFF_Uy(i,k,j,jrow,n)
     &                      - am3(jrow)*u(i,k,j,n,taum1) - am4(jrow,n)*
     &           (u(i+1,k,j,3-n,taum1) - u(i-1,k,j,3-n,taum1))*dxmetr(i)
     &                        )*umask(i,k,j)
            enddo
          enddo
          call setbcx (del2(1,1,j,n), imt, km)

!         set -del**2 = 0 on southern wall

          if (jrow .eq. 2) then
            do k=1,km
              do i=is-1,ie+1
                del2(i,k,j-1,n) = c0
              enddo
            enddo
          endif

!         set -del**2 = 0 on northern wall

          if (jrow .eq. jmt-1) then
            do k=1,km
              do i=is-1,ie+1
                del2(i,k,j+1,n) = c0
              enddo
            enddo
          endif
        enddo
      enddo
#endif

      return
      end
