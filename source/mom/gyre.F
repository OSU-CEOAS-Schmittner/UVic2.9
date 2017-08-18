      subroutine gyre (joff, js, je, is, ie, n)

#if defined O_mom && defined O_gyre_components
!-----------------------------------------------------------------------
!     compute the northward transport components of each tracer

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, is, ie, n
      integer mask, ll

      real small, totdxn, totdxs, vbr, tbrs, tbrn, tempdiff_fn
      real tempadv_fn, factor, totz, vbrz, tbrz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "cregin.h"
      include "diag.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
# if defined O_isopycmix
      include "isopyc.h"
# endif

      do j=js,je
        jrow = j + joff
        if (jrow .lt. jmtm1) then
          small = 1.e-10
          do k=1,km
            totdxn = small
            totdxs = small
            vbr    = c0
            tbrs   = c0
            tbrn   = c0
            do i=is,ie
              totdxn = totdxn + dxt(i)*tmask(i,k,j+1)
              totdxs = totdxs + dxt(i)*tmask(i,k,j)
              vbr    = vbr  + u(i,k,j,2,tau)*dxu(i)*csu(jrow)
              tbrn   = tbrn + t(i,k,j+1,n,tau)*tmask(i,k,j+1)*dxt(i)
              tbrs   = tbrs + t(i,k,j,n,tau)*tmask(i,k,j)*dxt(i)
            enddo
            tbrn          = tbrn/totdxn
            tbrs          = tbrs/totdxs
            ttn(1,jrow,n) = ttn(1,jrow,n) + vbr*p5*(tbrn+tbrs)*dzt(k)
            do i=is,ie
              tempdiff_fn =
# if defined O_consthmix && !defined O_biharmonic && !defined O_isopycmix
#  if defined O_bryan_lewis_horizontal
     &                      ahc_north(jrow,k)*cst(jrow)*dyt(jrow)*
     &                      (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))*
#  else
     &                      ahc_north(jrow)*cst(jrow)*dyt(jrow)*
     &                      (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))*
#  endif
# else
     &                      diff_fn(i,k,j)*
# endif
     &                           tmask(i,k,j+1)*tmask(i,k,j)*
     &                           dxt(i)*dzt(k)
              tempadv_fn       = p5*adv_vnt(i,k,j)*(t(i,k,j,n,tau) +
     &                           t(i,k,j+1,n,tau))*dxt(i)*dzt(k)
              ttn(6,jrow,n)    = ttn(6,jrow,n) + tempadv_fn
              ttn(7,jrow,n)    = ttn(7,jrow,n) - tempdiff_fn
              ttn2(6,jrow,n,0) = ttn2(6,jrow,n,0) + tempadv_fn
              ttn2(7,jrow,n,0) = ttn2(7,jrow,n,0) - tempdiff_fn
              if (mskhr(i,jrow) .ne. 0) then
                ttn2(6,jrow,n,mskhr(i,jrow)) =
     &                   ttn2(6,jrow,n,mskhr(i,jrow)) + tempadv_fn
                ttn2(7,jrow,n,mskhr(i,jrow)) =
     &                   ttn2(7,jrow,n,mskhr(i,jrow)) - tempdiff_fn
              endif
            enddo
          enddo

          do i=is,ie
            if (cori(i,jrow,1) .eq. c0 .and. jrow .gt. 1) then
              factor = c4*cori(i,jrow-1,1)
            else
              factor = c4*cori(i,jrow,1)
            endif
            totz = c0
            vbrz = c0
            tbrz = c0
            do k=1,km
              mask = tmask(i,k,j)*tmask(i,k,j+1)
              vbrz = vbrz + adv_vnt(i,k,j)*dxt(i)*dzt(k)
              tbrz = tbrz +mask*(t(i,k,j,n,tau)+t(i,k,j+1,n,tau))*dzt(k)
              totz = totz + mask*dzt(k)
            enddo
            if (totz .ne. c0) then
              tbrz = tbrz/totz
              ttn(3,jrow,n) = ttn(3,jrow,n) + vbrz*tbrz*p5
              ttn(5,jrow,n) = ttn(5,jrow,n) - (smf(i,j,1)*dxu(i) +
     &                        smf(i-1,j,1)*dxu(i-1))*(t(i,1,j,n,tau)
     &                        +t(i,1,j+1,n,tau)-tbrz)
     &         *csu(jrow)/factor
            endif
          enddo
          ttn(2,jrow,n) = ttn(6,jrow,n)-ttn(1,jrow,n)
          ttn(4,jrow,n) = ttn(6,jrow,n)-ttn(3,jrow,n)-ttn(5,jrow,n)
          ttn(8,jrow,n) = ttn(6,jrow,n)+ttn(7,jrow,n)
          do ll=0,nhreg
            ttn2(8,jrow,n,ll) = ttn2(6,jrow,n,ll)+ttn2(7,jrow,n,ll)
          enddo
        endif
      enddo
#endif

      return
      end
