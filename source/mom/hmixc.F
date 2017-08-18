      subroutine hmixc (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     set horizontal mixing coeffs on north and east face of "t" and
!     "u" cells.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer js ,je ,jrowstart ,joff ,jrowend ,jrow ,jm1 ,jp1 ,k ,ie
      integer is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
# if defined O_full_tensor
      include "vmixc.h"
# endif

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     set all horizontal mixing coefficients
!-----------------------------------------------------------------------

      jrowstart = js + joff
      jrowend   = min(je + joff + 2,jmt)

# if defined O_consthmix

!     for momentum... set coefficients for all latitudes

      if (first) then
#  if defined O_biharmonic
        visc_cnu = sqrt(abs(ambi))
        visc_ceu = sqrt(abs(ambi))
#  else
        visc_cnu  = am
        visc_ceu  = am
#  endif
        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          amc_north(jrow) = visc_cnu*cst(jp1)*dytr(jp1)
     &                              *csur(jrow)*dyur(jrow)
          amc_south(jrow) = visc_cnu*cst(jrow)*dytr(jrow)
     &                              *csur(jrow)*dyur(jrow)
        enddo
      endif
# endif
# if defined O_consthmix

!     for tracers... set coefficients for all latitudes

      if (first) then
#  if defined O_bryan_lewis_horizontal
        do k=1,km
          diff_cet(k) = Ahh(k)
          diff_cnt(k) = Ahh(k)
        enddo
        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          do k=1,km
            ahc_north(jrow,k) = diff_cnt(k)*csu(jrow)*dyur(jrow)
     &                                     *cstr(jrow)*dytr(jrow)
            ahc_south(jrow,k) = diff_cnt(k)*csu(jm1)*dyur(jm1)
     &                                     *cstr(jrow)*dytr(jrow)
          enddo
        enddo
#  else
#   if defined O_biharmonic
        diff_cet = sqrt(abs(ahbi))
        diff_cnt = sqrt(abs(ahbi))
#   else
        diff_cnt  = ah
        diff_cet  = ah
#   endif
#   if defined O_full_tensor
        diff_cnt = diff_cnt + kappa_h
        diff_cet = diff_cet + kappa_h
#   endif

        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          ahc_north(jrow) = diff_cnt*csu(jrow)*dyur(jrow)*cstr(jrow)
     &                              *dytr(jrow)
          ahc_south(jrow) = diff_cnt*csu(jm1)*dyur(jm1)*cstr(jrow)
     &                              *dytr(jrow)
        enddo
#  endif
      endif
# else
#  if defined O_smagnlmix && !defined O_consthmix

!-----------------------------------------------------------------------
!     limit the longitude indices and calculate deformation rates
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

      call smagnlc (joff, js, je, istrt, iend)
#  endif
# endif
#endif

      return
      end
