      subroutine setvbc (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     set momentum and tracer vertical boundary conditions

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer js, je, istrt, is, iend, ie, n, j, i, jrow, joff, kz

      real uvmag

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "grdvar.h"
      include "levind.h"
      include "scalar.h"
      include "mw.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

!----------------------------------------------------------------------
!       set no flux condition for all tracers at surface & bottom.
!----------------------------------------------------------------------

      do n=1,nt
        do j=js,je
          do i=istrt,iend
            stf(i,j,n) = c0
            btf(i,j,n) = c0
          enddo
        enddo
      enddo

!----------------------------------------------------------------------
!       apply surface tracer and momentum fluxes from the atmosphere
!       code is for 2 tracer and 2 momentum fluxes.
!----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do i=istrt,iend
          stf(i,j,itemp) = sbc(i,jrow,ihflx)*tmask(i,1,j)
          stf(i,j,isalt) = sbc(i,jrow,isflx)*tmask(i,1,j)
# if defined O_carbon
          stf(i,j,idic) = sbc(i,jrow,idicflx)*tmask(i,1,j)
#  if defined O_carbon_14
          stf(i,j,ic14) = sbc(i,jrow,ic14flx)*tmask(i,1,j)
#  endif
# endif
# if defined O_npzd_alk
          stf(i,j,ialk) = sbc(i,jrow,ialkflx)*tmask(i,1,j)
# endif
# if defined O_npzd_o2
          stf(i,j,io2) = sbc(i,jrow,io2flx)*tmask(i,1,j)
# endif
# if defined O_npzd
          stf(i,j,ipo4) = sbc(i,jrow,ipo4flx)*tmask(i,1,j)
#  if !defined O_npzd_no_vflux
          stf(i,j,iphyt) = sbc(i,jrow,iphytflx)*tmask(i,1,j)
          stf(i,j,izoop) = sbc(i,jrow,izoopflx)*tmask(i,1,j)
          stf(i,j,idetr) = sbc(i,jrow,idetrflx)*tmask(i,1,j)
#  endif
#  if defined O_npzd_nitrogen
          stf(i,j,ino3) = sbc(i,jrow,ino3flx)*tmask(i,1,j)
#   if !defined O_npzd_no_vflux
          stf(i,j,idiaz) = sbc(i,jrow,idiazflx)*tmask(i,1,j)
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
          stf(i,j,icfc11) = sbc(i,jrow,icfc11flx)*tmask(i,1,j)
          stf(i,j,icfc12) = sbc(i,jrow,icfc12flx)*tmask(i,1,j)
# endif
          smf(i,j,1) = sbc(i,jrow,itaux)*umask(i,1,j)
          smf(i,j,2) = sbc(i,jrow,itauy)*umask(i,1,j)
        enddo
      enddo

!----------------------------------------------------------------------
!       set bottom drag
!----------------------------------------------------------------------

      do n=1,2
        if (cdbot .eq. c0) then
          do j=js,je
            do i=istrt,iend
              bmf(i,j,n) = c0
            enddo
          enddo
        else
          do j=js,je
            jrow = j + joff
            do i=istrt,iend
              kz = kmu(i,jrow)
              if (kz .ne. 0) then
                uvmag    = sqrt(u(i,kz,j,1,taum1)**2 +
     &                          u(i,kz,j,2,taum1)**2)
                bmf(i,j,n) = cdbot*u(i,kz,j,n,taum1)*uvmag
              else
                bmf(i,j,n) = c0
              endif
            enddo
          enddo
        endif
      enddo

!----------------------------------------------------------------------
!     apply zonal boundary conditions
!----------------------------------------------------------------------

      do n=1,nt
        call setbcx (stf(1,js,n), imt, je-js+1)
        call setbcx (btf(1,js,n), imt, je-js+1)
      enddo
      do n=1,2
        call setbcx (smf(1,js,n), imt, je-js+1)
        call setbcx (bmf(1,js,n), imt, je-js+1)
      enddo

#endif

      return
      end
