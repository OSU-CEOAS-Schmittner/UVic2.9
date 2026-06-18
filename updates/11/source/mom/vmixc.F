      subroutine vmixc (joff, js, je, is, ie)

#if defined O_mom
!=======================================================================
!     set viscosity coefficient on bottom face of "u" cells
!     set diffusion coefficient on bottom face of "t" cells

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, istrt, is, iend, ie, jstrt
      integer jend, jrow, joff, ks, k1

      real zn2, hab, zkappa, qk1, qo1, edr, q2

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"
# if defined O_isopycmix || defined O_redi_diffusion
      include "isopyc.h"
#  if defined O_tidal_kv
      include "tidal_kv.h"
      include "diag.h"
      include "grdvar.h"
      include "levind.h"
#  endif
#  if defined O_save_kv
      include "timeavgs.h"
#  endif
# endif
# if defined O_constvmix && defined O_implicitvmix

      real temp(imt,km,jsmw:jmw)
# endif

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude and latitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)
      jstrt = max(2,js-1)
      jend  = je-1
# if defined O_constvmix

!-----------------------------------------------------------------------
!     constant vertical mixing coefficients
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = j + joff
        do i=istrt,iend
          if (abs(tlat(i,jrow)) .lt. 30.) then 
	    qk1 = 0.33
	    qo1 = 0.33
	  else
	    qk1 = 1.
	    qo1 = 1.
	  endif
          if (abs(tlat(i,jrow)) .lt. 70.) then 
            q2 = 0.33
	  else
            q2 = 1.
	  endif
c          q2=0.5
          do k=1,kmt(i,jrow)-1
            visc_cbu(i,k,j) = kappa_m
#  if defined O_tidal_kv && defined O_isopycmix

!           calculate N^2 = -g/rho drhodz on bottom of cell face
!           (where K33 and diff_cbt = kappa_h are defined). Note that
!           N2 is not guaranteed to be positive. If instability occurs,
!           convective adjustment will eliminate it.
!           drodzb is defined in isopyc.h

!           ZN2 is defined on T-cell bottom (zw pt)
            ZN2 = max(-gravrho0r*drodzb(i,k,j,0),1e-8)

!           subgrid-scale scheme: sum over all levels below k
            edr = 0.
            do k1=k+1,kmt(i,jrow)
!              height above bottom
               hab = zw(k) - zw(k1)
               edr = edr + (q2*(edrm2(i,k1,jrow)+edrs2(i,k1,jrow))
     &                   + qk1*edrk1(i,k1,jrow) + qo1*edro1(i,k1,jrow))
     &                   *exp(hab*zetar)/(1-exp(-zetar*zw(k1)))
            enddo
            zkappa = ogamma*edr/ZN2

! Andreas: high mixing in Southern Ocean below 500 m as in observations
c            zkappa = max(zkappa,tanh((zw(k)-50000.)/10000)*
c     &                          (1-tanh((tlat(i,j+joff)+40)/8))/2)
!           limit diff_cbt
#   if defined O_bryan_lewis_vertical
            diff_cbt(i,k,j) = max(Ahv(k), min(100., zkappa + Ahv(k)))
#   else
            diff_cbt(i,k,j) = max(kappa_h, min(100., zkappa + kappa_h))
#   endif
#  elif defined O_bryan_lewis_vertical
            diff_cbt(i,k,j) = Ahv(k)
#  else
            diff_cbt(i,k,j) = kappa_h
#  endif
          enddo
        enddo
      enddo
#  if defined O_implicitvmix

      do ks=1,2

!       find density

        call statec (t(1,1,1,1,taum1), t(1,1,1,2,taum1), temp(1,1,jsmw)
     &,                jstrt, jend, istrt, iend, ks)

!       set large diffusion coefficient between unstable layers
!       (note: viscosity coefficient is not limited but could be here)

        do j=jstrt,jend
          do k=ks,kmm1,2
            do i=istrt,iend
              if (temp(i,k,j) .gt. temp(i,k+1,j)) then
                diff_cbt(i,k,j) = diff_cbt_limit*tmask(i,k+1,j)
              endif
            enddo
          enddo
        enddo
      enddo
#  endif
# endif
# if defined O_ppvmix

!-----------------------------------------------------------------------
!     for momentum and tracers based on the pacanowski & philander
!     richardson mixing scheme (JPO vol 11, # 11, 1981).
!-----------------------------------------------------------------------

      call ppmix (joff, js, je, istrt, iend)
# endif
# if defined O_save_kv

!-----------------------------------------------------------------------
!     accumulate time average diapycnal (without K33) diffusivity
!-----------------------------------------------------------------------

      if (timavgperts .and. .not. euler2) then
        do j=jstrt,jend
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              ta_diff_cbt(i,k,jrow) = ta_diff_cbt(i,k,jrow)
     &                              + diff_cbt(i,k,j)
            enddo
          enddo
        enddo
      endif
# endif
# if defined O_isopycmix || defined O_redi_diffusion

!-----------------------------------------------------------------------
!     Add K33 component to vertical diffusion coefficient
!-----------------------------------------------------------------------

      do j=jstrt,jend
        do i=istrt,iend
          do k=1,km
            diff_cbt(i,k,j) = diff_cbt(i,k,j) + K33(i,k,j)
          enddo
        enddo
      enddo
# endif
#endif

      return
      end
