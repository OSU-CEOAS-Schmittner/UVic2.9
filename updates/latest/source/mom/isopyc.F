
      subroutine isopi (error, am, ah)

#if defined O_mom
# if defined O_isopycmix || defined O_gent_mcwilliams
!=======================================================================

!               Initialization for isopycnal mixing scheme

!      -Disopycmix gives either the full or small angle isopycnal
!                  mixing tensor

!      -Disopycmix -Dgent_mcwilliams gives the isopycnal mixing tensor
!       plus the advective velocity parameterization of Gent_McWilliams
!     input:
!       error  = logical to signal problems
!       am     = background mixing coeff for momentum
!       ah     = horizontaal mixing coeff for tracers
!       slmx   = max slope of isopycnals
!       ahisop = isopycnal tracer diffusivity(cm**2/sec)
!       athkdf = isopycnal thickness diffusivity (cm**2/sec)

!     output:
!       The above input can be reset via namelist
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name
      integer i, k, j, ip, kr, jq, io, i_delta, j_delta, k_delta, jrow
      integer iou, n, ib(10), ic(10)
      logical exists, inqvardef, error, debug
      real slmx, s_dm, ah, am, ft, delta1, delta2

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "iounit.h"
      include "isopyc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

      real tmpijk(imtm2,jmtm2,km)

c      namelist /isopyc/ slmx, ahisop, athkdf, del_dm, s_dm, 
c     &                  c_eden, gamma_eden

      write (stdout,'(/,20x,a,/,20x,a,/)')
     & 'I S O P Y C M I X    I N I T I A L I Z A T I O N'
#  if defined O_gent_mcwilliams
     &,'Isopycnal mixing tensor + GM thickness diffusion'
#  else
     &,'Isopycnal mixing tensor.'
#  endif

!-----------------------------------------------------------------------
!     USER INPUT
!     initialize variables (all mixing units are cm**2/sec.)
!-----------------------------------------------------------------------

!     maximum isopycnal slope

      slmx  = 0.01

!     isopycnal diffusion coefficient

      ahisop = 1.2e7

!     isopycnal thickness diffusion coefficient

      athkdf = 8.e6

!     transition for scaling diffusion coefficient

      del_dm = .4e-02

!     half width scaling for diffusion coefficient

      s_dm = .1e-2

!-----------------------------------------------------------------------
!     provide for namelist over-ride of above settings + documentation
!-----------------------------------------------------------------------

c      call getunit (io, 'control.in'
c     &,               'formatted sequential rewind')
c      read (io,isopyc,end=100)
c100   continue
c      write (stdout,isopyc)
c      call relunit (io)

!     set reciprocal of maximum isopycnal slope and scaling half width
      slmxr = c1/slmx
      s_dmr = c1/s_dm

!-----------------------------------------------------------------------
!     check for problems
!-----------------------------------------------------------------------

#  if defined O_gent_mcwilliams && !defined O_isopycmix
        write (stdout,*)
     & '==> Error: "isopycmix" must be enabled with "gent_mcwilliams"'
        error = .true.
#  endif
#  if defined O_full_tensor
        write (stdout,*)
     & '==> Note: full isopycnal tensor is being used.               '
#  else
        write (stdout,*)
     & '==> Note: small angle approximation to isopycnal tensor is'
     &,'             being used.                                     '
#   if defined O_dm_taper
        write (stdout,*)
     & '==> Note: Danabasoglu-McWilliams hyperbolic tangent taper '
     &,'             is being used to reduce mixing coeff within     '
     &,'             areas of steeply sloping isopycnals             '
#   else
        write (stdout,*)
     & '==> Note: Re-scaling of Gerdes et al is being used to        '
     &,'             reduce mixing coeff within areas of steeply     '
     &,'             sloping isopycnals                              '
#   endif
#  endif
#  if !defined O_constvmix
        write (stdout,*)
     & '==> Error: "isopycmix" only works with "constvmix" for now   '
        error = .true.
#  endif

!-----------------------------------------------------------------------
!     print out tracer diffusion coefficients in isopycnal case
!-----------------------------------------------------------------------

      write(stdout,9102) ahisop, athkdf, ah, am

9102  format(/' ahisop = ',e12.5,' along isopyncal tracer mixing '
     &,  '(cm**2/sec) '/' athkdf = ',e12.5,' isopycnal thickness '
     &,'diffusion (cm**2/sec) '/' ah = ',e12.5,' cm**2/sec for '
     &,' am = ',e12.5,' cm**2/sec for horizontal mixing of momentum'/)
      ft = c1/(4.0*ahisop*dtts)
      delta_iso = dxt(1)*cst(jmt/2)*dzt(1)*ft
      i_delta = 1
      j_delta = 1
      k_delta = 1
      do jrow=2,jmt-1
        do i=2,imt-1
          do k=1,km
            delta1 = dxt(i)*cst(jrow)*dzt(k)*ft
            delta2 = dyt(jrow)*dzt(k)*ft
            if (delta_iso .ge. delta1 .or. delta_iso .ge. delta2) then
              i_delta = i
              j_delta = jrow
              k_delta = k
              delta_iso = min(delta1,delta2)
            endif
          enddo
        enddo
      enddo
      if (delta_iso .lt. p5) then
        s_minus = (c1 - sqrt(c1 - 4.0*delta_iso**2))/(c2*delta_iso)
        s_plus  = c1/s_minus
      else
        s_minus = c0
        s_plus  = c0
      endif

      write(stdout,'(a)')
     &'------------------------------------------------'
      write(stdout,'(a,e14.7)')
     &'The isopycnal diffusion grid factor delta_iso =',delta_iso
      write(stdout,'(a)')
     &'was determined at the grid point'
      write(stdout,'(a,i4,a,e14.7)')
     &'dxt(',i_delta,') = ',dxt(i_delta)
      write(stdout,'(a,i4,a,e14.7)')
     &'dyt(',j_delta,') = ',dyt(j_delta)
      write(stdout,'(a,i4,a,e14.7)')
     &'dzt(',k_delta,') = ',dzt(k_delta)
      write(stdout,'(a,i4,a,e14.7)')
     &'cst(',j_delta,') = ',cst(j_delta)
      write(stdout,'(a,e14.7)')
     &'dtts             = ',dtts
      write(stdout,'(a,e14.7)')
     &'ahisop           = ',ahisop

#  if defined O_full_tensor
      if (delta_iso .lt. p5) then
        write(stdout,'(/a,a,/2(a,e14.7)/)')
     & ' The full isopycnal tensor will be rescaled'
     &,' when the slope is in the range:'
     &,' s(-) = ',s_minus,', s_(+) = ',s_plus
      else
        write(stdout,'(/a,a/)') 'Isopycnal slopes in the full tensor'
     &,                  'will not be rescaled since delta_iso > 0.5'
      endif
#  else
      write(stdout,'(/a,e14.7/)')
     & 'Maximum allowable isopycnal slope is specified as slmx = ',slmx
#  endif

      write(stdout,'(a)')
     &'------------------------------------------------'

!-----------------------------------------------------------------------
!     store the square root of the tracer timestep acceleration values
!     into variable "dtxsqr" for use in isopycnal mixing
!-----------------------------------------------------------------------

      do k=1,km
       dtxsqr(k) = sqrt(dtxcel(k))
      enddo

      write (stdout,'(a)') ' Acceleration with depth "dtxsqr(k)"='
      write (stdout,'(5(1x,e12.5))') (dtxsqr(k),k=1,km)

!-----------------------------------------------------------------------
!     read ocean diffusion factor
!-----------------------------------------------------------------------

      fisop(:,:,:) = 1.
      fname = new_file_name ("O_diffac.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        ic(3) = km
        call openfile (fname, iou)
        exists = inqvardef('O_diffac', iou)
        if (exists) then
          write(*,*) "reading ocean diffusion factor from: ",fname
          call getvara ('O_diffac', iou, imtm2*jmtm2*km, ib, ic, tmpijk
     &,     c1, c0)
          fisop(2:imtm1,2:jmtm1,1:km) = tmpijk(1:imtm2,1:jmtm2,1:km)
#  if defined O_anisotropic_zonal_mixing
c         Andreas read this in since it is constant in time
!-----------------------------------------------------------------------
!     include anisotropic zonal equatorial isopycnal mixing scheme
!     see Getzlaff and Dietze, 2013, GRL for more information
!     smooth values linearly between 5-10 deg N/S to default values
!----------------------------------------------------------------------
          write(*,*) "reading equ zonal isopycnal diffusivity: ",fname
          call getvara ('O_ISOP', iou, imtm2*jmtm2*km, ib, ic, tmpijk
     &,     c1, c0)
          do j=1,jmtm2
             do i=1,imtm2
                do k=1,km
                   addisop(i+1,k,j+1) = tmpijk(i,j,k)*1e4
                enddo
             enddo
          enddo
#  endif
        endif
      endif

!-----------------------------------------------------------------------
!     initialize arrays
!-----------------------------------------------------------------------
      do j=1,jmw
        do k=1,km
          do i=1,imt
            alphai(i,k,j) = c0
            betai(i,k,j)  = c0
          enddo
        enddo
      enddo

      do j=jsmw,jemw
        do n=1,2
          do k=1,km
            do i=1,imt
              ddxt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo

      do j=1,jemw
        do n=1,2
          do k=1,km
            do i=1,imt
              ddyt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo
      do j=1,jmw
        do n=1,2
          do k=0,km
            do i=1,imt
              ddzt(i,k,j,n) = c0
            enddo
          enddo
        enddo
      enddo

      do j=jsmw,jemw
        do k=1,km
          do i=1,imt
            do kr=0,1
              do ip=0,1
                Ai_ez(i,k,j,ip,kr)  = c0
              enddo
            enddo
            do kr=0,1
              do ip=0,1
                Ai_bx(i,k,j,ip,kr) = c0
              enddo
            enddo
            do kr=0,1
              do jq=0,1
                Ai_by(i,k,j,jq,kr) = c0
              enddo
            enddo
            K33(i,k,j) = c0
            K11(i,k,j) = c0
#  if defined O_gent_mcwilliams
            adv_vetiso(i,k,j) = c0
#  endif
          enddo
        enddo
      enddo

      do j=1,jemw
        do k=1,km
          do i=1,imt
            do kr=0,1
              do jq=0,1
                Ai_nz(i,k,j,jq,kr)  = c0
              enddo
            enddo
            K22(i,k,j)        = c0
#  if defined O_gent_mcwilliams
            adv_vntiso(i,k,j) = c0
#  endif
          enddo
        enddo
      enddo

#  if defined O_gent_mcwilliams
      do j=jsmw,jemw
        do k=0,km
          do i=1,imt
            adv_vbtiso(i,k,j) = c0
            adv_fbiso(i,k,j)  = c0
          enddo
        enddo
      enddo

#  endif

      return
      end

      subroutine elements (joff, js, je, is, ie)

!=======================================================================
!     Estimate alpha, beta, and normal gradients on faces of T cells
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, is, ie, n, kp1, jrow, joff

      real dens, tq, sq, drodt, drods, tprime, sprime

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "grdvar.h"
      include "levind.h"
      include "mw.h"
      include "state.h"
      include "isopyc.h"
      include "dens.h"

!-----------------------------------------------------------------------
!     alpha and beta at centers of T cells
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,km
          do i=is,ie
            tprime = t(i,k,j,1,taum1)-to(k)
            sprime = t(i,k,j,2,taum1)-so(k)
            alphai(i,k,j) = drodt (tprime, sprime, k)
            betai(i,k,j)  = drods (tprime, sprime, k)
          enddo
        enddo
        call setbcx (alphai(1,1,j), imt, km)
        call setbcx (betai(1,1,j),  imt, km)
      enddo

!-----------------------------------------------------------------------
!     gradients at bottom face of T cells
!-----------------------------------------------------------------------

      do j=js,je
        do n=1,2
          do k=1,km
            kp1 = min(k+1,km)
            do i=is,ie
              ddzt(i,k,j,n) = tmask(i,kp1,j)*dzwr(k)*
     &                        (t(i,k,j,n,taum1) - t(i,kp1,j,n,taum1))
            enddo
          enddo
          do i=is,ie
            ddzt(i,0,j,n) = c0
          enddo
          call setbcx (ddzt(1,0,j,n), imt, km+1)
        enddo
      enddo

!-----------------------------------------------------------------------
!     gradients at eastern face of T cells
!-----------------------------------------------------------------------

#  if defined O_full_tensor
      do j=js,je
#  else
      do j=max(js-1,2), je-1
#  endif
        jrow = j + joff
        do n=1,2
          do k=1,km
            do i=is,ie
              ddxt(i,k,j,n) = tmask(i,k,j)*tmask(i+1,k,j)*cstr(jrow)*
     &                   dxur(i)*(t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
          call setbcx (ddxt(1,1,j,n), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     gradients at northern face of T cells
!-----------------------------------------------------------------------

      do j=max(js-1,1), je-1
        jrow = j + joff
        do n=1,2
          do k=1,km
            kp1 = min(k+1,km)
            do i=is,ie
              ddyt(i,k,j,n) = tmask(i,k,j)*tmask(i,k,j+1)*dyur(jrow)*
     &                 (t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
          call setbcx (ddyt(1,1,j,n), imt, km)
        enddo
      enddo

      return
      end

      subroutine isopyc (joff, js, je, is, ie)

!=======================================================================

!     Compute isopycnal diffusion coefficients Ai_ez, Ai_nz, Ai_bx,
!     Ai_by, and the Gent/McWilliams advection velocities.

!     input:
!       joff = offset relating row "j" in the MW to latitude "jrow"
!       js   = starting row within the MW for calculations
!       je   = ending row within the MW for calculations
!       is   = starting index longitude within the MW
!       ie   = ending index longitude within the MW

!     output:
!       Ai_ez = diffusion coefficient centered on east face of T cells
!       Ai_nz = diffusion coefficient centered on north face of T cells
!       Ai_bx = diffusion coefficient centered on bottom face of T cells
!       Ai_by = diffusion coefficient centered on bottom face of T cells

#  if defined O_gent_mcwilliams
!       adv_vetiso = isopycnal advective vel on east face of "T" cell
!       adv_vntiso = isopycnal advective vel on north face of "T" cell
!               (Note: this includes the cosine factor as in "adv_vnt")
!       adv_vbtiso = isopycnal advective vel on bottom face of "T" cell
#  endif

!=======================================================================

      implicit none

      integer istrt, is, iend, ie, joff, js, je

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"

!-----------------------------------------------------------------------
!     set local constants
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     estimate alpha, beta, and gradients on sides of T cells
!-----------------------------------------------------------------------

      call elements (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!     compute Ai_ez centered on eastern face of T cells
!     1 <= MW < last: calculate for rows jsmw..jemw
!     last MW       : calculate for rows jsmw..<=jemw
!-----------------------------------------------------------------------

      call ai_east (joff, max(js-1,2), je-1, is, ie)

!-----------------------------------------------------------------------
!     compute Ai_nz centered on the northern face of T cells
!     first MW     : calculate for rows 1..jemw
!     1 < MW < last: calculate for rows jsmw..jemw
!     last MW      : calculate for rows jsmw..<=jemw
!-----------------------------------------------------------------------

      call ai_north (joff, max(js-1,1), je-1, is, ie)

!-----------------------------------------------------------------------
!     evaluate Ai_bx & Ai_by centered on bottom face of T cells
!     1 <= MW < last: calculate for rows jsmw..jemw
!     last MW       : calculate for rows jsmw..<=jemw
!-----------------------------------------------------------------------

      call ai_bottom (joff, max(js-1,2), je-1, is, ie)

#  if defined O_gent_mcwilliams

!-----------------------------------------------------------------------
!     compute isopycnal advective velocities for tracers
!     first MW     : calculate (viso)           for rows 1..jemw
!                              (uiso,wiso)      for rows 2..jemw
!     1 < MW < last: calculate (uiso,viso,wiso) for rows jsmw..jemw
!     last MW      : calculate (uiso,viso,wiso) for rows jsmw..<=jemw
!-----------------------------------------------------------------------

      call isopyc_adv (joff, max(js-1,1), je-1, is, ie)

#  endif

      return
      end

      subroutine ai_east (joff, js, je, is, ie)

!=======================================================================
!     compute "Ai_ez" & "K11" at the center of eastern face of T cells.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, ie, is

      real sc, dzt4r, ai0, sumz, sxe, sumy, facty 

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "hmixc.h"
      include "isopyc.h"
      include "mw.h"

!-----------------------------------------------------------------------
!     compute "Ai_ez" on east face of T cell. Note re-scaling factor
!     which reduces mixing coefficient "Ai" where slope "sxe"
!     exceeds the critical slope "sc" for the small slope approx. For
!     the full tensor, it is re-scaled if outside the stable range.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,km
          sc = c1/(slmxr*dtxsqr(k))
          dzt4r = p5*dzt2r(k)
          do i=2,imtm1
            Ai0 = .5 * (fisop(i,jrow,k) + fisop(i+1,jrow,k)) * 
#  if defined O_KGMvar
     & .5 * (ahisop_var(i,k,jrow,niso) + ahisop_var(i+1,k,jrow,niso)) + 
#  else
     &   ahisop +
#  endif
     &   addisop(i,k,jrow) 

            
#  if defined O_full_tensor
            sumz = c0
            do kr=0,1
              do ip=0,1
                sxe = abs(drodxe(i,k,j,ip)/(drodze(i,k,j,ip,kr)+epsln))
                if (delta_iso .lt. p5 .and.
     &              (s_minus .lt. sxe .and. sxe .lt. s_plus)) then
                  Ai_ez(i,k,j,ip,kr) = Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &                                    *delta_iso*(sxe + c1/sxe)
                else
                  Ai_ez(i,k,j,ip,kr) = Ai0*tmask(i,k,j)*tmask(i+1,k,j)
                endif
                sumz = sumz + dzw(k-1+kr)*Ai_ez(i,k,j,ip,kr)
     &               *drodze(i,k,j,ip,kr)**2/(drodxe(i,k,j,ip)**2
     &          +  p5*(drodye(i,k,j,ip,0)**2 + drodye(i,k,j,ip,1)**2)
     &               + drodze(i,k,j,ip,kr)**2 + epsln)
              enddo
            enddo

            sumy = c0
            do jq=0,1
              facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do ip=0,1
                sumy = sumy + facty*Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &               *drodye(i,k,j,ip,jq)**2/(drodxe(i,k,j,ip)**2
     &               + drodye(i,k,j,ip,jq)**2 + epsln
     &          +  p5*(drodze(i,k,j,ip,0)**2 + drodze(i,k,j,ip,1)**2))
              enddo
            enddo
            K11(i,k,j) = dzt4r*sumz + p25*cstdytr(jrow)*sumy
#  else
            sumz = c0
            do kr=0,1
              do ip=0,1
                sxe = abs(drodxe(i,k,j,ip)/(drodze(i,k,j,ip,kr)+epsln))
#   if defined O_dm_taper
                Ai_ez(i,k,j,ip,kr) = Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &                               *p5*(c1-tanh((sxe-del_dm)*s_dmr))
#   else
                if (sxe .gt. sc) then
                  Ai_ez(i,k,j,ip,kr)  =  Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &                                   *(sc/(sxe + epsln))**2
                else
                  Ai_ez(i,k,j,ip,kr)  =  Ai0*tmask(i,k,j)*tmask(i+1,k,j)
                endif
#   endif
                sumz = sumz + dzw(k-1+kr)*Ai_ez(i,k,j,ip,kr)
              enddo
            enddo
            K11(i,k,j) = dzt4r*sumz
#  endif
          enddo
        enddo
        call setbcx (Ai_ez(1,1,j,0,0), imt, km)
        call setbcx (Ai_ez(1,1,j,1,0), imt, km)
        call setbcx (Ai_ez(1,1,j,0,1), imt, km)
        call setbcx (Ai_ez(1,1,j,1,1), imt, km)
        call setbcx (K11(1,1,j), imt, km)
      enddo

      return
      end

      subroutine ai_north (joff, js, je, is, ie)

!=======================================================================
!     compute "Ai_nz" & "K22" at the center of northern face of T cells.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, ie, is

      real sc, dzt4r, ai0, sumz, syn, sumx

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "hmixc.h"
      include "isopyc.h"
      include "mw.h"

!-----------------------------------------------------------------------
!     compute "Ai_nz" on north face of T cell. Note re-scaling factor
!     which reduces mixing coefficient "Ai" where slope "syn"
!     exceeds the critical slope "sc" for the small slope approx. For
!     the full tensor, it is re-scaled if outside the stable range.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,km
          sc = c1/(slmxr*dtxsqr(k))
          dzt4r = p5*dzt2r(k)
          do i=2,imtm1
            Ai0 = p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k)) * 
#  if defined O_KGMvar
     &            p5*(ahisop_var(i,k,jrow,niso) + 
     &                ahisop_var(i,k,jrow+1,niso))
#  else
     &            ahisop
#  endif

#  if defined O_full_tensor
            sumz = c0
            do kr=0,1
              do jq=0,1
                syn=abs(drodyn(i,k,j,jq)/(drodzn(i,k,j,jq,kr)+epsln))
                if (delta_iso .lt. p5 .and.
     &            (s_minus .lt. syn .and. syn .lt. s_plus)) then
                  Ai_nz(i,k,j,jq,kr) = Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &                                    *delta_iso*(syn + c1/syn)
                else
                  Ai_nz(i,k,j,jq,kr) = Ai0*tmask(i,k,j)*tmask(i,k,j+1)
                endif
                sumz = sumz + dzw(k-1+kr)*Ai_nz(i,k,j,jq,kr)
     &               *drodzn(i,k,j,jq,kr)**2/(
     &          +  p5*(drodxn(i,k,j,0,jq)**2 + drodxn(i,k,j,1,jq)**2)
     &          + drodyn(i,k,j,jq)**2 + drodzn(i,k,j,jq,kr)**2 + epsln)
              enddo
            enddo

            sumx = c0
            do ip=0,1
              do jq=0,1
                  sumx =sumx+dxu(i-1+ip)*Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &               *drodxn(i,k,j,ip,jq)**2/(drodxn(i,k,j,ip,jq)**2
     &               + drodyn(i,k,j,jq)**2 + epsln
     &       +    p5*(drodzn(i,k,j,jq,0)**2 + drodzn(i,k,j,jq,1)**2))
              enddo
            enddo
            K22(i,k,j) = dzt4r*sumz + p25*dxtr(i)*sumx
#  else
            sumz = c0
            do kr=0,1
              do jq=0,1
                syn =abs(drodyn(i,k,j,jq)/(drodzn(i,k,j,jq,kr)+epsln))
#   if defined O_dm_taper
                Ai_nz(i,k,j,jq,kr) = Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &                               *p5*(c1-tanh((syn-del_dm)*s_dmr))
#   else
                if (syn .gt. sc) then
                  Ai_nz(i,k,j,jq,kr) = Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &                                   *(sc/(syn + epsln))**2
                else
                  Ai_nz(i,k,j,jq,kr) = Ai0*tmask(i,k,j)*tmask(i,k,j+1)
                endif
#   endif
                sumz = sumz + dzw(k-1+kr)*Ai_nz(i,k,j,jq,kr)
              enddo
            enddo
            K22(i,k,j) = dzt4r*sumz
#  endif
          enddo
        enddo
        call setbcx (Ai_nz(1,1,j,0,0), imt, km)
        call setbcx (Ai_nz(1,1,j,1,0), imt, km)
        call setbcx (Ai_nz(1,1,j,0,1), imt, km)
        call setbcx (Ai_nz(1,1,j,1,1), imt, km)
        call setbcx (K22(1,1,j), imt, km)
      enddo

      return
      end

      subroutine ai_bottom (joff, js, je, is, ie)

!=======================================================================
!     compute Ai_bx, Ai_by, and K33 at the center of the bottom face of
!     T cells.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, ie, is

      real sc, kp1, ai0, sumx, sxb, sumy, facty, syb

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "grdvar.h"
      include "hmixc.h"
      include "isopyc.h"
      include "mw.h"

!-----------------------------------------------------------------------
!     compute "Ai_bx", "Ai_by", & K33 on bottom face of T cell. Note
!     re-scaling factor to reduce mixing coefficient "Ai" where slopes
!     "sxb" and "syb" exceeds the critical slope "sc" for the small slope approx. For
!     the full tensor, it is re-scaled if outside the stable range.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,kmm1
          sc = c1/(slmxr*dtxsqr(k))
          kp1   = min(k+1,km)
          do i=2,imtm1
            Ai0 = p5*(fisop(i,jrow,k+1) + fisop(i,jrow,k)) * 
#  if defined O_KGMvar
     &        p5*(ahisop_var(i,k,jrow,niso)+ahisop_var(i,k+1,jrow,niso))
#  else
     &        ahisop
#  endif
#  if defined O_full_tensor

!           eastward slopes at the base of T cells
            sumx = c0
            do ip=0,1
              do kr=0,1
                sxb = abs(drodxb(i,k,j,ip,kr)/(drodzb(i,k,j,kr)+epsln))
                if (delta_iso .lt. p5 .and.
     &              (s_minus .lt. sxb .and. sxb .lt. s_plus)) then
                  Ai_bx(i,k,j,ip,kr)  =  Ai0*tmask(i,k+1,j)
     &                                   *delta_iso*(sxb + c1/sxb)
                else
                  Ai_bx(i,k,j,ip,kr)  =  Ai0*tmask(i,k+1,j)
                endif
                sumx = sumx + dxu(i-1+ip)*Ai_bx(i,k,j,ip,kr)
     &            *drodxb(i,k,j,ip,kr)**2/(drodxb(i,k,j,ip,kr)**2
     &       +  p5*(drodyb(i,k,j,0,kr)**2 + drodyb(i,k,j,1,kr)**2)
     &               + drodzb(i,k,j,kr)**2 + epsln)
              enddo
            enddo

!           northward slopes at the base of T cells

            sumy = c0
            jrow = j + joff
            do jq=0,1
              facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do kr=0,1
                syb = abs(drodyb(i,k,j,jq,kr)/(drodzb(i,k,j,kr)+epsln))
                if (delta_iso .lt. p5 .and.
     &              (s_minus .lt. syb .and. syb .lt. s_plus)) then
                  Ai_by(i,k,j,jq,kr)  =  Ai0*tmask(i,k+1,j)
     &                                   *delta_iso*(syb + c1/syb)
                else
                  Ai_by(i,k,j,jq,kr)  =  Ai0*tmask(i,k+1,j)
                endif
                sumy = sumy + facty*Ai_by(i,k,j,jq,kr)
     &            *drodyb(i,k,j,jq,kr)**2/(
     &           p5*(drodxb(i,k,j,0,kr)**2 + drodxb(i,k,j,1,kr)**2)
     &            + drodyb(i,k,j,jq,kr)**2
     &            + drodzb(i,k,j,kr)**2 + epsln)
              enddo
            enddo

            K33(i,k,j) = dxt4r(i)*sumx + dyt4r(jrow)*cstr(jrow)*sumy
#  else

!           eastward slopes at the base of T cells

            sumx = c0
            do ip=0,1
              do kr=0,1
                sxb = abs(drodxb(i,k,j,ip,kr)/(drodzb(i,k,j,kr)+epsln))
#   if defined O_dm_taper
                Ai_bx(i,k,j,ip,kr) = Ai0*tmask(i,k+1,j)
     &                               *p5*(c1-tanh((sxb-del_dm)*s_dmr))
#   else
                if (sxb .gt. sc) then
                  Ai_bx(i,k,j,ip,kr)  =  Ai0*tmask(i,k+1,j)
     &                                  *(sc/(sxb + epsln))**2
                else
                  Ai_bx(i,k,j,ip,kr)  =  Ai0*tmask(i,k+1,j)
                endif
#   endif
                sumx = sumx + dxu(i-1+ip)*Ai_bx(i,k,j,ip,kr)*sxb**2
              enddo
            enddo

!           northward slopes at the base of T cells

            sumy = c0
            do jq=0,1
              facty = csu(jrow-1+jq)*dyu(jrow-1+jq)
              do kr=0,1
                syb = abs(drodyb(i,k,j,jq,kr)/(drodzb(i,k,j,kr)+epsln))
#   if defined O_dm_taper
                Ai_by(i,k,j,jq,kr) = Ai0*tmask(i,k+1,j)
     &                               *p5*(c1-tanh((syb-del_dm)*s_dmr))
#   else
                if (syb .gt. sc) then
                  Ai_by(i,k,j,jq,kr)  =  Ai0*tmask(i,k+1,j)
     &                                  *(sc/(syb + epsln))**2
                else
                  Ai_by(i,k,j,jq,kr)  =  Ai0*tmask(i,k+1,j)
                endif
#   endif
                sumy = sumy + facty*Ai_by(i,k,j,jq,kr)*syb**2
              enddo
            enddo

            K33(i,k,j) = dxt4r(i)*sumx + dyt4r(jrow)*cstr(jrow)*sumy
#  endif
          enddo
        enddo
        call setbcx (Ai_bx(1,1,j,1,0), imt, km)
        call setbcx (Ai_bx(1,1,j,0,0), imt, km)
        call setbcx (Ai_bx(1,1,j,1,1), imt, km)
        call setbcx (Ai_bx(1,1,j,0,1), imt, km)
        call setbcx (Ai_by(1,1,j,1,0), imt, km)
        call setbcx (Ai_by(1,1,j,0,0), imt, km)
        call setbcx (Ai_by(1,1,j,1,1), imt, km)
        call setbcx (Ai_by(1,1,j,0,1), imt, km)
        call setbcx (K33(1,1,j), imt, km)
      enddo

      return
      end

      subroutine isoflux (joff, js, je, is, ie, n)

!=======================================================================
!     isopycnal diffusive tracer fluxes are computed.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, km1kr, kpkr
      integer n, ie, is

      real dzt4r, sumz, flux_x, ai0, sumy, facty, csu_dzt4r, flux_y
      real sumx, ck

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "grdvar.h"
      include "hmixc.h"
      include "isopyc.h"
      include "levind.h"
      include "mw.h"
      include "vmixc.h"

!-----------------------------------------------------------------------
!     construct total isopycnal tracer flux at east face of "T" cells
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,km
          dzt4r = p5*dzt2r(k)
          do i=2,imtm1
            sumz = c0
            do kr=0,1
              km1kr = max(k-1+kr,1)
              kpkr = min(k+kr,km)
              do ip=0,1
                sumz = sumz - Ai_ez(i,k,j,ip,kr)
     &                *(t(i+ip,km1kr,j,n,taum1)-t(i+ip,kpkr,j,n,taum1))
     &                *drodxe(i,k,j,ip)
#  if defined O_full_tensor
     &                *drodze(i,k,j,ip,kr)/(drodxe(i,k,j,ip)**2
     &          +  p5*(drodye(i,k,j,ip,0)**2 + drodye(i,k,j,ip,1)**2)
     &               + drodze(i,k,j,ip,kr)**2 + epsln)
#  else
     &               /(drodze(i,k,j,ip,kr) + epsln)
#  endif
              enddo
            enddo
            flux_x = dzt4r*sumz
#  if defined O_full_tensor
#   if defined O_anisotropic_zonal_mixing 
            ck = addisop(i,k,jrow)
            print*, "ADDISO", ck
#   endif
            Ai0 = (ahisop+addisop(i,k,jrow))*p5*(fisop(i,jrow,k)
     &           + fisop(i+1,jrow,k))
            sumy = c0
            do jq=0,1
              facty = csu(jrow-1+jq)
              do ip=0,1
                sumy = sumy - facty*Ai0*tmask(i,k,j)*tmask(i+1,k,j)
     &       *(t(i+ip,k,j+jq,n,taum1)-t(i+ip,k,j-1+jq,n,taum1))
     &       *drodye(i,k,j,ip,jq)*drodxe(i,k,j,ip)/(drodxe(i,k,j,ip)**2
     &               + drodye(i,k,j,ip,jq)**2 + epsln
     &          +  p5*(drodze(i,k,j,ip,0)**2 + drodze(i,k,j,ip,1)**2))
              enddo
            enddo
            flux_x = flux_x + p25*cstdytr(jrow)*sumy
#  endif
            diff_fe(i,k,j) = diff_fe(i,k,j) + K11(i,k,j)*cstdxur(i,j)*
     &                       (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
     &                     + flux_x
          enddo
        enddo
        call setbcx (diff_fe(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     construct total isopycnal tracer flux at north face of "T" cells
!-----------------------------------------------------------------------

      do j=js-1,je
        jrow = j + joff
        do k=1,km
          csu_dzt4r = csu(jrow)*p5*dzt2r(k)
          do i=2,imtm1
            sumz = c0
            do kr=0,1
              km1kr = max(k-1+kr,1)
              kpkr = min(k+kr,km)
              do jq=0,1
                sumz = sumz - Ai_nz(i,k,j,jq,kr)
     &                 *(t(i,km1kr,j+jq,n,taum1)-t(i,kpkr,j+jq,n,taum1))
     &                 *drodyn(i,k,j,jq)
#  if defined O_full_tensor
     &                 *drodzn(i,k,j,jq,kr)/(
     &          +  p5*(drodxn(i,k,j,0,jq)**2 + drodxn(i,k,j,1,jq)**2)
     &          + drodyn(i,k,j,jq)**2 + drodzn(i,k,j,jq,kr)**2 + epsln)
#  else
     &               /(drodzn(i,k,j,jq,kr) + epsln)
#  endif
              enddo
            enddo
            flux_y = csu_dzt4r*sumz
#  if defined O_full_tensor
            Ai0 = (ahisop+addisop(i,k,jrow))*p5*(fisop(i,jrow,k)
     &           + fisop(i,jrow+1,k))
            sumx = c0
            do ip=0,1
              do jq=0,1
                sumx = sumx - Ai0*tmask(i,k,j)*tmask(i,k,j+1)
     &               *(t(i+ip,k,j+jq,n,taum1)-t(i-1+ip,k,j+jq,n,taum1))
     &                 *cstr(jrow+jq)
     &               *drodxn(i,k,j,ip,jq)*drodyn(i,k,j,jq)/(
     &             drodxn(i,k,j,ip,jq)**2 + drodyn(i,k,j,jq)**2 + epsln
     &          +  p5*(drodzn(i,k,j,jq,0)**2 + drodzn(i,k,j,jq,1)**2))
              enddo
            enddo
            flux_y = flux_y + p25*csu(jrow)*dxtr(i)*sumx
#  endif
            diff_fn(i,k,j) = diff_fn(i,k,j)+K22(i,k,j)*csu_dyur(jrow)*
     &                       (t(i,k,j+1,n,taum1)-t(i,k,j,n,taum1))
     &                    + flux_y
          enddo
        enddo
        call setbcx (diff_fn(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     compute the vertical tracer flux "diff_fbiso" containing the K31
!     and K32 components which are to be solved explicitly. The K33
!     component will be treated implicitly. Note that there are some
!     cancellations of dxu(i-1+ip) and dyu(jrow-1+jq)
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,kmm1
          do i=2,imtm1
            sumx = c0
            do ip=0,1
              do kr=0,1
                sumx = sumx
     &            - Ai_bx(i,k,j,ip,kr)*cstr(jrow)*
     &             (t(i+ip,k+kr,j,n,taum1) - t(i-1+ip,k+kr,j,n,taum1))
     &             *drodxb(i,k,j,ip,kr)
#  if defined O_full_tensor
     &             *drodzb(i,k,j,kr)/(drodxb(i,k,j,ip,kr)**2
     &       +  p5*(drodyb(i,k,j,0,kr)**2 + drodyb(i,k,j,1,kr)**2)
     &            + drodzb(i,k,j,kr)**2 + epsln)
#  else
     &              /(drodzb(i,k,j,kr) + epsln)
#  endif
              enddo
            enddo
            sumy = c0
            do jq=0,1
              do kr=0,1
                sumy = sumy
     &           - Ai_by(i,k,j,jq,kr)*csu(jrow-1+jq)*
     &            (t(i,k+kr,j+jq,n,taum1)-t(i,k+kr,j-1+jq,n,taum1))
     &               *drodyb(i,k,j,jq,kr)
#  if defined O_full_tensor
     &               *drodzb(i,k,j,kr)/(
     &           p5*(drodxb(i,k,j,0,kr)**2 + drodxb(i,k,j,0,kr)**2)
     &            + drodyb(i,k,j,jq,kr)**2
     &            + drodzb(i,k,j,kr)**2 + epsln)
#  else
     &              /(drodzb(i,k,j,kr) + epsln)
#  endif
              enddo
            enddo
            diff_fbiso(i,k,j) = dxt4r(i)*sumx
     &                          + dyt4r(jrow)*cstr(jrow)*sumy
          enddo
        enddo
        do i=2,imtm1
          diff_fbiso(i,0,j)  = c0
          diff_fbiso(i,km,j) = c0
        enddo
        call setbcx (diff_fbiso(1,0,j), imt, km+1)
      enddo

#  if defined O_gent_mcwilliams

!-----------------------------------------------------------------------
!     compute advective tracer flux at the center of the bottom face of
!     the "T" cells
!-----------------------------------------------------------------------

      do j=js,je
        do k=1,kmm1
          do i=2,imt-1
            adv_fbiso(i,k,j) = adv_vbtiso(i,k,j)*
     &                         (t(i,k,j,n,taum1) + t(i,k+1,j,n,taum1))
          enddo
        enddo
      enddo

!     now consider the top and bottom boundaries

      do j=js,je
        do i=2,imt-1
          adv_fbiso(i,0,j)  = c0
          adv_fbiso(i,km,j) = c0
        enddo
      enddo
#  endif

      return
      end

#  if defined O_gent_mcwilliams
      subroutine isopyc_adv (joff, js, je, is, ie)

!=======================================================================
!     compute isopycnal transport velocities.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, km1, kp1
      integer jstrt, ie, is

      real sc, ath0, at, bt, stn, ab, bb, sbn, absstn, abssbn
      real ath_t, ath_b, ste, sbe, absste, abssbe

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "isopyc.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
#   if defined O_time_averages
      include "timeavgs.h"
#   endif

      real top_bc(km), bot_bc(km)

#   if defined O_KGMvar
      real sigma, max_cori
      Lm         = 50.e5
      eddy_min   = 3.0e6
      eddy_max   = 3.0e7
      coef       = 980./1.025
      pii        = 4.0 * atan(1.0)
#   endif

      do k=1,km
        top_bc(k) = c1
        bot_bc(k) = c1
      enddo
      top_bc(1)  = c0
      bot_bc(km) = c0

      do j=js,je
        jrow = j + joff
        do i=1,imtm1
          at =    p5 * (alphai(i,1,j) + alphai(i,1,j+1))
          bt =    p5 * (betai(i,1,j) + betai(i,1,j+1))
          drodytn(i,1,jrow) = at * ddyt(i,1,j,1)
     &                      + bt * ddyt(i,1,j,2)
          drodztn(i,1,jrow) = at * (ddzt(i,1,j,1) + ddzt(i,1,j+1,1))*p5
     &                      + bt * (ddzt(i,1,j,2) + ddzt(i,1,j+1,2))*p5
!         -- zonal --
          at =     p5 * (alphai(i,1,j) + alphai(i+1,1,j))
          bt =     p5 * (betai(i,1,j) + betai(i+1,1,j))
          drodxte(i,1,jrow) = at * ddxt(i,1,j,1)
     &                      + bt * ddxt(i,1,j,2) 
          drodzte(i,1,jrow) = at * (ddzt(i,1,j,1) + ddzt(i+1,1,j,1))*p5
     &                      + bt * (ddzt(i,1,j,2) + ddzt(i+1,1,j,2))*p5
          do k=1,km
            km1 = max(k-1,1)
            kp1 = min(k+1,km)
!     compute the meridional slopes of isopycnals at the north and top
!     or bottom faces of the "T" grid
            ab =     (alphai(i,k,j) + alphai(i,k,j+1) + alphai(i,kp1,j)
     &              + alphai(i,kp1,j+1)) * 0.25
            bb =     (betai(i,k,j) + betai(i,k,j+1) + betai(i,kp1,j)
     &              + betai(i,kp1,j+1)) * 0.25
            drodybn(i,k,jrow) = ab*p5*(ddyt(i,k,j,1) + ddyt(i,kp1,j,1))
     &                        + bb*p5*(ddyt(i,k,j,2) + ddyt(i,kp1,j,2))
            drodzbn(i,k,jrow) = ab*p5*(ddzt(i,k,j,1) + ddzt(i,k,j+1,1))
     &                        + bb*p5*(ddzt(i,k,j,2) + ddzt(i,k,j+1,2))

            ! Copy data to top face of present cell from bottom faces of
            ! cell above
            if (k.gt.1) then
              drodytn(i,k,jrow) = drodybn(i,km1,jrow)
              drodztn(i,k,jrow) = drodzbn(i,km1,jrow)
            endif
            
!     compute the zonal slopes of isopycnals at the east and top
!     or bottom faces of the "T" grid
            ab =     (alphai(i,k,j) + alphai(i+1,k,j) + alphai(i,kp1,j)
     &              + alphai(i+1,kp1,j)) * 0.25
            bb =     (betai(i,k,j) + betai(i+1,k,j) + betai(i,kp1,j)
     &              + betai(i+1,kp1,j)) * 0.25
            drodxbe(i,k,jrow) =  ab*p5*(ddxt(i,k,j,1) + ddxt(i,kp1,j,1))
     &                         + bb*p5*(ddxt(i,k,j,2) + ddxt(i,kp1,j,2))
            drodzbe(i,k,jrow) =  ab*p5*(ddzt(i,k,j,1) + ddzt(i+1,k,j,1))
     &                         + bb*p5*(ddzt(i,k,j,2) + ddzt(i+1,k,j,2))
     
            if (k.gt.1) then
              drodxte(i,k,jrow) = drodxbe(i,km1,jrow)
              drodzte(i,k,jrow) = drodzbe(i,km1,jrow)
            endif
            
          enddo  ! k
        enddo    ! i
      enddo      ! j

#   if defined O_KGMvar
!-----------------------------------------------------------------------
! Compute the 3 dimensional isopycnal diffiusivity coefficient K_GM
! as outlined in Oleg, Eden, et. al. See documentation for more detail.
! *ToDo: add more
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do i=2,imtm1
          Lr(i,jrow)    = 0. 
          LRhi(i,:,jrow)  = 0.
          sigma_ave     = 0.
          stratif_int   = 0.
          sum_zz        = 0.
          beta          = (cori(i,j+1,1) - cori(i,j-1,1)) /
     &                    (2.0 * dyt(j))

!         Only integrate baroclinicty below the mixed layer (it will be 
!         huge in the mixed layer, giving huge K_gm in shallow seas).
!         Also, can integrate downwards to about 2000m, not further, to
!         emphasize the mid-depth water where we have more faith in the
!         stratification. (To emphasize the first baroclinic mode??)
!          do k=5,min(kmt(i,jrow)-1,20)
          do k=1,kmt(i,jrow)  ! change hard value to a calculated
           ! stratif_int = vertical integral of stratification, on T
           !   cells
            stratif_int = stratif_int + dzw(k) *
     &                    sqrt(coef * p5 * abs(drodzb(i,k,jrow,0)
     &                                       + drodzb(i,k,jrow,1)))
          enddo
          do k=1,kmt(i,jrow)
            kgm(i,k,jrow,niso) = eddy_min
            ahisop_var(i,k,jrow,niso) = eddy_min

            km1 = max(k-1,1)
            countx = 2
            county = 2
            if (drodxte(i  ,k,jrow) .eq. 0) countx = countx - 1
            if (drodxte(i-1,k,jrow) .eq. 0) countx = countx - 1
            if (drodytn(i,k,jrow  ) .eq. 0) county = county - 1
            if (drodytn(i,k,jrow-1) .eq. 0) county = county - 1
            
            grd_rho_x = drodxte(i,k,jrow) + drodxte(i-1,k,jrow)
            grd_rho_y = drodytn(i,k,jrow) + drodytn(i,k,jrow-1)
            if (countx.eq.0 .and. county.eq.0) then
              abs_grd_rho2 = 0.
            elseif (countx.eq.0) then
              abs_grd_rho2 = (grd_rho_y / county) ** 2
            elseif (county.eq.0) then
              abs_grd_rho2 = (grd_rho_x / countx) ** 2
            else
              abs_grd_rho2 = (grd_rho_x / countx) ** 2
     &                     + (grd_rho_y / county) ** 2
            endif
            
            ! Vertical gradient of density at the top face of "T" cells.
            abs_drho_dz = p5 * abs(drodzb(i,km1,jrow,0) +
     &                             drodzb(i,km1,jrow,1))
            
            ! Limit N2 > 1.e-8 / coef
            abs_drho_dz = max(1.e-8/coef, abs_drho_dz)

            max_cori = max(abs(cori(i,j,1)),
     &                 (2.*beta*stratif_int/pii)**(1/2))             

            sigma     = sqrt(abs_grd_rho2) /
     &                       sqrt((abs_drho_dz/coef) +
     &                            (gamma_eden * abs_grd_rho2 / 
     &                             max_cori**2))
           
            sigma_ave     = sigma_ave + dzw(km1) * sigma
            sum_zz = sum_zz + dzw(km1)
           
c            kgm(i,k,jrow,niso) = sigma

c            ! Local Rhine's Scale
c            LRhi(i,k,jrow) = sigma / beta
            LRhi(i,k,jrow) = 50.e5 ! 50 km

          enddo ! k

          ! Now get the average
          if (sum_zz.ne.0.) then
            sigma_ave = sigma_ave / sum_zz
          endif
          kgm(i,:,jrow,niso) = sigma_ave

          ! Rhine's Scale. Averaged over depth
c          LRhi(i,jrow) = sigma_ave / beta

          ! Retrieve a minimum of Rossby Radius. See Eden 2009.
          Lr1 = stratif_int
     &                   / (pii * p5 * abs(cori(i,j,1) + cori(i,j-1,1)))
          Lr2 = sqrt(stratif_int / (2. * pii * beta))
          Lr(i,jrow) = min(Lr1, Lr2)
c          Lr(i,jrow) = 50.e5
        enddo  ! i
      enddo ! j=js,je
      
      ! Finally, derive Kgm using minimum of Rhine's scale and Rossby
      ! Radius as described in Eden 2009.
      do j=js,je
        jrow = j + joff
        do i=2,imtm1
c            if (Lr(i,jrow).lt.LRhi(i,jrow)) then
c              kgm(i,:,jrow,niso) = kgm(i,:,jrow,niso) * Lr(i,jrow)**2
c            else
c              kgm(i,:,jrow,niso) = kgm(i,:,jrow,niso) * LRhi(i,jrow)**2
c            endif
            kgm(i,:,jrow,niso) = c_eden*kgm(i,:,jrow,niso)
     &                           * Lr(i,jrow)*LRhi(i,:,jrow)
c     &                           * min(Lr(i,jrow),LRhi(i,:,jrow))**2
c     &                           * Lr(i,jrow)**2
            kgm(i,:,jrow,niso) = min(max(kgm(i,:,jrow,niso),eddy_min), 
     &                                                      eddy_max)
        enddo
      enddo 
      
c      ahisop_var(:,:,:,niso) = kgm(:,:,:,niso) * ( 3.0 / 2.0 )
      ahisop_var(:,:,:,niso) = ahisop
c      kgm(:,:,:,niso) = athkdf

c      print *, 'ahiso', ahisop_var(50,2,50,niso)

      do j=1,niso  ! steal index j
        call setbcx (kgm(1,1,1,j), imt, jmt)
        call setbcx (ahisop_var(1,1,1,j), imt, jmt)
      enddo
      call setbcx (Lr(1,1), imt, jmt)
      call setbcx (LRhi(1,1,1), imt, jmt)
#   endif  ! O_KGMvar

!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "T" grid cell.
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do k=1,km

          sc = c1/(slmxr*dtxsqr(k))
          km1 = max(k-1,1)
          kp1 = min(k+1,km)
          do i=1,imt

#   if defined O_KGMvar
            ! Select y component of kgm if it exists (if O_KGM_aniso)
            Ath0 = p5 * (kgm(i,k,jrow,niso) + kgm(i,k,jrow+1,niso)) * 
     &             p5 * (fisop(i,jrow,k) + fisop(i,jrow+1,k))          
c            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k)) 
#   else
            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k)) 
#   endif
c            print *,'Ath0',i,k,jrow, Ath0
            ! 0.125*epsln is done to match the original stn etc. code.
            ! (I divided my drod*'s before calculating slopes stn etc.)
            stn = - drodytn(i,k,jrow) / (drodztn(i,k,jrow) +0.125*epsln)
            sbn = - drodybn(i,k,jrow) / (drodzbn(i,k,jrow) +0.125*epsln)
            absstn = abs(stn)
            abssbn = abs(sbn)
#   if defined O_dm_taper
            ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
     &             *p5*(c1-tanh((absstn-del_dm)*s_dmr))
            ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
     &             *p5*(c1-tanh((abssbn-del_dm)*s_dmr))
#   else
            ! Gerdes et al 1991 taper. cf. MITGCM manual sec 6.4.1.5
            if (absstn .gt. sc) then
              ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
     &              *(sc/(absstn + epsln))**2
            else
              ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
            endif
            if (abssbn .gt. sc) then
              ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
     &              *(sc/(abssbn + epsln))**2
            else
              ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
            endif
#   endif
            adv_vntiso(i,k,j) = -(ath_t*stn*top_bc(k) -
     &                            ath_b*sbn*bot_bc(k))*dztr(k)*csu(jrow)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     compute the zonal component of the isopycnal mixing velocity
!     at the center of the eastern face of the "T" grid cell.
!-----------------------------------------------------------------------
      jstrt = max(js,jsmw)
      do j=jstrt,je
        jrow = j + joff
        do k=1,km
          sc = c1/(slmxr*dtxsqr(k))
          km1 = max(k-1,1)
          kp1 = min(k+1,km)
          do i=1,imtm1

#   if defined O_KGMvar

            Ath0 = p5 * (kgm(i,k,jrow,niso) + kgm(i+1,k,jrow,niso)) *
     &             p5 * (fisop(i,jrow,k) + fisop(i+1,jrow,k))
c            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k)) 
#   else
            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k)) 
#   endif
            ! 0.125*epsln is done to match the original stn etc. code.
            ! (I divided my drod*'s before calculating slopes stn etc.)
            ste = - drodxte(i,k,jrow) / (drodzte(i,k,jrow) +0.125*epsln)
            sbe = - drodxbe(i,k,jrow) / (drodzbe(i,k,jrow) +0.125*epsln)
            absste = abs(ste)
            abssbe = abs(sbe)
#   if defined O_dm_taper
            ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
     &             *p5*(c1-tanh((absste-del_dm)*s_dmr))
            ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
     &             *p5*(c1-tanh((abssbe-del_dm)*s_dmr))
#   else   
            ! Gerdes et al 1991 taper. cf. MITGCM manual sec 6.4.1.5
            if (absste .gt. sc) then
              ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
     &              *(sc/(absste + epsln))**2
            else
              ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
            endif
            if (abssbe .gt. sc) then
              ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
     &              *(sc/(abssbe + epsln))**2
            else
              ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
            endif
#   endif
            adv_vetiso(i,k,j) = -(ath_t*ste*top_bc(k) -
     &                            ath_b*sbe*bot_bc(k))*dztr(k)
          enddo
        enddo
      enddo

!     set the boundary conditions

      do j=jstrt,je
        call setbcx (adv_vetiso(1,1,j), imt, km)
      enddo

!-----------------------------------------------------------------------
!     compute the vertical component of the isopycnal mixing velocity
!     at the center of the bottom face of the "T" cells, using the
!     continuity equation for the isopycnal mixing velocities
!-----------------------------------------------------------------------

      do j=jstrt,je
        do i=1,imt
          adv_vbtiso(i,0,j) = c0
        enddo
      enddo

      do j=jstrt,je
        jrow = j + joff
        do k=1,kmm1
          do i=2,imt
            adv_vbtiso(i,k,j) = dzt(k)*cstr(jrow)*(
     &      (adv_vetiso(i,k,j) - adv_vetiso(i-1,k,j))*dxtr(i) +
     &      (adv_vntiso(i,k,j) - adv_vntiso(i,k,j-1))*dytr(jrow))
          enddo
        enddo
      enddo

      do j=jstrt,je
        do k=1,kmm1
          do i=2,imt
            adv_vbtiso(i,k,j) = adv_vbtiso(i,k,j) + adv_vbtiso(i,k-1,j)
          enddo
        enddo
      enddo

      do j=jstrt,je
        jrow = j + joff
        do i=2,imt
          adv_vbtiso(i,kmt(i,jrow),j) = c0
        enddo
      enddo

!     set the boundary conditions

      do j=jstrt,je
        call setbcx (adv_vbtiso(1,0,j), imt, km+1)
      enddo

#   if defined O_time_averages

!-----------------------------------------------------------------------
!     accumulate time average gm velocities
!-----------------------------------------------------------------------

      if (timavgperts .and. .not. euler2) then
        do j=jstrt,je
          jrow = j + joff
          do k=1,kmm1
            do i=2,imt
              ta_vetiso(i,k,jrow) = ta_vetiso(i,k,jrow) +
     &                              adv_vetiso(i,k,j)
              ta_vntiso(i,k,jrow) = ta_vntiso(i,k,jrow) +
     &                              adv_vntiso(i,k,j)
              ta_vbtiso(i,k,jrow) = ta_vbtiso(i,k,jrow) +
     &                              adv_vbtiso(i,k,j)
            enddo
          enddo
#    if defined O_KGMvar
          do k=1,km
             do i=1,imt
                ta_LRhi(i,k,jrow) = ta_LRhi(i,k,jrow) + LRhi(i,k,jrow)
                ta_kgm(i,k,jrow,niso) = ta_kgm(i,k,jrow,niso) + 
     &                                  kgm(i,k,jrow,niso)
                ta_isop(i,k,jrow,niso) = ta_isop(i,k,jrow,niso) + 
     &                                   ahisop_var(i,k,jrow,niso)
             enddo
          enddo
          do i=1, imt
             ta_Lr(i,jrow) = ta_Lr(i,jrow) + Lr(i,jrow)
          enddo
#    endif
        enddo
      endif
#   endif     
#  endif
# endif
#endif
  
      return
      end

