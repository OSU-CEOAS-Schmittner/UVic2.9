      subroutine isopi (error, am, ah)

#if defined O_mom
# if defined O_isopycmix
!=======================================================================

!               Initialization for isopycnal mixing scheme

!      -Disopycmix gives either the full or small angle isopycnal
!                  mixing tensor

!      -Disopycmix -Dgent_mcwilliams gives the isopycnal mixing tensor
!       plus the advective velocity parameterization of Gent_McWilliams

!     input:
!       error  = logical to signal problems
!       am     = background mixing coeff for momentum
!       ah     = horizontal mixing coeff for tracers
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

      logical exists, inqvardef, error

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
#  if !defined O_gent_mcwilliams
      real athkdf
#  endif

      namelist /isopyc/ slmx, ahisop, athkdf, del_dm, s_dm

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

      slmx  = 1.0/100.0

!     isopycnal diffusion coefficient

      ahisop = 1.e7

!     isopycnal thickness diffusion coefficient

      athkdf = 1.0e7

!     transition for scaling diffusion coefficient

      del_dm = 4.0/1000.0

!     half width scaling for diffusion coefficient

      s_dm = 1.0/1000.0

!-----------------------------------------------------------------------
!     provide for namelist over-ride of above settings + documentation
!-----------------------------------------------------------------------

      call getunit (io, 'control.in'
     &,               'formatted sequential rewind')
      read (io,isopyc,end=100)
100   continue
      write (stdout,isopyc)
      call relunit (io)

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

9102  format(/' ahisop = ',e12.6,' along isopyncal tracer mixing '
     &,  '(cm**2/sec) '/' athkdf = ',e12.6,' isopycnal thickness '
     &,'diffusion (cm**2/sec) '/' ah = ',e12.6,' cm**2/sec for '
     &,' am = ',e12.6,' cm**2/sec for horizontal mixing of momentum'/)
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
      write (stdout,'(5(1x,e12.6))') (dtxsqr(k),k=1,km)

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

      real sc, dzt4r, ai0, sumz, sxe, sumy, facty, ck  ! ck added by J. Getzlaff

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

# if defined O_increase_isopyc_diff
! This is the code and it's variations is based on
! Getzlaff and Dietze, 2013, GRL Vol. 40
! doi: 10.1002/grl.50419 and it increases zonal isopycnal diffusivity
!
# if defined O_increase_isopyc_diff_everywhere
      addisop=0;
      do j=48,53
       do i=2,imtm1
        do k=4,km
          addisop(i,k,j)= 5.0e8
        enddo
       enddo
      enddo
# else

# if defined O_increase_isopyc_diff_smooth
! This option was added by C. Somes to smooth out the decrease north
! and south of 5 degrees
      addisop(1:imtm1,js:je,1:km) = 0;
! Pacific
      
      do i=30,83
         do j=44,47
            do k=4,km
               addisop(i,k,j) = 1.25e8*j - 5.5e9
            enddo
         enddo
         do j=48,53
            do k=4,km
               addisop(i,k,j) = 5.0e8
            enddo
         enddo
         do j=54,57
            do k=4,km
               addisop(i,k,j) = -1.25e8*j + 7.125e9
            enddo
         enddo
      enddo     

! Indian
      
      do i=5,29
         do j=44,47
            do k=4,km
               addisop(i,k,j) = 1.25e8*j - 5.5e9
            enddo
         enddo
         do j=48,53
            do k=4,km
               addisop(i,k,j) = 5.0e8
            enddo
         enddo
         do j=54,57
            do k=4,km
               addisop(i,k,j) = -1.25e8*j + 7.125e9
            enddo
         enddo
      enddo     

! Atlantic

      do i=1,5
         do j=44,47
            do k=4,km
               addisop(i,k,j) = 5.e7*j - 2.2e9
            enddo
         enddo
         do j=48,53
            do k=4,km
               addisop(i,k,j) = 2.0e8
            enddo
         enddo
         do j=54,57
            do k=4,km
               addisop(i,k,j) = -5.e7*j + 2.85e9
            enddo
         enddo
      enddo     

      do i=84,imtm1
         do j=44,47
            do k=4,km
               addisop(i,k,j) = 5.e7*j - 2.2e9
            enddo
         enddo
         do j=48,53
            do k=4,km
               addisop(i,k,j) = 2.0e8
            enddo
         enddo
         do j=54,57
            do k=4,km
               addisop(i,k,j) = -5.e7*j + 2.85e9
            enddo
         enddo
      enddo

# else

      addisop(1:imtm1,js:je,1:km) = 0;
! Pacific

      do j=48,53
       do i=30,83
        do k=4,km
          addisop(i,k,j)= 3.0e8
        enddo
       enddo
      enddo

! Indian

      do j=48,53
       do i=5,29
        do k=4,km
          addisop(i,k,j)= 3.0e8
        enddo
       enddo
      enddo

! Atlantic

      do j=48,53
       do i=1,5
        do k=4,km
          addisop(i,k,j)= 2.0e8
        enddo
       enddo
      enddo

      do j=48,53
       do i=84,imtm1
        do k=4,km
          addisop(i,k,j)= 2.0e8
        enddo
       enddo
      enddo
# endif     
# endif
# endif

      do j=js,je
        jrow = j + joff
        do k=1,km
          sc = c1/(slmxr*dtxsqr(k))
          dzt4r = p5*dzt2r(k)
          do i=2,imtm1
# if defined O_increase_isopyc_diff
            ck = addisop(i,k,jrow)
            Ai0 =(ahisop+ck)*p5*(fisop(i,jrow,k)+fisop(i+1,jrow,k))
# else
            Ai0 = ahisop*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k))
# endif
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

      real sc, dzt4r, ai0, sumz, syn, sumx, ck

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
            Ai0 = ahisop*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k))
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

      real sc, kp1, ai0, sumx, sxb, sumy, facty, syb, ck

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
            Ai0 = ahisop*p5*(fisop(i,jrow,k+1) + fisop(i,jrow,k))
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
#  if defined O_increase_isopyc_diff
			ck = addisop(i,jrow,k)
			Ai0 = (ahisop+ck)*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k))
#  else
            Ai0 = ahisop*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k))
#  endif
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
#  if defined O_increase_isopyc_diff
            ck = addisop(i,jrow,k)
            Ai0 = (ahisop+ck)*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k))
#  else
            Ai0 = ahisop*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k))
#  endif
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

      subroutine isopyc_adv (joff, js, je, is, ie)

!=======================================================================
!     compute isopycnal transport velocities.
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, jrow, joff, km1, kp1
      integer jstrt, ie, is

#  if defined O_gent_mcwilliams
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

      do k=1,km
        top_bc(k) = c1
        bot_bc(k) = c1
      enddo
      top_bc(1)  = c0
      bot_bc(km) = c0

!-----------------------------------------------------------------------
!     compute the meridional component of the isopycnal mixing velocity
!     at the center of the northern face of the "T" cells.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff
        do k=1,km

          sc = c1/(slmxr*dtxsqr(k))
          km1 = max(k-1,1)
          kp1 = min(k+1,km)
          do i=1,imt
            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i,jrow+1,k))

            at =     (alphai(i,k,j) + alphai(i,k,j+1) + alphai(i,km1,j)
     &              + alphai(i,km1,j+1))
            bt =     (betai(i,k,j) + betai(i,k,j+1) + betai(i,km1,j)
     &              + betai(i,km1,j+1))
            stn = -(at*(ddyt(i,k,j,1) + ddyt(i,km1,j,1))
     &               + bt*(ddyt(i,k,j,2) + ddyt(i,km1,j,2))) /
     &                (at*(ddzt(i,km1,j,1) + ddzt(i,km1,j+1,1))
     &               + bt*(ddzt(i,km1,j,2) + ddzt(i,km1,j+1,2))+epsln)

            ab =     (alphai(i,k,j) + alphai(i,k,j+1) + alphai(i,kp1,j)
     &              + alphai(i,kp1,j+1))
            bb =     (betai(i,k,j) + betai(i,k,j+1) + betai(i,kp1,j)
     &              + betai(i,kp1,j+1))
            sbn = -(ab*(ddyt(i,k,j,1) + ddyt(i,kp1,j,1))
     &               + bb*(ddyt(i,k,j,2) + ddyt(i,kp1,j,2))) /
     &                (ab*(ddzt(i,k,j,1) + ddzt(i,k,j+1,1))
     &               + bb*(ddzt(i,k,j,2) + ddzt(i,k,j+1,2))+epsln)
            absstn = abs(stn)
            abssbn = abs(sbn)
#   if defined O_dm_taper
            ath_t = Ath0*tmask(i,k,j)*tmask(i,k,j+1)
     &             *p5*(c1-tanh((absstn-del_dm)*s_dmr))
            ath_b = Ath0*tmask(i,kp1,j)*tmask(i,kp1,j+1)
     &             *p5*(c1-tanh((abssbn-del_dm)*s_dmr))
#   else
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
            Ath0 = athkdf*p5*(fisop(i,jrow,k) + fisop(i+1,jrow,k))

            at =     (alphai(i,k,j) + alphai(i+1,k,j) + alphai(i,km1,j)
     &              + alphai(i+1,km1,j))
            bt =     (betai(i,k,j) + betai(i+1,k,j) + betai(i,km1,j)
     &              + betai(i+1,km1,j))
            ste =-(at*(ddxt(i,k,j,1) + ddxt(i,km1,j,1))
     &           + bt*(ddxt(i,k,j,2) + ddxt(i,km1,j,2)))
     &          / (at*(ddzt(i,km1,j,1) + ddzt(i+1,km1,j,1))
     &           + bt*(ddzt(i,km1,j,2) + ddzt(i+1,km1,j,2))+epsln)

            ab =     (alphai(i,k,j) + alphai(i+1,k,j) + alphai(i,kp1,j)
     &              + alphai(i+1,kp1,j))
            bb =     (betai(i,k,j) + betai(i+1,k,j) + betai(i,kp1,j)
     &              + betai(i+1,kp1,j))
            sbe =-(ab*(ddxt(i,k,j,1) + ddxt(i,kp1,j,1))
     &           + bb*(ddxt(i,k,j,2) + ddxt(i,kp1,j,2))) /
     &            (ab*(ddzt(i,k,j,1) + ddzt(i+1,k,j,1))
     &           + bb*(ddzt(i,k,j,2) + ddzt(i+1,k,j,2))+epsln)

            absste = abs(ste)
            abssbe = abs(sbe)
#   if defined O_dm_taper
            ath_t = Ath0*tmask(i,k,j)*tmask(i+1,k,j)
     &             *p5*(c1-tanh((absste-del_dm)*s_dmr))
            ath_b = Ath0*tmask(i,kp1,j)*tmask(i+1,kp1,j)
     &             *p5*(c1-tanh((abssbe-del_dm)*s_dmr))
#   else
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
!             remove "extra" csu in time averages of vntiso
              ta_vntiso(i,k,jrow) = ta_vntiso(i,k,jrow) +
     &                              adv_vntiso(i,k,j)/csu(jrow)
              ta_vbtiso(i,k,jrow) = ta_vbtiso(i,k,jrow) +
     &                              adv_vbtiso(i,k,j)
            enddo
          enddo
        enddo
      endif
#   endif
#  endif
# endif
#endif

      return
      end
