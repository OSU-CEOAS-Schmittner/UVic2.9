      subroutine setdata (is, ie, js, je)

!-----------------------------------------------------------------------
!     set up for data routine
!-----------------------------------------------------------------------

      implicit none

!     this is set up only for periodic monthly data.

!     before adding a surface boundary condition
!      1. check numsbc and ntdbc are large enough (in csbc.h & ctdbc.h)
!      2. add index definition to UVic_ESCM (eq. itaux)
!      3. add code below and to atmos.F

      integer i, ie, is, id, im, iou, j, je, js, k, m, n

      real c10, c100, p001, p035, realdays

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
# if defined O_embm
      include "atm.h"
# endif
#if defined O_mobi_silicon
      include "mobi.h"
#endif
      include "calendar.h"
      include "ctdbc.h"
      include "tmngr.h"
      include "switch.h"

      c10 = 10.
      c100 = 100.
      p001 = 0.001
      p035 = 0.035

      obc(:,:,:,:) = 0.

      do n=1,ntdbc
        ntdrec(n) = 12
        period(n) = .true.
        do m=1, ntdrec(n)
          k = m + 1
          aprec(m,n) = daypm(m)
          if (.not. eqyear .and. nint(aprec(m,n)) .eq. 28) then
            aprec(m,n) = aprec(m,n) + 0.2425
            print*, '=>Warning: adding 0.2425 days to feb for leap year'
          endif
!         create time stamp for the end of each month
          if (k .gt. 12) then
            call mkstmp (dstamp(m,n), 2, 1, 1, 0, 0, 0)
          else
            call mkstmp (dstamp(m,n), 1, k, 1, 0, 0, 0)
          endif
        enddo

        call timeinterpi (ntdrec(n), dstamp(1,n), aprec(1,n), tdrec(1,n)
     &,                  isbcstart(n), period(n))

        call addtime (initial, imodeltime, itemptime)
        call subtime (itemptime, isbcstart(n), itemptime)
        daysbc(n) = realdays(itemptime)
        iprevm(n) = 1
        inextm(n) = 2
        method    = 3
        call timeinterp (daysbc(n), n, tdrec(1,n), aprec(1,n), ntdrec(n)
     &,      period(n), method, inextd(n), iprevd(n), wprev(n)
     &,      rdtdbc(n), inextm(n), iprevm(n))
      enddo

      euler2 = .false.

!     load previous and next data for all boundary conditions
# if defined O_sbc_in_memory
!     this is set up only for monthly data.
      do k=1,12
# else
      do k=1,2
# endif
        n = 1
        do m=1,numsbc

          if (n .le. ntdbc) then

            rdtdbc(n) = .true.
# if defined O_sbc_in_memory
            id = k
            im = k
# else
            id = iprevd(n)
            im = iprevm(n)
            if (k .eq. 2) then
              id = inextd(n)
              im = inextm(n)
           endif
# endif

            if ( m .eq. itaux ) then
              if (k .eq. 1) print*, 'x component of windstress'
              call get_sbc (n, 'O_tau.nc', 'O_tauX', id, im, c10, c0)

            elseif ( m .eq. itauy ) then
              if (k .eq. 1) print*, 'y component of windstress'
              call get_sbc (n, 'O_tau.nc', 'O_tauY', id, im, c10, c0)

            elseif ( m .eq. iws ) then
              if (k .eq. 1) print*, 'surface wind speed'
              call get_sbc (n, 'A_windsur.nc', 'A_windspd', id, im
     &,         c100, c0)

            elseif ( m .eq. iaca ) then
              if (k .eq. 1) print*, 'atmospheric coalbedo'
              call get_sbc (n, 'A_calbatm.nc', 'A_calbatm', id, im
     &,         c1, c0)
# if defined O_mobi_iron         
            elseif ( m .eq. idfeadep ) then              
              if (k .eq. 1) print*, 'atmospheric iron deposition'
              call get_sbc (n, 'O_feflux.nc', 'O_feflux', id, im
     &,         c1, c0)  
# endif
# if defined O_mobi_silicon
            elseif ( m .eq. isilflx ) then              
              if (k .eq. 1) print*, 'atmospheric silica deposition'
              call get_sbc (n, 'O_sil_dep.nc', 'O_SIL_DEP', id, im
     &,         c100*sildustfluxfac, c0)  
# endif
# if defined O_read_my_stf
              if (k .eq. 1) print*, 'heat flux'
            elseif ( m .eq. ihflx ) then
              call get_sbc (n, 'F_heat.nc', 'F_heat', id, im, c1, c0)

            elseif ( m .eq. isflx ) then
              if (k .eq. 1) print*, 'salt flux'
              call get_sbc (n, 'F_salt.nc', 'F_salt', id, im, c1, c0)

# elif defined O_replacst
            elseif ( m .eq. isst ) then
              if (k .eq. 1) print*, 'sea surface temperature'
              call get_sbc (n, 'O_tempsur.nc', 'O_tempsur', id, im
     &,         c1, c0)

            elseif ( m .eq. isss ) then
              if (k .eq. 1) print*, 'sea surface salinity'
              call get_sbc (n, 'O_salsur.nc', 'O_salsur', id, im
     &,         p001, -p035)

# elif defined O_restorst
            elseif ( m .eq. ihflx ) then
              if (k .eq. 1) print*, 'data for restoring sst'
              call get_sbc (n, 'O_tempsur.nc', 'O_tempsur', id, im
     &,         c1, c0)

            elseif ( m .eq. isflx ) then
              if (k .eq. 1) print*, 'data for restoring sss'
              call get_sbc (n, 'O_salsur.nc', 'O_salsur', id, im
     &,         p001, -p035)

# endif
# if defined O_embm_awind &&  defined O_embm
            elseif ( m .eq. iwa ) then
              if (k .eq. 1) print*, 'surface wind angle'
              call get_sbc (n, 'A_windsur.nc', 'A_windang', id, im
     &,         c1, c0)

# endif
# if defined O_embm
            elseif ( m .eq. iwxq ) then
              if (k .eq. 1) print*, 'x component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windqX', id, im
     &,         c100, c0)

            elseif ( m .eq. iwyq ) then
              if (k .eq. 1) print*, 'y component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windqY', id, im
     &,         c100, c0)

            elseif ( m .eq. iwxt ) then
              if (k .eq. 1) print*, 'x component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windtX', id, im
     &,         c100, c0)

            elseif ( m .eq. iwyt ) then
              if (k .eq. 1) print*, 'y component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windtY', id, im
     &,         c100, c0)

#  if defined O_carbon_co2_2d
            elseif ( m .eq. iwxc ) then
              if (k .eq. 1) print*, 'x component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windcX', id, im
     &,         c100, c0)

            elseif ( m .eq. iwyc ) then
              if (k .eq. 1) print*, 'y component of advecting wind'
              call get_sbc (n, 'A_wind.nc', 'A_windcY', id, im
     &,         c100, c0)

#  endif
# endif
# if defined O_shortwave
            elseif ( m .eq. ipsw ) then
              if (k .eq. 1) print*, 'penetrating shortwave radiation'
              call get_sbc (n, 'O_penswr.nc', 'O_penswr', id, im
     &,         c1, c0)

# endif
# if defined O_mtlm
            elseif ( m .eq. idtr ) then
              if (k .eq. 1) print*, 'diurnal temperature range'
              call get_sbc (n, 'L_diurtemp.nc', 'L_diurtemp', id, im
     &,         c1, c0)

# endif
            endif

          endif

        enddo
      enddo

      return
      end

      subroutine get_sbc (n, file, name, id, im, scalar, offset)

      implicit none

      character (*) :: file, name
      character(120) :: fname, new_file_name, text

      integer id, im, iou, n, ib(10), ic(10)

      logical exists

      real offset, scalar, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "ctdbc.h"

      logical inqvardef

      real tmpij(imtm2,jmtm2)

      C2K = 273.15
      fname = new_file_name (file)
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'get_sbc in setatm.f'
      else
        ib(:) = 1
        ib(3) = id
        ic(:) = imtm2
        ic(2) = jmtm2
        ic(3) = 1
        call openfile (fname, iou)
        if (inqvardef(name, iou)) then
          call getvara (name, iou, imtm2*jmtm2, ib, ic, tmpij
     &,     scalar, offset)
          obc(2:imtm1,2:jmtm1,n,im) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (obc(:,:,n,im))
          text = "C"
          call getatttext (iou, name, 'units', text)
          if (name .ne. "L_diurtemp".and. trim(text) .eq. "K")
     &      obc(:,:,n,im) = obc(:,:,n,im) - C2K
          where (obc(:,:,n,im) .gt. 1.e30) obc(:,:,n,im) = 0.
        endif
      endif
      n = n + 1

      return
      end
