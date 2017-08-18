      subroutine reg1st (nunit, wrform, wrvol, wrk, wrmask, rmask)

#if defined O_mom
!=======================================================================

!     Subroutine reg1st is i/o routine for the user defined horizontal
!     and vertical regional masks used in MOM.
!     It is also the i/o routine that can be used to write volume
!     information for the horizontal regions used in calculating
!     volume weighted averages of tracers and surface tracer fluxes.
!     Both formatted and unformatted i/o is supported.
!     (see "cregin.h" for more details on variables).
!=======================================================================

      implicit none

!  nunit = unit to be written to or read from
!  wrform= true(false) switch for formatted(unformatted) writes
!  wrvol = switch to write volume & area information
!  wrk   = switch to write k-level volume information
!  wrmask= switch to write horizontal region masks field
!  rmask = switch to read horizontal region masks field from specified
!          unit
!  settop, setbot are used in defining level limits for vertical regions

      integer nunit, mask, k, i, linemx, linel, line, nwr, n, ia, ib
      integer jj, jjj, jr, kk

      logical  wrform, wrvol, wrk, wrmask, rmask
      logical  settop, setbot

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cregin.h"
      include "iounit.h"

      integer ncol(imt)

!-----------------------------------------------------------------------
!     write out regional volume information
!-----------------------------------------------------------------------

        if (wrvol) then
          if (wrform) then
            write(nunit,9000)
            if (wrk) then
              write(nunit,9001) (mask,mask=1,nhreg)
              do 700 k=1,km
                write(nunit,9002) k, volgk(k),
     $                           (volbk(mask,k),mask=1,nhreg)
700           continue
            endif
            write(nunit,9003) volgt, (volbt(mask),mask=1,nhreg)
            write(nunit,9004) areag, (areab(mask),mask=1,nhreg)
          else

            iotext ='read (iotavg) km, nhreg, volgt, areag'
            write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
            write (nunit) km, nhreg, volgt, areag

            iotext ='read (iotavg) (volgk(k),k=1,km)'
            write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
            call wrufio (nunit, volgk, km)

            iotext ='read (iotavg) ((volbk(l,k),l=1,nhreg),k=1,km)'
            write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
            call wrufio (nunit, volbk, nhreg*km)

            iotext ='read (iotavg) (volbt(l),l=1,nhreg)'
            write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
            call wrufio (nunit, volbt, nhreg)

            iotext ='read (iotavg) (areab(l),l=1,nhreg)'
            write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
            call wrufio (nunit, areab, nhreg)
          endif
        endif

!-----------------------------------------------------------------------
!       read in (write out) horizontal & vertical region masks
!       set linel to length of desired formatted printout line
!-----------------------------------------------------------------------

        if (wrmask .or. rmask) then

          if (wrform) then
            if (wrmask) write(nunit,9011)
            if (wrmask) write(nunit,9012)
     $         (' domain for hor mask # ',i,'=',hregnm(i),i=1,nhreg)
            if ( rmask) then
              read (nunit,9099)
              read (nunit,9099)
            endif
            if ( rmask) read (nunit,9013) (hregnm(i),i=1,nhreg)
            if (nunit .eq. stdout) then
              call iplot (mskhr, imt, imt, jmt)
            else
              linemx = 100
              linel  = 105
              line   = linel - 5
              if (line .gt. linemx) line = linemx
              nwr = (imt/line) + 1

              do 900 i=1,imt
                ncol(i) = mod(i,10)
900           continue

              do 1000 n=1,nwr
                ia = 1 + (line*(n-1))
                ib = ia + line - 1
                if (ib .gt. imt) ib = imt
                if (wrmask) write(nunit,9021) (ncol(i),i=ia,ib)
                if ( rmask) read (nunit,9099)
                do 990 jj=1,jmt
                  jjj = jmt - jj + 1
                  if (wrmask)write(nunit,9022)jjj,(mskhr(i,jjj),i=ia,ib)
                  if ( rmask) then
                    read (nunit,9022) jr , (mskhr(i,jjj),i=ia,ib)
                    if (jr .ne. jjj) then
                      write (stdout,999) nunit, jjj, jr
                      write (stderr,999) nunit, jjj, jr
                      stop '=>reg1st'
                    endif
                  endif
990             continue
1000          continue
            endif

            if (wrmask) write(nunit,9031)
            if (wrmask) write(nunit,9032)
     $         (' domain for ver mask # ',i,'=',vregnm(i),i=1,nvreg)
            if (wrmask) write(nunit,9034) mskvr

            if ( rmask) then
               read (nunit,9099)
               read (nunit,9099)
            endif
            if ( rmask) read (nunit,9033) (vregnm(i),i=1,nvreg)
            if ( rmask) read (nunit,9034) mskvr

          else
            if (wrmask) then
              iotext = ' read (nunit) mskhr, mskvr'
              write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
              write(nunit) mskhr, mskvr

              iotext = ' read (nunit) hregnm, vregnm'
              write (nunit) 'no stamp available in reg1st    ', iotext
     $,                   expnam
              write(nunit) hregnm, vregnm
            endif
            if (rmask) then
              read (nunit)
              read (nunit) mskhr, mskvr
              read (nunit)
              read (nunit) hregnm, vregnm
            endif
          endif
        endif

!       if vertical masks were read in, set level limits for defining
!       vertical regions in term balance calculations

        if (rmask) then
          do 1100 i=1,nvreg
            settop = .false.
            setbot = .false.
            do 1090 k=1,km
              kk = km-k+1
              if (mskvr(k)  .eq. i .and. .not. settop) then
                llvreg(i,1) = k
                settop = .true.
              endif
              if (mskvr(kk) .eq. i .and. .not. setbot) then
                llvreg(i,2) = kk
                setbot = .true.
              endif
1090        continue
1100      continue
        endif

999   format(/' error => bad j-row when reading regionmasks from unit ',
     $       i3,/'   expected',i4,'    read in',i4)
9000  format(/' Horizontal regional volumes [cubic m] and areas [sq m]')
9001  format('    k','  All Regions ',9(1x,i7,5x))
9002  format(1x,i4,10(1x,e12.6))
9003  format('  SUM',10(1x,e12.6))
9004  format(/' AREA',10(1x,e12.6))
9011  format(/' Horizontal regional mask names & domains:')
9012  format(a22,i4,a1,a40)
9013  format(27x,a40)
9021  format('  i=>',100(i1))
9022  format(1x,i3,1x,100(i1))
9031  format(/' Vertical regional mask names & domains:')
9032  format(a22,i4,a1,a20)
9033  format(27x,a20)
9034  format (1x, 42i3)
9099  format(1x)
#endif

      return
      end
