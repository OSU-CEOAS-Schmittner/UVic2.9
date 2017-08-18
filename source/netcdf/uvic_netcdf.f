
      subroutine opentime (fname, relyr, ntrec, ncid)
!=======================================================================
!     open file for reading or writing at record of relyr

!     input:
!       fname = file name to be opened
!       relyr = relative year
!     output:
!       ntrec = record number of relyr (next if relyr is not found)
!       ncid  = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: fname
      character(120) :: name

      integer ncid, ntrec
      integer i, id, iv, ln, n

      real, intent(in) :: relyr
      real(kind=8), allocatable :: time(:)

      logical exists, notopen

      name = fname
      inquire (file=trim(name), exist=exists)
      if (.not. exists) then
        call opennew (name, ncid)
        return
      endif
      ntrec = 0
      call openchk (fname, ncid, notopen)
      if (notopen) then
        i = nf_open (trim(name), nf_write, ncid)
        if (i .ne. nf_noerr) then
          i = nf_open (trim(name), nf_nowrite, ncid)
        endif
        call checkerror (i,'opentime nf_open '//trim(name))
        call openadd (fname, ncid)
      endif

      i = nf_inq_varid (ncid, 'time', iv)
!     return if no time variable
      if (i .ne. nf_noerr) return
      i = nf_inq_vardimid (ncid, iv, id)
      call checkerror (i,'opentime  nf_inq_vardimid'//trim(name))
      i = nf_inq_dimlen (ncid, id, ln)
      call checkerror (i,'opentime  nf_inq_dimlen'//trim(name))
      if (ln .le. 0) return
      allocate (time(ln))
!     find ntrec of relyr or next if relyr is not found
      i = nf_get_var_double (ncid, iv, time)
      call checkerror (i,'opentime nf_get_vara_double time')
      ntrec = ln + 1
      do n=1,ln
       if (abs(time(n) - relyr) .lt. 1.e-6) ntrec = n
      enddo
      deallocate (time)

      return
      end

      subroutine opennext (fname, relyr, ntrec, ncid)
!=======================================================================
!     open file for reading or writing at next record number

!     input:
!       fname = file name to be opened
!       relyr = relative year
!     output:
!       ntrec = next record number (or last if last is relyr)
!       ncid  = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: fname
      character(120) :: name

      integer ncid, ntrec
      integer i, id, iv, ln

      real, intent(in) :: relyr
      real(kind=8) time

      logical exists, notopen

      name = fname
      inquire (file=trim(name), exist=exists)
      if (.not. exists) then
        call opennew (name, ncid)
        return
      endif

      ntrec = 0
      call openchk (fname, ncid, notopen)
      if (notopen) then
        i = nf_open (trim(name), nf_write, ncid)
        if (i .ne. nf_noerr) then
          i = nf_open (trim(name), nf_nowrite, ncid)
        endif
        call checkerror (i,'opennext nf_open '//trim(name))
        call openadd (fname, ncid)
      endif

      i = nf_inq_varid (ncid, 'time', iv)
!     return if no time variable
      if (i .ne. nf_noerr) return
      i = nf_inq_vardimid (ncid, iv, id)
      call checkerror (i,'opennext  nf_inq_vardimid'//trim(name))
      i = nf_inq_dimlen (ncid, id, ln)
      call checkerror (i,'opennext  nf_inq_dimlen'//trim(name))
      if (ln .le. 0) return
!     get next record or last if last is relyr
      i = nf_get_vara_double (ncid, iv, ln, 1, time)
      call checkerror (i,'opennext nf_get_vara_double time')
      ntrec = ln + 1
      if (abs(time - relyr) .lt. 1.e-6) ntrec = ln

      return
      end

      subroutine openfile (fname, ncid)
!=======================================================================
!     open file for reading or writing

!     input:
!       fname = file name to be opened
!     output:
!       ncid  = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: fname
      character(120) :: name

      integer ncid
      integer i

      logical exists, notopen

      name = fname
      inquire (file=trim(name), exist=exists)
      if (.not. exists) then
        call opennew (name, ncid)
        return
      endif

      call openchk (fname, ncid, notopen)
      if (notopen) then
        i = nf_open (trim(name), nf_write, ncid)
        if (i .ne. nf_noerr) then
          i = nf_open (trim(name), nf_nowrite, ncid)
        endif
        call checkerror (i,'openfile nf_open '//trim(name))
        call openadd (fname, ncid)
      endif

      return
      end

      subroutine opennew (fname, ncid)
!=======================================================================
!     open file for reading or writing

!     input:
!       fname = file name to be opened
!     output:
!       ncid  = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: fname
      character(120) :: name

      integer ncid
      integer i

      logical notopen

      name = fname
      call openchk (fname, ncid, notopen)
      if (notopen) then
        i = nf_create (trim(name), nf_clobber, ncid)
        call checkerror (i,'openfile '//trim(name))
        i = nf_enddef (ncid)
        call checkerror (i,'openfile nf_enddef')
        call openadd (fname, ncid)
      endif

      return
      end

      subroutine openchk (fname, ncid, notopen)
!=======================================================================
!     check if file is open and add to list if not
!     input:
!       fname = file name to be opened
!     output:
!       notopen  = true if file is not open
!=======================================================================

      implicit none

      include 'netcdf.inc'
      include 'uvic_netcdf.h'

      integer ncid
      integer n

      character(*), intent(in) :: fname

      logical notopen

      notopen = .true.
      ncid = 0
      do n=1,num
        if (fname .eq. list_names(n)) then
          notopen = .false.
          ncid = list_ncid(n)
        endif
      enddo

      return
      end

      subroutine openadd (fname, ncid)
!=======================================================================
!     add new file to open list
!     input:
!       fname = file name to be opened
!       ncid = id of file to be opened
!=======================================================================

      implicit none

      include 'netcdf.inc'
      include 'uvic_netcdf.h'

      character(*), intent(in) :: fname

      integer, intent(in) :: ncid

      num = num + 1
      if (num .gt. max_num_files) then
        print*, "=> Error: increase max_num_files in uvic_netcf.f"
        stop
      endif
      list_names(num) = fname
      list_ncid(num) = ncid

      return
      end

      subroutine closeall
!=======================================================================
!     close all netcdf files
!=======================================================================

      implicit none

      include 'netcdf.inc'
      include 'uvic_netcdf.h'

      integer i, n

      do n=1,num
        i = nf_close (list_ncid(n))
        call checkerror (i,'closeall nf_close')
        list_names(n) = " "
        list_ncid(n) = 0
      enddo
      num = 0

      return
      end

      subroutine closefile (ncid)
!=======================================================================
!     close file

!     input:
!       ncid = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'
      include 'uvic_netcdf.h'

      integer, intent(in) :: ncid
      integer i, m, n

      m = 1
      do n=1,num
        if (ncid .eq. list_ncid(n)) then
          i = nf_close (ncid)
          call checkerror (i,'closefile nf_close')
          m = n + 1
          list_names(n) = list_names(m)
          list_ncid(n) = list_ncid(m)
        endif
      enddo
      num =  2*num - m

      return
      end

      subroutine redef (ncid)
!=======================================================================
!     redefine

!     input:
!       ncid = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      integer, intent(in) :: ncid
      integer i

      i = nf_redef(ncid)
      call checkerror (i,'redef nf_redef')

      return
      end

      subroutine enddef (ncid)
!=======================================================================
!     end definitions

!     input:
!       ncid = iou unit
!=======================================================================

      implicit none

      include 'netcdf.inc'

      integer, intent(in) :: ncid
      integer i

      i = nf_enddef (ncid)
      call checkerror (i,' enddef nf_enddef')

      return
      end

      subroutine checkerror(i, trace)
!=======================================================================
!     check for any netcdf errors

!     input:
!       i     = netcdf error index
!       trace = trace string
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: trace

      integer, intent(in) :: i

      if (i .ne. nf_noerr) then
        print*, 'netcdf error: ', nf_strerror(i)
        print*, 'trace string: ', trace
        stop
      endif

      return
      end

      function inqvardef (name, ncid)
!=======================================================================
!     check if a variable is defined

!     input:
!       name = variable name
!       ncid = iou unit
!     output:
!        inqvardef = (true, false) = (defined, not defined)
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ncid
      integer i, iv

      logical inqvardef

      inqvardef = .false.
      i = nf_inq_varid (ncid, name, iv)
      if (i .eq. nf_noerr) inqvardef = .true.

      return
      end

      subroutine putatttext (ncid, var, name, text)
!=======================================================================
!     put text attribute

!     input:
!       ncid = iou unit
!       var  = variable name ("global" for a global attribute)
!       name = text name
!       text = text
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: var, name, text

      integer, intent(in) :: ncid
      integer i, iv

      if (var .eq. "global" .or. var .eq. "Global") then
        i = nf_put_att_text (ncid, nf_global, trim(name)                                                                             &
     &,   len(trim(text)), trim(text))
        call checkerror(i,'putatttext global '//trim(name))
      else
        i = nf_inq_varid (ncid, var, iv)
        call checkerror (i,'putatttext nf_inq_varid '//trim(var))
        i = nf_put_att_text (ncid, iv, trim(name)                                                                                    &
     &,   len(trim(text)), trim(text))
        call checkerror(i,'putatttext '//trim(var)//' '//trim(name))
      endif

      return
      end

      subroutine getatttext (ncid, var, name, text)
!=======================================================================
!     get text attribute

!     input:
!       ncid = iou unit
!       var  = variable name ("global" for a global attribute)
!       name = text name
!     output
!       text = text
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: var, name
      character(*) text

      integer, intent(in) :: ncid
      integer i, iv

      if (var .eq. "global" .or. var .eq. "Global") then
        i = nf_get_att_text (ncid, nf_global, trim(name), text)
        if (i .ne. nf_noerr) then
          print*,'getatttext: global ',trim(name),' not found'
        endif
      else
        i = nf_inq_varid (ncid, var, iv)
        call checkerror (i,'getatttext nf_inq_varid '//trim(var))
        i = nf_get_att_text (ncid, iv, trim(name), text)
        if (i .ne. nf_noerr) then
          print*,'getatttext: ',trim(var),' ',trim(name),' not found'
        endif
      endif

      return
      end

      subroutine defdim (name, ncid, ln, id)
!=======================================================================
!     define dimension

!     input:
!       name = name of variable to be defined
!       ncid = iou unit
!       ln   = length of axis (0 = unlimited)
!       id   = dimension id
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: id, ln, ncid
      integer i

      i = nf_inq_dimid (ncid, name, id)
!     if dimension is already defined, return
      if (i .eq. nf_noerr) return

      if (ln .gt. 0) then
        i = nf_def_dim (ncid, name, ln, id)
      else
        i = nf_def_dim (ncid, name, nf_unlimited, id)
      endif
      call checkerror (i, 'defdim '//trim(name))

      return
      end

      subroutine getdimlen (name, ncid, ln)
!=======================================================================
!     define dimension

!     input:
!       name = name of dimension
!       ncid = iou unit
!     output:
!       ln   = length of dimension
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ncid
      integer ln
      integer i, id

      i = nf_inq_dimid (ncid, name, id)
      call checkerror (i,'getdimlen nf_inq_dimid '//name)
      i = nf_inq_dimlen (ncid, id, ln)
      call checkerror (i,'getdimlen nf_inq_dimlen '//name)

      return
      end

      subroutine getaxis (name, ncid, ln, dout, is, ie, s, o)
!=======================================================================
!     read a netcdf axis
!     the first value of the axis to be read must be found within the
!     global axis provided or the axis will be redefined. if the axis is
!     redefined, it will be centred in the global axis and padded with
!     nf_fill_double. if the read axis is larger than the global axis or
!     is defined outside of the global axis a stop error is generated.

!     input:
!       name = name of variable to be defined
!       ncid = unit
!       ln   = length of axis
!       s    = data scalar
!       o    = data offset
!     output:
!       dout = global axis
!       is   = starting index in global axis
!       ie   = ending index in global array
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ln, ncid
      integer ie, is
      integer i, id, iv, len

      real, intent(in) :: o, s
      real dout(ln)
      real rs
      real(kind=8), allocatable :: din(:)

      i = nf_inq_varid (ncid, name, iv)
      call checkerror (i,'getaxis nf_inq_varid '//name)
      i = nf_inq_vardimid (ncid, iv, id)
      call checkerror (i,'getaxis nf_inq_varndimid '//name)
      i = nf_inq_dimlen (ncid, id, len)
      call checkerror (i,'getaxis nf_inq_dimlen '//name)
      allocate (din(len))
      i = nf_get_vara_double (ncid, iv, 1, len, din)
      call checkerror(i,'getaxis nf_get_vara_double '//name)
      is = 0
      do i=ln,1,-1
        if (abs(dout(i)-din(1)) .lt. 1.e-5) is = i
      enddo
      if (is .eq. 0) then
        dout(:) = nf_fill_double
        is = 1
        if (len .lt. ln) is = is + (ln - len)/2
      endif
      if (len + is - 1 .gt. ln) then
        print*, 'error in getaxis => read axis not within global axis'
        stop
      endif
      ie = is - 1 + len
      rs = 0.0
      if (s .ne. 0.) rs = 1.0/s
      do i=1,len
        dout(i + is - 1) = (din(i) - o)*rs
      enddo
      deallocate (din)

      return
      end

      subroutine defvar (name, ncid, nd, id, rmin, rmax, axis                                                                        &
     &,                  type, lname, sname, units)
!=======================================================================
!     define data

!     input:
!       name  = name of variable to be defined
!       ncid  = unit
!       nd    = number dimensions of data
!       id    = data id
!       rmin  = minimum range (default real)
!       rmax  = maximum range (default real)
!       axis  = axis type
!       type  = data type (D=double,F=float,I=integer,Tn=char*n)
!       lname = long name
!       sname = standard name
!       units = data units
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name, axis, lname, sname, type, units

      integer, intent(in) :: nd, id(nd), ncid
      integer i, idt(nd+1), iv, ln, ivar(2)

      real, intent(in) :: rmax, rmin
      real(kind=4) fvar(2)
      real(kind=8) dvar(2)

      i = nf_inq_varid (ncid, name, iv)
!     if variable is already defined, return
      if (i .eq. nf_noerr) return

      if (type .eq. 'D') then
        i = nf_def_var (ncid, name, nf_double, nd, id, iv)
        call checkerror (i,'defvar double '//trim(name))
        if (rmin .ne. rmax) then
          dvar(1) = dble(rmin)
          dvar(2) = dble(rmax)
          i = nf_put_att_double(ncid,iv,'valid_range',nf_double,2,dvar)
          call checkerror(i,'defvar valid_range double '//trim(name))
        endif
        i = nf_put_att_double (ncid,iv,'FillValue',nf_double,1                                                                      &
     &,   nf_fill_double)
        call checkerror (i,'defvar FillValue double '//trim(name))
        call checkerror (i,'defvar missing_value double '//trim(name))
        i = nf_put_att_double (ncid,iv,'missing_value',nf_double,1                                                                   &
     &,   nf_fill_double)

      elseif (type .eq. 'F') then
        i = nf_def_var (ncid, name, nf_real, nd, id, iv)
        call checkerror (i,'defvar real '//name)
        if (rmin .ne. rmax) then
          fvar(1) = real(rmin)
          fvar(2) = real(rmax)
          i = nf_put_att_real (ncid,iv,'valid_range',nf_real,2,fvar)
          call checkerror (i,'defvar valid_range real '//trim(name))
        endif
        i = nf_put_att_double (ncid,iv,'FillValue',nf_real,1                                                                        &
     &,   nf_fill_double)
        call checkerror (i,'defvar FillValue real '//trim(name))
        i = nf_put_att_double (ncid,iv,'missing_value',nf_real,1                                                                     &
     &,   nf_fill_double)
        call checkerror (i,'defvar missing_value real '//trim(name))

      elseif (type .eq. 'I') then
        i = nf_def_var (ncid, name, nf_int, nd, id, iv)
        call checkerror (i,'defvar integer '//trim(name))
        if (rmin .ne. rmax) then
          ivar(1) = int(rmin)
          ivar(2) = int(rmax)
          i = nf_put_att_int (ncid,iv,'valid_range',nf_int,2,ivar)
          call checkerror (i,'defvar valid_range integer '//trim(name))
        endif
        i = nf_put_att_int (ncid,iv,'FillValue',nf_int,1                                                                            &
     &,   nf_fill_int)
        call checkerror (i,'defvar FillValue integer '//trim(name))
        i = nf_put_att_int (ncid,iv,'missing_value',nf_int,1                                                                         &
     &,   nf_fill_int)
        call checkerror (i,'defvar missing_value integer '//trim(name))

      elseif (type(1:1) .eq. 'T') then
        ln = 0
        do i=2,len(type)
          ln = ln*10.0 +  ichar(type(i:i)) - 48
        enddo
        if (ln .le. 0 .or. ln .ge. 1000) ln = 80
        do i=1,nd
         idt(i+1) = id(i)
        enddo
        call defdim (type, ncid, ln, idt(1))
        i = nf_def_var (ncid, name, nf_char, 2, idt, iv)
        call checkerror (i,'defvar text '//trim(name))
      endif

      if (axis .ne. ' ') then
        i = nf_put_att_text (ncid,iv,'axis',len(axis),axis)
        call checkerror (i,'defvar axis '//trim(name))
        if (axis .ne. 'T') then
          i = max(len(trim(name))-5,1)
          ln = len(name)-len(trim(name))
          if (name(i-ln:len(trim(name))-ln) .ne. '_edges') then
            i = nf_put_att_text (ncid,iv,'edges'
     &,       len(trim(name)//"_edges"), trim(name)//"_edges")
            call checkerror (i,'defvar edges '//trim(name))
          endif
        endif
      endif
      if (lname .ne. ' ') then
        i = nf_put_att_text (ncid,iv,'long_name',len(lname),lname)
        call checkerror (i,'defvar long_name '//trim(name))
      endif
      if (sname .ne. ' ') then
        i = nf_put_att_text (ncid,iv,'standard_name',len(sname),sname)
        call checkerror(i,'defvar standard_name '//trim(name))
      endif
      if (units .ne. ' ') then
        i = nf_put_att_text (ncid,iv,'units',len(units),units)
        call checkerror(i,'defvar units '//trim(name))
      endif

      return
      end

      subroutine putvaramsk (name, ncid, ln, is, ic, din, dm, s, o)
!=======================================================================
!     write data

!     input:
!       name = name of variable to be written
!       ncid = iou unit
!       ln   = length of data
!       is   = starting point for write in each dimension
!       ic   = count (or length) for write in each dimension
!       din  = data to be written (default real)
!       dm   = data mask
!       s    = data scalar
!       o    = data offset
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ic(10), is(10), ln, ncid
      integer i, iv, nd

      real, intent(in) :: din(ln), dm(ln), o, s
      real rs
      real(kind=8) dout(ln)

      i = nf_inq_varid (ncid, name, iv)
      call checkerror (i,'putvaramsk nf_inq_varid '//name)

      rs = 0.0
      if (s .ne. 0.) rs = 1.0/s
      do i=1,ln
        if (dm(i) .ge. 0.5) then
          dout(i) = (din(i) - o)*rs
        else
          dout(i) = nf_fill_double
        endif
      enddo
      i = nf_inq_varid (ncid, name, iv)
      call checkerror (i,'putvaramsk nf_inq_varid '//name)
      i = nf_inq_varndims(ncid, iv, nd)
      call checkerror (i,'putvaramsk nf_inq_varndims '//name)
      i = nf_put_vara_double (ncid, iv, is(1:nd), ic(1:nd), dout)
      call checkerror(i,'putvaramsk '//name)

      return
      end

      subroutine putvara (name, ncid, ln, is, ic, din, s, o)
!=======================================================================
!     write data

!     input:
!       name = name of variable to be written
!       ncid = iou unit
!       ln   = length of data
!       is   = starting point for write in each dimension
!       ic   = count (or length) for write in each dimension
!       din  = data to be written (default real)
!       s    = data scalar
!       o    = data offset
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: is(10), ic(10), ln, ncid
      integer i, iv, nd

      real, intent(in) :: din(ln), o, s
      real rs
      real(kind=8) dout(ln)

      rs = 0.0
      if (s .ne. 0.) rs = 1.0/s
      do i=1,ln
        dout(i) = (din(i) - o)*rs
      enddo
      i = nf_inq_varid (ncid, name, iv)
      call checkerror (i,'putvara nf_inq_varid '//name)
      i = nf_inq_varndims(ncid, iv, nd)
      call checkerror (i,'putvara nf_inq_varndims '//name)
      i = nf_put_vara_double (ncid, iv, is(1:nd), ic(1:nd), dout)
      call checkerror(i,'putvara '//name)

      return
      end

      subroutine getvara (name, ncid, ln, is, ic, dout, s, o)
!=======================================================================
!     read data

!     input:
!       name = name of variable to be written
!       ncid = iou unit
!       ln   = length of data
!       is   = starting point for read in each dimension
!       ic   = count (or length) for read in each dimension
!       s    = data scalar
!       o    = data offset
!     output:
!       dout = data (default real)
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: is(10), ic(10), ln, ncid
      integer i, iv, nd

      real, intent(in) :: o, s
      real dout(ln)
      real(kind=8) din(ln), offset, scale

      i = nf_inq_varid (ncid, name, iv)
      if (i .ne. nf_noerr) then
        print*, '==> Warning: netcdf variable ',trim(name), ' not found'
        return
      endif
      scale = 1.0
      offset = 0.0
      i = nf_inq_varndims(ncid, iv, nd)
      call checkerror (i,'getvara nf_inq_varndims '//name)
      i = nf_get_att_double (ncid, iv, 'add_offset', offset)
      i = nf_get_att_double (ncid, iv, 'scale_factor', scale)
      i = nf_get_vara_double (ncid, iv, is(1:nd), ic(1:nd), din)
      call checkerror(i,'getvara '//name)
      dout(1:ln) = (din(1:ln)*scale + offset)*s + o

      return
      end

      subroutine putvars (name, ncid, is, din, s, o)
!=======================================================================
!     write scalar data

!     input:
!       name = name of variable to be written
!       ncid = iou unit
!       is   = starting point for write
!       din  = data to be written (default real)
!       s    = data scalar
!       o    = data offset
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ncid, is
      integer i, iv

      real, intent(in) :: din, o, s
      real rs
      real(kind=8) dout

      rs = 0.0
      if (s .ne. 0.) rs = 1.0/s
      dout = (din - o)*rs
      i = nf_inq_varid (ncid, name, iv)
      call checkerror (i,'putvars nf_inq_varid '//name)
      i = nf_put_vara_double (ncid, iv, is, 1, dout)
      call checkerror(i,'putvars '//name)

      return
      end

      subroutine getvars (name, ncid, is, dout, s, o)
!=======================================================================
!     read scalar data

!     input:
!       name = name of variable to be read
!       ncid = iou unit
!       is   = starting point for read
!       s    = data scalar
!       o    = data offset
!     output:
!       dout = data (default real)
!=======================================================================

      implicit none

      include 'netcdf.inc'

      character(*), intent(in) :: name

      integer, intent(in) :: ncid, is
      integer i, iv

      real, intent(in) :: o, s
      real dout
      real(kind=8) din, offset, scale

      i = nf_inq_varid (ncid, name, iv)
!     return zero for data if variable is not found
      if (i .ne. nf_noerr) then
        print*, '==> Warning: netcdf variable ',trim(name), ' not found'
        return
      endif
      scale = 1.0
      offset = 0.0
      i = nf_get_att_double (ncid, iv, 'add_offset', offset)
      i = nf_get_att_double (ncid, iv, 'scale_factor', scale)
      i = nf_get_vara_double (ncid, iv, is, 1, din)
      call checkerror(i,'getvara '//name)
      dout = (din*scale + offset)*s + o

      return
      end

      subroutine edge_maker (it, edges, xt, dxt, xu, dxu, imt)
!=======================================================================
!     make edges for grid cells

!     input:
!       it  = flag for grid (t=1, u=2)
!       xt  = t grid position array
!       dxt = t grid width
!       xu  = u grid position array
!       dxu = u grid width
!       imt = array size
!     output:
!       edges = edge array
!=======================================================================

      implicit none

      integer, intent(in) :: imt, it
      integer i
      real, intent(in) :: xt(imt), dxt(imt), xu(imt), dxu(imt)
      real edges(0:imt)

      if (it .eq. 1) then
!       make edges for T cells
        edges(0) = xu(1) - dxt(1)
        do i=1,imt
          edges(i) = xu(i)
        enddo
      elseif (it .eq. 2) then
!       make edges for U cells
        edges(imt) = xt(imt) + dxu(imt)
        do i=1,imt
          edges(i-1) = xt(i)
        enddo
      else
        write (*,*) 'Error:  it = ',it, ' in edge_maker'
        stop
      endif

      return
      end
