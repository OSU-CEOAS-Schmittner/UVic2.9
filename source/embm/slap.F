      subroutine slap_sslugm(n, b, x, nelt, ia, ja, a, isym, nsave,
     $     itol, tol, itmax, iter, err, ierr, iunit, rwork, lenw,
     $     iwork, leniw )

#if defined O_embm && defined O_embm_slap

      implicit none

      integer locrb, locib
      parameter (locrb=1, locib=11)
      integer n, nelt, ia(nelt), ja(nelt), isym, nsave, itol
      integer itmax, iter, ierr, iunit, lenw, iwork(leniw), leniw
      integer lociw, locl, locnc, locju, locnr, myitol
      integer locw, locrgw, locdin, locu, locjl, lociu, locil
      integer j, locigw, jbgn, jend, icol, nl, nu

      real b(n), x(n), a(nelt), tol, err, rwork(lenw)

!         change the slap input matrix ia, ja, a to slap-column format.
      ierr = 0
      err  = 0.0
      if ( nsave.le.1 ) then
         ierr = 3
         return
      endif
!      write(*,*) 'sslugm:  changing data structure'
      call ss2y( n, nelt, ia, ja, a, isym )

!         count number of non-zero elements preconditioner ilu matrix.
!         then set up the work arrays.  we assume maxl=kmp=nsave.
!      write(*,*) 'sslugm:  counting non-zeros'
      nl = 0
      nu = 0
      do 20 icol = 1, n
!         don't count diagonal.
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if ( jbgn.le.jend ) then

            do 10 j = jbgn, jend
               if ( ia(j).gt.icol ) then
                  nl = nl + 1
                  if ( isym.ne.0 ) nu = nu + 1
               else
                  nu = nu + 1
               endif
 10         continue
         endif
 20   continue

!      write(*,*) 'nl,nu =', nl,nu
      locigw = locib
      locil = locigw + 20
      locjl = locil + n+1
      lociu = locjl + nl
      locju = lociu + nu
      locnr = locju + n+1
      locnc = locnr + n
      lociw = locnc + n

      locl = locrb
      locdin = locl + nl
      locu = locdin + n
      locrgw = locu + nu
      locw = locrgw + 1+n*(nsave+6)+nsave*(nsave+3)

      iwork(1) = locil
      iwork(2) = locjl
      iwork(3) = lociu
      iwork(4) = locju
      iwork(5) = locl
      iwork(6) = locdin
      iwork(7) = locu
      iwork(9) = lociw
      iwork(10) = locw

!         compute the incomplete lu decomposition.
!      write(*,*) 'sslugm:  computing ilu factors'
      call ssilus( n, nelt, ia, ja, a, isym, nl, iwork(locil),
     $     iwork(locjl), rwork(locl), rwork(locdin), nu,
     $           iwork(lociu),
     $     iwork(locju), rwork(locu), iwork(locnr), iwork(locnc) )

!         perform the incomplet lu preconditioned generalized minimum
!         residual iteration algorithm.  the following sgmres
!         defaults are used maxl = kmp = nsave, jscal = 0,
!         jpre = -1, nrmax = itmax/nsave
      iwork(locigw  ) = nsave
      iwork(locigw+1) = nsave
      iwork(locigw+2) = 0
      iwork(locigw+3) = -1
      iwork(locigw+4) = itmax/nsave
      myitol = 0

      call sgmres( n, b, x, nelt, ia, ja, a, isym, myitol,
     $     tol, itmax, iter, err, ierr, iunit, rwork, rwork,
     $     rwork(locrgw), lenw-locrgw, iwork(locigw), 20,
     $     rwork, iwork )

      if ( iter .gt.itmax ) ierr = 2

      return
      end

      subroutine ss2y(n, nelt, ia, ja, a, isym )

!***purpose  slap triad to slap column format converter.
!            routine to convert from the slap triad to slap column
!            format.
! *arguments:
! n      :in       integer
!         order of the matrix.
! nelt   :in       integer.
!         number of non-zeros stored in a.
! ia     :inout    integer ia(nelt).
! ja     :inout    integer ja(nelt).
! a      :inout    real a(nelt).
!         these arrays should hold the matrix a in either the slap
!         triad format or the slap column format.  see "long
!         description", below.  if the slap triad format is used
!         this format is translated to the slap column format by
!         this routine.
! isym   :in       integer.
!         flag to indicate symmetric storage format.
!         if isym=0, all nonzero entries of the matrix are stored.
!         if isym=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.

!       the sparse linear algebra package (slap) utilizes two matrix
!       data structures: 1) the  slap triad  format or  2)  the slap
!       column format.  the user can hand this routine either of the
!       of these data structures.  if the slap triad format is give
!       as input then this routine transforms it into slap column
!       format.  the way this routine tells which format is given as
!       input is to look at ja(n+1).  if ja(n+1) = nelt+1 then we
!       have the slap column format.  if that equality does not hold
!       then it is assumed that the ia, ja, a arrays contain the
!       slap triad format.

      implicit none

      integer n, nelt, ia(nelt), ja(nelt), isym
      integer i, itemp, ibgn, iend, icol, j

      real a(nelt), temp

!         check to see if the (ia,ja,a) arrays are in slap column
!         format.  if it's not then transform from slap triad.
      if ( ja(n+1).eq.nelt+1 ) return

!         sort into ascending order by column (on the ja array).
!         this will line up the columns.

      call qs2i1r( ja, ia, a, nelt, 1 )

!         loop over each column to see where the column indicies change
!         in the column index array ja.  this marks the beginning of the
!         next column.

      ja(1) = 1
      do 20 icol = 1, n-1
         do 10 j = ja(icol)+1, nelt
            if ( ja(j).ne.icol ) then
               ja(icol+1) = j
               goto 20
            endif
 10      continue
 20   continue
      ja(n+1) = nelt+1

!         mark the n+2 element so that future calls to a slap routine
!         utilizing the ysmp-column storage format will be able to tell.

      ja(n+2) = 0

!         now loop thru the ia(i) array making sure that the diagonal
!         matrix element appears first in the column.  then sort the
!         rest of the column in ascending order.

      do 70 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
         do 30 i = ibgn, iend
            if ( ia(i).eq.icol ) then
!         swap the diag element with the first element in the column.
               itemp = ia(i)
               ia(i) = ia(ibgn)
               ia(ibgn) = itemp
               temp = a(i)
               a(i) = a(ibgn)
               a(ibgn) = temp
               goto 40
            endif
 30      continue
 40      ibgn = ibgn + 1
         if ( ibgn.lt.iend ) then
            do 60 i = ibgn, iend
               do 50 j = i+1, iend
                  if ( ia(i).gt.ia(j) ) then
                     itemp = ia(i)
                     ia(i) = ia(j)
                     ia(j) = itemp
                     temp = a(i)
                     a(i) = a(j)
                     a(j) = temp
                  endif
 50            continue
 60         continue
         endif
 70   continue

      return
      end

      subroutine sslui(n, b, x, nelt, ia, ja, a, isym, rwork, iwork)

      implicit none

      integer n, nelt, ia(nelt), ja(nelt), isym, iwork(*)
      integer locl, locdin, locu, locju, locil, locjl, lociu

      real b(n), x(n), a(nelt), rwork(*)

!         pull out the locations of the arrays holding the ilu
!         factorization.
      locil = iwork(1)
      locjl = iwork(2)
      lociu = iwork(3)
      locju = iwork(4)
      locl = iwork(5)
      locdin = iwork(6)
      locu = iwork(7)

!         solve the system lux = b
      call sslui2(n, b, x, iwork(locil), iwork(locjl), rwork(locl),
     $     rwork(locdin), iwork(lociu), iwork(locju), rwork(locu) )

      return
      end

      subroutine sslui2(n, b, x, il, jl, l, dinv, iu, ju, u )

!***purpose  slap back solve for ldu factorization.
!            routine  to  solve a system of the form  l*d*u x  =  b,
!            where l is a unit  lower  triangular  matrix,  d  is  a
!            diagonal matrix, and u is a unit upper triangular matrix.
! *arguments:
! n      :in       integer
!         order of the matrix.
! b      :in       real b(n).
!         right hand side.
! x      :out      real x(n).
!         solution of l*d*u x = b.
! nel    :in       integer.
!         number of non-zeros in the el array.
! il     :in       integer il(n+1).
! jl     :in       integer jl(nl).
!  l     :in       real     l(nl).
!         il, jl, l contain the unit  lower triangular factor of the
!         incomplete decomposition of some matrix stored in slap row
!         format.  the diagonal of ones *is* stored.  this structure
!         can   be   set up  by   the  ssilus routine.   see
!         "description", below  for more   details about   the  slap
!         format.
! dinv   :in       real dinv(n).
!         inverse of the diagonal matrix d.
! nu     :in       integer.
!         number of non-zeros in the u array.
! iu     :in       integer iu(n+1).
! ju     :in       integer ju(nu).
! u      :in       real     u(nu).
!         iu, ju, u contain the unit upper triangular factor  of the
!         incomplete decomposition  of  some  matrix stored in  slap
!         column format.   the diagonal of ones  *is* stored.   this
!         structure can be set up  by the ssilus routine.  see
!         "description", below   for  more   details about  the slap
!         format.

! *description:
!       this routine is supplied with  the slap package as a routine
!       to  perform  the  sslui operation  in   the  sir and   sbcg
!       iteration routines for  the  drivers ssilur and sslubc.   it
!       must  be called  via   the  slap  sslui  calling   sequence
!       convention interface routine sslui.
!         **** this routine itself does not conform to the ****
!               **** slap sslui calling convention ****

!       il, jl, l should contain the unit lower triangular factor of
!       the incomplete decomposition of the a matrix  stored in slap
!       row format.  iu, ju, u should contain  the unit upper factor
!       of the  incomplete decomposition of  the a matrix  stored in
!       slap column format this ilu factorization can be computed by
!       the ssilus routine.  the diagonals (which is all one's) are
!       stored.

!       =================== s l a p column format ==================
!       this routine  requires that  the matrix a  be stored in  the
!       slap column format.  in this format the non-zeros are stored
!       counting down columns (except  for the diagonal entry, which
!       must appear first in each  "column") and  are stored  in the
!       real array a.  in other words, for each column in the matrix
!       put the diagonal entry in a.  then put in the other non-zero
!       elements going down   the  column (except  the diagonal)  in
!       order.  the ia array holds the row  index for each non-zero.
!       the ja array holds the offsets into the ia, a arrays for the
!       beginning of   each    column.    that  is,    ia(ja(icol)),
!       a(ja(icol)) points to the beginning of the icol-th column in
!       ia and  a.  ia(ja(icol+1)-1),  a(ja(icol+1)-1) points to the
!       end  of   the icol-th  column.  note   that  we  always have
!       ja(n+1) = nelt+1, where  n  is the number of columns in  the
!       matrix and  nelt   is the number of non-zeros in the matrix.

!       here is an example of the  slap column  storage format for a
!       5x5 matrix (in the a and ia arrays '|'  denotes the end of a
!       column):

!           5x5 matrix      slap column format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   a: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
!       |21 22  0  0  0|  ia:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  ja:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|

      implicit none

      integer n, il(*), jl(*), iu(*), ju(*)
      integer i, icol, irow, j, jbgn, jend

      real b(n), x(n), l(*), dinv(n), u(*)

!         solve  l*y = b,  storing result in x, l stored by rows.
      do 10 i = 1, n
         x(i) = b(i)
 10   continue
      do 30 irow = 2, n
         jbgn = il(irow)
         jend = il(irow+1)-1
         if ( jbgn.le.jend ) then
            do 20 j = jbgn, jend
               x(irow) = x(irow) - l(j)*x(jl(j))
 20         continue
         endif
 30   continue

!         solve  d*z = y,  storing result in x.
      do 40 i=1,n
         x(i) = x(i)*dinv(i)
 40   continue

!         solve  u*x = z, u stored by columns.
      do 60 icol = n, 2, -1
         jbgn = ju(icol)
         jend = ju(icol+1)-1
         if ( jbgn.le.jend ) then
            do 50 j = jbgn, jend
               x(iu(j)) = x(iu(j)) - u(j)*x(icol)
 50         continue
         endif
 60   continue

      return
      end

      subroutine ssilus(n, nelt, ia, ja, a, isym, nl, il, jl,
     $     l, dinv, nu, iu, ju, u, nrow, ncol)

!***purpose  incomplete lu decomposition preconditioner slap set up.
!            routine to generate the incomplete ldu decomposition of a
!            matrix.  the  unit lower triangular factor l is stored by
!            rows and the  unit upper triangular factor u is stored by
!            columns.  the inverse of the diagonal matrix d is stored.
!            no fill in is allowed.
! *arguments:
! n      :in       integer

!         order of the matrix.
! nelt   :in       integer.
!         number of elements in arrays ia, ja, and a.
! ia     :in       integer ia(nelt).
! ja     :in       integer ja(nelt).
! a      :in       real a(nelt).
!         these arrays should hold the matrix a in the slap column
!         format.  see "description", below.
! isym   :in       integer.
!         flag to indicate symmetric storage format.
!         if isym=0, all nonzero entries of the matrix are stored.
!         if isym=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! nl     :out      integer.
!         number of non-zeros in the el array.
! il     :out      integer il(n+1).
! jl     :out      integer jl(nl).
! l      :out      real     l(nl).
!         il, jl, l  contain the unit ower  triangular factor of  the
!         incomplete decomposition  of some  matrix stored  in   slap
!         row format.     the   diagonal  of ones  *is*  stored.  see
!         "description", below for more details about the slap format.
! nu     :out      integer.
!         number of non-zeros in the u array.
! iu     :out      integer iu(n+1).
! ju     :out      integer ju(nu).
! u      :out      real     u(nu).
!         iu, ju, u contain   the unit upper triangular factor of the
!         incomplete  decomposition    of some matrix  stored in slap
!         column  format.   the diagonal of ones   *is*  stored.  see
!         "description", below  for  more  details  about  the   slap
!         format.
! nrow   :work     integer nrow(n).
!         nrow(i) is the number of non-zero elements in the i-th row
!         of l.
! ncol   :work     integer ncol(n).
!         ncol(i) is the number of non-zero elements in the i-th
!         column of u.

! *description
!       il, jl, l should contain the unit  lower triangular factor of
!       the incomplete decomposition of the a matrix  stored in slap
!       row format.  iu, ju, u should contain  the unit upper factor
!       of the  incomplete decomposition of  the a matrix  stored in
!       slap column format this ilu factorization can be computed by
!       the ssilus routine.  the diagonals (which is all one's) are
!       stored.

      implicit none

      integer n, nelt, ia(nelt), ja(nelt), isym, nl, il(nl), jl(nl)
      integer nu, iu(nu), ju(nu), nrow(n), ncol(n)
      integer i, ibgn, iend, j, jbgn, jend, icol, irow
      integer itemp , jtemp, indx, indx1, indx2, indxr1, indxr2
      integer indxc2, indxc1, k, kc, kr

      real a(nelt), l(nl), dinv(n), u(nu), temp

!         count number of elements in each row of the lower triangle.
      do 10 i=1,n
         nrow(i) = 0
         ncol(i) = 0
 10   continue
      do 30 icol = 1, n
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if ( jbgn.le.jend ) then
            do 20 j = jbgn, jend
               if ( ia(j).lt.icol ) then
                  ncol(icol) = ncol(icol) + 1
               else
                  nrow(ia(j)) = nrow(ia(j)) + 1
                  if ( isym.ne.0 ) ncol(ia(j)) = ncol(ia(j)) + 1
               endif
 20         continue
         endif
 30   continue
      ju(1) = 1
      il(1) = 1
      do 40 icol = 1, n
         il(icol+1) = il(icol) + nrow(icol)
         ju(icol+1) = ju(icol) + ncol(icol)
         nrow(icol) = il(icol)
         ncol(icol) = ju(icol)
 40   continue

!         copy the matrix a into the l and u structures.
      do 60 icol = 1, n
         dinv(icol) = a(ja(icol))
         jbgn = ja(icol)+1
         jend = ja(icol+1)-1
         if ( jbgn.le.jend ) then
            do 50 j = jbgn, jend
               irow = ia(j)
               if ( irow.lt.icol ) then
!         part of the upper triangle.
                  iu(ncol(icol)) = irow
                  u(ncol(icol)) = a(j)
                  ncol(icol) = ncol(icol) + 1
               else
!         part of the lower triangle (stored by row).
                  jl(nrow(irow)) = icol
                  l(nrow(irow)) = a(j)
                  nrow(irow) = nrow(irow) + 1
                  if ( isym.ne.0 ) then
!         symmetric...copy lower triangle into upper triangle as well.
                     iu(ncol(irow)) = icol
                     u(ncol(irow)) = a(j)
                     ncol(irow) = ncol(irow) + 1
                  endif
               endif
 50         continue
         endif
 60   continue

!         sort the rows of l and the columns of u.
      do 110 k = 2, n
         jbgn = ju(k)
         jend = ju(k+1)-1
         if ( jbgn.lt.jend ) then
            do 80 j = jbgn, jend-1
               do 70 i = j+1, jend
                  if ( iu(j).gt.iu(i) ) then
                     itemp = iu(j)
                     iu(j) = iu(i)
                     iu(i) = itemp
                     temp = u(j)
                     u(j) = u(i)
                     u(i) = temp
                  endif
 70            continue
 80         continue
         endif
         ibgn = il(k)
         iend = il(k+1)-1
         if ( ibgn.lt.iend ) then
            do 100 i = ibgn, iend-1
               do 90 j = i+1, iend
                  if ( jl(i).gt.jl(j) ) then
                     jtemp = ju(i)
                     ju(i) = ju(j)
                     ju(j) = jtemp
                     temp = l(i)
                     l(i) = l(j)
                     l(j) = temp
                  endif
 90            continue
 100        continue
         endif
 110  continue

!         perform the incomplete ldu decomposition.
      do 300 i=2,n

!           i-th row of l
         indx1 = il(i)
         indx2 = il(i+1) - 1
         if (indx1 .gt. indx2) go to 200
         do 190 indx=indx1,indx2
            if (indx .eq. indx1) go to 180
            indxr1 = indx1
            indxr2 = indx - 1
            indxc1 = ju(jl(indx))
            indxc2 = ju(jl(indx)+1) - 1
            if (indxc1 .gt. indxc2) go to 180
 160        kr = jl(indxr1)
 170        kc = iu(indxc1)
            if (kr .gt. kc) then
               indxc1 = indxc1 + 1
               if (indxc1 .le. indxc2) go to 170
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if (indxr1 .le. indxr2) go to 160
            elseif(kr .eq. kc) then
               l(indx) = l(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if (indxr1.le.indxr2 .and. indxc1.le.indxc2) go to 160
            endif
 180        l(indx) = l(indx)/dinv(jl(indx))
 190     continue

!         ith column of u
 200     indx1 = ju(i)
         indx2 = ju(i+1) - 1
         if (indx1 .gt. indx2) go to 260
         do 250 indx=indx1,indx2
            if (indx .eq. indx1) go to 240
            indxc1 = indx1
            indxc2 = indx - 1
            indxr1 = il(iu(indx))
            indxr2 = il(iu(indx)+1) - 1
            if (indxr1 .gt. indxr2) go to 240
 210        kr = jl(indxr1)
 220        kc = iu(indxc1)
            if (kr .gt. kc) then
               indxc1 = indxc1 + 1
               if (indxc1 .le. indxc2) go to 220
            elseif(kr .lt. kc) then
               indxr1 = indxr1 + 1
               if (indxr1 .le. indxr2) go to 210
            elseif(kr .eq. kc) then
               u(indx) = u(indx) - l(indxr1)*dinv(kc)*u(indxc1)
               indxr1 = indxr1 + 1
               indxc1 = indxc1 + 1
               if (indxr1.le.indxr2 .and. indxc1.le.indxc2) go to 210
            endif
 240        u(indx) = u(indx)/dinv(iu(indx))
 250     continue

!         ith diagonal element
 260     indxr1 = il(i)
         indxr2 = il(i+1) - 1
         if (indxr1 .gt. indxr2) go to 300
         indxc1 = ju(i)
         indxc2 = ju(i+1) - 1
         if (indxc1 .gt. indxc2) go to 300
 270     kr = jl(indxr1)
 280     kc = iu(indxc1)
         if (kr .gt. kc) then
            indxc1 = indxc1 + 1
            if (indxc1 .le. indxc2) go to 280
         elseif(kr .lt. kc) then
            indxr1 = indxr1 + 1
            if (indxr1 .le. indxr2) go to 270
         elseif(kr .eq. kc) then
            dinv(i) = dinv(i) - l(indxr1)*dinv(kc)*u(indxc1)
            indxr1 = indxr1 + 1
            indxc1 = indxc1 + 1
            if (indxr1 .le. indxr2 .and. indxc1 .le. indxc2) go to 270
         endif

 300  continue

!         replace diagonal lts by their inverses.
      do 430 i=1,n
         dinv(i) = 1./dinv(i)
 430  continue

      return
      end

      subroutine sgmres(n, b, x, nelt, ia, ja, a, isym,
     $     itol, tol, itmax, iter, err, ierr, iunit, sb, sx,
     $     rgwk, lrgw, igwk, ligw, rwork, iwork )

      implicit none

      integer n, nelt, ia(nelt), ja(nelt), isym, itol, itmax, iter
      integer iunit, lrgw, ligw, igwk(ligw), iwork(*)
      integer jpre, kmp, maxl, nms, maxlp1, nmsl, nrsts, nrmax
      integer i, iflag, lr, ldl, lhes, lgmr, lq, lv, lw
      integer lz, ierr, lzm1, lxl, jscal

      real b(n), x(n), tol, err, sb(n), sx(n), a(nelt)
      real rgwk(lrgw), rwork(*)
      real bnrm, rhol, sum, slap_snrm2

      ierr = 0
      rhol = 0.0
!   ---------------------------------------------------------------
!         load method parameters with user values or defaults.
!   --------------------------------------------------------------
      maxl = igwk(1)
      if (maxl .eq. 0) maxl = 10
      if (maxl .gt. n) maxl = n
      kmp = igwk(2)
      if (kmp .eq. 0) kmp = maxl
      if (kmp .gt. maxl) kmp = maxl
      jscal = igwk(3)
      jpre = igwk(4)
!         check for consistent values of itol and jpre.
      if ( itol.eq.1 .and. jpre.lt.0 ) goto 650
      if ( itol.eq.2 .and. jpre.ge.0 ) goto 650
      nrmax = igwk(5)
      if ( nrmax.eq.0 ) nrmax = 10
!         if nrmax .eq. -1, then set nrmax = 0 to
!                turn off restarting.
      if ( nrmax.eq.-1 ) nrmax = 0

!         initialize counters.
      iter = 0
      nms = 0
      nrsts = 0
!   ---------------------------------------------------------------
!         form work array segment pointers.
!   --------------------------------------------------------------
      maxlp1 = maxl + 1
      lv = 1
      lr = lv + n*maxlp1
      lhes = lr + n + 1
      lq = lhes + maxl*maxlp1
      ldl = lq + 2*maxl
      lw = ldl + n
      lxl = lw + n
      lz = lxl + n

!     load igwk(6) with required minimum length of the rgwk array.
      igwk(6) = lz + n - 1
      if ( lz+n-1.gt.lrgw ) goto 640
!   --------------------------------------------------------------
!         calculate scaled-preconditioned norm of rhs vector b.
!   --------------------------------------------------------------
      if (jpre .lt. 0) then
         call sslui(n, b, rgwk(lr), nelt, ia, ja, a, isym,
     $        rwork, iwork)
         nms = nms + 1
      else
         call slap_scopy(n, b, 1, rgwk(lr), 1)
      endif
      if ( jscal.eq.2 .or. jscal.eq.3 ) then
         sum = 0.e0
         do 10 i = 1,n
            sum = sum + (rgwk(lr-1+i)*sb(i))**2
 10      continue
         bnrm = sqrt(sum)
      else
         bnrm = slap_snrm2(n,rgwk(lr),1)
      endif
!   -------------------------------------------------------------
!         calculate initial residual.
!   -------------------------------------------------------------
      call ssmv(n, x, rgwk(lr), nelt, ia, ja, a, isym)
      do 50 i = 1,n
         rgwk(lr-1+i) = b(i) - rgwk(lr-1+i)
 50   continue
!   -------------------------------------------------------------
!         if performing restarting, then load the residual into the
!         correct location in the rgwk array.
!   --------------------------------------------------------------
 100  continue
      if ( nrsts.gt.nrmax ) goto 610
      if ( nrsts.gt.0 ) then
!    copy the curr residual to different loc in the rgwk array.
         call slap_scopy(n, rgwk(ldl), 1, rgwk(lr), 1)
      endif
!   ------------------------------------------------------------
!         use the spigmr algorithm to solve the linear system a*z = r.
!   -------------------------------------------------------------
      call spigmr(n, rgwk(lr), sb, sx, jscal, maxl, maxlp1, kmp,
     $       nrsts, jpre, nmsl, rgwk(lz), rgwk(lv),
     $       rgwk(lhes), rgwk(lq), lgmr, rwork, iwork, rgwk(lw),
     $       rgwk(ldl), rhol, nrmax, b, bnrm, x, rgwk(lxl), itol,
     $       tol, nelt, ia, ja, a, isym, iunit, iflag, err)
      iter = iter + lgmr
      nms = nms + nmsl

!         increment x by the current approximate solution z of a*z = r.

      lzm1 = lz - 1
      do 110 i = 1,n
         x(i) = x(i) + rgwk(lzm1+i)
 110  continue
      if ( iflag.eq.0 ) goto 600
      if ( iflag.eq.1 ) then
         nrsts = nrsts + 1
         goto 100
      endif
      if ( iflag.eq.2 ) goto 620
!   ------------------------------------------------------------------
!         all returns are made through this section.
!   ------------------------------------------------------------------
!         the iteration has converged.

 600  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 0
      return

!         max number((nrmax+1)*maxl) of linear iterations performed.
 610  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 1
      return

!         gmres failed to reduce last residual in maxl iterations.
!         the iteration has stalled.
 620  continue
      igwk(7) = nms
      rgwk(1) = rhol
      ierr = 2
      return
!         error return.  insufficient length for rgwk array.
 640  continue
      err = tol
      ierr = -1
      return
!         error return.  inconsistent itol and jpre values.
 650  continue
      err = tol
      ierr = -2

      return
      end

      subroutine slap_scopy(n,sx,incx,sy,incy)

!     copies a vector, x, to a vector, y.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.

      implicit none

      integer i, incx, incy, ix, iy, m, mp1, n

      real sx(*), sy(*)

      if (n.le.0)return
      if (incx.eq.1.and.incy.eq.1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

!        code for both increments equal to 1

!        clean-up loop

   20 m = mod(n,7)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if ( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue

      return
      end

      subroutine qs2i1r( ia, ja, a, n, kflag )

!***purpose  sort an integer array also moving an integer and real array
!            this routine sorts the integer array ia and makes the same
!            interchanges in the integer array ja and the real array a.
!            the array ia may be sorted in increasing order or decreas-
!            ing order. a slightly modified quicksort algorithm is used.

!     description of parameters
!        ia - integer array of values to be sorted.
!        ja - integer array to be carried along.
!         a - real array to be carried along.
!         n - number of values in integer array ia to be sorted.
!     kflag - control parameter
!           = 1 means sort ia in increasing order.
!           =-1 means sort ia in decreasing order.

      implicit none

      integer il(21), iu(21)
      integer ia(n), ja(n), it, iit, jt, jjt, k
      integer kflag, ij, l, j, nn, n, m, i, kk

      real a(n), ta, tta, r

      nn = n
      if ( n.eq.1 ) return
      kk = iabs(kflag)

!     alter array ia to get decreasing order if needed.

      if ( kflag.lt.1 ) then
         do 20 i=1,nn
            ia(i) = -ia(i)
 20      continue
      endif

!     sort ia and carry ja and a along.
!     and now...just a little black magic...
      m = 1
      i = 1
      j = nn
      r = .375
 210  if ( r.le.0.5898437 ) then
         r = r + 3.90625e-2
      else
         r = r-.21875
      endif
 225  k = i

!     select a central element of the array and save it in location
!     it, jt, at.

      ij = i + ifix (float (j-i) *r)
      it = ia(ij)
      jt = ja(ij)
      ta = a(ij)

!     if first element of array is greater than it, interchange with it.

      if ( ia(i).gt.it ) then
         ia(ij) = ia(i)
         ia(i)  = it
         it     = ia(ij)
         ja(ij) = ja(i)
         ja(i)  = jt
         jt     = ja(ij)
         a(ij)  = a(i)
         a(i)   = ta
         ta     = a(ij)
      endif
      l=j

!     if last element of array is less than it, swap with it.

      if ( ia(j).lt.it ) then
         ia(ij) = ia(j)
         ia(j)  = it
         it     = ia(ij)
         ja(ij) = ja(j)
         ja(j)  = jt
         jt     = ja(ij)
         a(ij)  = a(j)
         a(j)   = ta
         ta     = a(ij)

!     if first element of array is greater than it, swap with it.

         if ( ia(i).gt.it ) then
            ia(ij) = ia(i)
            ia(i)  = it
            it     = ia(ij)
            ja(ij) = ja(i)
            ja(i)  = jt
            jt     = ja(ij)
            a(ij)  = a(i)
            a(i)   = ta
            ta     = a(ij)
         endif
      endif

!     find an element in the second half of the array which is
!     smaller than it.

  240 l=l-1
      if ( ia(l).gt.it ) go to 240

!     find an element in the first half of the array which is
!     greater than it.

  245 k=k+1
      if ( ia(k).lt.it ) go to 245

!     interchange these elements.

      if ( k.le.l ) then
         iit   = ia(l)
         ia(l) = ia(k)
         ia(k) = iit
         jjt   = ja(l)
         ja(l) = ja(k)
         ja(k) = jjt
         tta   = a(l)
         a(l)  = a(k)
         a(k)  = tta
         goto 240
      endif

!     save upper and lower subscripts of the array yet to be sorted.

      if ( l-i.gt.j-k ) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 260

!     begin again on another portion of the unsorted array.

  255 m = m-1
      if ( m.eq.0 ) go to 300
      i = il(m)
      j = iu(m)
  260 if ( j-i.ge.1 ) go to 225
      if ( i.eq.j ) go to 255
      if ( i.eq.1 ) go to 210
      i = i-1
  265 i = i+1
      if ( i.eq.j ) go to 255
      it = ia(i+1)
      jt = ja(i+1)
      ta =  a(i+1)
      if ( ia(i).le.it ) go to 265
      k=i
  270 ia(k+1) = ia(k)
      ja(k+1) = ja(k)
      a(k+1)  =  a(k)
      k = k-1
      if ( it.lt.ia(k) ) go to 270
      ia(k+1) = it
      ja(k+1) = jt
      a(k+1)  = ta
      go to 265

!     clean up, if necessary.

  300 if ( kflag.lt.1 ) then
         do 310 i=1,nn
            ia(i) = -ia(i)
 310     continue
      endif

      return
      end

      subroutine spigmr(n, r0, sr, sz, jscal, maxl, maxlp1, kmp,
     $     nrsts, jpre, nmsl, z, v, hes, q, lgmr,
     $     rpar, ipar, wk, dl, rhol, nrmax, b, bnrm, x, xl,
     $     itol, tol, nelt, ia, ja, a, isym, iunit, iflag, err)

      implicit none

      integer n, maxl, maxlp1, kmp, jpre, nmsl, lgmr, iflag, jscal
      integer nrsts, nrmax, itol, nelt, isym
      integer ipar(*), ia(*), ja(*)

      real rhol,bnrm,tol
      real r0(*), sr(*), sz(*), z(*), v(n,*), hes(maxlp1,*)
      real q(*), rpar(*), wk(*), dl(*)
      real a(*), b(*), x(*), xl(*)
      real slap_snrm2

!         local variables.

      integer    i, info, ip1, i2, j, k, ll, llp1
      real       r0nrm,c,dlnrm,prod,rho,s,snormw,tem, err
      integer         iter, iunit, itmax, issgmr

!         zero out the z array.
      do 5 i = 1,n
         z(i) = 0.0e0
 5    continue

      iflag = 0
      lgmr = 0
      nmsl = 0
!         load itmax, the maximum number of iterations.
      itmax =(nrmax+1)*maxl
!   -------------------------------------------------------------------
!         the initial residual is the vector r0.
!         apply left precon. if jpre < 0 and this is not a restart.
!         apply scaling to r0 if jscal = 2 or 3.
!   -------------------------------------------------------------------
      if ((jpre .lt. 0) .and.(nrsts .eq. 0)) then
         call slap_scopy(n, r0, 1, wk, 1)
         call sslui(n, wk, r0, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif
      if (((jscal.eq.2) .or.(jscal.eq.3)) .and.(nrsts.eq.0)) then
         do 10 i = 1,n
            v(i,1) = r0(i)*sr(i)
 10      continue
      else
         do 20 i = 1,n
            v(i,1) = r0(i)
 20      continue
      endif
      r0nrm = slap_snrm2(n, v, 1)
      iter = nrsts*maxl

!         call stopping routine issgmr.

      if (issgmr(n, b, x, xl, nelt, ia, ja, a, isym,
     $    nmsl, itol, tol, itmax, iter, err, iunit, v(1,1), z, wk,
     $    rpar, ipar, r0nrm, bnrm, sr, sz, jscal,
     $    kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $    hes, jpre) .ne. 0) return
      tem = 1.0e0/(r0nrm + 1.e-20)
      call slap_sscal(n, tem, v(1,1), 1)

!         zero out the hes array.

      do 50 j = 1,maxl
         do 40 i = 1,maxlp1
            hes(i,j) = 0.0e0
 40      continue
 50   continue
!   -------------------------------------------------------------------
!         main loop to compute the vectors v(*,2) to v(*,maxl).
!         the running product prod is needed for the convergence test.
!   -------------------------------------------------------------------
      prod = 1.0e0
      do 90 ll = 1,maxl
         lgmr = ll
!   -------------------------------------------------------------------
!        unscale  the  current v(ll)  and store  in wk.  call routine
!        sslui    to   compute(m-inverse)*wk,   where    m   is  the
!        preconditioner matrix.  save the answer in z.   call routine
!        ssmv to compute  vnew  = a*z,  where  a is  the the system
!        matrix.  save the answer in  v(ll+1).  scale v(ll+1).   call
!        routine sorth  to  orthogonalize the    new vector vnew   =
!        v(*,ll+1).  call routine sheqr to update the factors of hes.
!   -------------------------------------------------------------------
        if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
           do 60 i = 1,n
              wk(i) = v(i,ll)/sz(i)
 60        continue
        else
           call slap_scopy(n, v(1,ll), 1, wk, 1)
        endif
        if (jpre .gt. 0) then
           call sslui(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
           nmsl = nmsl + 1
           call ssmv(n, z, v(1,ll+1), nelt, ia, ja, a, isym)
        else
           call ssmv(n, wk, v(1,ll+1), nelt, ia, ja, a, isym)
        endif
        if (jpre .lt. 0) then
           call slap_scopy(n, v(1,ll+1), 1, wk, 1)
           call sslui(n,wk,v(1,ll+1),nelt,ia,ja,a,isym,rpar,ipar)
           nmsl = nmsl + 1
        endif
        if ((jscal .eq. 2) .or.(jscal .eq. 3)) then
           do 65 i = 1,n
              v(i,ll+1) = v(i,ll+1)*sr(i)
 65        continue
        endif
        call sorth(v(1,ll+1), v, hes, n, ll, maxlp1, kmp, snormw)
        hes(ll+1,ll) = snormw
        call sheqr(hes, maxlp1, ll, q, info, ll)
        if (info .eq. ll) go to 120
!   -------------------------------------------------------------------
!         update rho, the estimate of the norm of the residual r0-a*zl.
!         if kmp <  maxl, then the vectors v(*,1),...,v(*,ll+1) are not
!         necessarily orthogonal for ll > kmp.  the vector dl must then
!         be computed, and its norm used in the calculation of rho.
!   -------------------------------------------------------------------
        prod = prod*q(2*ll)
        rho = abs(prod*r0nrm)
        if ((ll.gt.kmp) .and.(kmp.lt.maxl)) then
           if (ll .eq. kmp+1) then
              call slap_scopy(n, v(1,1), 1, dl, 1)
              do 75 i = 1,kmp
                 ip1 = i + 1
                 i2 = i*2
                 s = q(i2)
                 c = q(i2-1)
                 do 70 k = 1,n
                    dl(k) = s*dl(k) + c*v(k,ip1)
 70              continue
 75           continue
           endif
           s = q(2*ll)
           c = q(2*ll-1)/snormw
           llp1 = ll + 1
           do 80 k = 1,n
              dl(k) = s*dl(k) + c*v(k,llp1)
 80        continue
           dlnrm = slap_snrm2(n, dl, 1)
           rho = rho*dlnrm
        endif
        rhol = rho
!   -------------------------------------------------------------------
!         test for convergence.  if passed, compute approximation zl.
!         if failed and ll < maxl, then continue iterating.
!   -------------------------------------------------------------------
        iter = nrsts*maxl + lgmr
        if (issgmr(n, b, x, xl, nelt, ia, ja, a, isym,
     $      nmsl, itol, tol, itmax, iter, err, iunit, dl, z, wk,
     $      rpar, ipar, rhol, bnrm, sr, sz, jscal,
     $      kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $      hes, jpre) .ne. 0) go to 200
        if (ll .eq. maxl) go to 100
!   -------------------------------------------------------------------
!         rescale so that the norm of v(1,ll+1) is one.
!   -------------------------------------------------------------------
        tem = 1.0e0/snormw
        call slap_sscal(n, tem, v(1,ll+1), 1)
 90   continue
 100  continue
      if (rho .lt. r0nrm) go to 150
 120  continue
      iflag = 2

!         load approximate solution with zero.

      do 130 i = 1,n
         z(i) = 0.e0
 130  continue
      return
 150  iflag = 1

!         tolerance not met, but residual norm reduced.

      if (nrmax .gt. 0) then

!        if performing restarting (nrmax > 0)  calculate the residual
!        vector rl and  store it in the dl  array.  if the incomplete
!        version is being used (kmp < maxl) then dl has  already been
!        calculated up to a scaling factor.   use srlcal to calculate
!        the scaled residual vector.

         call srlcal(n, kmp, maxl, maxl, v, q, dl, snormw, prod,
     $        r0nrm)
      endif
!   -------------------------------------------------------------------
!         compute the approximation zl to the solution.  since the
!         vector z was used as work space, and the initial guess
!         of the linear iteration is zero, z must be reset to zero.
!   -------------------------------------------------------------------
 200  continue
      ll = lgmr
      llp1 = ll + 1
      do 210 k = 1,llp1
         r0(k) = 0.0e0
 210  continue
      r0(1) = r0nrm
      call shels(hes, maxlp1, ll, q, r0)
      do 220 k = 1,n
         z(k) = 0.0e0
 220  continue
      do 230 i = 1,ll
         call slap_saxpy(n, r0(i), v(1,i), 1, z, 1)
 230  continue
      if ((jscal .eq. 1) .or.(jscal .eq. 3)) then
         do 240 i = 1,n
            z(i) = z(i)/sz(i)
 240     continue
      endif
      if (jpre .gt. 0) then
         call slap_scopy(n, z, 1, wk, 1)
         call sslui(n, wk, z, nelt, ia, ja, a, isym, rpar, ipar)
         nmsl = nmsl + 1
      endif

      return
      end

      subroutine slap_sscal(n,sa,sx,incx)

!     scales a vector by a constant.
!     uses unrolled loops for increment equal to 1.
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increments, 9/29/88.

      implicit none

      real sa, sx(*)

      integer i, ix, incx, m, mp1, n

      if (n.le.0)return
      if (incx.eq.1)go to 20

!        code for increment not equal to 1

      ix = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        sx(ix) = sa*sx(ix)
        ix = ix + incx
   10 continue
      return

!        code for increment equal to 1

!        clean-up loop

   20 m = mod(n,5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if ( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue

      return
      end

      subroutine sorth(vnew, v, hes, n, ll, ldhes, kmp, snormw)

      implicit none

      integer n, ll, ldhes, kmp

      real vnew(*), v(n,*), hes(ldhes,*), snormw
      real slap_snrm2, slap_sdot

!         internal variables.

      integer i, i0

      real arg, sumdsq, tem, vnrm

!         get norm of unaltered vnew for later use.
      vnrm = slap_snrm2(n, vnew, 1)
!   -------------------------------------------------------------------
!         perform the modified gram-schmidt procedure on vnew =a*v(ll).
!         scaled inner products give new column of hes.
!         projections of earlier vectors are subtracted from vnew.
!   -------------------------------------------------------------------
      i0 = max0(1,ll-kmp+1)
      do 10 i = i0,ll
         hes(i,ll) = slap_sdot(n, v(1,i), 1, vnew, 1)
         tem = -hes(i,ll)
         call slap_saxpy(n, tem, v(1,i), 1, vnew, 1)
 10   continue
!   -------------------------------------------------------------------
!         compute snormw = norm of vnew.  if vnew is small compared
!         to its input value (in norm), then reorthogonalize vnew to
!         v(*,1) through v(*,ll).  correct if relative correction
!         exceeds 1000*(unit roundoff).  finally, correct snormw using
!         the dot products involved.
!   -------------------------------------------------------------------
      snormw = slap_snrm2(n, vnew, 1)
      if (vnrm + 0.001e0*snormw .ne. vnrm) return
      sumdsq = 0.0e0
      do 30 i = i0,ll
         tem = -slap_sdot(n, v(1,i), 1, vnew, 1)
         if (hes(i,ll) + 0.001e0*tem .eq. hes(i,ll)) go to 30
         hes(i,ll) = hes(i,ll) - tem
         call slap_saxpy(n, tem, v(1,i), 1, vnew, 1)
         sumdsq = sumdsq + tem**2
 30   continue
      if (sumdsq .eq. 0.0e0) return
      arg = amax1(0.0e0,snormw**2 - sumdsq)
      snormw = sqrt(arg)

      return
      end

      subroutine sheqr(a, lda, n, q, info, ijob)

      implicit none

      integer lda, n, info, ijob

      real a(lda,*), q(*)

!         local variables.

      integer i, iq, j, k, km1, kp1, nm1

      real c, s, t, t1, t2

      if (ijob .gt. 1) go to 70
!   -------------------------------------------------------------------
!         a new factorization is desired.
!   -------------------------------------------------------------------
!         qr decomposition without pivoting.

      info = 0
      do 60 k = 1, n
         km1 = k - 1
         kp1 = k + 1

!           compute k-th column of r.
!           first, multiply the k-th column of a by the previous
!           k-1 givens rotations.

         if (km1 .lt. 1) go to 20
         do 10 j = 1, km1
            i = 2*(j-1) + 1
            t1 = a(j,k)
            t2 = a(j+1,k)
            c = q(i)
            s = q(i+1)
            a(j,k) = c*t1 - s*t2
            a(j+1,k) = s*t1 + c*t2
 10      continue

!         compute givens components c and s.

 20      continue
         iq = 2*km1 + 1
         t1 = a(k,k)
         t2 = a(kp1,k)
         if ( t2.eq.0.0e0 ) then
            c = 1.0e0
            s = 0.0e0
         elseif( abs(t2).ge.abs(t1) ) then
            t = t1/t2
            s = -1.0e0/sqrt(1.0e0+t*t)
            c = -s*t
         else
            t = t2/t1
            c = 1.0e0/sqrt(1.0e0+t*t)
            s = -c*t
         endif
         q(iq) = c
         q(iq+1) = s
         a(k,k) = c*t1 - s*t2
         if ( a(k,k).eq.0.0e0 ) info = k
 60   continue
      return
!   -------------------------------------------------------------------
!         the old factorization of a will be updated.  a row and a
!         column has been added to the matrix a.  n by n-1 is now
!         the old size of the matrix.
!   -------------------------------------------------------------------
 70   continue
      nm1 = n - 1
!   -------------------------------------------------------------------
!         multiply the new column by the n previous givens rotations.
!   -------------------------------------------------------------------
      do 100 k = 1,nm1
         i = 2*(k-1) + 1
         t1 = a(k,n)
         t2 = a(k+1,n)
         c = q(i)
         s = q(i+1)
         a(k,n) = c*t1 - s*t2
         a(k+1,n) = s*t1 + c*t2
 100  continue
!   -------------------------------------------------------------------
!         complete update of decomposition by forming last givens
!         rotation, and multiplying it times the column
!         vector(a(n,n),a(np1,n)).
!   -------------------------------------------------------------------
      info = 0
      t1 = a(n,n)
      t2 = a(n+1,n)
      if ( t2.eq.0.0e0 ) then
         c = 1.0e0
         s = 0.0e0
      elseif( abs(t2).ge.abs(t1) ) then
         t = t1/t2
         s = -1.0e0/sqrt(1.0e0+t*t)
         c = -s*t
      else
         t = t2/t1
         c = 1.0e0/sqrt(1.0e0+t*t)
         s = -c*t
      endif
      iq = 2*n - 1
      q(iq) = c
      q(iq+1) = s
      a(n,n) = c*t1 - s*t2
      if (a(n,n) .eq. 0.0e0) info = n

      return
      end

       subroutine srlcal(n, kmp, ll, maxl, v, q, rl, snormw, prod,
     $     r0nrm)

      implicit none

      integer n, kmp, ll, maxl

      real v(n,*), q(*), rl(n), snormw

!         internal variables.

      integer i, ip1, i2, k, llp1, llm1

      real c, s, tem, prod, r0nrm

      if (kmp .eq. maxl) then

!         calculate rl.  start by copying v(*,1) into rl.

         call slap_scopy(n, v(1,1), 1, rl, 1)
         llm1 = ll - 1
         do 20 i = 1,llm1
            ip1 = i + 1
            i2 = i*2
            s = q(i2)
            c = q(i2-1)
            do 10 k = 1,n
               rl(k) = s*rl(k) + c*v(k,ip1)
 10         continue
 20      continue
         s = q(2*ll)
         c = q(2*ll-1)/snormw
         llp1 = ll + 1
         do 30 k = 1,n
            rl(k) = s*rl(k) + c*v(k,llp1)
 30      continue
      endif

!         when kmp < maxl, rl vector already partially calculated.
!         scale rl by r0nrm*prod to obtain the residual rl.

      tem = r0nrm*prod
      call slap_sscal(n, tem, rl, 1)

      return
      end

      subroutine shels(a, lda, n, q, b)

      implicit none

      integer lda, n

      real a(lda,*), b(*), q(*)

!         local variables.

      integer iq, k, kb, kp1
      real c, s, t, t1, t2

!         minimize(b-a*x,b-a*x).  first form q*b.

      do 20 k = 1, n
         kp1 = k + 1
         iq = 2*(k-1) + 1
         c = q(iq)
         s = q(iq+1)
         t1 = b(k)
         t2 = b(kp1)
         b(k) = c*t1 - s*t2
         b(kp1) = s*t1 + c*t2
 20   continue

!         now solve  r*x = q*b.

      do 40 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call slap_saxpy(k-1, t, a(1,k), 1, b(1), 1)
 40   continue

      return
      end

      subroutine slap_saxpy(n,sa,sx,incx,sy,incy)

!     constant times a vector plus a vector.
!     uses unrolled loop for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

      implicit none

      integer i, incx, incy, ix, iy, m, mp1, n

      real sx(*), sy(*), sa

      if (n.le.0)return
      if (sa .eq. 0.0) return
      if (incx.eq.1.and.incy.eq.1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return

!        code for both increments equal to 1

!        clean-up loop

   20 m = mod(n,4)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if ( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue

      return
      end

      subroutine ssmv( n, x, y, nelt, ia, ja, a, isym )

!***purpose  slap column format sparse matrix vector product.
!            routine to calculate the sparse matrix vector product:
!            y = a*x.
! *arguments:
! n      :in       integer.
!         order of the matrix.
! x      :in       real x(n).
!         the vector that should be multiplied by the matrix.
! y      :out      real y(n).
!         the product of the matrix and the vector.
! nelt   :in       integer.
!         number of non-zeros stored in a.
! ia     :in       integer ia(nelt).
! ja     :in       integer ja(n+1).
! a      :in       integer a(nelt).
!         these arrays should hold the matrix a in the slap column
!         format.  see "description", below.
! isym   :in       integer.
!         flag to indicate symmetric storage format.
!         if isym=0, all nonzero entries of the matrix are stored.
!         if isym=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.

! *cautions:
!     this   routine   assumes  that  the matrix a is stored in slap
!     column format.  it does not check  for  this (for  speed)  and
!     evil, ugly, ornery and nasty things  will happen if the matrix
!     data  structure  is,  in fact, not slap column.  beware of the
!     wrong data structure

      implicit none

      integer n, nelt, ia(nelt), ja(nelt), isym
      integer i, ibgn, iend, j, jbgn, jend, irow, icol

      real a(nelt), x(n), y(n)

!         zero out the result vector.
      do 10 i = 1, n
         y(i) = 0.0
 10   continue

!         multiply by a.

      do 30 icol = 1, n
         ibgn = ja(icol)
         iend = ja(icol+1)-1
         do 20 i = ibgn, iend
            y(ia(i)) = y(ia(i)) + a(i)*x(icol)
 20      continue
 30   continue

      if ( isym.eq.1 ) then

!         the matrix is non-symmetric.  need to get the other half in...
!         this loops assumes that the diagonal is the first entry in
!         each column.

         do 50 irow = 1, n
            jbgn = ja(irow)+1
            jend = ja(irow+1)-1
            if ( jbgn.gt.jend ) goto 50
            do 40 j = jbgn, jend
               y(irow) = y(irow) + a(j)*x(ia(j))
 40         continue
 50      continue
      endif

      return
      end

      real function slap_sdot(n,sx,incx,sy,incy)

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

      implicit none

      real sx(*),sy(*),stemp

      integer i,incx,incy,ix,iy,m,mp1,n

      stemp = 0.0e0
      slap_sdot = 0.0e0
      if (n.le.0)return
      if (incx.eq.1.and.incy.eq.1)go to 20

!        code for unequal increments or equal increments
!          not equal to 1

      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      slap_sdot = stemp
      return

!        code for both increments equal to 1

!        clean-up loop

   20 m = mod(n,5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if ( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 slap_sdot = stemp

      return
      end

      real function slap_snrm2 ( n, sx, incx)

      implicit none

      integer next, i, j, n, nn, incx

      real sx(*),  cutlo, cuthi, hitest, sum, xmax, zero, one

      data zero, one /0.0e0, 1.0e0/

!     euclidean norm of the n-vector stored in sx() with storage
!     increment incx .
!     if    n .le. 0 return with result = 0.
!     if n .ge. 1 then incx must be .ge. 1

!     four phase method     using two built-in constants that are
!     hopefully applicable to all machines.
!         cutlo = maximum of  sqrt(u/eps)  over all known machines.
!         cuthi = minimum of  sqrt(v)      over all known machines.
!     where
!         eps = smallest no. such that eps + 1. .gt. 1.
!         u   = smallest positive no.   (underflow limit)
!         v   = largest  no.            (overflow  limit)

!     brief outline of algorithm..

!     phase 1    scans zero components.
!     move to phase 2 when a component is nonzero and .le. cutlo
!     move to phase 3 when a component is .gt. cutlo
!     move to phase 4 when a component is .ge. cuthi/m
!     where m = n for x() real and m = 2*n for complex.

      data cutlo, cuthi / 4.441e-16,  1.304e19 /

      if (n .gt. 0) go to 10
         slap_snrm2  = zero
         go to 300

   10 next = 30
      sum = zero
      nn = n * incx
!                                                 begin main loop
      i = 1
   20 if (next.eq.30) goto 30
      if (next.eq.50) goto 50
      if (next.eq.70) goto 70
      if (next.eq.110) goto 110

   30 if ( abs(sx(i)) .gt. cutlo) go to 85
      next = 50
      xmax = zero

!                        phase 1.  sum is zero

   50 if ( sx(i) .eq. zero) go to 200
      if ( abs(sx(i)) .gt. cutlo) go to 85

!                                prepare for phase 2.
      next = 70
      go to 105

!                                prepare for phase 4.

  100 i = j
      next = 110
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115

!                   phase 2.  sum is small.
!                             scale to avoid destructive underflow.

   70 if ( abs(sx(i)) .gt. cutlo ) go to 75

!                     common code for phases 2 and 4.
!                     in phase 4 sum is large.  scale to avoid overflow.

  110 if ( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200

  115 sum = sum + (sx(i)/xmax)**2
      go to 200

!                  prepare for phase 3.

   75 sum = (sum * xmax) * xmax

!     for real or d.p. set hitest = cuthi/n
!     for complex      set hitest = cuthi/(2*n)

   85 hitest = cuthi/float( n )

!                   phase 3.  sum is mid-range.  no scaling.

      do 95 j =i,nn,incx
        if (abs(sx(j)) .ge. hitest) go to 100
        sum = sum + sx(j)**2
   95 continue
      slap_snrm2 = sqrt( sum )
      go to 300

  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20

!              end of main loop.

!              compute square root and adjust for scaling.

      slap_snrm2 = xmax * sqrt(sum)
  300 continue

      return
      end

       function issgmr(n, b, x, xl, nelt, ia, ja, a, isym,
     $     nmsl, itol, tol, itmax, iter, err, iunit, r, z, dz,
     $     rwork, iwork, rnrm, bnrm, sb, sx, jscal,
     $     kmp, lgmr, maxl, maxlp1, v, q, snormw, prod, r0nrm,
     $     hes, jpre)

      implicit none

      integer kmp, lgmr, maxl, maxlp1, jpre, nmsl, n, nelt
      integer   ia(nelt), ja(nelt), iwork(*)

      real dxnrm, rnrm, r0nrm, snormw, solnrm, prod
      real slap_snrm2, sb(*), sx(*), q(*), v(n,*)
      real b(*), x(*), r(*), z(*), dz(*), rwork(*)
      real hes(maxlp1,*), xl(*), a(nelt)

!         local variables.

      integer i, itol, iter, itmax, issgmr, iunit, ielmax
      integer jscal, isym

      real tol, err, bnrm, fuzz, ratmax, rat
      real c,dlnrm,rho,s,tem

!      common /solblk/ soln(1)
!      save solnrm

      issgmr = 0
      if ( itol.eq.0 ) then

!       use input from spigmr to determine if stop conditions are met.

         err = rnrm/(bnrm + 1.e-20)
      endif
      if ( (itol.gt.0) .and. (itol.le.3) ) then

!       use srlcal to calculate the scaled residual vector.
!       store answer in r.

         if ( lgmr.ne.0 ) call srlcal(n, kmp, lgmr, maxl, v, q, r,
     $                                snormw, prod, r0nrm)
         if ( itol.le.2 ) then
!         err = ||residual||/||righthandside||(2-norms).
            err = slap_snrm2(n, r, 1)/(bnrm + 1.e-20)

!         unscale r by r0nrm*prod when kmp < maxl.

            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call slap_sscal(n, tem, r, 1)
            endif
         elseif ( itol.eq.3 ) then
!         err = max |(minv*residual)(i)/x(i)|
!         when jpre .lt. 0, r already contains minv*residual.
            if ( jpre.gt.0 ) then
               call sslui(n, r, dz, nelt, ia, ja, a, isym, rwork,
     $              iwork)
               nmsl = nmsl + 1
            endif

!         unscale r by r0nrm*prod when kmp < maxl.

            if ( (kmp.lt.maxl) .and. (lgmr.ne.0) ) then
               tem = 1.0e0/(r0nrm*prod)
               call slap_sscal(n, tem, r, 1)
            endif

            fuzz = 1.0e-20
            ielmax = 1
            ratmax = abs(dz(1))/amax1(abs(x(1)),fuzz)
            do 25 i = 2, n
               rat = abs(dz(i))/amax1(abs(x(i)),fuzz)
               if ( rat.gt.ratmax ) then
                  ielmax = i
                  ratmax = rat
               endif
 25         continue
            err = ratmax
            if ( ratmax.le.tol ) issgmr = 1
            if ( iunit.gt.0 ) write(iunit,1020) iter, ielmax, ratmax
            return
         endif
      endif

      if ( iunit.ne.0 ) then
         if ( iter.eq.0 ) then
            write(iunit,1000) maxl, kmp, n, itol
         endif
!         write(6,*) 'itol, iter =', itol, iter
!         write(6,*) 'rnrm, bnrm =', rnrm, bnrm
         write(iunit,1010) iter, rnrm/(bnrm + 1.e-20), err
      endif
      if ( err.le.tol ) issgmr = 1

 1000 format(' generalized minimum residual(',i3,i3,') for ',
     $     'n, itol = ',i6, i6,
     $     /' iter','   natral err est','   error estimate')
 1010 format(1x,i4,1x,e16.7,1x,e16.7)
 1020 format(1x,' iter = ',i5, ' ielmax = ',i5,
     $     ' |r(ielmax)/x(ielmax)| = ',e12.5)

#endif

      return
      end
