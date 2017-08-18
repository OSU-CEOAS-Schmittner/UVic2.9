      subroutine eqstate (zt, km, ro0_den, to_den, so_den, c_den
     &,                   tmink, tmaxk, smink, smaxk)

#if defined O_mom
!=======================================================================

!             E Q U A T I O N   O F   S T A T E   M O D U L E

!     Calculate polynomial coefficients for density computations in MOM.

!     To generate the coefficients:

!        1) set up the grid in "grids.F" module

!        2) compile and run this module by setting the options in and
!           executing the "run_denscoef" script

!     To install the coefficients in MOM:

!        3) follow the directions at the end of the output from 2)

!     This program calculates the 9 coefficients of a third order
!     polynomial approximation to the equation of state for sea water.
!       The program yields coefficients that will compute density as a
!     function of temperature, and salinity, at predetermined depths,
!     as used in the MOM subroutine "state".
!     More specifically, the densities calculated from the polynomial
!     formula are in the form of sigma anomalies.  The method is
!     taken from that described by Bryan & Cox (1972).
!       By default, the program uses the equation of state set by the
!     Joint Panel on Oceanographic Tables & Standards (UNESCO, 1981)
!     an described by Gill (1982).  An option exists to use the older
!     Knudsen-Ekman equation of state, as described by Fofonoff (1962),
!     if the user prefers.
!       Subroutine "lsqsl2" performs the iterative least-squares
!     polynomial fitting for the overdetermined system.  The algorithm
!     is outlined by Hanson and Lawson (1969), and the code looks as if
!     it has not be touched since that time.

!     references:
!        Bryan, K. & M. Cox, An approximate equation of state
!          for numerical models of ocean circulation, J. Phys.
!          Oceanogr., 2, 510-514, 1972.
!        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
!          New York, 1962, pp 3-30.
!        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical
!          Series No. 30. Academic Press, London, 1982, pp 599-600.
!        Hanson, R., & C. Lawson, Extensions and applications of the
!          Householder algorithm for solving linear least squares
!          problems. Math. Comput., 23, 787-812, 1969.
!        UNESCO, 10th report of the joint panel on oceanographic tables
!          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,
!          Paris, 1981.

!    ifdef options:
!       "knudsen"
!       To over-ride the default of using the UNESCO equation of state
!     and to instead employ the Knudsen-Ekman formula.
!     The default assumption is that potential temperatures will be used
!     (asin the ocean model code).
!-----------------------------------------------------------------------

      implicit none

      include "stdunits.h"
      character(10) :: fname

      integer kmax, kx, kxx, kk, ksdim, krdim, km, k, io, i, j, ka
      integer ndim, nrow, ncol, in, itmax, it, ieq, irank, nn, nhdim
      parameter (kmax=200)
      parameter (kx = 5, kxx = 2*kx, kk = kx*kxx )
      parameter (ksdim = kk+72, krdim = kk+36 )

      real (kind=8) c0, c1, c2, cmtocm, fkx, fi, t1, s1, tot, th1
      real (kind=8) fkk, d, s, t, densit, theta, tot1, tanom, sanom
      real (kind=8) eps, enorm

      real (kind=8) a(kk,9), sigma(kk), sigman(kk), c(kk,9), x(9)
      real (kind=8) sb(ksdim), r(krdim)
      real (kind=8) tmin(kmax), smin(kmax), tmax(kmax), smax(kmax)
      real (kind=8) z(kmax), dd(kmax), ss(kmax), ab(13,kmax), ts(33,4)
      real (kind=8) ta(kxx), sa(kxx), tp(kk), sp(kk), th(kk)

      real ro0_den(km), to_den(km), so_den(km), c_den(km,9), zt(km)
      real tmink(km), tmaxk(km), smink(km), smaxk(km), realz

      data fname /'dncoef.new'/

!  enter bounds for polynomial fit: at 33 levels from sfc to 8000 m.
!           ts(k,1)=lower bnd of t at z=(k-1)*250 meters
!           ts(k,2)=upper bnd of t
!           ts(k,3)=lower bnd of s
!           ts(k,4)=upper bnd of s
!  The user should review the appropriateness of the "ts" values set
!  below, and modify them if the intended modelling application could
!  be expected to yield temperature and salinity values outside of the
!  "ts" ranges set by default.

      data (ts(k,1),k=1,33) / 4*-2.0, 15*-1.0, 14*0.0 /
      data (ts(k,2),k=1,33) / 29.0, 19.0, 14.0, 11.0, 9.0, 28*7.0 /
      data (ts(k,3),k=1,33) / 28.5, 33.7, 34.0, 34.1, 34.2, 34.4,
     &                        2*34.5, 15*34.6, 10*34.7 /
      data (ts(k,4),k=1,33) / 37.0, 36.6, 35.8, 35.7, 35.3, 2*35.1,
     &                        26*35.0 /

!     z       = model levels (midpoint of model layers)
!     tmin, tmax, smin, smax = minimum and maximum in situ temperature
!               and salinity values which define the ranges to be used
!               when computing the polynomials at each model level
!     dd, ds  = increment between temperature and salinity values at
!               each model level to be used in constructing array of
!               temperature, salinity and density for curve fitting
!     ta, sa  = in situ temperature and salinity values available for
!               constructing array of data for curve fitting at each
!               model level
!     tp, sp  = in situ temperature and salinity values constructed from
!               all combinations of ta & sa
!     th      = potential temperature values associated with "tp" at a
!               given level and salinity
!     t1, s1, tot1, th1 = level mean insitu temp., salinity, density,
!               and potential temp. used in polynomial fitting
!     tot     = density (in sigma units) calculate from t1 and s1 at a
!               given model level
!     sigma   = insitu densities (in sigma units) calculated from "tp"
!               and "sp" values
!     sigman  = insitu density anomalies at a given level (formed by
!               subtracting "tot" from sigma)
!     tanom, sanom = temperature and salinity anomalies used in loading
!               array "a" for use in lsqsl2 curve fitting
!     x       = the 9 polynomial coefficients
!     r, sb   = used only in lsqsl2

!=======================================================================

      if (km .gt. kmax) then
        write (stdout,*) '=>Error: increase "kmax" > ',km,' in eqstate'
        stop
      endif

!     set some constants

      c0 = 0.0
      c1 = 1.0
      c2 = 2.0

!     construct depths (meters) from surface to midpoint of levels

      cmtocm = 1.0d-2
      do k=1,km
        z(k) = zt(k) * cmtocm
        if (z(k) .gt. 8000.0) then
          write (stdout,*) '=>Error:depth can`t exceed 8000m in eqstate'
          stop
        endif
      enddo

!     set the temperature and salinity ranges to be used for each
!     model level when performing the polynomial fitting

      do k=1,km
        realz = z(k)/250.0
        i = ifix (realz) + 1
        tmin(k) = ts(i,1)
        tmax(k) = ts(i,2)
        smin(k) = ts(i,3)
        smax(k) = ts(i,4)
      enddo

!  set temperature and salinity increments to be used in creating
!  curve fitting array at each level (twice as many temperature values
!  than salinity values)

      fkx = kx
      do k=1,km
        dd(k) = (tmax(k)-tmin(k)) / (c2*fkx-c1)
        ss(k) = (smax(k)-smin(k)) / (fkx-c1)
      enddo

!     loop over all model levels

      do k=1,km

        do i=1,kxx
          fi = i
          ta(i) = tmin(k) + (fi-c1)*dd(k)
          sa(i) = smin(k) + (fi-c1)*ss(k)
        enddo

!       load the "kxx" combinations of the 2*"kx" insitu temp. and "kx"
!       salinity values into "tp" and "sp"

        do i=1,kxx
          do j=1,kx
            ka = kx*i + j - kx
            tp(ka) = ta(i)
            sp(ka) = sa(j)
          enddo
        enddo

        t1  = c0
        s1  = c0
        tot = c0
        th1 = c0
        fkk = kk

!       calculate insitu density "sigma" for each t,s combination at
!       this depth "d"

        do ka=1,kk
          d = z(k)
          s = sp(ka)
          t = tp(ka)

# if defined O_knudsen
!         "knuekm" returns density (in sigma units) from insitu
!         temperature, salinity, & depth (pressure) using the
!         Knudsen-Ekman formula

          call knuekm(t,s,d,densit)

          sigma(ka) = densit
# else
!         "unesco" returns density (kg per m**3) from insitu
!         temperature, salinity, & depth (pressure) using the UNESCO
!         equation of state

          call unesco(t,s,d,densit)

          sigma(ka) = densit - 1.0d3 + 2.5d-2
# endif

!         "potem" returns potential temp. from from insitu temperature,
!         salinity, & depth (pressure)

          call potem(t,s,d,theta)

          th(ka) = theta
          t1 = t1 + tp(ka)
          s1 = s1 + sp(ka)
          tot = tot + sigma(ka)
          th1 = th1 + th(ka)
        enddo

!       form layer averages "t1", "s1", "th1", and "tot1", and compute
!       reference density "tot" from "t1" and "s1" at this depth "d"

        t1 = t1/fkk
        s1 = s1/fkk
        th1 = th1/fkk
        tot1 = tot/fkk
# if defined O_knudsen

!       "knuekm" returns density from insitu temp., salinity, & depth
!       (pressure) using the Knudsen-Ekman formula

        call knuekm (t1, s1, d, densit)

        tot = densit
# else

!       "unesco" returns density from insitu temp., salinity, & depth
!       (pressure) using the UNESCO equation of state

        call unesco (t1, s1, d, densit)
        tot = densit - 1.0d3 + 2.5d-2

# endif

!       define insitu
        t1 = th1

!       begin loading "ab" array with level averages

        ab(1,k) = z(k)
        ab(2,k) = tot
        ab(3,k) = t1
        ab(4,k) = s1

        do ka=1,kk

!         define insitu
          tp(ka) = th(ka)

!         create anomalies for temperature, salinity & density and
!         load work array "a" with the anomalies and their products

          tanom = tp(ka) - t1
          sanom = sp(ka) - s1
          sigman(ka) = sigma(ka) - tot
          a(ka,1) = tanom
          a(ka,2) = sanom
          a(ka,3) = tanom * tanom
          a(ka,4) = tanom * sanom
          a(ka,5) = sanom * sanom
          a(ka,6) = a(ka,3) * tanom
          a(ka,7) = a(ka,5) * tanom
          a(ka,8) = a(ka,3) * sanom
          a(ka,9) = a(ka,5) * sanom
        enddo

!       set the arguments used in call to "lsqsl2"
!       ndim = first dimension of array a
!       nrow =number of rows of array a
!       ncol = number of columns of array a
!       in = option number of lsqsl2
!       itmax = number of iterations

        ndim = 50
        nrow = kk
        ncol = 9
        in = 1
        itmax = 4

        it = 0
        ieq = 2
        irank = 0
        eps = 1.0e-7
        nhdim = 9

!       LSQL2 is  a Jet Propulsion Laboratory subroutine that
!       computes the least squares fit in an iterative manner for
!       overdetermined systems.

        call lsqsl2 (ndim, a, nrow, ncol, sigman, x, irank, in, itmax,
     &               it, ieq, enorm, eps, nhdim, c, r, sb)

        do i=1,ncol
          ab(i+4,k) = x(i)
        enddo

!       end loop on model levels

      enddo

      nn = ncol + 4

      do k=1,km
        ab(2,k)  = 1.e-3 * ab(2,k)
        ab(4,k)  = 1.e-3 * ab(4,k) - 0.035
        ab(5,k)  = 1.e-3 * ab(5,k)
        ab(7,k)  = 1.e-3 * ab(7,k)
        ab(10,k) = 1.e-3 * ab(10,k)
        ab( 9,k) = 1.e+3 * ab( 9,k)
        ab(11,k) = 1.e+3 * ab(11,k)
        ab(13,k) = 1.e+6 * ab(13,k)
      enddo

!     save this data for a model execution converting insitu
!     temperatures to potential temperatures

      do k=1,km
        ro0_den(k) = ab(2,k)
        to_den(k)  = ab(3,k)
        so_den(k)  = ab(4,k)
        d = z(k)
        s = smin(k)
        t = tmin(k)
        call potem(t,s,d,theta)
        tmink(k)   = theta
        s = smax(k)
        t = tmax(k)
        call potem(t,s,d,theta)
        tmaxk(k)   = theta
        smink(k)   = smin(k)
        smaxk(k)   = smax(k)
        do i=5,13
          c_den(k,i-4) = ab(i,k)
        enddo
      enddo

# if defined O_linearized_density

!     Density a linear function of temp. No pressure effects.
!     Take following numbers from Appendix 3 of Gill (1982)
!     with T=10C, S=35psu, p=0atm.
!     model salinity is S(model) = (S(psu)-35)/1000

      do k=1,km
        c_den(k,1) = -2.55e-4
        ab(1+4,k)   = c_den(k,1)
        ro0_den(k) = 0.025022
        ab(2,k)    = ro0_den(k)
        to_den(k)  = 19.0
        ab(3,k)    = to_den(k)
        so_den(k)  = 0.0
        ab(4,k)    = so_den(k)
        do i=2,9
           c_den(k,i) = c0
           ab(i+4,k) = c_den(k,i)
        enddo
      enddo
# endif

      write (stdout,9531)
      write (stdout,9532)
     & (i,z(i),tmin(i),tmax(i),smin(i),smax(i),tmink(i),tmaxk(i),i=1,km)

# if defined O_knudsen
      write (stdout,'(/a/a/)')
     & 'Note: density coefficients were calculated using the "knudsen"'
     &,'      formulation.'
# else
      write (stdout,'(/a/a/)')
     & 'Note: density coefficients were calculated using the UNESCO'
     &,'      formulation.'
# endif

 9531 format('    density coefficients were calculated  ',
# if defined O_knudsen
     &       /'     employing the Knudsen-Ekman equation of state ',
# else
     &       /'    (employing the UNESCO equation of state)',
# endif
     &       /'    and are valid for the following depths and',
     &       ' both INSITU and POTENTIAL T and S ranges'
     &       /' ',t7,'k',t14,'depth',t27,'tmin',t37,
     &       'tmax',t52,'smin',t62,'smax',t70,'Pot tmin',t80,'Pot tmax')
 9532 format(' ',t5,i3,t12,f7.2,'e2',t25,f7.3,t35,f7.3,t50,f7.4,
     &       t60,f7.4,t70,f7.3,t80,f7.3)

      return
      end

      subroutine knuekm (t, s, d, rho)
!=======================================================================
!     this subroutine calculates the density of seawater using the
!     Knudsen-Ekman equation of state.

!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [per mil]
!       depth (d): [meters of depth, to approximate pressure]
!     output [units]:
!       density (rho): sigma units

!     reference:
!        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
!          New York, 1962, pp 3-30.

!-----------------------------------------------------------------------

      implicit none

      real (kind=8) eps2, t2, t, t3, s2, s, s3, f1, f2, f3, fs, sigma
      real (kind=8) a, d, b1, b2, b, co, alpha, rho

      data eps2/1.d-16/

      t2 = t*t
      t3 = t2*t
      s2 = s*s
      s3 = s2*s
      f1 = -1.0d0 * (t - 3.98d0)**2 * (t + 2.83d2) /
     &     (5.0357d2*(t + 6.726d1))
      f2 = t3*1.0843d-6 - t2*9.8185d-5 + t*4.786d-3
      f3 = t3*1.6670d-8 - t2*8.1640d-7 + t*1.803d-5
      fs = s3*6.76786136d-6 - s2*4.8249614d-4 + s*8.14876577d-1

      sigma= f1 + (fs + 3.895414d-2)*
     &      (1.0d0 - f2 + f3*(fs - 2.2584586d-1))

      a= d*1.0d-4*(1.055d2 + t*9.50d0 - t2*1.58d-1 - d*t*1.5d-4)  -
     &   (2.27d2 + t*2.833d1 - t2*5.51d-1 + t3*4.0d-3)
      b1 = (fs - 2.81324d1)*1.d-1
      b2 = b1 * b1
      b  = -b1* (1.473d2 - t*2.72d0 + t2*4.0d-2 - d*1.0d-4*
     &     (3.24d1 - 0.87d0*t + 2.0d-2*t2))
      b  = b + b2*(4.5d0 - 1.0d-1*t - d*1.0d-4*(1.8d0 - 6.0d-2*t))
      co = 4.886d3/(1.0d0 + 1.83d-5*d)

      alpha = d*1.0d-6*(co + a + b)

      rho = (sigma + alpha)/(1.d0 - 1.0d-3*alpha)

      return
      end

      subroutine lsqsl2 (ndim,a,d,w,b,x,irank,in,itmax,it,ieq,enorm,eps1
     &,nhdim,aa,r,s)

!     this routine is a modification of lsqsol. march,1968. r. hanson.
!     linear least squares solution

!     this routine finds x such that the euclidean length of
!     (*) ax-b is a minimum.

!     here a has k rows and n columns, while b is a column vector with
!     k components.

!     an orthogonal matrix q is found so that qa is zero below
!     the main diagonal.
!     suppose that rank (a)=r
!     an orthogonal matrix s is found such that
!     qas=t is an r x n upper triangular matrix whose last n-r columns
!     are zero.
!     the system tz=c (c the first r components of qb) is then
!     solved. with w=sz, the solution may be expressed
!     as x = w + sy, where w is the solution of (*) of minimum euclid-
!     ean length and y is any solution to (qas)y=ty=0.

!     iterative improvements are calculated using residuals and
!     the above procedures with b replaced by b-ax, where x is an
!     approximate solution.

      implicit none

      integer nhdim, k, n, it, isw, l, m, j1, j2, j3, j4, j5, j6, j7
      integer j8, j9, lm, irank, in, j, i, ieq, n2, n7, n1, ip, km
      integer n8, ipm1, ii, n3, ipp1, n4, irp1, irm1, n6, itmax, ns
      integer k1, d, w, ndim

      logical erm

      real (kind=8) top1, top, eps2, am, a1, eps1, sp, enorm, top2
      real (kind=8) enm1, a2, sj, dp, up, bp, aj

!     in=1 for first entry.
!                   a is decomposed and saved. ax-b is solved.
!     in = 2 for subsequent entries with a new vector b.
!     in=3 to restore a from the previous entry.
!     in=4 to continue the iterative improvement for this system.
!     in = 5 to calculate solutions to ax=0, then store in the array h.
!     in  =  6   do not store a  in aa.  obtain  t = qas, where t is
!     min(k,n) x min(k,n) and upper triangular. now return.do not obtain
!     a solution.
!     no scaling or column interchanges are performed.
!     in  =  7   same as with  in = 6  except that soln. of min. length
!                is placed into x. no iterative refinement.  now return.
!     column interchanges are performed. no scaling is performed.
!     in  = 8    set addresses. now return.

!     options for computing  a matrix product   y*h  or  h*y are
!     available with the use of the entry points  myh and mhy.
!     use of these options in these entry points allow a great saving in
!     storage required.

      include "stdunits.h"

      real (kind=8) a(d,w), b(d), aa(d,w), s(d+72), x(w), h(nhdim,nhdim)
      real (kind=8) r(d+36)
# if !defined O_read_denscoef
!     d = depth of matrix.
!     w = width of matrix.
      k=d
      n=w
      erm=.true.
      top1 = 0.0
      top = 0.0
      eps2 = 0.0

!     if it=0 on entry, the possible error message will be suppressed.

      if (it.eq.0) erm=.false.

!     ieq = 2      if column scaling by least max. column length is
!     to be performed.

!     ieq = 1       if scaling of all components is to be done with
!     the scalar max(abs(aij))/k*n.

!     ieq = 3 if column scaling as with in =2 will be retained in
!     rank deficient cases.

!     the array s must contain at least max(k,n) + 4n + 4min(k,n) cells
!        the   array r must contain k+4n s.p. cells.

!     the last card controls desired relative accuracy.
!     eps1  controls  (eps) rank.

      isw=1
      l=min0(k,n)
      m=max0(k,n)
      j1=m
      j2=n+j1
      j3=j2+n
      j4=j3+l
      j5=j4+l
      j6=j5+l
      j7=j6+l
      j8=j7+n
      j9=j8+n
      lm=l
      if (irank.ge.1.and.irank.le.l) lm=irank
      if (in.eq.6) lm=l
      if (in.eq.8) return

!     return after setting addresses when in=8.

      if (in.eq.1) goto 10
      if (in.eq.2) goto 360
      if (in.eq.3) goto 810
      if (in.eq.4) goto 390
      if (in.eq.5) goto 830
      if (in.eq.6) goto 10
      if (in.eq.7) goto 10

!     equilibrate columns of a (1)-(2).

!     (1)

   10 continue

!     save data when in = 1.

      if (in.gt.5) go to 30
        do j=1,n
          do i=1,k
            aa(i,j)=a(i,j)
          enddo
        enddo
   30 continue
      if (ieq.eq.1) go to 60
        do j=1,n
!          am=0.e0
          am=1.e-30
          do i=1,k
            am = max(am,abs(a(i,j)))
          enddo

!         s(m+n+1)-s(m+2n) contains scaling for output variables.

          n2=j2+j
          if (in.eq.6) am=1.d0
          s(n2)=1.d0/am
          do i=1,k
            a(i,j)=a(i,j)*s(n2)
          enddo
        enddo
        go to 100
60      continue
!        am=0.d0
          am=1.e-30
        do j=1,n
          do i=1,k
            am= max(am,abs(a(i,j)))
          enddo
        enddo
        am=am/float(k*n)
        if (in.eq.6) am=1.d0
        do j=1,n
          n2=j2+j
          s(n2)=1.d0/am
        enddo
        do j=1,n
          n2=j2+j
          do i=1,k
            a(i,j)=a(i,j)*s(n2)
          enddo
        enddo
!       compute column lengths with d.p. sums finally rounded to s.p.

!     (2)

100     continue
        do j=1,n
          n7=j7+j
          n2=j2+j
          s(n7)=s(n2)
        enddo

!       s(m+1)-s(m+ n) contains variable permutations.

!       set permutation to identity.

        do j=1,n
          n1=j1+j
          s(n1)=j
        enddo

!       begin elimination on the matrix a with orthogonal matrices .

!       ip=pivot row

        do ip=1,lm
          dp=0.d0
          km=ip
          do j=ip,n
            sj=0.d0
            do i=ip,k
              sj=sj+a(i,j)**2
            enddo
            if (dp.gt.sj) go to 140
              dp=sj
              km=j
              if (in.eq.6) go to 160
140         continue
          enddo

!         maximize (sigma)**2 by column interchange.

!         suppress column interchanges when in=6.

!         exchange columns if necessary.

          if (km.eq.ip) go to 160
            do i=1,k
              a1=a(i,ip)
              a(i,ip)=a(i,km)
              a(i,km)=a1
            enddo

!           record permutation and exchange squares of column lengths.

            n1=j1+km
            a1=s(n1)
            n2=j1+ip
            s(n1)=s(n2)
            s(n2)=a1
            n7=j7+km
            n8=j7+ip
            a1=s(n7)
            s(n7)=s(n8)
            s(n8)=a1
160       continue
          if (ip.eq.1) go to 180
            a1=0.d0
            ipm1=ip-1
            do i=1,ipm1
              a1=a1+a(i,ip)**2
            enddo
            if (a1.gt.0.d0) go to 190
180       continue
          if (dp.gt.0.d0) go to 200

!           test for rank deficiency.

190         continue
            if (sqrt(dp/a1).gt.eps1) go to 200
            if (in.eq.6) go to 200
            ii=ip-1
            if (erm) write (stdout,1140) irank,eps1,ii,ii
            irank=ip-1
            erm=.false.
            go to 260

!           (eps1) rank is deficient.

200       continue
          sp=sqrt(dp)

!         begin front elimination on column ip.

!         sp=sqroot(sigma**2).

          bp=1.d0/(dp+sp*abs(a(ip,ip)))

!         store beta in s(3n+1)-s(3n+l).

          if (ip.eq.k) bp=0.d0
          n3=k+2*n+ip
          r(n3)=bp
          up=sign(dble(sp)+abs(a(ip,ip)),dble(a(ip,ip)))
          if (ip.ge.k) go to 250
            ipp1=ip+1
            if (ip.ge.n) go to 240
              do j=ipp1,n
                sj=0.d0
                do i=ipp1,k
                  sj=sj+a(i,j)*a(i,ip)
                enddo
                sj=sj+up*a(ip,j)
                sj=bp*sj

!               sj=yj now

                do i=ipp1,k
                  a(i,j)=a(i,j)-a(i,ip)*sj
                enddo
                a(ip,j)=a(ip,j)-sj*up
              enddo
240         continue
            a(ip,ip)=-sign(sp,a(ip,ip))

            n4=k+3*n+ip
            r(n4)=up
250       continue
        enddo
        irank=lm
260     continue
        irp1=irank+1
        irm1=irank-1
        if (irank.eq.0.or.irank.eq.n) go to 360
        if (ieq.eq.3) go to 290

!       begin back processing for rank deficiency case
!       if irank is less than n.

        do j=1,n
          n2=j2+j
          n7=j7+j
          l=min0(j,irank)

!         unscale columns for rank deficient matrices when ieq.ne.3.

          do i=1,l
            a(i,j)=a(i,j)/s(n7)
          enddo
          s(n7)=1.d0
          s(n2)=1.d0
        enddo
290     continue
        ip=irank
300     continue
        sj=0.d0
        do j=irp1,n
          sj=sj+a(ip,j)**2
        enddo
        sj=sj+a(ip,ip)**2
        aj=sqrt(sj)
        up=sign(dble(aj)+abs(a(ip,ip)),dble(a(ip,ip)))

!       ip th element of u vector calculated.

        bp=1.d0/(sj+abs(a(ip,ip))*aj)

!       bp = 2/length of u squared.

        ipm1=ip-1
        if (ipm1.le.0) go to 340
        do i=1,ipm1
          dp=a(i,ip)*up
          do j=irp1,n
            dp=dp+a(i,j)*a(ip,j)
          enddo
          dp=dp/(sj+abs(a(ip,ip))*aj)

!         calc. (aj,u), where aj=jth row of a

          a(i,ip)=a(i,ip)-up*dp

!         modify array a.

          do j=irp1,n
            a(i,j)=a(i,j)-a(ip,j)*dp
          enddo
        enddo
340     continue
        a(ip,ip)=-sign(dble(aj),dble(a(ip,ip)))

!       calc. modified pivot.

!       save beta and ip th element of u vector in r array.

        n6=k+ip
        n7=k+n+ip
        r(n6)=bp
        r(n7)=up

!       test for end of back processing.

        if (ip-1.gt.0) goto 350
        if (ip-1.le.0) goto 360
350     continue
        ip=ip-1
        go to 300
360     continue
        if (in.eq.6) return
        do j=1,k
          r(j)=b(j)
        enddo
        it=0

!       set initial x vector to zero.

        do j=1,n
          x(j)=0.d0
        enddo
        if (irank.eq.0) go to 690

!       apply q to rt. hand side.

390     continue
        do ip=1,irank
          n4=k+3*n+ip
          sj=r(n4)*r(ip)
          ipp1=ip+1
          if (ipp1.gt.k) go to 410
            do i=ipp1,k
              sj=sj+a(i,ip)*r(i)
            enddo
410       continue
          n3=k+2*n+ip
          bp=r(n3)
          if (ipp1.gt.k) go to 430
            do i=ipp1,k
              r(i)=r(i)-bp*a(i,ip)*sj
            enddo
430       continue
          r(ip)=r(ip)-bp*r(n4)*sj
        enddo
        do j=1,irank
          s(j)=r(j)
        enddo
        enorm=0.d0
        if (irp1.gt.k) go to 510
          do j=irp1,k
            enorm=enorm+r(j)**2
          enddo
          enorm=sqrt(enorm)
          go to 510
460       continue
          do j=1,n
            sj=0.d0
            n1=j1+j
            ip=s(n1)
            do i=1,k
              sj=sj+r(i)*aa(i,ip)
            enddo

!           apply at to rt. hand side.
!           apply scaling.

            n7=j2+ip
            n1=k+n+j
            r(n1)=sj*s(n7)
          enddo
          n1=k+n
          s(1)=r(n1+1)/a(1,1)
          if (n.eq.1) go to 510
          do j=2,n
            n1=j-1
            sj=0.d0
            do i=1,n1
              sj=sj+a(i,j)*s(i)
            enddo
            n2=k+j+n
            s(j)=(r(n2)-sj)/a(j,j)
          enddo

!         entry to continue iterating.  solves tz = c = 1st irank
!         components of qb .

510     continue
        s(irank)=s(irank)/a(irank,irank)
        if (irm1.eq.0) go to 540
          do j=1,irm1
            n1=irank-j
            n2=n1+1
            sj=0.
            do i=n2,irank
              sj=sj+a(n1,i)*s(i)
            enddo
            s(n1)=(s(n1)-sj)/a(n1,n1)
          enddo

!         z calculated.  compute x = sz.

540     continue
        if (irank.eq.n) go to 590
          do j=irp1,n
            s(j)=0.d0
          enddo
          do i=1,irank
            n7=k+n+i
            sj=r(n7)*s(i)
            do j=irp1,n
              sj=sj+a(i,j)*s(j)
            enddo
            n6=k+i
            do j=irp1,n
              s(j)=s(j)-a(i,j)*r(n6)*sj
            enddo
            s(i)=s(i)-r(n6)*r(n7)*sj
          enddo

!         increment for x of minimal length calculated.

590     continue
        do i=1,n
          x(i)=x(i)+s(i)
        enddo
        if (in.eq.7) go to 750

!         calc. sup norm of increment and residuals

          top1=0.d0
          do j=1,n
            n2=j7+j
            top1= max(top1,abs(s(j))*s(n2))
          enddo
          do i=1,k
            sj=0.d0
            do j=1,n
              n1=j1+j
              ip=s(n1)
              n7=j2+ip
              sj=sj+aa(i,ip)*x(j)*s(n7)
            enddo
            r(i)=b(i)-sj
          enddo
        if (itmax.le.0) go to 750

!         calc. sup norm of x.

          top=0.d0
          do j=1,n
            n2=j7+j
            top= max(top,abs(x(j))*s(n2))
          enddo

!         compare relative change in x with tolerance eps .

          if (top1-top*eps2.lt.0.) goto 690
          if (top1-top*eps2.ge.0.) goto 650
650   continue
      if (it-itmax.lt.0) goto 660
      if (it-itmax.ge.0) goto 680
660   continue
      it=it+1
      if (it.eq.1) go to 670
      if (top1.gt..25*top2) go to 690
670   continue
      top2=top1
      if (isw.eq.1) go to 390
      if (isw.eq.2) go to 460
680   continue
      it=0
690   continue
      sj=0.d0
      do j=1,k
        sj=sj+r(j)**2
      enddo
      enorm=sqrt(sj)
      if (irank.eq.n.and.isw.eq.1) go to 710
      go to 730
710   continue
      enm1=enorm

!     save x array.

      do j=1,n
        n1=k+j
        r(n1)=x(j)
      enddo
      isw=2
      it=0
      go to 460

!     choose best solution

730   continue
      if (irank.lt.n) go to 750
      if (enorm.le.enm1) go to 750
      do j=1,n
        n1=k+j
        x(j)=r(n1)
      enddo
      enorm=enm1

!     norm of ax - b located in the cell enorm .

!     rearrange variables.

750   continue
      do j=1,n
        n1=j1+j
        s(j)=s(n1)
      enddo
      do j=1,n
        do i=j,n
          ip=s(i)
          if (j.eq.ip) go to 780
        enddo
780     continue
        s(i)=s(j)
        s(j)=j
        sj=x(j)
        x(j)=x(i)
        x(i)=sj
      enddo

!     scale variables.

      do j=1,n
        n2=j2+j
        x(j)=x(j)*s(n2)
      enddo
      return

!     restore a.

810   continue
      do j=1,n
        n2=j2+j
        do i=1,k
          a(i,j)=aa(i,j)
        enddo
      enddo
      return

!     generate solutions to the homogeneous equation ax = 0.

  830 if (irank.eq.n) return
      ns=n-irank
      do i=1,n
        do j=1,ns
          h(i,j)=0.d0
        enddo
      enddo
      do j=1,ns
        n2=irank+j
        h(n2,j)=1.d0
      enddo
      if (irank.eq.0) return
      do j=1,irank
        do i=1,ns
          n7=k+n+j
          sj=r(n7)*h(j,i)

!         this part of the code should not be executed. However some
!         compilers generate warning msgs that "irp1" is not defined
!         so it is arbitrarily set here to fool the compilers

          irp1 = 1
          if (i .gt. 0) then
            write (stdout,*) ' Error in lsqsl2: search for this msg'
            stop
          endif
          do k1=irp1,n
            sj=sj+h(k1,i)*a(j,k1)
          enddo
          n6=k+j
          bp=r(n6)
          dp=bp*r(n7)*sj
          a1=dp
          a2=dp-a1
          h(j,i)=h(j,i)-(a1+2.*a2)
          do k1=irp1,n
            dp=bp*a(j,k1)*sj
            a1=dp
            a2=dp-a1
            h(k1,i)=h(k1,i)-(a1+2.*a2)
          enddo
        enddo
      enddo

!     rearrange rows of solution matrix.

      do j=1,n
        n1=j1+j
        s(j)=s(n1)
      enddo
      do j=1,n
        do i=j,n
          ip=s(i)
          if (j.eq.ip) go to 900
        enddo
900     continue
        s(i)=s(j)
        s(j)=j
        do k1=1,ns
          a1=h(j,k1)
          h(j,k1)=h(i,k1)
          h(i,k1)=a1
        enddo
      enddo
# endif
      return

 1140 format (/'warning. irank has been set to',i4,'  but(',1pe10.3,
     1 ') rank is',i4,'.  irank is now taken as ',i4)
      end

      subroutine potem (t, s, p, theta)

!=======================================================================
!     this subroutine calculates potential temperature as a function
!     of in-situ temperature, salinity, and pressure.

!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [per mil]
!       pressure (p): [decibars, approx. as meters of depth]
!     output [units]:
!       potential temperature (theta): [degrees centigrade]

!     references:
!        based on Fofonoff and Froese (1958) as shown in ...
!        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
!          New York, 1962, page 17, table iv.

!-----------------------------------------------------------------------

      implicit none

      real (kind=8) b1, p, b2, t, t2, t3, b3, b4, b5, s, b6, s2
      real (kind=8) p2, b7, b8, b9, b10, b11, potmp, theta

      b1    = -1.60d-5*p
      b2    = 1.014d-5*p*t
      t2    = t*t
      t3    = t2*t
      b3    = -1.27d-7*p*t2
      b4    = 2.7d-9*p*t3
      b5    = 1.322d-6*p*s
      b6    = -2.62d-8*p*s*t
      s2    = s*s
      p2    = p*p
      b7    = 4.1d-9*p*s2
      b8    = 9.14d-9*p2
      b9    = -2.77d-10*p2*t
      b10   = 9.5d-13*p2*t2
      b11   = -1.557d-13*p2*p
      potmp = b1+b2+b3+b4+b5+b6+b7+b8+b9+b10+b11
      theta = t-potmp

      return
      end

      subroutine unesco (t, s, pin, rho)

!=======================================================================
!     this subroutine calculates the density of seawater using the
!     standard equation of state recommended by unesco(1981).

!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [practical salinity units]
!       pressure (pin): [decibars, approx. as meters of depth]
!     output [units]:
!       density(rho): kilograms per cubic meter

!     references:
!        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical
!         Series No. 30. Academic Press, London, 1982, pp 599-600.
!        UNESCO, 10th report of the joint panel on oceanographic tables
!          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,
!          Paris, 1981.

!-----------------------------------------------------------------------

      implicit none

      real (kind=8) c1p5, p, pin, rw, t, rsto, s, xkw, xksto, xkstp, rho

      c1p5 = 1.5d0

!  convert from depth [m] (decibars) to bars
      p = pin * 1.0d-1

      rw =     9.99842594d2 + 6.793952d-2*t - 9.095290d-3*t**2
     &        + 1.001685d-4*t**3 - 1.120083d-6*t**4 + 6.536332d-9*t**5

      rsto =   rw + (8.24493d-1 - 4.0899d-3*t + 7.6438d-5*t**2
     &        - 8.2467d-7*t**3 + 5.3875d-9*t**4) * s
     &       + (-5.72466d-3 + 1.0227d-4*t - 1.6546d-6*t**2) * s**c1p5
     &       + 4.8314d-4 * s**2

      xkw =     1.965221d4 + 1.484206d2*t - 2.327105d0*t**2 +
     &         1.360477d-2*t**3 - 5.155288d-5*t**4

      xksto =   xkw + (5.46746d1 - 6.03459d-1*t + 1.09987d-2*t**2
     &        - 6.1670d-5*t**3) * s
     &       + (7.944d-2 + 1.6483d-2*t - 5.3009d-4*t**2) * s**c1p5

      xkstp =   xksto + (3.239908d0 + 1.43713d-3*t + 1.16092d-4*t**2
     &        - 5.77905d-7*t**3) * p
     &       + (2.2838d-3 - 1.0981d-5*t - 1.6078d-6*t**2) * p * s
     &       + 1.91075d-4 * p * s**c1p5
     &       + (8.50935d-5 - 6.12293d-6*t + 5.2787d-8*t**2) * p**2
     &       + (-9.9348d-7 + 2.0816d-8*t + 9.1697d-10*t**2) * p**2 * s

      rho =    rsto / (1.0d0 - p/xkstp)

#endif
      return
      end
