      subroutine filtr (s, im, mm, n, iss)

#if defined O_fourfil
!=======================================================================
!                                                                    ===
!  filter fourier analyses the arrays of various                     ===
!         physical quantities, then truncates the series and         ===
!         resynthesizes the filtered quantities where:               ===
!             s  =the string to be filtered                          ===
!             im =the length of s                                    ===
!             mm =1 (cosine series, deriv at bndry pts=0)            ===
!                =2 (  sine series,          bndry pts=0)            ===
!                =3 (full series, cyclic)                            ===
!             n  =number of waves to keep                            ===
!             iss=0 (cant use fourier coefs from previous call)      ===
!             iss>0 (can  use fourier coefs from previous call)      ===
!=======================================================================

      implicit none

!---------------------------------------------------------------------
!     define global data
!---------------------------------------------------------------------

      integer imtx2, ni, imtd2, lqmsum, lhsum, imtx4, imtx8, imtimt
      integer imp1x2, im, mm, n, iss, imsave, i, imm1, imqc, nmax
      integer nmaxp1, lcy, lh, lhm1, lqm, l2cy, lcym1, lcyp1, imx2
      integer imx4, imx8, nprint, maxind, ncyc, maxndx, npwr, np, j
      integer ioff1, ioff2, joff, ioff, jbase, ibase

      real ssm, fimr, cc1, cc2, fnorm, ssum, fim, stemp, fact1
      real fact2, genadj

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "scalar.h"
      include "switch.h"

!---------------------------------------------------------------------
!     define local data and dimension argument arrays
!---------------------------------------------------------------------

      parameter (imtx2=imt*2,ni=imt)
      parameter (imtd2=imt/2,lqmsum=imtd2*(imt-imtd2),lhsum=imt*imtp1/2)
      parameter (imtx4=imt*4,imtx8=imt*8,imtimt=imt*imt)
      parameter (imp1x2=imtp1*2)

      integer icbase(imtp1), idbase(imtp1), ind(imtx8), indx(imtx8)

      common /cfilt_i/ ind, idbase, icbase
      common /cfilt_i/ imsave, jbase, ibase

      real*8 cossav(lqmsum)
      real temp(imtx4), denmsv(lhsum), cosnpi(imt)
      real cof(imtx8), cosine(imtx8), ftarr(imtimt)
      real denom(imtx4), s(imt), sprime(imt), circle(4)

      common /cfilt_d/ cossav
      common /cfilt_r/ denmsv, cosnpi, ftarr, circle

!     data circle /0.,-1.,0.,1./
      circle(1) = 0.0
      circle(2) = -1.0
      circle(3) = 0.0
      circle(4) = 1.0

!---------------------------------------------------------------------
!  begin executable code
!---------------------------------------------------------------------

      if (im.lt.1 .or. mm.lt.1 .or. mm.gt.3 .or. n.lt.0 .or. iss.lt.0)
     $  then
        write (stdout,99) im, mm, n, iss
        write (stderr,99) im, mm, n, iss
        stop ' filtr 1'
      endif

      if (first) then

!       this section sets up tables for filter; it must be called once
!       per execution of ocean

!       note: lqmsum is the sum of (im-1)/2 for im=1,imtp1
!             lhsum is the sum of im-1 for im=1,imtp1

        imsave = im

!       assemble index array

        do 100 i=1,imtx8
          ind(i) = i
100     continue

!       calculate and save all cosines which will be needed

        ibase = 0
        jbase = 0

        do 200 im=1,imtp1
          fimr = c1/float(im)
          imm1 = im-1
          if (imm1.eq.0) goto 181
          do 180 i=1,imm1
            denmsv(ibase+i) = c1/(c1-cos(pi*float(i)*fimr))
180       continue
181       continue
          idbase(im) = ibase
          ibase = ibase + imm1
          imqc = (im-1)/2
          if (imqc .eq. 0) goto 191
          do 190 i=1,imqc
            cossav(jbase+i) = cos(pi*float(i)*fimr)
190       continue
191       continue
          icbase(im) = jbase
          jbase = jbase + imqc
200     continue

!       calculate adjustments for general fourier case if im=2*n

        do 300 im=1,imt
          cosnpi(im) = circle(mod(im-1,4)+1)
300     continue

        im = imsave
      endif

!     calculate some useful constants

      if (mm.eq.2 .and. n.eq.0) then

        do 400 i=1,im
          s(i) = c0
400     continue

        goto 3201
      endif

      if (mm .eq. 1) then
        nmax = n - 1
      else
        nmax = n
      endif

      nmaxp1 = nmax + 1
      cc1 = p5*float(nmax) + p25
      cc2 = float(nmax) + p5

      if (mm .eq. 2) then
        lcy   = 2*(im + 1)
        fnorm = c2/float(im + 1)
      else
        lcy   = 2*im
        fnorm = c2/float(im)
      endif

      lh    = lcy/2
      lhm1  = lh - 1
      lqm   = (lh - 1)/2
      l2cy  = 2*lcy
      lcym1 = lcy - 1
      lcyp1 = lcy + 1
      imx2  = im*2
      imx4  = im*4
      imx8  = im*8

!     average incoming array

      ssum = c0

      do 500 i=1,im
        ssum = ssum + s(i)
500   continue

!     mm = 1  derivative must be zero at boundaries (cosine)
!     mm = 2  value must be zero at boundaries (sine)
!     mm = 3  cyclic boundary conditions (general fourier series)

      fim   = float(im)
      fimr  = c1/fim
      stemp = ssum*fimr

      if (n.gt.1 .or. mm.ne.1) goto 601

      do 600 i=1,im
        s(i)=stemp
600   continue

      go to 3201
601   continue

      if (mm .ne. 2) then

        do 700 i=1,im
          s(i) = s(i) - stemp
700     continue

      endif

      if (iss .gt. 0) goto 2501

!     assemble appropriate 1-cycle (2*pi) cosine array

!     use stored 1/4 cycle to calculate first 1/2 cycle

      jbase = icbase(lh)

      do 800 i=1,lqm
        cosine(i) = cossav(jbase+i)
800   continue

      do 900 i=1,lqm
        cosine(lh-i) = -cossav(jbase+i)
900   continue

!     fill in cos(pi/2) if lh is even

      if (2*(lqm+1) .eq. lh) cosine(lqm+1) = c0

!     fill in cos(pi) in any case

      cosine(lh) = -c1

!     fill in rest of cycle

      do 1000 i=1,lh
        cosine(lh+i) = -cosine(i)
1000  continue

!     assemble denominator array

      ibase = idbase(lh)

      do 1100 i=1,lhm1
        denom(i) = p25*denmsv(ibase+i)
1100  continue

      denom(lh) = 0.125

      do 1200 i=1,lhm1
        temp(i) = denom(lh-i)
1200  continue

      do 1300 i=1,lhm1
        denom(lh+i) = temp(i)
1300  continue

      nprint = 0
      denom(lcy) = c0

      do 1400 i=lcyp1,imx4
        denom(i) = denom(i-lcy)
1400  continue

!     assemble appropriate subscript arrays

!     calculate needed indices

      if (mm.eq.3) then
        fact1 = 2*nmax
        fact2 = 2*nmaxp1
      else
        fact1 = nmax
        fact2 = nmaxp1
      endif

      do 1500 i=1,imx4
        indx(i) = ind(i)*fact1
1500  continue

      do 1600 i=1,imx4
        indx(imx4+i) = ind(i)*fact2
1600  continue

!     calculate parameters for reducing indices

      maxind = imx4*fact2
      ncyc   = (maxind-1)/lcy + 1
      maxndx = lcy
      if (maxndx .ge. maxind) goto 1801

      do 1700 npwr=1,ncyc+2
        maxndx = 2*maxndx
        if (maxndx .ge. maxind) goto 1701
1700  continue

      write (stdout,999)
      write (stderr,999)
      stop ' filtr 2'

1701  continue

      do 1800 np=1,npwr
        maxndx = maxndx/2
        do 1790 i=1,imx8
          if (indx(i) .gt. maxndx) indx(i) = indx(i) - maxndx
1790    continue
1800  continue

1801  continue

!     gather coefficients

      do 1900 j=1,imx8
        cof(j) = cosine(indx(j))
1900  continue

!     assemble transformation array which will filter s

      if (mm.eq.1) then

!       cosine transform

        ioff1 = lcy
        ioff2 = lcy + imx4

        do 2000 j=1,im
          joff = (j-1)*imt
          do 1990 i=1,im
            ftarr(joff+i) =
     $         (cof(i-j+ioff1) - cof(i-j+ioff2)) *denom(i-j+ioff1) +
     $         (cof(i + j - 1) - cof(imx4+i+j-1))*denom(i+j-1) - p5
1990      continue
2000    continue

        do 2100 j=1,im
          ftarr(j*imtp1-imt) = ftarr(j*imtp1-imt) + cc1
2100    continue

      elseif (mm .eq. 2) then

!       sine transform

        ioff1 = lcy
        ioff2 = lcy + imx4

        do 2200 j=1,im
          joff = (j-1)*imt
          do 2190 i=1,im
            ftarr(joff+i) =
     $         (cof(i-j+ioff1) - cof(i-j+ioff2))*denom(i-j+ioff1) -
     $         (cof(i + j)     - cof(imx4+i+j)) *denom(i+j)
2190      continue
2200    continue

        do 2300 j=1,im
          ftarr(j*imtp1-imt) = ftarr(j*imtp1-imt) + cc1
2300   continue

      elseif (mm.eq.3) then

!       general fourier transform

        if (2*n .eq. im) then
          genadj = p5
        else
          genadj = c0
        endif

        ioff1 = lcy
        ioff2 = lcy + imx4

        do 2400 j=1,im
          joff = (j-1)*imt
          do 2390 i=1,im
            ftarr(joff+i) = (c2*(cof(i-j+ioff1) - cof(i-j+ioff2)))
     $          *denom(2*i-2*j+ioff1) - p5 - genadj*cosnpi(i)*cosnpi(j)
2390      continue
2400    continue

        do 2500 j=1,im
          ftarr(j*imtp1-imt) = ftarr(j*imtp1-imt) + cc2
2500    continue

      endif

!     filter s

2501  continue

      do 2600 i=1,im
        sprime(i) = c0
2600  continue

!     note that ftarr(j,i)=ftarr(i,j), so following is legal

      do 2700 i=1,im
        ioff = (i-1)*imt
        do 2690 j=1,im
          sprime(j) = sprime(j) + s(i)*ftarr(ioff+j)
2690    continue
2700  continue

      do 2800 i=1,im
        sprime(i) = fnorm*sprime(i)
2800  continue

      if (mm.eq.2) then

        do 2900 i=1,im
          s(i) = sprime(i)
2900    continue

        goto 3201
      endif

      ssm = c0

      do 3100 i=1,im
        ssm = ssm + sprime(i)
3100  continue

      ssm = (ssum-ssm)*fimr

      do 3200 i=1,im
        s(i) = ssm+sprime(i)
3200  continue

3201  continue

   99 format (/' error => bad argument(s) in call to filtr'
     $       /' im,mm,n,iss = ',4i10)
  999 format (/' error => can not calculate parameters for reducing',
     $        ' indices in filtr')
#endif

      return
      end
