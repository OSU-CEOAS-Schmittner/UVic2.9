      function geoheatflux(lon,lat)
#if defined O_gthflx
!     written by Alex Hoffman
!     put into 2.9 by Andreas Schmittner
!     input arguments: latitude lat, longitude lon in degrees!
!     p is longitude and q is colatitude
      real*8 anm(0:12,0:12), bnm(0:12,0:12), p, con, gcon
      real*8 q, x, krt, sum, pprime, pp, qq, pi, qtest
      real*8 lon, lat, xt, yt, geoheatflux
      integer n, h, m, t, i, j
    
    
      anm(0,0)=86.674
      anm(1,0)=-12.999
      anm(1,1)=-2.689
      bnm(1,1)=-10.417
      anm(2,0)=-1.917
      anm(2,1)=4.578
      bnm(2,1)=1.022
      anm(2,2)=-14.076
      bnm(2,2)=6.507
      anm(3,0)=7.122
      anm(3,1)=-2.934
      bnm(3,1)=3.555
      anm(3,2)=7.232
      bnm(3,2)=-3.295
      anm(3,3)=10.299
      bnm(3,3)=4.646
      anm(4,0)=-3.511
      anm(4,1)=2.778
      bnm(4,1)=-1.873
      anm(4,2)=1.728
      bnm(4,2)=-2.546
      anm(4,3)=-4.822
      bnm(4,3)=.486
      anm(4,4)=4.408
      bnm(4,4)=-17.946
      anm(5,0)=5.316
      anm(5,1)=-1.984
      bnm(5,1)=-2.642
      anm(5,2)=2.167
      bnm(5,2)=3.835
      anm(5,3)=4.570
      bnm(5,3)=-6.087
      anm(5,4)=-8.353
      bnm(5,4)=10.283
      anm(5,5)=-6.896
      bnm(5,5)=-4.199
      anm(6,0)=-5.204
      anm(6,1)=2.795
      bnm(6,1)=3.162
      anm(6,2)=2.065
      bnm(6,2)=-2.889
      anm(6,3)=-2.740
      bnm(6,3)=-.252
      anm(6,4)=-.012
      bnm(6,4)=-1.897
      anm(6,5)=.637
      bnm(6,5)=.476
      anm(6,6)=3.739
      bnm(6,6)=7.849
      anm(7,0)=2.010
      anm(7,1)=.912
      bnm(7,1)=.116
      anm(7,2)=-6.044
      bnm(7,2)=-.179
      anm(7,3)=4.999
      bnm(7,3)=-.123
      anm(7,4)=-1.605
      bnm(7,4)=-3.721
      anm(7,5)=-.334
      bnm(7,5)=3.466
      anm(7,6)=-4.111
      bnm(7,6)=-.639
      anm(7,7)=4.126
      bnm(7,7)=-1.659
      anm(8,0)=2.621
      anm(8,1)=-1.376
      bnm(8,1)=1.795
      anm(8,2)=7.201
      bnm(8,2)=1.436
      anm(8,3)=-1.947
      bnm(8,3)=.679
      anm(8,4)=.204
      bnm(8,4)=1.171
      anm(8,5)=1.851
      bnm(8,5)=1.771
      anm(8,6)=3.579
      bnm(8,6)=-.250
      anm(8,7)=1.886
      bnm(8,7)=4.903
      anm(8,8)=-5.285
      bnm(8,8)=-4.412
      anm(9,0)=-.211
      anm(9,1)=3.140
      bnm(9,1)=.886
      anm(9,2)=-.360
      bnm(9,2)=-3.894
      anm(9,3)=-3.004
      bnm(9,3)=-2.056
      anm(9,4)=1.947
      bnm(9,4)=-2.511
      anm(9,5)=.328
      bnm(9,5)=-3.064
      anm(9,6)=1.030
      bnm(9,6)=-.745
      anm(9,7)=-4.117
      bnm(9,7)=-3.888
      anm(9,8)=6.529
      bnm(9,8)=3.889
      anm(9,9)=-4.084
      bnm(9,9)=-.082
      anm(10,0)=2.735
      anm(10,1)=-1.624
      bnm(10,1)=-1.998
      anm(10,2)=-1.309
      bnm(10,2)=1.333
      anm(10,3)=4.576
      bnm(10,3)=.641
      anm(10,4)=-4.506
      bnm(10,4)=.927
      anm(10,5)=-.363
      bnm(10,5)=-.927
      anm(10,6)=-4.528
      bnm(10,6)=-1.353
      anm(10,7)=-.952
      bnm(10,7)=1.810
      anm(10,8)=-1.104
      bnm(10,8)=-.739
      anm(10,9)=.129
      bnm(10,9)=.644
      anm(10,10)=4.164
      bnm(10,10)=-3.463
      anm(11,0)=-1.708
      anm(11,1)=.429
      bnm(11,1)=2.902
      anm(11,2)=2.106
      bnm(11,2)=.915
      anm(11,3)=-5.078
      bnm(11,3)=.595
      anm(11,4)=3.441
      bnm(11,4)=.907
      anm(11,5)=.784
      bnm(11,5)=2.762
      anm(11,6)=.158
      bnm(11,6)=.782
      anm(11,7)=-.377
      bnm(11,7)=-.355
      anm(11,8)=-.818
      bnm(11,8)=1.851
      anm(11,9)=3.654
      bnm(11,9)=1.336
      anm(11,10)=-1.765
      bnm(11,10)=4.245
      anm(11,11)=-.505
      bnm(11,11)=-3.520
      anm(12,0)=1.003
      anm(12,1)=-.689
      bnm(12,1)=-1.476
      anm(12,2)=-2.359
      bnm(12,2)=-.066
      anm(12,3)=3.863
      bnm(12,3)=.504
      anm(12,4)=.793
      bnm(12,4)=-1.034
      anm(12,5)=-1.761
      bnm(12,5)=-.267
      anm(12,6)=2.439
      bnm(12,6)=-2.484
      anm(12,7)=-2.080
      bnm(12,7)=3.714
      anm(12,8)=2.237
      bnm(12,8)=.809
      anm(12,9)=.289
      bnm(12,9)=-.838
      anm(12,10)=1.516
      bnm(12,10)=-4.821
      anm(12,11)=4.114
      bnm(12,11)=-.533
      anm(12,12)=-3.033
      bnm(12,12)=2.175


! identify constants
      pi = acos(-1.)
      con = 1./41840000
      gcon = 360./102

! Functions from Appendix 1 of Hamza et al. 2007
! Begin Loops
!CHANGE pi and qj to depend on xt and yt (i am assuming the array is built-in)
! Andreas Sep 8, 2015 the following line used to be yt=lat*pi/180. bug fixed
! thanks for Juan Muglia for finding it
      yt=(90.-lat)*pi/180.
      xt=lon*pi/180.
      qq=0
      qtest=0
      do n=0,12
         do m=0,n
            sum=0
            do t=0,int((n-m)/2)

               sum = ((((-1)**t)*factorial(2*n-2*t))/(factorial(t)*
     &         factorial(n-t)*factorial(n-m-2*t))*(cos(yt)**(n-m-2*t)
     &              ) )+sum
 

            enddo
            pprime = ((sin(yt)**m)/(2**n))*sum
            if (m.eq.0) then
               h=1
            else
               h=2
            endif
            krt = ((factorial(n+m)/factorial(n-m))/(h*(2*n+1.)))**(.5)
            pp = (pprime)/(krt)
         
            qq = (((anm(n,m)*cos(m*xt))+(bnm(n,m)*
     &           sin(m*xt)))*(pp))+qq

         enddo
      enddo
     
      geoheatflux = qq*con
  
      end
! End loops
! Define factorial function!     
      real function factorial(n)
      integer n
      real f
      f=1.
      do i=1,n
         f=i*f
      enddo
      factorial = f
#endif
      return
      end
! End Spatially Variable Heat Flux (spavar) !

     
