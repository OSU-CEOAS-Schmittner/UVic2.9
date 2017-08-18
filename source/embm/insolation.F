      subroutine zenith (points, t, dt, daylen, lat, lon, sindec, cosz)

#if defined O_embm
!-----------------------------------------------------------------------
! Calculate the cosine of the average zenith angle over a time period
!-----------------------------------------------------------------------

      implicit none

! points = number of points
! t      = start time (sec from GMT)
! dt     = time step (forced to be between 1 msec and 1 day - 1 msec)
! daylen = day length (sec)
! lat    = latitude (degrees)
! lon    = longitude (degrees)
! sindec = sine of solar declination
! cosz   = mean cosine of solar zenith angle

      integer, intent(in) :: points
      real, intent(in) :: t, dt, daylen
      real, intent(in) :: lat(points), lon(points)
      real, intent(in) :: sindec
      real, intent(out) :: cosz(points)

! sinlat = sine of latitude
! sinsin = product of the sines of solar declination and of latitude.
! coscos = product of the cosines of solar declination and of latitude.
! hld    = half length of the day in radians
! coshld = cosine of hld
! hat    = local hour angle at the start time.
! omegab = hour angle at beginning of timestep (radians from sunrise)
! omegae = hour angle at end of timestep (radians from sunrise)
! omegae = hour angle of end of sunset (radians from sunrise)
! omega1 = hour angle at beginning of integration (radians from sunrise)
! omega2 = hour angle at end of integration (radians from sunrise)
! omegas = hour angle at sunset (radians from sunrise)
! difsin = an intermediate difference of sines
! diftim = an intermediate difference of times
! trad   = start time (radians after midday GMT)
! dtrad  = timestep (radians after midday GMT)
! pi     = pi
! twopi  = 2*pi
! rtwopi = reciprocal of 2*pi
! secrad = seconds to radians converter
! degrad = degrees to radians converter
! msec   = seconds in a msec (used to restrict the range of dt)

      integer j

      real (kind=8) sinlat, sinsin, coscos, hld, coshld, hat
      real (kind=8) omegab, omegae, omega1, omega2, omegas
      real (kind=8) difsin, diftim, trad, dtrad
      real (kind=8) pi, twopi, rtwopi, secrad, degrad, msec

      parameter (msec=0.001)

      pi = 4.0*atan(1.0)
      twopi = 2.0*pi
      rtwopi = 0.5/pi
      secrad = twopi/daylen
      degrad = pi/180.

      trad = t*secrad - pi
      dtrad = dt
      dtrad = max(msec,min(daylen-msec,dtrad))*secrad

      do j=1,points

        sinlat = sin(lat(j)*degrad)
        sinsin = sindec*sinlat
        coscos = sqrt((1. - sindec**2)*(1. - sinlat**2))
        coshld = sinsin/(coscos + 1.e-20)

        if (coshld .lt. -1.) then
!         perpetual night (all longitudes)
          cosz(j) = 0.

        elseif (coshld .gt. 1.) then
!         perpetual day (all longitudes)
          hat = lon(j)*degrad + trad
          cosz(j) = (sin(hat + dtrad) - sin(hat))*coscos/dtrad + sinsin

        else
!         day and night (for at least some longitudes)
          hat = lon(j)*degrad + trad
          hld = acos(-coshld)
          omegab = hat + hld
!         times relative to sunrise but they must be kept in the range
!         0 to 2pi for the tests on their orders to work
          if (omegab.lt.0.)   omegab = omegab + twopi
          if (omegab.ge.twopi) omegab = omegab - twopi
          if (omegab.ge.twopi) omegab = omegab - twopi
          omegae = omegab + dtrad
          if (omegae.gt.twopi) omegae = omegae - twopi
          omegas = 2.*hld
          if (omegab .le. omegas .or. omegab .lt. omegae) then
             omega1 = omegab - hld
           else
             omega1 = - hld
          endif
          if (omegae .le. omegas) then
             omega2 = omegae - hld
           else
             omega2 = omegas - hld
          endif
!         case were the sun is not up during the timestep
          if (omegae.gt.omegab .and. omegab.gt.omegas) omega2 = omega1
          difsin = sin(omega2) - sin(omega1)
          diftim = omega2 - omega1
          if (diftim .lt. 0.) then
!           case where the sun sets and then rises again within the
!           timestep so the original integration is backward.
            difsin = difsin + 2.*sqrt(1.-coshld**2)
            diftim = diftim + 2.*hld
          endif
          if (diftim .eq. 0.) then
            cosz(j) = 0.
          else
            cosz(j) = (difsin*coscos/diftim + sinsin)*(diftim/dtrad)
          endif

        endif

      enddo

      return
      end

      subroutine decl (day, eccen, obliq, mvelp, lambm0, sindec, eccf)

!-----------------------------------------------------------------------
! Compute earth/orbit parameters using Berger, Andre.  1978
! "A Simple Algorithm to Compute Long-Term Variations of Daily
! Insolation".  Contribution 18, Institute of Astronomy and Geophysics,
! Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!-----------------------------------------------------------------------

      implicit none

! day    = Calendar day, including fraction
! eccen  = orbital eccentricity
! obliq  = obliquity (degrees)
! mvelp  = moving vernal equinox longitude of perihelion (degrees)
! lambm0 = Mean long of perihelion at vernal equinox (radians)
! sindec = sine of solar declination angle in rad
! eccf   = Earth-sun distance factor (ie. (1/r)**2)

      real, intent(in) :: day
      real, intent(in) :: eccen, obliq, mvelp, lambm0
      real, intent(out) :: sindec, eccf

! lambm  = Lambda m, mean long of perihelion (rad)
! lmm    = Intermediate argument involving lambm
! lamb   = Lambda, the earths long of perihelion
! invrho = Inverse normalized sun/earth distance
! sinl   = Sine of lmm
! dayspy = days per year
! ve     = day of vernal equinox assumes Jan 1 = day 1
! pi     = pi
! degrad = degree to radian converter
! obliqr = Earths obliquity (radians)
! mvelpp = moving vernal equinox long of perihelion plus pi (radians)

      real (kind=8) lambm, lmm, lamb, invrho, sinl, dayspy, ve
      real (kind=8) pi, degrad, obliqr, mvelpp

!     parameter (dayspy=365.0, ve=80.5)
!     to reproduce berger you need the following parameter settings:
      parameter (dayspy=365.25, ve=80.)

      pi = 4.*atan(1.)
      degrad = pi/180.
      obliqr = obliq*degrad
      mvelpp = (mvelp + 180.)*degrad

! Compute eccentricity factor and solar declination using
! day value where a round day (such as 213.0) refers to 0z at
! Greenwich longitude.

! Use formulae from Berger, Andre 1978: Long-Term Variations of Daily
! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci.
! 35:2362-2367.

! To get the earths true longitude (position in orbit; lambda in Berger
! 1978) which is necessary to find the eccentricity factor and
! declination, must first calculate the mean longitude (lambda m in
! Berger 1978) at the present day.  This is done by adding to lambm0
! (the mean longitude at the vernal equinox, set as March 21 at noon,
! when lambda=0; in radians) an increment (delta lambda m in Berger
! 1978) that is the number of days past or before (a negative increment)
! the vernal equinox divided by the days in a model year times the 2*pi
! radians in a complete orbit.

      lambm = lambm0 + (day - ve)*2.*pi/dayspy
      lmm   = lambm  - mvelpp

! The earths true longitude, in radians, is then found from
! the formula in Berger 1978:

      sinl  = sin(lmm)
      lamb  = lambm  + eccen*(2.*sinl + eccen*(1.25*sin(2.*lmm)
     &      + eccen*((13.0/12.0)*sin(3.*lmm) - 0.25*sinl)))

! Using the obliquity, eccentricity, moving vernal equinox longitude of
! perihelion (plus), and earths true longitude, the declination (delta)
! and the normalized earth/sun distance (rho in Berger 1978; actually
! inverse rho will be used), and thus the eccentricity factor (eccf),
! can be calculated from formulae given in Berger 1978.

      invrho = (1. + eccen*cos(lamb - mvelpp)) / (1. - eccen*eccen)

! Set solar declination and eccentricity factor

!     delta  = asin(sin(obliqr)*sin(lamb))
      sindec = sin(obliqr)*sin(lamb)
      eccf   = invrho*invrho

      return
      end

      subroutine orbit (year, eccen, obliq, mvelp, lambm0)

!-----------------------------------------------------------------------
! Calculate earths orbital parameters using Berger, Andre.  1978
! "A Simple Algorithm to Compute Long-Term Variations of Daily
! Insolation".  Contribution 18, Institute of Astronomy and Geophysics,
! Universite Catholique de Louvain, Louvain-la-Neuve, Belgium
!-----------------------------------------------------------------------

      implicit none

! year   = calender year of orbit (-ve => B.C. +ve => A.D.)
! eccen  = orbital eccentricity
! obliq  = obliquity in degrees
! mvelp  = moving vernal equinox longitude
! lambm0 = Mean long of perihelion at vernal equinox (radians)

      real, intent(in) :: year
      real, intent(inout) :: eccen, obliq, mvelp
      real, intent(out) :: lambm0

! obamp  = amplitudes for obliquity cos series (arc seconds)
! obrate = rates for obliquity cosine series (arc seconds/year)
! obphas = phases for obliquity cosine series (degrees)

      real (kind=8) obamp(47), obrate(47), obphas(47)

      data obamp(1:3)   /-2462.2214466, -857.3232075, -629.3231835/
      data obamp(4:6)   / -414.2804924, -311.7632587,  308.9408604/
      data obamp(7:9)   / -162.5533601, -116.1077911,  101.1189923/
      data obamp(10:12) /  -67.6856209,   24.9079067,   22.5811241/
      data obamp(13:15) /  -21.1648355,  -15.6549876,   15.3936813/
      data obamp(16:18) /   14.6660938,  -11.7273029,   10.2742696/
      data obamp(19:21) /    6.4914588,    5.8539148,   -5.4872205/
      data obamp(22:24) /   -5.4290191,    5.1609570,    5.0786314/
      data obamp(25:27) /   -4.0735782,    3.7227167,    3.3971932/
      data obamp(28:30) /   -2.8347004,   -2.6550721,   -2.5717867/
      data obamp(31:33) /   -2.4712188,    2.4625410,    2.2464112/
      data obamp(34:36) /   -2.0755511,   -1.9713669,   -1.8813061/
      data obamp(37:39) /   -1.8468785,    1.8186742,    1.7601888/
      data obamp(40:42) /   -1.5428851,    1.4738838,   -1.4593669/
      data obamp(43:45) /    1.4192259,   -1.1818980,    1.1756474/
      data obamp(46:47) /   -1.1316126,    1.0896928/

      data obrate(1:3)   /  31.609974,    32.620504,    24.172203/
      data obrate(4:6)   /  31.983787,    44.828336,    30.973257/
      data obrate(7:9)   /  43.668246,    32.246691,    30.599444/
      data obrate(10:12) /  42.681324,    43.836462,    47.439436/
      data obrate(13:15) /  63.219948,    64.230478,     1.010530/
      data obrate(16:18) /   7.437771,    55.782177,     0.373813/
      data obrate(19:21) /  13.218362,    62.583231,    63.593761/
      data obrate(22:24) /  76.438310,    45.815258,     8.448301/
      data obrate(25:27) /  56.792707,    49.747842,    12.058272/
      data obrate(28:30) /  75.278220,    65.241008,    64.604291/
      data obrate(31:33) /   1.647247,     7.811584,    12.207832/
      data obrate(34:36) /  63.856665,    56.155990,    77.448840/
      data obrate(37:39) /   6.801054,    62.209418,    20.656133/
      data obrate(40:42) /  48.344406,    55.145460,    69.000539/
      data obrate(43:45) /  11.071350,    74.291298,    11.047742/
      data obrate(46:47) /   0.636717,    12.844549/

      data obphas(1:3)   / 251.9025,     280.8325,     128.3057/
      data obphas(4:6)   / 292.7252,      15.3747,     263.7951/
      data obphas(7:9)   / 308.4258,     240.0099,     222.9725/
      data obphas(10:12) / 268.7809,     316.7998,     319.6024/
      data obphas(13:15) / 143.8050,     172.7351,      28.9300/
      data obphas(16:18) / 123.5968,      20.2082,      40.8226/
      data obphas(19:21) / 123.4722,     155.6977,     184.6277/
      data obphas(22:24) / 267.2772,      55.0196,     152.5268/
      data obphas(25:27) /  49.1382,     204.6609,      56.5233/
      data obphas(28:30) / 200.3284,     201.6651,     213.5577/
      data obphas(31:33) /  17.0374,     164.4194,      94.5422/
      data obphas(34:36) / 131.9124,      61.0309,     296.2073/
      data obphas(37:39) / 135.4894,     114.8750,     247.0691/
      data obphas(40:42) / 256.6114,      32.1008,     143.6804/
      data obphas(43:45) /  16.8784,     160.6835,      27.5932/
      data obphas(46:47) / 348.1074,      82.6496/

! ecamp  = ampl for eccen/fvelp cos/sin series (arc seconds)
! ecrate = rates for eccen/fvelp cos/sin series (arc seconds/year)
! ecphas = phases for eccen/fvelp cos/sin series (degrees)

      real (kind=8) ecamp(19), ecrate(19), ecphas(19)

      data ecamp(1:3)    /   0.01860798,  0.01627522, -0.01300660/
      data ecamp(4:6)    /   0.00988829, -0.00336700,  0.00333077/
      data ecamp(7:9)    /  -0.00235400,  0.00140015,  0.00100700/
      data ecamp(10:12)  /   0.00085700,  0.00064990,  0.00059900/
      data ecamp(13:15)  /   0.00037800, -0.00033700,  0.00027600/
      data ecamp(16:18)  /   0.00018200, -0.00017400, -0.00012400/
      data ecamp(19)     /   0.00001250/

      data ecrate(1:3)   /   4.2072050,   7.3460910,  17.8572630/
      data ecrate(4:6)   /  17.2205460,  16.8467330,   5.1990790/
      data ecrate(7:9)   /  18.2310760,  26.2167580,   6.3591690/
      data ecrate(10:12) /  16.2100160,   3.0651810,  16.5838290/
      data ecrate(13:15) /  18.4939800,   6.1909530,  18.8677930/
      data ecrate(16:18) /  17.4255670,   6.1860010,  18.4174410/
      data ecrate(19)    /   0.6678630/

      data ecphas(1:3)   /  28.620089,  193.788772,  308.307024/
      data ecphas(4:6)   / 320.199637,  279.376984,   87.195000/
      data ecphas(7:9)   / 349.129677,  128.443387,  154.143880/
      data ecphas(10:12) / 291.269597,  114.860583,  332.092251/
      data ecphas(13:15) / 296.414411,  145.769910,  337.237063/
      data ecphas(16:18) / 152.092288,  126.839891,  210.667199/
      data ecphas(19)    /  72.108838/

! mvamp  = amplitudes for mvelp sine series  (arc seconds)
! mvrate = rates for mvelp sine series (arc seconds/year)
! mvphas = phases for mvelp sine series (degrees)

      real (kind=8) mvamp(78), mvrate(78), mvphas(78)

      data mvamp(1:3)    / 7391.0225890, 2555.1526947, 2022.7629188/
      data mvamp(4:6)    /-1973.6517951, 1240.2321818,  953.8679112/
      data mvamp(7:9)    / -931.7537108,  872.3795383,  606.3544732/
      data mvamp(10:12)  / -496.0274038,  456.9608039,  346.9462320/
      data mvamp(13:15)  / -305.8412902,  249.6173246, -199.1027200/
      data mvamp(16:18)  /  191.0560889, -175.2936572,  165.9068833/
      data mvamp(19:21)  /  161.1285917,  139.7878093, -133.5228399/
      data mvamp(22:24)  /  117.0673811,  104.6907281,   95.3227476/
      data mvamp(25:27)  /   86.7824524,   86.0857729,   70.5893698/
      data mvamp(28:30)  /  -69.9719343,  -62.5817473,   61.5450059/
      data mvamp(31:33)  /  -57.9364011,   57.1899832,  -57.0236109/
      data mvamp(34:36)  /  -54.2119253,   53.2834147,   52.1223575/
      data mvamp(37:39)  /  -49.0059908,  -48.3118757,  -45.4191685/
      data mvamp(40:42)  /  -42.2357920,  -34.7971099,   34.4623613/
      data mvamp(43:45)  /  -33.8356643,   33.6689362,  -31.2521586/
      data mvamp(46:48)  /  -30.8798701,   28.4640769,  -27.1960802/
      data mvamp(49:51)  /   27.0860736,  -26.3437456,   24.7253740/
      data mvamp(52:54)  /   24.6732126,   24.4272733,   24.0127327/
      data mvamp(55:57)  /   21.7150294,  -21.5375347,   18.1148363/
      data mvamp(58:60)  /  -16.9603104,  -16.1765215,   15.5567653/
      data mvamp(61:63)  /   15.4846529,   15.2150632,   14.5047426/
      data mvamp(64:66)  /  -14.3873316,   13.1351419,   12.8776311/
      data mvamp(67:69)  /   11.9867234,   11.9385578,   11.7030822/
      data mvamp(70:72)  /   11.6018181,  -11.2617293,  -10.4664199/
      data mvamp(73:75)  /   10.4333970,  -10.2377466,   10.1934446/
      data mvamp(76:78)  /  -10.1280191,   10.0289441,  -10.0034259/

      data mvrate(1:3)   /   31.609974,    32.620504,    24.172203/
      data mvrate(4:6)   /    0.636717,    31.983787,     3.138886/
      data mvrate(7:9)   /   30.973257,    44.828336,     0.991874/
      data mvrate(10:12) /    0.373813,    43.668246,    32.246691/
      data mvrate(13:15) /   30.599444,     2.147012,    10.511172/
      data mvrate(16:18) /   42.681324,    13.650058,     0.986922/
      data mvrate(19:21) /    9.874455,    13.013341,     0.262904/
      data mvrate(22:24) /    0.004952,     1.142024,    63.219948/
      data mvrate(25:27) /    0.205021,     2.151964,    64.230478/
      data mvrate(28:30) /   43.836462,    47.439436,     1.384343/
      data mvrate(31:33) /    7.437771,    18.829299,     9.500642/
      data mvrate(34:36) /    0.431696,     1.160090,    55.782177/
      data mvrate(37:39) /   12.639528,     1.155138,     0.168216/
      data mvrate(40:42) /    1.647247,    10.884985,     5.610937/
      data mvrate(43:45) /   12.658184,     1.010530,     1.983748/
      data mvrate(46:48) /   14.023871,     0.560178,     1.273434/
      data mvrate(49:51) /   12.021467,    62.583231,    63.593761/
      data mvrate(52:54) /   76.438310,     4.280910,    13.218362/
      data mvrate(55:57) /   17.818769,     8.359495,    56.792707/
      data mvrate(58:60) /    8.448301,     1.978796,     8.863925/
      data mvrate(61:63) /    0.186365,     8.996212,     6.771027/
      data mvrate(64:66) /   45.815258,    12.002811,    75.278220/
      data mvrate(67:69) /   65.241008,    18.870667,    22.009553/
      data mvrate(70:72) /   64.604291,    11.498094,     0.578834/
      data mvrate(73:75) /    9.237738,    49.747842,     2.147012/
      data mvrate(76:78) /    1.196895,     2.133898,     0.173168/

      data mvphas(1:3)   /  251.9025,     280.8325,     128.3057/
      data mvphas(4:6)   /  348.1074,     292.7252,     165.1686/
      data mvphas(7:9)   /  263.7951,      15.3747,      58.5749/
      data mvphas(10:12) /   40.8226,     308.4258,     240.0099/
      data mvphas(13:15) /  222.9725,     106.5937,     114.5182/
      data mvphas(16:18) /  268.7809,     279.6869,      39.6448/
      data mvphas(19:21) /  126.4108,     291.5795,     307.2848/
      data mvphas(22:24) /   18.9300,     273.7596,     143.8050/
      data mvphas(25:27) /  191.8927,     125.5237,     172.7351/
      data mvphas(28:30) /  316.7998,     319.6024,      69.7526/
      data mvphas(31:33) /  123.5968,     217.6432,      85.5882/
      data mvphas(34:36) /  156.2147,      66.9489,      20.2082/
      data mvphas(37:39) /  250.7568,      48.0188,       8.3739/
      data mvphas(40:42) /   17.0374,     155.3409,      94.1709/
      data mvphas(43:45) /  221.1120,      28.9300,     117.1498/
      data mvphas(46:48) /  320.5095,     262.3602,     336.2148/
      data mvphas(49:51) /  233.0046,     155.6977,     184.6277/
      data mvphas(52:54) /  267.2772,      78.9281,     123.4722/
      data mvphas(55:57) /  188.7132,     180.1364,      49.1382/
      data mvphas(58:60) /  152.5268,      98.2198,      97.4808/
      data mvphas(61:63) /  221.5376,     168.2438,     161.1199/
      data mvphas(64:66) /   55.0196,     262.6495,     200.3284/
      data mvphas(67:69) /  201.6651,     294.6547,      99.8233/
      data mvphas(70:72) /  213.5577,     154.1631,     232.7153/
      data mvphas(73:75) /  138.3034,     204.6609,     106.5938/
      data mvphas(76:78) /  250.4676,     332.3345,      27.3039/

! i       = Index for series summations
! obsum   = Obliquity series summation
! cossum  = cos series summation for eccentricity/fvelp
! sinsum  = Sin series summation for eccentricity/fvelp
! fvelp   = Fixed vernal equinox long of perihelion
! mvsum   = mvelp series summation
! beta    = Intermediate argument for lambm0
! ryear   = year relative to 1950 (1950 = 0)
! eccen2  = eccentricity squared
! eccen3  = eccentricity cubed
! mvelpp  = moving vernal equinox long of perihelion plus pi (radians)
! pi      = pi
! psecdeg = arc sec to deg conversion
! degrad  = degree to radian conversion factor

      integer i

      real (kind=8) obsum, cossum, sinsum, fvelp, mvsum, beta, ryear
      real (kind=8) eccen2, eccen3, mvelpp, pi, psecdeg, degrad

      pi = 4.*atan(1.)
      psecdeg = 1./3600.
      degrad = pi/180.

# if !defined O_orbit_user
! The following calculates the earths obliquity, orbital eccentricity
! and vernal equinox mean longitude of perihelion using constants
! given in the program of:

! Berger, Andre.  1978  A Simple Algorithm to Compute Long-Term
! Variations of Daily Insolation.  Contribution 18, Institute of
! Astronomy and Geophysics, Universite Catholique de Louvain,
! Louvain-la-Neuve, Belgium.

! and formulae given in the paper (where less precise constants are also
! given):

! Berger, Andre.  1978.  Long-Term Variations of Daily Insolation and
! Quaternary Climatic Changes.  J. of the Atmo. Sci. 35:2362-2367

! The algorithm is valid only to 1,000,000 years past or hence.
! For a solution valid to 5-10 million years past see the above author.
! Algorithm below is better for years closer to present than is the
! 5-10 million year solution.

! In the summations below, cosine or sine arguments, which end up in
! degrees, must be converted to radians via multiplication by degrad.

! Summation of cosine series for obliquity (epsilon in Berger 1978) in
! degrees. Convert the amplitudes and rates, which are in arc secs,
! into degrees via multiplication by psecdeg (arc seconds to degrees
! conversion factor).  For obliq, first term is Berger 1978 epsilon
! star; second term is series summation in degrees.

      ryear = year - 1950.
      if ( abs(ryear) .gt. 1000000.0 )then
        print*, 'Error ==> orbit only valid for ryear+-1000000'
        stop
      endif

      obsum = 0.0
      do i=1,47
        obsum = obsum + obamp(i)*psecdeg*cos((obrate(i)*psecdeg*ryear
     &        + obphas(i))*degrad)
      enddo
      obliq = 23.320556 + obsum

! Summation of cosine and sine series for computation of eccentricity
! (eccen; e in Berger 1978) and fixed vernal equinox longitude of
! perihelion (fvelp; pi in Berger 1978), which is used for computation
! of moving vernal equinox longitude of perihelion.  Convert the rates,
! which are in arc seconds, into degrees via multiplication by psecdeg.

      cossum = 0.0
      do i=1,19
        cossum = cossum + ecamp(i)*cos((ecrate(i)*psecdeg*ryear
     &         + ecphas(i))*degrad)
      enddo

      sinsum = 0.0
      do i=1,19
        sinsum = sinsum+ecamp(i)*sin((ecrate(i)*psecdeg*ryear
     &         + ecphas(i))*degrad)
      enddo

! Use summations to calculate eccentricity

      eccen2 = cossum*cossum + sinsum*sinsum
      eccen  = sqrt(eccen2)
      eccen3 = eccen2*eccen

! A series of cases for fvelp, which is in radians.

      if (abs(cossum) .le. 1.0e-8) then
        if (sinsum .eq. 0.0) then
          fvelp = 0.0
        elseif (sinsum .lt. 0.0) then
          fvelp = 1.5*pi
        elseif (sinsum .gt. 0.0) then
          fvelp = .5*pi
        endif
      elseif (cossum .lt. 0.0) then
        fvelp = atan(sinsum/cossum) + pi
      elseif (cossum .gt. 0.0) then
        if (sinsum .lt. 0.0) then
          fvelp = atan(sinsum/cossum) + 2.*pi
        else
          fvelp = atan(sinsum/cossum)
        endif
      endif

! Summation of sin series for computation of moving vernal equinox long
! of perihelion (mvelp; omega bar in Berger 1978) in degrees.  For
! mvelp, first term is fvelp in degrees; second term is Berger 1978 psi
! bar times years and in degrees; third term is Berger 1978 zeta; fourth
! term is series summation in degrees.  Convert the amplitudes and
! rates, which are in arc seconds, into degrees via multiplication by
! psecdeg. Series summation plus second and third terms constitute
! Berger 1978 psi, which is the general precession.

      mvsum = 0.0
      do i=1,78
        mvsum = mvsum + mvamp(i)*psecdeg*sin((mvrate(i)*psecdeg*ryear
     &        + mvphas(i))*degrad)
      enddo
      mvelp = fvelp/degrad + 50.439273*psecdeg*ryear + 3.392506 + mvsum

      ! Cases to make sure mvelp is between 0 and 360.

      do while (mvelp .lt. 0.0)
        mvelp = mvelp + 360.0
      enddo
      do while (mvelp .ge. 360.0)
        mvelp = mvelp - 360.0
      enddo

# else
      eccen2 = eccen*eccen
      eccen3 = eccen2*eccen

# endif
! 180 degrees must be added to mvelp since observations are made from
! the earth and the sun is considered (wrongly for the algorithm) to go
! around the earth. For a more graphic explanation see Appendix B in:

! A. Berger, M. Loutre and C. Tricot. 1993.  Insolation and Earth
! Orbital Periods.  J. of Geophysical Research 98:10,341-10,362.

! So mvelp becomes mvelpp (mvelp plus pi)

      mvelpp = (mvelp + 180.)*degrad

! Set up an argument used several times in lambm0 calculation ahead.

      beta = sqrt(1. - eccen2)

! The mean longitude at the vernal equinox (lambda m nought in Berger
! 1978; in radians) is calculated from the following formula given in
! Berger 1978.  At the vernal equinox the true longitude (lambda in
! Berger 1978) is 0.

      lambm0 = 2.*((.5*eccen + .125*eccen3)*(1. + beta)*sin(mvelpp)
     &       - .250*eccen2*(.5    + beta)*sin(2.*mvelpp)
     &       + .125*eccen3*(1./3. + beta)*sin(3.*mvelpp))
# endif

      return
      end
