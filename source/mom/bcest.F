      subroutine bcest (tlat, ulat, bcval)

#if defined O_mom
!=======================================================================
!     this routine linearly interpolates global, zonal mean values of
!     ocean model surface boundary conditions (SST, salinity, WSX, WSY)
!     having 4.5 deg spacing, to the requested latitude.
!=======================================================================

      implicit none

      integer nbc, nolat, nolatp, n, nn, i
      parameter (nbc=4, nolat=40, nolatp=nolat+1)

      real c1, c0, p5, phione, phitwo, phithree, phifour, pi, radulat
      real ulat, radtlat, tlat, dolat, rdolat, olatt, olatv, ylatt, d
      real ylatv
      parameter (c1=1.0, c0=0.0, p5=0.5)

      common /cbcest/ olatt(nolat), olatv(nolatp), dolat, rdolat
      real sstobs(nolat), salobs(nolat)
      real wsxobs(nolatp), wsyobs(nolatp)
      real bcval(nbc)

!     bcval  = estimated boundary condition values ( wsx, wsy,t,s)
!       bcval(1) and bcval(2) units = dynes per square centimetre
!       bcval(3) units = degrees C
!       bcval(4) units = parts per thousand
!     nbc    = number of boundary conditions
!     olatt  = latitude points for observed data
!     olatv  = latitude points for observed data
!     dolat  = latitude spacing for observed data
!     ylatt  = latitude where t,s boundary conditions are desired
!     ylatv  = latitude where windstress boundary conditions are desired

!     "observed" temperature and salinity data are based on global,
!     annual mean zonally averaged values from the Levitus Atlas (1982).
!     "observed" windstress data are based on global, annual mean,
!     zonally averaged values from Hellerman and Rosenstein (1981).
!     some smoothing was done.

!     references:
!       Hellerman, S, and M. Rosenstein, normal monthly wind stress
!     over the world ocean with error estimates, J. Phys, Oceanogr., 13,
!     1093-1104,1983.
!       Levitus, S., Climatological atlas of the world ocean, NOAA
!     Prof. Paper 13, US Gov`t printing Office, Washington, DC, 1982.

      data sstobs / -1.75, -1.75, -1.50, -1.50, -1.28,
     &              -0.55,  0.90,  2.92,  5.45,  8.62,
     &              12.27, 15.49, 18.30, 20.67, 22.64,
     &              24.14, 25.27, 26.37, 26.52, 26.16,
     &              26.85, 27.27, 26.82, 26.42, 25.53,
     &              24.03, 22.07, 19.73, 17.02, 12.77,
     &               8.93,  7.25,  6.22,  4.67,  4.57,
     &               3.03, -0.01, -1.05, -1.75, -1.75/

      data salobs / 34.30, 34.30, 34.30, 34.13, 33.98,
     &              33.97, 33.97, 33.98, 34.03, 34.24,
     &              34.61, 35.02, 35.37, 35.61, 35.72,
     &              35.68, 35.51, 35.22, 35.05, 35.12,
     &              34.80, 34.56, 34.71, 34.90, 35.27,
     &              35.67, 35.56, 35.49, 35.23, 34.28,
     &              33.57, 33.57, 33.60, 33.80, 34.04,
     &              34.05, 32.65, 32.30, 32.10, 32.00/

      data wsxobs /  0.00,
     &               0.00,  0.00, -0.02,  0.15,  0.31,
     &               0.50,  0.82,  1.08,  1.23,  1.16,
     &               0.84,  0.41,  0.02, -0.35, -0.55,
     &              -0.67, -0.64, -0.46, -0.29, -0.19,
     &              -0.16, -0.33, -0.52, -0.59, -0.55,
     &              -0.32,  0.09,  0.42,  0.56,  0.76,
     &               0.81,  0.65,  0.29,  0.06, -0.10,
     &              -0.05, -0.03,  0.05,  0.10,  0.01/

      data wsyobs /  .000,
     &               .000,  .009,  .032,  .005, -.023,
     &              -.075, -.155, -.202, -.230, -.179,
     &              -.049,  .093,  .214,  .294,  .344,
     &               .383,  .364,  .269,  .189,  .178,
     &               .125, -.122, -.213, -.251, -.259,
     &              -.202, -.189, -.179, -.183, -.009,
     &               .023,  .053, -.048, -.185, -.225,
     &              -.097, -.050, -.023, -.006,  .000/

!---------------------------------------------------------------------
!     set latitudes of sst and salinity observations
!     and set latitudes of windstress observations
!---------------------------------------------------------------------

      dolat = 180.0/nolat
      rdolat = c1/dolat
      do n=1,nolat
        olatt(n) = -90.0 + (n-p5)*dolat
        olatv(n) = -90.0 + (n-1.0)*dolat
      enddo
      olatv(nolatp) = -90.0 + (nolat)*dolat

!---------------------------------------------------------------------
!   use linear interpolation to produce the estimated surface boundary
!   condition values for temperature and salinity at t,s row j
!---------------------------------------------------------------------

      ylatt = tlat

      if (ylatt .le. olatt(1)) then
        nn = 1
        d = c0
      elseif (ylatt .ge. olatt(nolat)) then
        nn = nolat-1
        d = dolat
      else
        do i=2,nolat
          if (ylatt .le. olatt(i)) then
            nn = i - 1
            d  = ylatt - olatt(nn)
            goto 201
          endif
        enddo
      endif

201   continue
      bcval(3) = (sstobs(nn)*(dolat - d) + sstobs(nn+1)*d)*rdolat
      bcval(4) = (salobs(nn)*(dolat - d) + salobs(nn+1)*d)*rdolat

!---------------------------------------------------------------------
!   use linear interpolation to produce the estimated surface boundary
!   condition values for wind stress components at u,v row j
!---------------------------------------------------------------------

      ylatv = ulat

      if (ylatv .le. olatv(1)) then
        nn = 1
        d = c0
      elseif (ylatv .ge. olatv(nolatp)) then
        nn = nolatp - 1
        d = dolat
      else
        do i=2,nolatp
          if (ylatv .le. olatv(i)) then
            nn = i - 1
            d  = ylatv - olatv(nn)
            goto 301
          endif
        enddo
      endif
301   continue
      bcval(1) = (wsxobs(nn)*(dolat - d) + wsxobs(nn+1)*d) *rdolat
      bcval(2) = (wsyobs(nn)*(dolat - d) + wsyobs(nn+1)*d) *rdolat
#endif

      return
      end
