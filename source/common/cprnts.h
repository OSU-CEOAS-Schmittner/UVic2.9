!====================== include file "cprnts.h" =======================

!     variables used for controlling matrix printouts during diagnostic
!     timesteps:

!     nlatpr = maximum number of latitudes for matrix printouts
!              on diagnostic time steps
!     prlat  = latitudes (deg) at which (x,z) printouts are desired

!     start & end coordinates for matrix printouts of (x,z) sections

!     prslon = starting longitudes (deg)
!     prelon = ending longitudes (deg)
!     prsdpt = starting depths  (cm)
!     predpt = ending depths  (cm)

!     start & end coordinates for matrix printouts of (x,y) sections

!     slonxy = starting longitude (deg)
!     elonxy = ending longitude (deg)
!     slatxy = starting latitude (deg)
!     elatxy = ending latitude (deg)

!     matrix printouts of (y,z) sections will use above coordinates

      real prlat, prslon, prelon, prsdpt, predpt, slatxy, elatxy
      real slonxy, elonxy

      common /cprnts/ prlat(nlatpr), prslon(nlatpr), prelon(nlatpr)
      common /cprnts/ prsdpt(nlatpr), predpt(nlatpr)
      common /cprnts/ slatxy, elatxy, slonxy, elonxy
