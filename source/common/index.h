!======================= include file "index.h" ========================

!    various starting & ending indices for controlling where quantities
!    are filtered. this removes the time step restriction due to
!    convergence of meridians.

#if defined O_fourfil || defined O_firfil
!    istf  = starting index for filtering "t" grid points
!    ietf  = ending index for filtering "t" grid points
!    isuf  = starting index for filtering "u,v" grid points
!    ieuf  = ending index for filtering "u,v" grid points
!    iszf  = starting index for filtering "vorticity" grid points
!            (also for divergence when using implicit free surface)
!    iezf  = ending index for filtering "vorticity" grid points
!            (also for divergence when using implicit free surface)

!     define latitudinal domain for filtering

!     filter t to yield equiv dx at lat rjft0 from lat rjfrst to rjft1
!       and lat rjft2 to yt(jmtm1)
!     filter u to yield equiv dx at lat rjfu0 from lat rjfrst to rjfu1
!       and lat rjfu2 to yt(jmtm2)

!     lsegf  = max number of longitudinal strips per filtering latitude
!     jmtfil = max number of latitudes to be filtered
!     numflt = specifies the number of filter applications to use
!              for the tracers, vorticity, or divergence when using the
!               finite impulse response filter
!     numflu = specifies the number of filter applications to use
!              for velocities when using the finite impulse response
!               filter

      integer lsegf, jmtfil
      parameter (lsegf=20, jmtfil=50)

      integer jfrst, jft0, jft1, jft2, jfu0, jfu1, jfu2
      integer jskpt, jskpu, njtbft, njtbfu, numflt, numflu
      integer istf, ietf, isuf, ieuf, iszf, iezf

      common /index_i/ jfrst, jft0, jft1, jft2, jfu0, jfu1, jfu2
      common /index_i/ jskpt, jskpu, njtbft, njtbfu
# if defined O_firfil
      common /index_i/ numflt(jmtfil), numflu(jmtfil)
# endif
      common /index_i/ istf(jmtfil,lsegf,km), ietf(jmtfil,lsegf,km)
      common /index_i/ isuf(jmtfil,lsegf,km), ieuf(jmtfil,lsegf,km)
      common /index_i/ iszf(jmtfil,lsegf),    iezf(jmtfil,lsegf)
#endif

      real rjfrst, rjft0, rjft1, rjft2, rjfu0, rjfu1, rjfu2
      common /index_r/ rjfrst, rjft0, rjft1, rjft2, rjfu0, rjfu1, rjfu2
