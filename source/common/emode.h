!======================= include file "emode.h" ========================

!     variables for  external mode

!     psi   = stream function (,,1) is for tau; (,,2) is for tau-1
!     zu    = vertically averaged forcing from momentum equations
!             (,,1) is zonal and (,,2) is meridional component
!     ztd   = curl of "zu" for the stream function equation
!     ptd   = time change of stream function
!     h     = depth over "u" points
!     hr    = reciprocal depth over "u" points
!     res   = residual from elliptic solver

!     map   = land mass map distinguishing, ocean, land, and perimeters

!     mxscan  = max number of allowable scans for Poisson solvers
!     mscan   = actual number of scans taken by Poisson solvers
!     sor     = successive over-relaxation constant
!     tolrsf  = tolerance for stream function calculation.
!               the solution is halted when it is within "tolrsf"
!               of the "true" solution assuming geometric convergence.
!     tolrsp  = tolerance for surface pressure calculation
!               the solution is halted when it is within "tolrsp"
!               of the "true" solution assuming geometric convergence.
!     tolrfs  = tolerance for implicit free surface calculation
!               the solution is halted when it is within "tolrfs"
!               of the "true" solution assuming geometric convergence.
!     esterr  = estimated maximum error in elliptic solver assuming
!               geometric convergence

!     nisle = number of land masses
!     nippts= number of land mass perimeter points
!     iperm = "i" coordinate for the land mass perimeter point
!     jperm = "j" coordinate for the land mass perimeter point
!     iofs  = offset for indexing into the land mass perimeter points
!     imask = controls whether calculations get done on perimeters
!     set mask for land mass perimeters on which to perform calculations
!     imask(-n) = .false.  [no equations ever on dry land mass n]
!     imask(0)  = .true.   [equations at all mid ocean points]
!     imask(n)  = .true./.false [controls whether there will be
!                                equations on the ocean perimeter of
!                                land mass n]
!     note: land mass 1 is the northwest-most land mass
!     for the numbering of the other landmasses, see generated map(i,j)
#if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface

!     surface pressure/free surface time centring parameters:

!     alph,gam = parameters to centre the Coriolis and surface pressure
!                gradients in time (leapfrog)

!     theta    =parameters to centre the Coriolis and surface pressure
!                gradients in time (mixing ts )

!     apgr     = is = to alpha/theta leapfrog/mixing ts

# if defined O_implicit_free_surface

!     to set free surface time centring parameters...
!     note: for proper time centring of Coriolis
!     and pressure gradients alph, gam should
!     satisfy gam = 1 -2*alph.  for stability
!     alph should be > 1/4.  recommended values
!     are alph = gam = 1/3.
# endif

!     ps    = surface pressure (,,1) is for tau; (,,2) is for tau-1
!     divf  = barotropic divergence of uncorrected ubar & vbar (rhs of
!             surface pressure eqn)
!     ubar  = barotropic velocity defined on "u" point for "tau"
!             (,,1) is zonal and (,,2) is meridional velocity
!     ubarm1= barotropic velocity defined on "u" point for "tau-1"
!             (,,1) is zonal and (,,2) is meridional velocity
#endif

      character(16) :: variable
      common /emode_c/ variable

      logical converged, imask
      common /emode_l/ converged, imask (-mnisle:mnisle)

      integer mxscan, mscan, map, nippts, iofs, iperm
      integer jperm, nisle, imain

      common /emode_i/ mxscan, mscan
      common /emode_i/ map(imt,jmt)
      common /emode_i/ nippts(mnisle), iofs(mnisle), iperm(maxipp)
      common /emode_i/ jperm(maxipp), nisle, imain

      real tolrsf, tolrsp, tolrfs, sor, esterr, ptd, res, hr
      real h, zu, psi, ztd, alph, gam, theta, apgr, uhat, pguess
      real ps, divf, ubar, ubarm1

      common /emode_r/ tolrsf, tolrsp, tolrfs, sor, esterr
      common /emode_r/ ptd(imt,jmt), res(imt,jmt), hr(imt,jmt)
      common /emode_r/ h(imt,jmt), zu(imt,jmt,2)
#if defined O_stream_function
      common /emode_r/ psi(imt,jmt,2), ztd(imt,jmt)
#endif
#if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface
      common /emode_r/ alph, gam, theta, apgr
      common /emode_r/ uhat(imt,jmt,2), pguess(imt,jmt)
      common /emode_r/ ps(imt,jmt,2)
      common /emode_r/ divf(imt,jmt), ubar(imt,jmt,2), ubarm1(imt,jmt,2)
#endif
