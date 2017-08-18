!======================= include file "accel.h" ========================

!     depth dependent tracer timestep acceleration multipliers used to
!     hasten the convergence to equilibrium of the deeper portions of
!     ocean-climate models.

!     accelerate abyssal processes by varying the length of the tracer
!     timestep with depth.  by using longer timesteps at depth, one can
!     in effect reduce the heat capacity of the deeper levels and speed
!     convergence to equilibrium.
!     note:
!     by applying this method, one is assuming that there is a single
!     steady-state solution to the model being considered.
!     also, since the diagnostic timestep calculations of "termbt" do
!     not attempt to account for depth variant timestep lengths, the
!     truncation error reported will increase, because it will include
!     the tracer changes due to variations in "dtxcel".

!     reference:
!       Bryan, K., 1984: accelerating the convergence to equilibrium
!     of ocean climate models, J. Phys. Oceanogr., 14, 666-673.
!     ("dtxcel" here is the same as 1/gamma in the above reference)
!     set "dtxcel" to 1.0 at the surface and for upper levels not
!     to be accelerated
!     set "dtxcel" to values greater than 1.0 at deeper levels to
!     accelerate convergence if above requirements are met

!     dtxcel   = model level dependent tracer timestep multipliers
#if defined O_isopycmix
!     dtxsqr   = square root of "dtxcel" (used in computation of
!                maximum slope constraint for isopycnal mixing)
#endif
#if !defined O_implicitvmix || defined O_isopycmix
!     dztxcl   = layer thickness divided by the timestep multiplier
!                (needed for convection code)
!     dzwxcl   = multiplication factor relating to the vertical
!                distance between ts points, scaled according
!                to timestep multipliers for use in convection code
#endif

      real dtxcel, dtxsqr, dztxcl, dzwxcl

      common /accel/ dtxcel(km)
#if defined O_isopycmix
      common /accel/ dtxsqr(km)
#endif
      common /accel/ dztxcl(km), dzwxcl(km)
