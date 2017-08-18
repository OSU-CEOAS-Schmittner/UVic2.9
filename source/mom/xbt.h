!====================== include file "xbt.h" ===========================

!     Each XBT station is located at latitude "xbtlat" and longitude
!     "xbtlon". Data is collected at each grid point from the first
!     level down through the nearest model level corresponding to a
!     depth of "xbtdpt" cm. Actually, all coordinates are converted to
!     the nearest model temperature grid point.

!     All basic quantities as well as all terms in the momentum,
!     temperature, and salinity equations are averaged over the time
!     period specified by "xbtint".

!     The maximum number of XBTs may be increased by changing parameter
!     "maxxbt" below.

!     inputs:

!     maxxbt  = maximum number of XBTs allowed.
!     kmxbt   = maximum number of levels from surface downward (<=km)
!               set kmxbt  < km to save space
!     xbtlat  = real latitude of XBTs in degrees
!     xbtlon  = real longitude of XBTs in degrees
!     xbtdpt  = real depths of XBTs in cm
!     items   = number of items in the XBT
!     xname   = character*12 names of XBT quantities

!     outputs:

!     numxbt  = actual number of XBTs used
!     nxbtts  = current number of time steps in accumulated XTB data

!     ixbt    = longitude index of nearest model temperature grid point
!               corresponding to "xbtlon"
!     jxbt    = latitude index of nearest model temperature grid point
!               corresponding to "xbtlat"
!     kxbt    = depth index of nearest model temperature grid point
!               corresponding to "xbtdpt"
!     nsxbt   = starting number for the XBTs on each latitude
!     nexbt   = ending number for the XBTs on each latitude

!     txbt    = accumulator array for time rate of change of
!                tracers. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (10) where each
!                term has the units of "tracer units/sec" using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the tracer
!                (2)  = change due to zonal nonlinear term (UT)x
!                (3)  = change due to meridional nonlinear term (VT)y
!                (4)  = change due to vertical nonlinear term (WT)z
!                (5)  = change due to zonal diffusion: Ah*Txx
!                (6)  = change due to meridional diffusion: Ah*Tyy
!                (7)  = change due to vertical diffusion:  kappa_h*Tzz
!                (8)  = change due to source term
!                (9)  = change due to explicit convection
!                (10) = change due to filtering
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of tracer "A" is -U(A)x = A(Ux) - (UA)x
!                (11) = zonal advection U(Ax)
!                (12) = meridional advection V(Ay)
!                (13) = vertical advection W(Az)
!                (14) = change of tracer variance
!                (15) = average tracer within volume (tracer units)
!     uxbt    = accumulator array for time rate of change of
!                momentum. the total time rate of change
!                is broken down into components as follows:
!                the form is d( )/dt = terms (2) ... (13) where each
!                term has the units of "cm/sec**2" and "Q" is the
!                momentum component {zonal or meridional} using
!                schematic terms for illustration.
!                (1)  = total time rate of change for the momentum
!                (2)  = change due to the pressure gradient: grad_p
!                       without the surface pressure gradients
!                       (i.e., for computing the internal modes)
!                (3)  = change due to zonal nonlinear term: (UQ)x
!                (4)  = change due to meridional nonlinear term: (VQ)y
!                (5)  = change due to vertical nonlinear term: (wQ)z
!                (6)  = change due to zonal viscosity: Am*Qxx
!                (7)  = change due to meridional viscosity: Am*Qyy
!                (8)  = change due to vertical viscosity: kappa_m*Qzz
!                (9)  = change due to metric terms
!                (10) = change due to coriolis terms: fQ
!                (11) = change due to source terms
!                (12) = change due to surface pressure gradient
!                       this is obtained after solving the external mode
!                       in the stream function technique. It is solved
!                       directly from the elliptic equation for the
!                       prognostic surface pressure technique
!                (13) = change due to metric advection
!     the nonlinear terms can be broken into two parts: advection and a
!     continuity part: The physically meaningful part is advection.
!     eg: Zonal advection of vel component "Q" is -U(Q)x = Q(U)x - (UQ)x
!                (14) = zonal advection U(Qx)
!                (15) = meridional advection V(Qy)
!                (16) = vertical advection W(Qz)
!                (17) = average velocity component
!     xbtw    = accumulator array for vertical velocity. (cm/sec)
!               this is the average of adv_vbu at top and bottom of cell
!     txbtsf  = accumulator array for tracer surface flux terms.
!                tracer (#1,#2) units = (cal/cm**2/sec, gm/cm**2/sec)
!     uxbtsf  = accumulator array for wind stress terms. (dynes/cm**2)

!     ntxbt   = number of terms for tracers
!     nuxbt   = number of terms for velocity

      integer maxxbt, kmxbt, ntxbt, nuxbt
      parameter (maxxbt=3, kmxbt=4)
      parameter (ntxbt=15, nuxbt=17)

      integer nxbtts, numxbt, nsxbt, nexbt, ixbt, jxbt, kxbt
      common /cxbt_i/ nxbtts, numxbt, nsxbt(jmt), nexbt(jmt)
      common /cxbt_i/ ixbt(maxxbt), jxbt(maxxbt), kxbt(maxxbt)

      character(12) :: xnamet, xnameu, xnamex
      common /cxbt_c/ xnamet(ntxbt), xnameu(nuxbt,2), xnamex(4)

      real xbtlat, xbtlon, xbtdpt, txbt, txbtsf, uxbt, uxbtsf, xbtw
      common /cxbt_r/ xbtlat(maxxbt), xbtlon(maxxbt), xbtdpt(maxxbt)
      common /cxbt_r/ txbt(kmxbt,ntxbt,nt,maxxbt), txbtsf(nt,maxxbt)
      common /cxbt_r/ uxbt(kmxbt,nuxbt,2,maxxbt), uxbtsf(2,maxxbt)
      common /cxbt_r/ xbtw(kmxbt,maxxbt)
