      subroutine mom

#if defined O_mom
!=======================================================================

!                    GFDL Modular Ocean Model

!      A primitive equation ocean model developed by researchers at
!      the Geophysical Fluid Dynamics Laboratory /NOAA in
!      Princeton, NJ. 08542.

!      The model is based on the pioneering work of

!      Kirk Bryan: A numerical method for the study of the  of the
!      circulation world ocean: 1969, J. Computat. Phys 4 347-376

!                              and

!      the invaluable work of Mike Cox & Bert Semtner on earlier
!                    fortran implementations.

!      The GFDL Modular Ocean Model (acronym MOM) is a three
!      dimensional primitive equation ocean model intended  to be
!      a flexible tool useful for ocean and coupled air-sea modeling
!      applications over a wide range of space & time scales.
!      It is also intended to run efficiently on scalar and vector
!      architectures. The programming approach is modular and
!      additions to this model are encouraged to follow this
!      approach. Additional modules will be added with time and
!      new versions will be released when ready.

!      Documentation:

!      For documentation refer to the postscript manual which is
!      included along with this code.

!      Requirements:

!      Standard fortran 77 is used (except for namelist which is
!      fortran 90 compliant, do enddo, use of "max" function in
!      parameter statements and variable names > than 6 characters)
!      The preprocessor "cpp" (available on systems using "c" or UNIX)

!      Disclaimer:

!      MOM is an ocean modeling research tool developed at GFDL.
!      Others may use it freely but we assume no responsibility
!      for problems or incorrect use of MOM. It is left to the user to
!      satisfy (him/her)self that a particular configuration is
!      working correctly. To this end, many of the included
!      diagnostics will be helpful.

!=======================================================================

      implicit none

      character (120) :: fname

      integer jrow, i, num_mw, mw, js, joff, je, is, ie, jscalc
      integer jecalc, jstrac, jetrac, ntaux

      logical first_mw

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "emode.h"
      include "iounit.h"
      include "mw.h"
      include "csbc.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"

!-----------------------------------------------------------------------
!     integrate one time step
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     update timestep counter, set corresponding model time, and set
!     time dependent logical switches which determine program flow.
!-----------------------------------------------------------------------

      itt = itt + 1
      call increment_time (dtts)
      call set_time_switches

!-----------------------------------------------------------------------
!       initialize diagnostic variables
!-----------------------------------------------------------------------

      call diagi

!-----------------------------------------------------------------------
!     adjust various quantities for leapfrog/mixing timesteps

!     leapfrog----------> h(tau+1) = h(tau-1) + 2dt*F(tau)

!     forward-----------> tau-1 <= tau
!                         h(tau+1) = h(tau-1) + dt*F(tau)

!     euler backward:     tau-1 <= tau
!       euler1----------> h(tau` ) = h(tau-1) + dt*F(tau)
!       euler2----------> h(tau+1) = h(tau-1) + dt*F(tau`)
!-----------------------------------------------------------------------

      if (leapfrog) then

!       normal leapfrog time step

        nots = nots + 1
        euler1  = .false.
        euler2  = .false.
        forward = .false.
        eots    = .true.

        c2dtts  = c2*dtts
        c2dtuv  = c2*dtuv
        c2dtsf  = c2*dtsf
      else
        nots = 1
!       mixing time step (forward step or euler backward step)

        if (eb) then
          euler1  = .true.
          euler2  = .false.
          forward = .false.
          eots    = .false.
        else
          euler1  = .false.
          euler2  = .false.
          forward = .true.
          eots    = .true.
        endif

        c2dtts = dtts
        c2dtuv = dtuv
        c2dtsf = dtsf
# if defined O_stream_function

        do jrow=1,jmt
          do i=1,imt
            psi(i,jrow,2) = psi(i,jrow,1)
          enddo
        enddo
# endif
      endif
# if defined O_sed

!     determine if this is a flux summing time step
      addflxo = .true.
#  if defined O_even_fluxes
      if (mod(nots,2) .ne. 0) addflxo = .false.
#  endif
# endif
# if defined O_rigid_lid_surface_pressure || defined O_implicit_free_surface

      if (leapfrog) then
        apgr = alph
      else
        apgr = theta
        do jrow=1,jmt
          do i=1,imt
            ubarm1(i,jrow,1) = ubar(i,jrow,1)
            ubarm1(i,jrow,2) = ubar(i,jrow,2)
          enddo
        enddo
      endif
# endif
# if defined O_implicit_free_surface

      if (leapfrog) then
        if (acor .ne. c0) acor = alph
      else
        if (acor .ne. c0) acor = theta
      endif
# endif

!-----------------------------------------------------------------------
!     set time centering "gcor" for coriolis term
!-----------------------------------------------------------------------

      if (acor .eq. c0) then
        gcor = c1
      elseif (acor .ne. c0) then
        gcor = c0
      endif
# if defined O_implicit_free_surface

!     set "gcor" for leapfrog or mixing time steps

      if (leapfrog) then
        if (acor .eq. c0) then
          gcor = c1
        else
          gcor = gam
        endif
      else
        gcor = c0
      endif
# endif

!-----------------------------------------------------------------------
!     update pointers to tau-1, tau, & tau+1 data on disk.
!     for latitude rows they point to latdisk
!     for 2D fields they point to records on kflds
!-----------------------------------------------------------------------

      taum1disk = mod(itt+1,2) + 1
      taudisk   = mod(itt  ,2) + 1
      taup1disk = taum1disk

!-----------------------------------------------------------------------
!     update pointers (indices) to tau-1, tau, & tau+1 data in the MW
!-----------------------------------------------------------------------

      if (wide_open_mw) then

!       rotate time levels instead of moving data

        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      else

!       they are being held constant in time.

      endif

!=======================================================================

!             SOLVE THE BAROCLINIC AND TRACER EQUATIONS

!     Since all latitude rows may not fit into central memory, a
!     flexible MW (memory window) approach is used. The minimum MW
!     holds 3 latitude rows and the maximum MW holds "jmt" latitude
!     rows in central memory. Choose the size to fit into available
!     central memory. The MW is loaded with variables from disk as
!     many times as needed to solve latitude rows 2 through "jmt-2".

!     Example using a MW with 3 rows (jmw=3)

!     "loadmw" loads variables from the first 3 latitude jrows into
!     rows js=1 through je=3 in the 1st MW (mw=1). Equations are
!     computed for j=2 in the MW (corresponding to latitude jrow=2)
!     then written to disk. For the second MW (mw=2), "loadmw"
!     first copies variables from j=2 to j=1, then variables from
!     j=3 to j=2 in the MW, before loading latitude jrow=4 variables
!     into row js=je=3 in the MW. Equations are computed for j=2 in
!     the MW (corresponding to latitude jrow=3) then written to disk. The
!     process continues until latitude jrows 2 through jmt-1 are
!     computed.

!     Example using a MW with 5 rows (jmw=5)

!     "loadmw" loads variables from the first 5 latitude jrows into
!     rows js=1 through je=5 in the 1st MW (mw=1). Equations are
!     computed for j=2,3,4 in the MW (latitude jrows=2,3,4)
!     and written to disk. For the second MW (mw=2),
!     "loadmw" first moves variables from j=2 to j=1, then moves
!     variables from j=3 to j=2 in the MW, before loading jrow 6,7,8
!     variables into rows js=3 to je=5 in the MW. Equations
!     are computed for j=2,3,4 in the MW (latitude jrows=5,6,7)
!     then written to disk. The process continues until latitude
!     jrows 2 through jmt-1 are computed. Note that "je" for the last
!     MW may be less than "jmw" (depending on the size of
!     "jmw" and "jmt") and there may be fewer than 3 calculated rows.
!     On the last MW, "je" will correspond to latitude row "jmt".

!     Note:

!     When the MW is fully opened (jmw=jmt), all latitude rows
!     reside in the MW (none on disk). Instead of reading/writing
!     from MW to disk, data is moved between time levels within the
!     MW.

!=======================================================================

1000    continue

!-----------------------------------------------------------------------
!     Solve equations for rows within each MW
!     MW size is controlled by parameter "jmw" in file size.h
!-----------------------------------------------------------------------

!     num_mw  = number of MW`s needed to solve rows 2 -> jmt-1
      num_mw = (jmt-2)/ncrows + (jmt-3)/(ncrows*((jmt-2)/ncrows))
      do mw = 1,num_mw
        first_mw = (mw .eq. 1)

!       define starting and ending rows for each MW

        if (first_mw) then
          js = 1
        else
          js = jmw - ncrows + 1
        endif

        joff = (mw-1)*ncrows
        je = min(jmw,jmt-joff)

        is = 2
        ie = imt - 1

!-----------------------------------------------------------------------
!       load prognostic and related variables into the MW

!       joff = offset relating row "j" in the MW to latitude "jrow"
!       js   = starting row within the MW for LOADING latitude rows
!       je   = ending row within the MW for LOADING latitude rows
!       is   = starting index for longitude
!       ie   = ending index for longitude

!       typically, for a 2nd order window (schemes needing to access
!       one cell in all directions).

!       first MW     : load latitude data into rows js=1 ... je=jmw
!       1 < MW < last: load latitude data into rows js=3 ... je=jmw
!       last MW      : load latitude data into rows js=3 ... je<=jmw
!                      On last MW, row "je" corresponds to latitude
!                      row "jmt"
!-----------------------------------------------------------------------

        call loadmw (joff, js, je, is, ie, latdisk(taum1disk)
     &,                latdisk(taudisk), first_mw)

!-----------------------------------------------------------------------
!       calculate advection velocities for momentum and tracers
!-----------------------------------------------------------------------

        call adv_vel (joff, js, je, is, ie)
# if defined O_isopycmix

!-----------------------------------------------------------------------
!       calculate isopycnal diffusion tensor components (and
!       gent_mcwilliams advective velocities) for use with tracers
!-----------------------------------------------------------------------

        call isopyc (joff, js, je, is, ie)
# endif

!-----------------------------------------------------------------------
!       set vertical mixing coefficients for momentum and tracers
!-----------------------------------------------------------------------

        call vmixc (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!       set horizontal mixing coefficients for momentum and tracers
!-----------------------------------------------------------------------

        call hmixc (joff, js, je, is, ie)
# if defined O_biharmonic

!-----------------------------------------------------------------------
!       calculate del**2 quantities for biharmonic mixing
!-----------------------------------------------------------------------

        call delsq (joff, js, je, is, ie)
# endif

!-----------------------------------------------------------------------
!       set vertical boundary conditions for momentum and tracers
!-----------------------------------------------------------------------

        call setvbc (joff, js, je, is, ie)

!-----------------------------------------------------------------------
!       set which MW rows to calculate: "jscalc" through "jecalc"
!-----------------------------------------------------------------------

        jscalc = 2
        jecalc = min(jsmw+ncrows-1,jmt-1-joff)
# if defined O_pressure_gradient_average

!       compute one extra row of tracers
        jstrac = jscalc + 1 - 1/mw
        jetrac = min(jsmw+ncrows,jmt-1-joff)
# else
        jstrac = jscalc
        jetrac = jecalc
# endif

!-----------------------------------------------------------------------
!       compute tracers and internal mode velocities
!-----------------------------------------------------------------------

        call tracer (joff, jstrac, jetrac, is, ie)
        call clinic (joff, jscalc, jecalc, is, ie)

!-----------------------------------------------------------------------
!       calculate diagnostics
!-----------------------------------------------------------------------

        call diag (joff, jscalc, jecalc, is, ie)

!-----------------------------------------------------------------------
!     write prognostic variables from the MW to disk "tau+1"
!-----------------------------------------------------------------------

        if (wide_open_mw) then
!         do nothing since variables are already in "tau+1" MW
        else
          call putmw (joff, jscalc, jecalc, latdisk(taup1disk))
        endif

      enddo

!=======================================================================

!     SOLVE THE BAROTROPIC EQUATION

!=======================================================================

      call tropic (c2dtsf, acor, cori, itt, dtts)

!-----------------------------------------------------------------------
!     if this is the 1st pass of an euler backward timestep, set the
!     disk pointers so the proper time levels are read on the 2nd pass
!     and go back to do the 2nd pass.
!-----------------------------------------------------------------------

      if (euler1) then
        eots      = .true.
        euler1    = .false.
        euler2    = .true.
        ntaux     = taum1disk
        taum1disk = taudisk
        taudisk   = taup1disk
        taup1disk = ntaux
        go to 1000
      endif
      if (wide_open_mw .and. euler2) then

!       shuffle "tau" and "tau+1" after euler backward to
!       insure data is in the right place for the next timestep

        call euler_shuffle

!       re-establish correct pointers for this timestep

        taum1 = mod(itt+0,3) - 1
        tau   = mod(itt+1,3) - 1
        taup1 = mod(itt+2,3) - 1
      endif

!-----------------------------------------------------------------------
!     output all remaining diagnostics
!-----------------------------------------------------------------------

      call diago

!-----------------------------------------------------------------------
!     save restart
!-----------------------------------------------------------------------

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_mom (0, fname)
          call mom_rest_out (fname, 1, imt, 1, jmt)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_mom (1, fname)
          call mom_rest_out (fname, 1, imt, 1, jmt)
        endif
      endif

!-----------------------------------------------------------------------
!     if it`s the last timestep then clean things up otherwise return
!-----------------------------------------------------------------------

      if (eorun) write (stdout,'(1x,a)') 'MOMdone'
#endif

      return
      end
