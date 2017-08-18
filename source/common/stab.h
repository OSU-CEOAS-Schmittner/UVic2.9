!====================== include file "stab.h" ==========================

!     CFL and other stability criteria information

!     cflons  = starting longitude (degrees) for stability tests
!     cflone  = ending longitude (degrees) for stability tests
!     cflats  = starting latitude (degrees) for stability tests
!     cflate  = ending latitude (degrees) for stability tests
!     cfldps  = starting depth (cm) for stability tests
!     cfldpe  = ending depth (cm) for stability tests

!     iscfl   = index corresponding to "cflons"
!     iecfl   = index corresponding to "cflone"
!     jscfl   = index corresponding to "cflats"
!     jecfl   = index corresponding to "cflate"
!     kscfl   = index corresponding to "cfldps"
!     kecfl   = index corresponding to "cfldpe"

!     cflcrt  = factor by which the cfl criteria must be exceeded in
!               order to show local values  (see blkdta.F)
!     maxcfl  = maximum number of times the "cflcrt" factor can be
!               exceeded before stopping.
!     numcfl  = counter for number of times the "cflcrt" factor was
!               exceeded.
!     cflum   = zonal velocity which comes closest to its cfl criteria
!     cflup   = percent of cfl criteria reached by "cflum"
!     icflu   = "i" coordinate of "cflum"
!     jcflu   = "j" coordinate of "cflum"
!     kcflu   = "k" coordinate of "cflum"

!     cflvm   = meridional velocity which comes closest to its cfl
!               criteria
!     cflvp   = percent of cfl criteria reached by "cflvm"
!     icflv   = "i" coordinate of "cflvm"
!     jcflv   = "j" coordinate of "cflvm"
!     kcflv   = "k" coordinate of "cflvm"

!     cflwtm  = vertical velocity at "t" box bottom
!                which comes closest to its cfl criteria
!     cflwtp  = percent of cfl criteria reached by "cflwtm"
!     icflwt  = "i" coordinate of "cflwtm"
!     jcflwt  = "j" coordinate of "cflwtm"
!     kcflwt  = "k" coordinate of "cflwtm"

!     cflwum  = vertical velocity at "u,v" box bottom
!                which comes closest to its cfl criteria
!     cflwup  = percent of cfl criteria reached by "cflwum"
!     icflwu  = "i" coordinate of "cflwum"
!     jcflwu  = "j" coordinate of "cflwum"
!     kcflwu  = "k" coordinate of "cflwum"

!     reynx   = maximum reynolds number in the zonal direction
!     ireynx  = "i" coordinate of "reynx"
!     jreynx  = "j" coordinate of "reynx"
!     kreynx  = "k" coordinate of "reynx"
!     reynu   = "u" for computing "reynx"
!     reynmu  = zonal mixing of momentum for computing "reynx"

!     reyny   = maximum reynolds number in the meridional direction
!     ireyny  = "i" coordinate of "reyny"
!     jreyny  = "j" coordinate of "reyny"
!     kreyny  = "k" coordinate of "reyny"
!     reynv   = "v" for computing "reyny"
!     reynmv  = meridional mixing of momentum for computing "reyny"

!     reynz   = maximum reynolds number in the vertical direction
!     ireynz  = "i" coordinate of "reynz"
!     jreynz  = "j" coordinate of "reynz"
!     kreynz  = "k" coordinate of "reynz"
!     reynw   = "w" for computing "reynz"
!     reynmw  = vertical mixing of momentum for computing "reynz"

!     peclx   = maximum peclet number in the zonal direction
!     ipeclx  = "i" coordinate of "peclx"
!     jpeclx  = "j" coordinate of "peclx"
!     kpeclx  = "k" coordinate of "peclx"
!     peclu   = "u" for computing "peclx"
!     peclmu  = zonal mixing of tracer for computing "peclx"

!     pecly   = maximum peclet number in the meridional direction
!     ipecly  = "i" coordinate of "pecly"
!     jpecly  = "j" coordinate of "pecly"
!     kpecly  = "k" coordinate of "pecly"
!     peclv   = "v" for computing "pecly"
!     peclmv  = meridional mixing of tracer for computing "pecly"

!     peclz   = maximum peclet number in the vertical direction
!     ipeclz  = "i" coordinate of "peclz"
!     jpeclz  = "j" coordinate of "peclz"
!     kpeclz  = "k" coordinate of "peclz"
!     peclw   = "w" for computing "peclz"
!     peclmw  = vertical mixing of tracer for computing "peclz"

!     tdig    = factor by which local tracer extremum must be exceeded
!               before showing ficticious creation of tracer

      integer iscfl,  iecfl,  jscfl,  jecfl,  kscfl,  kecfl
      integer icflu, jcflu, kcflu, icflv, jcflv, kcflv
      integer icflwt, jcflwt, kcflwt, icflwu, jcflwu, kcflwu
      integer numcfl, maxcfl
      integer ireynx, jreynx, kreynx, ireyny, jreyny, kreyny
      integer ireynz, jreynz, kreynz, ipeclx, jpeclx, kpeclx
      integer ipecly, jpecly, kpecly, ipeclz, jpeclz, kpeclz

      common /stab_i/ iscfl,  iecfl,  jscfl,  jecfl,  kscfl,  kecfl
      common /stab_i/ icflu, jcflu, kcflu, icflv, jcflv, kcflv
      common /stab_i/ icflwt, jcflwt, kcflwt, icflwu, jcflwu, kcflwu
      common /stab_i/ numcfl, maxcfl
      common /stab_i/ ireynx, jreynx, kreynx, ireyny, jreyny, kreyny
      common /stab_i/ ireynz, jreynz, kreynz, ipeclx, jpeclx, kpeclx
      common /stab_i/ ipecly, jpecly, kpecly, ipeclz, jpeclz, kpeclz

      real cflons, cflone, cflats, cflate, cfldps, cfldpe
      real cflup, cflum, cflvp, cflvm, cflwtp, cflwtm
      real cflwup, cflwum, cflcrt, tdig
      real reynx, reynu, reynmu, reyny, reynv, reynmv
      real reynz, reynw, reynmw, peclx, peclu, peclmu
      real pecly, peclv, peclmv, peclz, peclw, peclmw

      common /stab_r/ cflons, cflone, cflats, cflate, cfldps, cfldpe
      common /stab_r/ cflup, cflum, cflvp, cflvm, cflwtp, cflwtm
      common /stab_r/ cflwup, cflwum, cflcrt, tdig
      common /stab_r/ reynx, reynu, reynmu, reyny, reynv, reynmv
      common /stab_r/ reynz, reynw, reynmw, peclx, peclu, peclmu
      common /stab_r/ pecly, peclv, peclmv, peclz, peclw, peclmw
