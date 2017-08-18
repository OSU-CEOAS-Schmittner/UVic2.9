!====================== include file "cshort.h" ========================

!     pen    = double exponential penetration function defined at
!              bottoms of "T" boxes

!              note: pen(0) is set 0.0 instead of 1.0 to compensate for
!                    the shortwave part of the total surface flux "stf"

!     divpen = divergence of penetration defined at the center of
!              "T" boxes

      real pen, divpen
      common /cshort_r/ pen(0:km), divpen(km)
