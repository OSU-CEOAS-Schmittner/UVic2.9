!====================== include file "cpolar.h" =========================

!     polar transform coefficients used to transform velocities near
!     poles before filtering

      real spsin, spcos
      common /cpolar_r/ spsin(imt), spcos(imt)
