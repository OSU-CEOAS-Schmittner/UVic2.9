!====================== include file "matrix.h" ==========================
!       variables for extracting the matrix

      integer numtiles
      parameter (numtiles = km)
      character (10) :: tilename

      real tile(imt,km,jmw,numtiles)
      real Aexp(imt,km,jmw,numtiles)
      real Aimp(imt,km,jmw,numtiles)
      real dtmat, dtmatrix(km)
      common /MATRIX_FIELDS/ tile, Aexp, Aimp, dtmat, dtmatrix

      integer iTile(numtiles), ntilestart, ntileend
      integer stepCount
      integer expMatrixCounter, impMatrixCounter
      integer expMatrixWriteCount, impMatrixWriteCount
      real expMatrixWriteTime, impMatrixWriteTime
      common /MATRIX_PARAMS_I/iTile, ntilestart, ntileend, 
     &   stepCount, 
     &   expMatrixCounter, impMatrixCounter,
     &   expMatrixWriteCount, impMatrixWriteCount
      common /MATRIX_PARAMS_R/
     &   expMatrixWriteTime, impMatrixWriteTime
      