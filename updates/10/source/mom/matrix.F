      SUBROUTINE MATRIX_INIT( )

      IMPLICIT NONE

# if defined O_matrix
      include "size.h"
      include "pconst.h"
      include "matrix.h"
      
      integer ntile, n, k
      character(120) :: fn, fname, new_file_name
      logical exists

!     set time step for computation of implicit matrix
      do k=1,km
        dtmatrix(k) = dtmat
      enddo

C     Initialize counters
      stepCount = 0
      expMatrixWriteCount = 1
      impMatrixWriteCount = 1
	  expMatrixCounter=0
	  impMatrixCounter=0      

      Aexp(:,:,:,:) = c0
      Aimp(:,:,:,:) = c0
            
!     read tracer tiles/initial conditions
      do ntile=1,numtiles
		write(fn,'(A6,I2.2,A4)') 'O_tile',ntile,'.bin'
		fname = new_file_name (fn)
		write(*,*)'TMM: reading tile from ',fname
		inquire (file=trim(fname), exist=exists)
		if (exists) then			  
		  call read_tracer_array(fname,tile(1,1,1,ntile),1)
		else
          print*, "Warning => Can not find ",trim(fname)
		endif				
	  enddo  

# endif

      RETURN
      END

      SUBROUTINE MATRIX_STORE_EXPLICIT(
     I                  ntile,tracer,twodt )

      IMPLICIT NONE

      include "size.h"
      include "matrix.h"
      
      integer ntile
      real tracer(imt,km,jmw), twodt(km)

      character(120) :: fname, new_file_name
      integer k

!     compute and accumulate explicit matrix
	  do k=1,km
		Aexp(:,k,:,ntile) = Aexp(:,k,:,ntile) + 
     &       ((tracer(:,k,:)-tile(:,k,:,ntile))/twodt(k))
      enddo
	  if (ntile.eq.numtiles) then
	    stepCount = stepCount + 1
	    expMatrixCounter = expMatrixCounter + 1
	    write(*,*)'DEBUG3: accumulating explicit matrix ',
     &   expMatrixCounter, stepCount, twodt(1),twodt(2)
      endif
	  
!     reset tracer to tiles/initial condition
	  tracer(:,:,:)=tile(:,:,:,ntile)		  

      RETURN
      END

      SUBROUTINE MATRIX_STORE_IMPLICIT(
     I                  ntile,tracer )

      IMPLICIT NONE

      include "size.h"
      include "matrix.h"
      
      integer ntile
      real tracer(imt,km,jmw)

      character(120) :: fname, new_file_name

!     compute and accumulate implicit matrix
	  Aimp(:,:,:,ntile) = Aimp(:,:,:,ntile) + tracer(:,:,:)
	  
	  if (ntile.eq.numtiles) then
	    impMatrixCounter = impMatrixCounter + 1
	    write(*,*)'DEBUG4: accumulating implicit matrix ',
     &   impMatrixCounter
      endif
      
      RETURN
      END

      SUBROUTINE MATRIX_WRITE(writeExpMatrix, writeImpMatrix)

      IMPLICIT NONE

      include "size.h"
      include "pconst.h"      
      include "matrix.h"

      logical writeExpMatrix, writeImpMatrix
      
      integer ntile
      character(120) :: fn, fname, new_file_name
      logical exists
      real recipExpMatrixCounter, recipImpMatrixCounter

C      writeExpMatrix = .TRUE.
C      writeImpMatrix = .TRUE.

      IF (writeExpMatrix) THEN
		recipExpMatrixCounter = c1/expMatrixCounter
		write(*,*)'DEBUG1 ',recipExpMatrixCounter,expMatrixCounter
        Aexp(:,:,:,:) = Aexp(:,:,:,:)*recipExpMatrixCounter
        do ntile=1,numtiles
		  write(*,*) 'Writing explicit matrix :',ntile,
     &        expMatrixWriteCount, expMatrixCounter        
		  write(fn,'(A9,I2.2,A4)') 'MATRIXEXP',ntile,'.bin'
		  fname = new_file_name (fn)
		  call write_tracer_array(fname,Aexp(1,1,1,ntile),
     &	                          expMatrixWriteCount)
        enddo
		expMatrixCounter=0        
		expMatrixWriteCount = expMatrixWriteCount + 1

C       Reset explicit matrix
        Aexp(:,:,:,:) = c0        
      ENDIF

      IF (writeImpMatrix) THEN
		recipImpMatrixCounter = c1/impMatrixCounter
		write(*,*)'DEBUG2 ',recipImpMatrixCounter,impMatrixCounter
        Aimp(:,:,:,:) = Aimp(:,:,:,:)*recipImpMatrixCounter
        do ntile=1,numtiles
		  write(*,*) 'Writing implicit matrix :',ntile,
     &        impMatrixWriteCount, impMatrixCounter                
		  write(fn,'(A9,I2.2,A4)') 'MATRIXIMP',ntile,'.bin'
		  fname = new_file_name (fn)
		  call write_tracer_array(fname,Aimp(1,1,1,ntile),
     &	                          impMatrixWriteCount)
        enddo
		impMatrixCounter=0        
		impMatrixWriteCount = impMatrixWriteCount + 1

C       Reset implicit matrix
        Aimp(:,:,:,:) = c0        
      ENDIF


      RETURN
      END

C======================================================
C This subroutine reads a model tracer array from file
C======================================================
      SUBROUTINE READ_TRACER_ARRAY(filename,fld,iSlice)

      implicit none

      include "size.h"
      include "param.h"

      integer iSlice      
      real*8 fld(1:imt,1:km,1:jmt)
      character*(*) filename

      real*8 tmp(imtm2,jmtm2,km)       
      integer j, k

      call read_r8_field(imtm2*jmtm2,km,iSlice,tmp,filename)
      
      do j=2,jmtm1
        do k=1,km
          fld(2:imtm1,k,j)=tmp(1:imtm2,j-1,k)
        enddo
      enddo
      fld(2:imtm1,1:km,1)=fld(2:imtm1,1:km,2)
      fld(2:imtm1,1:km,jmt)=fld(2:imtm1,1:km,jmtm1)
	  do j=1,jmt
		call setbcx (fld(1:imt,1:km,j),imt,km)
	  enddo  
      
      
      END

C======================================================
C This subroutine writes a model tracer array from file
C======================================================
      SUBROUTINE WRITE_TRACER_ARRAY(filename,fld,iSlice)

      implicit none

      include "size.h"
      include "param.h"
      
      integer iSlice      
      real*8 fld(1:imt,1:km,1:jmt)
      character*(*) filename

      real*8 tmp(imtm2,jmtm2,km)       
      integer j, k

      do j=2,jmtm1
        do k=1,km
          tmp(1:imtm2,j-1,k)=fld(2:imtm1,k,j)
        enddo
      enddo  

      call write_r8_field(imtm2*jmtm2,km,iSlice,tmp,filename)
            
      END
