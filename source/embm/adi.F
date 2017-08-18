      subroutine adi (x, an, ans, as, ae, aew, aw, b, id, jd)

#if defined O_embm && defined O_embm_adi && !defined O_embm_mgrid
!=======================================================================
!                    UVic ADI Solver

!     implementation of an Alternating-Direction-Implicit solver
!=======================================================================

!=======================================================================
!     subroutine to solve
!   input:
!     an,ans,as: active coef. in north-south dir.
!     ae,aew,aw: active coef. in east-west dir.
!     b: accumulated fixed source term
!     id,jd:  array dimensions
!   output:
!     x: solution
!=======================================================================

      implicit none

      integer id, jd

      real an(jd,id), ans(jd,id), as(jd,id)
      real ae(id,jd), aew(id,jd), aw(id,jd)
      real x(id,jd), b(id,jd)

      integer i, j
      real y(jd,id), c(jd,id)

!     factor and solve in x (periodic)

      do j=1,jd
         call ctdma(x(1,j),aw(1,j),aew(1,j),ae(1,j),b(1,j),id)
      enddo

!     transfer x to c (transposing)

      do j=1,jd
         do i=1,id
            c(j,i) = x(i,j)
         enddo
      enddo

!     factor and solve in y (adiabatic)

      do i=1,id
         call tdma(y(1,i),as(1,i),ans(1,i),an(1,i),c(1,i),jd)
      enddo

!     transfer y to x (transposing)

      do j=1,jd
         do i=1,jd
            x(i,j) = y(j,i)
         enddo
      enddo

      return
      end

      subroutine ctdma (x, aw, ap, ae, b, id)
!=======================================================================
!   subroutine to do a cyclic tridiagonal matrix solve
!   input:
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     id:  array dimension
!   output:
!     x: solution
!=======================================================================

      implicit none

      integer id, i

      real x(id), aw(id), ap(id), ae(id), b(id)
      real factor, alpha(id), beta(id), theta(id)

      alpha(1) = 2*ap(1)
      do i=2,id-1
        alpha(i) = ap(i)
      enddo
      alpha(id) = ap(id) + ae(id)*aw(1)/ap(1)

      call tdma (x, aw, alpha, ae, b, id)

      beta(1) = -ap(1)
      do i=2,id-1
        beta(i) = 0.0
      enddo
      beta(id) = -ae(id)

      call tdma (theta, aw, alpha, ae, beta, id)

      factor = (x(1) + aw(1)/ap(1)*x(id))/
     &         (1. + theta(1) + aw(1)/ap(1)*theta(id))
      do i=1,id
        x(i) = x(i) - factor*theta(i)
      enddo

      return
      end

      subroutine tdma (x, aw, ap, ae, b, id)
!=======================================================================
!   subroutine to do a tridiagonal matrix solve
!   input:
!     ap,aw,ae: active coefficients for p,w,e nodes
!     b: accumulated fixed source term
!     id:  array dimension
!   output:
!     x: solution
!=======================================================================

      implicit none

      integer id, i

      real x(id), aw(id), ap(id), ae(id), b(id), alpha(id), beta

!     forward sweep

      beta = ap(1)
      x(1) = b(1)/beta
      do i=2,id
        alpha(i) = -ae(i-1)/beta
        beta = ap(i) + aw(i)*alpha(i)
        x(i) = (b(i) + aw(i)*x(i-1))/beta
      enddo

!     backward sweep

      do i=id-1,1,-1
        x(i) = x(i) - alpha(i+1)*x(i+1)
      enddo
#endif
      return
      end
