module Exact_Solver
  implicit none
contains
  !---------------------------------------------------------------------------!
  ! LU Solver for Poisson equation
  ! The RHS and output Phi are two dimensional matrix (Which is reshaped to calculate )
  ! Reshape(RHS) => B
  ! Solve AX = B
  ! Reshape(X)   => Phi
  !---------------------------------------------------------------------------!

  subroutine LU_Solver(Nx,Ny,dx,dy,RHS,Phi)
    implicit none
    integer,intent(in) :: Nx,Ny
    real*8, intent(in) :: dx,dy
    real*8,dimension(0:Nx,0:Ny),intent(in)        :: RHS
    real*8,dimension(0:Nx,0:Ny),intent(inout)     :: Phi
    real*8,dimension((Nx-1)*(Ny-1),(Nx-1)*(Ny-1)) :: A
    real*8,dimension((Nx-1)*(Ny-1),(Nx-1)*(Ny-1)) :: L,U
    real*8,dimension((Nx-1)*(Ny-1))               :: B,X,bp
    integer :: i,j,k,Np
    real*8  :: sum

    !---------------------------------------------------------------------------!
    ! Initialize
    Np=(Nx-1)*(Ny-1)

    do j=0,Ny
    do i=0,Nx
      Phi(i,j)   = 0.0d0
    enddo
    enddo

    do j=1,Np
    do i=1,Np
      A(i,j)   = 0.0d0
      L(i,j)   = 0.0d0
      U(i,j)   = 0.0d0
    enddo
    enddo

    do i=1,Np
      bp(i) = 0.0d0
      B(i)  = 0.0d0
      X(i)  = 0.0d0
    enddo

    ! Initialize A Matrix
    A(1,1:2)=(/-4,1/)
    A(Np,Np-1:Np)=(/1,-4/)

    do i=2,Np-1
      A(i,i-1:i+1)=(/1,-4,1/)
    end do

    do i=1,Np-Nx+1
      A(i+Nx-1,i)=1
      A(i,i+Nx-1)=1
    enddo

    ! Correct Boundary
    do i=Nx-1,Np-1,Nx-1
      A(i,i+1)=0
    enddo
    do i=Nx,Np-1,Nx-1
      A(i,i-1)=0
    enddo


    !---------------------------------------------------------------------------!
    ! Reshape RHS Matrix
    k=0
    do j=1,Ny-1
    do i=1,Nx-1
      k = k + 1
      B(k) = RHS(i,j)*(dx*dx)
      if (i .eq. 1) then
        B(k) = B(k) - Phi(0,j)
      endif
      if (i .eq. Nx-1) then
        B(k) = B(k) - Phi(Nx,j)
      endif
      if (j .eq. 1) then
        B(k) = B(k) - Phi(i,0)
      endif
      if (j .eq. Ny-1) then
        B(k) = B(k) - Phi(i,Ny)
      endif
    enddo
    enddo

    !---------------------------------------------------------------------------!
    ! LU decomposition
    do i=1,Np
      L(i,1)=A(i,1)
      U(1,i)=A(1,i)/L(1,1)
    end do

    do j=2,Np
      do i=j,Np
        sum=dble(0)
        do k=1,j-1
          sum=sum+L(i,k)*U(k,j)
        end do
        L(i,j)=A(i,j)-sum
      end do
      U(j,j)=dble(1)
      do i=j+1,Np
        sum=dble(0)
        do k=1,j-1
          sum=sum+L(j,k)*U(k,i)
        end do
        U(j,i)=(A(j,i)-sum)/L(j,j)
      end do
    end do

    ! Backward substitution
    bp(1)=B(1)/L(1,1)
    do i=2,Np
      sum=dble(0)
      do k=1,i-1
        sum=sum+L(i,k)*bp(k)
      enddo
      bp(i)=(b(i)-sum)/L(i,i)
    end do

    X=bp
    do j=Np-1,1,-1
      sum=dble(0)
      do k=j+1,Np
        sum=sum+U(j,k)*X(k)
      end do
      X(j)=bp(j)-sum
    end do

    !---------------------------------------------------------------------------!
    ! Reshape RHS Matrix
    k=0
    do j=1,Ny-1
    do i=1,Nx-1
      k = k + 1
      Phi(i,j) = X(k)
    enddo
    enddo


  end subroutine

  !---------------------------------------------------------------------------!
  ! CG Solver for Poisson equation
  ! The RHS and output Phi are two dimensional matrix (Which is reshaped to calculate )
  ! Reshape(RHS) => B
  ! Solve AX = B
  ! Reshape(X)   => Phi
  !---------------------------------------------------------------------------!

  subroutine CG_Solver(Nx,Ny,dx,dy,RHS,Phi)
    implicit none
    integer,intent(in) :: Nx,Ny
    real*8, intent(in) :: dx,dy
    real*8,dimension(0:Nx,0:Ny),intent(in)        :: RHS
    real*8,dimension(0:Nx,0:Ny),intent(inout)     :: Phi
    real*8,dimension((Nx-1)*(Ny-1),(Nx-1)*(Ny-1)) :: A
    real*8,dimension((Nx-1)*(Ny-1))               :: B,X,bp
    real*8,dimension((Nx-1)*(Ny-1))               :: r,p,temp_p
    integer :: i,j,k,Np,niter,iter
    real*8  :: Allow_iter_error,sum
    real*8  :: ak,ak_denominator,ak_numerator,bk,bk_numerator,r_max


    !---------------------------------------------------------------------------!
    ! Initialize
    Allow_iter_error = 1d-12
    niter = 200000
    Np=(Nx-1)*(Ny-1)

    do j=0,Ny
    do i=0,Nx
      Phi(i,j)   = 0.0d0
    enddo
    enddo

    do j=1,Np
    do i=1,Np
      A(i,j)   = 0.0d0
    enddo
    enddo

    do i=1,Np
      bp(i) = 0.0d0
      B(i)  = 0.0d0
      X(i)  = 0.0d0
    enddo

    ! Initialize A Matrix
    A(1,1:2)=(/-4,1/)
    A(Np,Np-1:Np)=(/1,-4/)

    do i=2,Np-1
      A(i,i-1:i+1)=(/1,-4,1/)
    end do

    do i=1,Np-Nx+1
      A(i+Nx-1,i)=1
      A(i,i+Nx-1)=1
    enddo

    ! Correct Boundary
    do i=Nx-1,Np-1,Nx-1
      A(i,i+1)=0
    enddo
    do i=Nx,Np-1,Nx-1
      A(i,i-1)=0
    enddo


    !---------------------------------------------------------------------------!
    ! Reshape RHS Matrix
    k=0
    do j=1,Ny-1
    do i=1,Nx-1
      k = k + 1
      B(k) = RHS(i,j)*(dx*dx)
      if (i .eq. 1) then
        B(k) = B(k) - Phi(0,j)
      endif
      if (i .eq. Nx-1) then
        B(k) = B(k) - Phi(Nx,j)
      endif
      if (j .eq. 1) then
        B(k) = B(k) - Phi(i,0)
      endif
      if (j .eq. Ny-1) then
        B(k) = B(k) - Phi(i,Ny)
      endif
    enddo
    enddo


    !---------------------------------------------------------------------------!
    ! CG Solver
    iter=0
    r=B
    p=r

    do k=1,niter
      temp_p=dble(0)
      ak_numerator=dble(0)
      ak_denominator=dble(0)
      bk_numerator=dble(0)
      r_max=dble(0)

      do i=1,Np
        do j=1,Np
          temp_p(i)=p(j)*A(i,j)+temp_p(i)
        end do
        ak_numerator=r(i)*r(i)+ak_numerator
      enddo

      do i=1,Np
        ak_denominator=temp_p(i)*p(i)+ak_denominator
      end do
      ak=ak_numerator/ak_denominator

      do i=1,Np
        X(i)=X(i)+ak*p(i)
        r(i)=r(i)-ak*temp_p(i)
        r_max=r(i)*r(i)+r_max
      end do

      if ( dsqrt(r_max) .lt. Allow_iter_error) then
        iter=k
        exit
      endif

      do i=1,Np
        bk_numerator=r(i)*r(i)+bk_numerator
      end do

      bk=bk_numerator/ak_numerator

      do i=1,Np
        p(i)=r(i)+bk*p(i)
      end do

      iter=k

    enddo

    !---------------------------------------------------------------------------!
    ! Reshape RHS Matrix
    k=0
    do j=1,Ny-1
    do i=1,Nx-1
      k = k + 1
      Phi(i,j) = X(k)
    enddo
    enddo

  end subroutine
end module
