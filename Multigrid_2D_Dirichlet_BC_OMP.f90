module Multigrid_2D_Dirichlet_BC_OMP
  use omp_lib
  use Exact_Solver
  implicit none
contains
  !---------------------------------------------------------------------------!
  ! Multigrid V cycle scheme (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_MG_Vcycle(Nx,Ny,dx,dy,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Ny,Level_num,exact_solver
    real*8 ,intent(in) :: dx,dy,tol
    real*8,dimension(0:Nx,0:Ny), intent(in)     :: RHS
    real*8,dimension(0:Nx,0:Ny), intent(inout)  :: U
    integer,dimension(Level_num) :: Level_Nx,Level_Ny
    real*8 ,dimension(Level_num) :: Level_dx,Level_dy
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,j,k,Level,iter
    integer, intent(inout) :: I_cycle

    !Defined all data space
    type space
      real*8,allocatable :: data(:,:)
    end type space

    type (space) UL(Level_num),R(Level_num),F(Level_num),P(Level_num)

    allocate(UL(1)%data(0:Nx,0:Ny),R(1)%data(0:Nx,0:Ny),F(1)%data(0:Nx,0:Ny),P(1)%data(0:Nx,0:Ny))

    Level_Nx(1) = Nx
    Level_Ny(1) = Ny
    Level_dx(1) = dx
    Level_dy(1) = dy

    !$omp parallel do private(i,k)
    do i=2,Level_num
      k=2**(i-1)
      Level_Nx(i) = Nx / k
      Level_Ny(i) = Ny / k
      Level_dx(i) = dx * dble(k)
      Level_dy(i) = dy * dble(k)
      allocate(UL(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( R(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( F(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( P(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
    enddo
    !$omp end parallel do


    Max_Iter   = 100000     ! Allowed maximum number of outer iteration
    Relax_Iter = 2          ! Number of relaxation for restriction in V-cycle

    if (I_cycle .eq. -1) then
      Max_Iter = 1
    endif

    !Check the coarsest grid
    if (Level_Nx(Level_num).le.3) then
    write(*,*) Level_num," level is high for this grid.."
    stop
    end if

    !$omp parallel do private(i,j)
    do j=0,Level_Ny(1)
    do i=0,Level_Nx(1)
      F(1)%data(i,j)  = RHS(i,j)
      UL(1)%data(i,j) = U(i,j)
      R(1)%data(i,j)  = 0.0d0
    enddo
    enddo
    !$omp end parallel do

    !Compute initial resitual:
    call P_Residual(Nx,Ny,dx,dy,F(1)%data,UL(1)%data,R(1)%data)
    call P_L2norm(Nx,Ny,R(1)%data,rms0)

    !open(50,file='residual_All_PMG5V2.plt')
    !write(50,*) 'variables ="k","iter","rms","rms/rms0"'

    iter=0

    do I_cycle=1,Max_Iter

    !---------------------------------------------------------------------------!
      do Level=1,Level_num-1

        ! Relax
        do i=1,Relax_Iter
          !iter=iter+1
          call P_Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
          !call P_Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)
          !call P_L2norm(Level_Nx(Level),Level_Ny(Level),R(Level)%data,rms)
          !write(50,*) I_cycle,iter,rms,rms/rms0
        end do

        ! Compute residual
        call P_Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)

        ! Check for convergence on finest grid
        call P_L2norm(Level_Nx(Level),Level_Ny(Level),R(Level)%data,rms)
        if (rms/rms0.le.tol .and. Level .eq. 1) goto 10

        ! Restriction
        call P_Restriction(Level_Nx(Level),Level_Ny(Level),Level_Nx(Level+1),Level_Ny(Level+1),R(Level)%data,F(Level+1)%data)

        !$omp parallel do private(i,j)
        do j=0,Level_Ny(Level+1)
        do i=0,Level_Nx(Level+1)
          UL(Level+1)%data(i,j) = 0.0d0
        end do
        end do
        !$omp end parallel do


      end do
    !---------------------------------------------------------------------------!

      ! Compute residual on coarsest grid
      call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
      call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rmsc)

      ! Solve exact solution on coarsest grid
      if (exact_solver .eq. 1) then
        do while (rms/rmsc .gt. tol)
          iter=iter+1
          call P_Relax(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
          call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
          ! Check for convergence on smallest grid
          call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
          !write(50,*) I_cycle,iter,rms,rms/rms0
        end do
      else if (exact_solver .eq. 2) then
        call LU_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
        !call Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
        !call L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
        !write(50,*) I_cycle,iter,rms,rms/rms0
      else
        call CG_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)

      endif

    !---------------------------------------------------------------------------!

      do Level=Level_num-1,1,-1

        ! Prolongation
        call P_Prolongation(Level_Nx(Level+1),Level_Ny(Level+1),Level_Nx(Level),Level_Ny(Level),UL(Level+1)%data,P(Level)%data)

        ! Correct
        !$omp parallel do private(i,j)
        do j=1,Level_Ny(Level)-1
        do i=1,Level_Nx(Level)-1
          UL(Level)%data(i,j) = UL(Level)%data(i,j) + P(Level)%data(i,j)
        end do
        end do
        !$omp end parallel do

        ! Relax
        do i=1,Relax_Iter
          !iter=iter+1
          call P_Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
          !call P_Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)
          !call P_L2norm(Level_Nx(Level),Level_Ny(Level),R(Level)%data,rms)
          !write(50,*) I_cycle,iter,rms,rms/rms0
        end do

      end do


    end do  ! Outer iteration loop

    10 continue

    !$omp parallel do private(i,j)
    do j=0,Ny
    do i=0,Nx
      U(i,j) = UL(1)%data(i,j)
    end do
    end do
    !$omp end parallel do

    !close(50)

    do i=1,Level_num
      deallocate(UL(i)%data,R(i)%data,F(i)%data,P(i)%data)
    enddo


    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Multigrid Full cycle scheme (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_MG_Fcycle(Nx,Ny,dx,dy,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Ny,Level_num,exact_solver
    real*8 ,intent(in) :: dx,dy,tol
    real*8,dimension(0:Nx,0:Ny), intent(in)     :: RHS
    real*8,dimension(0:Nx,0:Ny), intent(inout)  :: U
    integer,dimension(Level_num) :: Level_Nx,Level_Ny
    real*8 ,dimension(Level_num) :: Level_dx,Level_dy
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,j,k,Level,Inner_Level,iter,inner_cycle
    integer, intent(out) :: I_cycle

    !Defined all data space
    type space
      real*8,allocatable :: data(:,:)
    end type space

    type (space) UL(Level_num),R(Level_num),F(Level_num),P(Level_num)

    allocate(UL(1)%data(0:Nx,0:Ny),R(1)%data(0:Nx,0:Ny),F(1)%data(0:Nx,0:Ny),P(1)%data(0:Nx,0:Ny))

    Level_Nx(1) = Nx
    Level_Ny(1) = Ny
    Level_dx(1) = dx
    Level_dy(1) = dy

    !$omp parallel do private(i,k)
    do i=2,Level_num
      k=2**(i-1)
      Level_Nx(i) = Nx / k
      Level_Ny(i) = Ny / k
      Level_dx(i) = dx * dble(k)
      Level_dy(i) = dy * dble(k)
      allocate(UL(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( R(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( F(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
      allocate( P(i)%data(0:Level_Nx(i),0:Level_Ny(i)))
    enddo
    !$omp end parallel do


    Max_Iter   = 100000     ! Allowed maximum number of outer iteration
    Relax_Iter = 2          ! Number of relaxation for restriction in V-cycle

    !Check the coarsest grid
    if (Level_Nx(Level_num).le.3) then
      write(*,*) Level_num," level is high for ",Nx," grid "
    stop
    end if

    !$omp parallel do private(i,j)
    do j=0,Level_Ny(1)
    do i=0,Level_Nx(1)
      F(1)%data(i,j)  = RHS(i,j)
      UL(1)%data(i,j) = U(i,j)
      R(1)%data(i,j)  = 0.0d0
    enddo
    enddo
    !$omp end parallel do

    !Compute initial resitual:
    call P_Residual(Nx,Ny,dx,dy,F(1)%data,UL(1)%data,R(1)%data)
    call P_L2norm(Nx,Ny,R(1)%data,rms0)

    !open(70,file='residual_All_PMG5F2.plt')
    !write(70,*) 'variables ="k","iter","rms","rms/rms0"'

    iter=0

    do I_cycle=1,Max_Iter

    !---------------------------------------------------------------------------!
      do Level=1,Level_num-1

        ! Compute residual
        call P_Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)

        ! Check for convergence on finest grid
        call P_L2norm(Level_Nx(Level),Level_Ny(Level),R(Level)%data,rms)
        if (rms/rms0.le.tol .and. Level .eq. 1) goto 10

        ! Restriction
        call P_Restriction(Level_Nx(Level),Level_Ny(Level),Level_Nx(Level+1),Level_Ny(Level+1),R(Level)%data,F(Level+1)%data)

        !$omp parallel do private(i,j)
        do j=0,Level_Ny(Level+1)
        do i=0,Level_Nx(Level+1)
          UL(Level+1)%data(i,j) = 0.0d0
        end do
        end do
        !$omp end parallel do


      end do
    !---------------------------------------------------------------------------!

      ! Compute residual on coarsest grid
      call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
      call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rmsc)

      ! Solve exact solution on coarsest grid
      if (exact_solver .eq. 1) then
        do while (rms/rmsc .gt. tol)
          iter=iter+1
          call P_Relax(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
          call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
          ! Check for convergence on smallest grid
          call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
          !write(70,*) I_cycle,iter,rms,rms/rms0
        end do
      else if (exact_solver .eq. 2) then
        call LU_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
        !call Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
        !call L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
        !write(70,*) I_cycle,iter,rms,rms/rms0
      else
        call CG_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)

      endif


    !---------------------------------------------------------------------------!

      do Level=Level_num-1,1,-1

        ! Prolongation
        call P_Prolongation(Level_Nx(Level+1),Level_Ny(Level+1),Level_Nx(Level),Level_Ny(Level),UL(Level+1)%data,P(Level)%data)

        ! Correct
        !$omp parallel do private(i,j)
        do j=1,Level_Ny(Level)-1
        do i=1,Level_Nx(Level)-1
          UL(Level)%data(i,j) = UL(Level)%data(i,j) + P(Level)%data(i,j)
        end do
        end do
        !$omp end parallel do


        ! Start Inner V Cycle
        !---------------------------------------------------------------------------!
        do Inner_Level=Level,Level_num-1

          ! Relax
          do i=1,Relax_Iter
            call P_Relax(Level_Nx(Inner_Level),Level_Ny(Inner_Level),Level_dx(Inner_Level),Level_dy(Inner_Level),F(Inner_Level)%data,UL(Inner_Level)%data)
          end do

          ! Compute residual
          call P_Residual(Level_Nx(Inner_Level),Level_Ny(Inner_Level),Level_dx(Inner_Level),Level_dy(Inner_Level),F(Inner_Level)%data,UL(Inner_Level)%data,R(Inner_Level)%data)

          ! Restriction
          call P_Restriction(Level_Nx(Inner_Level),Level_Ny(Inner_Level),Level_Nx(Inner_Level+1),Level_Ny(Inner_Level+1),R(Inner_Level)%data,F(Inner_Level+1)%data)

          !$omp parallel do private(i,j)
          do j=0,Level_Ny(Inner_Level+1)
          do i=0,Level_Nx(Inner_Level+1)
            UL(Inner_Level+1)%data(i,j) = 0.0d0
          end do
          end do
          !$omp end parallel do


        end do
        !---------------------------------------------------------------------------!

        ! Compute residual on coarsest grid
        call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
        call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rmsc)

        ! Solve exact solution on coarsest grid
        if (exact_solver .eq. 1) then
          do while (rms/rmsc .gt. tol)
            call P_Relax(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
            call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
            ! Check for convergence on smallest grid
            call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
          end do
        else if (exact_solver .eq. 2) then
          call LU_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
        else
          call CG_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)

        endif

        !---------------------------------------------------------------------------!
        do Inner_Level=Level_num-1,Level,-1

          ! Prolongation
          call P_Prolongation(Level_Nx(Inner_Level+1),Level_Ny(Inner_Level+1),Level_Nx(Inner_Level),Level_Ny(Inner_Level),UL(Inner_Level+1)%data,P(Inner_Level)%data)

          ! Correct
          !$omp parallel do private(i,j)
          do j=1,Level_Ny(Inner_Level)-1
          do i=1,Level_Nx(Inner_Level)-1
            UL(Inner_Level)%data(i,j) = UL(Inner_Level)%data(i,j) + P(Inner_Level)%data(i,j)
          end do
          end do
          !$omp end parallel do

          ! Relax
          do i=1,Relax_Iter
            call P_Relax(Level_Nx(Inner_Level),Level_Ny(Inner_Level),Level_dx(Inner_Level),Level_dy(Inner_Level),F(Inner_Level)%data,UL(Inner_Level)%data)
          end do

        end do
        ! End Inner V Cycle
        !---------------------------------------------------------------------------!

      end do


      ! Start V Cycle
      !---------------------------------------------------------------------------!

      do Level=1,Level_num-1

        ! Relax
        do i=1,Relax_Iter
          call P_Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
        end do

        ! Compute residual
        call P_Residual(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data,R(Level)%data)

        ! Restriction
        call P_Restriction(Level_Nx(Level),Level_Ny(Level),Level_Nx(Level+1),Level_Ny(Level+1),R(Level)%data,F(Level+1)%data)

        !$omp parallel do private(i,j)
        do j=0,Level_Ny(Level+1)
        do i=0,Level_Nx(Level+1)
          UL(Level+1)%data(i,j) = 0.0d0
        end do
        end do
        !$omp end parallel do


      end do
    !---------------------------------------------------------------------------!

      ! Compute residual on coarsest grid
      call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
      call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rmsc)

      ! Solve exact solution on coarsest grid
      if (exact_solver .eq. 1) then
        do while (rms/rmsc .gt. tol)
          call P_Relax(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
          call P_Residual(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data,R(Level_num)%data)
          ! Check for convergence on smallest grid
          call P_L2norm(Level_Nx(Level_num),Level_Ny(Level_num),R(Level_num)%data,rms)
        end do
      else if (exact_solver .eq. 2) then
        call LU_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)
      else
        call CG_Solver(Level_Nx(Level_num),Level_Ny(Level_num),Level_dx(Level_num),Level_dy(Level_num),F(Level_num)%data,UL(Level_num)%data)

      endif

    !---------------------------------------------------------------------------!

      do Level=Level_num-1,1,-1

        ! Prolongation
        call P_Prolongation(Level_Nx(Level+1),Level_Ny(Level+1),Level_Nx(Level),Level_Ny(Level),UL(Level+1)%data,P(Level)%data)

        ! Correct
        !$omp parallel do private(i,j)
        do j=1,Level_Ny(Level)-1
        do i=1,Level_Nx(Level)-1
          UL(Level)%data(i,j) = UL(Level)%data(i,j) + P(Level)%data(i,j)
        end do
        end do
        !$omp end parallel do


        ! Relax
        do i=1,Relax_Iter
          call P_Relax(Level_Nx(Level),Level_Ny(Level),Level_dx(Level),Level_dy(Level),F(Level)%data,UL(Level)%data)
        end do

      end do

      ! End V Cycle
      !---------------------------------------------------------------------------!


    end do  ! Outer iteration loop

    10 continue

    !$omp parallel do private(i,j)
    do j=0,Ny
    do i=0,Nx
      U(i,j) = UL(1)%data(i,j)
    end do
    end do
    !$omp end parallel do

    !close(70)

    do i=1,Level_num
      deallocate(UL(i)%data,R(i)%data,F(i)%data,P(i)%data)
    enddo


    return
  end subroutine


  !---------------------------------------------------------------------------!
  ! Relaxation formula for Poisson equation (Parallel Version)
  ! Uses GS relaxation
  !---------------------------------------------------------------------------!
  subroutine P_Relax(Nx,Ny,dx,dy,F,U)
    implicit none
    integer ,intent(in) :: Nx,Ny
    real*8  ,intent(in) :: dx,dy
    real*8  ,dimension(0:Nx,0:Ny),intent(in)    :: F
    real*8  ,dimension(0:Nx,0:Ny),intent(inout) :: U

    real*8  :: a
    integer :: i,j

    a = -2.0d0/(dx*dx) - 2.0d0/(dy*dy)

    !$omp parallel do private(i,j)
    do j=1,Ny-1
    do i=1,Nx-1
      U(i,j) = (1.0d0/a)*(F(i,j) &
              &  - (U(i+1,j)+U(i-1,j))/(dx*dx) &
              &  - (U(i,j+1)+U(i,j-1))/(dy*dy) )
    end do
    end do
    !$omp end parallel do

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Residual formula for Poisson equation (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_Residual(Nx,Ny,dx,dy,F,U,R)
    implicit none
    integer,intent(in) :: Nx,Ny
    real*8 ,intent(in) :: dx,dy
    real*8, dimension(0:Nx,0:Ny),intent(in)  :: U,F
    real*8, dimension(0:Nx,0:Ny),intent(out) :: R

    integer :: i,j

    !$omp parallel
    !$omp do private(i,j)
    do j=1,Ny-1
    do i=1,Nx-1
      R(i,j) = F(i,j) - (U(i+1,j) - 2.0d0*U(i,j) + U(i-1,j))/(dx*dx) &
                   &  - (U(i,j+1) - 2.0d0*U(i,j) + U(i,j-1))/(dy*dy)
    end do
    end do
    !$omp end do

    !Boundary conditions for residuals
    !$omp do private(i)
    do i=0,Nx
      R(i,0)  = 0.0d0
      R(i,Ny) = 0.0d0
    end do
    !$omp end do

    !$omp do private(j)
    do j=0,Ny
      R(0,j)  = 0.0d0
      R(Nx,j) = 0.0d0
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Compute L2-norm for an array (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_L2norm(Nx,Ny,R,RMS)
    implicit none
    integer ,intent(in) :: Nx,Ny
    real*8, dimension(0:Nx,0:Ny),intent(in) :: R
    real*8  ,intent(out):: RMS
    integer :: i,j

    RMS=0.0d0

    !$omp parallel
    !$omp do private(i,j) reduction(+:RMS)
    do j=1,Ny-1
    do i=1,Nx-1
      RMS = RMS + R(i,j)*R(i,j)
    end do
    end do
    !$omp end do
    !$omp end parallel

    RMS= dsqrt(RMS/dfloat((Nx-1)*(Ny-1)))

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Restriction operators (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_Restriction(Nxf,Nyf,Nxh,Nyh,R,F)
    implicit none
    integer ,intent(in) :: Nxf,Nyf,Nxh,Nyh
    real*8, dimension(0:Nxf,0:Nyf),intent(in)  :: R       !on higher grid
    real*8, dimension(0:Nxh,0:Nyh),intent(out) :: F       !on lower grid
    integer :: i,j

    !$omp parallel
    !$omp do private(i,j)
    do j=1,Nyh-1
    do i=1,Nxh-1
      F(i,j) = 1.0d0/16.0d0 * (4.0d0*R(2*i,2*j) &
            &      + 2.0d0 * (R(2*i+1,2*j)+R(2*i-1,2*j)+R(2*i,2*j+1)+R(2*i,2*j-1))        &
            &      + 1.0d0 * (R(2*i+1,2*j+1)+R(2*i-1,2*j-1)+R(2*i-1,2*j+1)+R(2*i+1,2*j-1)))
    end do
    end do
    !$omp end do

    !Boundary conditions
    !$omp do private(i)
    do i=0,Nxh
      F(i,0)   = R(2*i,0)
      F(i,Nyh) = R(2*i,Nyf)
    end do
    !$omp end do

    !$omp do private(j)
    do j=0,Nyh
      F(0,j)   = R(0,2*j)
      F(Nxh,j) = R(Nxf,2*j)
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine

  !---------------------------------------------------------------------------!
  ! Prolongation operator (Parallel Version)
  !---------------------------------------------------------------------------!
  subroutine P_Prolongation(Nxh,Nyh,Nxf,Nyf,U,P)
    implicit none
    integer ,intent(in) :: Nxf,Nyf,Nxh,Nyh
    real*8, dimension(0:Nxh,0:Nyh) ,intent(in)  :: U       !on lower grid
    real*8, dimension(0:Nxf,0:Nyf) ,intent(out) :: P       !on higher grid
    integer :: i,j

    !$omp parallel
    !$omp do private(i,j)
    do j=0,Nyh-1
    do i=0,Nxh-1
      P(2*i,2*j)     = U(i,j)
      P(2*i+1,2*j)   = 1.0d0/2.0d0*(U(i,j)+U(i+1,j))
      P(2*i,2*j+1)   = 1.0d0/2.0d0*(U(i,j)+U(i,j+1))
      P(2*i+1,2*j+1) = 1.0d0/4.0d0*(U(i,j)+U(i,j+1)+U(i+1,j)+U(i+1,j+1))
    end do
    end do
    !$omp end do

    !Boundary conditions
    !$omp do private(j)
    do j=0,Nyh
      P(Nxf,2*j)     = U(Nxh,j)
    end do
    !$omp end do

    !$omp do private(i)
    do i=0,Nxh
      P(2*i,Nyf)     = U(i,Nyh)
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine

end module
