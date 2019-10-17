!-----------------------------------------------------------------------------!
!  Multigrid method for solving 2D Poisson equation
!     d2u/dx2 + d2u/dy2 = f(x,y)
!     Drichlet b.c.
!-----------------------------------------------------------------------------!
program Poisson_2D
  use omp_lib
  use Multigrid_2D_Dirichlet_BC
  use Multigrid_2D_Dirichlet_BC_OMP
  implicit none
  integer::i,j,k,m,nx,ny,exact_solver,Level_num,I_cycle
  character(len=7) :: mode
  character(len=1) :: Level,version
  real*8,dimension(:,:),allocatable ::u,f,ue,e
  real*8,dimension(:),allocatable ::x,y
  real*8 ::dx,dy,tol,rms,x0,xL,y0,yL,start,finish


  !Domain
  x0 =-1.0d0 !left
  xL = 1.0d0 !right

  y0 =-1.0d0 !bottom
  yL = 1.0d0 !up

  !Tolerance
  tol= 1.0d-6

  !---------------------------------------------------------------------------!
  ! Exact Solver
  !     1. Gauss Seidel Solver
  !     2. LU decomposition
  !     3. CG Solver 
  !---------------------------------------------------------------------------!
  exact_solver = 1  !
  !---------------------------------------------------------------------------!

  open(66,file='Result.txt')
  write(66,*)"-------------------------------------------------------"
  write(66,*)"Solver     Np     Iter        CPU Time        L2-norm"
  write(66,*)"-------------------------------------------------------"

  write(*,*)"-------------------------------------------------------"
  write(*,*)"Solver     Np     Iter        CPU Time        L2-norm"
  write(*,*)"-------------------------------------------------------"

  do k=8,11

    ! Multigrid Level
    do Level_num=6,7

    !number of points
    nx = 2**k !number of grid points in x
    ny = nx   !number of grid points in y

    !grid spacing (spatial)
    dx = (xL-x0)/dfloat(nx)
    dy = (yL-y0)/dfloat(ny)

    !spatial coordinates
    allocate(x(0:nx))
    do i=0,nx
    x(i) = x0 + dfloat(i)*dx
    end do

    allocate(y(0:ny))
    do j=0,ny
    y(j) = y0 + dfloat(j)*dy
    end do

    allocate(u(0:nx,0:ny))
    allocate(f(0:nx,0:ny))
    allocate(e(0:nx,0:ny))
    allocate(ue(0:nx,0:ny))

    !---------------------------------------------!
    ! Exact solution
    !---------------------------------------------!
    do j=0,ny
    do i=0,nx
    f(i,j) =-2.0d0*(2.0d0-x(i)*x(i)-y(j)*y(j))
    ue(i,j)= (x(i)*x(i)-1.0d0)*(y(j)*y(j)-1.0d0)
    end do
    end do


    !Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    !Boundary conditions has to satisfy exact solution
    do i=0,nx
    u(i,0)  = ue(i,0)
    u(i,ny) = ue(i,ny)
    end do

    do j=0,ny
    u(0,j)  = ue(0,j)
    u(nx,j) = ue(nx,j)
    end do

    !---------------------------------------------------------------------------!
    ! Solver:
    !---------------------------------------------------------------------------!

    !Level_num = 6
    write(Level,'(I1.1)') Level_num
    write(version,'(I1.1)') exact_solver
    mode = 'MG'//Level//'-V'//version

    ! Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    start=omp_get_wtime()
    call MG_Vcycle(Nx,Ny,dx,dy,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do j=0,nx
    do i=0,ny
    e(i,j) = dabs(u(i,j)-ue(i,j))
    end do
    end do

    !L-2 Norm:
    call L2norm(nx,ny,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    mode = 'PMG'//Level//'-V'//version

    ! Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    start=omp_get_wtime()
    call P_MG_Vcycle(Nx,Ny,dx,dy,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do j=0,nx
    do i=0,ny
    e(i,j) = dabs(u(i,j)-ue(i,j))
    end do
    end do

    call L2norm(nx,ny,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    mode = 'MG'//Level//'-F'//version
    !Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    start=omp_get_wtime()
    call MG_Fcycle(Nx,Ny,dx,dy,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do j=0,nx
    do i=0,ny
    e(i,j) = dabs(u(i,j)-ue(i,j))
    end do
    end do

    call L2norm(nx,ny,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    mode = 'PMG'//Level//'-F'//version
    !Numerical solution:
    do j=0,nx
    do i=0,ny
    u(i,j)=0.0d0
    end do
    end do

    start=omp_get_wtime()
    call P_MG_Fcycle(Nx,Ny,dx,dy,F,U,Level_num,tol,exact_solver,I_cycle)
    finish=omp_get_wtime()

    !----------------------!
    ! Error analysis:
    !----------------------!
    do j=0,nx
    do i=0,ny
    e(i,j) = dabs(u(i,j)-ue(i,j))
    end do
    end do

    call L2norm(nx,ny,e,rms)

    write(*,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms
    write(66,'(A7,i8,i8,f18.12,f15.10)') mode,nx,I_cycle,finish-start,rms

    !---------------------------------------------------------------------------!
    write(*,*)"-------------------------------------------------------"
    write(66,*)"-------------------------------------------------------"

    deallocate(u,f,e,ue,x,y)
    enddo

  enddo

  close(66)


end program
