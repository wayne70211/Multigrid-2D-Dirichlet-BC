# 2D Poisson Multigrid Solver

### 2D Poisson Equation Problem
<br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;\nabla^2&space;u=f(x,y)=2x^2&plus;2y^2-4" title="\large \nabla^2 u=f(x,y)=2x^2+2y^2-4" />
</p>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;x,y\in&space;[-1,1]" title="\large x,y\in [-1,1]" /> 
</p>

with *Dirichlet Boundary Condition*  <br>
<br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;u(-1,y)=u(1,y)=u(x,-1)=u(x,1)=0" title="\large u(-1,y)=u(1,y)=u(x,-1)=u(x,1)=0" /><br>
</p>

The exact solution is <br>
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;u(x,y)=(x^2-1)(y^2-1)" title="\large u(x,y)=(x^2-1)(y^2-1)" /> <br>
</p>

The residual is <br>

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\large&space;r_{i,j}=f_{i,j}-\frac{u_{i&plus;1,j}-2u_{i,j}&plus;u_{i-1,j}}{\Delta&space;x^2}-\frac{u_{i,j&plus;1}-2u_{i,j}&plus;u_{i,j-1}}{\Delta&space;y^2}" title="\large r_{i,j}=f_{i,j}-\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{\Delta x^2}-\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{\Delta y^2}" /><br>
</p>

---
### Exact Solver on Coarsest Level
The coarsest grid should be solved to exact solution<br>
The program provides two solvers
1. **Gauss-Seidel** with residual L2-norm less than tolerance
2. **LU decomposition**

``` fortran
tol = 1.0d-6
exact_solver = 1 
```

###  Set Level

Set arbitrary level. The coarsest grid should greater than 3

``` fortran
Level_num = 6
if (Level_Nx(Level_num).le.3) then
    write(*,*) Level_num," level is high for ",Nx," grid "
    stop
end if
```

###  Set Multigrid Solver

**Gauss-Seidel** scheme is selected to be `smoother` <br>

**V Cycle Scheme**
``` fortran
subroutine MG_Vcycle(Nx,Ny,dx,dy,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Ny,Level_num,exact_solver
    real*8 ,intent(in) :: dx,dy,tol
    real*8,dimension(0:Nx,0:Ny), intent(in)    :: RHS
    real*8,dimension(0:Nx,0:Ny), intent(inout) :: U
    integer,dimension(Level_num) :: Level_Nx,Level_Ny
    real*8 ,dimension(Level_num) :: Level_dx,Level_dy
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,j,k,Level,iter
    integer, intent(inout) :: I_cycle
```
**Full Cycle Scheme**
``` fortran
subroutine MG_Fcycle(Nx,Ny,dx,dy,RHS,U,Level_num,tol,exact_solver,I_cycle)
    implicit none
    integer,intent(in) :: Nx,Ny,Level_num,exact_solver
    real*8 ,intent(in) :: dx,dy,tol
    real*8,dimension(0:Nx,0:Ny), intent(in)    :: RHS
    real*8,dimension(0:Nx,0:Ny), intent(inout) :: U
    integer,dimension(Level_num) :: Level_Nx,Level_Ny
    real*8 ,dimension(Level_num) :: Level_dx,Level_dy
    real*8  :: rms0,rms,rmsc
    integer :: Max_Iter,Relax_Iter
    integer :: i,j,k,Level,iter,inner_cycle
    integer, intent(out) :: I_cycle
```
---
## Multigrid V cycle Result
The following tables are the result of **Multigrid V cycle** with 5 and 6 levels <br>
The **P-** is parallel version implemented by `OpenMP` (Use 8 threads for test)  <br>


|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |  
| MG5       | 1024  |   6     |     0.8668  | 
| P-MG5     | 1024  |   6     |     0.2620  | 
| MG6       | 1024  |   6     |     0.3376  | 
| P-MG6     | 1024  |   6     |     0.1236  |

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |
| MG5       | 2048  |   6     |     9.7714  | 
| P-MG5     | 2048  |   6     |     1.9416  | 
| MG6       | 2048  |   6     |     1.7439  | 
| P-MG6     | 2048  |   6     |     0.6016  |
| MG7       | 2048  |   6     |     1.2188  | 
| P-MG7     | 2048  |   6     |     0.4571  |

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |
| MG5       | 4096  |   5     |   123.5220  | 
| P-MG5     | 4096  |   6     |    21.8231  | 
| MG6       | 4096  |   6     |    13.2918  | 
| P-MG6     | 4096  |   6     |     3.2305  | 
| MG7       | 4096  |   6     |     5.2741  | 
| P-MG7     | 4096  |   6     |     1.8986  | 

|   Solver  |  Np   |  Cycle  |   CPU Time  |  
| :---:     | :---: | :---:   |   :---:     |
| MG6       | 8192  |   5     |   135.0223  | 
| P-MG6     | 8192  |   6     |    26.9059  | 
| MG7       | 8192  |   6     |    27.4078  | 
| P-MG7     | 8192  |   6     |     8.2924  | 

## Compile and Run 
Use [`PGI Compiler`](https://www.pgroup.com/products/community.htm)

```shell
make
```

Run
```shell
./Run
```
## Reference
See [wikipedia](https://en.wikipedia.org/wiki/Multigrid_method)
