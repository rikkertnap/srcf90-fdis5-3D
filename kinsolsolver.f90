!------------------------------------------------------------------------------| 
!     kinsolver.f                                                              |  
!     module containing functions and subroutines for the usage of             |  
!     the numerical kinsol library routine                                     |
!------------------------------------------------------------------------------| 

module kinsolvars

  implicit none
  !     .. variables
  
  real*8, dimension(:), allocatable ::   pp ! pre-processor variable

end module kinsolvars


!     The routine fkpset is the preconditioner setup routine. It must have
!     that specific name to be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

  use globals, only : neq
  use kinsolvars
  implicit none
  
  integer :: ier
  integer*8 :: i
  double precision :: udata(*), uscale(*), fdata(*), fscale(*)
  double precision :: vtemp1(*), vtemp2(*)
  
  do i = 1, neq
     pp(i) = 0.5d0 / (udata(i) + 5.0d0)
  enddo
  ier = 0
  
end subroutine fkpset

!     The routine fkpsol is the preconditioner solve routine. It must have
!     that specific name to be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
  
  use globals, only : neq
  use kinsolvars
  implicit none
  
  integer :: ier
  integer*8 :: i
  double precision :: udata(*), uscale(*), fdata(*), fscale(*)
  double precision :: vv(*), ftem(*)
  
  do  i = 1, neq
     vv(i) = vv(i) * pp(i)
  enddo
  ier = 0
  
  return
end subroutine fkpsol



!     Initialiazaition and run of the kinsolver                                 
!     Scale Preconditioned GMRES solver                                        
!     pre:  input vector x and guess of x (xguess) and n (number of equations) 
!            and fnorm                                                         
!     post: solution of SCMFT stored in x and  fnorm residual error            


subroutine kinsol_gmres_solver(x, xguess, n, error, fnorm)
  
  use globals, only : nsize,neq
  use kinsolvars

  implicit none
  
  !     .. variables and constant declaractions 
  
  !     .. array arguments     
 
  real*8 :: x(neq)
  real*8 :: xguess(neq)

  !      real*8, dimension(:) :: x    
  !      real*8, dimension(:) :: xguess
  
  !     .. scalar arguments
  real*8 :: fnorm 
  real*8 :: error
  integer :: n                 !  number of equations
  
  !     .. local arguments 
  
  integer*4 :: ier             ! Kinsol error flag
!  integer*8 :: neq             ! Kinsol number of equations
  
  integer*8 :: iout(15)        ! Kinsol additional output information
  real*8 :: rout(2)            ! Kinsol additional out information
  
  integer :: i                 ! dummy index 
  integer*8 :: msbpre
  integer*8 :: maxniter
  real*8 :: fnormtol, scsteptol
  real*8 :: scale(neq)
  real*8 :: constr(neq)
  
  integer*4 :: globalstrat, maxl, maxlrst
  
  
  !     .. executable statements 
  
  !     .. init of kinsol variables 

  ! neq = n                                
  msbpre  = 5               ! maximum number of iterations without prec. setup 
  fnormtol = error          ! Function-norm stopping tolerance
  scsteptol = error         ! Function-norm stopping tolerance
  maxl = neq                ! maximum Krylov subspace dimension 
  maxlrst = 20               ! maximum number of restarts
  globalstrat = 0           ! inexact Newton  
  maxniter =1000            ! maximum of nonlinear iterations default 200  


  allocate(pp(neq))

  
  do i = 1, neq             
     constr(i) = 0.0        ! constraint vector  
     scale(i) = 1.0         ! scaling vector  
     x(i) = xguess(i)       ! initial guess
  enddo
  
  
  call fnvinits(3, neq, ier) ! inits NVECTOR module
    
  if (ier .ne. 0) then      ! 3 for Kinsol, neq number of equations, ier error flag (must be 0 on output)
     print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
     stop
  endif
  
  call fkinmalloc(iout, rout, ier) ! Allocates memory and output additional information
  
  if (ier .ne. 0) then
     print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
     stop
  endif
  
  call fkinsetiin('MAX_SETUPS', msbpre, ier) ! Additional input information
  call fkinsetiin('MAX_NITER', maxniter, ier) 
  call fkinsetrin('FNORM_TOL', fnormtol, ier)
  call fkinsetrin('SSTEP_TOL', scsteptol, ier)
  call fkinsetvin('CONSTR_VEC', constr, ier) 
  
  !     .. init Scale Preconditioned GMRES solver 
  
  call fkinspgmr(maxl, maxlrst, ier)   
  if (ier .ne. 0) then
     print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
     call fkinfree          ! free memory
     stop
  endif
  
  !     .. preconditioner   
  call fkinspilssetprec(0, ier) 

  !     .. call solver
  
  call fkinsol(x, globalstrat, scale, scale, ier) 
  fnorm=rout(1) 
  
  if (ier .lt. 0) then
     print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
     print*, 'Linear Solver returned IER = ', iout(9)
     call fkinfree
     stop
  endif
  
  print*,'Found solution: fnorm = ',fnorm
  
  call fkinfree             ! free memory
  deallocate(pp)
  return
  
end subroutine kinsol_gmres_solver


!     .. wrapper function 

subroutine fkfun(x,f,ier)

  use globals,  only  : neq
  use fcnpointer

  implicit none

  real*8, dimension(neq) :: x  ! expliciet size array                                                                                                                                        
  real*8, dimension(neq) :: f

  integer*4 ier

  call fcnptr(x,f,neq)

end subroutine fkfun


















