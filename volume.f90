!     makes volume elements 
!     for spherical coordinates 
module volume      
  implicit none
  
  
  !     .. variables
   
  real*8  :: delta              ! delta  spacing of lattice site in z-direction
  integer :: nz                 ! nz number of lattice sites in z-direction  nz <= nsize
  integer :: nzmax              ! nzmax  maximum number of lattice sites in z-direction
  integer :: nzmin              ! nzmin minimumal number of lattice sites in z-direction  
  integer :: nzstep             ! nzstep number of lattice sites stepped over or reduced 
  real*8, dimension(:), allocatable :: zc ! z-coordinate
  real*8, dimension(:), allocatable :: G ! geometrical factor 
  real*8, dimension(:), allocatable :: deltaG ! geometrical factor
  real*8, dimension(:), allocatable :: Fplus ! factor in Poisson Eq  
  real*8, dimension(:), allocatable :: Fmin
  
contains
  
  subroutine allocate_geometry(N)
    implicit none
    integer, intent(in) :: N
    allocate(zc(N))
    allocate(G(N))
    allocate(deltaG(N+1))
!    allocate(Fplus(N))
!    allocate(Fmin(N))
    
  end subroutine allocate_geometry
  
  
  
  subroutine  make_geometry()
    
    use globals
    
    implicit none
    
    integer ::  i
    real*8  ::  vol
    real*8  ::  vtest
    real*8, parameter :: vtol=1.0d-5

  
    vol=0.0d0
    
    do i=1,nz
       zc(i)= (i-0.5d0) * delta  ! z-coordinate 
       G(i) =  1.0d0  ! geometrical factor 
       deltaG(i)=1.0d0  
!G(i) +(delta*delta)/(12.0d0*radius*radius) ! delta G(i)= (1/delta) \int dr G(r) 
!       Fplus(i)=1.0d0+ delta/rc(i)
!       Fmin(i) = 2.0d0-Fplus(i)      ! factors in Poisson Equation
       vol=vol+ deltaG(i)
    enddo
    
    vtest=nz
    vtest=vtest
    if(dabs(vtest-vol)>=vtol) then 
       print*,"Warning vtest=",Vtest," not equal to vol=",vol
    end if
    
    return
      
    end subroutine make_geometry
    
  end module volume
  
