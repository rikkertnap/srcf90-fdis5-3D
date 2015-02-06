module field
  
  !     .. variables

!  implicit none
    
  real*8, dimension(:), allocatable :: xpolAB  ! total volume fraction of polymer on sphere
  real*8, dimension(:), allocatable :: xpolC   ! total volume fraction of polymer on sphere: hydrocarbon chain
  real*8, dimension(:), allocatable :: rhopolA ! density A monomer of polymer on sphere
  real*8, dimension(:), allocatable :: rhopolB ! density B monomer of polymer on sphere
  real*8, dimension(:), allocatable :: rhopolC ! density C monomer of polymer on sphere
  real*8, dimension(:), allocatable :: xsol    ! volume fraction solvent
  real*8, dimension(:), allocatable :: psi     ! electrostatic potential 
  real*8, dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
  real*8, dimension(:), allocatable :: xK      ! volume fraction of positive K+ ion
  real*8, dimension(:), allocatable :: xCa     ! volume fraction of positive Ca2+ ion
  real*8, dimension(:), allocatable :: xNaCl   ! volume fraction of NaCl ion pair
  real*8, dimension(:), allocatable :: xKCl    ! volume fraction of KCl  ion pair
  real*8, dimension(:), allocatable :: xCl     ! volume fraction of negative ion
  real*8, dimension(:), allocatable :: xHplus  ! volume fraction of Hplus
  real*8, dimension(:), allocatable :: xOHmin  ! volume fraction of OHmin 
  real*8, dimension(:), allocatable :: rhoq    ! total charge density in units of vsol
  real*8, dimension(:), allocatable :: qpol    ! charge density of polymer
  real*8, dimension(:,:), allocatable :: fdisA   ! degree of dissociation 
  real*8, dimension(:,:), allocatable :: fdisB   ! degree of dissociation
  
  real*8 :: qAB              ! normalization partion fnc polymer 
  real*8 :: qC              ! normalization partion fnc polymer 
  
contains

  subroutine allocate_field(N)
    implicit none
    
    integer, intent(in) :: N
    
    allocate(xpolAB(N))
    allocate(xpolC(N))
    allocate(rhopolA(N))
    allocate(rhopolB(N))
    allocate(rhopolC(N))
    allocate(xsol(N))
    allocate(psi(N+1))
    allocate(xNa(N))
    allocate(xK(N))
    allocate(xCa(N))
    allocate(xNaCl(N)) 
    allocate(xKCl(N)) 
    allocate(xCl(N)) 
    allocate(xHplus(N))
    allocate(xOHmin(N))
    allocate(rhoq(N))
    allocate(qpol(N))
    allocate(fdisA(5,N))
    allocate(fdisB(5,N))
    
    
  end subroutine allocate_field
  

  subroutine deallocate_field()
    implicit none
    
    
    deallocate(xpolAB)
    deallocate(xpolC)
    deallocate(rhopolA)
    deallocate(rhopolB)
    deallocate(rhopolC)
    deallocate(xsol)
    deallocate(psi)
    deallocate(xNa)
    deallocate(xK)
    deallocate(xCa)
    deallocate(xNaCl) 
    deallocate(xKCl) 
    deallocate(xCl) 
    deallocate(xHplus)
    deallocate(xOHmin)
    deallocate(rhoq)
    deallocate(qpol)
    deallocate(fdisA)
    deallocate(fdisB)
    
    
  end subroutine deallocate_field
  

  
  !     .. compute average height of tethered layer 
  !     .. first moment of density profile 
  
  subroutine average_height()
    
    use globals
    use parameters
    use volume
    
    implicit none
    
    integer :: i
    
    real*8 :: zero,first
    
    first=0.0d0               ! first moment 
    zero=0.0d0                ! zero moment  
    do i=1,nz
       zero=zero+xpolAB(i)*deltaG(i)
       first=first+xpolAB(i)*zc(i)*deltaG(i)
    enddo
    
    heightAB=first/zero

    first=0.0d0               ! first moment 
    zero=0.0d0                ! zero moment  
    do i=1,nz
        zero=zero+xpolC(i)*deltaG(i)
        first=first+xpolC(i)*zc(i)*deltaG(i)
    enddo

    heightC=first/zero


  end subroutine average_height
  
  subroutine charge_polymer()
    
    use globals
    use volume
    use parameters

    implicit none

    integer :: i
    real*8  :: npolA,npolB
    
    
    qpolA=0.0d0
    qpolB=0.0d0
    
    
    do i=1,nz
       qpolA=qpolA+(zpolA(1)*fdisA(1,i)*rhopolA(i)+&  
            zpolA(4)*fdisA(4,i)*rhopolA(i))*deltaG(i) 
       qpolB=qpolB+(zpolB(1)*fdisB(1,i)*rhopolB(i)+&
            zpolB(4)*fdisB(4,i)*rhopolB(i))*deltaG(i)
    enddo
    
    qpolA=qpolA*delta
    qpolB=qpolB*delta
    qpol_tot=qpolA+qpolB
    
  end subroutine charge_polymer

  ! .. post : return average charge of state of 
  !   of polymers

  subroutine average_charge_polymer()
        
    use globals
    use volume
    use parameters
    use chains

    implicit none 
   
    integer :: i,s,k
    real*8  :: npolA,npolB
    
    
    ! .. number of A and B monomors 
    npolA=0
    do s=1,nsegAB
       if(isAmonomer(s).eqv..true.) then
          npolA=npolA+1
       endif
    enddo
    npolB=nsegAB-npolA
    
    if(npolA.ne.0) then
       do k=1,5
          avfdisA(k)=0.0d0
          do i=1,nz
             avfdisA(k)=avfdisA(k)+fdisA(k,i)*rhopolA(i)*deltaG(i) 
          enddo
          avfdisA(k)=avfdisA(k)*delta/(sigmaAB*delta*npolA)
       enddo
    else
       do k=1,5
          avfdisA(k)=0.0d0
       enddo
    endif
    
    if(npolB.ne.0) then
       do k=1,5
          avfdisB(k)=0.0d0
          do i=1,nz
             avfdisB(k)=avfdisB(k)+fdisB(k,i)*rhopolB(i)*deltaG(i) 
          enddo
          avfdisB(k)=avfdisB(k)*delta/(sigmaAB*delta*npolB)
       enddo
    else
       do k=1,5
          avfdisB(k)=0.0d0
       enddo
    endif
    
  end subroutine average_charge_polymer
  

  logical function isNaN(x)
    implicit none
    real*8 :: x
    if (x /= x) then
        isNaN=.true.
    else
      isNaN=.false.
    endif 

  end function isNaN 
  
end module field

