!---------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers       |
! coated onto a spherical surface,                              | 
!---------------------------------------------------------------|
      
program brushweakpolyelectrolyte 
  
  !     .. variable and constant declaractions 

  use globals       ! parameters definitions 
  use physconst
  use mathconst
  use volume
  use random
  use field
  use parameters
  use matrices
  use energy
  use chains
  use VdW

  implicit none  
 
  real*8,  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
  real*8,  dimension(:), allocatable :: xguess    ! guess fraction  solvent 
  real*8,  dimension(:), allocatable :: fvec     
  
  integer :: i,c             ! dummy indices    
     
  real*8  :: pKdCa(20)             
     
  !     .. executable statements 
  !     .. init 

  call read_inputfile()
  call init_constants()
  call init_matrices()      ! init matrices for chain generation
  call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)
  call make_chains(chainmethod) ! generate polymer configurations 
  call make_sequence_chain(period,chaintype)
  call allocate_geometry(nsize)
  call make_geometry()      ! generate volume elements lattice 
  call allocate_field(nsize) 
  call read_VdWCoeff()      

  

  allocate(x(neq)) ! 
  allocate(xguess(4*nsize))
  allocate(fvec(4*nsize))

  !     .. computation starts

  !     .. pKaCad binding strength
  pKdCa(1)=1.0d0
  pKdCa(2)=1.5d0
  pKdCa(3)=2.0d0
  pKdCa(4)=2.5d0
  pKdCa(5)=3.0d0
  pKdCa(6)=3.5d0
  pKdCa(7)=4.0d0
  pKdCa(8)=4.5d0
  pKdCa(9)=5.0d0
  pKdCa(10)=5.5d0
  pKdCa(11)=6.0d0
  pKdCa(12)=6.5d0
  pKdCa(13)=7.0d0
  pKdCa(14)=7.5d0
  pKdCa(15)=8.0d0
  pKdCa(16)=8.5d0
  pKdCa(17)=9.0d0
  pKdCa(18)=9.5d0
  pKdCa(19)=10.0d0
  
  neq = 5*nr                  ! size of lattice        
  iter = 0                  ! iteration counter 
  count= 1                  ! file counter 
  
  do c=1,19                 ! loop over pKa Ca 
     
       pKa(4)=pKdCa(c)
   
       call init_expmu()
       call make_guess(x,xguess,pKa(4),pKdCa(1))

       ! .. solver  SCFMT eqs 
!       call fcnelect(x,fvec,n)
       call solver(x, xguess, n, error, fnorm) 
       
       
       do i=1,nr           
          xsol(i)=x(i)        ! solvent density=volume fraction   
          psi(i) =x(i+nr)     ! potential           
       enddo
       
       call fcnenergy()       ! free energy 
       !         call endpoint(count)   ! compute endpoint pdf
       call average_height()      
       call charge_polymer()
       call average_charge_polymer()
       call output(count)     ! writing of output
       
       count = count+1        ! next  
       iter  = 0              ! reset of iteration counter 
          
    enddo
    
    stop
    
  end program brushweakpolyelectrolyte
  
