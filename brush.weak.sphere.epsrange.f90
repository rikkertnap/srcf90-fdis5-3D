!---------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers       |
! coated onto a spherical surface,                              | 
! input: see myio.f90                                           | 
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
    use listfcn
    use initxvector
    
    implicit none  
    
    real*8,  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
    real*8,  dimension(:), allocatable :: xguess    ! guess fraction  solvent 
    real*8,  dimension(:), allocatable :: fvec     
    
    integer :: i,c              ! dummy indices    
    integer :: countfile        ! file counter        
    real*8  :: eps(20)          ! VdW strength 
    
    !     .. executable statements 
    !     .. init 
    
    call read_inputfile()
    call init_constants()
    call init_matrices()        ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)
    call make_chains(chainmethod) ! generate polymer configurations 
    call make_sequence_chain(period,chaintype)
    call allocate_geometry(nsize)
    call make_geometry()        ! generate volume elements lattice 
    call allocate_field(nsize) 
    call read_VdWCoeff()
    print*,"main: sysflag=",sysflag
    call set_size_neq()         ! number of non-linear equation neq
    
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(fvec(neq))
    
    !     .. computation starts
    
    !     .. salt concentrations
    
    eps(1)=0.00d0
    eps(2)=0.25d0
    eps(3)=0.50d0
    eps(4)=0.75d0
    eps(5)=1.00d0
    eps(6)=1.25d0
    eps(7)=1.50d0
    eps(8)=1.75d0
    eps(9)=2.00d0
    eps(10)=2.5d0
    eps(11)=3.0d0
    eps(16)=3.5d0
    eps(17)=4.0d0
    eps(18)=4.5d0
    eps(19)=5.0d0
    
    countfile= 1    
 
    do c=1,17               ! loop over salt concentration
     
        VdWepsB = eps(c)
        call init_expmu()
        
        call make_guess(x,xguess,eps(c),eps(1))
!        call fcnelect(x,fvec,neq)
        
        ! solver  SCFMT eqs 
        
        call solver(x, xguess, error, fnorm) 
        
        do i=1,nr           
            xsol(i)=x(i)        ! solvent density=volume fraction   
            psi(i) =x(i+nr)     ! potential           
        enddo
       
        call fcnenergy()       ! free energy 
        !         call endpoint(count)   ! compute endpoint pdf
        call average_height()      
        call charge_polymer()
        call average_charge_polymer()
        call output(countfile)     ! writing of output
        
        countfile = countfile+1        ! next  
        iter  = 0              ! reset of iteration counter 
        
    enddo
     
    stop
     
end program brushweakpolyelectrolyte
  
