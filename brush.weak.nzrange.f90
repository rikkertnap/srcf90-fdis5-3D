! ---------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers        |
! coated onto a planar surface,                                  | 
! input/output: see myio.f90                                     | 
! ---------------------------------------------------------------|
      
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
    use surface
    
    implicit none  
    
    real*8,  dimension(:), allocatable :: x         ! iteration vector 
    real*8,  dimension(:), allocatable :: xguess    ! guess iteration vector
    real*8,  dimension(:), allocatable :: xstored   ! stored iteration vector
    real*8,  dimension(:), allocatable :: fvec     
  
    integer :: i,c,eps          ! dummy indices
    integer, parameter :: cmax=20
    integer, parameter :: epsmax=19
    real*8  :: cs(cmax)         ! salt concentrations     
    real*8  :: VdWeps(epsmax)   ! VdW eps values
    integer :: countfile        ! file counter        
    logical :: use_xstored       
    logical :: isfirstguess   
   
    ! .. executable statements 
    ! .. init 

    call read_inputfile()
    call init_constants()
    call init_matrices()        ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)
    call make_chains(chainmethod) ! generate polymer configurations 
    call make_sequence_chain(period,chaintype)
    call allocate_geometry(nsize)
    call make_geometry()        ! generate volume elements lattice 
    call allocate_field(nsize) 
    !    call read_VdWCoeff() ! THIS NEED TO BE CHANGED !!!!!!!!!!!!!!!!!!!
    call set_size_neq()
    call init_expmu()
    call init_surface(bcflag)

    !  .. computation starts
    
    nz=nzmax                    ! first distance 
    call set_size_neq()         ! number of non-linear equation neq
    neqmax = neq
    allocate(xstored(neq))
    isfirstguess = .true.    
    use_xstored = .false.
    countfile = 1                
    iter = 0

    do while (nz.ge.nzmin)    
 
        call set_size_neq() 
        allocate(x(neq))
        allocate(xguess(neq))
        call chain_filter()    
        call make_guess(x, xguess, isfirstguess, use_xstored, xstored)
        call solver(x, xguess, error, fnorm) 
        call fcnenergy()         ! free energy
        call average_height()      
        call charge_polymer()
        call average_charge_polymer()
        call output(countfile)    ! writing of output

        isfirstguess =.false.    
        use_xstored = .true.
        countfile = countfile+1  ! next  
        iter = 0                 ! reset of iteration counter 
        nz = nz-nzstep             ! reduce distance 
        do i=1,neq
            xstored(i)=x(i)
        enddo

        deallocate(x)   
        deallocate(xguess)
    enddo   

    deallocate(xstored)

    call deallocate_field()

end program brushweakpolyelectrolyte
