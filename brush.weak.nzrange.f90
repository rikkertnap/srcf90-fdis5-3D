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
    use myio
    use myutils
    use chaingenerator
    
    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector

    integer :: i             ! dummy indices       
    logical :: use_xstored       
    logical :: isfirstguess   
    character(len=lenText) :: text

    ! .. executable statements 
    ! .. init 

    LogName='status.log'
    call open_logfile(LogUnit,LogName) 
    text='program begins'
    call print_to_log(LogUnit,text)

    call read_inputfile()
    call init_constants()
    call init_matrices()            ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)  
    call make_sequence_chain(period,chaintype)
    call set_properties_chain(period,chaintype)  
    call make_chains(chainmethod)   ! generate polymer configurations 
    call allocate_geometry(nsize)
    call make_geometry()            ! generate volume elements lattice 
    call allocate_field(nsize) 

    !    call read_VdWCoeff() ! THIS NEED TO BE CHANGED !!!
    
    call set_size_neq()             ! number of non-linear equation neq    
    call init_expmu()
    call init_surface(bcflag)

    !  .. computation starts
    
    nz=nzmax                    
    neqmax = neq    
    allocate(xstored(neq))
    isfirstguess = .true.    
    use_xstored = .false.             
    iter = 0

    do while (nz>=nzmin)        ! loop distances

        call set_size_neq()  
        
        allocate(x(neq))
        allocate(xguess(neq))

        call chain_filter()    
        call make_guess(x, xguess, isfirstguess, use_xstored, xstored)
        call solver(x, xguess, error, fnorm)
        call fcnenergy()        
        call average_height()      
        call charge_polymer()
        call average_charge_polymer()
        call output()           ! writing of output

        isfirstguess =.false.    
        use_xstored = .true.
        iter = 0                ! reset of iteration counter 
        nz = nz-nzstep          ! reduce distance 
        do i=1,neq
            xstored(i)=x(i)
        enddo

        deallocate(x)   
        deallocate(xguess)
    enddo   

    deallocate(xstored)
    call deallocate_field()

    text="program end"
    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program brushweakpolyelectrolyte
