! ---------------------------------------------------------------|
! Solves the SCMFT eqs for WEAK polyelectrolytes polymers        |
! coated onto a planar surface,                                  | 
! input/output: see myio.f90                                     | 
! ---------------------------------------------------------------|
      
program brushweakpolyelectrolyte 
    
    !     .. variable and constant declaractions

    use mpivars 
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
    use chaingenerator
!    use VdW
    use listfcn 
    use fcnpointer
    use initxvector
    use surface
    use myio
    use myutils

    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec  

    integer :: i             ! dummy indices       
    logical :: use_xstored       
    logical :: isfirstguess
    integer :: info
    character(len=lenText) :: text, istr, rstr
    character(len=20) :: fname
    integer :: iend    

    ! .. executable statements 

    ! .. mpi

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

    ! .. logfile 
    write(istr,'(I2)')rank
    fname='status.'//trim(adjustl(istr))//'.log'
    call open_logfile(logUnit,fname) 
    write(istr,'(I4)')rank
    text='program begins : rank '//istr
    call print_to_log(LogUnit,text)


    ! .. init
    call read_inputfile(info)
    if(info/=0) then
        write(istr,'(I3)')info
        text="Error in input file: info = "//istr//" : end program."
        call print_to_log(LogUnit,text)
        print*,text
        call MPI_FINALIZE(ierr)
        stop
    endif


    call init_constants()

    call allocate_geometry(nx,ny,nz)
    call make_geometry()            ! generate volume elements lattice 
    call init_matrices()            ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC,ngr_node)  
    call make_chains(chainmethod)   ! generate polymer configurations 


    call make_sequence_chain(period,chaintype)
    call set_properties_chain(period,chaintype)  

    print*,"main: ierr=",ierr," rank=",rank
    

   
    call make_sequence_chain(period,chaintype)
    call set_properties_chain(period,chaintype)  

   
    ! call allocate_geometry(nx,ny,nz)
    !call make_geometry()            ! generate volume elements lattice 
    call allocate_field(nx,ny,nz) 
    call allocate_part_fnc(ngr)
    call set_size_neq()             ! number of non-linear equation neq  
    call set_fcn()     
    call init_surface(bcflag,nsurf)

    !  .. computation starts
    
    print*,"main: cuantasAB=",cuantasAB


    nz = nzmax                    
    neqmax = neq    
    allocate(xstored(neq))
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(fvec(neq))

    isfirstguess = .true.    
    use_xstored = .false.         ! with both flags set false make_guess will set xguess equal to x    
    
    pH%val=pH%min   
    call init_expmu()             ! set chemical potenitals  
    
    iter = 0
    print*,"main : nz=",nz,"nzmin=",nzmin," rank=",rank,"size=",size
    do while (nz>=nzmin)        ! loop distances

        print*,"main : inside while loop rank=",rank
        call set_size_neq()  
        
        if(.not.allocated(x)) allocate(x(neq))
        if(.not.allocated(xguess)) allocate(xguess(neq))
        if(.not.allocated(fvec)) allocate(fvec(neq))

        !call init_expmu   
        call init_vars_input()  ! sets up chem potenitals 
        call set_fcn()    
        
        flag_solver = 0
           
        if(rank.eq.0) then     ! node rank=0

            print*,"main: solver called rank=",rank," neqint=",neqint            
            call make_guess(x, xguess, isfirstguess, use_xstored, xstored)  
            call fcnptr(x, fvec, neq)  
            !  call solver(x, xguess, error, fnorm)

            flag_solver = 0   ! stop nodes
                
            do i = 1, size-1
                dest =i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
            enddo
        
        endif
  
        print*,"main: rank=",rank

        if(rank.ne.0) then 

            print*,"main: rank=",rank," neqint=",neqint  
            flag_solver = 1
 
            do while(flag_solver.eq.1) 
 
                flag_solver = 0
                source = 0 
                call MPI_RECV(flag_solver, 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
                if(flag_solver.eq.1) then
                    print*,"main : flag_solver=",flag_solver
                    call MPI_RECV(x, neqint, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)                        
                    call fcnptr(x, fvec, neq)  
                     print*,"main : fcnptr neq=",neq," neqint=",neqint
                endif    

            enddo  
        endif


        if(rank==0) then
            ! call copy_solution(x)
            !call compute_vars_and_output()
            call output()          

            isfirstguess =.false.    
            use_xstored = .true.
            iter = 0                ! reset of iteration counter 
            nz = nz-nzstep          ! reduce distance 
            do i=1,neq
                xstored(i)=x(i)
            enddo
            ! commucitate new values of nz (loop%val and loop%stepsize) from master to slave nodes
            ! to advance while loop on slave nodes
            do i = 1, size-1
                dest = i
                call MPI_SEND(nz, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
            !   call MPI_SEND(loop%val     , 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
            !   call MPI_SEND(loop%stepsize, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
            enddo
        else 
            ! receive values 
            source = 0
            call MPI_RECV(nz , 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
        endif
            
        iter  = 0              ! reset of iteration counter 
 
    enddo ! end while loop 

    call MPI_FINALIZE(ierr)  



    deallocate(x)   
    deallocate(xguess)
    deallocate(xstored)
    deallocate(fvec)
!    call deallocate_field()

    text="program end"

    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program brushweakpolyelectrolyte
