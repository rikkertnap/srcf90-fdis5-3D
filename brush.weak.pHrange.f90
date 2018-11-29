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
    use conform_entropy
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
    logical :: issolution
    integer :: info
    character(len=lenText) :: text, istr, rstr
    character(len=20) :: fname, conffilename
    integer :: iend , un_conf   
    type (looplist), pointer :: loop
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
    write(istr,'(A40)')VERSION
    text='program version     = '//istr
    call print_to_log(LogUnit,text)
    if(rank==0) print*,text

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
    ! call allocate_geometry(nx,ny,nz)
    call make_geometry()            ! generate volume elements lattice 
    call init_matrices()            ! init matrices for chain generation
    call allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC,ngr_node)  
    call make_sequence_chain(period,chaintype)
    call set_properties_chain(period,chaintype)      
    call make_chains(chainmethod)   ! generate polymer configurations    
    call allocate_field(nx,ny,nz) 
    call allocate_part_fnc(ngr)
    call init_field(Nx,Ny,Nz)
    call set_size_neq()             ! number of non-linear equation neq  
    call set_fcn()   
    call init_surface(bcflag,nsurf)
    call init_sigma()

    !  .. computation starts
       
    allocate(xstored(neq))
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(fvec(neq))

    if(runflag=="rangedist") then ! loop over distances
         
        nz = nzmax                    
        neqmax = neq   

        isfirstguess = .true.    
        use_xstored = .false.         ! with both flags set false make_guess will set xguess equal to x    

        pH%val=pH%min   
        call init_expmu()             ! set chemical potenitals     
        iter = 0
       
        do while (nz>=nzmin)        ! loop distances

            call set_size_neq()  
            
            if(.not.allocated(x)) allocate(x(neq))
            if(.not.allocated(xguess)) allocate(xguess(neq))
            if(.not.allocated(fvec)) allocate(fvec(neq))

            call init_vars_input()  ! sets up chem potenitals 
            call chain_filter()
            call set_fcn()! why 
            
            flag_solver = 0
               
            if(rank.eq.0) then     ! node rank=0    
                call make_guess(x, xguess, isfirstguess, use_xstored, xstored)  
                call solver(x, xguess, error, fnorm, issolution)
                flag_solver = 0   ! stop nodes
                do i = 1, size-1
                    dest =i
                    call MPI_SEND(flag_solver, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else
                flag_solver = 1
                do while(flag_solver.eq.1) 
                    flag_solver = 0
                    source = 0 
                    call MPI_RECV(flag_solver, 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
                    if(flag_solver.eq.1) then
                        call MPI_RECV(x, neqint, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)      
                        call fcnptr(x, fvec, neq)          
                    endif    
                enddo  
            endif

            call FEconf_entropy(FEconfAB,FEconfC) ! parrallel computation of conf entropy

            if(rank==0) then
                
                call compute_vars_and_output()
                
                isfirstguess =.false.    
                use_xstored = .true.
                iter = 0                ! reset of iteration counter 
                nz = nz-nzstep          ! reduce distance 
                do i=1,neq
                    xstored(i)=x(i)
                enddo
                ! communicate new values of nz from master to  compute  nodes to advance while loop on compute nodes
                do i = 1, size-1
                    dest = i
                    call MPI_SEND(nz, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else ! receive values 
                source = 0
                call MPI_RECV(nz , 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
            endif
                
            iter  = 0              ! reset of iteration counter 
            deallocate(x)   
            deallocate(xguess)
            deallocate(fvec)
     
        enddo ! end while loop 


    else if(runflag=="rangepH") then  ! loop over pH values 
        
        ! select pH variable with which loop is associated
        loop => pH 

        nz = nzmax                    
        neqmax = neq   
        isfirstguess = .true.    
        use_xstored = .false.         ! with both flags set false make_guess will set xguess equal to x          
        iter = 0                    ! iteration counter      

        if(loop%stepsize>0) then 
            loop%val=loop%min
        else
            loop%val=loop%max
        endif    

        do while (loop%min<=loop%val.and.loop%val<=loop%max.and.&
                (abs(loop%stepsize)>=loop%delta))

            call init_vars_input()  ! sets up chem potenitals 
            call chain_filter()
            call set_fcn()

            flag_solver = 0
           
            if(rank==0) then     ! node rank=0    
                call make_guess(x, xguess, isfirstguess)  
                call solver(x, xguess, error, fnorm, issolution)
                flag_solver = 0   ! stop nodes
                do i = 1, size-1
                    dest =i
                    call MPI_SEND(flag_solver, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else
                flag_solver = 1
                do while(flag_solver.eq.1) 
                    flag_solver = 0
                    source = 0 
                    call MPI_RECV(flag_solver, 1, MPI_INTEGER, source, tag,MPI_COMM_WORLD,stat, ierr)
                    if(flag_solver.eq.1) then
                        call MPI_RECV(x, neqint, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)      
                        call fcnptr(x, fvec, neq)          
                    endif    
                enddo  
            endif

            call FEconf_entropy(FEconfAB,FEconfC) ! parrallel computation of conf FEconf_entropy
            
            if(rank==0) then
            
                if(isSolution) then 
                    call compute_vars_and_output()
                    isfirstguess =.false.    
                    loop%val=loop%val+loop%stepsize   

                else
                    text="no solution: backstep" 
                    call print_to_log(LogUnit,text)
                    loop%stepsize=loop%stepsize/2.0d0   ! decrease increment 
                    loop%val=loop%val-loop%stepsize     ! step back
                    do i=1,neq
                        x(i)=xguess(i)
                    enddo       
                endif

                ! commucitate new values of loop from master to compute nodes to advance while loop on compute nodes
                do i = 1, size-1
                    dest = i
                    call MPI_SEND(loop%val, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD,ierr)
                enddo
            else ! receive values 
                source = 0
                call MPI_RECV(loop%val, 1, MPI_DOUBLE_PRECISION, source, tag,MPI_COMM_WORLD,stat, ierr)
            endif
                
            iter  = 0              ! reset of iteration counter 
            
         enddo ! end while loop 

    endif    

    call MPI_FINALIZE(ierr)  

    deallocate(xstored)
    call deallocate_field()

    text="program end"

    call print_to_log(LogUnit,text)
    call close_logfile(LogUnit)

end program brushweakpolyelectrolyte
