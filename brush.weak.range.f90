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
    
    implicit none  
    
    real(dp),  dimension(:), allocatable :: x         ! iteration vector 
    real(dp),  dimension(:), allocatable :: xguess    ! guess iteration vector
    real(dp),  dimension(:), allocatable :: xstored   ! stored iteration vector
    real(dp),  dimension(:), allocatable :: fvec     
  
    integer :: i,c,eps          ! dummy indices
    integer, parameter :: cmax=20
    integer, parameter :: epsmax=19
    real(dp)  :: cs(cmax)         ! salt concentrations     
    real(dp)  :: VdWeps(epsmax)   ! VdW eps values
    integer :: countfile        ! file counter        
    logical :: use_xstored       
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
    call set_size_neq()         ! number of non-linear equation neq
    
    allocate(x(neq))
    allocate(xguess(neq))
    allocate(xstored(neq))
    allocate(fvec(neq))
    
    !     .. computation starts
    
    cs(13)=0.01d0
    cs(12)=0.05d0
    cs(11)=0.10d0
    cs(10)=0.15d0
    cs(9)=0.20d0
    cs(8)=0.25d0
    cs(7)=0.50d0
    cs(6)=0.75d0
    cs(5)=1.00d0
    cs(4)=1.25d0
    cs(3)=1.50d0
    cs(2)=1.75d0
    cs(1)=2.00d0
        
    VdWeps(1)=0.00d0
    VdWeps(2)=0.05d0
    VdWeps(3)=0.10d0
    VdWeps(4)=0.15d0
    VdWeps(5)=0.20d0
    VdWeps(6)=0.25d0
    VdWeps(7)=0.50d0
    VdWeps(8)=0.75d0
    VdWeps(9)=1.00d0
    VdWeps(10)=1.25d0
    VdWeps(11)=1.50d0
    VdWeps(12)=1.75d0
    VdWeps(13)=2.00d0
    VdWeps(14)=2.5d0
    VdWeps(15)=3.0d0
    VdWeps(16)=3.5d0
    VdWeps(17)=4.0d0
    VdWeps(18)=4.5d0
    VdWeps(19)=5.0d0


    if(sysflag.eq."elect") then
        ! .. salt concentrations
    
        countfile= 1    
        use_xstored=.false.

        do c=1,1                ! loop over salt concentration
            do  eps=1,epsmax               ! loop over eps values
         
                cNaCl=cs(c)
                VdWepsC=VdWeps(eps)                
                call init_expmu()
                call make_guess(x,xguess,VdWeps(eps),VdWeps(1),use_xstored,xstored)
                call solver(x, xguess, error, fnorm) 
                do i=1,nz                
                    xsol(i)=x(i)        ! solvent volume fraction   
                    psi(i) =x(i+nz)     ! potential           
                enddo
           
                call fcnenergy()        ! free energy
                call average_height()      
                call charge_polymer()
                call average_charge_polymer()
                call output(countfile)    ! writing of output
            
                countfile = countfile+1  ! next  
                iter  = 0                ! reset of iteration counter 
               
                if(eps==1) then
                    do i=1,neq
                        xstored(i)=x(i)
                    enddo
                endif
                if(eps==epsmax) then 
                    use_xstored=.true.
                else
                    use_xstored=.false.
                endif     
            enddo
        enddo

    elseif (sysflag.eq."neutral") then

        countfile= 1    

        do eps=1,epsmax               ! loop over salt concentration
         
            VdWepsB=VdWeps(eps)

            call init_expmu()
            call make_guess(x,xguess,VdWeps(eps),VdWeps(1))

            call solver(x, xguess, error, fnorm) 
            
            do i=1,nz           
                xsol(i)=x(i)        ! solvent volume fraction   
            enddo
           
            call fcnenergy()        ! free energy 
            call average_height()      
            call output(countfile)   ! writing of output
            
            countfile = countfile+1        ! next  
            iter  = 0              ! reset of iteration counter 

        enddo        
    else
        print*,"Wrong value sysflag : ", sysflag
        stop    
    endif   
    
    deallocate(x)
    deallocate(xguess)
    deallocate(xstored)
    deallocate(fvec)
    call deallocate_field()

    stop
end program brushweakpolyelectrolyte
