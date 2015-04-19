! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f:                                            |       
! generator chains on planar surface                            |
! --------------------------------------------------------------|



subroutine make_chains(chainmethod)

    implicit none

    character(len=8) :: chainmethod

    if(chainmethod.eq.'MC') then
        call make_chains_mc()
    elseif(chainmethod.eq.'FILE') then
        call make_chains_file()
    else
        print*,"chainmethod not equalt to MC or FILE"
        stop
    endif

end subroutine make_chains


!     purpose: init of cuantas polymer
!     configurations of polymer chain anchored onto a spherical surface 

subroutine make_chains_mc()
  
    use globals
    use chains
    use random
    use parameters
    use volume


    implicit none 

    !     .. variable and constant declaractions      

    integer :: j,s               ! dummy indices
    integer :: nchains            ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: conf              ! counts number of conformations
    real(dp)  :: chain(3,nsegAB,200) ! chains(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp)  :: x,y,z,r           ! coordinates
    character(len=50) :: text 

    !     .. executable statements 

    !     .. initializations of variables     
    !     .. init of  cuantas polymer configurations  
    !     .. anchoring of polymer chain onto spherical surface

    conf=0                    ! counter for conformations
    seed=435672               ! seed for random number generator 
    maxnchains=12

    do while (conf.le.max_conforAB)
        nchains= 0      ! init zero 
        !call cadenas(chain,nsegAB,lsegAB,nchain)  ! chain generator ! f77
        call cadenas(chain,nchains,maxnchains,nsegAB,lsegAB)  ! chain generator f90
        do j=1,nchains   
            conf=conf +1
            do s=1,nsegAB         !     transforming form real- to lattice coordinates
                z = chain(1,s,j)
                indexchainAB_init(conf,s)= int(z/delta)+1
!                print*,"indexchain(",conf,",",s,")=",indexchainAB_init(conf,s)  
            enddo
         enddo                  ! end j loop
    enddo                     ! end while loop
             
    !     .. end chains generation 
     
    write(text,'(A19)')'AB Chains generated'
    call print_to_log(LogUnit,text)

    conf=0                    ! counter for conformations                                                                      
    seed=43567                ! seed for random number generator 
    maxnchains=12

    do while (conf.le.max_conforC)
        nchains= 0      ! init zero 
        call cadenas(chain,nchains,maxnchains,nsegC,lsegC) ! chain generator                                                               
        do j=1,nchains
            conf=conf +1
            do s=1,nsegC         !     transforming form real- to lattice coordinates                                            
                z = chain(1,s,j)
                indexchainC_init(conf,s)= int(z/delta)+1
            enddo
        enddo                  ! end j loop  
    enddo                     ! end while loop  

    write(text,'(A18)')'C Chains generated'
    call print_to_log(LogUnit,text)

end subroutine make_chains_mc


subroutine make_chains_file()
  
    !     .. variable and constant declaractions                                                           
    use globals
    use chains
    use random
    use parameters
    use volume

    implicit none

    integer :: j,s,rot           ! dummy indices 
    integer :: conf,conffile     ! counts number of conformations
    integer :: nsegfile          ! nseg in input file
    integer :: cuantasfile       ! cuantas in input file
    real(dp)  :: chain(3,nsegAB+1)   ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp)  :: chains_rot(3,nsegAB+1) ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp)  :: x,y,z,r            ! coordinates
    real(dp)  ::  x0,y0,z0           ! origin coordinates

    character(len=10) :: fname
    integer :: ios
    logical :: rottest
    integer :: nchain,rotmax,maxattempts

    !     .. executable statements                                                                          
    !     .. reading in chains  from file  

    write(fname,'(A10)')'chains.dat'
    open(unit=1,file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        stop
    endif

    !     .. first line  contains nseg

    !      read(1,*)nsegfile
    !      read(1,*)cuantasfile

    nsegfile=50
    cuantasfile=1000000
    conf=0                    ! counter for conformations                   
    conffile=0                ! counter for conformations in file
    seed=435672               ! seed for random number generator 
    rotmax=12                 ! maximum number of rotation conf 
    maxattempts=72            ! maximum number of attempts to rotate conf

    if(nsegfile-1.ne.nsegAB) then
        print*,"nseg chain file not equal internal nseg : stop program"
        stop
    endif

    do while ((conf.le.max_conforAB).and.(conffile.le.cuantasfile))
     
        conffile=conffile +1
     
        read(1,*)x0,y0,z0      ! .. origin
     
        chain(1,1) = 0.0
        chain(2,1) = 0.0
        chain(3,1) = 0.0
     
        do s=2,nsegfile        ! .. read form file
            read(1,*)x,y,z
            chain(1,s) = x-x0
            chain(2,s) = y-y0
            chain(3,s) = z-z0
        enddo
         
        do rot=1,rotmax        !  ..  transforming form real- to lattice coordinates                                  
            rottest=.FALSE.
            nchain=1
            
            do while ((rottest.eqv. .false.).and.(nchain.lt.maxattempts)) 
                call rotation(chain,chains_rot,nsegAB,rottest,lsegAB)
                nchain=nchain+1
            enddo
            
            if(rottest) then 
                conf=conf+1 
                do s=1,nsegAB 
                    z = chains_rot(1,s+1)
                    indexchainAB_init(conf,s)= int(z/delta)+1
                enddo
            endif
        enddo                  ! end rot loop  
         
    enddo                     ! end while loop                

    !     .. end chains generation 

    if(conf.le.max_conforAB) then
        print*,"Something went wrong"
        print*,"conf     = ",conf," max_conforAB    = ",max_conforAB
        print*,"conffile = ",conffile," cuantasfile = ",cuantasfile
        stop
    else
        print*,"Chains generated"
        readinchains=conffile
    endif

    close(1)

end subroutine make_chains_file


subroutine make_sequence_chain(freq,chaintype)
  
    use globals
    use chains

    implicit none

    integer :: freq
    character(len=8)  :: chaintype

    !     .. local variables

    integer :: s

    if(chaintype.eq.'altA') then
        do s=1,nsegAB
            if(mod(s,freq).ne.0) then ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            endif
        enddo
    else if(chaintype.eq.'altB') then
        do s=1,nsegAB
            if(mod(s,freq).eq.0) then ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            endif
        enddo
    else if(chaintype.eq.'diblock') then
        do s=1,nsegAB
            if(s.le.freq) then   ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            endif
        enddo
    else
        print*,"Wrong chaintype: aborting program"
        stop
    endif
  
end subroutine make_sequence_chain



subroutine chain_filter()
    
    use  globals
    use  chains
    use  random
    use  parameters
    use  volume

    implicit none

    integer :: conf,s,allowed_confAB,allowed_confC
    integer :: flag
    integer(2) :: ind  ! temporary index of chain


    allowed_confAB=1            ! counts allowed conformations 
    do conf=1,max_conforAB  ! loop of all polymer conformations to filter out allowed ones 
        flag=0
        do s=1,nsegAB
            ind=indexchainAB_init(conf,s)
            if(ind<=nz) then   ! nz is in between plates  
                indexchainAB(allowed_confAB,s)=ind
                flag=flag+1
            endif
        enddo
        if (flag.eq.nsegAB) allowed_confAB= allowed_confAB+1 ! conformation  is allowed  
    enddo

    cuantasAB=allowed_confAB-1    ! the number of allowed conformations

    allowed_confC=1            ! counts allowed conformations 
    do conf=1,max_conforC  ! loop of all polymer conformations to filter out allowed ones 
        flag=0
        do s=1,nsegC
            ind=indexchainC_init(conf,s)
            if(ind<=nz) then   ! nz is in between plates  
                indexchainC(allowed_confC,s)=ind
                flag=flag+1
            endif
        enddo
        if (flag.eq.nsegC) allowed_confC= allowed_confC+1 ! conformation  is allowed  
    enddo

    cuantasC=allowed_confC-1    ! the number of allowed conformations

end subroutine  chain_filter


