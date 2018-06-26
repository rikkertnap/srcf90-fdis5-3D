! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f90:                                          |       
! generator chains on planar surface                            |
! --------------------------------------------------------------|


module chaingenerator 

    use precision_definition   
    implicit none

    private :: pbc 

contains


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
  
    use mpivars
    use globals
    use chains
    use random
    use parameters, only : geometry, lsegAB
    use volume, only : ngrx, ngry, nx, ny,  nz, delta, ngr, ngr_node, ngr_freq
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use myutils
    use cadenas_linear
    use cadenas_sequence

    implicit none

    !     .. variable and constant declaractions      

    real(dp), dimension(:), allocatable :: x_ngr  ! x coordinate of  grafting point
    real(dp), dimension(:), allocatable :: y_ngr  

    integer :: i,j,k,s,g,gn      ! dummy indices
    integer :: idx               ! index label
    integer :: ix,iy,idxtmp
    integer :: nchains           ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: conf              ! counts number of conformations
    integer :: allowedconfAB
    real(dp) :: chain(3,nsegAB,200) ! chain(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp) :: xp(nsegAB), yp(nsegAB), zp(nsegAB) ! coordinates
    real(dp) :: x(nsegAB), y(nsegAB), z(nsegAB)
    real(dp) :: Lx,Ly,Lz         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    character(len=lenText) :: text, istr
    integer  :: xc,yc,zc        

    !     .. executable statements
    !     .. initializations of variables     
    
    allocate(x_ngr(ngrx))
    allocate(y_ngr(ngry))

    do i=1,ngrx            ! location graft points 
        x_ngr(i)=(i-0.5_dp)*delta*ngr_freq
    enddo
    do i=1,ngry
        y_ngr(i)=(i-0.5_dp)*delta*ngr_freq  
    enddo

    conf=1                  ! counter for conformations
    seed=435672*(rank+1)    ! seed for random number generator  different on each node
    maxnchains=12

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 

    if(isHomopolymer.eqv..FALSE.) then 
        allocate(lsegseq(nsegAB))
        call make_lsegseq(lsegseq,nsegAB)
    endif    


    do while (conf.le.max_conforAB)
        nchains= 0      ! init zero 
        if(isHomopolymer) then 
            call make_linear_chains(chain,nchains,maxnchains,nsegAB,lsegAB) ! chain generator f90
        else
            call make_linear_seq_chains(chain,nchains,maxnchains,nsegAB) 
        endif  

        select case (geometry) 
        case("cubic")
            
            do j=1,nchains   
            
                do s=1,nsegAB                      !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j) 
                enddo

                do gn=1,ngr_node                   ! loop over grafted points per node
               
                    g = rank*ngr_node +gn          ! g real number grafted postion gn relative to rank node    
                  
                    idxtmp = g
                    ix     = mod(idxtmp-1,ngrx)+1  ! inverse of g
                    idxtmp = int((idxtmp-1)/ngrx)+1
                    iy = idxtmp
               
                    xpt = x_ngr(ix)                ! position of graft point
                    ypt = y_ngr(iy) 
                      
                    weightchainAB(gn,conf)=.TRUE.  ! init weight

                    do s=1,nsegAB

                        ! .. translation onto correct grafting aream translation in xy plane 
                        x(s)=xp(s)+xpt
                        y(s)=yp(s)+ypt

                        ! .. check z coordinate 
                        if((0>z(s)).or.(z(s)>Lz)) weightchainAB(gn,conf)=.FALSE. 

                        ! .. periodic boundary conditions in x-direction and y-direction  
                        x(s)=pbc(x(s),Lx)
                        y(s)=pbc(y(s),Ly)

                        ! .. transforming form real- to lattice coordinates                 
                        xc=int(x(s)/delta)+1
                        yc=int(y(s)/delta)+1
                        zc=int(zp(s)/delta)+1

                        if(weightchainAB(gn,conf).eqv..TRUE.) then
                            call linearIndexFromCoordinate(xc,yc,zc,idx)
                            indexchainAB(s,gn,conf) = idx
                            if(idx<=0) then
                                print*,"index=",idx, " xc=",xc," yc=",yc," zc=",zc, "conf=",conf,"s=",s 
                            endif
                        endif

                    enddo   
                enddo         ! end g loop

                conf=conf +1
            
            enddo             ! end j loop

        case("square")
        

            do j=1,nchains   
            
                do s=1,nsegAB         !     transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j) 
                enddo

                do gn=1,ngr_node ! loop over grafted points per node
               
                    g = rank*ngr_node +gn          ! g real number grafted postion gn relative to rank node    
                    idxtmp = g
                    ix     = mod(idxtmp-1,ngrx)+1  ! inverse of g
                    xpt = x_ngr(ix)                ! position of graft point

                    weightchainAB(gn,conf)=.TRUE.   ! init weight

                    do s=1,nsegAB

                        ! .. translation onto correct grafting aream translation in xy plane 
                        x(s)=xp(s)+xpt
                   
                        ! .. check z coordinate 
                        if((0>z(s)).or.(z(s)>Lz)) weightchainAB(gn,conf)=.FALSE. 

                        ! .. periodic boundary conditions in x-direction and y-direction  
                        x(s)=pbc(x(s),Lx)
                   
                        ! .. transforming form real- to lattice coordinates                 
                        xc=int(x(s)/delta)+1
                        zc=int(zp(s)/delta)+1

                        if(weightchainAB(gn,conf).eqv..TRUE.) then
                            call linearIndexFromCoordinate(xc,1,zc,idx)
                            indexchainAB(s,gn,conf) = idx
                            if(idx<=0) then
                                print*,"index=",idx, " xc=",xc," yc=",yc," zc=",zc, "conf=",conf,"s=",s 
                            endif
                        endif

                    enddo   
                enddo         ! end g loop

                conf=conf +1
            
            enddo             ! end j loop

        case("hexagon")     

        case default
           print*,"Error: geometry not cubic, square, or hexagon"
           print*,"stopping program"
           stop
        end select
     

    enddo                     ! end while loop
             
    !     .. end chains generation 
      
    write(istr,'(I4)')rank
    text='AB Chains generated on node '//istr
    call print_to_log(LogUnit,text)
    print*,text

    do gn=1,ngr_node
        allowedconfAB=0
        do k=1,cuantasAB 
            if(weightchainAB(gn,k).eqv..TRUE.)then
                allowedconfAB=allowedconfAB+1
            endif
        enddo
        print*,"allowedconf=",allowedconfAB,"gn=",gn,"rank=",rank
    enddo
      
    deallocate(x_ngr)
    deallocate(y_ngr)
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)


end subroutine make_chains_mc


subroutine make_chains_file()
  
    !     .. variable and constant declaractions                                                           
    use globals
    use chains
    use random
    use parameters
    use volume
    use chain_rotation, only : rotation

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

    ! !     .. first line  contains nseg

    ! !      read(1,*)nsegfile
    ! !      read(1,*)cuantasfile

    ! nsegfile=50
    ! cuantasfile=1000000
    ! conf=0                    ! counter for conformations                   
    ! conffile=0                ! counter for conformations in file
    ! seed=435672               ! seed for random number generator 
    ! rotmax=12                 ! maximum number of rotation conf 
    ! maxattempts=72            ! maximum number of attempts to rotate conf

    ! if(nsegfile-1.ne.nsegAB) then
    !     print*,"nseg chain file not equal internal nseg : stop program"
    !     stop
    ! endif

    ! do while ((conf.le.max_conforAB).and.(conffile.le.cuantasfile))
     
    !     conffile=conffile +1
     
    !     read(1,*)x0,y0,z0      ! .. origin
     
    !     chain(1,1) = 0.0
    !     chain(2,1) = 0.0
    !     chain(3,1) = 0.0
     
    !     do s=2,nsegfile        ! .. read form file
    !         read(1,*)x,y,z
    !         chain(1,s) = x-x0
    !         chain(2,s) = y-y0
    !         chain(3,s) = z-z0
    !     enddo
         
    !     do rot=1,rotmax        !  ..  transforming form real- to lattice coordinates                                  
    !         rottest=.FALSE.
    !         nchain=1
            
    !         do while ((rottest.eqv. .false.).and.(nchain.lt.maxattempts)) 
    !             rottest=rotation(chain,chains_rot,nsegAB)
    !             nchain=nchain+1
    !         enddo
            
    !         if(rottest) then 
    !             conf=conf+1 
    !             do s=1,nsegAB 
    !                 z = chains_rot(1,s+1)
    !                 indexchainAB_init(conf,s)= int(z/delta)+1
    !             enddo
    !         endif
    !     enddo                  ! end rot loop  
         
    ! enddo                     ! end while loop                

    ! !     .. end chains generation 

    ! if(conf.le.max_conforAB) then
    !     print*,"Something went wrong"
    !     print*,"conf     = ",conf," max_conforAB    = ",max_conforAB
    !     print*,"conffile = ",conffile," cuantasfile = ",cuantasfile
    !     stop
    ! else
    !     print*,"Chains generated"
    !     readinchains=conffile
    ! endif

    ! close(1)

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
    integer :: count_seg 
    integer :: indx  ! temporary index of chain
    integer :: ix,iy,iz,gn

    allowed_confAB=0            ! counts allowed conformations 
    do conf=1,max_conforAB      ! loop of all polymer conformations to filter out allowed ones 
        do gn=1,ngr_node        ! loop over grafted points per node
            count_seg=0
            do s=1,nsegAB
                indx=indexchainAB_init(s,gn,conf)
                call coordinateFromLinearIndex(indx,ix,iy,iz)
                if(iz<=nz) then  ! nz is distance in between plates  
                    indexchainAB(s,gn,allowed_confAB)=indx
                    count_seg=count_seg+1
                endif
            enddo
            if (count_seg.eq.nsegAB) allowed_confAB = allowed_confAB+1 ! conformation  is allowed  
        enddo
    enddo

    cuantasAB=allowed_confAB    ! the number of allowed conformations 

    allowed_confC=0             ! counts allowed conformations 
    do conf=1,max_conforC       ! loop of all polymer conformations to filter out allowed ones 
        do gn=1,ngr_node        ! loop over grafted points per node
            count_seg=0
            do s=1,nsegC
                indx=indexchainC_init(s,gn,conf)
                call coordinateFromLinearIndex(indx,ix,iy,iz)
                if(iz<=nz) then  ! nz is distance in between plates  
                    indexchainC(s,gn,allowed_confC)=indx
                    count_seg=count_seg+1
                endif
            enddo
            if (count_seg.eq.nsegC) allowed_confC = allowed_confC+1 ! conformation  is allowed  
        enddo
    enddo

    cuantasC=allowed_confC     ! the number of allowed conformations

end subroutine  chain_filter



! set isHomopolymer and set proper segment lenght homopolymer 

subroutine set_properties_chain(freq,chaintype)

    use globals
    use parameters, only : lsegAB, lsegPAA, lsegPAMPS
    use chains


    implicit none

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype

    if(freq>nsegAB) then 
        isHomopolymer=.TRUE.
    else 
        if(freq==0.and.chaintype.eq."diblock") then 
            isHomopolymer=.TRUE.
        else
            isHomopolymer=.FALSE.
        endif      
    endif    

    if(isHomopolymer) then 
        if (isAmonomer(1)) then ! it is a homopolymer so all segments the same 
            lsegAB=lsegPAA
        else
            lsegAB=lsegPAMPS
        endif   
    endif    

end subroutine set_properties_chain



! .. compute periodic boundary condition in z direction
! .. maps z onto interval [0,deltaz]

function pbc(z,deltaz) result(zpbc)
    implicit none 
    real(dp), intent(in) :: deltaz
    real(dp), intent(in) :: z
    real(dp) :: zpbc 
    if(z>0) then 
        zpbc=z-int(z/deltaz)*deltaz
    else
        zpbc=z-(int(z/deltaz)-1)*deltaz
    endif   
end function


end module chaingenerator