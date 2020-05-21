! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f90:                                          |       
! generator chains on planar surface                            |
! --------------------------------------------------------------|


module chaingenerator 

    use precision_definition   
    implicit none

    integer, parameter :: lenfname=40

    private :: lenfname
    private :: pbc 

contains


subroutine make_chains(chainmethod)

    use myutils,  only :  print_to_log, LogUnit, lenText

    character(len=15), intent(in)  :: chainmethod

    integer :: i
    character(len=lenText) :: text, istr


    select case (chainmethod)
    case ("MC")
        call make_chains_mc()
    case ("FILE_lammps_XYZ") 
        call read_chains_lammps_xyz  
    case ("FILE_TXT") 
        call make_chains_file()
    case default
        text="chainmethod not equalt to MC or FILE"
        call print_to_log(LogUnit,text)
        print*,text
        stop
    end select

end subroutine make_chains


!     purpose: init of cuantas polymer
!     configurations of polymer chain anchored onto a spherical surface 

subroutine make_chains_mc()
  
    use mpivars
    use globals
    use chains
    use random
    use parameters, only : geometry, lseg
    use volume, only : ngrx, ngry, nx, ny, nz, delta, ngr, ngr_node, ngr_freq 
    use volume, only : position_graft
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use volume, only : ut, vt
    use myutils
    use cadenas_linear
    use cadenas_sequence

    !     .. variable and constant declaractions      

    real(dp), dimension(:), allocatable :: x_ngr  ! x coordinate of  grafting point
    real(dp), dimension(:), allocatable :: y_ngr  

    integer :: i,j,k,s,g,gn      ! dummy indices
    integer :: idx               ! index label
    integer :: ix,iy,idxtmp,ntheta
    integer :: nchains           ! number of rotations
    integer :: maxnchains        ! number of rotations
    integer :: maxntheta         ! maximum number of rotation in xy-plane
    integer :: conf              ! counts number of conformations
    integer :: allowedconf
    real(dp) :: chain(3,nseg,200) ! chain(x,i,l)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp) :: x(nseg), y(nseg), z(nseg) ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates
    real(dp) :: xpp(nseg), ypp(nseg), zpp(nseg)  
    real(dp) :: Lx,Ly,Lz         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    character(len=lenText) :: text, istr
    integer  :: xc,yc,zc        

    !     .. executable statements
    !     .. initializations of variables     
       

    conf=1                  ! counter for conformations
    seed=435672             ! *(rank+1) ! seed for random number generator  different on each node
    maxnchains=12
    maxntheta =6            ! maximum number of rotation in xy-plane  

    theta_angle= 2.0_dp*pi/maxntheta
    
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 

    if(isHomopolymer.eqv..FALSE.) then 
        allocate(lsegseq(nseg))
        call make_lsegseq(lsegseq,nseg)
    endif    
  
    do while (conf<=max_confor)
        nchains= 0      ! init zero 
        if(isHomopolymer) then 
            call make_linear_chains(chain,nchains,maxnchains,nseg,lseg) ! chain generator f90
        else
            call make_linear_seq_chains(chain,nchains,maxnchains,nseg) 
        endif  

        if(geometry=="cubic") then
            
            do j=1,nchains   
            
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j) 
                enddo

                do ntheta=1,maxntheta                  ! rotation in xy-plane

                    theta= ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s)= xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s)=-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    

                    do gn=1,ngr_node                   ! loop over grafted points per node
                   
                        g = rank*ngr_node +gn          ! g real number grafted postion gn relative to rank node    
                      
                        idxtmp = g
                        ix = mod(idxtmp-1,ngrx)+1      ! inverse of g
                        iy = int((idxtmp-1)/ngrx)+1
                        
                        xpt =  position_graft(g,1)  !x_ngr(ix)                ! position of graft point
                        ypt =  position_graft(g,2)  !y_ngr(iy) 
                        
                        weightchain(gn,conf)=.TRUE.  ! init weight     

                        do s=1,nseg

                            ! .. translation onto correct grafting area translation in xy plane 
                            x(s)=xpp(s)+xpt
                            y(s)=ypp(s)+ypt

                            ! .. check z coordinate 
                            if((0>zp(s)).or.(zp(s)>Lz)) then 
                                weightchain(gn,conf)=.FALSE. 
                                print*,"conf=",conf,"s=",s,"zp=",zp(s)
                            endif   

                            ! .. periodic boundary conditions in x-direction and y-direction  
                            x(s)=pbc(x(s),Lx)
                            y(s)=pbc(y(s),Ly)

                            ! .. transforming form real- to lattice coordinates                 
                            xc=int(x(s)/delta)+1
                            yc=int(y(s)/delta)+1
                            zc=int(zp(s)/delta)+1

                            if(weightchain(gn,conf).eqv..TRUE.) then
                                call linearIndexFromCoordinate(xc,yc,zc,idx)
                                indexchain_init(s,gn,conf) = idx
                                if(idx<=0) then
                                    print*,"index=",idx, " xc=",xc," yc=",yc," zc=",zc, "conf=",conf,"s=",s 
                                endif
                            endif

                        enddo  
                                       
                    enddo     ! end loop over graft points

                    conf=conf +1

                enddo         ! end loop over rotations
            
            enddo             ! end loop over nchains

        else if(geometry=="prism") then
            
            do j=1,nchains   
            
                do s=1,nseg                     !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j) 
                enddo

                do ntheta=1,maxntheta                  ! rotation in xy-plane

                    theta= ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s)= xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s)=-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    


                    do gn=1,ngr_node                   ! loop over grafted points per node
                   
                        g = rank*ngr_node +gn          ! g real number grafted postion gn relative to rank node    
                         
                        xpt = position_graft(g,1)
                        ypt = position_graft(g,2)
                        
                        weightchain(gn,conf)=.TRUE.  ! init weight     

                        do s=1,nseg

                            ! .. translation onto correct grafting area translation in xy plane 
                            x(s)=ut(xpp(s),ypp(s))+xpt
                            y(s)=vt(xpp(s),ypp(s))+ypt

                            ! .. check z coordinate 
                            if((0>zp(s)).or.(zp(s)>Lz)) weightchain(gn,conf)=.FALSE. 

                            ! .. periodic boundary conditions in x-direction and y-direction  
                            x(s)=pbc(x(s),Lx)
                            y(s)=pbc(y(s),Ly)

                            ! .. transforming form real- to lattice coordinates                 
                            xc=int(x(s)/delta)+1
                            yc=int(y(s)/delta)+1
                            zc=int(zp(s)/delta)+1

                            if(weightchain(gn,conf).eqv..TRUE.) then
                                call linearIndexFromCoordinate(xc,yc,zc,idx)
                                indexchain_init(s,gn,conf) = idx
                                if(idx<=0) then
                                    print*,"index=",idx, " xc=",xc," yc=",yc," zc=",zc, "conf=",conf,"s=",s 
                                endif
                            endif

                        enddo  
                                       
                    enddo         ! end loop over graft points

                    conf=conf +1  ! end loop over rotations

                enddo     
            
            enddo                 ! end loop over nchains

        else 
           print*,"Error: in make_chains_mc: geometry not cubic, prism, or square"
           print*,"stopping program"
           stop
        endif
     
    enddo                     ! end while loop
            
    !     .. end chains generation 
      
    write(istr,'(I4)')rank
    text='AB Chains generated on node '//istr
    call print_to_log(LogUnit,text)
    !print*,text

    do gn=1,ngr_node
        allowedconf=0
        do k=1,cuantas
            if(weightchain(gn,k).eqv..TRUE.)then
                allowedconf=allowedconf+1
            endif
        enddo
      
        if(allowedconf/=cuantas) then 
            print*,"allowedconf=",allowedconf,"gn=",gn,"rank=",rank
            print*,"make_chains_mc: Not all conformations were allowed"    
        endif

    enddo
        
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)

end subroutine make_chains_mc




! Reads confomations from a file called traj.xyz
! Format lammps trajectory file in xyz format
! number of ATOMS much equal nseg+1  
! differeent processor rank get assigned different conformations

subroutine read_chains_lammps_XYZ(info)

    !     .. variable and constant declaractions  
    use mpivars                                                                                     
    use globals
    use chains
    use random
    use parameters
    use volume
    use chain_rotation
    use myio, only : myio_err_chainsfile, myio_err_energyfile
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit

    ! .. argument

    integer, intent(out),optional :: info

    ! .. local variables

    integer :: i,j,s,rot,g,gn      ! dummy indices
    integer :: idx                 ! index label
    integer :: ix,iy,idxtmp,ntheta
    integer :: nchains              ! number of rotations
    integer :: maxnchains           ! number of rotations
    integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg+1)     ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1  
    real(dp) :: chain_rot(3,nseg+1) 
    real(dp) :: xseg(3,nseg+1)
    real(dp) :: x(nseg), y(nseg), z(nseg) ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: Lx,Ly,Lz         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    real(dp) :: xc,yc,zc          
    integer  :: xi,yi,zi          
    real(dp) :: energy                                               
    character(len=8) :: fname
    integer :: ios 
    character(len=30) :: str
    integer :: nchain, rotmax, maxattempts, idatom
    integer :: un, un_ene ! unit number
    logical :: is_positive_rot, exist
    character(len=lenText) :: text,istr 


    ! .. executable statements                                                                                                     
    ! .. open file                                                                                              

    write(fname,'(A8)')'traj.xyz' ! trajectory coordinates
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj.xyz file does not exit'
        if (present(info)) info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    write(fname,'(A8)')'traj.ene' ! traectory energy
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
    else
        print*,'traj.ene file does not exit'
        if (present(info)) info = myio_err_energyfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_energyfile
        return
    endif

    conf=0                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file                                                                    
    seed=435672               ! seed for random number generator    
    maxnchains=12                                                                         
    maxattempts=72            ! maximum number of attempts to rotate conf                                                                   
    !rotmax=maxnchainsrotations  ! maximum number of rotation conf  
    rotmax=maxnchains
    ios=0

    maxntheta =6                ! maximum number of rotation in xy-plane  
    theta_angle= 2.0_dp*pi/maxntheta    
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 


    do while ((conf<cuantas).and.(ios==0))
    
        if(conf.ne.0) then ! skip lines
            read(un,*,iostat=ios)str
            read(un,*,iostat=ios)str
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            read(un,*,iostat=ios)
            if(nsegfile.ne.(nseg+1)) then ! nsegfile =nseg+1 . The chain has one 'segement' more: tethered point.
                print*,"nseg chain file not equal internal (nseg+1) : stop program"
                stop
            endif    
        endif
        
        do s=1,nseg+1              ! .. read form  trajecotory file
            read(un,*,iostat=ios)idatom,xc,yc,zc
            xseg(1,s) = xc/10.0_dp  ! .. convert from Angstrom to nm
            xseg(2,s) = yc/10.0_dp
            xseg(3,s) = zc/10.0_dp 
            !print*,"index=",idatom, " xc=",xc," yc=",yc," zc=",zc, "conf=",conf,"s=",s 
        enddo

        read(un_ene,*,iostat=ios)energy

        conffile=conffile +1 


        if(mod(conffile,size)==rank) then 
            ! .. assignes this conformation to process numberrank
            ! .. 'first' segment, sets origin 

            do s=1,nseg+1        
                chain(2,s) = (xseg(1,s)-xseg(1,1)) 
                chain(3,s) = (xseg(2,s)-xseg(2,1))
                chain(1,s) = (xseg(3,s)-xseg(3,1))
                !print*,"conf=",conf,"s=",s,(chain(i,s),i=1,3)
            enddo
 
            do rot=1,rotmax        
                 
                nchain=1
                is_positive_rot=.false.

                do while ((is_positive_rot.eqv..false.).and.(nchain.lt.maxattempts))
                    is_positive_rot=rotation(chain,chain_rot,nseg)
                    nchain=nchain+1
                enddo
                !print*,"conf=",conf,"is_positve=",is_positive_rot

                if(is_positive_rot) then
            
                    if(geometry=="cubic") then !  transforming form real- to lattice coordinates
            
                        do ntheta=1,maxntheta                  ! rotation in xy-plane

                            conf=conf+1
                            theta= ntheta * theta_angle

                            do s=1,nseg
                                xp(s)= chain_rot(2,s+1)*cos(theta)+chain_rot(3,s+1)*sin(theta) 
                                yp(s)=-chain_rot(2,s+1)*sin(theta)+chain_rot(3,s+1)*cos(theta)
                                zp(s)= chain_rot(1,s+1)   
                            enddo    

                            do gn=1,ngr_node                   ! loop over grafted points per node
                               
                                g = rank*ngr_node +gn          ! g real number grafted postion gn relative to rank node    
                               
                                xpt =  position_graft(g,1)     ! position of graft point
                                ypt =  position_graft(g,2)   
                                
                                weightchain(gn,conf)=.TRUE.  ! init weight     

                                do s=1,nseg

                                    ! .. translation onto correct grafting area translation in xy plane 
                                    x(s)=xp(s)+xpt
                                    y(s)=yp(s)+ypt
                                    ! .. check z coordinate 
                                    if((0>zp(s)).or.(zp(s)>Lz)) weightchain(gn,conf)=.FALSE. 
                                    ! .. periodic boundary conditions in x-direction and y-direction  
                                    x(s)=pbc(x(s),Lx)
                                    y(s)=pbc(y(s),Ly)

                                    ! .. transforming form real- to lattice coordinates                 
                                    xi=int(x(s)/delta)+1
                                    yi=int(y(s)/delta)+1
                                    zi=int(zp(s)/delta)+1

                                    if(weightchain(gn,conf).eqv..TRUE.) then
                                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                                        indexchain_init(s,gn,conf) = idx
                                        if(idx<=0) then
                                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                        endif
                                        exp_energychain(conf)=exp(-energy)
                                    endif
                                enddo                    
                            enddo     ! end loop over graft points

                        enddo         ! end loop over rotations
            
                    else if(geometry=="prism") then
                     ! to be done later 

                    else 
                        print*,"Error: in make_chains_mc: geometry not cubic, prism, or square"
                        print*,"stopping program"
                        stop
                    endif
                endif  ! if 
            enddo      ! rotation loop
        endif          ! if distrubtion if 
    enddo              ! end while loop                                                                                                          
    !  .. end chains generation    
    

    if(conf<=cuantas) then
        print*,"subroutine make_chains_lammps_XYZ : something went wrong"
        print*,"conf     = ",conf," cuantas     = ",cuantas
        print*,"conffile = ",conffile
        stop
    else
        text="Chains generated: subroutine make_chains_lammps_file"
        call print_to_log(LogUnit,text)
        readinchains=conffile
    endif

    close(un)
    close(un_ene)

end subroutine read_chains_lammps_XYZ

subroutine make_chains_file()
  
    !     .. variable and constant declaractions                                                           
    use globals
    use chains
    use random
    use parameters
    use volume
    use chain_rotation, only : rotation


    integer :: j,s,rot           ! dummy indices 
    integer :: conf,conffile     ! counts number of conformations
    integer :: nsegfile          ! nseg in input file
    integer :: cuantasfile       ! cuantas in input file
    real(dp)  :: chain(3,nseg+1)   ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
    real(dp)  :: chains_rot(3,nseg+1) ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1
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

end subroutine make_chains_file

! post: isAmonomer set 
! for chaintype==multi , type_of_monomer,type_of_monomer_char
! ismonomer_of_type are set.

subroutine make_sequence_chain(freq,chaintype)
  
    use globals
    use chains
    use parameters, only : typesfname

    integer, intent(in) :: freq
    character(len=8),intent(in)  :: chaintype

    !     .. local variables
    integer :: info
    integer :: s

    select case (chaintype ) 
    case ('altA')
        do s=1,nseg
            if(mod(s,freq).ne.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            endif
        enddo
    case('altB') 
        do s=1,nseg
            if(mod(s,freq).eq.0) then ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            endif
        enddo
    case('diblockA')
        do s=1,nseg
            if(s<=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            endif
        enddo
    case('diblockB') 
        do s=1,nseg
            if(s>=freq) then   ! A segment
                isAmonomer(s)=.TRUE.
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                isAmonomer(s)=.FALSE.
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            endif
        enddo  
    case('copolyAB')

        call read_sequence_copoly_from_file(info)

         do s=1,nseg
           if(isAmonomer(s)) then 
                type_of_monomer_char(s)="A"
                type_of_monomer(s)=1
            else
                type_of_monomer_char(s)="B"
                type_of_monomer(s)=2
            endif
        enddo

    case('multi')
    
        call read_type_of_monomer(type_of_monomer,type_of_monomer_char,typesfname, nseg) 
        call make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)    
        
        do s=1,nseg
            isAmonomer(s)=(type_of_monomer_char(s)=="A")
        enddo

    case default
        print*,"Wrong chaintype: aborting program"
        stop
    end select

end subroutine make_sequence_chain


subroutine read_sequence_copoly_from_file (info)

    use globals, only : nseg
    use chains,  only : isAmonomer
    use myutils, only : newunit

    integer, intent(out),optional :: info

    character(len=11) :: fname
    integer :: ios, un_seq, s
    character :: char

    if (present(info)) info = 0

    write(fname,'(A11)')'sequence.in'
    open(unit=newunit(un_seq),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file sequence.in : iostat =', ios
        stop
    endif
        
    s=0
    ios=0

    do while (s<nseg.and.ios==0)

        read(un_seq,*,iostat=ios)char

        if (ios==0) then
            s=s+1 
            if(char=="A") then   ! A segment
                isAmonomer(s)=.TRUE.
            else
                isAmonomer(s)=.FALSE.
            endif
        endif    
    enddo 
    if(s/=nseg) then 
        print*,"reached end of file before all elements read"
        info = 1
        stop "read sequence file failed"
    endif

end subroutine read_sequence_copoly_from_file


subroutine chain_filter()
    
    use  globals
    use  chains
    use  random
    use  parameters
    use  volume

    implicit none

    integer :: conf,s,allowed_conf
    integer :: count_seg 
    integer :: indx  ! temporary index of chain
    integer :: ix,iy,iz,gn

    do gn=1,ngr_node              ! loop over grafted points per node
        allowed_conf=1            ! counts allowed conformations 
        do conf=1,max_confor      ! loop of all polymer conformations to filter out allowed ones 
            count_seg=0
            do s=1,nseg
                indx=indexchain_init(s,gn,conf)
                call coordinateFromLinearIndex(indx,ix,iy,iz)
                if(iz<=nz) then  ! nz is distance in between plates  
                    indexchain(s,gn,allowed_conf)=indx
                    count_seg=count_seg+1
                endif
            enddo
            if (count_seg.eq.nseg) allowed_conf = allowed_conf+1 ! conformation  is allowed  
        enddo
    enddo

    cuantas=allowed_conf-1    ! the number of allowed conformations 

end subroutine  chain_filter


subroutine set_properties_chain(freq,chaintype)
    
    use globals
    use chains
    use parameters, only : lseg

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype

    call set_isHomopolymer(freq,chaintype)
    call set_lsegAA()

    if(isHomopolymer) lseg=set_lseg_Homopolymer()

end subroutine set_properties_chain

! set logical isHomoolymer 
! pre : type_of_monomer need to be initialized

subroutine set_isHomopolymer(freq,chaintype)

    use globals, only : nseg
    use chains, only : type_of_monomer, isHomopolymer

    integer, intent(in) :: freq
    character(len=8), intent(in)   :: chaintype
    
    integer :: s,  type_number

    if(chaintype=="multi") then
        type_number=type_of_monomer(1)
        do s=2,nseg
            isHomopolymer=(type_of_monomer(s)==type_number)
        enddo
    else if(freq>nseg) then 
        isHomopolymer=.TRUE.
    else 
        if(freq==0.and.(chaintype.eq."diblockA".or. chaintype=="diblockB")) then 
            isHomopolymer=.TRUE.
        else
            isHomopolymer=.FALSE.
        endif      
    endif    

end subroutine set_isHomopolymer


!  Assigns specific values to lsegAA initial set in init_lseg in parameters

subroutine set_lsegAA

    use globals
    use chains
    use parameters, only : lseg, lsegAA,lsegPAA, lsegPAMPS, lsegPEG
    use parameters, only : chainmethod 

    if(chainmethod=='MC') then  ! chain are not read in from file 
        
        select case (systype)
        case ("elect","electVdWAB","electdouble") ! diblock copolymer lseg determined in cadenas_sequence      
            ! assume A is present and type_of_monomer =1 
            lsegAA(1) = lsegPAA
            lsegAA(2) = lsegPAMPS        
        case ("electA")                 ! homopolymer weak polyacid VdW     
            lsegAA = lsegPAA
        case ("neutral")                ! homopolymer neutral
            lsegAA = lsegPEG
        case ("brush_mul","brush","brush_neq","brushvarelec","brushborn","brushssdna")
            ! lsegAA = lsegPAA !0.36287_dp
        case default
            print*,"Error: in set_lsegAA, systype=",systype
            print*,"stopping program"
            stop
        end select  
    endif   

end subroutine set_lsegAA

     
function set_lseg_Homopolymer()result(lengthseg)

    use parameters, only : lsegAA

    real(dp) :: lengthseg
    
    lengthseg=lsegAA(1)
    
end function set_lseg_Homopolymer



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



!     .. returns number of A monomer
!     .. this works only for binary AB polymers

function number_Amonomers(nseg)result(npolA)
  
    use chains, only : isAmonomer
    implicit none

    integer , intent(in) :: nseg
    
    integer :: npolA

    !     .. local variables
    integer :: s

    npolA=0
    do s=1,nseg
        if(isAmonomer(s).eqv..true.) npolA=npolA+1
    enddo
  

end function number_Amonomers



!  .. assign type_of_monomer and type_of_monomer_char from values in file named filename

subroutine read_type_of_monomer(type_of_monomer, type_of_monomer_char,filename, nseg)

    use mpivars
    use  myutils
    implicit none 
    
    !     .. arguments 
    integer, intent(inout) :: type_of_monomer(:)
    character(len=2), intent(inout) :: type_of_monomer_char(:)
    character(lenfname), intent(in) :: filename
    integer,  intent(in) :: nseg

    !      .. local variables
    integer :: ios, un  ! un = unit number
    integer :: s
    character(80) :: istr,str,letter

    !     .. reading in of variables from file
    open(unit=newunit(un),file=filename,iostat=ios,status='old')
    if(ios/=0 ) then
        write(istr,'(I2)')ios
        str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif
    
    s=0
    ios=0

    do while (s<nseg.and.ios==0)
        s=s+1
        read(un,*,iostat=ios)type_of_monomer(s),type_of_monomer_char(s)
    enddo

    
    if(s/=nseg.or.ios/=0) then 
        str="reached end of file before all elements read or ios error"
        print*,str
        str="read file "//trim(adjustl(filename))//" failed"
        print*,str
        call MPI_FINALIZE(ierr)
        stop
    endif

    close(un)

end subroutine read_type_of_monomer


! ismonomer_of_type is a table which row index is the segment  and column index correspond to the segment type  
! pre: type_of_monomer needs to be initialized
! post: table list of logicals indication is segement s is of type t

subroutine make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)
    
    implicit none

    logical, intent(out) :: ismonomer_of_type(:,:)
    integer, intent(in) :: type_of_monomer(:) 
    integer, intent(in) :: nseg
    integer, intent(in) :: nsegtypes

    ! local variable 
    integer :: s,t 

    do s=1,nseg 
        do t=1,nsegtypes   
            ismonomer_of_type(s,t)=.false.
        enddo
        !print*,"type_of_monomer(",s,")=",type_of_monomer(s)
        ismonomer_of_type(s,type_of_monomer(s))=.true.    
    enddo

end subroutine make_type_table




! routine determines is segment type t is chargeable
! pre: zpol needs to be initialized
! post: ismonomer_chargeable list of logicals

subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
    
    implicit none
 
    logical, intent(out) :: ismonomer_chargeable(:)
    integer, intent(in) :: zpol(:,:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t 

    do t=1,nsegtypes
        if(zpol(t,1)==0.and.zpol(t,2)==0) then 
            ismonomer_chargeable(t)=.false.
        else
            ismonomer_chargeable(t)=.true.
        endif
    enddo

end subroutine make_charge_table


logical function is_polymer_neutral(ismonomer_chargeable, nsegtypes)
    
    implicit none
 
    logical, intent(in) :: ismonomer_chargeable(:)
    integer, intent(in) :: nsegtypes
    
    ! local variable
    integer :: t,flag 

    flag=0
    do t=1,nsegtypes
        if(.not.ismonomer_chargeable(t)) flag=flag+1
    enddo

    if(flag==nsegtypes) then 
        is_polymer_neutral=.true. 
    else
        is_polymer_neutral=.false. 
    endif

end  function is_polymer_neutral


end module chaingenerator
