! --------------------------------------------------------------|
!                                                               | 
! chainsgenerator.f90:                                          |       
! generator chains in center of box                             |
! --------------------------------------------------------------|


module chaingenerator 

    use precision_definition   
    implicit none

    integer, parameter :: lenfname=40
    integer :: conf_write

    private :: lenfname, conf_write
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
    case ("FILE_lammps_xyz") 
        call read_chains_lammps_XYZ
    case ("FILE_lammps_trj") 
        call read_chains_lammps_trj  
    case default
        text="chainmethod not equal to MC, FILEE_lammps_XYZ, or FILE_lammps_trj"
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
    use parameters, only : geometry, lseg, write_mc_chains, isVdW, isVdWintEne
    use volume, only : nx, ny, nz, delta
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use volume, only : ut, vt
    use myutils
    use cadenas_linear
    use cadenas_sequence

    !     .. variable and constant declaractions      

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
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    character(len=lenText) :: text, istr
    integer  :: xi,yi,zi ,un_trj, un_ene
    real(dp), dimension(:), allocatable :: energy        
   

    !     .. executable statements
    !     .. initializations of variables     
       
    conf = 1                 ! counter for conformations
    seed = 435672*(rank+1)   ! seed for random number generator  different on each node
    maxnchains = 12
    maxntheta = 6            ! maximum number of rotation in xy-plane  
    theta_angle = 2.0_dp*pi/maxntheta
    Lz = nz*delta            ! maximum height box 
    Lx = nx*delta            ! maximum width box 
    Ly = ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp           ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp
   
    allocate(energy(maxnchains))
    energy=0.0_dp
            
    if(write_mc_chains) then 
        conf_write=0
        un_trj= open_chain_lammps_trj()
        un_ene= open_chain_energy()
    endif   

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

        if(isVdWintEne) energy= VdWpotentialenergy(chain,nchains)

        if(write_mc_chains)  then
            call write_chain_lammps_trj(un_trj,chain,nchains)    
            do j=1,nchains
                write(un_ene,*)energy(j)    
            enddo    
        endif
            
        if(geometry=="cubic") then
            
            do j=1,nchains   
            
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)-chain(1,1,j)
                    xp(s) = chain(2,s,j)-chain(2,1,j)
                    yp(s) = chain(3,s,j)-chain(3,1,j)
                enddo

                do ntheta=1,maxntheta                  ! rotation in xy-plane

                    theta = ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s) = xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s) =-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    

                    weightchain(conf)=.TRUE.

                    do s=1,nseg

                        ! .. translation onto center box 
                        x(s) = xpp(s) + xcm
                        y(s) = ypp(s) + ycm
                        z(s) = zp(s)  + zcm

                        ! .. check z coordinate 
                        if((z(s)>Lz)) weightchain(conf)=.FALSE.

                        ! .. periodic boundary conditions in x-direction and y-direction an z-direction 
                        x(s) = pbc(x(s),Lx)
                        y(s) = pbc(y(s),Ly)
                        z(s) = pbc(z(s),Lz)

                        ! .. transforming form real- to lattice coordinates                 
                        xi = int(x(s)/delta)+1
                        yi = int(y(s)/delta)+1
                        zi = int(z(s)/delta)+1
                        
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        indexchain_init(s,conf) = idx
                        if(idx<=0) then
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        endif
                    enddo            
            
                    energychain(conf)=energy(j) 
                    conf = conf +1

                enddo         ! end loop over rotations
            
            enddo             ! end loop over nchains

        else if(geometry=="prism") then
            
            do j=1,nchains   
            
                do s=1,nseg                     !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)-chain(1,1,j)
                    xp(s) = chain(2,s,j)-chain(2,1,j)
                    yp(s) = chain(3,s,j)-chain(3,1,j)
                enddo

                do ntheta=1,maxntheta                  ! rotation in xy-plane

                    theta = ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s)= xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s)=-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    

                    do s=1,nseg

                        xpp(s)=xpp(s) + xcm
                        ypp(s)=ypp(s) + ycm
                        z(s)  = zp(s) + zcm
                            
                        ! .. check z coordinate 
                        if((z(s)>Lz)) weightchain(conf)=.FALSE. 

                        ! .. translation onto correct grafting area translation in xy plane 
                        x(s) = ut(xpp(s),ypp(s))
                        y(s) = vt(xpp(s),ypp(s))
                        
                        ! .. periodic boundary conditions in x-direction and y-direction  
                        x(s) = pbc(x(s),Lx)
                        y(s) = pbc(y(s),Ly)
                        z(s) = pbc(z(s),Ly)

                        ! .. transforming form real- to lattice coordinates                 
                        xi=int(x(s)/delta)+1
                        yi=int(y(s)/delta)+1
                        zi=int(z(s)/delta)+1

                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        indexchain_init(s,conf) = idx
                        if(idx<=0) then
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        endif
                        
                    enddo         ! end loop over graft points

                    energychain(conf)=energy(j)
                    conf = conf +1  ! end loop over rotations

                enddo     
            
            enddo                 ! end loop over nchains

        else 
           print*,"Error: in make_chains_mc: geometry not cubic or prism"
           print*,"stopping program"
           stop
        endif
     
    enddo                     ! end while loop
    
    !  .. end chains generation 
      
    write(istr,'(I4)')rank
    text='AB Chains generated on node '//istr
    call print_to_log(LogUnit,text)

    ! .. find lowest global chainenergy and substract from all 
    call global_minimum_chainenergy()   

    allowedconf=0
    do k=1,cuantas
        if(weightchain(k).eqv..TRUE.) allowedconf = allowedconf+1
    enddo
      
    if(allowedconf/=cuantas) then 
        print*,"allowedconf=",allowedconf,"rank=",rank
        print*,"make_chains_mc: Not all conformations met constraint z(s)<Lz)"    
    endif
    
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)

    if(write_mc_chains) then
        close(un_trj)
        close(un_ene)
    endif   

    deallocate(energy)  

end subroutine make_chains_mc



! Reads confomations from a file called traj.xyz
! Format lammps trajectory file in xyz format
! number of ATOMS much equal nseg+1  
! differeent processor rank get assigned different conformations

subroutine read_chains_lammps_XYZ(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank                                                                                     
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
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle     
    integer  :: xi,yi,zi
    real(dp) :: xc,yc,zc          
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

    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.xyz'
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
    xcm= Lx/2.0_dp           ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    do while ((conf<cuantas).and.(ios==0))
    
        if(conf.ne.0) then ! skip lines
            read(un,*,iostat=ios)str
            read(un,*,iostat=ios)str
        else               ! read preamble
            read(un,*,iostat=ios)nsegfile
            read(un,*,iostat=ios)
            if(nsegfile.ne.nseg) then ! nsegfile =nseg+1 . The chain has one 'segement' more: tethered point.
                print*,"nseg chain file not equal internal nseg : stop program"
                stop
            endif    
        endif
        
        do s=1,nseg              ! .. read form  trajecotory file
            read(un,*,iostat=ios)idatom,xc,yc,zc
            xseg(1,s) = xc/10.0_dp  ! .. convert from Angstrom to nm
            xseg(2,s) = yc/10.0_dp
            xseg(3,s) = zc/10.0_dp 
        enddo

        if(ios==0) then ! read was succesfull 

             conffile=conffile +1 

            ! .. 'first' segment, sets origin 
            do s=1,nseg        
                chain(2,s) = (xseg(1,s)-xseg(1,1)) 
                chain(3,s) = (xseg(2,s)-xseg(2,1))
                chain(1,s) = (xseg(3,s)-xseg(3,1))
            enddo

            ! .. compute chain energy 
            energy= VdWpotentialenergy_lammps(chain)
        
            do rot=1,rotmax        
                 
                nchain=1
                is_positive_rot=.true.

                do while (nchain.lt.maxattempts)
                    is_positive_rot=rotation(chain,chain_rot,nseg-1)
                    nchain=nchain+1
                enddo
                !print*,"conf=",conf,"is_positve=",is_positive_rot

                if(is_positive_rot) then
            
                    if(geometry=="cubic") then !  transforming form real- to lattice coordinates
            
                        do ntheta=1,maxntheta                  ! rotation in xy-plane

                            conf=conf+1
                            theta= ntheta * theta_angle

                            do s=1,nseg

                                xp(s)= chain_rot(2,s)*cos(theta)+chain_rot(3,s)*sin(theta) 
                                yp(s)=-chain_rot(2,s)*sin(theta)+chain_rot(3,s)*cos(theta)
                                zp(s)= chain_rot(1,s)  

                            enddo 

                            do s=1,nseg

                                ! .. translation onto correct volume translation to cneter xyz 
                                x(s) = xp(s) + xcm
                                y(s) = yp(s) + ycm
                                z(s) = zp(s) + zcm

                                ! .. periodic boundary conditions in x-direction and y-direction and z-direction 
                                x(s) = pbc(x(s),Lx)
                                y(s) = pbc(y(s),Ly)
                                z(s) = pbc(z(s),Lz)

                                ! .. transforming form real- to lattice coordinates                 
                                xi = int(x(s)/delta)+1
                                yi = int(y(s)/delta)+1
                                zi = int(z(s)/delta)+1
                                
                                call linearIndexFromCoordinate(xi,yi,zi,idx) 

                                indexchain_init(s,conf) = idx
                                if(idx<=0) then
                                    print*,"index=",idx, " xc=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                endif
                            
                            enddo
                                
                            energychain(conf)=energy
                            
                        enddo         ! end loop over rotations
            
                    else if(geometry=="prism") then
                     ! to be done later 

                    else 
                        print*,"Error: in make_chains_lammps_xyz: geometry not cubic or prism"
                        print*,"stopping program"
                        stop
                    endif
                endif  ! if 
            enddo      ! rotation loop
        endif     
    enddo          ! end while loop                                                                                                          
    
    !  .. end chains generation     

    if(conf<=cuantas) then
        print*,"subroutine make_chains_lammps_XYZ :"
        print*,"conf     = ",conf," less then cuantas     = ",cuantas
        print*,"conffile = ",conffile
        cuantas = conf
    else
        text="Chains generated: subroutine make_chains_lammps_file"
        call print_to_log(LogUnit,text)
        readinchains=conffile
    endif

    ! .. find lowest global chainenergy and substract from all 
    call global_minimum_chainenergy()   

    close(un)

end subroutine read_chains_lammps_XYZ


! Reads confomations from a file called traj.xyz
! Format lammps trajectory file in xyz format
! number of ATOMS much equal nseg+1  
! differeent processor rank get assigned different conformations

subroutine read_chains_lammps_trj(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank                                                                                    
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
    integer :: ix,iy,iz,idxtmp,ntheta
    integer :: nchains              ! number of rotations
    integer :: maxnchains           ! number of rotations
    integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg+1)     ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1  
    real(dp) :: chain_rot(3,nseg+1) 
    real(dp) :: xseg(3,nseg+1)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(dp) :: theta, theta_angle
    real(dp) :: xc,yc,zc          
    integer  :: xi,yi,zi         
    real(dp) :: energy                                               
    character(len=25) :: fname
    integer :: ios 
    character(len=30) :: str
    real(dp) :: xbox0,xbox1,scalefactor
    integer :: nchain, rotmax, maxattempts, idatom, item, moltype
    integer :: un ! unit number
    logical :: is_positive_rot, exist
    character(len=lenText) :: text,istr

    ! .. executable statements                                                                                                     
    ! .. open file                                                                                              
    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'chains.lammostrj file does not exit'
        if (present(info)) info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
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

    maxntheta =12                ! maximum number of rotation in xy-plane  
    theta_angle= 2.0_dp*pi/maxntheta    
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= Lz/2.0_dp

    do while ((conf<cuantas).and.(ios==0))
    
        if(conf.ne.0) then ! skip preamble 
            do i=1,9
                read(un,*,iostat=ios)str
            enddo    
        else               ! read preamble
            do i=1,3
                read(un,*,iostat=ios)str
            enddo   
            read(un,*,iostat=ios)nsegfile
            read(un,*,iostat=ios)str
            read(un,*,iostat=ios)xbox0,xbox1
            read(un,*,iostat=ios)str
            read(un,*,iostat=ios)str 
            read(un,*,iostat=ios)str 
            if(nsegfile.ne.nseg) then ! nsegfile =nseg+1 . The chain has one 'segement' more: tethered point.
                print*,"nseg chain file not equal internal  nseg : stop program"
                stop
            endif    
            scalefactor=(xbox1-xbox0)*unit_conv
        endif

    
        do s=1,nseg              ! .. read form  trajecotory file
            read(un,*,iostat=ios)item,idatom,moltype,xc,yc,zc,ix,iy,iz
            xseg(1,item) = xc*scalefactor 
            xseg(2,item) = yc*scalefactor
            xseg(3,item) = zc*scalefactor
        enddo
     
       
        if(ios==0) then ! read was succesfull 

            conffile=conffile +1 

            ! .. 'first' segment, sets origin 

            do s=1,nseg        
                chain(2,s) = (xseg(1,s)-xseg(1,1)) 
                chain(3,s) = (xseg(2,s)-xseg(2,1))
                chain(1,s) = (xseg(3,s)-xseg(3,1))
            enddo

            ! .. compute chain energy 
            energy= VdWpotentialenergy_lammps(chain)

            do rot=1,rotmax        
                 
                is_positive_rot=rotation(chain,chain_rot,nseg-1)
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    zp(s) = chain_rot(1,s) !1->z
                    xp(s) = chain_rot(2,s) !2->x 
                    yp(s) = chain_rot(3,s) !3->y
                enddo
                

                if(geometry=="cubic") then !  transforming form real- to lattice coordinates
        
                    do ntheta=1,maxntheta                  ! rotation in xy-plane

                        conf=conf+1
                        theta = ntheta * theta_angle          
                
                        do s=1,nseg
                            xpp(s)= xp(s)*cos(theta)+yp(s)*sin(theta) 
                            ypp(s)=-xp(s)*sin(theta)+yp(s)*cos(theta)   
                        enddo    

                        do s=1,nseg

                        ! .. translation onto correct volume translation to cneter xyz 
                            x(s) = xpp(s) + xcm
                            y(s) = ypp(s) + ycm
                            z(s) =  zp(s) + zcm

                            ! .. periodic boundary conditions in x-direction and y-direction an z-direction 
                            x(s) = pbc(x(s),Lx)
                            y(s) = pbc(y(s),Ly)
                            z(s) = pbc(z(s),Lz)

                            ! .. transforming form real- to lattice coordinates                 
                            xi = int(x(s)/delta)+1
                            yi = int(y(s)/delta)+1
                            zi = int(z(s)/delta)+1
                            
                            call linearIndexFromCoordinate(xi,yi,zi,idx)
                            
                            indexchain_init(s,conf) = idx

                            if(idx<=0.or.idx>nsize) then
                                print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            endif
                        
                        enddo
                            
                        energychain(conf)=energy
                        
                    enddo         ! end loop over rotations
        
                else if(geometry=="prism") then
                 
                ! to be done later 

                else 
                    print*,"Error: in make_chains_lajmmps_trj: geometry not cubic or prism"
                    print*,"stopping program"
                    stop
                endif
         
            enddo      ! rotation loop
        endif     
    enddo          ! end while loop                                                                                                          
    !  .. end chains generation    
    

    if(conf<=cuantas) then
        print*,"subroutine make_chains_lammps_trj :" 
        print*,"conf     = ",conf," less then imposed cuantas     = ",cuantas
        print*,"conffile = ",conffile
        cuantas=conf           
    else
        text="Chains generated: subroutine make_chains_lammps_trj"
        call print_to_log(LogUnit,text)
        readinchains=conffile
    endif

    ! .. find lowest global chainenergy and substract from all 
    call global_minimum_chainenergy()   


    close(un)

end subroutine read_chains_lammps_trj


! post: isAmonomer set 
! for chaintype==multi , type_of_monomer,type_of_monomer_char
! ismonomer_type are set.

subroutine make_sequence_chain(freq,chaintype)
  
    use globals, only : nseg, nsegtypes
    use chains, only :  isAmonomer,type_of_monomer_char, type_of_monomer,ismonomer_of_type
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
    
    use  globals, only : nseg, cuantas
    use  chains, only : indexchain,indexchain_init,max_confor
    use  volume, only : nz, coordinateFromLinearIndex,linearIndexFromCoordinate

    integer :: conf,s 
    integer :: indx  ! temporary index of chain
    integer :: ix,iy,iz,izpbc

    do conf=1,cuantas        
        do s=1,nseg
            indx=indexchain_init(s,conf)
            call coordinateFromLinearIndex(indx,ix,iy,iz)
            if(iz<=nz) then  ! nz is distance in between plates  
                indexchain(s,conf)=indx
            else
                izpbc=ipbc(iz,nz)
                print*,"iz=",iz,"izpbc=",izpbc
                call linearIndexFromCoordinate(ix,iy,izpbc,indx)
                indexchain(s,conf)=indx
            endif
        enddo
    enddo

end subroutine  chain_filter

function minimum_chainenergy() result(min_chainenergy)

    use  mpivars
    use  globals, only :cuantas
    use  chains, only : energychain

    real(dp) :: min_chainenergy

    integer :: conf, unitno
    
    unitno=10*rank+100
    min_chainenergy =0.0_dp
    do conf=1,cuantas        
        if(min_chainenergy > energychain(conf)) min_chainenergy=energychain(conf)
    enddo
   
end function


subroutine global_minimum_chainenergy()

    use  mpivars
    use  globals, only :cuantas
    use  chains, only : energychain

    real(dp) :: localmin(2), globalmin(2)
    integer :: i

    localmin(1)=minimum_chainenergy()
    localmin(2)=rank   
    print*,"rank ",rank ," has local lowest value of", localmin(1)
    
    call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
    call MPI_ALLREDUCE(localmin, globalmin, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
    
    if (rank == 0) then  
        print*,"Rank ",globalmin(2)," has lowest value of", globalmin(1)  
    endif

    call shift_energy_chain(globalmin(1))

end subroutine



! post: substract minimum internal energy from all internal chain energies
subroutine shift_energy_chain(min_chainenergy)
    
    use  globals, only :cuantas
    use  chains, only : energychain
    
    real(dp), intent(in) :: min_chainenergy

    ! local
    integer :: conf
    do conf=1,cuantas        
        energychain(conf)=energychain(conf)-min_chainenergy
    enddo

end subroutine shift_energy_chain


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
        case ("neutral","neutralnoVdW")                ! homopolymer neutral
            !lsegAA = lsegPEG
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


function ipbc(ival,imax) result(intpbc)
    implicit none 
    integer, intent(in) :: ival
    integer, intent(in) :: imax
    integer :: intpbc

    if(ival>0) then
        intpbc=ival-int((ival-1)/imax)*imax
    else
        intpbc=ival-(int((ival-1)/imax)-1)*imax
    endif

end function


!     .. returns number of A monomer
!     .. this works only for binary AB polymers

function number_Amonomers(nseg)result(npolA)
  
    use chains, only : isAmonomer

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
        ismonomer_of_type(s,type_of_monomer(s))=.true.    
    enddo

end subroutine make_type_table


! routine determines is segment type t is chargeable
! pre: zpol needs to be initialized
! post: ismonomer_chargeable list of logicals

subroutine make_charge_table(ismonomer_chargeable,zpol,nsegtypes)
 
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


! make a lammps trajectory file based on  indexchain_init

subroutine write_indexchain_lammps_trj(info)

    use mpivars, only : rank
    use globals, only : nseg, cuantas
    use myutils, only : newunit, lenText
    use volume, only : delta, nz, coordinateFromLinearIndex
    use chains, only : indexchain_init, type_of_monomer
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    character(len=lenText) :: text, istr
    character(len=25) :: fname
    integer :: ios, un_trj 
    real(dp):: xs, ys, zs
    integer :: ix, iy, iz, i, j, k, idx
    real(dp) :: xbox0, xbox1
    integer :: idatom, item, moltype, conf

    ! open file 

    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

    ix=0
    iy=0
    iz=0
    xbox0=0.0_dp
    xbox1=nz*delta
    
    do conf=1,cuantas
        ! write preamble 
        write(un_trj,*)'ITEM: TIMESTEP' 
        write(un_trj,*)conf 
        write(un_trj,*)'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg 
        write(un_trj,*)'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)'ITEM: ATOMS id mol type xs ys zs ix iy iz'
        !  determine xs, ys, zs,
        do item=1,nseg
            idx=indexchain_init(item,conf)
            call coordinateFromLinearIndex(idx, i, j, k)
            xs=(i-0.5_dp)*delta/xbox1
            ys=(j-0.5_dp)*delta/xbox1
            zs=(k-0.5_dp)*delta/xbox1
            idatom=1
            moltype=type_of_monomer(item) 
            write(un_trj,*)item,idatom,moltype,xs,ys,zs,ix,iy,iz
        enddo    
    enddo
        
    close(un_trj)    


end subroutine write_indexchain_lammps_trj


function open_chain_lammps_trj(info)result(un_trj)

    use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    integer :: un_trj
    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios
    logical :: exist
    
    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    inquire(file=fname,exist=exist)
    if(.not.exist) then
        open(unit=newunit(un_trj),file=fname,status='new',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    else
        open(unit=newunit(un_trj),file=fname,status='old',iostat=ios)
        if(ios >0 ) then
            print*, 'Error opening : ',fname,' file : iostat =', ios
            if (present(info)) info = myio_err_chainsfile
            return
        endif
    endif    
end function


subroutine write_chain_lammps_trj(un_trj,chain,nchains)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nz, coordinateFromLinearIndex
    use chains, only :  type_of_monomer

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains
    integer , intent(in) :: un_trj
    
    character(len=lenText) :: istr
    real(dp) :: xs, ys, zs
    integer :: ix, iy, iz, i, j, k
    real(dp) :: xbox0, xbox1
    integer :: idatom, item, moltype, conf

    ix=0
    iy=0
    iz=0
    xbox0=0.0_dp
    xbox1=nz*delta
    do j=1,nchains   
        
        conf_write=conf_write+1 
        ! write preamble 
        write(un_trj,*)'ITEM: TIMESTEP' 
        write(un_trj,*)conf_write 
        write(un_trj,*)'ITEM: NUMBER OF ATOMS' 
        write(un_trj,*)nseg 
        write(un_trj,*)'ITEM: BOX BOUNDS ff ff ff'
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)xbox0,xbox1
        write(un_trj,*)'ITEM: ATOMS id mol type xs ys zs ix iy iz'
        do item=1,nseg
            zs = chain(1,item,j)/xbox1
            xs = chain(2,item,j)/xbox1
            ys = chain(3,item,j)/xbox1 
            idatom=1
            moltype=type_of_monomer(item) 
            write(un_trj,*)item,idatom,moltype,xs,ys,zs,ix,iy,iz
        enddo    
    enddo
        
end subroutine write_chain_lammps_trj


function open_chain_energy(info)result(un_ene)

    use mpivars, only : rank
    use myutils, only : newunit, lenText
    use myio, only : myio_err_chainsfile

    integer, optional, intent(inout) :: info

    integer :: un_ene
    ! local
    character(len=lenText) :: istr
    character(len=25) :: fname
    integer :: ios

    write(istr,'(I4)')rank
    fname='traj.'//trim(adjustl(istr))//'.ene'
    open(unit=newunit(un_ene),file=fname,status='replace',iostat=ios)
    if(ios >0 ) then
        print*, 'Error opening : ',fname,' file : iostat =', ios
        if (present(info)) info = myio_err_chainsfile
        return
    endif

end function

function VdWpotentialenergy(chain,nchains)result(Energy)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in) :: chain(:,:,:)
    integer, intent(in) :: nchains

    real(dp) :: Energy(nchains)
    real(dp) :: Ene,sqrlseg,sqrdist
    integer :: k,i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
    real(dp) :: Lz,Ly,Lx

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 

    do k=1,nchains
        Ene=0.0_dp      
        do i=1,nseg
            do j=i+1,nseg

                s=type_of_monomer(i)
                t=type_of_monomer(j)
               
                zi = pbc(chain(1,i,k),Lz)
                xi = pbc(chain(2,i,k),Lx)
                yi = pbc(chain(3,i,k),Ly) 
                zj = pbc(chain(1,j,k),Lz)
                xj = pbc(chain(2,j,k),Lx)
                yj = pbc(chain(3,j,k),Ly) 

                sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
                sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

                Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3

            enddo
        enddo 
        Energy(k) = Ene            
    enddo

 end function VdWpotentialenergy
 
function VdWpotentialenergy_lammps(chain)result(Energy)

    use globals, only : nseg
    use myutils, only : newunit, lenText
    use volume, only : delta, nx, ny, nz, coordinateFromLinearIndex
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in) :: chain(:,:)
   
    real(dp) :: Energy
    real(dp) :: Ene,sqrlseg,sqrdist
    integer :: k,i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
    real(dp) :: Lz,Ly,Lx

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    
    Ene=0.0_dp      

    do i=1,nseg
        do j=i+1,nseg

            s=type_of_monomer(i)
            t=type_of_monomer(j)

            zi = pbc(chain(1,i),Lz)
            xi = pbc(chain(2,i),Lx)
            yi = pbc(chain(3,i),Ly) 
            zj = pbc(chain(1,j),Lz)
            xj = pbc(chain(2,j),Lx)
            yj = pbc(chain(3,j),Ly) 

            sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
            sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

            Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3

        enddo
    enddo 
    Energy = Ene            

end function VdWpotentialenergy_lammps
   
        
end module chaingenerator
