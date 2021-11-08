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
    real(dp) :: xgraftloop(3,2), xgraftlinear(3)

    private :: lenfname, conf_write
    private :: pbc 
    private :: xgraftloop, xgraftlinear
   
contains


subroutine make_chains(chainmethod)

    use mpivars, only : ierr
    use myutils, only : print_to_log, LogUnit, lenText
    use myio, only : myio_err_chainmethod

    character(len=15), intent(in)  :: chainmethod

    integer :: i, info
    character(len=lenText) :: text, istr

    info=0

    select case (chainmethod)
    case ("MC")
        call make_chains_mc()
    case ("FILE_lammps_xyz") 
        call read_chains_lammps_XYZ(info)
    case ("FILE_lammps_trj") 
        call read_chains_lammps_trj(info)
    case ("FILE_XYZ")
        call read_chains_xyz(info)  
    case default
        text="chainmethod not equal to MC, FILE_lammps_XYZ, FILE_lammps_trj or FILKE_XYZ"
        call print_to_log(LogUnit,text)
        print*,text
        info=myio_err_chainmethod
    end select

    if(info/=0) then
        write(istr,'(I3)')info
        text="Error make_chains: chain generation failed: info = "//istr//" : end program."
        call print_to_log(LogUnit,text)
        print*,text
        call MPI_FINALIZE(ierr)
        stop
    endif    

end subroutine make_chains


!  init of cuantas polymer configurations of polymer chain anchored onto a flat surface 

subroutine make_chains_mc()
  
    use mpivars
    use globals
    use chains
    use random
    use parameters, only : geometry, lseg, write_mc_chains, isVdW, isVdWintEne
    use parameters, only : maxnchainsrotations, maxnchainsrotationsxy
    use volume, only : nx, ny, nz, delta
    use volume, only : coordinateFromLinearIndex, linearIndexFromCoordinate
    use volume, only : ut, vt
    use volume, only : position_graft, nset_per_graft
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
    real(dp) :: chain_rot(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg) ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates
    real(dp) :: xpp(nseg), ypp(nseg), zpp(nseg)  
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: xpt,ypt          ! coordinates
    real(dp) :: theta, theta_angle
    character(len=lenText) :: text, istr
    integer  :: xi,yi,zi ,un_trj, un_ene, segcenter
    real(dp) :: energy   
    logical :: saw     
   
    !     .. executable statements
    !     .. initializations of variables     
       
    conf = 1                 ! counter for conformations
    seed = 435672*(rank+1)   ! seed for random number generator  different on each node
    maxnchains = maxnchainsrotations
    maxntheta = maxnchainsrotationsxy         ! maximum number of rotation in xy-plane  
    theta_angle = 2.0_dp*pi/maxntheta
    Lz = nz*delta            ! maximum height box 
    Lx = nx*delta            ! maximum width box 
    Ly = ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp           ! center box
    ycm= Ly/2.0_dp
    zcm= 0.0_dp
   
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

        if(write_mc_chains) then 
            energy=0.0_dp
            write(un_ene,*)energy  
            call write_chain_lammps_trj(un_trj,chain,nchains)  
        endif    
   
            
        if(geometry=="cubic") then

            g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

            xpt =  position_graft(g,1)    ! position of graft point
            ypt =  position_graft(g,2)  

            do j=1,nchains   
            
                do s=1,nseg                          !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j)
                enddo

                do ntheta=1,maxntheta                  ! rotation in xy-plane

                    theta = ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s) = xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s) =-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    

                    do s=1,nseg

                        ! .. translation onto center box 
                        x(s) = xpp(s) + xpt
                        y(s) = ypp(s) + ypt
                        z(s) = zp(s)  

                        ! .. periodic boundary conditions in x-direction and y-direction 
                        chain_rot(2,s) = pbc(x(s),Lx)
                        chain_rot(3,s) = pbc(y(s),Ly)
                        chain_rot(1,s) = z(s)     !  no pbc in z-direction    

                        ! .. transforming form real- to lattice coordinates                 
                        xi = int(chain_rot(2,s)/delta)+1
                        yi = int(chain_rot(3,s)/delta)+1
                        zi = int(chain_rot(1,s)/delta)+1
                        
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        indexchain_init(s,conf) = idx
                        if(idx<=0) then
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        endif
                    enddo            
            
                    conf = conf +1 

                enddo         ! end loop over rotations
            
            enddo             ! end loop over nchains

        else if(geometry=="prism") then
            
            do j=1,nchains   
            
                do s=1,nseg                     !  transforming form real- to lattice coordinates
                    zp(s) = chain(1,s,j)
                    xp(s) = chain(2,s,j)
                    yp(s) = chain(3,s,j)
                enddo

                do ntheta=1,maxntheta            ! rotation in xy-plane

                    theta = ntheta * theta_angle          
                    do s=1,nseg
                        xpp(s)= xp(s)*cos(theta)+yp(s)*sin(theta) 
                        ypp(s)=-xp(s)*sin(theta)+yp(s)*cos(theta)   
                    enddo    

                    do s=1,nseg

                        xpp(s)=xpp(s) + xcm
                        ypp(s)=ypp(s) + ycm
                        z(s)  = zp(s) + zcm
                            
                        ! .. translation onto correct grafting area translation in xy plane 
                        x(s) = ut(xpp(s),ypp(s))
                        y(s) = vt(xpp(s),ypp(s))
                        
                        ! .. periodic boundary conditions in x-direction and y-direction an z-direction 
                        chain_rot(2,s) = pbc(x(s),Lx)
                        chain_rot(3,s) = pbc(y(s),Ly)
                        chain_rot(1,s) = z(s)        ! .. no pbc in z-direction    

                        ! .. transforming form real- to lattice coordinates                 
                        xi = int(chain_rot(2,s)/delta)+1
                        yi = int(chain_rot(3,s)/delta)+1
                        zi = int(chain_rot(1,s)/delta)+1

                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                        indexchain_init(s,conf) = idx
                        if(idx<=0) then
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                        endif
                        
                    enddo         ! end loop over graft points

                    conf = conf +1 
                
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

  
    if(isHomopolymer.eqv..FALSE.) deallocate(lsegseq)

    if(write_mc_chains) then
        close(un_trj)
        close(un_ene)
    endif   

    energychain_init=0.0_dp ! no internal energy 

end subroutine make_chains_mc


subroutine read_chains_XYZ(info)

    use parameters, only : chaintopol
    
    ! .. argument
    integer, intent(out) :: info

    
    if(chaintopol=="loop") then 
        call read_chains_XYZ_loop(info)
    else if(chaintopol=="linear") then
        call read_chains_XYZ_linear(info)
    else
        print*,"Error: in read_chains_XYZ: chaintop not loop or linear"
        print*,"stopping program"
        stop    
    endif

end subroutine


! Reads confomations from a file called traj.xyz
! Format lammps trajectory file in xyz format
! number of ATOMS much equal nseg+1  
! differeent processor rank get assigned different conformations

subroutine read_chains_lammps_XYZ(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank, numproc                                                                                     
    use globals
    use chains
    use random
    use parameters
    use volume, only : sgraft, nset_per_graft,position_graft
    use chain_rotation
    use myio, only : myio_err_chainsfile, myio_err_energyfile,myio_err_index,myio_err_conf
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,rot,g         ! dummy indices
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
    real(dp) :: xpp(nseg), ypp(nseg)
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm         ! sizes box
    real(dp) :: theta, theta_angle     
    integer  :: xi,yi,zi
    real(dp) :: xc,yc,zc          
    real(dp) :: energy                                               
    character(len=8) :: fname
    integer :: ios ,rankfile
    character(len=30) :: str
    integer :: nchain, rotmax, maxattempts, idatom
    integer :: un    ! unit number
    logical :: is_positive_rot, exist ,saw
    character(len=lenText) :: text,istr 


    ! .. executable statements  

    info=0

     ! get location graft points loop
    call read_graftpts_lammps_trj(info)
    if(info/=0) return

    ! .. open file 
    rankfile=mod(rank,nset_per_graft)                                                                                                 
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj.xyz file does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file                                                                    
    seed=435672               ! seed for random number generator                                                                             
    maxattempts=72            ! maximum number of attempts to rotate conf   
    maxnchains= maxnchainsrotations  ! maximum number of rotation conf                                                                  
    rotmax=maxnchains
    maxntheta = maxnchainsrotationsxy! maximum number of rotation in xy-plane  
    theta_angle= 2.0_dp*pi/maxntheta ! angle of rotation 

    ios=0

    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center box
    ycm= Ly/2.0_dp
    zcm= 0.0_dp

    do while ((conf<max_confor).and.(ios==0))
    
        if(conf.ne.1) then ! skip lines
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
                chain(1,s) = xseg(1,s)-xgraftloop(1,1) 
                chain(2,s) = xseg(2,s)-xgraftloop(2,1) 
                chain(3,s) = xseg(3,s)-xgraftloop(3,1) 
            enddo

            do rot=1,rotmax        
                 
                nchain=1
        
                is_positive_rot=.false.

                do while (nchain.lt.maxattempts.and..not.is_positive_rot)
                    is_positive_rot=rotationXaxis(chain,chain_rot,nseg-1)
                    nchain=nchain+1
                enddo

                if(is_positive_rot) then
            
                    if(geometry=="cubic") then !  transforming form real- to lattice coordinates
            
                       do ntheta=1,maxntheta      ! .. rotation in xy-plane and translation to center of xy-plane

                            theta= ntheta * theta_angle

                            do s=1,nseg
                                xp(s)=  chain_rot(1,s)*cos(theta)+chain_rot(2,s)*sin(theta) + xcm
                                yp(s)= -chain_rot(1,s)*sin(theta)+chain_rot(2,s)*cos(theta) + ycm
                                zp(s)=  chain_rot(3,s)   !+ zcm =  0                                 
                            enddo 
                       
                            do s=1,nseg

                                x(s) = pbc(xp(s),Lx) ! .. periodic boundary conditions in x and y direction
                                y(s) = pbc(yp(s),Ly)
                                z(s) = zp(s)       ! not in z- direction 

                                ! .. transforming form real- to lattice coordinates                 
                                xi = int(x(s)/delta)+1
                                yi = int(y(s)/delta)+1
                                zi = int(z(s)/delta)+1
                                    
                                call linearIndexFromCoordinate(xi,yi,zi,idx)
                                    
                                indexchain_init(s,conf) = idx

                                if(idx<=0.or.idx>nsize) then    
                                    text="Conformation outside box:"
                                    call print_to_log(LogUnit,text)  
                                    print*,text                          
                                    print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                    info= myio_err_index
                                    return
                                endif
                                
                            enddo
                            
                            call VdWpotentialenergy(chain_rot,energy)  
                            energychain_init(conf)=energy

                            conf=conf+1   
                        
                        enddo ! .. rotation  
            
                    else if(geometry=="prism") then

                        do ntheta=1,maxntheta      ! .. rotation in xy-plane and translation to center of xy-plane

                            theta = ntheta * theta_angle

                            do s=1,nseg
                                xp(s)=  chain_rot(1,s)*cos(theta)+chain_rot(2,s)*sin(theta) + xcm
                                yp(s)= -chain_rot(1,s)*sin(theta)+chain_rot(2,s)*cos(theta) + ycm
                                zp(s)=  chain_rot(3,s)   !+ zcm =  0                                 
                            enddo 
                       
                            do s=1,nseg

                                ! .. transformation to prism coordinates 
                                xpp(s) = ut(xp(s),yp(s))
                                ypp(s) = vt(xp(s),yp(s))

                                x(s) = pbc(xpp(s),Lx) ! .. periodic boundary conditions in x and y direction
                                y(s) = pbc(ypp(s),Ly)
                                z(s) = zp(s)          ! .. no pbc  in z- direction 

                                ! .. transforming form real- to lattice coordinates                 
                                xi = int(x(s)/delta)+1
                                yi = int(y(s)/delta)+1
                                zi = int(z(s)/delta)+1
                                    
                                call linearIndexFromCoordinate(xi,yi,zi,idx)
                                    
                                indexchain_init(s,conf) = idx

                                if(idx<=0.or.idx>nsize) then    
                                    text="Conformation outside box:"
                                    call print_to_log(LogUnit,text)  
                                    print*,text                          
                                    print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                    info= myio_err_index
                                    return
                                endif
                                
                            enddo
                            
                            call VdWpotentialenergy(chain_rot,energy)  
                            energychain_init(conf)=energy

                            conf=conf+1   
                        
                        enddo ! .. rotation  

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

    if(conf<=max_confor) then
        print*,"subroutine make_chains_lammps_XYZ :"
        print*,"conf     = ",conf," less then max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas = conf
    else
        text="Chains generated: subroutine make_chains_lammps_file"
        call print_to_log(LogUnit,text)
        readinchains=conffile
    endif

    if(.not.(isVdWintEne))energychain_init=0.0_dp


    close(un)

end subroutine read_chains_lammps_XYZ


subroutine read_graftpts_lammps_trj(info)

    use mpivars, only : rank, numproc
    use parameters, only : unit_conv
    use myio, only : myio_err_chainsfile, myio_err_graft
    use myutils,  only : newunit
    use volume, only : sgraft,nset_per_graft,position_graft

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    character(len=25) :: fname
    integer :: ios, rankfile
    real(dp) :: xc,yc,zc          
    integer  :: ix,iy,iz,un,s, i,t
    integer  :: item,moltype,nsegfile,idatom
    character(len=30) :: istr,str
    real(dp) :: xbox0,xbox1,scalefactor
    logical :: exist, isGraftItem

    ! .. executable statements 
    info=0

    !. . open file     
    rankfile=mod(rank,nset_per_graft)                                                                                               
    write(istr,'(I4)')rankfile
    fname='traj-graft.'//trim(adjustl(istr))//'.lammpstrj'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj-graft.rank.lammpstrj file does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    ! read preamble/header
    do i=1,3
        read(un,*,iostat=ios)str
    enddo   
    read(un,*,iostat=ios)nsegfile
    read(un,*,iostat=ios)str
    read(un,*,iostat=ios)xbox0,xbox1
    read(un,*,iostat=ios)str
    read(un,*,iostat=ios)str 
    read(un,*,iostat=ios)str 
   
    scalefactor=(xbox1-xbox0)*unit_conv
    isGraftItem=.false.
    do s=1,2
        read(un,*,iostat=ios)item,idatom,moltype,xc,yc,zc,ix,iy,iz
        if(item==sgraft) then
            t=1
            isGraftItem=.true.
        else
            t=2
        endif    
        xgraftloop(1,t)=xc*scalefactor
        xgraftloop(2,t)=yc*scalefactor
        xgraftloop(3,t)=zc*scalefactor    
    enddo
    
    close(un)

    if(.not.isGraftItem) then
        print*,'Error read_graftpts_lammps_trj: sgraft not in traj file.'
        info = myio_err_graft
    endif
        
end subroutine



! Reads confomations from a file called traj."rank".lammpstrj
! Format lammps trajectory file in xyz format
! number of ATOMS much equal nseg+1  
! differeent processor rank get assigned different conformations

subroutine read_chains_lammps_trj(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank, numproc                                                                                
    use globals
    use chains
    use random
    use parameters
    use volume,  only: position_graft, sgraft, nx, ny,nz, delta, nset_per_graft
    use chain_rotation, only : rotationXaxis
    use myio, only : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only : myio_err_conf, myio_err_nseg, myio_err_geometry
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit

    ! .. argument

    integer, intent(out) :: info

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
    real(dp) :: energy, energy2                                              
    character(len=25) :: fname
    integer :: ios, rankfile
    character(len=30) :: str
    real(dp) :: xbox0,xbox1,scalefactor
    integer :: nchain, rotmax, maxattempts, idatom, item, moltype, segcenter
    integer :: un,unw! unit number
    logical :: is_positive_rot, exist
    character(len=lenText) :: text,istr
    logical :: saw

    ! .. executable statements   

    info=0

    ! get location graft points loop
    call read_graftpts_lammps_trj(info)
    if(info/=0) return

    ! .. open file   
    !rankfile=int(rank*nset_per_graft/size) 
    rankfile=mod(rank,nset_per_graft)                                                                                         
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.lammpstrj'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj.rank.lammostrj file does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file                                                                    
    seed=435672               ! seed for random number generator                                                                             
    maxattempts=72            ! maximum number of attempts to rotate conf   
    maxnchains= maxnchainsrotations  ! maximum number of rotation conf                                                                  
    rotmax=maxnchains
    maxntheta = maxnchainsrotationsxy! maximum number of rotation in xy-plane  
    theta_angle= 2.0_dp*pi/maxntheta ! angle of rotation 

    ios=0
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center x-y plane
    ycm= Ly/2.0_dp
    zcm= 0.0_dp


    do while ((conf<max_confor).and.(ios==0))
    
        if(conf.ne.1) then ! skip preamble 
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
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            endif    
            scalefactor=(xbox1-xbox0)*unit_conv
        endif

    
        do s=1,nseg              ! .. read form  trajecotory file
            read(un,*,iostat=ios)item,idatom,moltype,xc,yc,zc,ix,iy,iz
            xseg(1,item) = xc*scalefactor 
            xseg(2,item) = yc*scalefactor  !
            xseg(3,item) = zc*scalefactor  ! permutated y and z 
        enddo
     
        if(ios==0) then ! read was succesfull 

            conffile=conffile +1 

            ! .. 'first' segment, sets origin 
           
            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xgraftloop(1,1) 
                chain(2,s) = xseg(2,s)-xgraftloop(2,1)  
                chain(3,s) = xseg(3,s)-xgraftloop(3,1) 
            enddo
  
            do rot=1,rotmax        
                 
                nchain=1
                is_positive_rot=.false.

                do while (nchain.lt.maxattempts.and..not.is_positive_rot)
                    is_positive_rot=rotationXaxis(chain,chain_rot,nseg-1)
                    nchain=nchain+1
                enddo

                if(is_positive_rot) then

                    
                    if(geometry=="cubic") then 

                        g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                        xpt =  position_graft(g,1)  ! position of graft point
                        ypt =  position_graft(g,2)  

                        do ntheta=1,maxntheta      ! .. rotation in xy-plane and translation to center of xy-plane

                            theta= ntheta * theta_angle

                            do s=1,nseg
                                xp(s)=  chain_rot(1,s)*cos(theta)+chain_rot(2,s)*sin(theta) + xpt
                                yp(s)= -chain_rot(1,s)*sin(theta)+chain_rot(2,s)*cos(theta) + ypt
                                zp(s)=  chain_rot(3,s)   !+ zcm =  0                                 
                            enddo 
                       
                            do s=1,nseg

                                x(s) = pbc(xp(s),Lx) ! .. periodic boundary conditions in x and y direction
                                y(s) = pbc(yp(s),Ly)
                                z(s) = zp(s)       ! not in z- direction 

                                ! .. transforming form real- to lattice coordinates                 
                                xi = int(x(s)/delta)+1
                                yi = int(y(s)/delta)+1
                                zi = int(z(s)/delta)+1
                                    
                                call linearIndexFromCoordinate(xi,yi,zi,idx)
                                    
                                indexchain_init(s,conf) = idx

                                if(idx<=0.or.idx>nsize) then    
                                    text="Conformation outside box:"
                                    call print_to_log(LogUnit,text)  
                                    print*,text                          
                                    print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                    info= myio_err_index
                                    return
                                endif
                                
                            enddo
                            
                            call VdWpotentialenergy(chain_rot,energy)  
                            energychain_init(conf)=energy

                            conf=conf+1   
                        
                        enddo ! .. rotation 
                            
                    else if(geometry=="prism") then
                    
                        g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                        xpt =  position_graft(g,1)  ! position of graft point
                        ypt =  position_graft(g,2)  


                        do ntheta=1,maxntheta      ! .. rotation in xy-plane and translation to center of xy-plane

                            theta= ntheta * theta_angle

                            do s=1,nseg
                                xp(s)=  chain_rot(1,s)*cos(theta)+chain_rot(2,s)*sin(theta) + xpt
                                yp(s)= -chain_rot(1,s)*sin(theta)+chain_rot(2,s)*cos(theta) + ypt
                                zp(s)=  chain_rot(3,s)   !+ zcm =  0                                 
                            enddo 
                       
                            do s=1,nseg

                                ! .. transformation to prism coordinates 
                                xpp(s) = ut(xp(s),yp(s))
                                ypp(s) = vt(xp(s),yp(s))

                                x(s) = pbc(xpp(s),Lx) ! .. periodic boundary conditions in x and y direction
                                y(s) = pbc(ypp(s),Ly)
                                z(s) = zp(s)          ! .. no pbc  in z- direction 

                                ! .. transforming form real- to lattice coordinates                 
                                xi = int(x(s)/delta)+1
                                yi = int(y(s)/delta)+1
                                zi = int(z(s)/delta)+1
                                    
                                call linearIndexFromCoordinate(xi,yi,zi,idx)
                                    
                                indexchain_init(s,conf) = idx

                                if(idx<=0.or.idx>nsize) then    
                                    text="Conformation outside box:"
                                    call print_to_log(LogUnit,text)  
                                    print*,text                          
                                    print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                                    info= myio_err_index
                                    return
                                endif
                                
                            enddo
                            
                            call VdWpotentialenergy(chain_rot,energy)  
                            energychain_init(conf)=energy

                            conf=conf+1   
                        
                        enddo ! .. rotation     

                    else 
                        text="Error: in make_chains_lajmmps_trj: geometry not cubic or prism: stopping program"
                        call print_to_log(LogUnit,text)
                        info = myio_err_geometry
                        return 
                    endif
                
                endif

            enddo      ! rotation loop
        endif    

    enddo          ! end while loop                                                                                                          
    !  .. end chains generation    
    

    if(conf<=max_confor) then
        print*,"subroutine make_chains_lammps_trj :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_lammps_trj"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isVdWintEne))energychain_init=0.0_dp

    close(un)

end subroutine read_chains_lammps_trj


! Reads confomations from a file called traj.xyz
! Format repeated lammps trajectory file 
! number of ATOMS much equal nseg 
! conformation is a loop molecule  

subroutine read_chains_XYZ_loop(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank, numproc                                                                                 
    use globals
    use chains
    use random
    use parameters
    use volume, only : position_graft, sgraft, nx, ny,nz, delta, nset_per_graft
    use volume, only : init_loop_rot_angle  
    use chain_rotation, only : rotationXaxis
    use myio, only : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only : myio_err_conf, myio_err_nseg, myio_err_geometry
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit


    ! .. argument

    integer, intent(out) :: info

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
    real(dp) :: chain(3,nseg)       ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1  
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    integer  :: xi,yi,zi
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(dp) :: theta 
    real(dp), allocatable, dimension(:,:) :: theta_array
    real(dp) :: xc,yc,zc               
    real(dp) :: energy                                             
    character(len=25) :: fname
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr

    ! .. executable statements   

    info=0

    ! get location graft points loop

    call read_graftpts_xyz_loop(info)
    if(info/=0) return

    ! .. open file   

    rankfile=mod(rank,nset_per_graft)                                                                                     
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='traj.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        endif
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        endif
    endif    


    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file    
    ios=0
    scalefactor=unit_conv
    energy=0.0_dp

    seed=435672               ! seed for random number generator                                                                               
    maxnchains= maxnchainsrotations  ! maximum number of rotation conf                                                                  
    maxntheta = maxnchainsrotationsxy! maximum number of rotation in xy-plane
    call init_loop_rot_angle(maxntheta,ngr,theta_array) ! array of angles of rotation 

    ios=0
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 
    xcm= Lx/2.0_dp          ! center x-y plane
    ycm= Ly/2.0_dp
    zcm= 0.0_dp

    do while ((conf<=max_confor).and.(ios==0))
    
        read(un,*,iostat=ios)nsegfile
        read(un,*,iostat=ios)str    
        
        if(conf==1) then
            if(nsegfile.ne.nseg) then 
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            endif    
        endif
    
        do s=1,nseg              ! .. read form  trajecotory file
            read(un,*,iostat=ios)xc,yc,zc
            xseg(1,s) = xc*scalefactor 
            xseg(2,s) = yc*scalefactor  
            xseg(3,s) = zc*scalefactor  
        enddo
     
        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        if(ios==0) then ! read was succesfull 

            conffile=conffile +1 
           
            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xgraftloop(1,1) 
                chain(2,s) = xseg(2,s)-xgraftloop(2,1)  
                chain(3,s) = xseg(3,s)-xgraftloop(3,1) 
            enddo
  
            select case (geometry)
            case ("cubic")

                g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                xpt =  position_graft(g,1)    ! position of graft point
                ypt =  position_graft(g,2)  

                do ntheta=1,maxntheta         ! rotation in xy-plane and translation to graft location of xy-plane

                    theta=theta_array(ntheta,g)

                    do s=1,nseg
                        xp(s)=  chain(1,s)*cos(theta)+chain(2,s)*sin(theta) + xpt
                        yp(s)= -chain(1,s)*sin(theta)+chain(2,s)*cos(theta) + ypt
                        zp(s)=  chain(3,s)                                  
                    enddo 
               
                    do s=1,nseg

                        x(s) = pbc(xp(s),Lx) ! periodic boundary conditions in x and y direction
                        y(s) = pbc(yp(s),Ly)
                        z(s) = zp(s)         ! no pbc in z- direction 

                        ! transforming form real- to lattice coordinates                 
                        xi = int(x(s)/delta)+1
                        yi = int(y(s)/delta)+1
                        zi = int(z(s)/delta)+1
                            
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                            
                        indexchain_init(s,conf) = idx

                        if(idx<=0.or.idx>nsize) then   

                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        endif
                        
                    enddo
                    
                    energychain_init(conf)=energy

                    conf=conf+1   
                
                end do ! .. rotation 
                        
            case("prism") 
                    
                g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                xpt =  position_graft(g,1)    ! position of graft point
                ypt =  position_graft(g,2)  

                do ntheta=1,maxntheta         ! rotation in xy-plane and translation to center of xy-plane

                    theta=theta_array(ntheta,g)

                    do s=1,nseg
                        xp(s)=  chain(1,s)*cos(theta)+chain(2,s)*sin(theta) + xpt
                        yp(s)= -chain(1,s)*sin(theta)+chain(2,s)*cos(theta) + ypt
                        zp(s)=  chain(3,s)   !+ zcm =  0                                 
                    end do 
               
                    do s=1,nseg

                        ! .. transformation to prism coordinates 
                        xpp(s) = ut(xp(s),yp(s))
                        ypp(s) = vt(xp(s),yp(s))

                        x(s) = pbc(xpp(s),Lx) ! .. periodic boundary conditions in x and y direction
                        y(s) = pbc(ypp(s),Ly)
                        z(s) = zp(s)          ! .. no pbc  in z- direction 

                        ! .. transforming form real- to lattice coordinates                 
                        xi = int(x(s)/delta)+1
                        yi = int(y(s)/delta)+1
                        zi = int(z(s)/delta)+1
                            
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                            
                        indexchain_init(s,conf) = idx

                        if(idx<=0.or.idx>nsize) then    
                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        end if
                        
                    end do
                    
                    energychain_init(conf)=energy

                    conf=conf+1   
                
                end do ! .. rotation     

            case default
                text="Error: in make_chains_XYZ_loop geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        end if   ! read 

    end do       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine make_chains_XYZ_loop :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_XYZ_loop"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain_init=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)
    
    deallocate(theta_array)

end subroutine read_chains_XYZ_loop


subroutine read_graftpts_xyz_loop(info)

    use mpivars, only : rank
    use parameters, only : unit_conv
    use myio, only : myio_err_chainsfile, myio_err_graft
    use myutils,  only : newunit
    use volume, only : sgraft,nset_per_graft  

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    character(len=25) :: fname
    integer :: ios 
    real(dp) :: xc,yc,zc          
    integer :: un,s, t
    integer :: rankfile
    integer :: item,moltype,nsegfile,idatom
    character(len=30) :: istr,str
    real(dp) :: scalefactor
    logical :: exist, isGraftItem

    ! .. executable statements 
    info=0

    !. . open file    
    rankfile=mod(rank,nset_per_graft)                                                                                            
    write(istr,'(I4)')rankfile
    fname='traj-graft.'//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj-graft.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    end if
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    end if

    ! read preamble/header
   
    scalefactor=unit_conv
    isGraftItem=.false.
    do s=1,2
        read(un,*,iostat=ios)item,xc,yc,zc
        if(item==sgraft) then
            t=1
            isGraftItem=.true.
        else
            t=2
        end if    
        xgraftloop(1,t)=xc*scalefactor
        xgraftloop(2,t)=yc*scalefactor
        xgraftloop(3,t)=zc*scalefactor    
    end do
    
    close(un)

    if(.not.isGraftItem) then
        print*,'Error read_graftpts_xyz_loop: sgraft not in traj file.'
        info = myio_err_graft
    end if
        
end subroutine



! Reads confomations from a file called traj.xyz
! Format repeated lammps trajectory file 
! number of ATOMS much equal nseg
! conformation is a linear molecule 

subroutine read_chains_XYZ_linear(info)

    !     .. variable and constant declaractions  
    use mpivars, only : rank, numproc                                                                                 
    use globals
    use chains
    use random
    use parameters
    use volume, only : position_graft, sgraft, nx, ny,nz, delta, nset_per_graft
    use volume, only : init_loop_rot_angle  
    use chain_rotation, only : rotationXaxis,rotationZcorr3
    use myio, only : myio_err_chainsfile, myio_err_energyfile, myio_err_index
    use myio, only : myio_err_conf, myio_err_nseg, myio_err_geometry
    use myutils,  only :  print_to_log, LogUnit, lenText, newunit


    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    integer :: i,j,s,rot,g,gn ,k     ! dummy indices
    integer :: idx                 ! index label
    integer :: ix,iy,iz,idxtmp,ntheta
    integer :: nchains              ! number of rotations
    integer :: maxnchains           ! number of rotations
    integer :: maxntheta            ! maximum number of rotation in xy-plane
    integer :: nchain, rotmax,  maxattempts
    integer :: conf,conffile        ! counts number of conformations  
    integer :: nsegfile             ! nseg in chain file      
    integer :: cuantasfile          ! cuantas in chain file                                              
    real(dp) :: chain(3,nseg)       ! chains(x,i)= coordinate x of segement i ,x=2 y=3,z=1   
    real(dp) :: chain_rot(3,nseg) ! 
    real(dp) :: xseg(3,nseg)
    real(dp) :: x(nseg), y(nseg), z(nseg)    ! coordinates
    real(dp) :: xp(nseg), yp(nseg), zp(nseg) ! coordinates 
    real(dp) :: xpp(nseg),ypp(nseg)
    integer  :: xi,yi,zi
    real(dp) :: Lx,Ly,Lz,xcm,ycm,zcm ! sizes box and center of mass box
    real(dp) :: xpt,ypt              ! coordinates
    real(dp) :: theta 
    real(dp), allocatable, dimension(:,:) :: theta_array
    real(dp) :: xc,yc,zc               
    real(dp) :: energy                                             
    character(len=25) :: fname
    integer :: ios, rankfile, iosene
    character(len=30) :: str
    real(dp) :: scalefactor
    integer :: un,unw,un_ene ! unit number
    logical :: exist
    character(len=lenText) :: text,istr
    logical :: is_positive_rot

    ! .. executable statements   

    info=0

    ! get location graft points linear

    call read_graftpts_xyz_linear(info)
    if(info/=0) return

    ! .. open file   

    rankfile=mod(rank,nset_per_graft)  
    
    write(istr,'(I4)')rankfile
    fname='traj.'//trim(adjustl(istr))//'.xyz'
    
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    end if
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    end if

    if(isChainEnergyFile) then
        write(istr,'(I4)')rankfile
        fname='energy.'//trim(adjustl(istr))//'.ene'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un_ene),file=fname,status='old',iostat=ios)
        else
            text='traj.'//trim(adjustl(istr))//'.ene file does not exit'
            print*,text
            info = myio_err_energyfile
            return
        end if
        if(ios >0 ) then
            print*, 'Error opening file : iostat =', ios
            info = myio_err_energyfile
            return
        end if
    end if    


    conf=1                    ! counter for conformations                                                           
    conffile=0                ! counter for conformations in file    
    ios=0
    scalefactor=unit_conv
    energy=0.0_dp

    seed=435672               ! seed for random number generator                                                                               
    maxnchains= maxnchainsrotations  ! maximum number of rotation conf                                                                  
    rotmax=maxnchains
    maxntheta = maxnchainsrotationsxy! maximum number of rotation in xy-plane
    call init_loop_rot_angle(maxntheta,ngr,theta_array) ! array of angles of rotation 

    ios=0
    Lz= nz*delta            ! maximum height box 
    Lx= nx*delta            ! maximum width box 
    Ly= ny*delta            ! maximum depth box 

    do while ((conf<=max_confor).and.(ios==0))
    
        read(un,*,iostat=ios)nsegfile
        read(un,*,iostat=ios)str    
        
        if(conf==1) then
            if(nsegfile.ne.nseg) then 
                text="nseg chain file not equal internal nseg : stop program"
                call print_to_log(LogUnit,text)
                info=myio_err_nseg 
                return
            end if    
        end if
    
        do s=1,nseg              ! .. read form  trajecotory file
            read(un,*,iostat=ios)xc,yc,zc
            xseg(1,s) = xc*scalefactor 
            xseg(2,s) = yc*scalefactor  
            xseg(3,s) = zc*scalefactor  
        end do
     
        if(isChainEnergyFile) read(un_ene,*,iostat=ios)energy

        if(ios==0) then ! read was succesfull 

            conffile=conffile +1 
           
            do s=1,nseg        
                chain(1,s) = xseg(1,s)-xgraftlinear(1) 
                chain(2,s) = xseg(2,s)-xgraftlinear(2)  
                chain(3,s) = xseg(3,s)-xgraftlinear(3) 
            end do
            
            select case (geometry)
            case ("cubic")

                g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                xpt =  position_graft(g,1)    ! position of graft point
                ypt =  position_graft(g,2)  

                do ntheta=1,maxntheta         ! rotation in xy-plane and translation to graft location of xy-plane

                    theta=theta_array(ntheta,g)

                    do s=1,nseg
                        xp(s)=  chain_rot(1,s)*cos(theta)+chain_rot(2,s)*sin(theta) + xpt
                        yp(s)= -chain_rot(1,s)*sin(theta)+chain_rot(2,s)*cos(theta) + ypt
                        zp(s)=  chain_rot(3,s)                                  
                    enddo 
               
                    do s=1,nseg

                        x(s) = pbc(xp(s),Lx) ! periodic boundary conditions in x and y direction
                        y(s) = pbc(yp(s),Ly)
                        z(s) = zp(s)         ! no pbc in z- direction 

                        ! transforming form real- to lattice coordinates                 
                        xi = int(x(s)/delta)+1
                        yi = int(y(s)/delta)+1
                        zi = int(z(s)/delta)+1
                            
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                            
                        indexchain_init(s,conf) = idx

                        if(idx<=0.or.idx>nsize) then   

                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return 
                        endif
                        
                    enddo
                    
                    energychain_init(conf)=energy

                    conf=conf+1   
                

                enddo   ! .. rotation 
                        
            case("prism") 
                    
                g=int(rank/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)

                xpt =  position_graft(g,1)    ! position of graft point
                ypt =  position_graft(g,2) 


                do ntheta=1,maxntheta         ! rotation in xy-plane and translation to center of xy-plane

                    theta=theta_array(ntheta,g)

                    do s=1,nseg
                        xp(s)=  chain(1,s)*cos(theta)+chain(2,s)*sin(theta) + xpt
                        yp(s)= -chain(1,s)*sin(theta)+chain(2,s)*cos(theta) + ypt
                        zp(s)=  chain(3,s)   !+ zcm =  0                                 
                    enddo 
               
                    do s=1,nseg

                        ! .. transformation to prism coordinates 
                        xpp(s) = ut(xp(s),yp(s))
                        ypp(s) = vt(xp(s),yp(s))

                        x(s) = pbc(xpp(s),Lx) ! .. periodic boundary conditions in x and y direction
                        y(s) = pbc(ypp(s),Ly)
                        z(s) = zp(s)          ! .. no pbc  in z- direction 

                        ! .. transforming form real- to lattice coordinates                 
                        xi = int(x(s)/delta)+1
                        yi = int(y(s)/delta)+1
                        zi = int(z(s)/delta)+1
                            
                        call linearIndexFromCoordinate(xi,yi,zi,idx)
                            
                        indexchain_init(s,conf) = idx

                        if(idx<=0.or.idx>nsize) then    
                            text="Conformation outside box:"
                            call print_to_log(LogUnit,text)  
                            print*,text                          
                            print*,"index=",idx, " xi=",xi," yi=",yi," zi=",zi, "conf=",conf,"s=",s 
                            info= myio_err_index
                            return
                        endif
                        
                    enddo
                    
                    energychain_init(conf)=energy

                    conf=conf+1   
                
                enddo ! .. rotation
                    

            case default

                text="Error: in make_chains_XYZ_linear geometry not cubic or prism: stopping program"
                call print_to_log(LogUnit,text)
                info = myio_err_geometry
                return 
                    
            end select

        endif   ! read 

    enddo       ! end while loop                                                                                                          
    ! end chains generation    
    
    conf=conf-1  ! lower by one  

    if(conf<max_confor) then
        print*,"subroutine make_chains_XYZ_loop :" 
        print*,"conf     = ",conf," less then imposed max cuantas     = ",max_confor
        print*,"conffile = ",conffile
        cuantas=conf   
        info=myio_err_conf        
    else
        text="Chains generated: subroutine make_chains_XYZ_loop"
        call print_to_log(LogUnit,text)
        readinchains=conffile
        info=0
    endif


    write(istr,'(L2)')isVdWintEne
    text="isVdWintEne = "//trim(adjustl(istr))
    call print_to_log(LogUnit,text)

    if(.not.(isChainEnergyFile)) energychain_init=0.0_dp

    close(un) 
    if(isChainEnergyFile) close(un_ene)
    
    deallocate(theta_array)

end subroutine read_chains_XYZ_linear

subroutine read_graftpts_XYZ_linear(info)

    use mpivars, only : rank
    use parameters, only : unit_conv
    use myio, only : myio_err_chainsfile, myio_err_graft
    use myutils,  only : newunit
    use volume, only : sgraft,nset_per_graft  

    ! .. argument

    integer, intent(out) :: info

    ! .. local variables

    character(len=25) :: fname
    integer :: ios, un, item
    real(dp) :: xc,yc,zc          
    integer :: rankfile
    character(len=30) :: istr,str
    real(dp) :: scalefactor
    logical :: exist

    ! .. executable statements 
    info=0

    !. . open file    
    rankfile=mod(rank,nset_per_graft)                                                                                            
    write(istr,'(I4)')rankfile
    fname='traj-graft.'//trim(adjustl(istr))//'.xyz'
    inquire(file=fname,exist=exist)
    if(exist) then
        open(unit=newunit(un),file=fname,status='old',iostat=ios)
    else
        print*,'traj-graft.rank.xyz file does not exit'
        info = myio_err_chainsfile
        return
    endif
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        info = myio_err_chainsfile
        return
    endif

    ! read preamble/header
   
    scalefactor=unit_conv
    
    read(un,*,iostat=ios)item,xc,yc,zc
    xgraftlinear(1)=xc*scalefactor
    xgraftlinear(2)=yc*scalefactor
    xgraftlinear(3)=zc*scalefactor    
    

    close(un)

    if(item/=1) then
        print*,'Error read_graftpts_xyz_linear: segement number not 1 in traj file.'
        info = myio_err_graft
    endif
        
end subroutine read_graftpts_XYZ_linear


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
            end if
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
            end if
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
            end if
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
            end if
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
            end if
        enddo

    case('multi')
    
        call read_type_of_monomer(type_of_monomer,type_of_monomer_char,typesfname, nseg) 
        call make_type_table(ismonomer_of_type,type_of_monomer,nseg,nsegtypes)    
        
        do s=1,nseg
            isAmonomer(s)=(type_of_monomer_char(s)=="A")
        end do

    case default
        print*,"Wrong chaintype: aborting program"
        stop
    end select

end subroutine make_sequence_chain


subroutine read_sequence_copoly_from_file(info)

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
    end if
        
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
            end if
        end if    
    end do 
    if(s/=nseg) then 
        print*,"reached end of file before all elements read"
        info = 1
        stop "read sequence file failed"
    end if

end subroutine read_sequence_copoly_from_file


! select indexchain and energy that have a weightchain=.true.
! return actual number of conformations
! rational: need identical set of confors for loop of distance /volumesizes
! Allows moding of brush compressesd by second surface at z=nz
 
subroutine chain_filter()
    
    use globals, only : nseg, cuantas, max_confor
    use chains, only : indexchain,indexchain_init,energychain,energychain_init
    use volume, only : nz,coordinateFromLinearIndex

    integer :: conf, c, s,  count_seg
    integer :: indx, ix, iy, iz 

    c=0           ! counts allowed conformations 
   
    do conf=1,max_confor      ! loop of all polymer conformations to filter out allowed ones 
        count_seg=0
        do s=1,nseg
            indx=indexchain_init(s,conf)
            call coordinateFromLinearIndex(indx,ix,iy,iz)
            if(iz<=nz) count_seg=count_seg+1
        end do

        if (count_seg.eq.nseg) then
            c= c+1 ! conformation  is allowed  
            energychain(c)=energychain_init(conf)
            do s=1,nseg
                indexchain(s,c)=indexchain_init(s,conf)
            end do
        end if    
    end do    

    cuantas=c ! actual number of conformation   

    call normed_weightchains()


end subroutine  chain_filter


! Compute weight chain w=e^E/Tre^E and normalize
! layout conformations : there are nset_per_graft nset 
! each graft point has nset of confomations thus total numberf of nodes= ngr*nset_per_graft
! the confomation on different graft point are the same excpet for translation. 
! Each graft point is idivual normalizes thus we need only normalize one set 

subroutine normed_weightchains()

    use mpivars
    use globals, only : cuantas
    use chains, only : energychain, logweightchain
    use volume, only : nset_per_graft
   
    integer :: un, c, k
    real(dp) :: localsum, totalsum, logtotalsum

        
    !    assymetric only use the first graft point to find normalization
    localsum=0.0_dp    
    do k=0,nset_per_graft-1   
        if(k==rank) then    ! 
            do c=1,cuantas
                localsum=localsum+exp(energychain(c)) 
            enddo   
        endif
    enddo    
  
    call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
    call MPI_ALLREDUCE(localsum, totalsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD,ierr) 
    
    ! normalize
    logtotalsum=log(totalsum)
    !un=rank+100
    !write(un,*)"rank=",rank," ",localsum,totalsum,logtotalsum
    do c=1,cuantas
        logweightchain(c)=energychain(c)-logtotalsum
        !write(un,*)rank,c,energychain(c),logweightchain(c)
    enddo    

    !call make_histogram(400)

end subroutine



function minimum_chainenergy() result(min_chainenergy)

    use  mpivars
    use  globals, only :cuantas
    use  chains, only : energychain

    real(dp) :: min_chainenergy

    integer :: conf
    
    min_chainenergy =0.0_dp
    do conf=1,cuantas      
        if(min_chainenergy > energychain(conf)) min_chainenergy=energychain(conf)
    enddo
   
end function


subroutine global_minimum_chainenergy()

    use  mpivars
    use  globals, only : cuantas
    use  chains, only : energychain, energychain_min
    use  parameters, only: isEnergyShift

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
    
    energychain_min =globalmin(1)

    if(isEnergyShift) then        
        call shift_energy_chain(globalmin(1))
        print*,"Shift minimum by :",energychain_min 
    endif
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
        case ("brush_mul","brush_mulnoVdW","brush","brush_neq","brushvarelec","brushborn","brushdna")
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

function VdWpotentialenergy_MC(chain,nchains)result(Energy)

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

    Lz= nz*delta            ! maximum height box s
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

end function VdWpotentialenergy_MC
 
! pre : chain conformation 
! post: VdW of conformation and if conformation is saw or not

subroutine VdWpotentialenergySaw(chain,Energy,saw)

    use globals, only : nseg, nsegtypes
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in)  :: chain(:,:)
    real(dp), intent(out) :: Energy
    logical, intent(out) :: saw

    real(dp) :: Ene,sqrlseg,sqrdist
    integer :: i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
   
    Ene=0.0_dp      
    saw=.true.

    do i=1,nseg
        do j=i+1,nseg

            s=type_of_monomer(i)
            t=type_of_monomer(j)

            zi = chain(1,i)
            xi = chain(2,i)
            yi = chain(3,i)

            zj = chain(1,j)
            xj = chain(2,j)
            yj = chain(3,j) 

            sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
            sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

            Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
            
            if(j/=(i+1).and.sqrdist<sqrlseg) then 
                saw=.false.
             !   print*,"overlap occured for i= ",i," and j= ",j 
            endif    
        enddo
    enddo 

    !print*,"total ",Ene," saw ",saw,' ',(lsegAA(t),t=1,nsegtypes),VdWeps 
    Energy = Ene            

end subroutine VdWpotentialenergySaw
   

 
! pre : chain conformation 
! post: VdW of conformation 
! Warning: need to to add VdWcutoff

subroutine VdWpotentialenergy(chain,Energy)

    use globals, only : nseg, nsegtypes
    use chains, only : type_of_monomer
    use parameters, only :  lsegAA,VdWeps

    real(dp), intent(in)  :: chain(:,:)
    real(dp), intent(out) :: Energy
    
    real(dp) :: Ene,sqrlseg,sqrdist !,maxdist
    integer ::  i,j,s,t
    real(dp) :: xi,xj,yi,yj,zi,zj
    
    Ene=0.0_dp 

    !maxdist=VdWcutoff*delta 
   
    do i=1,nseg
        do j=i+1,nseg

            s=type_of_monomer(i)
            t=type_of_monomer(j)

            zi = chain(1,i)
            xi = chain(2,i)
            yi = chain(3,i)

            zj = chain(1,j)
            xj = chain(2,j)
            yj = chain(3,j) 

            sqrdist=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
            sqrlseg=((lsegAA(t)+lsegAA(s))/2.0_dp)**2

            if(sqrdist<sqrlseg) sqrdist=sqrlseg

            Ene=Ene - VdWeps(s,t)*(sqrlseg/sqrdist)**3
            
        enddo
    enddo 

    Energy = Ene            

end subroutine VdWpotentialenergy

        
end module chaingenerator
