!     makes volume elements 
!     for spherical coordinates 
module volume      

    use precision_definition
    ! use mpivars

    implicit none
  
    !     .. variables
   
    real(dp) :: delta               ! delta  spacing of lattice site in x-, y- and z-direction  
    integer :: nz                   ! nz number of lattice sites in z-direction 
    integer :: nx                   ! nx number of lattice sites in x-direction 
    integer :: ny                   ! ny number of lattice sites in y-direction 
    integer :: nsurf                ! nsurf numer of lattice siteß  at z=0 ore z=nz*delta
    integer :: nzmax                ! nzmax  maximum number of lattice sites in z-direction
    integer :: nzmin                ! nzmin minimumal number of lattice sites in z-direction  
    integer :: nzstep               ! nzstep number of lattice sites stepped over or reduced
    real(dp) :: volcell             ! volcell=delta**3
    real(dp) :: areacell            ! areacell=delta**2
    real(dp) :: areasurf            ! area surfaces spanned by x and y direction
    real(dp) :: gamma               ! angle between oblique basis vectors u and v: cubic = gamma=90=pi/2 hexagonl gamma=60=2pi/3                          
    real(dp) :: beta                ! related beta = (pi/2- gamma)/2, angle between basis vector u and x and v and y
    real(dp) :: cos_two_beta        ! sqrt(cos(beta)**2 - sin(beta)**2)=cos(2beta) scales u and v coordinates 
    real(dp) :: sin_two_beta        ! sin(2beta)  
    character(len=11) :: geometry

    ! variable for grafting position
    integer :: ngr                  ! total number of graft points    
    integer :: ngr_node             ! number of grafted areas assigned to an individual node
    integer :: ngr_freq             ! frequence spacing in terms of delta 
    integer :: ngrx                 ! total number of graft points along x-axis 
    integer :: ngry                 ! total number of graft points along y-axis                                                       
    logical :: isRandom_pos_graft   ! true random position, false regual pattern
    integer :: seed_graft           ! seed for graft points 
    real(dp) :: scale_ran_step      ! scale random postition,minimum value =1.0 =half ngr_freq
    integer :: sgraft               ! item number of graft point to which loop is attached 
    integer :: nset_per_graft       ! number of confomation set to read in graft point
    real(dp), dimension(:,:), allocatable :: position_graft
   
    ! variable for rotation loop chain
    logical  :: isRandom_rot_loop   ! true random rotation/oewinetation loop chain, false regaualr rotation 
    integer  :: seed_rot_loop       ! seed for rotation loops 

    ! hash table lattice 
    integer, dimension(:,:,:), allocatable :: coordtoindex 
    integer, dimension(:,:),  allocatable  :: indextocoord
    
    ! hash table lattice neighbors
    integer, dimension(:,:), allocatable   :: indexlatneighbor
    integer ::  maxlatneigh

    private :: beta
    private :: ipbc

contains
    
    subroutine init_lattice

        use globals, only : nsize, nsizepsi, DEBUG
        use mathconst
        use mpivars, only : numproc

        implicit none 

        if(geometry=="cubic") then 
            gamma = pi/2.0_dp
        else
            gamma = gamma * pi /180.0_dp ! convert from degree to radians
        endif        

        beta = (pi/2.0_dp - gamma) / 2.0_dp       ! used by ut, vt, xt and yt functions
        cos_two_beta = cos(2*beta) ! sqrt(cos(beta)**2 - sin(beta)**2)  ! scaling of u and v coordinates
        sin_two_beta = sin(2*beta)

        ! cubic lattice  or prism surface in x-y direction at z=0 and z=nz  
        nz = nzmax
        nsize = nx*ny*nz                ! total number of cells or layers
        nsizepsi = nsize + 2 * nx * ny  ! total number of cells for el potential  
        volcell = delta*delta*delta*1.0_dp     ! volume of one latice volume 
        areacell = delta*delta
        nsurf = nx*ny  
        areasurf=nsurf*delta*delta
        ngrx=int(nx/ngr_freq)
        ngry=int(ny/ngr_freq)
        ngr = int(nx/ngr_freq)*int(ny/ngr_freq) ! number of surface elements to be end-grafted with chains

        ! check 
        
        if(.not.(mod(nx,ngr_freq).eq.0)) then 
             print*,"ngr test failed in x-direction: exiting"
             print*,"nx= ",nx," ngr_freq = ",ngr_freq
             stop
        endif    
        
        if(.not.(mod(ny,ngr_freq).eq.0))then 
          print*,"ngr test failed in y-direction: exiting"
             print*,"ny= ",ny," ngr_freq = ",ngr_freq
             stop
        endif    
        
        if(ngr*nset_per_graft/=numproc) then
            print*,"nset_per_graft test failed: exiting"
            print*,"nset_per_graft=",nset_per_graft,"ngr=",ngr,"numproc=",numproc
            stop
        endif

        allocate(position_graft(ngr,2)) ! only after ngr has been established position_graft can be allocated

        call init_graftpoints()

        !   use of hash table        
        call allocate_hashtable(nx,ny,nz)
        call make_hashtable()

    end subroutine init_lattice
         

    subroutine  make_geometry()
    
        use globals
        use myutils

        implicit none
    
        character(len=lenText) :: text, str 

        call init_lattice

        write(str,'(A20)')geometry
        text="geometry="//trim(adjustl(str)) 
        call print_to_log(LogUnit,text) 


    end subroutine make_geometry

    subroutine linearIndexFromCoordinate(x,y,z,idx)
     
        implicit none 
      
        integer, intent(in)  :: x,y,z
        integer, intent(out) :: idx
      
        integer :: a,b,c,d

        a = 1
        b = nx 
        c = nx * ny 
        d = 1-a-b-c
        idx = a*x + b*y + c*z + d

    end subroutine linearIndexFromCoordinate

    subroutine linearIndexFromCoordinate_general(x,y,z,nx,ny,idx)

        implicit none

        integer, intent(in)  :: x,y,z,nx,ny
        integer, intent(out) :: idx

        integer :: a,b,c,d

        a = 1
        b = nx
        c = nx * ny
        d = 1-a-b-c
        idx = a*x + b*y + c*z + d

    end subroutine linearIndexFromCoordinate_general


    subroutine coordinateFromLinearIndex(idx, x,y, z)

        implicit none

        integer, intent(out)  :: x,y,z
        integer, intent(in)   :: idx
        integer :: idxtmp
   
        idxtmp=idx
        x =  mod(idxtmp-1,nx)+1
        idxtmp =int((idxtmp-1)/nx)+1
        y = mod(idxtmp-1,ny)+1
        idxtmp = int((idxtmp-1)/ny)+1
        z = idxtmp
    
    end subroutine coordinateFromLinearIndex

    subroutine  coordinateFromLinearIndex_general(idx,x,y,z,nx,ny)

        implicit none

        integer, intent(out)  :: x,y,z
        integer, intent(in)   :: idx, nx,ny
        integer :: idxtmp

        idxtmp=idx
        x =  mod(idxtmp-1,nx)+1
        idxtmp =int((idxtmp-1)/nx)+1
        y = mod(idxtmp-1,ny)+1
        idxtmp = int((idxtmp-1)/ny)+1
        z = idxtmp

    end subroutine coordinateFromLinearIndex_general

     
    function mirror_index(idx,nz) result(idx_mirror)

        implicit none 
      
        integer, intent(in)  :: idx, nz
        integer :: idx_mirror

        integer :: x,y,z

        call coordinateFromLinearIndex(idx,x,y,z)
        call LinearIndexFromCoordinate(x,y,nz+1-z,idx_mirror)

    end function mirror_index

    ! coordinate transformation from cartesian x,y to prism or oblique u,v coordinates
    ! u = s(cos(gamma)*x -sin(gamma)*y); s=1/srqt(cos(gamma)^2-sin(gamma)^2)
    ! v = s(-sin(gamma)*x +sin(gamma)*y); 
    ! metric ds^2= s^2 du^2 + s^2 dv^2 + ss^2 sin(2 gamma) du dv +dz^2

    function ut(x, y) result(ut_val)
       
        real(dp), intent(in) :: x, y
        real(dp) :: ut_val
        
        ut_val = cos(beta)*x - sin(beta)*y
        ut_val = ut_val / cos_two_beta
    
    end function ut


    ! coordinate transformation from cartesian x,y to prism u,v  coordinates

    function vt(x, y) result(vt_val)
    
        real(dp), intent(in) :: x, y
        real(dp) :: vt_val
       
        vt_val = -sin(beta)*x + cos(beta)*y
        vt_val = vt_val / cos_two_beta
    
    end function vt



    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates

    function xt(u, v) result(xt_val)
        
        real(dp), intent(in) :: u, v
        real(dp) :: xt_val
    
        xt_val = cos(beta)*u + sin(beta)*v
        xt_val = xt_val / cos_two_beta
    
    end function xt

    ! inverse coordinate transformation from prism u,v to  cartesian x,y  coordinates
    
    function yt(u, v) result(yt_val)

        real(dp), intent(in) :: u, v
        real(dp) :: yt_val
    
        yt_val = sin(beta)*u + cos(beta)*v
        yt_val = yt_val / cos_two_beta
    
    end function yt
 
    function isEven(number) result(val) 
        integer, intent(in) :: number
        logical :: val

        if ( mod(number, 2) ==0) then
            val=.True. 
        else
            val=.False.
        endif 
   
    end function

    subroutine allocate_hashtable(nx,ny,nz)

        use globals, only : nsize
        integer, intent(in) :: nx,ny,nz

        allocate(coordtoindex(nx,ny,nz))
        allocate(indextocoord(nsize,3))
        
    end subroutine allocate_hashtable

    subroutine make_hashtable

        use globals, only : nsize

        integer :: idx, ix, iy, iz
            
        do idx=1,nsize
            
            call coordinateFromLinearIndex(idx, ix, iy,iz)

            coordtoindex(ix,iy,iz)=idx
            
            indextocoord(idx,1)=ix
            indextocoord(idx,2)=iy
            indextocoord(idx,3)=iz

        enddo
        
    end subroutine make_hashtable


    subroutine init_graftpoints()

        use globals,  only : DEBUG
        use random 
        use myutils, only : newunit,lenText
        use mpivars, only : rank


        integer :: i, j, ig
        real(dp) :: rnd
        real(dp) :: u, v
        character(len=lenText) :: fname, istr
        integer :: ios, un_pgpt

        if(.not.isRandom_pos_graft) then

            ig = 1   

            do i=1,ngrx            ! location graft points 
                do j=1,ngry
                    position_graft(ig,1) = (i-0.5_dp)*delta*ngr_freq
                    position_graft(ig,2) = (j-0.5_dp)*delta*ngr_freq
                    ig = ig + 1
                enddo
            enddo
        
        else 
        
            ig = 1   

            do i=1,ngrx            ! location regular graft points 
                do j=1,ngry
                    position_graft(ig,1) = (i-0.5_dp)*delta*ngr_freq
                    position_graft(ig,2) = (j-0.5_dp)*delta*ngr_freq
                    ig = ig + 1
                enddo
            enddo

            
            seed=seed_graft         ! randomize

            do ig=1,ngr          
                rnd = (rands(seed)-0.5_dp)*ngr_freq/scale_ran_step
                position_graft(ig,1) = position_graft(ig,1) + delta*rnd
                rnd = (rands(seed)-0.5_dp)*ngr_freq/scale_ran_step
                position_graft(ig,2) = position_graft(ig,2) + delta*rnd
            enddo
        
        endif  
      
        ! position_graft is on a regular square pattern applicable to geometry ="cubic" 
        ! For hexagonal/oblique geometry one needs to tranform (u,v) to (x,y),
        ! then above square pattern becomes a hexagonal/oblique pattern. 

        if(geometry=="prism") then 

            do ig=1,ngr   
                u = position_graft(ig,1)
                v = position_graft(ig,2)  
                position_graft(ig,1) = xt(u, v)  
                position_graft(ig,2) = yt(u, v)
            enddo    
        
        endif   


        ! output
        if(DEBUG.or.rank==0) then
            write(fname,'(A18)')'positiongraft-rank'
            write(istr,'(I4)')rank
            fname=trim(fname)//trim(adjustl(istr))//'.dat'
            open(unit=newunit(un_pgpt),file=fname,iostat=ios)
            if(ios >0 ) then
                 print*, 'Error opening positiongraftpt.dat file : iostat =', ios
            else

                do ig=1,ngr          
                    write(un_pgpt,*)position_graft(ig,1),position_graft(ig,2) 
                enddo
                close(un_pgpt)
            endif
        endif    

    end subroutine init_graftpoints

    subroutine init_loop_rot_angle(maxntheta,ngr,theta_angle)

        use globals, only : DEBUG
        use mathconst
        use random, only : seed, rands
        use myutils, only : newunit,lenText
        use mpivars, only : rank 


        integer, intent(in) :: maxntheta,ngr
        real(dp), allocatable, dimension(:,:), intent(inout) :: theta_angle

        ! local variables

        integer :: i, g
        real(dp) :: delta_theta_angle, rnd  
        character(len=lenText) :: fname, istr
        integer :: ios, un_rotgpt

        ! maxntheta = maxnchainsrotationsxy  maximum number of rotation in xy-plane 

        allocate(theta_angle(maxntheta,ngr))

        delta_theta_angle = 2.0_dp*pi/maxntheta ! angle of rotation   
        
        do g=1,ngr
            do i=1,maxntheta           ! location graft points 
                theta_angle(i,g) = i * delta_theta_angle
            end do
        end do    

        if(isRandom_rot_loop) then   
            
            seed=seed_rot_loop     ! randomize
            do g =1,ngr
                do i=1,maxntheta          
                    rnd = (rands(seed)-0.5_dp)*delta_theta_angle/2.0_dp 
                    theta_angle(i,g) = theta_angle(i,g) + rnd
                end do
            end do              
        end if        
    

        ! output
        if(DEBUG.or.rank==0) then
            write(fname,'(A18)')'rotationgraft-rank'
            write(istr,'(I4)')rank
            fname=trim(fname)//trim(adjustl(istr))//'.dat'
            open(unit=newunit(un_rotgpt),file=fname,iostat=ios)
            if(ios >0 ) then
                 print*, 'Error opening rotationgraftpt.dat file : iostat =', ios
            else
                do g=1,ngr
                    do i=1,maxntheta          
                        write(un_rotgpt,*)theta_angle(i,g) 
                    end do
                end do
                close(un_rotgpt)
            end if
        end if    

    end subroutine init_loop_rot_angle



    subroutine allocate_indexlatneighbor(nsize,maxlatneigh)

        integer, intent(in) :: nsize
        integer, intent(in) :: maxlatneigh

        allocate(indexlatneighbor(nsize,maxlatneigh))

    end subroutine allocate_indexlatneighbor



    ! Calculate indexlatneighbor(idx,k) 
    ! indexlatneighbor(idx,k) = return lattice cell index of (lattice) neighbor number k  of lattice cell index idx
    !                         = return index of  kth lattice neigbor asccoated with lattice cell index idx                   
    ! input real(dp) : distphoscutoff : maximum distance or separation of a phosphate pair
    ! output integer,   indexlatneighbor(idx,k) 
            

    subroutine make_table_index_neighbors(distphoscutoff)

        use globals, only : nsize

        real(dp), intent(in) :: distphoscutoff
        
        integer :: idx, i, j, k, ix, iy, iz, deltaix, deltaiy, deltaiz
        integer :: latneighbornumber, idxneigh
        integer :: ip, jp, kp
        integer :: sqrdist, rangecutoff
        real(dp) :: sqrphoscutoff_delta


        rangecutoff=int(distphoscutoff/delta)+2  ! 2 just to be sure !!

        sqrphoscutoff_delta= (distphoscutoff/delta)**2

        do idx=1,nsize

            ! use hash table to get coordiantes 
            ix=indextocoord(idx,1)
            iy=indextocoord(idx,2)
            iz=indextocoord(idx,3)

            latneighbornumber=0 ! 

            do deltaix=-rangecutoff,rangecutoff
                i=ix+deltaix
                do deltaiy=-rangecutoff,rangecutoff
                    j=iy+deltaiy
                    do deltaiz=-rangecutoff,rangecutoff
                        k=iz+deltaiz
                        sqrdist = deltaix**2 + deltaiy**2 +deltaiz**2  ! square distance in unit of delta squared   
                        if (sqrdist <= sqrphoscutoff_delta) then 
                            ! apply pbc 
                            ip=ipbc(i,nx)
                            jp=ipbc(j,ny)
                            kp=ipbc(k,nz)

                            latneighbornumber = latneighbornumber+1
                            idxneigh = coordtoindex(ip,jp,kp) ! index of neighbour
                            indexlatneighbor(idx,latneighbornumber) = idxneigh
                        end if     

                    end do
                end do
            end do

        end do    

    end subroutine make_table_index_neighbors


    ! Calculate maxlatneightbor : total number of lattice cells neighbor any lattice cell that is with a distance  smaller of 
    ! equal to maximum allow  distance or separation of a phosphate pair : distphoscutoff 
    ! input real(dp) : distphoscutoff : maximum distance or separation of a phosphate pair
    ! input real(dp) :: delta <=> defined in  global
       ! output integer maxlatneigh
    
    function set_maxlatneigh(delta, distphoscutoff) result(maxlatneigh)

        real(dp), intent(in) :: delta
        real(dp), intent(in) :: distphoscutoff
        integer :: maxlatneigh

        ! local variables 

        integer :: deltaix, deltaiy, deltaiz
        integer :: latneighbornumber
        integer :: sqrdist, rangecutoff
        real(dp) :: sqrphoscutoff_delta

        rangecutoff=int(distphoscutoff/delta)+2  ! 2 just to be sure !!
        sqrphoscutoff_delta= (distphoscutoff/delta)**2

        latneighbornumber = 0 ! 

        do deltaix=-rangecutoff,rangecutoff
            do deltaiy=-rangecutoff,rangecutoff
                do deltaiz=-rangecutoff,rangecutoff
                    sqrdist = deltaix**2 + deltaiy**2 +deltaiz**2  ! square distance in unit of delta squared   
                    if (sqrdist<=sqrphoscutoff_delta) latneighbornumber = latneighbornumber+1
                    
                end do
            end do
        end do    

        maxlatneigh = latneighbornumber

    end function  set_maxlatneigh


    ! compute periodic boundary condition in integer units 
    ! real(dp) version  : pbbc is located in chaingenerator  

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



end module volume
  
