!     makes volume elements 
!     for spherical coordinates 
module volume      

    use precision_definition
    use mpivars

    implicit none
  
    !     .. variables
   
    real(dp) :: delta               ! delta  spacing of lattice site in x-, y- and z-direction  
    integer :: nz                   ! nz number of lattice sites in z-direction 
    integer :: nx                   ! nx number of lattice sites in x-direction 
    integer :: ny                   ! ny number of lattice sites in y-direction 
    integer :: nsurf                ! nsurf numer of lattice site  at z=0 ore z=nz*delta
    integer :: nzmax                ! nzmax  maximum number of lattice sites in z-direction
    integer :: nzmin                ! nzmin minimumal number of lattice sites in z-direction  
    integer :: nzstep               ! nzstep number of lattice sites stepped over or reduced 
    integer :: ngr                  ! total number of graft points    
    integer :: ngr_node             ! number of grafted areas assigned to an individual node
    integer :: ngr_freq             ! frequence spacing in terms of delta 
    integer :: ngrx                 ! total number of graft points along x-axis 
    integer :: ngry                 ! total number of graft points along y-axis                                                       
    logical :: isRandom_pos_graft   ! true random position, false regual pattern   
    real(dp) :: volcell             ! volcell=delta**3
    real(dp) :: areacell            ! areacell=delta**2
    real(dp) :: areasurf            ! area surfaces spanned by x and y direction
    real(dp) :: gamma               ! angle between oblique basis vectors u and v: cubic = gamma=90=pi/2 hexagonl gamma=60=2pi/3                          
    real(dp) :: beta                ! related beta = (pi/2- gamma)/2, angle between basis vector u and x and v and y
    real(dp) :: cos_two_beta        ! sqrt(cos(beta)**2 - sin(beta)**2)=cos(2beta) scales u and v coordinates 
    real(dp) :: sin_two_beta        ! sin(2beta)      
    integer  :: seed_graft          ! seed for graft points 
    integer  :: sgraft              ! item number of graft point to which loop is attached 
    integer  :: nset_per_graft      ! number of confomation set to read in graft point
    
    real(dp), dimension(:,:), allocatable :: position_graft

    character(len=11) :: geometry
  
    private :: beta

contains
    
    subroutine init_lattice

        use globals, only : nsize, DEBUG
        use mathconst

        implicit none 

        if(geometry=="cubic") then 
            gamma = pi/2.0_dp
        else
            gamma = gamma * pi /180.0_dp ! convert from degree to radians
        endif        

        beta = (pi/2.0_dp - gamma) / 2.0_dp       ! used by ut, vt, xt and yt functions
        cos_two_beta=cos(2*beta) ! sqrt(cos(beta)**2 - sin(beta)**2)  ! scaling of u and v coordinates
        sin_two_beta=sin(2*beta)

        ! cubic lattice  or prism surface in x-y direction at z=0 and z=nz  
        nz=nzmax
        nsize = nx*ny*nz                ! total number of cells or layers
        volcell = delta*delta*delta*1.0_dp     ! volume of one latice volume 
        areacell = delta*delta
        nsurf = nx*ny  
        areasurf=nsurf*delta*delta
        ngrx=int(nx/ngr_freq)
        ngry=int(ny/ngr_freq)

        ngr = int(nx/ngr_freq)*int(ny/ngr_freq) ! number of surface elements to be end-grafted with chains

        ! check 
        if(.not.((mod(nx,ngr_freq).eq.0).and.(mod(ny,ngr_freq).eq.0))) then 
             print*,"ngr test failed: exiting"
             print*,"nx= ",nx," ny= ",ny," ngr_freq = ",ngr_freq
             stop
        endif    
        
        !nset_per_graft = int(size/ngr)
        !..testing
        if(ngr*nset_per_graft/=size) then
            print*,"nset_per_graft test failed: exiting"
            print*,"nset_per_graft=",nset_per_graft,"ngr=",ngr,"size=",size
            stop
        endif

        allocate(position_graft(ngr,2)) ! only after ngr has been established position_graft can be allocated

        call init_graftpoints()

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


    subroutine init_graftpoints()

        use random 
        use myutils, only : newunit,lenText


        integer :: i, j, ig
        real(dp) :: rnd
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
            ! regulat lattice

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
                rnd = (rands(seed)-0.5_dp)*1.5_dp
                position_graft(ig,1) = position_graft(ig,1) + delta*rnd
                rnd = (rands(seed)-0.5_dp)*1.5_dp
                position_graft(ig,2) = position_graft(ig,2) + delta*rnd
            enddo
        
        endif        
        
        write(fname,'(A18)')'positiongraft-rank'
        write(istr,'(I4)')rank
        fname=trim(fname)//trim(adjustl(istr))//'.dat'
        print*,fname
        print*,istr
        open(unit=newunit(un_pgpt),file=fname,iostat=ios)
        if(ios >0 ) then
             print*, 'Error opening postiongraftpt.dat file : iostat =', ios
        else

            do ig=1,ngr          
                write(un_pgpt,*)position_graft(ig,1),position_graft(ig,2) 
            enddo
            close(un_pgpt)
        endif
    
    end subroutine init_graftpoints


end module volume
  
