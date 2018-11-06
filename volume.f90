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
    real(dp) :: volcell             ! volcell=delta**3
    integer :: nzmax                ! nzmax  maximum number of lattice sites in z-direction
    integer :: nzmin                ! nzmin minimumal number of lattice sites in z-direction  
    integer :: nzstep               ! nzstep number of lattice sites stepped over or reduced 
    integer :: ngr                  ! total number of graft points    
    integer :: ngr_node             ! number of grafted areas assigned to an individual node
    integer :: ngr_freq             ! frequence spacing in terms of delta 
    integer :: ngrx                 ! total number of graft points along x-axis 
    integer :: ngry                 ! total number of graft points along y-axis 
 
    real(dp) :: gamma               ! angle oblique lattice: cubic = gamma=90 degree 
                                    !                        hexagonal gamma=45 degree     


    real(dp), dimension(:), allocatable :: zc ! z-coordinate
    real(dp), dimension(:), allocatable :: xc ! z-coordinate
    real(dp), dimension(:), allocatable :: yc ! z-coordinate
   
    real(dp), dimension(:), allocatable :: deltaG ! geometrical factor
    real(dp), dimension(:), allocatable :: Fplus ! factor in Poisson Eq  
    real(dp), dimension(:), allocatable :: Fmin

    character(len=11) :: geometry
  
contains
  
    subroutine allocate_geometry(nx,ny,nz)
    
        implicit none
        integer, intent(in) :: nx, ny, nz
        
        integer :: ntot

        ntot=nx*ny*nz    
        
    !    print*," allocate_geometry : ntot=",ntot,"rank=",rank

        allocate(xc(nx))
        allocate(yc(ny))
        allocate(zc(nz))
        allocate(deltaG(ntot))

!    allocate(Fplus(N))
!    allocate(Fmin(N))
    
    end subroutine allocate_geometry
  
    !     computes geometical factors: deltaG
    !     deltaG(p) = the finite volume element p(i,j,k)
    !     subspended by volume element [i-1][j-1][k-1]
    !     divided by volcell. For volume elemetn outside and not intersecting 
    !     with cylinder deltaG=1
      
    subroutine  volume_elements_cubic(nx,ny,nz)

        implicit none
        
        integer, intent(in) :: nx, ny, nz

        !     .. local variables
        real(dp) :: vtol                 ! tolerance test volume integral 
        parameter (vtol=1.0E-10_dp)
        integer :: ix, iy, iz, ntot, idx
        real(dp) :: vol, vtest

        vol=0.0_dp
    
        do iz=1,nz
            zc(iz)= (iz-0.5_dp) * delta  ! x-coordinate 
        enddo
        do iy=1,ny
            yc(iy)= (iy-0.5_dp) * delta  ! y-coordinate 
        enddo
        do ix=1,nx
            xc(ix)= (ix-0.5_dp) * delta  ! z-coordinate 
        enddo

        vol=0.0_dp
        
        ntot=nx*ny*nz
  
        do idx=1,ntot ! nsize=nx*ny*nz 
            deltaG(idx)=1.0_dp  
            vol=vol+ deltaG(idx)
        enddo
    
        vtest=ntot
    
        if(abs(vtest-vol)>=vtol) then 
            print*,"Warning vtest=",vtest," not equal to vol=",vol
        endif

      
    end subroutine
  

    subroutine init_lattice

        use globals, only : nsize, DEBUG
        use mathconst

        implicit none 

        select case (geometry) 

        case("cubic")   ! cubic lattice  surfac in x-y direction at z=0 and z=nz  
           
            nsize = nx*ny*nz                ! total number of cells or layers
            volcell = delta*delta*delta     ! volume of one latice volume 
            nsurf = nx*ny  
            ngrx=int(nx/ngr_freq)
            ngry=int(ny/ngr_freq)

            ngr = int(nx/ngr_freq)*int(ny/ngr_freq) ! number of surface elements to be end-grafted with chains
            ! check 
            if(.not.((mod(nx,ngr_freq).eq.0).and.(mod(ny,ngr_freq).eq.0))) then 
                print*,"ngr test failed: exiting"
                print*,"nx= ",nx," ny= ",ny," ngr_freq = ",ngr_freq
                stop
            endif    
            !     .. compute ngr_node = the number of grafted areas assigned to one node
            !     .. part of the parallelization of program
            ngr_node = int(ngr/size)
            !     ..testing
            if(ngr/=(ngr_node*size)) then
                print*,"ngr_node test failed: exiting"
                print*,"ngrnode=",ngr_node,"ngr=",ngr,"size=",size
                stop
            endif
           
        case("hexagonal") ! hexagonal lattica 
            nsize = nx*ny*nz                ! total number of cells or layers
            !volcell = delta*delta*delta     ! volume of one latice volume 
            nsurf = nx*ny
            ngr = int(nx/ngr_freq)*int(ny/ngr_freq) ! number of surface elements to be end-grafted with chains
            ngrx=int(nx/ngr_freq)
            ngry=int(ny/ngr_freq)
            ! check 
            if(.not.((mod(nx,ngr_freq).eq.0).and.(mod(ny,ngr_freq).eq.0))) then 
                print*,"ngr test failed: exiting"
                print*,"nx= ",nx," ny= ",ny," ngr_freq = ",ngr_freq
                stop
            endif    

            !     .. compute ngr_node = the number of grafted areas assigned to one node
            !     .. part of the parallelization of program
            ngr_node = int(ngr/size)
            !     ..testing
            if(ngr/=(ngr_node*size)) then
                print*,"ngr_node test failed: exiting"
                print*,"ngrnode=",ngr_node,"ngr=",ngr,"size=",size
                stop
            endif
        
        case("square") ! square lattice in x-z axes 
            ny=1      ! one layer in y-direction 
            nsize = nx*nz                ! total number of cells or layers
            !volcell = delta*delta*delta     ! volume of one latice volume 
            nsurf = nx

            ngr = int(nx/ngr_freq)! number of surface elements to be end-grafted with chains

            ! check 
            if(.not.(mod(nx,ngr_freq).eq.0)) then 
                print*,"ngr test failed: exiting"
                print*,"nx= ",nx,"ngr_freq = ",ngr_freq
                stop
            endif    

            !     .. compute ngr_node = the number of grafted areas assigned to one node
            !     .. part of the parallelization of program
            ngr_node = int(ngr/size)
            !     ..testing
            if(ngr/=(ngr_node*size)) then
                print*,"ngr_node test failed: exiting"
                print*,"ngrnode=",ngr_node,"ngr=",ngr,"size=",size
                stop
            endif
    
        case default
            print*,"Error: init_lattice: geometry not CUBIC, or PLANAR"
            print*,"stopping program"
            stop
        end select
            
    end subroutine init_lattice
         

    subroutine  make_geometry()
    
        use globals
        use myutils

        implicit none
    
        integer ::  i
        character(len=lenText) :: text, str 

        call init_lattice

        select case (geometry) 
        case("cubic")
            call volume_elements_cubic(nx,ny,nz)   
        case default
            print*,"Error: geometry not CUBIC"
            print*,"stopping program"
            stop
        end select

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


    ! coordinate transformation from cartesian x,y to oblique u,v  coordinates

    function ut(x, y, gamma) result(ut_val)
    
        use mathconst

        real(dp), intent(in) :: x, y, gamma
        real(dp) :: ut_val
        real(dp) :: beta 

        beta = (pi/2.0_dp - gamma) / 2.0_dp
        ut_val = -sin(beta)*x + cos(beta)*y
        ut_val = ut_val / sqrt(cos(beta)**2 - sin(beta)**2)
    
    end function ut

   ! coordinate transformation from cartesian x,y to oblique u,v  coordinates

    function vt(x, y, gamma) result(vt_val)
       
        use mathconst

        real(dp), intent(in) :: x, y, gamma
        real(dp) :: vt_val
        real(dp) :: beta 

        beta = (pi/2.0_dp - gamma) / 2.0_dp
        vt_val = cos(beta)*x - sin(beta)*y
        vt_val = vt_val / sqrt(cos(beta)**2 - sin(beta)**2)
    
    end function vt

    ! inverse coordinate transformation from oblique u,v to  cartesian x,y  coordinates

    function xt(v, u, gamma) result(xt_val)
        
        use mathconst

        real(dp), intent(in) :: u, v, gamma
        real(dp) :: xt_val
        real(dp) :: beta
    
        beta = (pi/2.0_dp - gamma) / 2.0_dp
        xt_val = cos(beta)*v + sin(beta)*u
        xt_val = xt_val / sqrt(cos(beta)**2 - sin(beta)**2)
    end function xt

    ! inverse coordinate transformation from oblique u,v to  cartesian x,y  coordinates

    function yt(v, u, gamma) result(yt_val)
        
        use mathconst

        real(dp), intent(in) :: u, v, gamma
        real(dp) :: yt_val
        real(dp) :: beta
    
        beta = (pi/2.0_dp - gamma) / 2.0_dp

        yt_val = sin(beta)*v + cos(beta)*u
        yt_val = yt_val / sqrt(cos(beta)**2 - sin(beta)**2)
    
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






end module volume
  
