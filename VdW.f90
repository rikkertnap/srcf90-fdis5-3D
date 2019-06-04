!    .. module file for VdWcoeff
!    .. module file for computation of Van der Waals interaction 
!       Van der Waals energy :
!       EvdW = - 1/2 \sum{a,b} \esplison_{ab}\int dr \int dr' \rho_a(r) V_{a,b}(|r-r'|) \rho_b(r')
!       with       V_ab(r) =  (l_ab/r)^6  for l_ab<r < alpha l_ab 
!                  V_ab(r) =  0 otherwise !              
!    .. module comnputes V_ab using Monter Carlo simulation


module VdW
    
    use precision_definition
    implicit none

    real(dp), dimension(:,:,:), allocatable :: VdWcoeff

    real(dp), dimension(:,:,:), allocatable :: VdWcoeffAA
    real(dp), dimension(:,:,:), allocatable :: VdWcoeffAB
    real(dp), dimension(:,:,:), allocatable :: VdWcoeffBB

    real(dp), dimension(:,:,:), allocatable :: rhopoltmp  
    real(dp), dimension(:,:,:), allocatable :: rhopolAtmp  
    real(dp), dimension(:,:,:), allocatable :: rhopolBtmp
    
    integer, parameter :: range = 2 
    integer, parameter :: MCsteps = 100000000
    
    integer, parameter :: VdW_err_allocation = 1
    integer, parameter :: VdW_err_vdwcoeff   = 2
    integer, parameter :: Vdw_err_systype    = 3

    private
    public :: make_VdWcoeff, VdW_contribution_exp,VdW_contribution_exp_diblock
    public :: VdW_energy_diblock,VdW_energy

contains

  
  
subroutine allocate_VdWcoeff(info)
    
    use globals, only : systype

    integer :: ier
    logical :: alloc_fail

    integer,  intent(out), optional :: info
    
    alloc_fail=.FALSE.    

    if(systype=="electVdWAB") then 
    
        if (.not. allocated(VdWcoeffAA))  then 
            allocate(VdWcoeffAA(-range:range, -range:range, -range:range),stat=ier)
            if( ier/=0) alloc_fail=.true.
        endif
    
        if (.not. allocated(VdWcoeffAA))  then 
            allocate(VdWcoeffAA(-range:range, -range:range, -range:range),stat=ier)
            if( ier/=0) alloc_fail=.true.
        endif
    
    else if(systype=="electA".or.systype=="dipolarweakA") then 
        
        if (.not. allocated(VdWcoeff))  then 
            allocate(VdWcoeff(-range:range, -range:range, -range:range),stat=ier)
            if( ier/=0) alloc_fail=.true.
        endif 
    else 
        print*, 'Allocate_VdWcoeff: wrong systype =', systype
        if(present(info)) info= VdW_err_systype
    endif      


    if(alloc_fail) then 
        print*,'Allocation warning: allocate_VdWcoeff failed'
        if(present(info)) info= VdW_err_allocation
    endif   

    
end subroutine allocate_VdWcoeff

  
subroutine allocate_auxdensity(info)

    use globals, only : systype
    use volume, only : nx, ny, nz 
   
    integer,  intent(out), optional :: info

    integer :: ier
    
    if(systype=="electVdWAB") then 

        if (.not. allocated(rhopoltmp))  then 
           allocate(rhopoltmp(nx,ny,nz),stat=ier)
           
            if( ier.ne.0 ) then
                print*, 'Allocation error in auxilary polymer density: stat =', ier
                if(present(info)) info= VdW_err_allocation
                return
           endif
        else
           print*,'Allocation warning: auxilary polymer density already allocated'
        endif
    else if(systype=="electA".or.systype=="dipolarweakA") then 

        if (.not. allocated(rhopolAtmp))  then 
           allocate(rhopolAtmp(nx,ny,nz),stat=ier)
           
            if( ier.ne.0 ) then
                print*, 'Allocation error in auxilary polymer density: stat =', ier
                if(present(info)) info= VdW_err_allocation
                return
            endif
        else
           print*,'Allocation warning: auxilary polymer density already allocated'
        endif

        if (.not. allocated(rhopolAtmp))  then 
           allocate(rhopolBtmp(nx,ny,nz),stat=ier)
           
           if( ier.ne.0 ) then
                print*, 'Allocation error in auxilary polymer density: stat =', ier
                if(present(info)) info= VdW_err_allocation
                return
           endif
        else
           print*,'Allocation warning: auxilary polymer density already allocated'
        endif

    else 
        print*, 'Allocate_auxdensity: wrong systype =', systype
        if(present(info)) info= VdW_err_systype
    endif    

end subroutine allocate_auxdensity
 


subroutine make_VdWcoeff(info)

    use globals, only : systype
    use volume, only : geometry
    use parameters, only : lsegPAA,lsegPS
    use myutils

    integer,  intent(out), optional :: info
    
    real(dp) :: lsegAB
    character(len=lenText) :: text
    integer :: info_allocate


    if (present(info)) info = 0

    call allocate_VdWcoeff(info_allocate)
    call allocate_auxdensity(info_allocate)

    if (info_allocate/=0) then 

         text="make_VdWcoeff allocation failure"
         call print_to_log(LogUnit,text)
         if(present(info)) info= VdW_err_allocation
        
    else   

        if(systype=="electVdWAB") then
      
             lsegAB=(lsegPAA+lsegPS)/2.0_dp

            call MC_VdWcoeff(lsegPAA, VdWcoeffAA)
            call MC_VdWcoeff(lsegPS , VdWcoeffBB)
            call MC_VdWcoeff(lsegAB , VdWcoeffAB)
        
        else if (systype=="electA") then 

            call MC_VdWcoeff(lsegPAA , VdWcoeff)

        else

            text="make_VdWcoeff called with wrong systype"
            call print_to_log(LogUnit,text)
            if(present(info)) info= VdW_err_systype
        
        endif      

    endif    

end subroutine




! Monte Carlo simulation to deterimine VdW coefficients

subroutine MC_VdWcoeff(lseg,VdWcoeff)

    use mpivars
    use mathconst
    use volume, only : xt, yt, ut, vt
    use volume, only : delta
    use random 
    use myutils

    real(dp), intent(in)   :: lseg ! size  segment 
    real(dp), intent(out) :: VdWcoeff(-range:range, -range:range, -range:range)
    

    integer :: ix, iy , iz
    real(dp) :: x,y,z, radius, u, v
    real(dp) :: rn
    integer :: limit, plimit
    parameter (limit = range+1) 

    real(dp) :: matriz(-limit:limit, -limit:limit, -limit:limit) ! matrix for chi

    integer :: i
    real(dp) :: sum
    character(len=lenText) :: text, rstr


    if(rank==0) then 
        text="VdWcoefficient calculation"
        call print_to_log(LogUnit,text)
    endif    

    call allocate_VdWcoeff()
    call allocate_auxdensity()

    sum = 0.0_dp

    do ix = -limit, limit
        do iy = -limit, limit
            do iz = -limit, limit
                matriz(ix, iy, iz) = 0.0_dp
            enddo
        enddo
    enddo

    seed = 1010
      
    do i = 1, MCsteps

        x = 3.0*(rands(seed)-0.5)*delta ! random number between -1.5 * delta and 1.5 * delta
        y = 3.0*(rands(seed)-0.5)*delta 
        z = 3.0*(rands(seed)-0.5)*delta 

        u = ut(x, y)
        v = vt(x, y)

        radius = sqrt(x**2 + y**2 + z**2) ! real space
 
        if(radius<=(1.5*delta)) then  ! It is not inside the cut-off sphere
            if(radius>=lseg) then     ! is within the sphere of the segment

                ! cell

                ix = anint(v/delta)   ! space grid
                iy = anint(u/delta) 
                iz = anint(z/delta)

                matriz(ix, iy, iz) = matriz(ix, iy, iz) + (lseg/radius)**6
            endif
        endif
                
    enddo

    do ix = -2, 2
        do iy = -2, 2
            do iz = -2, 2

                VdWcoeff(ix, iy, iz) = matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
                sum = sum +  matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
        
            enddo
        enddo
    enddo

    if(rank.eq.0)then 
        write(rstr,'(F5.3)') sum
        text="VdWcoefficient calculation : Sum 5x5 = "//trim(adjustl(rstr))
        call print_to_log(LogUnit,text)
    endif    

    sum = 0.0_dp

    do ix = -limit, limit
        do iy = -limit, limit
            do iz = -limit, limit
                sum = sum +  matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
            enddo
        enddo
    enddo

    if(rank.eq.0) then 
        write(rstr,'(F5.3)') sum
        text="VdWcoefficient calculation : Sum Total = "//trim(adjustl(rstr))
        call print_to_log(LogUnit,text)
    endif    
 
end subroutine


! this compute contribution to Palpha

subroutine VdW_contribution_exp(rhopol,exppi)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWeps

!    real(dp), intent(in) :: VdWcoeff(:,:,:)
    real(dp), intent(in) :: rhopol(:)
    real(dp), intent(inout) :: exppi(:)

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protemp
    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz) = rhopol(id)   
    enddo
    
    !VdWstr = VdWeps

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                protemp = 0.0_dp
                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            if(1<=jz.and.jz<=nz) then 
                                protemp=protemp + VdWcoeff(ax,ay,az)*rhopoltmp(jx, jy, jz)
                            endif    
                        enddo
                    enddo
                enddo
                
                call linearIndexFromCoordinate(ix,iy,iz,id)
                exppi(id) = exppi(id)*exp(VdWeps%val*protemp)

            enddo
        enddo
    enddo

end subroutine



! this compute contribution to Palpha

subroutine VdW_contribution_exp_diblock(rhopolA,rhopolB,exppiA,exppiB)


    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    use parameters, only : VdWepsAA, VdWepsAB, VdWepsBB

    real(dp), intent(in) :: rhopolA(:),rhopolB(:)
    real(dp), intent(inout) :: exppiA(:),exppiB(:)

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protempAA,protempAB,protempBA,protempBB
    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopolAtmp(ix,iy,iz) = rhopolA(id)  
        rhopolBtmp(ix,iy,iz) = rhopolB(id)    
    enddo
    
    !VdWstr = VdWeps

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                protempAA = 0.0_dp
                protempAB = 0.0_dp
                protempBA = 0.0_dp
                protempBB = 0.0_dp

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            if(1<=jz.and.jz<=nz) then 
                                protempAA=protempAA + VdWcoeffAA(ax,ay,az)*rhopolAtmp(jx, jy, jz)
                                protempAB=protempAB + VdWcoeffAB(ax,ay,az)*rhopolBtmp(jx, jy, jz)
                                protempBA=protempBA + VdWcoeffAB(ax,ay,az)*rhopolAtmp(jx, jy, jz)
                                protempBB=protempBB + VdWcoeffBB(ax,ay,az)*rhopolBtmp(jx, jy, jz)
                            endif    
                        enddo
                    enddo
                enddo
                
                call linearIndexFromCoordinate(ix,iy,iz,id)

                exppiA(id) = exppiA(id)*exp(VdWepsAA%val*protempAA+VdWepsAB%val*protempAB)
                exppiB(id) = exppiB(id)*exp(VdWepsAB%val*protempBA+VdWepsBB%val*protempAB)
            enddo
        enddo
    enddo

end subroutine


! this compute contribution to Palpha

function VdW_energy(rhopol)result(EVdW)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate,volcell
    use parameters, only : VdWeps


    real(dp), intent(in) :: rhopol(:)

    real(dp) :: EVdW
    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: Etemp

    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz) = rhopol(id)   
    enddo
    
    Etemp =0.0_dp

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            if(1<=jz.and.jz<=nz) then 
                                Etemp=Etemp + VdWcoeff(ax,ay,az)*rhopoltmp(jx, jy, jz)*rhopoltmp(ix, iy, iz)
                            endif    
                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo

    EVdW = VdWeps%val*Etemp*volcell/2.0_dp

end function



function VdW_energy_diblock(rhopolA,rhopolB)result(EVdW)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate,volcell
    use parameters, only : VdWepsAA, VdWepsBB,VdWepsAB


    real(dp), intent(in) :: rhopolA(:),rhopolB(:)

    real(dp) :: EVdW
    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: EtempAA,EtempAB,EtempBB

    
    !  allocate rhopol tmp 
    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopolAtmp(ix,iy,iz) = rhopolA(id)
        rhopolBtmp(ix,iy,iz) = rhopolB(id)   
    enddo
    
    EtempAA =0.0_dp
    EtempAB =0.0_dp
    EtempBB =0.0_dp

    do ix=1,nx
        do iy=1,ny
            do iz=1,nz

                do ax = -2,2 
                
                    jx = ix+ax
                    jx = ipbc(jx,nx) ! mod(jx-1+5*dimx, dimx) + 1
                
                    do ay = -2,2
                        jy = iy+ay
                        jy = ipbc(jy,ny) ! mod(jy-1+5*dimy, dimy) + 1
                
                        do az = -2,2
                            jz = iz+az
                            if(1<=jz.and.jz<=nz) then 
                                EtempAA=EtempAA + VdWcoeffAA(ax,ay,az)*rhopolAtmp(jx, jy, jz)*rhopolAtmp(ix, iy, iz)
                                EtempAB=EtempAB + VdWcoeffAB(ax,ay,az)*rhopolAtmp(jx, jy, jz)*rhopolBtmp(ix, iy, iz)
                                EtempBB=EtempBB + VdWcoeffBB(ax,ay,az)*rhopolBtmp(jx, jy, jz)*rhopolBtmp(ix, iy, iz)
                            endif    
                        enddo
                    enddo
                enddo

            enddo
        enddo
    enddo

    EVdW = (VdWepsAA%val*EtempAA+2.0_dp*VdWepsAB%val*EtempAB+VdWepsBB%val*EtempBB)*volcell/2.0_dp

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


  
end module VdW

      
