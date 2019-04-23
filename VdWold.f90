!    .. module file for VdWcoeff
module VdW
    use precision_definition
    implicit none
  
    real(dp), dimension(:,:,:), allocatable :: VdWcoeff
    real(dp), dimension(:,:,:), allocatable :: rhopoltmp

    real(dp), dimension(:,:), allocatable :: chis  ! still needed for defenition in energy  module

    private :: ipbc, rhopoltmp, allocate_VdWcoeff
    public :: MC_VdWcoeff, VdW_contribution_exp 

contains

  
  
subroutine allocate_VdWcoeff
    implicit none

    integer :: range
    integer :: ier
    
    range=2

    if (.not. allocated(VdWcoeff))  then 
       allocate(VdWcoeff(-range:range, -range:range, -range:range),stat=ier)
       
       if( ier.ne.0 ) then
          print*, 'Allocation error in matrix VdWcoeff: stat =', ier
          stop
       endif
    else
       print*,'Allocation warning: matrix already allocated'
    endif
    
end subroutine allocate_VdWcoeff

  
subroutine allocate_auxdensity

    use volume, only : nx, ny, nz 
    implicit none

    integer :: ier
    
    if (.not. allocated(rhopoltmp))  then 
       allocate(rhopoltmp(nx,ny,nz),stat=ier)
       
       if( ier.ne.0 ) then
          print*, 'Allocation error in auxilary polymer density: stat =', ier
          stop
       endif
    else
       print*,'Allocation warning: auxilary polymer density already allocated'
    endif
    
end subroutine allocate_auxdensity
 
! read in VdW coeff from file
! needs to be updates   
! subroutine read_VdWcoeff()     
! end subroutine read_VdWcoeff


! Monte Carlo simulation to deterimin VdW coefficients


subroutine MC_VdWcoeff(lseg)

    use mpivars
    use mathconst
    use volume, only : xt, yt, ut, vt
    use volume, only : delta
    use random 
    use myutils

    !real(dp), intent(inout) :: VdWcoeff(-2:2, -2:2, -2:2)
    real(dp), intent(in)   :: lseg ! size  segment 
    
 
    integer :: MCsteps ! number of MC steps 
    integer :: ix, iy , iz
    real(dp) :: x,y,z, radius, u, v
    real(dp) :: rn
    integer :: limit, plimit
    parameter (limit = 3)

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
      
!      MCsteps = 10
    MCsteps = 100000000


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
        endd
    enddo

    if(rank.eq.0) then 
        write(rstr,'(F5.3)') sum
        text="VdWcoefficient calculattion : Sum Total = "//trim(adjustl(rstr))
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
    !real(dp), dimension(:,:,:), allocatable ::rhopoltmp

    
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

      
