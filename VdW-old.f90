!    .. module file for VdWcoeff
module VdW
    use precision_definition
    implicit none
  
    real(dp), dimension(:,:), allocatable :: chis  ! matrix of interaction coeff
    integer, dimension(:), allocatable :: minj 
    integer, dimension(:), allocatable :: maxj
  

    private :: ipbc 

contains
  

  
  
subroutine allocate_matrix(M,n)
    implicit none

    real(dp), dimension(:,:), allocatable :: M
    integer, intent(in) :: n
    integer :: ier
    
    if (.not. allocated(M))  then 
       allocate(M(n,n),stat=ier)
       
       if( ier.ne.0 ) then
          print*, 'Allocation error in matrix M: stat =', ier
          stop
       endif
    else
       print*,'Allocation warning: matrix already allocated'
    endif
    
end subroutine allocate_matrix

 
! read in VdW coeff from file
! needs to be updates   
subroutine read_VdWcoeff()
    
    use globals
    use parameters
    use volume 
    
    implicit none
    
    
    
end subroutine read_VdWcoeff


! Monte Carlo simulation to deterimin VdW coefficients


subroutine MC_VdWcoeff(VdWcoeff,lseg)

    use mpivars
    use mathconst
    use volume, only : xt, yt, ut, vt
    use volume, only : delta
    use random 

    real(dp), intent(inout) :: VdWcoeff(-2:2, -2:2, -2:2)
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


    if(rank.eq.0) print*,'chi calculation'

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

    if(rank.eq.0)print*, 'Sum 5x5', sum

    sum = 0.0_dp

    do ix = -limit, limit
        do iy = -limit, limit
            do iz = -limit, limit
                sum = sum +  matriz(ix, iy, iz)/MCsteps*((3.0_dp*delta)**3)
            enddo
        enddo
    enddo

    if(rank.eq.0) print*, 'Sum Total', sum
 
end subroutine


! this compute contribution to Palpha

subroutine VdW_contribution_exp(VdWcoeff,rhopol,exppi)

    use globals, only : nsize
    use volume, only : nx,ny,nz,coordinateFromLinearIndex,linearIndexFromCoordinate
    
    real(dp), intent(in) :: VdWcoeff(:,:,:)
    real(dp), intent(in) :: rhopol(:)
    real(dp), intent(inout) :: exppi(:)

    integer :: id,ix,iy,iz, jx,jy,jz, ax,ay,az
    real(dp) :: protemp

    real(dp), dimension(:,:,:), allocatable ::rhopoltmp

    

    do id=1,nsize
        call coordinateFromLinearIndex(id, ix, iy, iz)
        rhopoltmp(ix,iy,iz) = rhopol(id)   
    enddo
    
    VdWstr = VdWeps

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
                exppi(id) = exppi(id)*exp(VdWstr*protemp)

            enddo
        enddo
    enddo

end subroutine


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

      
