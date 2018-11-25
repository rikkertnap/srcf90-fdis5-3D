!    .. module file for VdWcoeff
module VdW
  use precision_definition
  implicit none
  
  real(dp), dimension(:,:), allocatable :: chis  ! matrix of interaction coeff
  integer, dimension(:), allocatable :: minj 
  integer, dimension(:), allocatable :: maxj
  
contains
  
  integer function minrange(i)
    
    use parameters, only :  VdWcutoffdelta
    implicit none
    
    integer, intent(in) :: i

    if (i-VdWcutoffdelta>1) then
       minrange=i-VdWcutoffdelta
    else 
       minrange=1
    endif
    
  end function minrange
  
  integer function maxrange(i,n)

    use parameters, only :  VdWcutoffdelta
    implicit none

    integer, intent(in) :: i,n
   
    if(i+VdWcutoffdelta<=n) then
       maxrange=i+VdWcutoffdelta
    else
       maxrange=n
    endif

  end function maxrange

  
  
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

  
    subroutine read_VdWcoeff()
        
        use globals
        use parameters
        use volume 
        
        implicit none
        
        character(len=9) :: fname
        integer :: ios
        
        !     .. local variables 
        integer :: nsize_file
        real(dp)  :: delta_file
        real(dp)  :: chibulk_file
        integer :: cutoff_file
        integer :: i,j !,nradius
        integer :: ival,jval
        
        real(dp), dimension(:,:), allocatable :: VdWcoeff  ! matrix of interaction coeff

        !     .. reading in of variables from file                                              
        
        if(sysflag=="elect") then 
            write(fname,'(A9)')'chisC.dat'
            VdWCutoffdelta=int((VdWcutoff*lsegC)/delta)+2
        endif   
        if(sysflag=="neutral") then 
            write(fname,'(A9)')'chisB.dat'
            VdWCutoffdelta=int((VdWcutoff*lsegB)/delta)+2 
        endif

        open(unit=1,file=fname,iostat=ios,status='old')
        if(ios >0 ) then
           print*, 'Error opening file : iostat =', ios
           stop
        endif
        
        read(1,*)nsize_file
        read(1,*)delta_file
        read(1,*)cutoff_file
        read(1,*)chibulk_file
        
        chibulk=chibulk_file
        
          

        if(delta_file /= delta) then 
           print*,'delta in VdWcoeff file not equal delta'
           stop
        endif
        if(cutoff_file /= VdWcutoff) then 
           print*,'cuttoff in VdWcoeff file not equal imposed cutoff'
           stop
        endif
        
        if(nsize_file <= nsize+VdWcutoff) then
           print*,'Size in VdWcoeff file not large enough'
           print*,'nsize_file        =',nsize_file
           stop
        endif
          
        call allocate_matrix(chis,nsize)
        call allocate_matrix(VdWcoeff,nsize_file)
        
        ! read file
        do i=1,nsize_file
           do j=1,nsize_file
              read(1,*)ival,jval,VdWcoeff(i,j)
           enddo
        enddo

        ! transform VdW coefficients to chi values
        do i=1,nsize
           do j=1,nsize
              chis(i,j)=-1.0_dp*VdWcoeff(i,j)
           enddo
        enddo
        
        deallocate(VdWcoeff)

        close(1)
        
        ! make range maxj and minj to be used in fcnCa.f90
      
        allocate(minj(nsize))
        allocate(maxj(nsize))

        do i=1,nsize
           minj(i)=minrange(i)
           maxj(i)=maxrange(i,nsize)
        enddo
        
    end subroutine read_VdWcoeff



    subroutine VdWcoeff

        use mpivars
        use mathconst
        use volume, only : xt, yt, ut, vt
        use volume, only : delta
        use random 


        real(dp) :: lseg ! largo del segmento
        real(dp) :: l ! medio lseg, radio del segmento

        real(dp) :: Xu(-2:2, -2:2, -2:2)
     
        integer :: MCsteps ! numero de steps de MC
        integer :: ix, iy , iz
        real(dp) :: x,y,z, radius, u, v
        real(dp) :: rn
        integer :: limit, plimit
        parameter (limit = 3)

        real(dp) :: matriz(-limit:limit, -limit:limit, -limit:limit) ! matriz de kai

        integer :: i
        real(dp) :: suma


        if(rank.eq.0)print*,'chi calculation'

        suma = 0.0

        do ix = -limit, limit
            do iy = -limit, limit
                do iz = -limit, limit
                    matriz(ix, iy, iz) = 0.0
                enddo
            enddo
        enddo

        seed = 1010
          
    !      MCsteps = 10
        MCsteps = 100000000

        lseg=0.35
        l = lseg ! OJO!!!

        do i = 1, MCsteps
            x = 3.0*(rands(seed)-0.5)*delta ! random number between -1.5 * delta and 1.5 * delta
            y = 3.0*(rands(seed)-0.5)*delta 
            z = 3.0*(rands(seed)-0.5)*delta 

            u = ut(x, y)
            v = vt(x, y)

            radius = sqrt(x**2 + y**2 + z**2) ! real space
     
            if(radius<=(1.5*delta)) then  ! It is not inside the cut-off sphere
                if(radius>=l) then  ! is within the sphere of the segment

                    ! cell

                    ix = anint(v/delta)   ! space grid
                    iy = anint(u/delta) 
                    iz = anint(z/delta) 

                    matriz(ix, iy, iz) = matriz(ix, iy, iz) + (l/radius)**6
                endif
            endif
                    
        enddo

        do ix = -2, 2
            do iy = -2, 2
                do iz = -2, 2

                    Xu(ix, iy, iz) = matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)
                    suma = suma +  matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)

                enddo
            enddo
        enddo

        if(rank.eq.0)print*, 'Suma 5x5', suma

        suma = 0.0

        do ix = -limit, limit
            do iy = -limit, limit
                do iz = -limit, limit
                    suma = suma +  matriz(ix, iy, iz)/MCsteps*((3.0*delta)**3)
                enddo
            enddo
        enddo

        if(rank.eq.0)print*, 'Suma Total', suma
     
    end

  
end module VdW

      
