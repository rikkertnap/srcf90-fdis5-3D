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
  
end module VdW

      
