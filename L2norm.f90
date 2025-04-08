
module vectornorm

! computes L2norm
    use precision_definition
    implicit none

contains

    function l2norm (f,n)result(norm)
      
        implicit none

        integer, intent(in)  :: n 
        real(dp), intent(in) :: f(n)
        real(dp)             :: norm ! output

        integer :: i ! dummy index

        norm=0.0_dp
        do i=1,n
            norm = norm + f(i)*f(i)
        enddo
        norm=sqrt(norm)
        
    end function l2norm


    ! computes L2norm of vector f

    function l2norm_f90(f)result(norm)
      
        implicit none

        real(dp), intent(in) :: f(:)
        real(dp)             :: norm ! output
        integer              :: i ,n

        n=size(f)
        norm=0.0_dp
        do i=1,n
            norm = norm + f(i)**2
        enddo
        norm=sqrt(norm)
        
    end function l2norm_f90


end module vectornorm
      
