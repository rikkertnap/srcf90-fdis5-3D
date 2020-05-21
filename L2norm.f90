
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

end module vectornorm
      
