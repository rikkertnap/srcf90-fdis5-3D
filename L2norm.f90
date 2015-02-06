module vectornorm
! computes L2norm
 
    implicit none

    contains

    real*8 function l2norm (f,n)
      
        implicit none

        integer, intent(in)  :: n 
        real*8, intent(in) :: f(n)
  
        integer :: i
        real*8 :: norm

        norm=0.0d0

        do i=1,n
            norm = norm + f(i)*f(i)
        enddo
  
        l2norm=dsqrt(norm)
        return

    end function l2norm

end module vectornorm      
