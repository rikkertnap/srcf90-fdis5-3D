!     module of mathconst

module precision_definition 
   implicit none
   integer,  parameter :: dp=kind(1.0D0)
!  integer , parameter :: sp=selected_real(6,37)
!  integer , parameter :: dp=selected_real(15,307)
!  integer , parameter :: qp=selected_real(33,4931)

end module precision_definition


module mathconst

  use precision_definition  
  implicit none
  
  real(dp) :: pi
  
end module mathconst
