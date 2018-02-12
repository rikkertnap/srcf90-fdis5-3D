module loopvar

    use precision_definition
    implicit none
    
    type looplist
       real(dp) :: val
       real(dp) :: min
       real(dp) :: max
       real(dp) :: stepsize
       real(dp) :: delta
    end type looplist
    
  contains

end module
