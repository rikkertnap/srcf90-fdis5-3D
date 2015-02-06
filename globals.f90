!     .. module file of global variables 

module globals
  
    use mathconst
  
    implicit none
  
    !     .. variables

    integer  :: nsize         ! size lattice, numer of layers
    integer*8  :: neq         ! number of non-linear equations
    integer*8  :: neqmax      ! maximum number of non-linear equations
    integer  :: nsegAB        ! length of AB polymer 
    integer  :: nsegC         ! length of C polymer 
    integer  :: cuantasAB     ! number of polymer configurations
    integer  :: cuantasC      ! number of polymer configurations
    integer  :: max_conforAB  ! maximum number of polymer configurations
    integer  :: max_conforC   ! maximum number of polymer configurations
    character(len=15) :: sysflag   ! sysflag selects fcn    
    character(len=2)  :: bcflag(2) ! bcflag selects bc surface 

    integer, parameter :: LEFT = 1
    integer, parameter :: RIGHT = 2

end module globals

