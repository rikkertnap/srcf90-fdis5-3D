!  module file of phyiscal constants
module physconst

    use precision_definition

    implicit none  

    real(dp), parameter :: Na          = 6.0221417930d23         ! Avogadro's number unit (Na=6.02d23)
    real(dp), parameter :: kBoltzmann  = 1.3806504d-23           ! Boltzmann constant unit J K^-1
    real(dp), parameter :: elemcharge  = 1.60217648740d-19       ! elementary charge unit C 
    real(dp), parameter :: dielect0    = 8.854187817d-12         ! dielectric constant of vacuum unit C 
   
end module physconst
