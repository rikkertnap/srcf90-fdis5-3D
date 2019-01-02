!  module file of phyiscal constants
!  constants from NIST : CODATA values 2014

module physconst

    use precision_definition

    implicit none  
!     CODATA values

    real(dp), parameter :: Na          = 6.022140857e23_dp          ! Avogadro's number unit 
    real(dp), parameter :: kBoltzmann  = 1.38064852e-23_dp          ! Boltzmann constant unit J K^-1 
    real(dp), parameter :: elemcharge  = 1.6021766208e-19_dp        ! elementary charge unit C 
    real(dp), parameter :: dielect0    = 8.854187817e-12_dp         ! dielectric constant of vacuum unit C 
   
!     old values 
!    real(dp), parameter :: Na          = 6.0221417930e23_dp         ! Avogadro's number unit (Na=6.02d23)
!    real(dp), parameter :: kBoltzmann  = 1.3806504e-23_dp           ! Boltzmann constant unit J K^-1
!    real(dp), parameter :: elemcharge  = 1.60217648740e-19_dp       ! elementary charge unit C 
!    real(dp), parameter :: dielect0    = 8.854187817e-12_dp         ! dielectric constant of vacuum unit C 

end module physconst


