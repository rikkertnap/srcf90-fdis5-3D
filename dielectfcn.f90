! module dielectfcn.f90:                                                 | 
! computes the effectice dielectic function epsfcn and          |
! the derivate of the dielectric function  with respect to      |
! the polymer volume fraction                                   |
! both are in units of the dielectric constant of ionic solvent |
! pre : phi = 1-xsol -xpos -xneg xHplus- xOHmin                 |
! post: espfcn and Despfcn                                      |

module dielectric_const

    use precision_definition  
    implicit none
    
    abstract interface
        subroutine dielectfunct(phi,epsfcn,Depsfcn,dielectP,dielectW,n) 
            use precision_definition
            implicit none
    
            integer, intent(in)  :: n  
            real(dp), intent(in) ::  dielectP, dielectW 
            real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
            real(dp), intent(in) :: phi(:)

        end subroutine dielectfunct
    end interface

    procedure(dielectfunct), pointer :: dielectfcn => null()


    private                    ! default all routines in this module private 
    public  ::  dielectfcn,born,set_dielect_fcn

contains



subroutine set_dielect_fcn(dielect_env)

        character(len=15), intent(in) :: dielect_env

        select case (dielect_env)
            case ("constant") 
                dielectfcn=> dielectfcnConst 
            case ("linear") 
                dielectfcn => dielectfcnAV
            case ("MaxwellGarnett") 
                dielectfcn => dielectfcnMG
            case default
                print*,"Error in set_dielect_fcn"    
                print*,"Wrong value dielect_env : ",dielect_env
                stop
        end select  
    
end subroutine set_dielect_fcn



subroutine dielectfcnConst(phi,epsfcn,Depsfcn,dielectP,dielectW,n) 

    integer, intent(in)  :: n  
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
    real(dp), intent(in) :: phi(:)
    
    integer :: i  

    do i=1,n  
        epsfcn(i)  = 1.0d0 ! dielectric function
        Depsfcn(i) = 0.0d0 ! derivative dielectric function    
    enddo
                                
end subroutine

! volume fraction weighted average of dielectric constant

subroutine dielectfcnAV(phi,epsfcn,Depsfcn,dielectP,dielectW ,n) 

    integer, intent(in)  :: n
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp), intent(inout) :: epsfcn(:),Depsfcn(:)
    real(dp), intent(in) :: phi(:)
    

    integer :: i
    real(dp)  :: ratioeps

    ratioeps = dielectP/dielectW 
   
    do i=1,n  
        epsfcn(i)  = 1.0_dp-phi(i) + ratioeps * phi(i) ! dieletric function
        Depsfcn(i) = -1.0_dp+ratioeps ! derivative dieletric function    
    enddo
                                
end subroutine

! Maxwell-Garnett mixing rule for dielectic constant
! valid only for low phi<10-5 :

subroutine dielectfcnMG(phi,epsfcn,Depsfcn,dielectP, dielectW, n)  
   
    use mathconst

    integer, intent(in) :: n
    real(dp), intent(in) ::  dielectP, dielectW 
    real(dp) , intent(inout):: epsfcn(:),Depsfcn(:)
    real(dp), intent(in) :: phi(:)

    !     .. local variables

    real(dp) :: Kalpha, Kmossotti, gamma
    real(dp) :: epspol,epssol
    integer ::i
  
    epssol =dielectW          ! permittivty ionic solution 
    epspol =dielectP          ! permittivty polymer  
   
    Kmossotti = (epssol-epspol)/(2.0_dp*epssol+epspol)  ! Mossotti factor
    gamma=1.5_dp         ! ratio between electric radius and matter radius
    Kalpha= Kmossotti * 4.0_dp* pi /( 3.0_dp* gamma**3)
      
    do i=1,n  
        epsfcn(i)= 1.0_dp-3.0_dp*Kalpha*phi(i)/(1.0_dp+Kalpha*phi(i)) ! dieletric function
        Depsfcn(i)= -(3.0_dp*Kalpha/((1.0_dp+Kalpha*phi(i))**2))      ! derivative dieletric function
    enddo
                                
end subroutine

                                                
! computes the born energy for a given  Bjerrum lenght and size 
! U_B = (beta z^2e^2)/(8 pi eps eps0 a )= lB z /(2 *a)                       


function born(lB, radius, z) result(born_energy)

    real(dp), intent(in) :: lB, radius
    integer, intent(in) :: z
    real(dp) :: born_energy

    born_energy=(lB*z*z/(2.0_dp*radius))

end function
      


end module

