! --------------------------------------------------------------|
! fcnMg-excpl.f90:                                              |
! Constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for                |
! polyelctolyte  with ionbinding and  Mg binding                |
! using phosphate pairs.                                        |
! --------------------------------------------------------------|

module modfcnMgexpl
   
    use mpivars
    use precision_definition

    implicit none


contains
   
     subroutine compute_fdisPP(fdisPP,fdisP2Mg,position1,position2)


        real(dp), intent(inout), dimension(:,:) :: fdisPP
        real(dp), intent(inout) :: fdisP2Mg
        integer , intent(in) :: position1, position2
                        
    end subroutine compute_fdisPP


    ! polylectrolyte with phosphate acid groups
    ! with ion chargeable group being on one acid (tA) with counterion binding etc 
    ! distribute volume of neighboring cells

    subroutine fcn_Mg_expl(x,f,nn)

        use mpivars
        use precision_definition
        use globals, only : neq 

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        

    end subroutine fcn_Mg_expl



    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_average_charge_PP_expl(avfdisP2Mg,avfdisPP)

        use mpivars
        use precision_definition

        real(dp), intent(inout) :: avfdisP2Mg
        real(dp), intent(inout) :: avfdisPP(5,5)

        !     .. local variables


    end subroutine compute_average_charge_PP_expl


    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_FEchem_react_PP_expl(FEchemPP)

        use mpivars
        use precision_definition

        real(dp), intent(inout) :: FEchemPP

    end subroutine compute_FEchem_react_PP_expl

   

end module modfcnMgexpl

   

