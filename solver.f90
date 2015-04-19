subroutine solver(x, xguess, accuracy, residual)
  
    use globals
    use parameters
    use listfcn
    use fcnpointer

    implicit none
  
    !     .. variables and constant declaractions

    !     .. array arguments
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)  
    !     .. scalar arguments
    real(dp) :: accuracy
    real(dp) :: residual
    character(len=80) :: fcnname
   
    write(fcnname,'(A6)')'solver'

    call set_size_neq()
    call check_value_sysflag(fcnname)
    
    select case (sysflag)
        case ("elect") 
            fcnptr => fcnelect
        case ("electdouble")  
            fcnptr => fcnelectdouble
        case ("electnopoly") 
            fcnptr => fcnelectNoPoly 
        case ("electHC") 
            fcnptr => fcnelectHC
        case ("neutral") 
            fcnptr => fcnneutral
        case ("bulk water") 
             fcnptr => fcnbulk
        case default
            print*,"Error in call to solver subroutine"    
            print*,"Wrong value sysflag : ", sysflag
            stop
    end select  


    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, neq, accuracy, residual)
     
     !  elseif (method.eq."zspow") then 
     !     call zspow_solver(x, xguess, neq, accuracy, residual)
     
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
