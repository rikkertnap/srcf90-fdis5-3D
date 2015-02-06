subroutine solver(x, xguess, accuracy, residual)
  
    use globals
    use parameters
    use listfcn
    use fcnpointer

    implicit none
  
    !     .. variables and constant declaractions

    !     .. array arguments
    real*8 :: x(neq)
    real*8 :: xguess(neq)  
    !     .. scalar arguments
    real*8 :: accuracy
    real*8 :: residual
 
    call set_size_neq()
 
    if(sysflag=="electdouble") then 
        fcnptr => fcnelectdouble
    elseif(sysflag=="elect") then 
        fcnptr => fcnelect
    elseif(sysflag=="electHC") then
        fcnptr => fcnelectHC  
    elseif(sysflag=="bulk water") then 
        fcnptr => fcnbulk
    elseif(sysflag=="neutral") then 
        fcnptr => fcnneutral
    else
        print*,"Wrong value sysflag : ", sysflag
        stop
    endif

    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, neq, accuracy, residual)
     
     !  elseif (method.eq."zspow") then 
     !     call zspow_solver(x, xguess, neq, accuracy, residual)
     
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
