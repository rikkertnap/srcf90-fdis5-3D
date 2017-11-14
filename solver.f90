subroutine solver(x, xguess, accuracy, residual)
  
    use globals
    use parameters
    use listfcn
    use fcnpointer

    implicit none
  
    !     .. arguments
    real(dp) :: x(neq)
    real(dp) :: xguess(neq)  
    real(dp) :: accuracy
    real(dp) :: residual


!    character(len=80) :: fcnname
   
!    write(fcnname,'(A6)')'solver'
!    call set_size_neq()
!    call check_value_sysflag(fcnname)
    

    call set_size_neq  
    call set_fcn

    if(method.eq."kinsol") then
     
        call kinsol_gmres_solver(x, xguess, neq, accuracy, residual)
     
    !  elseif (method.eq."zspow") then 
    !     call zspow_solver(x, xguess, neq, accuracy, residual)
     
    else  
        print*,"Solver method incorrect"
        stop
    endif
  
end subroutine solver
