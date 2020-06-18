module initxvector
    
    use precision_definition, only : dp
    
    implicit none

    private 
    public  :: make_guess

contains

subroutine make_guess(x, xguess, isfirstguess, flagstored, xstored)
  
    use globals, only : neq,neqmax,systype,bcflag,LEFT,RIGHT
    use volume, only : nsurf 

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
    logical, intent(in) :: isfirstguess     ! first guess   
    logical, optional, intent(in) :: flagstored
    real(dp), optional, intent(in) :: xstored(:)
    
    !  ..local variables 
    integer :: i ,neq_bc

    if(present(flagstored)) then
        if(present(xstored)) then
            if(flagstored) then  
                call make_guess_from_xstored(xguess,xstored)
            else if(isfirstguess) then       ! first guess
                call init_guess(x,xguess)
            else  
                do i=1,neq
                    xguess(i)=x(i)      
                enddo
            endif
        else
            print*,"Error: argument xstored not present, while flagstored present"
            stop 
        endif 
    else if(isfirstguess) then       ! first guess
        call init_guess(x,xguess)
    else     
        do i=1,neq
            xguess(i)=x(i)     
        enddo
    endif

end subroutine make_guess



subroutine init_guess(x, xguess)
    
    use globals, only : systype

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
       
    select case (systype)
        case ("elect")   
            call init_guess_elect(x,xguess)    
        case ("neutral")  
            call init_guess_neutral(x,xguess)
        case ("neutralnoVdW")  
            call init_guess_neutralnoVdW(x,xguess)
        case ("brush_mul")  
            call init_guess_multi(x,xguess)
        case ("brushssdna")  
            call init_guess_multi(x,xguess)
        case ("brushborn") 
            print*,"not implemented yet."  
        case default   
            print*,"Init_guess: Wrong value systype : ", systype
    end select 

end subroutine init_guess

         




!     purpose: initalize x and xguess

subroutine init_guess_elect(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
    integer, parameter :: A=1, B=2   
  
    ! .. init guess all xbulk     

    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.0_dp
        x(i+2*nsize)=0.0_dp
        x(i+3*nsize)=0.0_dp
    enddo
   
    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'rhoA.in'
        write(fname(4),'(A7)')'rhoB.in'
     
        do i=1,4 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif            
        do i=1,nsize
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! degree of complexation A
            read(un_file(3),*)rhopol(i,A) ! degree of complexation A
            read(un_file(4),*)rhopol(i,B)   ! degree of complexation A
            x(i)         = xsol(i)    ! placing xsol  in vector x
            x(i+nsize)   = psi(i)     ! placing xsol  in vector x
            x(i+2*nsize) = rhopol(i,A) ! placing xsol  in vector x
            x(i+3*nsize) = rhopol(i,B)   ! placing xsol  in vector x
        enddo
    
        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
        do i=1,4
            close(un_file(i))
        enddo

    endif

    !     .. end init from file 
    do i=1,neq
        xguess(i)=x(i)
    enddo

    
end subroutine init_guess_elect


subroutine init_guess_neutral(x, xguess)
    
    use globals, only : neqint,nsize,nsegtypes
    use volume, only : nz
    use field, only : xsol,rhopol
    use parameters, only : xbulk, infile
    use myutils, only : newunit
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, t
    character(len=8) :: fname(2)
    integer :: ios,un_file(2)
  
    !     .. init guess all xbulk      

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo


    if (infile.eq.1) then   ! infile is read in from file/stdio  
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'rho.in'
        do i=1,2 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                print*, 'file number =',un_file,' file name =',fname
                print*, 'Error opening file : iostat =', ios                
                stop
            endif
        enddo
        do i=1,nsize
            read(un_file(1),*)xsol(i) ! solvent
            read(un_file(2),*)(rhopol(i,t),t=1,nsegtypes)
            x(i) = xsol(i)            ! placing xsol  in vector x
            do t=1,nsegtypes          ! placing rhopol(:,t) in vector x
                x(i+t*nsize)=rhopol(i,t)
            enddo      
        enddo     
        close(un_file(1))
        close(un_file(2))
    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutral


subroutine init_guess_neutralnoVdW(x, xguess)
    
    use globals, only : neqint,nsize,nsegtypes
    use volume, only : nz
    use field, only : xsol,rhopol
    use parameters, only : xbulk, infile
    use myutils, only : newunit
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, t
    character(len=8) :: fname
    integer :: ios,un_file
  
    !     .. init guess all xbulk      

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo


    if (infile.eq.1) then   ! infile is read in from file/stdio  
        write(fname,'(A7)')'xsol.in'
        open(unit=newunit(un_file),file=fname,iostat=ios,status='old')
        if(ios >0 ) then
            print*, 'file number =',un_file,' file name =',fname
            print*, 'Error opening file : iostat =', ios                
            stop
        endif
        do i=1,nsize
            read(un_file,*)xsol(i) ! solvent
            x(i) = xsol(i)            ! placing xsol  in vector x   
        enddo     
        close(un_file)
    endif
    !     .. end init from file 
  
    do i=1,neqint
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutralnoVdW


subroutine init_guess_multi(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize,neqint,nsegtypes
    use volume, only : nsurf
    use field, only : xsol,psi,rhopol
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i, t
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
  
    ! .. init guess all xbulk     

    do i=1,neqint
        x(i)=0.0_dp    
    enddo

    do i=1,nsize
        x(i)=xbulk%sol
    enddo

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A6)')'rho.in'
     
        do i=1,3 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)
            enddo
        endif            
        do i=1,nsize
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! potential
            read(un_file(3),*)(rhopol(i,t),t=1,nsegtypes)

            x(i)         = xsol(i)    ! placing xsol in vector x
            x(i+nsize)   = psi(i)     ! placing psi in vector x
            do t=1,nsegtypes          ! placing rhopol(:,t) in vector x
                x(i+(t+1)*nsize)=rhopol(i,t)
            enddo        
        enddo
    
        if(bcflag(RIGHT)/="cc") then
            do i=1,nsurf 
                read(un_file(2),*)psisurfR(i)
            enddo
        endif            
       
         do i=1,3
            close(un_file(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_multi



! .. copy solution of previous solution ( distance ) to create new guess
! .. data x=(pi,psi) and pi and psi order and split into  blocks


subroutine make_guess_from_xstored(xguess,xstored)

    use globals, only : neq, neqmax, systype, bcflag, LEFT, RIGHT, nsize
    use volume, only : nz, nzstep, nsurf

    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potentia
    real(dp), intent(in) :: xstored(:)

    !   .. local variables
    integer :: i, neq_bc

    neq_bc=0  
      
    if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+nsurf
    if(bcflag(LEFT)/="cc") neq_bc=neq_bc+nsurf 

    if(systype=="elect") then 
        
        do i=1,nsize
            xguess(i)=xstored(i)                                ! volume fraction solvent 
            xguess(i+  nsize)=xstored(i+   nsize+nsurf*nzstep)       ! potential
            xguess(i+2*nsize)=xstored(i+2*(nsize+nsurf*nzstep))  
            xguess(i+3*nsize)=xstored(i+3*(nsize+nsurf*nzstep))  
        enddo
        
        do i=1,neq_bc
            xguess(4*nsize+i)=xstored(4*(nsize+nsurf*nzstep)+i) 
        enddo

    elseif(systype=="electdouble") then 
        ! assume xstored xsol,psi symetric xsol(1)=xsol(nsize) etc     
        do i=1,nsize/2
            xguess(i)=xstored(i)                                ! volume fraction solvent 
            xguess(i+nsize)=xstored(i+nsize+nsurf*nzstep)       ! potential
            xguess(i+2*nsize)=xstored(i+2*(nsize+nsurf*nzstep))  
            xguess(i+3*nsize)=xstored(i+3*(nsize+nsurf*nzstep))  
        enddo
        do i=nsize/2+1,nsize  ! shift by nzstep
            xguess(i)=xstored(i+nsurf*nzstep)    
            xguess(i+  nsize)=xstored(i+  nsize+2*nsurf*nzstep)      
            xguess(i+2*nsize)=xstored(i+2*nsize+3*nsurf*nzstep)  
            xguess(i+3*nsize)=xstored(i+3*nsize+4*nsurf*nzstep)  
        enddo       

        do i=1,neq_bc
            xguess(4*nsize+i)=xstored(4*(nsize+nsurf*nzstep)+i) 
        enddo
    else
        print*,"Error : make_guess_from_xstored wrong systype"
        print*,"systype",systype
        stop
    endif    

end subroutine make_guess_from_xstored



end module initxvector
