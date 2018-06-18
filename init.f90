module initxvector
    
    use precision_definition, only : dp
    
    implicit none

    private 
    public  :: make_guess

contains

subroutine make_guess(x, xguess, isfirstguess, flagstored, xstored)
  
    use globals, only : neq,neqmax,sysflag,bcflag,LEFT,RIGHT

    implicit none

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
    logical, intent(in) :: isfirstguess     ! first guess   
    logical, optional, intent(in) :: flagstored
    real(dp), optional, intent(in) :: xstored(:)

    !  ..local variables 
    integer :: i,neq_bc
  
    neq_bc=0    

    if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+1
    if(bcflag(LEFT)/="cc") neq_bc=neq_bc+1 

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
    
    use globals, only : sysflag
    
    implicit none

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
        
    if(sysflag=="elect") then 
        call init_guess_elect(x,xguess)
    else if(sysflag=="electdouble") then 
        call init_guess_electdouble(x,xguess)  
    else if(sysflag=="electnopoly") then 
        call init_guess_neutral(x,xguess)
    else if(sysflag=="neutral") then 
        call init_guess_electnopoly(x,xguess)
    else     
        print*,"Wrong value sysflag : ", sysflag
    endif

end subroutine init_guess



subroutine init_guess_electdouble(x, xguess)
      
    use globals, only : neq,bcflag,LEFT,RIGHT
    use volume, only : nz
    use field, only : xsol,psi,fdisA,rhopolA,rhopolB
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit

    implicit none
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
  
    ! .. init guess all xbulk     
    do i=1,nz
        x(i)=xbulk%sol
        x(i+nz)=0.000_dp
        x(i+2*nz)=0.000001_dp
        x(i+3*nz)=0.0_dp
    enddo
  
    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'rhoA.in'
        write(fname(4),'(A7)')'rhoB.in'
        
        do i=1,4 ! loop files
        
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file number =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag(LEFT)/="cc") read(un_file(2),*)psisurfL     ! degree of complexation A
        do i=1,nz
            read(un_file(1),*)xsol(i)       ! solvent
            read(un_file(2),*)psi(i)        ! degree of complexation A
            read(un_file(3),*)rhopolA(i)    ! degree of complexation A
            read(un_file(4),*)rhopolB(i)    ! degree of complexation A
            x(i)      = xsol(i)       ! placing xsol  in vector x
            x(i+nz)   = psi(i)        ! placing xsol  in vector x
            x(i+2*nz) = rhopolA(i)    ! placing xsol  in vector x
            x(i+3*nz) = rhopolB(i)    ! placing xsol  in vector x
        enddo
        if(bcflag(RIGHT)/="cc") read(un_file(2),*)psisurfR     ! degree of complexation A 
       
         do i=1,4
            close(un_file(i))
        enddo

    endif  !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_electdouble

subroutine init_guess_electnopoly(x, xguess)
    
    use globals, only : neq,bcflag,LEFT,RIGHT
    use volume, only : nz
    use field, only : xsol,psi
    use surface, only : psisurfL, psisurfR   
    use parameters, only : xbulk ,infile
    use myutils, only : newunit

    implicit none
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(2)
    integer :: ios,un_file(4)
    integer :: neq_bc 
  
  
    ! .. init guess all xbulk     

    do i=1,nz
        x(i)=xbulk%sol
        x(i+nz)=0.000_dp
    enddo

    neq_bc=0
    if(bcflag(LEFT)/="cc") then
        neq_bc=neq_bc+1 
        x(2*nz+neq_bc)=0.00_dp
    endif      
    if(bcflag(RIGHT)/="cc") then 
        neq_bc=neq_bc+1 
        x(2*nz+neq_bc)=0.00_dp
    endif     

    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
     
        do i=1,2 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file number =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag(LEFT)/="cc") read(200,*)psisurfL     ! degree of complexation A
        do i=1,nz
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! degree of complexation A
            x(i)      = xsol(i)    ! placing xsol  in vector x
            x(i+nz)   = psi(i)     ! placing xsol  in vector x
        enddo
        if(bcflag(RIGHT)/="cc") read(un_file(2),*)psisurfR     ! degree of complexation A 
       
        do i=1,2
            close(un_file(i))
        enddo

    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_electnopoly
!     purpose: initalize x and xguess

subroutine init_guess_elect(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT
    use volume, only : nz
    use field, only : xsol,psi,fdisA,xpolC
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit

    implicit none
  
    real(dp) :: x(neq)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(neq)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(4)
    integer :: ios,un_file(4)
  
  
    ! .. init guess all xbulk     

    do i=1,nz
        x(i)=xbulk%sol
        x(i+nz)=0.000_dp
        x(i+2*nz)=0.000001_dp
        x(i+3*nz)=0.00001_dp
!        x(i+4*nz)=0.00_dp
    enddo
  
    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A5)')'gA.in'
        write(fname(4),'(A5)')'xpolC.in'
     
        do i=1,4 ! loop files
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then    
                print*, 'file num ber =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
        if(bcflag(LEFT)/="cc") read(200,*)psisurfL    
        do i=1,nz
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! degree of complexation A
            read(un_file(3),*)fdisA(5,i) ! degree of complexation A
            read(un_file(4),*)xpolC(i)   ! degree of complexation A
            x(i)      = xsol(i)    ! placing xsol  in vector x
            x(i+nz)   = psi(i)     ! placing xsol  in vector x
            x(i+2*nz) = fdisA(5,i) ! placing xsol  in vector x
            x(i+3*nz) = xpolC(i)   ! placing xsol  in vector x
        enddo
    
        if(bcflag(RIGHT)/="cc") read(200,*)psisurfR    
       
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
    
    use globals, only : neq
    use volume, only : nz
    use field, only : xsol,rhopolB 
    use parameters, only : xbulk, infile

    implicit none
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(2)
    integer :: ios,nfile(2)
  
    !     .. init guess all xbulk      

    do i=1,nz
        x(i)=xbulk%sol
        x(i+nz)=0.000_dp
    enddo
  
  
    if (infile.eq.1) then   ! infile is read in from file/stdio  
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'rhopolB.in'
     
        nfile(1)=100
        nfile(2)=200
     
        do i=1,2
            open(unit=nfile(i),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                print*, 'file number =',nfile(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        do i=1,nz
            read(100,*)xsol(i)       ! solvent
            read(200,*)rhopolB(i)    ! density polymer B
            x(i)      = xsol(i)       ! placing xsol  in vector x
            x(i+nz)   = rhopolB(i)    ! placing rhopolB  in vector x
        enddo
        
        close(100)
        close(200)
    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutral


!     .. copy solution of previous solution ( distance ) to create new guess
!     .. data x=(pi,psi) and pi and psi order and split into  blocks


subroutine make_guess_from_xstored(xguess,xstored)

    use globals, only : neq, neqmax, sysflag, bcflag, LEFT, RIGHT
    use volume, only : nz, nzstep

    implicit none

    real(dp), intent(out) :: xguess(neq)    ! guess volume fraction solvent and potentia
    real(dp), intent(in) :: xstored(neqmax)

    !   .. local variables
    integer :: i,neq_bc

    neq_bc=0    
    if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+1
    if(bcflag(LEFT)/="cc") neq_bc=neq_bc+1 

    if(sysflag=="elect".or.sysflag=="electdouble") then 
        do i=1,nz/2
            xguess(i)=xstored(i)                    ! volume fraction solvent 
            xguess(i+nz)=xstored(i+nz+nzstep)       ! potential
            xguess(i+2*nz)=xstored(i+2*(nz+nzstep))  
            xguess(i+3*nz)=xstored(i+3*(nz+nzstep))  
        enddo
        do i=nz/2+1,nz  ! shift by nzstep
            xguess(i)=xstored(i+nzstep)    
            xguess(i+nz)=xstored(i+nz+2*nzstep)      
            xguess(i+2*nz)=xstored(i+2*nz+3*nzstep)  
            xguess(i+3*nz)=xstored(i+3*nz+4*nzstep)  
        enddo       

        do i=1,neq_bc
            xguess(4*nz+i)=xstored(4*(nz+nzstep)+i) 
        enddo   
    elseif (sysflag=="electnopoly") then 
        do i=1,nz/2
            xguess(i)=xstored(i)                    ! volume fraction solvent 
            xguess(i+nz)=xstored(i+nz+nzstep)       ! potential
        enddo
        do i=nz/2+1,nz  ! shift by nzstep
            xguess(i)=xstored(i+nzstep)    
            xguess(i+nz)=xstored(i+nz+2*nzstep)      
        enddo       

        do i=1,neq_bc
            xguess(2*nz+i)=xstored(2*(nz+nzstep)+i) 
        enddo   
    else
        print*,"Error : make_guess_from_xstored wrong sysflag"
        print*,"sysflag",sysflag
        stop
    endif    

end subroutine make_guess_from_xstored



end module initxvector
