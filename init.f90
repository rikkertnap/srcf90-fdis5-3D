module initxvector
    
    use precision_definition, only : dp
    
    implicit none

    private 
    public  :: make_guess

contains

subroutine make_guess(x, xguess, isfirstguess, flagstored, xstored)
  
    use globals, only : neq,neqmax,sysflag,bcflag,LEFT,RIGHT
    use volume, only : nsurf 

    implicit none

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
    
    use globals, only : sysflag
    
    implicit none

    real(dp), intent(in) :: x(:)          ! iteration vector 
    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potential 
        
    if(sysflag=="elect") then 
        call init_guess_elect(x,xguess)    
    else if(sysflag=="dipolarstrong") then 
        call init_guess_dipolar(x,xguess)
    else if(sysflag=="dipolarweak") then 
        call init_guess_dipolar(x,xguess)    
    else if(sysflag=="electdouble") then 
        call init_guess_electdouble(x,xguess)  
    else if(sysflag=="electnopoly") then 
        call init_guess_electnopoly(x,xguess)
    else if(sysflag=="neutral") then 
        call init_guess_neutral(x,xguess)
    else if(sysflag=="electA") then 
        call init_guess_electA(x,xguess)
    else     
        print*,"Wrong value sysflag : ", sysflag
    endif

end subroutine init_guess

         

subroutine init_guess_dipolar(x, xguess)
      
    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
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
    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.0_dp
    enddo
  
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
     
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf 
                read(un_file(2),*)psisurfL(i)     ! surface potetial 
            enddo
        endif       

        do i=1,nsize
            read(un_file(1),*)xsol(i)       ! solvent
            read(un_file(2),*)psi(i)        ! potential
            x(i)         = xsol(i)       ! placing xsol  in vector x
            x(i+nsize)   = psi(i)        ! placing xsol  in vector x
        enddo

        if(bcflag(RIGHT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfR(i)     ! degree of complexation A 
            enddo
        endif    
        
        do i=1,2
            close(un_file(i))
        enddo

    endif  !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_dipolar

subroutine init_guess_electdouble(x, xguess)
      
    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
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
    do i=1,nsize
        x(i)=xbulk%sol
        x(i+  nsize)=0.0_dp
        x(i+2*nsize)=0.0_dp        ! 0.000001_dp
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
                print*, 'file number =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        if(bcflag(LEFT)/="cc") then 
            do i=1,nsurf 
                read(un_file(2),*)psisurfL(i)     ! surface potetial 
            enddo
        endif        
        do i=1,nsize
            read(un_file(1),*)xsol(i)       ! solvent
            read(un_file(2),*)psi(i)        ! degree of complexation A
            read(un_file(3),*)rhopolA(i)    ! degree of complexation A
            read(un_file(4),*)rhopolB(i)    ! degree of complexation A
            x(i)         = xsol(i)       ! placing xsol  in vector x
            x(i+nsize)   = psi(i)        ! placing xsol  in vector x
            x(i+2*nsize) = rhopolA(i)    ! placing xsol  in vector x
            x(i+3*nsize) = rhopolB(i)    ! placing xsol  in vector x
        enddo
        if(bcflag(RIGHT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfR(i)     ! degree of complexation A 
            enddo
        endif    
        
        do i=1,4
            close(un_file(i))
        enddo

    endif  !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo

end subroutine init_guess_electdouble

subroutine init_guess_electnopoly(x, xguess)
    
    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
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

    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.0_dp
    enddo

    neq_bc=0
    if(bcflag(LEFT)/="cc") then
        neq_bc=neq_bc+nsurf
        do i=1,nsurf 
            x(2*nsize+i)=0.0_dp
        enddo  
    endif
    if(bcflag(RIGHT)/="cc") then 
        do i=1,nsurf
            x(2*nsize+neq_bc+i)=0.0_dp
        enddo
        neq_bc=neq_bc+nsurf
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
     
        if(bcflag(LEFT)/="cc") then
            do i=1,nsurf
                read(un_file(2),*)psisurfL(i)     ! degree of complexation A
            enddo   
        endif     
        do i=1,nsize
            read(un_file(1),*)xsol(i)    ! solvent
            read(un_file(2),*)psi(i)     ! degree of complexation A
            x(i)      = xsol(i)     ! placing xsol  in vector x
            x(i+nsize) = psi(i)     ! placing xsol  in vector x
        enddo
        if(bcflag(RIGHT)/="cc") then 
            do i=1,nsurf
                read(un_file(2),*)psisurfR(i)     ! degree of complexation A 
            enddo
        endif        
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

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi,rhopolA,rhopolB
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
            read(un_file(3),*)rhopolA(i) ! degree of complexation A
            read(un_file(4),*)rhopolB(i)   ! degree of complexation A
            x(i)         = xsol(i)    ! placing xsol  in vector x
            x(i+nsize)   = psi(i)     ! placing xsol  in vector x
            x(i+2*nsize) = rhopolA(i) ! placing xsol  in vector x
            x(i+3*nsize) = rhopolB(i)   ! placing xsol  in vector x
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


subroutine init_guess_electA(x, xguess)

    use globals, only : neq,bcflag,LEFT,RIGHT,nsize
    use volume, only : nsurf
    use field, only : xsol,psi,rhopolA
    use surface, only : psisurfL, psisurfR 
    use parameters, only : xbulk, infile
    use myutils, only : newunit

    implicit none
  
    real(dp) :: x(neq)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(neq)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(3)
    integer :: ios,un_file(3)
  
    ! .. init guess all xbulk     

    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.0_dp
        x(i+2*nsize)=0.0_dp
    enddo
  
    if (infile.eq.1) then   ! infile is read in from file/stdio  
    
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'psi.in'
        write(fname(3),'(A7)')'rhoA.in'
        
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
            read(un_file(2),*)psi(i)     ! degree of complexation A
            read(un_file(3),*)rhopolA(i) ! degree of complexation A
            x(i)         = xsol(i)    ! placing xsol  in vector x
            x(i+nsize)   = psi(i)     ! placing xsol  in vector x
            x(i+2*nsize) = rhopolA(i) ! placing xsol  in vector x
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

end subroutine init_guess_electA


subroutine init_guess_neutral(x, xguess)
    
    use globals, only : neq,nsize
    use volume, only : nz
    use field, only : xsol,rhopolB 
    use parameters, only : xbulk, infile
    use myutils, only : newunit

    implicit none
  
    real(dp) :: x(:)       ! volume fraction solvent iteration vector 
    real(dp) :: xguess(:)  ! guess fraction  solvent 
  
    !     ..local variables 
    integer :: n, i
    character(len=8) :: fname(2)
    integer :: ios,un_file(2)
  
    !     .. init guess all xbulk      

    do i=1,nsize
        x(i)=xbulk%sol
        x(i+nsize)=0.000_dp
    enddo
  
  
    if (infile.eq.1) then   ! infile is read in from file/stdio  
        write(fname(1),'(A7)')'xsol.in'
        write(fname(2),'(A6)')'rhopolB.in'
     
        do i=1,2
            open(unit=newunit(un_file(i)),file=fname(i),iostat=ios,status='old')
            if(ios >0 ) then
                print*, 'file number =',un_file(i),' file name =',fname(i)
                print*, 'Error opening file : iostat =', ios
                stop
            endif
        enddo
     
        do i=1,nsize
            read(un_file(1),*)xsol(i)       ! solvent
            read(un_file(2),*)rhopolB(i)    ! density polymer B
            x(i)      = xsol(i)       ! placing xsol  in vector x
            x(i+nz)   = rhopolB(i)    ! placing rhopolB  in vector x
        enddo
        
        close(un_file(1))
        close(un_file(2))
    endif
    !     .. end init from file 
  
    do i=1,neq
        xguess(i)=x(i)
    enddo
    
end subroutine init_guess_neutral


! .. copy solution of previous solution ( distance ) to create new guess
! .. data x=(pi,psi) and pi and psi order and split into  blocks


subroutine make_guess_from_xstored(xguess,xstored)

    use globals, only : neq, neqmax, sysflag, bcflag, LEFT, RIGHT, nsize
    use volume, only : nz, nzstep, nsurf

    implicit none

    real(dp), intent(out) :: xguess(:)    ! guess volume fraction solvent and potentia
    real(dp), intent(in) :: xstored(:)

    !   .. local variables
    integer :: i, neq_bc

    neq_bc=0  
      
    if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+nsurf
    if(bcflag(LEFT)/="cc") neq_bc=neq_bc+nsurf 

    if(sysflag=="elect") then 
        
        do i=1,nsize
            xguess(i)=xstored(i)                                ! volume fraction solvent 
            xguess(i+  nsize)=xstored(i+   nsize+nsurf*nzstep)       ! potential
            xguess(i+2*nsize)=xstored(i+2*(nsize+nsurf*nzstep))  
            xguess(i+3*nsize)=xstored(i+3*(nsize+nsurf*nzstep))  
        enddo
        
        do i=1,neq_bc
            xguess(4*nsize+i)=xstored(4*(nsize+nsurf*nzstep)+i) 
        enddo

    elseif(sysflag=="electdouble") then 
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

    elseif (sysflag=="electnopoly".or.sysflag=="dipolarstrong".or.sysflag=="dipolarweak") then 
        do i=1,nsize/2
            xguess(i)=xstored(i)                    ! volume fraction solvent 
            xguess(i+nsize)=xstored(i+nsize+nsurf*nzstep)       ! potential
        enddo
        do i=nsize/2+1,nsize  ! shift by nzstep
            xguess(i)=xstored(i+nsurf*nzstep)    
            xguess(i+nsize)=xstored(i+nsize+2*nsurf*nzstep)      
        enddo       

        do i=1,neq_bc
            xguess(2*nsize+i)=xstored(2*(nsize+nsurf*nzstep)+i) 
        enddo   
    else
        print*,"Error : make_guess_from_xstored wrong sysflag"
        print*,"sysflag",sysflag
        stop
    endif    

end subroutine make_guess_from_xstored



end module initxvector
