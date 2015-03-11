! --------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes on a coated onto a spherical surface             |
! --------------------------------------------------------------|


module fcnpointer

    implicit none

    abstract interface
        subroutine fcn(x,f,n)
            implicit none
            real*8, dimension(n), intent(in) :: x
            real*8, dimension(n), intent(out) :: f
            integer*8, intent(in) :: n
        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer


module listfcn

    implicit none

    contains

    ! pre:  vector x=xsol+psi+rhopolA+rhoolB+xpolC                   
    ! post: vector f=xpolAB+xpolC+xsol+Sum_i x_i -1,poisson equation,rhoA,rhoB,xpolC              

    subroutine fcnelectHC(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface
        use vectornorm
 
        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out) :: f(neq)
        integer*8, intent(in) :: nn

        !     .. declare local variables

        real*8 :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real*8 :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real*8 :: xA(3),xB(3),sumxA,sumxB
        real*8 :: constA,constB
        real*8 :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real*8 :: tmp,expVdW 
        real*8 :: norm
        integer :: conf              ! counts number of conformations
        real*8 :: cn                 ! auxilary variable for Poisson Eq


        real*8, parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nz                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
            xpolCin(i)=x(i+4*n)    ! volume fraction C-polymer
        enddo

        psiSurfL = psi(1)          ! surface potentail

        do i=1,n                  ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolA(i) = 0.0d0     ! A polymer density 
            rhopolB(i) = 0.0d0
            rhopolC(i) = 0.0d0
        
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))     ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))  ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))  ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0d0*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0d0/(1.0d0+sumxA)
                fdisA(5,i)=0.0d0
            else
                fdisA(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constA/((sumxA+1.0d0)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0d0)/(2.0d0*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0d0*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0d0/(1.0d0+sumxB)
                fdisB(5,i)=0.0d0
            else
                fdisB(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constB/((sumxB+1.0d0)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0d0)/(2.0d0*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

    !        exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
            ! Na condensed ANa reference state
            exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

       
            !     .. VdW interaction   
            tmp = 0.0d0
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*xpolCin(j)
                enddo
            endif
            expVdW=dexp(-VdWepsC*tmp)
            exppiC(i)=(xsol(i)**vpolC)*expVdW ! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0d0                 ! init q   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0d0)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0d0)
       
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 

        sigmaqSurfL=0.0d0 ! charge regulating surface charge 
        psi(n+1)=0.0d0   ! bulk potential

        !    .. Poisson Eq 

        f(n+1)= -0.5d0*((psi(2)-psi(1)) + sigmaqSurfL +rhoq(1)*constqW)      !     boundary

        do i=2,n
            f(n+i)= -0.5d0*(psi(i+1)-2.0d0*psi(i) + psi(i-1) +rhoq(i)*constqW)
        enddo

        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
            f(4*n+i)=xpolC(i)-xpolCin(i)
        enddo

        norm=l2norm(f,5*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelectHC

    subroutine fcnelectNoPoly(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out) :: f(neq)
        integer*8, intent(in) :: nn


        !     .. declare local variables

        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real*8 :: norm
        integer :: conf              ! counts number of conformations
!        real*8 :: cn                 ! auxilary variable for Poisson Eq
        integer :: neq_bc           

        !     .. executable statements 

        n=nz                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
        enddo
        
        neq_bc=0

        if(bcflag(RIGHT)/="cc") then
            neq_bc=neq_bc+1 
            psiSurfR =x(2*n+neq_bc)          ! surface potentail
        endif   
        if(bcflag(LEFT)/="cc") then 
            neq_bc=neq_bc+1
            psiSurfL =x(2*n+neq_bc)          ! surface potentail
        endif    

        neq_bc=0
    
        do i=1,n                  ! init volume fractions 
        
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
    
        enddo

        !   .. construction of fcn 
        

        do i=1,n

            f(i)=xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)
       
            !   ..  total charge density in units of vsol
        enddo 

        ! .. electrostatics 
        ! .. charge regulating surface charge 
        sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
        sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

        ! .. Poisson Eq 
        
        if(n/=1) then 
            f(n+1)= -0.5d0*( psi(2)-psi(1)  + sigmaqSurfL +rhoq(1)*constqW )      !     boundary
            f(2*n)= -0.5d0*( sigmaqSurfR- (psi(n)-psi(n-1)) +rhoq(n)*constqW )
            do i=2,n-1
                f(n+i)= -0.5d0*( psi(i+1)-2.0d0*psi(i) + psi(i-1) +rhoq(i)*constqW) 
            enddo
        else
             f(2)= -0.5d0*( sigmaqSurfR  + sigmaqSurfL +rhoq(1)*constqW )
        endif       

        ! self consistent boundary conditions
        neq_bc=0
        if(bcflag(RIGHT)/='cc') then 
            neq_bc=neq_bc+1
            f(2*n+neq_bc)=psisurfR-psi(n)-sigmaqSurfR/2.0d0
        else
            psisurfR=psi(n)+sigmaqSurfR/2.0d0
        endif   

        if(bcflag(LEFT)/='cc') then 
            neq_bc=neq_bc+1
            f(2*n+neq_bc)=psi(1)-psisurfL+sigmaqSurfL/2.0d0
        else    
            psisurfL=psi(1)+sigmaqSurfL/2.0d0
        endif   

        norm=l2norm(f,2*n+neq_bc)
        
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelectNoPoly


    subroutine fcnelect(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out) :: f(neq)
        integer*8, intent(in) :: nn


        !     .. declare local variables

        real*8 :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real*8 :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real*8 :: xA(3),xB(3),sumxA,sumxB
        real*8 :: constA,constB
        real*8 :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real*8 :: tmp,expVdW 
        real*8 :: norm
        integer :: conf              ! counts number of conformations
        real*8 :: cn                 ! auxilary variable for Poisson Eq
        integer :: neq_bc           

        real*8, parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nz                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
        neq_bc=0

        if(bcflag(RIGHT)/="cc") then
            neq_bc=neq_bc+1 
            psiSurfR =x(4*n+neq_bc)          ! surface potentail
        endif   
        if(bcflag(LEFT)/="cc") then 
            neq_bc=neq_bc+1
            psiSurfL =x(4*n+neq_bc)          ! surface potentail
        endif    
        neq_bc=0
    
        do i=1,n                  ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolA(i) = 0.0d0     ! A polymer density 
            rhopolB(i) = 0.0d0

            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))     ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))  ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))  ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0d0*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0d0/(1.0d0+sumxA)
                fdisA(5,i)=0.0d0
            else
                fdisA(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constA/((sumxA+1.0d0)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0d0)/(2.0d0*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0d0*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0d0/(1.0d0+sumxB)
                fdisB(5,i)=0.0d0
            else
                fdisB(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constB/((sumxB+1.0d0)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0d0)/(2.0d0*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
            ! Na condensed ANa reference state
    !        exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

       
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo


        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
       

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
           
       
            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0d0)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0d0)
       
            f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 
        ! .. charge regulating surface charge 
        sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
        sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

        ! .. Poisson Eq 

        f(n+1)= -0.5d0*((psi(2)-psi(1)) + sigmaqSurfL +rhoq(1)*constqW)      !     boundary
        f(2*n)= -0.5d0*(sigmaqSurfR- (psi(n)-psi(n-1)) +rhoq(n)*constqW)
        do i=2,n-1
            f(n+i)= -0.5d0*(psi(i+1)-2.0d0*psi(i) + psi(i-1) +rhoq(i)*constqW)
        enddo

        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
        enddo

        ! self consistent boundary conditions
        neq_bc=0
        if(bcflag(RIGHT)/='cc') then 
            neq_bc=neq_bc+1
            f(4*n+neq_bc)=psisurfR-psi(n)-sigmaqSurfR/2.0d0
        else
            psisurfR=psi(n)+sigmaqSurfR/2.0d0
!            print*,'psisurfR=',psisurfR,' sigmaqSurfR=',sigmaqSurfR
        endif   

        if(bcflag(LEFT)/='cc') then 
            neq_bc=neq_bc+1
            f(4*n+neq_bc)=psi(1)-psisurfL+sigmaqSurfL/2.0d0
        else    
            psisurfL=psi(1)+sigmaqSurfL/2.0d0
!            print*,'psisurfL=',psisurfL,' sigmaqSurfL=',sigmaqSurfL
        endif   

        norm=l2norm(f,4*n+neq_bc)
        
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelect

    subroutine fcnelectdouble(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out) :: f(neq)
        integer*8, intent(in) :: nn

        !     .. declare local variables

        real*8 :: exppiA(nsize),exppiB(nsize)  ! auxilairy variable for computing P(\alpha) 
        real*8 :: rhopolAin(nsize),rhopolBin(nsize)
        real*8 :: xA(3),xB(3),sumxA,sumxB
        real*8 :: constA,constB
        real*8 :: proL,rhopolABL0,proR,rhopolABR0
        integer :: n                 ! half of n
        integer :: i,j,k,kL,kR,c,s         ! dummy indices
        real*8 :: norm
        integer :: conf              ! counts number of conformations
        real*8 :: cn                 ! auxilary variable for Poisson Eq

        real*8, parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 

        !  .. executable statements 

        n=nz                        ! size vector neq=5*nz

        do i=1,n                    ! init x 
            xsol(i)= x(i)           ! solvent volume fraction 
            psi(i) = x(i+n)         ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo

        psiSurfR = psi(1)          ! surface potentail
        psiSurfL = psi(n)
   
        do i=1,n                   ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            rhopolAL(i) = 0.0d0     ! A polymer density 
            rhopolBL(i) = 0.0d0     ! B polymer density  
            rhopolAR(i) = 0.0d0     ! A polymer density 
            rhopolBR(i) = 0.0d0     ! B polymer density 

            xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))      ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))   ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))   ! ACa+/A-
       
            sumxA=xA(1)+xA(2)+xA(3)
            constA=(2.0d0*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0d0/(1.0d0+sumxA)
                fdisA(5,i)=0.0d0
            else
                fdisA(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constA/((sumxA+1.0d0)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0d0)/(2.0d0*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                       ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                       ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                       ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))     ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))   ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))   ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0d0*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0d0/(1.0d0+sumxB)
                fdisB(5,i)=0.0d0
            else
                fdisB(1,i)= (-1.0d0+dsqrt(1.0d0+4.0d0*constB/((sumxB+1.0d0)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0d0)/(2.0d0*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       
            ! A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

    !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
    !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
            ! Na condensed ANa reference state
    !        exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
    !        exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

        enddo

        !  .. computation polymer volume fraction 

        qABL = 0.0d0                  ! init q
        qABR = 0.0d0

        do c=1,cuantasAB              ! loop over cuantas
            proL=1.0d0                ! initial weight conformation 
            proR=1.0d0
            do s=1,nsegAB             ! loop over segments 
                kL=indexchainAB(c,s)
                kR=nz+1-kL
                if(isAmonomer(s)) then ! A segment 
                    proL = proL*exppiA(kL)
                    proR = proR*exppiA(kR)
                else
                    proL = proL*exppiB(kL)
                    proR = proR*exppiB(kR)
                endif
            enddo

            qABL = qABL+proL
            qABR = qABR+proR

            do s=1,nsegAB
                kL=indexchainAB(c,s)
                kR=nz+1-kL
                if(isAmonomer(s)) then ! A segment 
                    rhopolAL(kL)=rhopolAL(kL)+proL
                    rhopolAR(kR)=rhopolAR(kR)+proR
                else
                    rhopolBL(kL)=rhopolBL(kL)+proL
                    rhopolBR(kR)=rhopolBR(kR)+proR
                endif
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolABL0=sigmaABL/qABL
        rhopolABR0=sigmaABR/qABR
!        print*,"sigmaABR=",sigmaABR,"sigmaABL=",sigmaABL
!        print*,"qABL=",qABL,"qABR=",qABR

        do i=1,n

            rhopolA(i)= (rhopolABL0*rhopolAL(i)+rhopolABR0*rhopolAR(i))/deltaG(i)
            rhopolB(i)= (rhopolABL0*rhopolBL(i)+rhopolABR0*rhopolBR(i))/deltaG(i)
            
!            print*,i,rhopolB(i),rhopolBL(i),rhopolBR(i)

            do k=1,4               ! polymer volume fraction
                xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                    +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo    
            xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0d0)
            xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0d0)
       
       
            f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0d0
       
            rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
            !   ..  total charge density in units of vsol
        enddo  !  .. end computation polymer density and charge density  

        ! .. electrostatics 
        ! .. no charge regulating surface charge 
        psisurfR = psi(n)+sigmaqSurfR/2.0d0 
        psisurfL = psi(1)+sigmaqSurfL/2.0d0

        ! .. Poisson Eq 

        f(n+1)  = -0.5d0*( (psi(2)-psi(1)) + sigmaqSurfL  + rhoq(1)*constqW)      !     boundary
        f(2*n)= -0.5d0*( sigmaqSurfR -(psi(n)-psi(n-1)) + rhoq(1)*constqW)      !     boundary

        do i=2,n-1
            f(n+i)= -0.5d0*( psi(i+1)-2.0d0*psi(i) + psi(i-1) +rhoq(i)*constqW)
        enddo

        do i=1,n
            f(2*n+i)=rhopolA(i)-rhopolAin(i)
            f(3*n+i)=rhopolB(i)-rhopolBin(i)
        enddo

        norm=l2norm(f,4*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnelectdouble

    subroutine fcnneutral(x,f,nn)

    !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use vectornorm

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out) :: f(neq)
        integer*8, intent(in) :: nn

        !     .. declare local variables

        real*8 :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real*8 :: xpolBin(nsize)
        real*8 :: xA(3),xB(3),sumxA,sumxB
        real*8 :: constA,constB
        real*8 :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real*8 :: tmp,expVdW 
        real*8 :: norm
        integer :: conf              ! counts number of conformations
        real*8 :: cn                 ! auxilary variable for Poisson Eq


        real*8, parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        n=nz                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                  ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            xpolBin(i)=x(i+n)    ! volume fraction B-polymer 
        enddo

        do i=1,n                  ! init volume fractions 
            xpolAB(i)  = 0.0d0     ! AB polymer volume fraction 
            xpolC(i)   = 0.0d0     ! C polymer volume fraction 
            rhopolA(i) = 0.0d0     ! A polymer density 
            rhopolB(i) = 0.0d0
            rhopolC(i) = 0.0d0
        
            exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
            exppiC(i)=(xsol(i)**vpolC)
       
            !     .. VdW interaction   
            tmp = 0.0d0
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*xpolBin(j)
                enddo
            endif
            expVdW=dexp(-VdWepsB*tmp)
            exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0d0                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0d0                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0d0                 ! init q   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(c,s)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

        rhopolAB0=sigmaAB/qAB
        rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            ! polymer volume fraction A and B all in state 3 ANa and BNa
            xpolAB(i)=(rhopolA(i)*vpolA(3)+rhopolB(i)*vpolB(3))*vsol
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)-1.0d0
            f(n+i)=rhopolB(i)*vpolB(3)*vsol-xpolBin(i)
        enddo  !  .. end computation polymer density 

        norm=l2norm(f,2*n)
        iter=iter+1

        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnneutral


    !     .. function solves for bulk volume fraction 

    subroutine fcnbulk(x,f,nn)   

        !     .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use physconst
        use vectornorm

        implicit none

        !     .. scalar arguments

        integer*8, intent(in) :: nn

        !     .. array arguments

        real*8, intent(in) :: x(neq)
        real*8, intent(out):: f(neq)


        !     .. local variables

        real*8 :: phiNaCl,phiNa,phiCl,phiK,phiKCl
        real*8 :: deltavolNaCl,deltavolKCl,norm

        !     .. executable statements 


        deltavolNaCl=(vNaCl-vNa-vCl)
        deltavolKCl=(vKCl-vNa-vCl)

        phiNa   =x(1)
        phiCl   =x(2)
        phiNaCl =x(3)
        phiK    =x(4)
        phiKCl  =x(5)

        !     .. condensation equilibrium eq NaCl

        f(1)= phiNaCl-(phiNa*phiCl*K0ionNa*vsol*vNaCl/(vNa*vCl*vsol))*      & 
         ((1.0d0 -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolNaCl)

        !     .. charge neutrality
        f(2)=phiNa/vNa-phiCl/vCl+2.0d0*xbulk%Ca/vCa+xbulk%Hplus-xbulk%OHmin+phiK/vK

        !     .. conservation of number NaCl

        f(3)=phiNa/vNa+phiCl/vCl+2.0d0*phiNaCl/vNaCl + phiKCL/vKCL &    
             -2.0d0*vsol*(Na/1.0d24)*cNaCl-dabs(xbulk%Hplus-xbulk%OHmin) &
            -2.0d0*xbulk%Ca/vCa-vsol*(Na/1.0d24)*cKCl

        !     .. condensation equilibrium eq KCl

        f(4)= phiKCl-(phiK*phiCl*K0ionK*vsol*vKCl/(vK*vCl*vsol))* &
             ((1.0d0 -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolKCl)

        !     .. conservation of number K
        f(5)= phiK/vK+phiKCl/vKCl-vsol*(Na/1.0d24)*cKCl

        norm=l2norm(f,5)
        iter=iter+1
     
        print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulk


end module listfcn
