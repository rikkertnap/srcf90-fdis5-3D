! --------------------------------------------------------------|
! fcnCa.f90:                                                    |
! constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for weak poly-     |
! electrolytes onto a tethered planar surface                   |
! --------------------------------------------------------------|


module fcnpointer

    implicit none

    abstract interface
        subroutine fcn(x,f,n)
            use precision_definition
            implicit none
            integer(8), intent(in) :: n
            real(dp), dimension(n), intent(in) :: x
            real(dp), dimension(n), intent(out) :: f
        end subroutine fcn
    end interface

    procedure(fcn), pointer :: fcnptr => null()

end module fcnpointer


module listfcn
   
    use mpivars
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
        integer(8), intent(in) :: nn
        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize),xpolCin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s,g,gn        ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn                 ! auxilary variable for Poisson Eq


        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


!         !     .. executable statements 

!         n=nz                      ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

!         do i=1,n                  ! init x 
!             xsol(i)= x(i)          ! solvent volume fraction 
!             psi(i) = x(i+n)        ! potential
!             rhopolAin(i)=x(i+2*n)
!             rhopolBin(i)=x(i+3*n)
!             xpolCin(i)=x(i+4*n)    ! volume fraction C-polymer
!         enddo

!         psiSurfL = psi(1)          ! surface potentail

!         do i=1,n                  ! init volume fractions 
!             xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
!             xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
!             rhopolA(i) = 0.0_dp     ! A polymer density 
!             rhopolB(i) = 0.0_dp
!             rhopolC(i) = 0.0_dp
        
!             xNa(i)   = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
!             xK(i)    = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
!             xCa(i)   = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
!             xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
!             xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
!             xCl(i)   = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
!             xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
!             xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
       
!             xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))     ! AH/A-
!             xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))  ! ANa/A-
!             xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))  ! ACa+/A-
       
!             sumxA=xA(1)+xA(2)+xA(3)
!             constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
!             if(constA<=tolconst) then 
!                 fdisA(1,i)=1.0_dp/(1.0_dp+sumxA)
!                 fdisA(5,i)=0.0_dp
!             else
!                 fdisA(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constA/((sumxA+1.0_dp)**2)))
!                 fdisA(1,i)= fdisA(1,i)*(sumxA+1.0_dp)/(2.0_dp*constA)
!                 fdisA(5,i)= (fdisA(1,i)**2)*constA
!             endif    
       
!             fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
!             fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
!             fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
!             xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
!             xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
!             xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
!             sumxB=xB(1)+xB(2)+xB(3)
!             constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
!             if(constB<=tolconst) then
!                 fdisB(1,i)=1.0_dp/(1.0_dp+sumxB)
!                 fdisB(5,i)=0.0_dp
!             else
!                 fdisB(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constB/((sumxB+1.0_dp)**2)))
!                 fdisB(1,i)= fdisB(1,i)*(sumxB+1.0_dp)/(2.0_dp*constB) !B^-
!                 fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
!             endif
       
!             fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
!             fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
!             fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       

!     !        exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
!     !        exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

!     !       exppiA(i)=(xsol(i)**vpolA(2))*dexp(-zpolA(2)*psi(i))/fdisA(2,i) ! auxiliary variable                                           
!     !       exppiB(i)=(xsol(i)**vpolB(2))*dexp(-zpolB(2)*psi(i))/fdisB(2,i) ! auxiliary variable   
!             ! Na condensed ANa reference state
!             exppiA(i)=(xsol(i)**vpolA(3))*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
!             exppiB(i)=(xsol(i)**vpolB(3))*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   

       
!             !     .. VdW interaction   
!             tmp = 0.0_dp
!             if((i+VdWcutoffdelta)<=nsize) then 
!                 do j=minrange(i),i+VdWcutoffdelta
!                     tmp = tmp + chis(i,j)*xpolCin(j)
!                 enddo
!             endif
!             expVdW=dexp(-VdWepsC*tmp)
!             exppiC(i)=(xsol(i)**vpolC)*expVdW ! auxiliary variable
!         enddo

!         !   .. computation polymer volume fraction 

!         qAB = 0.0_dp                 ! init q

!         do c=1,cuantasAB             ! loop over cuantas
!             pro=1.0_dp               ! initial weight conformation 
!             do s=1,nsegAB            ! loop over segments 
!                 k=indexchainAB(s,gn,c)
!                 if(isAmonomer(s)) then ! A segment 
!                     pro = pro*exppiA(k)
!                 else
!                     pro = pro*exppiB(k)
!                 endif
!             enddo

!             qAB = qAB+pro
!             do s=1,nsegAB
!                 k=indexchainAB(s,gn,c)
!                 if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
!                     rhopolA(k)=rhopolA(k)+pro
!                 else
!                     rhopolB(k)=rhopolB(k)+pro
!                 endif
!             enddo
!         enddo

!         qC = 0.0_dp                 ! init q   
!         do c=1,cuantasC             ! loop over cuantas                                                      
!             pro=1.0_dp              ! initial weight conformation                                                   
!             do s=1,nsegC            ! loop over segments                
!                 k=indexchainC(s,gn,c)
!                 pro = pro*exppiC(k)
!             enddo
!             qC = qC+pro
!             do s=1,nsegC
!                 k=indexchainC(s,gn,c)
!                 rhopolC(k)=rhopolC(k)+pro
!             enddo
!         enddo

!         !   .. construction of fcn and volume fraction polymer        

! !        rhopolAB0=sigmaAB/qAB
! !        rhopolC0=sigmaC/qC

!         do i=1,n
!             rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
!             rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
!             rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
!             do k=1,4               ! polymer volume fraction
!                 xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
!                     +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
!             enddo    
!             xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0_dp)
!             xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0_dp)
       
!             xpolC(i)=rhopolC(i)*vpolC*vsol
       
!             f(i)=xpolAB(i)+xpolC(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
!             rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
!                 zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
!                 zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
!             !   ..  total charge density in units of vsol
!         enddo  !  .. end computation polymer density and charge density  

!         ! .. electrostatics 

!         sigmaqSurfL=0.0_dp ! charge regulating surface charge 
!         psi(n+1)=0.0_dp   ! bulk potential

!         !    .. Poisson Eq 

!         f(n+1)= -0.5_dp*((psi(2)-psi(1)) + sigmaqSurfL +rhoq(1)*constqW)      !     boundary

!         do i=2,n
!             f(n+i)= -0.5_dp*(psi(i+1)-2.0_dp*psi(i) + psi(i-1) +rhoq(i)*constqW)
!         enddo

!         do i=1,n
!             f(2*n+i)=rhopolA(i)-rhopolAin(i)
!             f(3*n+i)=rhopolB(i)-rhopolBin(i)
!             f(4*n+i)=xpolC(i)-xpolCin(i)
!         enddo

!         iter=iter+1

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
        use myutils
        use Poisson

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn


        !     .. declare local variables

        integer :: n                 ! half of n
        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        integer :: neq_bc           
        character(len=lenText) :: text, istr, rstr
        !     .. executable statements 

        !     .. communication between processors: only done to keep fc call in line with fcnelect and fcnelectdouble

        if (rank.eq.0) then      ! 
        
            flag_solver = 1      !  continue program
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1,   MPI_INTEGER, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x,neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo    
        
            n=nsize                    ! size vector neq=2*nsize x=(pi,psi)

            do i=1,n                   ! init x 
                xsol(i)= x(i)          ! solvent volume fraction 
                psi(i) = x(i+n)        ! potential
            enddo
            
            neq_bc=0
            if(bcflag(RIGHT)/="cc") then
                neq_bc=nx*ny
                do i=1,neq_bc
                    psiSurfR(i) =x(2*n+i)                 ! surface potentail
                enddo
            endif   
            if(bcflag(LEFT)/="cc") then 
                do i=1,nx*ny
                    psiSurfL(i) =x(2*n+neq_bc+i)          ! surface potentail
                enddo
                neq_bc=neq_bc+nx*ny
            endif    


            do i=1,n                  ! init volume fractions 
                xNa(i)    = expmu%Na*(xsol(i)**vNa)*dexp(-psi(i)*zNa) ! ion plus volume fraction
                xK(i)     = expmu%K*(xsol(i)**vK)*dexp(-psi(i)*zK)    ! ion plus volume fraction
                xCa(i)    = expmu%Ca*(xsol(i)**vCa)*dexp(-psi(i)*zCa) ! ion divalent pos volume fraction
                xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
                xKCl(i)   = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
                xCl(i)    = expmu%Cl*(xsol(i)**vCl)*dexp(-psi(i)*zCl) ! ion neg volume fraction
                xHplus(i) = expmu%Hplus*(xsol(i))*dexp(-psi(i))      ! H+  volume fraction
                xOHmin(i) = expmu%OHmin*(xsol(i))*dexp(+psi(i))      ! OH-  volume fraction
            enddo

            !   .. construction of fcn 
            
            do i=1,n

                f(i)=xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
           
                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)
                !   ..  total charge density in units of vsol
            enddo 

            ! .. electrostatics 
            ! .. charge regulating surface charge 

            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)
            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

            ! .. selfconsistent boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)
                    
            norm=l2norm(f,2*n+neq_bc)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)

        endif

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
        use myutils
        use Poisson

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn


        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize)
        real(dp) :: rhopolA_tmp(nsize),rhopolB_tmp(nsize)
        real(dp) :: rhopolA_local(nsize),rhopolB_local(nsize)
        real(dp) :: qAB_local(ngr_node)


        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s,g,gn         ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn               ! auxilary variable for Poisson Eq
        integer :: neq_bc           
        character(len=lenText) :: text, istr, rstr
        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 
 
        !     .. executable statements 
        !     .. communication between processors

        if (rank.eq.0) then
            flag_solver = 1      !  continue program
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1,   MPI_INTEGER, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x,neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo    
        endif
            
        n=nsize                      ! size vector 4*n x=(pi,psi,rhopolA,rhopolB]
   
        do i=1,n                   ! init x 
            xsol(i)= x(i)          ! solvent volume fraction 
            psi(i) = x(i+n)        ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
        
        neq_bc=0
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(4*n+i) ! surface potentail
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(4*n+neq_bc+i)          ! surface potentail
            enddo
            neq_bc=neq_bc+nx*ny
        endif    
        
        do i=1,n                    ! init volume fractions 
            xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
            xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
            rhopolA(i) = 0.0_dp     ! A polymer density 
            rhopolB(i) = 0.0_dp

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
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) ! A2Ca/(A-)^2
            if(constA<=tolconst) then 
                fdisA(1,i)=1.0_dp/(1.0_dp+sumxA)
                fdisA(5,i)=0.0_dp
            else
                fdisA(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constA/((sumxA+1.0_dp)**2)))
                fdisA(1,i)= fdisA(1,i)*(sumxA+1.0_dp)/(2.0_dp*constA)
                fdisA(5,i)= (fdisA(1,i)**2)*constA
            endif    
       
            fdisA(2,i)  = fdisA(1,i)*xA(1)                      ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                      ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                      ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))    ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))  ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))  ! BCa+/B-
       
       
            sumxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            if(constB<=tolconst) then
                fdisB(1,i)=1.0_dp/(1.0_dp+sumxB)
                fdisB(5,i)=0.0_dp
            else
                fdisB(1,i)= (-1.0_dp+dsqrt(1.0_dp+4.0_dp*constB/((sumxB+1.0_dp)**2)))
                fdisB(1,i)= fdisB(1,i)*(sumxB+1.0_dp)/(2.0_dp*constB) !B^-
                fdisB(5,i)= (fdisB(1,i)**2)*constB                    ! B2Ca
            endif
       
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       
            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable
       
        enddo

        if(rank==0) then            ! global polymer density
            do i=1,n
                xpolAB(i)  = 0.0_dp 
                rhopolA(i) = 0.0_dp     ! A polymer density 
                rhopolB(i) = 0.0_dp     ! B polymer density 
            enddo
            do g=1,ngr
                qAB(g) = 0.0_dp
            enddo
        endif
        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qAB_local(gn)=0.0d0      ! init qB
        enddo 

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                pro=1.0_dp                ! initial weight conformation 
            
                if(weightchainAB(gn,c)) then ! initial weight conformation 

                    pro=1.0_dp
  
                    do s=1,nsegAB              ! loop over segments 
                        k=indexchainAB(s,gn,c)             
                        if(isAmonomer(s)) then ! A segment 
                            pro = pro*exppiA(k)
                        else
                            pro = pro*exppiB(k)
                        endif
                    enddo

                    qAB_local(gn) = qAB_local(gn)+pro
                    
                    do s=1,nsegAB
                        k=indexchainAB(s,gn,c)
                        if(isAmonomer(s)) then ! A segment 
                            rhopolA_tmp(k)=rhopolA_tmp(k)+pro
                        else
                            rhopolB_tmp(k)=rhopolB_tmp(k)+pro
                        endif
                    enddo
                
                endif
            
            enddo

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolA_local(i)=rhopolA_local(i)+rhopolA_tmp(i)/qAB_local(gn)
                rhopolA_tmp(i)=0.0_dp
                rhopolB_local(i)=rhopolB_local(i)+rhopolB_tmp(i)/qAB_local(gn)
                rhopolB_tmp(i)=0.0_dp
            enddo

        enddo    

        !     .. import results

        if (rank==0) then

            call MPI_REDUCE(rhopolA_local, rhopolA, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolB_local, rhopolB, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            
            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                qAB(g)=qAB_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(qAB_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    qAB(g)=qAB_local(gn)
                enddo
            enddo 
  
            !  .. construction of fcn and volume fraction polymer   

            rhopolAB0=1.0_dp/volcell ! volume polymer segment per volume cell
            
            do i=1,n

                rhopolA(i) = rhopolAB0*rhopolA(i)   ! /deltaG(i) not nessesary since deltaG=1
                rhopolB(i) = rhopolAB0*rhopolBL(i)
            
                do k=1,4               ! polymer volume fraction
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                        +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
                enddo    
                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0_dp)
                xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0_dp)
       
                f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                    zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
                    !  .. total charge density in units of vsol
            enddo   !  .. end computation polymer density and charge density  


            ! .. electrostatics 
            ! .. charge regulating surface charge 

            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
            ! ..selfconsistent boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

            do i=1,n
                f(2*n+i)=rhopolA(i)-rhopolAin(i)
                f(3*n+i)=rhopolB(i)-rhopolBin(i)
            enddo

            qABL=qAB   ! defined both qAB and qABL, communcicates value     

            norm=l2norm(f,4*n+neq_bc)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)


        else          ! Export results

            dest = 0
         
            call MPI_REDUCE(rhopolA_local, rhopolA, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolB_local, rhopolB, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qAB_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
    
        endif

    end subroutine fcnelect

    subroutine fcnelectdouble(x,f,nn)

        ! .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface
        use vectornorm
        use myutils
        use poisson

        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)  ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize)    
        real(dp) :: rhopolAL_tmp(nsize),rhopolBL_tmp(nsize)
        real(dp) :: rhopolAR_tmp(nsize),rhopolBR_tmp(nsize)
        real(dp) :: rhopolAL_local(nsize),rhopolBL_local(nsize)
        real(dp) :: rhopolAR_local(nsize),rhopolBR_local(nsize)
        real(dp) :: qABL_local(ngr_node), qABR_local(ngr_node)
        real(dp) :: xA(3),xB(3),sumxA,sumxB, sgxA, sgxB, qAD, qBD 
        real(dp) :: constA,constB
        real(dp) :: proL,rhopolABL0,proR,rhopolABR0
        integer :: n                  ! half of n
        integer :: i,j,k,kL,kR,c,s,g,gn   ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr

        ! quick fix for psiSurf 
        !real(dp) :: psiSurfR_local(nx*ny), psiSurfL_local(nx*ny)
        
        ! real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 
        
        !     .. executable statements 
        !     .. communication between processors

        print*,"fcn: hello rank=", rank

        
        if (rank.eq.0) then
            flag_solver = 1      !  continue program
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1,   MPI_INTEGER, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x,neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo    
        endif
        print*,"fcn: send rank=", rank
        n=nsize                      ! size vector 4*n x=(pi,psi,rhopolA,rhopolB]
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo

        !psiSurfR = psi(1)            ! surface potentail
        !psiSurfL = psi(n)
   
        do i=1,n                     ! init volume fractions 
            rhopolAL_local(i) = 0.0_dp     ! A polymer density 
            rhopolBL_local(i) = 0.0_dp     ! B polymer density  
            rhopolAR_local(i) = 0.0_dp     ! A polymer density 
            rhopolBR_local(i) = 0.0_dp     ! B polymer density 

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
       
            sgxA=1.0_dp+xA(1)+xA(2)+xA(3)                                                         
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) 
            qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

            fdisA(1,i)  = 1.0_dp/qAD   ! removed double minus sign 
            fdisA(5,i)  = (fdisA(1,i)**2)*constA
            fdisA(2,i)  = fdisA(1,i)*xA(1)                       ! AH 
            fdisA(3,i)  = fdisA(1,i)*xA(2)                       ! ANa 
            fdisA(4,i)  = fdisA(1,i)*xA(3)                       ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))     ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))   ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))   ! BCa+/B-
       
       
            sgxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            qBD = (sgxB+sqrt(sgxB*sgxB+4.0_dp*constB))/2.0_dp  ! remove minus

            fdisB(1,i)  = 1.0_dp/qBD
            fdisB(5,i)  = (fdisB(1,i)**2)*constB
            fdisB(2,i)  = fdisB(1,i)*xB(1)                      ! BH 
            fdisB(3,i)  = fdisB(1,i)*xB(2)                      ! BNa 
            fdisB(4,i)  = fdisB(1,i)*xB(3)                      ! BCa+ 
       
            ! use A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable

        enddo

         print*,"fcn: compute rank=", rank

        if(rank==0) then            ! global polymer density
            do i=1,n
                xpolAB(i) = 0.0_dp 
                rhopolAL(i) = 0.0_dp     ! A polymer density 
                rhopolBL(i) = 0.0_dp     ! B polymer density  
                rhopolAR(i) = 0.0_dp     ! A polymer density 
                rhopolBR(i) = 0.0_dp     ! B polymer density 
            enddo
            do g=1,ngr
                qABL(g)=0.0_dp
                qABR(g)=0.0_dp
            enddo
        endif

        print*,"fcn: compute polymer  rank=", rank
        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qABL_local(gn)=0.0d0       ! init qB
            qABR_local(gn)=0.0d0
        enddo 


        print*,"fcn: ngr_node=",ngr_node," cuantasAB=",cuantasAB

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=1.0_dp                ! initial weight conformation 
                proR=1.0_dp
            
                if(weightchainAB(gn,c)) then ! initial weight conformation 

                    proL=1.0_dp
                    proR=1.0_dp

                    do s=1,nsegAB              ! loop over segments 
                        kL=indexchainAB(s,gn,c)           
                        kR=mirror_index(kL,nz)    
                        if(isAmonomer(s)) then ! A segment 
                            proL = proL*exppiA(kL)
                            proR = proR*exppiA(kR)
                        else
                            proL = proL*exppiB(kL)
                            proR = proR*exppiB(kR)
                        endif
                    enddo

                    qABL_local(gn) = qABL_local(gn)+proL
                    qABR_local(gn) = qABR_local(gn)+proR

                    do s=1,nsegAB
                        kL=indexchainAB(s,gn,c)
                        kR=mirror_index(kL,nz)
                        if(isAmonomer(s)) then ! A segment 
                            rhopolAL_tmp(kL)=rhopolAL_tmp(kL)+proL
                            rhopolAR_tmp(kR)=rhopolAR_tmp(kR)+proR
                        else
                            rhopolBL_tmp(kL)=rhopolBL_tmp(kL)+proL
                            rhopolBR_tmp(kR)=rhopolBR_tmp(kR)+proR
                        endif
                    enddo
                
                endif
            
            enddo

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolAL_local(i)=rhopolAL_local(i)+rhopolAL_tmp(i)/qABL_local(gn)
                rhopolAL_tmp(i)=0.0_dp
                rhopolBL_local(i)=rhopolBL_local(i)+rhopolBL_tmp(i)/qABL_local(gn)
                rhopolBL_tmp(i)=0.0_dp
                rhopolAR_local(i)=rhopolAR_local(i)+rhopolAR_tmp(i)/qABR_local(gn)
                rhopolAR_tmp(i)=0.0_dp
                rhopolBR_local(i)=rhopolBR_local(i)+rhopolBR_tmp(i)/qABR_local(gn)
                rhopolBR_tmp(i)=0.0_dp
            enddo

        enddo    

        !     .. import results

        if (rank==0) then
            print*,"fcn: reduce"

            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolAR_local, rhopolAR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBR_local, rhopolBR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)

            print*,"fcn: reduce complete"

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                qABL(g)=qABL_local(gn)
                qABR(g)=qABR_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(qABR_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)                
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    qABL(g)=qABL_local(gn)
                    qABR(g)=qABR_local(gn)
                enddo
            enddo 

            !   .. construction of fcn and volume fraction polymer        
            rhopolABL0=1.0_dp/volcell ! volume polymer segment per volume cell
            rhopolABR0=1.0_dp/volcell 

            do i=1,n

                rhopolAL(i) = rhopolABL0*rhopolAL(i)   ! /deltaG(i) not nessesary since deltaG=1
                rhopolAR(i) = rhopolABR0*rhopolAR(i)
                rhopolA(i)  = rhopolAL(i) + rhopolAR(i)

                rhopolBL(i) = rhopolABL0*rhopolBL(i)
                rhopolBR(i) = rhopolABR0*rhopolBR(i)
                rhopolB(i)  = rhopolBL(i) + rhopolBR(i)
            
                do k=1,4               ! polymer volume fraction
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol  & 
                        +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
                enddo    

                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(5,i)*vpolA(5)*vsol/2.0_dp)
                xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(5,i)*vpolB(5)*vsol/2.0_dp)
       
                f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                    zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol         
       
                !   ..  total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

            ! .. electrostatics 
            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
 
            ! .. selfconsistent boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

            do i=1,n
                f(2*n+i)=rhopolA(i)-rhopolAin(i)
                f(3*n+i)=rhopolB(i)-rhopolBin(i)
            enddo

            norm=l2norm(f,4*n)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)


        else          ! Export results

            dest = 0
            print*,"fcn, rank=",rank
         
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolAR_local, rhopolAR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBR_local, rhopolBR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)

            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)

        endif



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

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 
        real(dp) :: xpolBin(nsize)
        real(dp) :: xA(3),xB(3),sumxA,sumxB
        real(dp) :: constA,constB
        real(dp) :: pro,rhopolAB0,rhopolC0
        integer :: n                 ! half of n
        integer :: i,j,k,c,s,gn         ! dummy indices
        real(dp) :: tmp,expVdW 
        real(dp) :: norm
        integer :: conf              ! counts number of conformations
        real(dp) :: cn               ! auxilary variable for Poisson Eq


        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 

        !     .. executable statements 

        n=nsize                        ! size vector neq=5*nz x=(pi,psi,rhopolA,rhopolB,xpolC)

        do i=1,n                    ! init x 
            xsol(i)= x(i)           ! solvent volume fraction 
            xpolBin(i)=x(i+n)       ! volume fraction B-polymer 
        enddo

        do i=1,n                    ! init volume fractions 
            xpolAB(i)  = 0.0_dp     ! AB polymer volume fraction 
            xpolC(i)   = 0.0_dp     ! C polymer volume fraction 
            rhopolA(i) = 0.0_dp     ! A polymer density 
            rhopolB(i) = 0.0_dp
            rhopolC(i) = 0.0_dp
        
            exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
            exppiC(i)=(xsol(i)**vpolC)
       
            !     .. VdW interaction   
            tmp = 0.0_dp
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*xpolBin(j)
                enddo
            endif
            expVdW=dexp(-VdWepsB*tmp)
            exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        enddo

        !   .. computation polymer volume fraction 

        qAB = 0.0_dp                 ! init q

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0_dp                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(s,gn,c)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

            qAB = qAB+pro
            do s=1,nsegAB
                k=indexchainAB(s,gn,c)
                if(isAmonomer(s)) then ! A segment  !        if(isAmonomer(s).eqv..TRUE.) then ! A segment 
                    rhopolA(k)=rhopolA(k)+pro
                else
                    rhopolB(k)=rhopolB(k)+pro
                endif
            enddo
        enddo

        qC = 0.0_dp                 ! init q   
        do c=1,cuantasC             ! loop over cuantas                                                      
            pro=1.0_dp              ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(s,gn,c)
                pro = pro*exppiC(k)
            enddo
            qC = qC+pro
            do s=1,nsegC
                k=indexchainC(s,gn,c)
                rhopolC(k)=rhopolC(k)+pro
            enddo
        enddo

        !   .. construction of fcn and volume fraction polymer        

!        rhopolAB0=sigmaAB/qAB
!       rhopolC0=sigmaC/qC

        do i=1,n
            rhopolA(i)= rhopolAB0*rhopolA(i)/deltaG(i)
            rhopolB(i)= rhopolAB0*rhopolB(i)/deltaG(i)
            rhopolC(i)= rhopolC0*rhopolC(i)/deltaG(i)
       
            ! polymer volume fraction A and B all in state 3 ANa and BNa
            xpolAB(i)=(rhopolA(i)*vpolA(3)+rhopolB(i)*vpolB(3))*vsol
            xpolC(i)=rhopolC(i)*vpolC*vsol
       
            f(i)=xpolAB(i)+xpolC(i)+xsol(i)-1.0_dp
            f(n+i)=rhopolB(i)*vpolB(3)*vsol-xpolBin(i)
        enddo  !  .. end computation polymer density 

    !    norm=l2norm(f,2*n)
        iter=iter+1

    !    print*,'iter=', iter ,'norm=',norm

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

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out):: f(neq)


        !     .. local variables

        real(dp) :: phiNaCl,phiNa,phiCl,phiK,phiKCl
        real(dp) :: deltavolNaCl,deltavolKCl,norm

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
         ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolNaCl)

        !     .. charge neutrality
        f(2)=phiNa/vNa-phiCl/vCl+2.0_dp*xbulk%Ca/vCa+xbulk%Hplus-xbulk%OHmin+phiK/vK

        !     .. conservation of number NaCl

        f(3)=phiNa/vNa+phiCl/vCl+2.0_dp*phiNaCl/vNaCl + phiKCL/vKCL &    
             -2.0_dp*vsol*(Na/1.0e24_dp)*cNaCl-dabs(xbulk%Hplus-xbulk%OHmin) &
            -2.0_dp*xbulk%Ca/vCa-vsol*(Na/1.0e24_dp)*cKCl

        !     .. condensation equilibrium eq KCl

        f(4)= phiKCl-(phiK*phiCl*K0ionK*vsol*vKCl/(vK*vCl*vsol))* &
             ((1.0_dp -phiNaCl-phiNa-phiCl-phiK-phiKCl-xbulk%Hplus-xbulk%OHmin-xbulk%Ca)**deltavolKCl)

        !     .. conservation of number K
        f(5)= phiK/vK+phiKCl/vKCl-vsol*(Na/1.0e24_dp)*cKCl

    !    norm=l2norm(f,5)
        iter=iter+1
     
    !    print*,'iter=', iter ,'norm=',norm

    end subroutine fcnbulk




 ! selects correct fcn function 

    subroutine set_fcn

        use globals
        use fcnpointer
        
        implicit none   

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
        
    end subroutine set_fcn





end module listfcn
