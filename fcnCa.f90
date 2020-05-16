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


    ! brush of multiblock copolymers

    subroutine fcnelectbrushmulti(x,f,nn)

        use mpivars
        use globals
        use parameters, Tlocal=>Tref 
        use volume
        use chains
        use field
        use vectornorm
        use VdW 
        use surface
        use Poisson

        implicit none
        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)
        real(dp) :: local_q(ngr_node)
        real(dp) :: exppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: rhopol_tmp(nsize,nsegtypes)
        real(dp) :: pro
        integer  :: n,i,j,k,l,c,s,ln,t,g,gn   ! dummy indices
        real(dp) :: norm
        real(dp) :: rhopol0 !integra_q
        integer  :: noffset
        


        !     .. executable statements 

        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, size-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        n=nsize
        
        ! read out x 
        k=(nsegtypes+1)*n
       
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            psi(i)  = x(i+k)      ! potential
        enddo           
        
        do t=1,nsegtypes
            k=t*n
            do i=1,n 
                rhopolin(i,t) = x(i+k) ! density 
            enddo    
        enddo
             
        !  .. assign global and local polymer density 

        do t=1,nsegtypes
            do i=1,n
                rhopol(i,t)=0.0_dp 
                local_rhopol(i,t)=0.0_dp
            enddo    
        enddo    
       
        do i=1,n                  ! init volume fractions
            xpol(i)    = 0.0_dp                                   ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction 
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
        enddo

        !  fdis(i,t) is assocaited with fraction of monomer of type t at i in state 2 
        !  acid : AH  <=> A- +H+
        !  base : BH+ <=> B + H+

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,n                                         
                    fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                    exppi(i,t) = (xsol(i)**vpol(t))*exp(-zpol(t,2)*psi(i) )/fdis(i,t)   ! auxilary variable palpha
                enddo  
            else
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    exppi(i,t) = xsol(i)**vpol(t)
                enddo  
            endif   
        enddo      
       
        ! Van der Waals   
        if(isVdW) then 
            call VdW_contribution_exp_diblock(rhopolin(:,1),rhopolin(:,2),exppi(:,1),exppi(:,2))
            !do t=1,nsegtypes  
            !    call VdW_contribution_exp_multi(rhopolin,exppi(:,t),t)
            !enddo
        endif 

        !  .. computation polymer volume fraction 

         do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node    
            local_q(gn) = 0.0_dp    ! init q
             
            do c=1,cuantas         ! loop over cuantas
                pro=1.0_dp         ! initial weight conformation (1 or 0)
                do s=1,nseg        ! loop over segments 
                    k=indexchain(s,gn,c)
                    t=type_of_monomer(s)                
                    pro = pro*exppi(k,t)
                enddo    
                local_q(gn) = local_q(gn)+pro
                do s=1,nseg
                    k=indexchain(s,gn,c) 
                    t=type_of_monomer(s)
                    rhopol_tmp(k,t)=rhopol_tmp(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
                enddo
            enddo

            do t=1,nsegtypes
                do i=1,n                   ! normalization with local_qi(g) 
                    local_rhopol(i,t)=local_rhopol(i,t)+rhopol_tmp(i,t)/local_q(gn)
                enddo
            enddo
                
        enddo    
              
        
        !   .. import results 

        if (rank==0) then 
          
            do t=1,nsegtypes
                call MPI_REDUCE(local_rhopol(:,t), rhopol(:,t), nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            enddo
            
             
            do gn=1,ngr_node
                g = gn+ rank*ngr_node
                q(g)=local_q(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(local_q, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    q(g)=local_q(gn)
                enddo
            enddo 
          
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=1.0_dp/volcell ! volume polymer segment per volume cell

            do t=1, nsegtypes
                do i=1,n
                    rhopol(i,t) = rhopol0 * rhopol(i,t)       ! density polymer of type t  
                    xpol(i)     = xpol(i) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    if(ismonomer_chargeable(t)) then 
                        rhoqpol(i)  = rhoqpol(i) + (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol ! charge density polymer
                    endif    
                    ! ... scf eq for density
                    f(i+t*n)    = rhopol(i,t) - rhopolin(i,t) 
                enddo
            enddo        


            do i=1,n
                f(i)    = xpol(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)-1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb ! total charge density in units of vsol  
            !   print*,i,rhoq(i)
            enddo
          
            !     .. end computation polymer density and charge density  

          
            ! .. electrostatics 

            ! .. surface charge  
            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

            ! noffset=nsegtypes*n

            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL) ! work on noffset 
 
            ! .. selfconsistent boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

            ! k=(nsegtypes+1)*n            
             ! f(k+1)=-0.5_dp*(gplus(1)*(psi(2)-psi(1))+ rhoq(1)*constqW)     

            !  .. end electrostatics   
         
            norm=l2norm(f,(nsegtypes+2)*n)
            iter=iter+1

            print*,'iter=', iter ,'norm=',norm
         
        else                      ! Export results 
         
            dest = 0 
           
            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t), nsize,MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo

            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
         
        endif


    end subroutine fcnelectbrushmulti


    subroutine fcnelectNoPoly(x,f,nn)

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
                xNa(i)    = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
                xK(i)     = expmu%K *(xsol(i)**vK )*exp(-psi(i)*zK ) ! ion plus volume fraction
                xRb(i)    = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! ion plus volume fraction
                xCa(i)    = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
                xMg(i)    = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! ion divalent pos volume fraction
                xNaCl(i)  = expmu%NaCl*(xsol(i)**vNaCl)              ! ion pair  volume fraction
                xKCl(i)   = expmu%KCl *(xsol(i)**vKCl)               ! ion pair  volume fraction
                xCl(i)    = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
                xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
                xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH-  volume fraction
            enddo

            !   .. construction of fcn 
            
            do i=1,n

                f(i)=xsol(i)+xNa(i)+xCl(i)+xRb(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xMg(i)+&
                    xHplus(i)+xOHmin(i)-1.0_dp
           
                rhoq(i)= zNa*xNa(i)/vNa+zRb*xRb(i)/vRb+zCa*xCa(i)/vCa + zMg*xMg(i)/vMg +zK*xK(i)/vK + &
                    zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)
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

        real(dp) :: exppiA(nsize),exppiB(nsize)  ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize)    
        real(dp) :: rhopolAL_tmp(nsize),rhopolBL_tmp(nsize)
        real(dp) :: rhopolAL_local(nsize),rhopolBL_local(nsize)
        real(dp) :: qABL_local(ngr_node)

        real(dp) :: xA(3),xB(3),sumxA,sumxB, sgxA, sgxB, qAD, qBD 
        real(dp) :: constA,constB
        real(dp) :: proL,rhopolABL0,proR,rhopolABR0
        integer :: n,ix,iy,iz,neq_bc                  
        integer :: i,j,k,kL,kR,c,s,g,gn        ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr

        real(dp):: checkintegralAB,sumrhopolAB
        
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
    
        n=nsize                      ! size vector neq=4*nsizeß x=(pi,psi,rhopolA,rhopolB]
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
   
        if(rank.eq.0) then 
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
        endif

        do i=1,n                     ! init volume fractions 
            rhopolAL_local(i) = 0.0_dp     ! A polymer density 
            rhopolBL_local(i) = 0.0_dp     ! B polymer density  
            rhopolAL_tmp(i)   = 0.0_dp     ! 
            rhopolBL_tmp(i)   = 0.0_dp
           
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))      ! OH-  volume fraction
       
            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))      ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))   ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))   ! ACa+/A-
       
            sgxA=1.0_dp+xA(1)+xA(2)+xA(3)                                                         
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) 
            qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

            fdisA(i,1)  = 1.0_dp/qAD   ! removed double minus sign 
            fdisA(i,5)  = (fdisA(i,1)**2)*constA
            fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
            fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
            fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
       
            xB(1)= xHplus(i)/(K0b(1)*(xsol(i) **deltavB(1)))     ! BH/B-
            xB(2)= (xNa(i)/vNa)/(K0b(2)*(xsol(i)**deltavB(2)))   ! BNa/B-
            xB(3)= (xCa(i)/vCa)/(K0b(3)*(xsol(i)**deltavB(3)))   ! BCa+/B-
       
            sgxB=xB(1)+xB(2)+xB(3)
            constB=(2.0_dp*(rhopolBin(i)*vsol)*(xCa(i)/vCa))/(K0b(4)*(xsol(i)**deltavB(4)))
            qBD = (sgxB+sqrt(sgxB*sgxB+4.0_dp*constB))/2.0_dp  ! remove minus

            fdisB(i,1)  = 1.0_dp/qBD
            fdisB(i,5)  = (fdisB(i,1)**2)*constB
            fdisB(i,2)  = fdisB(i,1)*xB(1)                      ! BH 
            fdisB(i,3)  = fdisB(i,1)*xB(2)                      ! BNa 
            fdisB(i,4)  = fdisB(i,1)*xB(3)                      ! BCa+ 
       
            ! use A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*exp(-zpolB(1)*psi(i))/fdisB(i,1) ! auxiliary variable

        enddo

        if(rank==0) then            ! global polymer density
            do i=1,n
                xpolAB(i) = 0.0_dp 
                rhopolAL(i) = 0.0_dp     ! A polymer density 
                rhopolBL(i) = 0.0_dp     ! B polymer density  
               
            enddo
            do g=1,ngr
                qABL(g)=0.0_dp
            enddo
        endif

        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qABL_local(gn)=0.0_dp       ! init qB
        enddo 

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=1.0_dp

                do s=1,nsegAB              ! loop over segments 
                    kL=indexchainAB(s,gn,c)           
                
                    if(isAmonomer(s)) then ! A segment 
                        proL = proL*exppiA(kL)
                    else
                        proL = proL*exppiB(kL)
                    endif
                enddo

                qABL_local(gn) = qABL_local(gn)+proL

                do s=1,nsegAB
                    kL=indexchainAB(s,gn,c)
        
                    if(isAmonomer(s)) then ! A segment 
                        rhopolAL_tmp(kL)=rhopolAL_tmp(kL)+proL
                    else
                        rhopolBL_tmp(kL)=rhopolBL_tmp(kL)+proL
                    endif
                enddo
          
            
            enddo   ! end cuantas loop 

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolAL_local(i)=rhopolAL_local(i)+rhopolAL_tmp(i)/qABL_local(gn)
                rhopolAL_tmp(i)=0.0_dp
                rhopolBL_local(i)=rhopolBL_local(i)+rhopolBL_tmp(i)/qABL_local(gn)
                rhopolBL_tmp(i)=0.0_dp
            enddo

        enddo    

        !     .. import results

        if (rank==0) then
          
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
    

            do gn=1,ngr_node
                g = gn+ rank*ngr_node
                qABL(g)=qABL_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    qABL(g)=qABL_local(gn)
                enddo
            enddo 

            !   .. construction of fcn and volume fraction polymer        
            rhopolABL0=1.0_dp/volcell ! volume polymer segment per volume cell

            do i=1,n

                rhopolAL(i) = rhopolABL0*rhopolAL(i)   ! /deltaG(i) not nessesary since deltaG=1
                rhopolA(i)  = rhopolAL(i) 

                rhopolBL(i) = rhopolABL0*rhopolBL(i)
                rhopolB(i)  = rhopolBL(i) 
            
                do k=1,4               ! polymer volume fraction
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(i,k)*vpolA(k)*vsol  & 
                        +rhopolB(i)*fdisB(i,k)*vpolB(k)*vsol
                enddo    

                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(i,5)*vpolA(5)*vsol/2.0_dp)
                xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(i,5)*vpolB(5)*vsol/2.0_dp)
       
                f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                !rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                !    zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                !    zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol  

                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    ((zpolA(1)*fdisA(i,1)+ zpolA(4)*fdisA(i,4))*rhopolA(i) + &
                     (zpolB(1)*fdisB(i,1)+ zpolB(4)*fdisB(i,4))*rhopolB(i) )*vsol         
       
                !   ..  total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

            ! .. electrostatics 

            ! .. surface charge  
            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

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
         
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)

        endif


    end subroutine fcnelect


    ! ... fcn for homopolymer A 
    !  .. A weak acid monomre with potential Ca binding and  VdW inteaction 

    subroutine fcnelectA(x,f,nn)

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW
        use surface 
        use vectornorm
        use myutils
        use surface 
        use Poisson


        implicit none

        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: exppiA(nsize) ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize)    
        real(dp) :: rhopolAL_tmp(nsize)
        real(dp) :: rhopolAL_local(nsize)
        real(dp) :: qABL_local(ngr_node)

        real(dp) :: xA(3),sumxA, sgxA, qAD
        real(dp) :: constA
        real(dp) :: proL,rhopolABL0
        integer :: n,ix,iy,iz,neq_bc                  
        integer :: i,j,k,kL,kR,c,s,g,gn        ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr

        real(dp):: checkintegralAB,sumrhopolAB
        
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
    
        n=nsize                      ! size vector neq=2*nsize x=(pi,psi)
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolAin(i)=x(i+2*n)
        enddo
   
        ! This part needs to be checked 
        if(rank.eq.0) then 
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
        endif

        ! end to be checked

        do i=1,n                     ! init volume fractions 
            
            rhopolAL_local(i) = 0.0_dp     ! A polymer density 
            rhopolAL_tmp(i)   = 0.0_dp     ! 
           
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))      ! OH-  volume fraction


            xA(1)= xHplus(i)/(K0a(1)*(xsol(i)**deltavA(1)))      ! AH/A-
            xA(2)= (xNa(i)/vNa)/(K0a(2)*(xsol(i)**deltavA(2)))   ! ANa/A-
            xA(3)= (xCa(i)/vCa)/(K0a(3)*(xsol(i)**deltavA(3)))   ! ACa+/A-
       
            sgxA=1.0_dp+xA(1)+xA(2)+xA(3)                                                         
            constA=(2.0_dp*(rhopolAin(i)*vsol)*(xCa(i)/vCa))/(K0a(4)*(xsol(i)**deltavA(4))) 
            qAD = (sgxA+sqrt(sgxA*sgxA+4.0_dp*constA))/2.0_dp  ! remove minus

            fdisA(i,1)  = 1.0_dp/qAD   ! removed double minus sign 
            fdisA(i,5)  = (fdisA(i,1)**2)*constA
            fdisA(i,2)  = fdisA(i,1)*xA(1)                       ! AH 
            fdisA(i,3)  = fdisA(i,1)*xA(2)                       ! ANa 
            fdisA(i,4)  = fdisA(i,1)*xA(3)                       ! ACa+ 
           
            ! use A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable

        enddo

        if(isVdW) call VdW_contribution_exp(rhopolAin,exppiA)

        if(rank==0) then            ! global polymer density
            do i=1,n
                xpolAB(i) = 0.0_dp 
                rhopolAL(i) = 0.0_dp     ! A polymer density 
            enddo
            do g=1,ngr
                qABL(g)=0.0_dp
            enddo
        endif

        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qABL_local(gn)=0.0_dp       ! init qB
        enddo 

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=1.0_dp

                do s=1,nsegAB              ! loop over segments 
                    kL=indexchainAB(s,gn,c)           
                    proL = proL*exppiA(kL)
                enddo

                qABL_local(gn) = qABL_local(gn)+proL

                do s=1,nsegAB
                    kL=indexchainAB(s,gn,c)
                    ! A segment 
                    rhopolAL_tmp(kL)=rhopolAL_tmp(kL)+proL
                enddo
          
            enddo   ! end cuantas loop 

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolAL_local(i)=rhopolAL_local(i)+rhopolAL_tmp(i)/qABL_local(gn)
                rhopolAL_tmp(i)=0.0_dp
            enddo

        enddo    

        !     .. import results

        if (rank==0) then
          
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
          
            do gn=1,ngr_node
                g = gn+ rank*ngr_node
                qABL(g)=qABL_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    qABL(g)=qABL_local(gn)
                enddo
            enddo 

            !   .. construction of fcn and volume fraction polymer        
            rhopolABL0=1.0_dp/volcell ! volume polymer segment per volume cell

            do i=1,n

                rhopolAL(i) = rhopolABL0*rhopolAL(i)   ! /deltaG(i) not nessesary since deltaG=1
                rhopolA(i)  = rhopolAL(i) 
  
                do k=1,4               ! polymer volume fraction
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(i,k)*vpolA(k)*vsol 
                enddo    
                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(i,5)*vpolA(5)*vsol/2.0_dp)
       
                f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    ((zpolA(1)*fdisA(i,1)+ zpolA(4)*fdisA(i,4))*rhopolA(i) )*vsol         
       
                !   ..  total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

            ! .. electrostatics 

            ! .. surface charge  
            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
 
            ! .. selfconsistent boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

            do i=1,n
                f(2*n+i)=rhopolA(i)-rhopolAin(i)
            enddo

            norm=l2norm(f,3*n)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)

        else          ! Export results

            dest = 0
         
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)

        endif

    end subroutine fcnelectA

    !  .. with calcium binding to A monomers 
    !  .. A weak acid monomre with potential Ca binding 
    !  .. B neutal monomer 
    !  .. VdW inteaction for A and B 

    subroutine fcnelectVdWAB(x,f,nn)

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

        real(dp) :: exppiA(nsize),exppiB(nsize)  ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolBin(nsize)    
        real(dp) :: rhopolAL_tmp(nsize),rhopolBL_tmp(nsize)
        real(dp) :: rhopolAL_local(nsize),rhopolBL_local(nsize)
        real(dp) :: qABL_local(ngr_node)

        real(dp) :: xA(3),xB(3),sumxA,sumxB, sgxA, sgxB, qAD, qBD 
        real(dp) :: constA,constB
        real(dp) :: proL,rhopolABL0,proR,rhopolABR0
        integer :: n,ix,iy,iz,neq_bc                  
        integer :: i,j,k,kL,kR,c,s,g,gn        ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr

        real(dp):: checkintegralAB,sumrhopolAB
        
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
    
        n=nsize                      ! size vector neq=4*nsizeß x=(pi,psi,rhopolA,rhopolB]
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
   
        if(rank.eq.0) then 
            neq_bc=0
            if(bcflag(RIGHT)/="cc") then
                neq_bc=nx*ny
                do i=1,neq_bc
                    !psiSurfR(i) =x(2*n+i)    ! surface potentail
                    psiSurfR(i) =x(4*n+i)    ! surface potentail
                    
                enddo
            endif   
            if(bcflag(LEFT)/="cc") then 
                do i=1,nx*ny
                    ! psiSurfL(i) =x(2*n+neq_bc+i) ! surface potentail
                    psiSurfL(i) =x(4*n+neq_bc+i) ! surface potentail
                enddo
                neq_bc=neq_bc+nx*ny
            endif    
        endif

        do i=1,n                     ! init volume fractions 
            rhopolAL_local(i) = 0.0_dp     ! A polymer density 
            rhopolBL_local(i) = 0.0_dp     ! B polymer density  
            rhopolAL_tmp(i)   = 0.0_dp     ! 
            rhopolBL_tmp(i)   = 0.0_dp
           
            xNa(i)   = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))      ! OH-  volume fraction
       
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
       
            ! use A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(2))  ! auxiliary variable

        enddo

        if(isVdW) call VdW_contribution_exp_diblock(rhopolAin,rhopolBin,exppiA,exppiB)


        if(rank==0) then            ! global polymer density
            do i=1,n
                xpolAB(i) = 0.0_dp 
                rhopolAL(i) = 0.0_dp     ! A polymer density 
                rhopolBL(i) = 0.0_dp     ! B polymer density  
            enddo
            do g=1,ngr
                qABL(g)=0.0_dp
            enddo
        endif

        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qABL_local(gn)=0.0_dp       ! init qB
        enddo 

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=1.0_dp

                do s=1,nsegAB              ! loop over segments 
                    kL=indexchainAB(s,gn,c)           
                
                    if(isAmonomer(s)) then ! A segment 
                        proL = proL*exppiA(kL)
                    else
                        proL = proL*exppiB(kL)
                    endif
                enddo

                qABL_local(gn) = qABL_local(gn)+proL

                do s=1,nsegAB
                    kL=indexchainAB(s,gn,c)
        
                    if(isAmonomer(s)) then ! A segment 
                        rhopolAL_tmp(kL)=rhopolAL_tmp(kL)+proL
                    else
                        rhopolBL_tmp(kL)=rhopolBL_tmp(kL)+proL
                    endif
                enddo
          
            
            enddo   ! end cuantas loop 

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolAL_local(i)=rhopolAL_local(i)+rhopolAL_tmp(i)/qABL_local(gn)
                rhopolAL_tmp(i)=0.0_dp
                rhopolBL_local(i)=rhopolBL_local(i)+rhopolBL_tmp(i)/qABL_local(gn)
                rhopolBL_tmp(i)=0.0_dp
            enddo

        enddo    

        !     .. import results

        if (rank==0) then
          
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
    

            do gn=1,ngr_node
                g = gn+ rank*ngr_node
                qABL(g)=qABL_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    qABL(g)=qABL_local(gn)
                enddo
            enddo 

            !   .. construction of fcn and volume fraction polymer        
            rhopolABL0=1.0_dp/volcell ! volume polymer segment per volume cell

            do i=1,n

                rhopolAL(i) = rhopolABL0*rhopolAL(i)   ! /deltaG(i) not nessesary since deltaG=1
                rhopolA(i)  = rhopolAL(i) 

                rhopolBL(i) = rhopolABL0*rhopolBL(i)
                rhopolB(i)  = rhopolBL(i) 
            
                do k=1,4               ! polymer volume fraction
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(i,k)*vpolA(k)*vsol  & 
                        +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
                enddo    

                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(i,5)*vpolA(5)*vsol/2.0_dp)
                xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(i,5)*vpolB(5)*vsol/2.0_dp)
       
                f(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    ((zpolA(1)*fdisA(i,1)+ zpolA(4)*fdisA(i,4))*rhopolA(i) + &
                     (zpolB(1)*fdisB(i,1)+ zpolB(4)*fdisB(i,4))*rhopolB(i) )*vsol         
       
                !   ..  total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

            ! .. electrostatics 

            ! .. surface charge  
            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

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
         
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)

        endif


    end subroutine fcnelectVdWAB


    subroutine fcnelectdouble(x,fvec,nn)

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
        real(dp), intent(out) :: fvec(neq)
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
        integer :: n,ix,iy,iz                  ! half of n
        integer :: i,j,k,kL,kR,c,s,g,gn   ! dummy indices
        real(dp) :: norm
        integer :: conf               ! counts number of conformations
        real(dp) :: cn                ! auxilary variable for Poisson Eq
        character(len=lenText) :: text, istr, rstr

        real(dp):: checkintegralAB,sumrhopolAB
        
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
    
        n=nsize                      ! size vector neq=4*nsizeß x=(pi,psi,rhopolA,rhopolB]
    
        do i=1,n                     ! init x 
            xsol(i)= x(i)            ! solvent volume fraction 
            psi(i) = x(i+n)          ! potential
            rhopolAin(i)=x(i+2*n)
            rhopolBin(i)=x(i+3*n)
        enddo
   
        do i=1,n                     ! init volume fractions 
            rhopolAL_local(i) = 0.0_dp     ! A polymer density 
            rhopolBL_local(i) = 0.0_dp     ! B polymer density  
            rhopolAR_local(i) = 0.0_dp     ! A polymer density 
            rhopolBR_local(i) = 0.0_dp     ! B polymer density 
            rhopolAL_tmp(i)   = 0.0_dp     ! 
            rhopolBL_tmp(i)   = 0.0_dp
            rhopolAR_tmp(i)   = 0.0_dp     ! 
            rhopolBR_tmp(i)   = 0.0_dp

            xNa(i)   = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! ion plus volume fraction
            xK(i)    = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! ion plus volume fraction
            xCa(i)   = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! ion divalent pos volume fraction
            xNaCl(i) = expmu%NaCl*(xsol(i)**vNaCl)               ! ion pair  volume fraction
            xKCl(i)  = expmu%KCl*(xsol(i)**vKCl)                 ! ion pair  volume fraction
            xCl(i)   = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! ion neg volume fraction
            xHplus(i) = expmu%Hplus*(xsol(i))*exp(-psi(i))      ! H+  volume fraction
            xOHmin(i) = expmu%OHmin*(xsol(i))*exp(+psi(i))      ! OH-  volume fraction
       
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

            fdisB(i,1)  = 1.0_dp/qBD
            fdisB(i,5)  = (fdisB(i,1)**2)*constB
            fdisB(i,2)  = fdisB(i,1)*xB(1)                      ! BH 
            fdisB(i,3)  = fdisB(i,1)*xB(2)                      ! BNa 
            fdisB(i,4)  = fdisB(i,1)*xB(3)                      ! BCa+ 
       
            ! use A^- reference state

            exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(i,1) ! auxiliary variable

        enddo

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

        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qABL_local(gn)=0.0_dp       ! init qB
            qABR_local(gn)=0.0_dp
        enddo 

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                proL=0.0_dp                ! initial weight conformation 
                proR=0.0_dp
            
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
            
            enddo   ! end cuantas loop 

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
          
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolAR_local, rhopolAR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBR_local, rhopolBR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)

            do gn=1,ngr_node
                g = gn+ rank*ngr_node
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
                    xpolAB(i)=xpolAB(i)+rhopolA(i)*fdisA(i,k)*vpolA(k)*vsol  & 
                        +rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
                enddo    

                xpolAB(i)=xpolAB(i)+rhopolA(i)*(fdisA(i,5)*vpolA(5)*vsol/2.0_dp)
                xpolAB(i)=xpolAB(i)+rhopolB(i)*(fdisB(i,5)*vpolB(5)*vsol/2.0_dp)
       
                fvec(i)=xpolAB(i)+xsol(i)+xNa(i)+xCl(i)+xNaCl(i)+xK(i)+xKCl(i)+xCa(i)+xHplus(i)+xOHmin(i)-1.0_dp
       
                !rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                !    zpolA(1)*fdisA(1,i)*rhopolA(i)*vsol+ zpolA(4)*fdisA(4,i)*rhopolA(i)*vsol+ &
                !    zpolB(1)*fdisB(1,i)*rhopolB(i)*vsol+ zpolB(4)*fdisB(4,i)*rhopolB(i)*vsol  

                rhoq(i)= zNa*xNa(i)/vNa + zCa*xCa(i)/vCa +zK*xK(i)/vK + zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    ((zpolA(1)*fdisA(i,1)+ zpolA(4)*fdisA(i,4))*rhopolA(i) + &
                     (zpolB(1)*fdisB(i,1)+ zpolB(4)*fdisB(i,4))*rhopolB(i) )*vsol         
       
                !   ..  total charge density in units of vsol
            enddo  !  .. end computation polymer density and charge density  

           
            !call check_integral_rholpolAB(sumrhopolAB, checkintegralAB)
            !print*,"sumrholAB=",sumrhopolAB," check= ",checkintegralAB
            !print*,"qABL=",qABL
            !print*,"qABR=",qABR

            ! .. electrostatics 

            ! .. surface charge  
            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)

            call Poisson_Equation(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
 
            ! .. selfconsistent boundary conditions
            call Poisson_Equation_Surface(fvec,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

            do i=1,n
                fvec(2*n+i)=rhopolA(i)-rhopolAin(i)
                fvec(3*n+i)=rhopolB(i)-rhopolBin(i)
            enddo

            norm=l2norm(fvec,4*n)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)


        else          ! Export results

            dest = 0
         
            call MPI_REDUCE(rhopolAL_local, rhopolAL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBL_local, rhopolBL, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolAR_local, rhopolAR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(rhopolBR_local, rhopolBR, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)

            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)

        endif

    end subroutine fcnelectdouble


    subroutine fcnneutral(x,f,nn)

        use globals
        use volume
        use chains
        use field
        use parameters
        use vectornorm
        use myutils

        implicit none
        !     .. scalar arguments
        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)
        integer(8), intent(in) :: nn

        !     .. declare local variables

        real(dp) :: exppiA(nsize) ! auxilairy variable for computing P(\alpha) 
        real(dp) :: rhopolAin(nsize),rhopolA_tmp(nsize),rhopolA_local(nsize)
        real(dp) :: qAB_local(ngr_node)
        real(dp) :: pro,rhopolAB0
        integer :: n,ix,iy,iz              
        integer :: i,j,k,c,s,g,gn        ! dummy indices
        real(dp) :: norm
        character(len=lenText) :: text, istr, rstr

        real(dp):: checkintegralAB,sumrhopolAB
        
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
    
        n=nsize                      ! size vector neq=nsize 
    
        do i=1,n                          ! init x 
            xsol(i)= x(i)                 ! solvent volume fraction 
            rhopolA_local(i) = 0.0_dp     ! A polymer density 
            rhopolA_tmp(i)   = 0.0_dp     ! 
            exppiA(i)=(xsol(i)**vpolA(2)) ! auxiliary variable    
        enddo

        if(rank==0) then                  ! global polymer density
            do i=1,n
                xpolAB(i) = 0.0_dp 
                rhopolA(i) = 0.0_dp       ! A polymer density 
            enddo
            do g=1,ngr
                qAB(g)=0.0_dp
            enddo
        endif

        !     .. computation polymer volume fraction 
        do gn=1,ngr_node
            qAB_local(gn)=0.0_dp           ! init qB
        enddo 

        do gn=1,ngr_node                   ! loop over grafted points <=>  grafted area on different nodes 
 
            g=gn+rank*ngr_node
     
            do c=1,cuantasAB               ! loop over cuantas
            
                pro=1.0_dp

                do s=1,nsegAB              ! loop over segments 
                    k=indexchainAB(s,gn,c)           
                    pro = pro*exppiA(k)
                enddo

                qAB_local(gn) = qAB_local(gn)+pro

                do s=1,nsegAB
                    k=indexchainAB(s,gn,c)
                    rhopolA_tmp(k)=rhopolA_tmp(k)+pro
                enddo
            enddo   ! end cuantas loop 

            do i=1,n                   ! normalization with local_qi(g) 
                rhopolA_local(i)=rhopolA_local(i)+rhopolA_tmp(i)/qAB_local(gn)
                rhopolA_tmp(i)=0.0_dp
            enddo
        enddo    

        !     .. import results

        if (rank==0) then
          
            call MPI_REDUCE(rhopolA_local, rhopolA, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            
            do gn=1,ngr_node
                g = gn+ rank*ngr_node
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

            !   .. construction of fcn and volume fraction polymer        
            rhopolAB0=1.0_dp/volcell ! volume polymer segment per volume cell

            do i=1,n
                rhopolA(i) = rhopolAB0*rhopolA(i)  
                xpolAB(i) = rhopolA(i)*vpolA(2)*vsol 
                f(i) = xpolAB(i)+xsol(i)-1.0_dp     
            enddo  

            norm=l2norm(f,n)
            iter=iter+1
        
            write(rstr,'(E25.16)')norm
            write(istr,'(I6)')iter
            text="iter = "//trim(istr)//" fnorm = "//trim(rstr)
            call print_to_log(LogUnit,text)

        else          ! Export results
            dest = 0         
            call MPI_REDUCE(rhopolA_local, rhopolA, nsize, MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qAB_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif

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


    !     set constrains on vector x  depending on systype value

    subroutine set_contraints(constr)
    
        use precision_definition
        use globals, only : systype, neq , nsize,  LEFT, RIGHT, bcflag
        use volume, only : nx,ny 


        implicit none
            
        real(dp), intent(inout):: constr(:)

        integer :: i, neqint ,neq_bc

        neqint=int(neq,kind(neqint))     ! explict conversion from integer(8) to integer
        
        neq_bc=0 
        if(bcflag(LEFT)/="cc") neq_bc=neq_bc+nx*ny
        if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+nx*ny
    
        select case (systype)
            case ("brush_mul")                 ! multi copolymer:

                do i=1,nsize                    
                    constr(i)=1.0_dp           ! solvent volume fraction   
                    constr(i+nsize)=0.0_dp     ! electrostatic potential
                    constr(i+2*nsize)=1.0_dp   ! number density A
                    constr(i+3*nsize)=1.0_dp   ! number density B 
                enddo  
                do i=1,neq_bc                  ! surface electrostatic potential if bcflag/=cc
                    constr(i+4*nsize)=0.0_dp
                enddo    
 

            case ("elect")                     ! AB copolymer: weak acid A weak acid B

                do i=1,nsize                    
                    constr(i)=1.0_dp           ! solvent volume fraction   
                    constr(i+nsize)=0.0_dp     ! electrostatic potential
                    constr(i+2*nsize)=1.0_dp   ! number density A
                    constr(i+3*nsize)=1.0_dp   ! number density B 
                enddo  
                do i=1,neq_bc                  ! surface electrostatic potential if bcflag/=cc
                    constr(i+4*nsize)=0.0_dp
                enddo    
 
            case ("electA")                    ! A homopolymer weak acid A plus VdW

                do i=1,nsize
                    constr(i)=1.0_dp
                    constr(i+nsize)=0.0_dp   
                    constr(i+2*nsize)=1.0_dp 
                enddo  
                do i=1,neq_bc
                    constr(i+3*nsize)=0.0_dp
                enddo    

             case ("electVdWAB")               ! AB copolymer: weak acid A and neutral B plus VdW A and B

                do i=1,nsize
                    constr(i)=1.0_dp
                    constr(i+nsize)=0.0_dp   
                    constr(i+2*nsize)=1.0_dp
                    constr(i+3*nsize)=1.0_dp 
                enddo  
                do i=1,neq_bc
                    constr(i+4*nsize)=0.0_dp
                enddo    

            
            case ("electdouble")   

                do i=1,nsize
                    constr(i)=1.0_dp
                    constr(i+nsize)=0.0_dp   
                    constr(i+2*nsize)=1.0_dp
                    constr(i+3*nsize)=1.0_dp 
                enddo  

                do i=1,neq_bc
                    constr(i+4*nsize)=0.0_dp
                enddo  

            case ("electnopoly")              ! no polymers only surface chgarges  

                do i=1,nsize
                    constr(i)=2.0_dp
                    constr(i+nsize)=0.0_dp
                enddo  
            
            case ("neutral")                 ! neutral polymers
            
                do i=1,nsize
                    constr(i)=1.0_dp
                enddo 
              
            case default
            
                do i=1,neqint
                    constr(i)=1.0_dp
                enddo 
        
        end select  

    end subroutine set_contraints


 ! selects correct fcn function 

    subroutine set_fcn

        use globals
        use fcnpointer
        
        implicit none   

        select case (systype) 
        case ("brush_mul")
            fcnptr => fcnelectbrushmulti   ! acid and base : no counterion binding VdW
        case ("elect")                  ! copolymer weak polyacid, no VdW
             fcnptr => fcnelect
        case ("electA")                 ! homopolymer weak polyacid VdW 
            fcnptr => fcnelectA         
        case ("electVdWAB")             ! copolymer weak polyacid, VdW
            fcnptr => fcnelectVdWAB        
        case ("electdouble")            ! copolymer weak polyacid, no VdW, both surfaces grafted with polymers		
            fcnptr => fcnelectdouble
        case ("electnopoly")            ! electolyte solution only 
            fcnptr => fcnelectNoPoly 
        case ("neutral")                ! homopolymer neutral
            fcnptr => fcnneutral
        case ("bulk water")             ! determines compositon bulk electrolyte solution
             fcnptr => fcnbulk
        case default
            print*,"Error in call to set_fcn subroutine"    
            print*,"Wrong value systype : ", systype
            stop
        end select  
        
    end subroutine set_fcn


end module listfcn
