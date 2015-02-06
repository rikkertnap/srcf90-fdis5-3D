!---------------------------------------------------------------|
! augustus  2009                                                | 
! fenergy.f:                                                    |
! calculates the free energy and bulk free energy               |
!---------------------------------------------------------------|


!     .. module file of energy variables

module energy 

    implicit none
  
    !     .. variables
  
    real*8 :: FE                  ! free energy
    real*8 :: FEbulk              ! free energybulk
    real*8 :: deltaFE             ! free energy difference
  
    !     .. auxiliary variable used in free energy computation  
    real*8 :: FEq                 ! partition function poly A and B 
    real*8 :: FEpi                ! sum over pi
    real*8 :: FErho               ! sum over densities
    real*8 :: FEel                ! electrostatics energ
    real*8 :: FEelsurf            ! electrostatics energy from  surface
    real*8 :: FEchem              ! chemical energy
    real*8 :: sumphiA             ! check integrale over phiA
    real*8 :: sumphiB             ! check integrale over phiB
    real*8 :: sumphiC             ! check integrale over phiC
    real*8 :: qres                ! charge charge
    real*8 :: checkphi            ! check integrate over phi
    real*8 :: FEbind              ! complexation contribution
    real*8 :: FEVdW               ! Van der Waals contribution
    real*8 :: FEVdWB
    real*8 :: FEVdWC
    !   
    real*8 :: FEconfAB
    real*8 :: FEConfC
    real*8 :: FEtranssol
    real*8 :: FEalt

contains

    subroutine fcnenergy()
 
        use globals
        implicit none
    
        if(sysflag.eq."elect") call fcnenergy_elect()
        if(sysflag.eq."neutral") call fcnenergy_neutral()
    
    end subroutine fcnenergy

   
    subroutine fcnenergy_elect()
    
        use globals
        use volume
        use parameters
        use field
        use VdW
        use surface
        implicit none

        !  .. local arguments 
    
        real*8 :: sigmaq0,psi0
        real*8 :: qsurf              ! total charge on surface 
        real*8 :: qsurfg             ! total charge on grafting surface 
        real*8 :: sigmaTOL           ! tolerance of surface coverage 
        integer :: i,j               ! dummy variables 
        real*8 :: volumelat          ! volume lattice 
        integer :: nzadius

        !  .. computation of free energy 
    
        FEpi  = 0.0d0
        FErho = 0.0d0
        FEel  = 0.0d0
        FEelsurf = 0.0d0
        sumphiA = 0.0d0
        sumphiB = 0.0d0
        sumphiC = 0.0d0

        FEq = 0.0d0
        FEbind = 0.0d0
        FEchem = 0.0d0
        FEVdWC = 0.0d0
        FEVdWB = 0.0d0     
        qres = 0.0d0

        do i=1,nz
            FEpi = FEpi  + deltaG(i)*dlog(xsol(i))
            FErho = FErho - deltaG(i)*(xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel = FEel  - deltaG(i) * rhoq(i) * psi(i)/2.0d0        
            FEbind = FEbind + fdisA(5,i)*rhopolA(i)+fdisB(5,i)*rhopolB(i)

            do j=1,nz 
                FEVdWC = FEVdWC + deltaG(i)*rhopolC(i)* rhopolC(j)*chis(i,j)
                FEVdWB = FEVdWB + deltaG(i)*rhopolB(i)* rhopolB(j)*chis(i,j)       
            enddo   

            qres = qres + deltaG(i) * rhoq(i)
            sumphiA = sumphiA + deltaG(i) * rhopolA(i)
            sumphiB = sumphiB + deltaG(i) * rhopolB(i)
            sumphiC = sumphiC + deltaG(i) * rhopolC(i)
        enddo
    
        FEel  = (delta/vsol)*FEel
        FEpi  = (delta/vsol)*FEpi
        FErho = (delta/vsol)*FErho
        qres = (delta/vsol)*qres
        sumphiA = delta*sumphiA
        sumphiB = delta*sumphiB
        sumphiC = delta*sumphiC

        FEVdWC  = delta*FEVdWC*VdWepsC*vpolC*vsol/2.0d0   
        FEVdWB  = delta*FEVdWB*VdWepsB*vpolB(3)*vsol/2.0d0   
    
        if (sysflag=="elect") then 
            FEVdW=FEVdWC
        elseif (sysflag=="neutral") then
            FEVdW=FEVdWB
        else 
            print*,"Wrong value sysflag : ", sysflag
            stop    
        endif   


        FEq = -(delta/(vsol))*(sigmaAB*dlog(qAB)+sigmaC*dlog(qC) )
    
        sigmaTOL = 0.00000001     ! tolerance of surface coverage below no polymers 
    
        if(sigmaAB <= sigmaTOL) then 
            FEq = -(delta/(vsol))*(sigmaC*dlog(qC) )
        elseif (sigmaC <= sigmaTOL) then 
            FEq = -(delta/(vsol))*(sigmaAB*dlog(qAB) )
        endif
    
        !     .. surface charge constribution  
       !     NEED WORK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        sigmaq0= sigmaqSurfL/(delta*4.0d0*pi*lb) ! dimensional charge density  
        psi0 = psi(1)+sigmaqSurfL/2.0d0 ! surface potential units !!
    
        FEelsurf = FEelsurf + sigmaq0 * psi0 
        FEelsurf = FEelsurf/2.0d0
    
        FEchem = 0.0d0 !(sigmaSurf/(4.0d0*pi*lb*delta))*dlog(1.0-fdisSurf) -2.0d0*FEelsurf

        !     .. total free energy per area of surface 

        FE = FEq + FEpi + FErho + FEel + FEelsurf + FEchem - FEVdW + FEbind
        
        print*,"FE = " ,FE
    
        qsurf  = sigmaqSurfL/(4.0d0*pi*lb*delta) ! total surface charge
  
        print*,"qsurf=",qsurf,"qres=",qres    

        qres = qres + qsurf  ! total residual charge 

        
        volumelat=nz*delta   ! volume lattice divide by area surface

        FEbulk   = dlog(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ & 
            xbulk%Na/vNa +xbulk%Ca/vCa +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_elect

    

    subroutine fcnenergy_neutral()

        !  .. variable and constant declaractions 
    
        use globals
        use volume
        use parameters
        use field
        use VdW

        implicit none

        !     .. local arguments 
    
        real*8 :: sigmaq0,psi0
        real*8 :: qsurf              ! total charge on surface 
        real*8 :: qsurfg             ! total charge on grafting surface 
        real*8 :: sigmaTOL           ! tolerance of surface coverage 
        integer :: i,j               ! dummy variables 
        real*8 :: volumelat          ! volume lattice 
        integer :: nzadius

        !     .. computation of free energy 
    
        FEpi  = 0.0d0
        FErho = 0.0d0
        FEel  = 0.0d0
        FEelsurf = 0.0d0
        sumphiA = 0.0d0
        sumphiB = 0.0d0
        sumphiC = 0.0d0

        FEq = 0.0d0
        FEbind = 0.0d0
        FEchem = 0.0d0
        FEVdW = 0.0d0 
        qres = 0.0d0

        do i=1,nz
            FEpi = FEpi  + deltaG(i)*dlog(xsol(i))
            FErho = FErho - deltaG(i)*xsol(i) 

            do j=1,nz 
                FEVdW = FEVdW + deltaG(i)*rhopolB(i)* rhopolB(j)*chis(i,j)       
            enddo   

            sumphiA = sumphiA + deltaG(i) * rhopolA(i)
            sumphiB = sumphiB + deltaG(i) * rhopolB(i)
            sumphiC = sumphiC + deltaG(i) * rhopolC(i)
        enddo
    
        FEpi  = (delta/vsol)*FEpi
        FErho = (delta/vsol)*FErho
    
        sumphiA = delta*sumphiA
        sumphiB = delta*sumphiB
        sumphiC = delta*sumphiC

        FEVdW  = delta*FEVdW*(VdWepsB*vpolB(3)*vsol)/2.0d0   

        FEq = -delta*(sigmaAB*dlog(qAB)+sigmaC*dlog(qC) )
    
        sigmaTOL = 0.00000001     ! tolerance of surface coverage below no polymers 
    
        if(sigmaAB <= sigmaTOL) then 
            FEq = -delta*sigmaC*dlog(qC)
        elseif (sigmaC <= sigmaTOL) then 
            FEq = -delta*sigmaAB*dlog(qAB) 
        endif
    
        FEel = 0.0d0    
        FEelsurf =0.0d0
        FEchem =0.0d0
        FEbind =0.0d0
        qres =0.0d0
        qsurf =0.0d0 

        !  .. total free energy per area of surface 

        FE = FEq + FEpi + FErho - FEVdW 
        
        print*,"FE =",FE
    
       
        volumelat=nz*delta  ! volume lattice divide by area surface
        
        print*,"volumelat=",volumelat
        FEbulk  = dlog(xbulk%sol)-xbulk%sol
        print*,"FEbulk=",FEbulk
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE  = FE -FEbulk
    
        ! altnative computation free energy

        call FEconf_neutral()
        FEtranssol=FEtrans(xsol,vsol)     
        FEalt= FEconfAB+FEconfC+FEtranssol+FEVdW
        print*,"FEalt = ",FEalt
    
    end subroutine fcnenergy_neutral

    ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral()

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW

        implicit none
        
        !     .. declare local variables

        real*8 :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 

        integer :: i,j,k,c,s         ! dummy indices
        real*8 :: pro,tmp,expVdW 
        integer :: conf              ! counts number of conformations
        

        real*8, parameter :: tolconst = 1.0d-9  ! tolerance for constA and constB 


        !     .. executable statements 

        do i=1,nz    
            exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
            exppiC(i)=(xsol(i)**vpolC)
       
            !     .. VdW interaction   
            tmp = 0.0d0
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*rhopolB(j)*vpolB(3)*vsol
                enddo
            endif
            expVdW=dexp(-VdWepsB*tmp)
            exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        enddo

        
        FEconfAB=0.0d0

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
            FEconfAB=FEconfAB+pro*dlog(pro)
        enddo

        ! normalize
        FEconfAB=(FEconfAB/qAB-dlog(qAB))*(sigmaAB*delta)    

        FEconfC = 0.0d0                   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0d0               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            FEconfC = FEconfC+pro*dlog(pro)
        enddo
        ! normalize
        FEconfC=(FEConfC/qC -dlog(qC))*(sigmaC*delta)    

    end subroutine FEconf_neutral

    real*8 function FEtrans(xvol,vol)
        use globals
        use field
        use parameters
        implicit none

        real*8, intent(in) :: xvol(nsize)
        real*8, intent(in) :: vol    

        integer :: i

        FEtrans=0.0d0

        do i=1,nz
            FEtrans=FEtrans + deltaG(i)*xvol(i)*(dlog(xvol(i))-1.0d0)
        enddo
        FEtrans= delta *FEtrans /vol
        return 

    end function FEtrans   

  
end module energy
