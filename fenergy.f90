!  .. module file for free energy variables /calculations

module energy 

    use molecules

    implicit none
  
    !     .. variables
    real*8 :: FE                  ! free energy
    real*8 :: FEbulk              ! free energybulk
    real*8 :: deltaFE             ! free energy difference delteFE=FE-FEbulk
  
    !     .. auxiliary variable used in free energy computation  

    real*8 :: FEq                 ! partition function poly A and B 
    real*8 :: FEpi                ! sum over pi
    real*8 :: FErho               ! sum over densities
    real*8 :: FEel                ! electrostatics energ
    real*8 :: FEelsurf(2)         ! electrostatics energy from  surface
    real*8 :: FEchemsurf(2)       ! chemical free energy surface
    real*8 :: FEchem
    real*8 :: FEbind              ! complexation contribution
    real*8 :: FEVdW               ! Van der Waals contribution
    real*8 :: FEVdWB
    real*8 :: FEVdWC 
    real*8 :: FEconfAB
    real*8 :: FEConfC

    real*8 :: FEalt               ! free energy
    real*8 :: FEbulkalt           ! free energybulk
    real*8 :: deltaFEalt          ! free energy difference delteFE=FE-FEbulk

    real*8 :: FEchemsurfalt(2)    ! chemical free energy surface
    real*8 :: diffFEchemsurf(2)   ! difference cheme

    type(moleclist) :: FEtrans,FEchempot,FEtransbulk,FEchempotbulk
    type(moleclist) :: deltaFEtrans,deltaFEchempot

    real*8 :: sumphiA             ! check integral over phiA
    real*8 :: sumphiB             ! check integral over phiB
    real*8 :: sumphiC             ! check integral over phiC
    real*8 :: qres                ! charge charge
    real*8 :: checkphi            ! check integrate over phi
    
    real*8, parameter :: sigmaTOL = 0.00000001     ! tolerance of surface coverage below no polymers 

    private :: sigmaTOL

contains

    subroutine fcnenergy()
 
        use globals
        implicit none
    
        if(sysflag.eq."elect") then 
            call fcnenergy_elect()
        elseif(sysflag.eq."electdouble") then 
            print*,"sysflag=electdouble"    
            call fcnenergy_elect()
        elseif(sysflag.eq."electnopoly") then 
            call fcnenergy_elect()
            call fcnenergy_elect_alternative()
        elseif(sysflag.eq."neutral") then 
            call fcnenergy_neutral()
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value sysflag : ",sysflag
        endif 
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
        real*8 :: qsurf(2)           ! total charge on surface 
        real*8 :: qsurfg             ! total charge on grafting surface  
        integer :: i,j               ! dummy variables 
        real*8 :: volumelat          ! volume lattice 
        integer :: nzadius
        real*8 :: sigmaSurf(2),sigmaqSurf(2),sigmaq0Surf(2),psiSurf(2)

        !  .. computation of free energy 
    
        FEpi  = 0.0d0
        FErho = 0.0d0
        FEel  = 0.0d0
!        FEelsurf = 0.0d0
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
            FEpi = FEpi  + dlog(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel = FEel  - rhoq(i) * psi(i)/2.0d0        
            FEbind = FEbind + fdisA(5,i)*rhopolA(i)+fdisB(5,i)*rhopolB(i)

!            do j=1,nz 
!                FEVdWC = FEVdWC + deltaG(i)*rhopolC(i)* rhopolC(j)*chis(i,j)
!                FEVdWB = FEVdWB + deltaG(i)*rhopolB(i)* rhopolB(j)*chis(i,j)       
!            enddo   

            qres = qres + rhoq(i)
            sumphiA = sumphiA +  rhopolA(i)
            sumphiB = sumphiB +  rhopolB(i)
            sumphiC = sumphiC +  rhopolC(i)
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
    
!        if (sysflag=="elect") then 
!            FEVdW=FEVdWC
!        elseif (sysflag=="neutral") then
!            FEVdW=FEVdWB
!        else 
!            print*,"Wrong value sysflag : ", sysflag
!            stop    
!        endif   

      
        if((sigmaABL > sigmaTOL).and.(sigmaABR > sigmaTOL).and.(sigmaC > sigmaTOL)) then
        
            FEq =-delta*(sigmaABL*dlog(qABL)+sigmaABR*dlog(qABR)+ sigmaABR*dlog(qABR) +sigmaC*dlog(qC) )
        
        elseif((sigmaABL <= sigmaTOL).and.(sigmaABR <= sigmaTOL).and.(sigmaC > sigmaTOL)) then
        
            FEq = -delta*(sigmaC*dlog(qC) )
        
        elseif((sigmaABL > sigmaTOL).and.(sigmaABR <= sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
            FEq = -delta*(sigmaABL*dlog(qABL) )
        
        elseif((sigmaABL > sigmaTOL).and.(sigmaABR > sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
            FEq = -delta*(sigmaABL*dlog(qABL) +sigmaABR*dlog(qABR))
        
        elseif((sigmaABL <= sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
            FEq = 0.0d0
        
        else
        
            print*,"Error in fcnerergy"
            print*,"Something went wrong in evaluating FEq"   
            print*,"sigmaABL=",sigmaABL
            print*,"sigmaABL=",sigmaABL
            print*,"sigmaC=",sigmaC
            stop    
        
        endif
    
        ! .. surface charge constribution 

        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL
        sigmaqSurf(RIGHT) = sigmaqSurfR
        sigmaqSurf(LEFT)  = sigmaqSurfL
        psiSurf(RIGHT)    = psiSurfR
        psiSurf(LEFT)     = psiSurfL
      
        do i = 1,2    
            sigmaq0Surf(i)=  sigmaqSurf(i)/(delta*4.0d0*pi*lb) ! dimensional charge density  
            FEelsurf(i) = sigmaq0Surf(i) * psiSurf(i) /2.0d0 
        enddo   

        if(bcflag(RIGHT)=='qu') then ! quartz

            FEchemSurf(RIGHT) = dlog(fdisS(2))*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay
        
            FEchemSurf(RIGHT) = (dlog(fdisS(2))+qS(2)*psiSurfR)*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurf(RIGHT) =(dlog(fdisS(2))+dlog(fdisS(5)))*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
        
            FEchemSurf(RIGHT)= (dlog(fdisTaR(2))*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb)) -2.0d0*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
        
            FEchemSurf(RIGHT)=0.0d0
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 

            FEchemSurf(LEFT)= dlog(fdisTaL(2))*sigmaSurf(LEFT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(LEFT)
        
        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurf(LEFT)=0.0d0
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 


        !     .. total free energy per area of surface 

        FE = FEq  + FEpi + FErho + FEel + FEelSurf(RIGHT) + FEelSurf(LEFT)
        FE = FE + FEchemSurf(RIGHT)+FEchemSurf(LEFT) - FEVdW + FEbind
        
!        print*,"FE = " ,FE
        
        do i=LEFT,RIGHT     
            qsurf(i) = sigmaqSurf(i)/(4.0d0*pi*lb*delta)
        enddo

!        print*,"qsurf(LEFT)=",qsurf(LEFT),"qsurf(RIGHT)=",qsurf(RIGHT),"qres=",qres    

        qres = qres + qsurf(RIGHT)+qsurf(LEFT)  ! total residual charge 

        
        volumelat=nz*delta   ! volume lattice

        FEbulk   = dlog(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ & 
            xbulk%Na/vNa +xbulk%Ca/vCa +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_elect

    subroutine fcnenergy_elect_alternative()
    
        use globals
        use volume
        use parameters
        use field
        use VdW
        use surface


        implicit none

        !  .. local arguments 
    
        real*8 :: sigmaq0,psi0
        real*8 :: qsurf(2)           ! total charge on surface 
        real*8 :: qsurfg             ! total charge on grafting surface  
        integer :: i,j               ! dummy variables 
        real*8 :: volumelat          ! volume lattice 
        integer :: nzadius
        real*8 :: sigmaSurf(2),sigmaqSurf(2),sigmaq0Surf(2),psiSurf(2)
        real*8 :: diffFEchemTa

        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL
        sigmaqSurf(RIGHT) = sigmaqSurfR
        sigmaqSurf(LEFT)  = sigmaqSurfL
        psiSurf(RIGHT)    = psiSurfR
        psiSurf(LEFT)     = psiSurfL


        !  .. computation of free energy 
    
        !  .. alternative computation free energy

!        call FEconf_neutral()

        ! .. translational entropy 

        FEtrans%sol   = FEtrans_entropy(xsol,xbulk%sol,vsol,"w")   
        FEtrans%Na    = FEtrans_entropy(xNa,xbulk%Na,vNa)
        FEtrans%Cl    = FEtrans_entropy(xCl,xbulk%Cl,vCl)
        FEtrans%Ca    = FEtrans_entropy(xCa,xbulk%Ca,vCa)
        FEtrans%K     = FEtrans_entropy(xK,xbulk%K,vK)
        FEtrans%KCl   = FEtrans_entropy(xKCl,xbulk%KCl,vKCl)
        FEtrans%NaCl  = FEtrans_entropy(xNaCl,xbulk%NaCl,vNaCl)
        FEtrans%Hplus = FEtrans_entropy(xHplus,xbulk%Hplus,vsol,"w")
        FEtrans%OHmin = FEtrans_entropy(xOHmin,xbulk%OHmin,vsol,"w")    


        ! .. chemical potential + standard chemical potential 

        FEchempot%sol   = 0.0d0 ! by construction  
        FEchempot%Na    = FEchem_pot(xNa,expmu%Na,vNa)
        FEchempot%Cl    = FEchem_pot(xCl,expmu%Cl,vCl)
        FEchempot%Ca    = FEchem_pot(xCa,expmu%Ca,vCa)
        FEchempot%K     = FEchem_pot(xK,expmu%K,vK) 
        FEchempot%KCl   = FEchem_pot(xKCl,expmu%KCl,vKCl)
        FEchempot%NaCl  = FEchem_pot(xNaCl,expmu%NaCl,vNaCl)
        FEchempot%Hplus = FEchem_pot(xHplus,expmu%Hplus,vsol,"w")
        FEchempot%OHmin = FEchem_pot(xOHmin,expmu%OHmin,vsol,"w")


        ! .. surface chemical contribution

        if(bcflag(RIGHT)=='qu') then ! quartz
            FEchemSurfalt(RIGHT) = (dlog(fdisS(1))+qS(1)*psiSurfR)*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay        
            FEchemSurfalt(RIGHT) = (dlog(fdisS(1))+qS(1)*psiSurfR)*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
            FEchemSurfalt(RIGHT) =(dlog(fdisS(2))+dlog(fdisS(5)))*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(RIGHT)
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
            FEchemSurfalt(RIGHT)= ((dlog(fdisTaR(1))+qTA(1)*psiSurfR)*sigmaSurf(RIGHT)/(delta*4.0d0*pi*lb)) -2.0d0*FEelsurf(RIGHT)
        elseif(bcflag(RIGHT)=="cc") then  
            FEchemSurfalt(RIGHT)=0.0d0
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 
            FEchemSurfalt(LEFT)= (dlog(fdisTaL(1))+qTA(1)*psiSurfL)*sigmaSurf(LEFT)/(delta*4.0d0*pi*lb) -2.0d0*FEelsurf(LEFT)
        elseif(bcflag(LEFT)=="cc") then  
            FEchemSurfalt(LEFT)=0.0d0
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 



        ! .. summing all contrubutions
        
!        FEalt= FEconfAB+FEconfC+FEVdW
        FEalt = FEtrans%sol +FEtrans%Na+ FEtrans%Cl +FEtrans%NaCl+FEtrans%Ca 
        FEalt = FEalt+FEtrans%OHmin +FEtrans%Hplus +FEtrans%K +FEtrans%KCl
        FEalt = FEalt+FEchempot%sol +FEchempot%Na+ FEchempot%Cl +FEchempot%NaCl+FEchempot%Ca 
        FEalt = FEalt+FEchempot%OHmin +FEchempot%Hplus+ FEchempot%K +FEchempot%K+FEchempot%KCl
        ! be vary carefull FE = -1/2 \int dz rho_q(z) psi(z)

        FEalt = FEalt- FEel + FEelSurf(RIGHT)+FEelSurf(LEFT)+FEchemSurfalt(RIGHT)+FEchemSurfalt(LEFT) 

        ! .. delta translational entropy

        FEtransbulk%sol   = FEtrans_entropy_bulk(xbulk%sol,vsol,"w")   
        FEtransbulk%Na    = FEtrans_entropy_bulk(xbulk%Na,vNa)
        FEtransbulk%Cl    = FEtrans_entropy_bulk(xbulk%Cl,vCl)
        FEtransbulk%Ca    = FEtrans_entropy_bulk(xbulk%Ca,vCa)
        FEtransbulk%K     = FEtrans_entropy_bulk(xbulk%K,vK)
        FEtransbulk%KCl   = FEtrans_entropy_bulk(xbulk%KCl,vKCl)
        FEtransbulk%NaCl  = FEtrans_entropy_bulk(xbulk%NaCl,vNaCl)
        FEtransbulk%Hplus = FEtrans_entropy_bulk(xbulk%Hplus,vsol,"w")
        FEtransbulk%OHmin = FEtrans_entropy_bulk(xbulk%OHmin,vsol,"w")    

        ! .. delta chemical potential + standard chemical potential 

        FEchempotbulk%sol   = 0.0d0 ! by construction  
        FEchempotbulk%Na    = FEchem_pot_bulk(xbulk%Na,expmu%Na,vNa)
        FEchempotbulk%Cl    = FEchem_pot_bulk(xbulk%Cl,expmu%Cl,vCl)
        FEchempotbulk%Ca    = FEchem_pot_bulk(xbulk%Ca,expmu%Ca,vCa)
        FEchempotbulk%K     = FEchem_pot_bulk(xbulk%K,expmu%K,vK) 
        FEchempotbulk%KCl   = FEchem_pot_bulk(xbulk%KCl,expmu%KCl,vKCl)
        FEchempotbulk%NaCl  = FEchem_pot_bulk(xbulk%NaCl,expmu%NaCl,vNaCl)
        FEchempotbulk%Hplus = FEchem_pot_bulk(xbulk%Hplus,expmu%Hplus,vsol,"w")
        FEchempotbulk%OHmin = FEchem_pot_bulk(xbulk%OHmin,expmu%OHmin,vsol,"w")

        
        ! .. bulk free energy

        volumelat=nz*delta   ! volume lattice divide by area surface
        FEbulkalt = FEtransbulk%sol +FEtransbulk%Na+ FEtransbulk%Cl +FEtransbulk%NaCl+FEtransbulk%Ca 
        FEbulkalt = FEbulkalt+FEtransbulk%OHmin +FEtransbulk%Hplus +FEtransbulk%K +FEtransbulk%KCl
        FEbulkalt = FEbulkalt+FEchempotbulk%sol +FEchempotbulk%Na+FEchempotbulk%Cl +FEchempotbulk%NaCl+FEchempotbulk%Ca 
        FEbulkalt = FEbulkalt+FEchempotbulk%OHmin +FEchempotbulk%Hplus -FEchempotbulk%K -FEchempotbulk%KCl
      
        FEbulkalt = volumelat*FEbulkalt

        ! .. delta

        deltaFEtrans%sol   = FEtrans%sol  - FEtransbulk%sol * volumelat
        deltaFEtrans%Na    = FEtrans%Na   - FEtransbulk%Na * volumelat
        deltaFEtrans%Cl    = FEtrans%Cl   - FEtransbulk%Cl * volumelat
        deltaFEtrans%Ca    = FEtrans%Ca   - FEtransbulk%Ca * volumelat
        deltaFEtrans%K     = FEtrans%K    - FEtransbulk%K * volumelat
        deltaFEtrans%KCl   = FEtrans%KCl  - FEtransbulk%KCl * volumelat
        deltaFEtrans%NaCl  = FEtrans%NaCl - FEtransbulk%NaCl * volumelat
        deltaFEtrans%Hplus = FEtrans%Hplus- FEtransbulk%Hplus * volumelat
        deltaFEtrans%OHmin = FEtrans%OHmin- FEtransbulk%OHmin * volumelat
         
        deltaFEchempot%sol   = FEchempot%sol  - FEchempotbulk%sol * volumelat
        deltaFEchempot%Na    = FEchempot%Na   - FEchempotbulk%Na * volumelat
        deltaFEchempot%Cl    = FEchempot%Cl   - FEchempotbulk%Cl * volumelat
        deltaFEchempot%Ca    = FEchempot%Ca   - FEchempotbulk%Ca * volumelat
        deltaFEchempot%K     = FEchempot%K    - FEchempotbulk%K * volumelat
        deltaFEchempot%KCl   = FEchempot%KCl  - FEchempotbulk%KCl * volumelat
        deltaFEchempot%NaCl  = FEchempot%NaCl - FEchempotbulk%NaCl * volumelat
        deltaFEchempot%Hplus = FEchempot%Hplus- FEchempotbulk%Hplus * volumelat
        deltaFEchempot%OHmin = FEchempot%OHmin- FEchempotbulk%OHmin * volumelat


        ! .. differences

        deltaFEalt = FEalt - FEbulkalt


!        print*,"delta FEchemsurfalt(RIGHT)=",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)

        
        if(bcflag(LEFT)=="ta" ) then 
            diffFEchemsurf(LEFT)= (sigmaSurf(LEFT)/(4.0d0*pi*lb*delta))*( dlog(K0Ta(1)) -dlog(expmu%Hplus))
        else
            diffFEchemsurf(LEFT)=0.0d0
        endif

        if(bcflag(RIGHT)=='qu') then ! quartz
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0d0*pi*lb*delta))*( dlog(K0S(1)) -dlog(expmu%Hplus))
        elseif(bcflag(RIGHT)=='cl') then ! quartz
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0d0*pi*lb*delta))*( dlog(K0S(1)) -dlog(expmu%Hplus))
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine       
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0d0*pi*lb*delta))*( dlog(K0Ta(1)) -dlog(expmu%Hplus))
        else
            diffFEchemsurf(RIGHT)=0.0d0
            ! not yet implemented 
        endif 


       
!        print*,"expmu%Hplus=",expmu%Hplus
!        print*,"sigmaSurf(RIGHT)=",sigmaSurf(RIGHT)
!        print*,"dlog(expmu%Hplus))=",dlog(expmu%Hplus)
!        print*,"diffFEchemTa",diffFEchemTa

!        print*,"difference(LEFT)  =", FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
!        print*,"difference(RIGHT) =", FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)


    end subroutine fcnenergy_elect_alternative
    

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

        FEtrans%sol=FEtrans_entropy(xsol,xbulk%sol,vsol,"w")     
        
        FEalt= FEconfAB+FEconfC+FEtrans%sol+FEVdW
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

    real*8 function FEtrans_entropy(xvol,xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real*8, intent(in) :: xvol(nsize)
        real*8, intent(in) :: xvolbulk 
        real*8, intent(in) :: vol    
        character(len=1), optional :: flag    

        integer :: i


        if(xvolbulk==0.0d0) then 
            FEtrans_entropy=0.0d0
        else
            FEtrans_entropy=0.0d0
            if(present(flag)) then
            ! water special case because vsol treated diffetent then vi  
                do i=1,nz
                    FEtrans_entropy=FEtrans_entropy + xvol(i)*(dlog(xvol(i))-1.0d0)
                enddo 
                FEtrans_entropy = delta*FEtrans_entropy/vol
            else 
                do i=1,nz
                    FEtrans_entropy = FEtrans_entropy + xvol(i)*(log(xvol(i)/vol)-1.0d0)
                enddo
                FEtrans_entropy = delta*FEtrans_entropy/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy   

    real*8 function FEchem_pot(xvol,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real*8, intent(in) :: xvol(nsize)
        real*8, intent(in) :: expchempot 
        real*8, intent(in) :: vol    
        character(len=1), optional :: flag    

        ! .. local 
        integer :: i
        real*8 :: chempot ! chemical potential difference 
        real*8 :: sumdens 

        if(expchempot==0.0d0) then 
            FEchem_pot=0.0d0
        else     
            sumdens=0.0d0
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                chempot = -dlog(expchempot)    
                do i=1,nz
                    sumdens=sumdens +xvol(i)
                enddo
                FEchem_pot=delta*chempot*sumdens/vol            
            else
                chempot = -dlog(expchempot/vol)    
                do i=1,nz
                    sumdens=sumdens +xvol(i)
                enddo
                FEchem_pot=delta*chempot*sumdens/(vol*vsol)               
            endif
        endif

    end function FEchem_pot   

    real*8 function FEtrans_entropy_bulk(xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real*8, intent(in) :: xvolbulk 
        real*8, intent(in) :: vol    
        character(len=1), optional :: flag    

        integer :: i


        if(xvolbulk==0.0d0) then 
            FEtrans_entropy_bulk=0.0d0
        else
            FEtrans_entropy_bulk=0.0d0
            if(present(flag)) then
                ! water special case because vsol treated diffetent then vi  
                FEtrans_entropy_bulk=xvolbulk*(dlog(xvolbulk)-1.0d0)/vol
            else 
                FEtrans_entropy_bulk=xvolbulk*(dlog(xvolbulk/vol)-1.0d0)/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy_bulk   

    real*8 function FEchem_pot_bulk(xvolbulk,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real*8, intent(in) :: xvolbulk
        real*8, intent(in) :: expchempot 
        real*8, intent(in) :: vol    
        character(len=1), optional :: flag    

        ! .. local 
        integer :: i
        real*8 :: chempot ! chemical potential difference 
        real*8 :: sumdens 

        if(expchempot==0.0d0) then 
            FEchem_pot_bulk=0.0d0
        else     
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                FEchem_pot_bulk=-dlog(expchempot)*xvolbulk/vol            
            else
                FEchem_pot_bulk=-dlog(expchempot/vol)*xvolbulk/(vol*vsol)               
            endif
        endif

    end function FEchem_pot_bulk   

  
end module energy
