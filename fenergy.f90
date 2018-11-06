
!  .. module file for free energy variables /calculations

module energy 

    use precision_definition
    use molecules
    
    implicit none
  
    !     .. variables
    real(dp) :: FE                  ! free energy
    real(dp) :: FEbulk              ! free energybulk
    real(dp) :: deltaFE             ! free energy difference delteFE=FE-FEbulk
  
    !     .. auxiliary variable used in free energy computation  

    real(dp) :: FEq                 ! partition function poly A and B 
    real(dp) :: FEpi                ! sum over pi
    real(dp) :: FErho               ! sum over densities
    real(dp) :: FEel                ! electrostatics energ
    real(dp) :: FEelsurf(2)         ! electrostatics energy from  surface
    real(dp) :: FEchemsurf(2)       ! chemical free energy surface
    real(dp) :: FEchem
    real(dp) :: FEbind              ! complexation contribution
    real(dp) :: FEVdW               ! Van der Waals contribution
    real(dp) :: FEVdWB
    real(dp) :: FEVdWC 
    real(dp) :: FEconfAB
    real(dp) :: FEConfC

    real(dp) :: FEalt               ! free energy
    real(dp) :: FEbulkalt           ! free energybulk
    real(dp) :: deltaFEalt          ! free energy difference delteFE=FE-FEbulk

    real(dp) :: FEchemsurfalt(2)    ! chemical free energy surface
    real(dp) :: diffFEchemsurf(2)   ! difference cheme

    type(moleclist) :: FEtrans,FEchempot,FEtransbulk,FEchempotbulk
    type(moleclist) :: deltaFEtrans,deltaFEchempot

    real(dp) :: sumphiA             ! check integral over phiA
    real(dp) :: sumphiB             ! check integral over phiB
    real(dp) :: sumphiC             ! check integral over phiC
    real(dp) :: qres                ! charge charge
    real(dp) :: checkphi            ! check integrate over phi
    
    real(dp), parameter :: sigmaTOL = 0.00000001_dp     ! tolerance of surface coverage below no polymers 

    private :: sigmaTOL

contains

    subroutine fcnenergy()
 
        use globals
        implicit none
    
        if(sysflag.eq."elect") then 
            call fcnenergy_elect()
        elseif(sysflag.eq."dipolarstrong") then
            !call fcnenergy_elect()
            print*,"Free energy: warning not yet implemented "
        elseif(sysflag.eq."dipolarweak") then
            !call fcnenergy_elect()
            print*,"Free energy: warning not yet implemented "
        elseif(sysflag.eq."electdouble") then 
        !    print*,"fcnenergy sysflag=electdouble"    
            call fcnenergy_elect()
            call fcnenergy_elect_alternative()
        elseif(sysflag.eq."electnopoly") then 
            call fcnenergy_elect()
            !call fcnenergy_elect_alternative()
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
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf(2)           ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface  
        integer  :: i,j,s               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: nzadius
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx),sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: FEchemSurftmp

        !  .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
!        FEelsurf = 0.0_dp
        sumphiA = 0.0_dp
        sumphiB = 0.0_dp
        sumphiC = 0.0_dp

        FEq = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdWC = 0.0_dp
        FEVdWB = 0.0_dp     
        qres = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel = FEel  - rhoq(i) * psi(i)/2.0_dp        
            FEbind = FEbind + fdisA(5,i)*rhopolA(i)+fdisB(5,i)*rhopolB(i)

            qres = qres + rhoq(i)
            sumphiA = sumphiA +  rhopolA(i)
            sumphiB = sumphiB +  rhopolB(i)
            sumphiC = sumphiC +  rhopolC(i)


!            do j=1,nz 
!                FEVdWC = FEVdWC + deltaG(i)*rhopolC(i)* rhopolC(j)*chis(i,j)
!                FEVdWB = FEVdWB + deltaG(i)*rhopolB(i)* rhopolB(j)*chis(i,j)       
!            enddo   

        enddo
    
        FEel  = (volcell/vsol)*FEel
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        
        FEbind = volcell*FEbind/2.0_dp !  check this 

        qres = (volcell/vsol)*qres
        sumphiA = volcell*sumphiA
        sumphiB = volcell*sumphiB
        sumphiC = volcell*sumphiC

!        FEVdWC  = deltavoll*FEVdWC*VdWepsC*vpolC*vsol/2.0_dp   
!        FEVdWB  = deltavoll*FEVdWB*VdWepsB*vpolB(3)*vsol/2.0_dp   
    
!        if (sysflag=="elect") then 
!            FEVdW=FEVdWC
!        elseif (sysflag=="neutral") then
!            FEVdW=FEVdWB
!        else 
!            print*,"Wrong value sysflag : ", sysflag
!            stop    
!        endif   

      
        ! if((sigmaABL > sigmaTOL).and.(sigmaABR > sigmaTOL).and.(sigmaC > sigmaTOL)) then
        
        !     FEq =-delta*(sigmaABL*log(qABL)+ sigmaABR*log(qABR) +sigmaC*log(qC) )
       
        ! elseif((sigmaABL <= sigmaTOL).and.(sigmaABR <= sigmaTOL).and.(sigmaC > sigmaTOL)) then
        
        !     FEq = -delta*(sigmaC*log(qC) )
            
        ! elseif((sigmaABL > sigmaTOL).and.(sigmaABR <= sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
        !     FEq = -delta*(sigmaABL*log(qABL) )
       
        ! elseif((sigmaABL > sigmaTOL).and.(sigmaABR > sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
        !     FEq = -delta*(sigmaABL*log(qABL) +sigmaABR*log(qABR))
       
        ! elseif((sigmaABL <= sigmaTOL).and.(sigmaC <= sigmaTOL)) then
        
        !     FEq = 0.0_dp
        
        ! else
        
        !     print*,"Error in fcnerergy"
        !     print*,"Something went wrong in evaluating FEq"   
        !     print*,"sigmaABL=",sigmaABL
        !     print*,"sigmaABL=",sigmaABL
        !     print*,"sigmaC=",sigmaC
        !     stop    
        
        ! endif
    
        ! .. surface charge constribution 

        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL

        do s=1,nx*ny
            sigmaqSurf(RIGHT,s) = sigmaqSurfR(s)
            sigmaqSurf(LEFT,s)  = sigmaqSurfL(s)
            psiSurf(RIGHT,s)    = psiSurfR(s)
            psiSurf(LEFT,s)     = psiSurfL(s)
        enddo    
      
        do i = 1,2
            FEelsurf(i)=0.0_dp
            do s = 1, nx*ny    
                sigmaq0Surf(i,s)=  sigmaqSurf(i,s)/(delta*4.0_dp*pi*lb) ! dimensional charge density  
                FEelsurf(i) = FEelsurf(i)+sigmaq0Surf(i,s) * psiSurf(i,s) /2.0_dp 
            enddo   
        enddo    

        if(bcflag(RIGHT)=='qu') then ! quartz

            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(2)) ! area surface integration measure 
            enddo    
            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay
        
            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+(log(fdisS(2))+qS(2)*psiSurfR(s)) 
            enddo    

            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurf(RIGHT) =(log(fdisS(2))+dlog(fdisS(5)))*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
        
            FEchemSurf(RIGHT)= (log(fdisTaR(2))*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
        
            FEchemSurf(RIGHT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 

            FEchemSurf(LEFT)= log(fdisTaL(2))*sigmaSurf(LEFT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        
        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurf(LEFT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 


        !     .. total free energy per area of surface 

        FE = FEq  + FEpi + FErho + FEel + FEelSurf(RIGHT) + FEelSurf(LEFT)
        FE = FE + FEchemSurf(RIGHT)+FEchemSurf(LEFT) - FEVdW + FEbind
        
!        print*,"FE = " ,FE
        
        do i=LEFT,RIGHT
            qsurf(i)=0.0_dp
            do s=1,nx*ny     
                qsurf(i) = qsurf(i)+sigmaqSurf(i,s)
            enddo
            qsurf(i)=(qsurf(i)/(4.0_dp*pi*lb*delta))*delta*delta  ! delta*delta=area size one surface element, 4*pi*lb*delta make correct dimensional full unit    
        enddo

!        print*,"qsurf(LEFT)=",qsurf(LEFT),"qsurf(RIGHT)=",qsurf(RIGHT),"qres=",qres    

        qres = qres + (qsurf(RIGHT)+qsurf(LEFT))  ! total residual charge 

        
        volumelat= volcell*nsize !nz*delta   ! volume lattice

        FEbulk   = log(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ & 
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
        use conform_entropy


        implicit none

        !  .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf(2)           ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface  
        integer :: i,j,s               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer :: nzadius
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,nx*ny),sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: diffFEchemTa, FEchemSurftmp

        sigmaSurf(RIGHT)  = sigmaSurfR 
        sigmaSurf(LEFT)   = sigmaSurfL
        
        do s=1,nx*ny
            sigmaqSurf(RIGHT,s) = sigmaqSurfR(s)
            sigmaqSurf(LEFT,s)  = sigmaqSurfL(s)
            psiSurf(RIGHT,s)    = psiSurfR(s)
            psiSurf(LEFT,s)     = psiSurfL(s)
        enddo    


        !  .. computation of free energy 
    
        !  .. alternative computation free energy

        call FEconf_entropy(FEconfAB,FEconfC)


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

        FEchempot%sol   = 0.0_dp ! by construction  
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
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(1))+qS(1)*psiSurfR(s)
            enddo

            FEchemSurfalt(RIGHT) = (FEchemSurftmp)*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)

        elseif(bcflag(RIGHT)=="cl" ) then  ! clay        
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(1))+qS(1)*psiSurfR(s)
            enddo

            FEchemSurfalt(RIGHT) = (FEchemSurftmp)*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)

        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurfalt(RIGHT) =(log(fdisS(2))+log(fdisS(5)))*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisTaR(1))+qTA(1)*psiSurfR(s)
            enddo

            FEchemSurfalt(RIGHT)= (FEchemSurftmp*sigmaSurf(RIGHT)/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
            
            FEchemSurfalt(RIGHT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisTaL(1))+qTA(1)*psiSurfL(s)
            enddo
            FEchemSurfalt(LEFT)= FEchemSurftmp*sigmaSurf(LEFT)/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        

        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurfalt(LEFT)=0.0_dp
        
        else
        
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        
        endif 



        ! .. summing all contrubutions
        
        FEalt = FEtrans%sol +FEtrans%Na+ FEtrans%Cl +FEtrans%NaCl+FEtrans%Ca 
        FEalt = FEalt+FEtrans%OHmin +FEtrans%Hplus +FEtrans%K +FEtrans%KCl
        FEalt = FEalt+FEchempot%sol +FEchempot%Na+ FEchempot%Cl +FEchempot%NaCl+FEchempot%Ca 
        FEalt = FEalt+FEchempot%OHmin +FEchempot%Hplus+ FEchempot%K +FEchempot%K+FEchempot%KCl
        ! be vary carefull FE = -1/2 \int dz rho_q(z) psi(z)

         ! .. chemical and binding contribution
        FEchem = FEchem_react()
        
        FEalt = FEalt + FEconfAB + FEchem
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

        FEchempotbulk%sol   = 0.0_dp ! by construction  
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
            diffFEchemsurf(LEFT)= (sigmaSurf(LEFT)/(4.0_dp*pi*lb*delta))*( dlog(K0Ta(1)) -dlog(expmu%Hplus))
        else
            diffFEchemsurf(LEFT)=0.0_dp
        endif

        if(bcflag(RIGHT)=='qu') then ! quartz
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0_dp*pi*lb*delta))*( dlog(K0S(1)) -dlog(expmu%Hplus))
        elseif(bcflag(RIGHT)=='cl') then ! quartz
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0_dp*pi*lb*delta))*( dlog(K0S(1)) -dlog(expmu%Hplus))
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine       
            diffFEchemsurf(RIGHT)= (sigmaSurf(RIGHT)/(4.0_dp*pi*lb*delta))*( dlog(K0Ta(1)) -dlog(expmu%Hplus))
        else
            diffFEchemsurf(RIGHT)=0.0_dp
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
    
        use globals, only : nsize
        use volume, only : volcell
        use parameters
        use field
        use VdW
        use conform_entropy

        implicit none

        !     .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf              ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface 
        integer :: i,j               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
    !    integer :: nzadius

        !     .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        sumphiA = 0.0_dp
        sumphiB = 0.0_dp
        sumphiC = 0.0_dp

        FEq = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdW = 0.0_dp 
        qres = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + deltaG(i)*dlog(xsol(i))
            FErho = FErho - deltaG(i)*xsol(i) 

            do j=1,nz 
                FEVdW = FEVdW + deltaG(i)*rhopolB(i)* rhopolB(j)*chis(i,j)       
            enddo   

            sumphiA = sumphiA + deltaG(i) * rhopolA(i)
            sumphiB = sumphiB + deltaG(i) * rhopolB(i)
            sumphiC = sumphiC + deltaG(i) * rhopolC(i)
        enddo
    
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
    
        sumphiA = volcell*sumphiA
        sumphiB = volcell*sumphiB
        sumphiC = volcell*sumphiC

        FEVdW  = volcell*FEVdW*(VdWepsB*vpolB(3)*vsol)/2.0_dp   

!        FEq = -delta*(sigmaAB*dlog(qAB)+sigmaC*dlog(qC) )
    
        if(sigmaAB <= sigmaTOL) then 
!            FEq = -delta*sigmaC*dlog(qC)
        elseif (sigmaC <= sigmaTOL) then 
!            FEq = -delta*sigmaAB*dlog(qAB) 
        endif
    
        FEel = 0.0_dp    
        FEelsurf =0.0_dp
        FEchem =0.0_dp
        FEbind =0.0_dp
        qres =0.0_dp
        qsurf =0.0_dp 

        !  .. total free energy per area of surface 

        FE = FEq + FEpi + FErho - FEVdW 
        
        print*,"FE =",FE
    
       
        volumelat= volcell*nsize  ! nz*delta  ! volume lattice divide by area surface
        
        print*,"volumelat=",volumelat
        FEbulk  = log(xbulk%sol)-xbulk%sol
        print*,"FEbulk=",FEbulk
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE  = FE -FEbulk
    
        ! altnative computation free energy

        call FEconf_entropy(FEconfAB,FEconfC)

        FEtrans%sol=FEtrans_entropy(xsol,xbulk%sol,vsol,"w")     
        
        FEalt= FEconfAB+FEconfC+FEtrans%sol+FEVdW
        print*,"FEalt = ",FEalt
    
    end subroutine fcnenergy_neutral


    real(dp) function FEtrans_entropy(xvol,xvolbulk,vol,flag)
    
        use globals, only : nsize
        use parameters, only : vsol 
        use volume, only : volcell
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        integer :: i


        if(xvolbulk==0.0_dp) then 
            FEtrans_entropy=0.0_dp
        else
            FEtrans_entropy=0.0_dp
            if(present(flag)) then
            ! water special case because vsol treated different then vi  
                do i=1,nsize
                    FEtrans_entropy=FEtrans_entropy + xvol(i)*(log(xvol(i))-1.0_dp)
                enddo 
                FEtrans_entropy = volcell*FEtrans_entropy/vol
            else 
                do i=1,nsize
                    FEtrans_entropy = FEtrans_entropy + xvol(i)*(log(xvol(i)/vol)-1.0_dp)
                enddo
                FEtrans_entropy = volcell*FEtrans_entropy/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy


    real(dp) function FEchem_pot(xvol,expchempot,vol,flag)
    
        use globals, only : nsize
        use volume, only : volcell
        use parameters, only : vsol
        implicit none

        real(dp), intent(in) :: xvol(nsize)
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        ! .. local 
        integer :: i
        real(dp) :: chempot ! chemical potential difference 
        real(dp) :: sumdens 

        if(expchempot==0.0_dp) then 
            FEchem_pot=0.0_dp
        else     
            sumdens=0.0_dp
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                chempot = -log(expchempot)    
                do i=1,nsize
                    sumdens=sumdens +xvol(i)
                enddo
                FEchem_pot=volcell*chempot*sumdens/vol            
            else
                chempot = -log(expchempot/vol)    
                do i=1,nsize
                    sumdens=sumdens +xvol(i)
                enddo
                FEchem_pot=volcell*chempot*sumdens/(vol*vsol)               
            endif
        endif

    end function FEchem_pot


    real(dp) function FEtrans_entropy_bulk(xvolbulk,vol,flag)
    
        use globals
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        integer :: i


        if(xvolbulk==0.0_dp) then 
            FEtrans_entropy_bulk=0.0_dp
        else
            FEtrans_entropy_bulk=0.0_dp
            if(present(flag)) then
                ! water special case because vsol treated diffetent then vi  
                FEtrans_entropy_bulk=xvolbulk*(log(xvolbulk)-1.0_dp)/vol
            else 
                FEtrans_entropy_bulk=xvolbulk*(log(xvolbulk/vol)-1.0_dp)/(vol*vsol)
            endif
        endif

    end function FEtrans_entropy_bulk

    real(dp) function FEchem_pot_bulk(xvolbulk,expchempot,vol,flag)
    
        use globals
        use field
        use parameters
        implicit none

        real(dp), intent(in) :: xvolbulk
        real(dp), intent(in) :: expchempot 
        real(dp), intent(in) :: vol    
        character(len=1), optional :: flag    

        ! .. local 
        integer :: i
        real(dp) :: chempot ! chemical potential difference 
        real(dp) :: sumdens 

        if(expchempot==0.0_dp) then 
            FEchem_pot_bulk=0.0_dp
        else     
            if(present(flag)) then  ! water special case because vsol treated diffetent then vi
                FEchem_pot_bulk=-log(expchempot)*xvolbulk/vol            
            else
                FEchem_pot_bulk=-log(expchempot/vol)*xvolbulk/(vol*vsol)               
            endif
        endif

    end function FEchem_pot_bulk  


    real(dp) function FEchem_react()

        use globals, only : sysflag, nsize
        use field
        use volume
        use parameters, only : vpolA,vsol,zpolA,vpolB,zpolB

        implicit none

        integer :: i, k
        real(dp) :: lambdaA, lambdaB, rhopolAq, rhopolBq, xpolA, xpolB
        real(dp) :: betapi

        FEchem_react = 0.0_dp

        do i=1,nsize

            betapi=-log(xsol(i))/vsol

            lambdaA=-log(fdisA(1,i)) -psi(i)*zpolA(1)-betapi*vpolA(1)*vsol
            lambdaB=-log(fdisB(1,i)) -psi(i)*zpolB(1)-betapi*vpolB(1)*vsol

            rhopolAq = 0.0_dp
            rhopolBq = 0.0_dp
            xpolA=0.0_dp
            xpolB=0.0_dp

            do k=1,4
                rhopolAq=rhopolAq+ zpolA(k)*fdisA(k,i)*rhopolA(i)
                xpolA   =xpolA + rhopolA(i)*fdisA(k,i)*vpolA(k)*vsol
                rhopolBq=rhopolBq+ zpolB(k)*fdisB(k,i)*rhopolB(i)
                xpolB   =xpolB + rhopolB(i)*fdisB(k,i)*vpolB(k)*vsol
            enddo   
            xpolA = xpolA +rhopolA(i)*fdisA(5,i)*vpolA(5)*vsol/2.0_dp
            xpolB = xpolB +rhopolB(i)*fdisB(5,i)*vpolB(5)*vsol/2.0_dp
            
            FEchem_react = FEchem_react + (- rhopolA(i)*lambdaA -psi(i)*rhopolAq -betapi*xpolA &
                +fdisA(5,i)*rhopolA(i)/2.0_dp)

            FEchem_react = FEchem_react + (- rhopolB(i)*lambdaB -psi(i)*rhopolBq -betapi*xpolB &
                +fdisB(5,i)*rhopolB(i)/2.0_dp)

        enddo

        FEchem_react=volcell*FEChem_react    
        
        if(sysflag=="neutral") FEchem_react=0.0_dp


    end function FEchem_react


end module energy
