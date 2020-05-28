
!  .. module file for free energy variables /calculations

module energy 

    use mpivars    
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
    real(dp) :: FEel                ! electrostatics energy
!    real(dp) :: FEelb               ! contrubution of bound charges to electrostatic energy

    real(dp) :: FEelsurf(2)         ! electrostatics energy from  surface
    real(dp) :: FEchemsurf(2)       ! chemical free energy surface
    real(dp) :: FEchem
    real(dp) :: FEbind,FEbindA,FEbindB    ! complexation contribution
    real(dp) :: FEVdW,FEVdWB,FEVdWC       ! Van der Waals contribution
    real(dp) :: FEconf
    real(dp) :: Econf
    real(dp) :: FEalt               ! free energy
    real(dp) :: FEbulkalt           ! free energybulk
    real(dp) :: deltaFEalt          ! free energy difference delteFE=FE-FEbulk

    real(dp) :: FEchemsurfalt(2)    ! chemical free energy surface
    real(dp) :: diffFEchemsurf(2)   ! difference cheme

    type(moleclist) :: FEtrans,FEchempot,FEtransbulk,FEchempotbulk
    type(moleclist) :: deltaFEtrans,deltaFEchempot

    real(dp), dimension(:), allocatable :: sumphi  ! check integral over phiA 
    real(dp) :: sumphiA             ! check integral over phiA
    real(dp) :: sumphiB             ! check integral over phiB
    real(dp) :: qres                ! charge charge
    real(dp) :: checkphi            ! check integrate over phi
    
!   real(dp), parameter :: sigmaTOL = 0.00000001_dp     ! tolerance of surface coverage below no polymers 

!    private :: sigmaTOL

contains

    subroutine fcnenergy()
 
        use globals, only : systype 
        use myutils

        implicit none
        
        character(len=lenText) :: text 

        select case (systype) 
        case ("brush_mul","brushssdna")
        
            call fcnenergy_electbrush_mul() 
            call fcnenergy_elect_alternative()

        case("elect","electdouble","electnopoly","electA")  
            
            call fcnenergy_elect()
            call fcnenergy_elect_alternative()
        
        case("neutral")
            
            call fcnenergy_neutral()
        
        case default  

            text="fcnenergy: wrong systype: "//systype//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        
        end select 
  
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
        integer  :: i,j,s,g               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: nzadius
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx),sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: FEchemSurftmp
        integer, parameter :: A=1, B=2    

        !  .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        sumphiA = 0.0_dp
        sumphiB = 0.0_dp
        FEq     = 0.0_dp
        FEbindA = 0.0_dp
        FEbindB = 0.0_dp
        FEchem = 0.0_dp
        FEVdWB = 0.0_dp     
        qres   = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel  = FEel  - rhoq(i) * psi(i)/2.0_dp      
            FEbindA = FEbindA + fdisA(i,5)*rhopol(i,A)
            FEbindB = FEbindB + fdisB(i,5)*rhopol(i,B)

            qres = qres + rhoq(i)
            sumphiA = sumphiA +  rhopol(i,A)
            sumphiB = sumphiB +  rhopol(i,B)
        enddo

        ! .. calcualtion of FEVdW
        if(isVdW) then 
            if(systype=="electA") then
                FEVdW=-VdW_energy_homo(rhopol(:,A))
            else if(systype=="electVdWAB") then
                FEVdW =-VdW_energy_diblock(rhopol(:,A),rhopol(:,B))
            endif
        else   
            FEVdW=0.0_dp
        endif

        FEel  = (volcell/vsol)*FEel
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        FEbindA = volcell*FEbindA/2.0_dp !  
        FEbindB = volcell*FEbindB/2.0_dp !
        qres = (volcell/vsol)*qres
        sumphiA = volcell*sumphiA
        sumphiB = volcell*sumphiB
        checkphi= nseg*ngr-sumphiA-sumphiB


        if(systype=="nopoly".or.systype=="neutral") then 
            FEbind = 0.0_dp 
        else if(systype=="electVdWpos") then 
            FEbind = FEbindA
        else 
            FEbind = FEbindA+FEbindB
        endif     
        
        ! .. calcualtion of FEq
        if(systype=="electnopoly".or.systype=="dipolarnopoly") then 
            FEq=0.0_dp
        elseif(systype=="elect") then 
            FEq=0.0_dp
            do g=1,ngr
                FEq=FEq-log(q(g))    
            enddo
        elseif(systype=="electA".or.systype=="electVdWAB") then 
            FEq=0.0_dp
            do g=1,ngr
                FEq=FEq-log(qABL(g))    
            enddo
        elseif(systype=="electdouble") then
            FEq=0.0_dp
            do g=1,ngr
                FEq=FEq-log(qABL(g)) -log(qABR(g))   
            enddo
        elseif(systype=="neutral") then
            FEq=0.0_dp
            do g=1,ngr
                FEq=FEq-log(qAB(g))  
            enddo
        else
            print*,"Error in fcnenergy"
            print*,"Calculation FEq failed, wrong systype : ",systype 
            stop   
        endif
            
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
            FEelsurf(i)=FEelsurf(i)*areacell
        enddo    

        if(bcflag(RIGHT)=='qu') then ! quartz

            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(2)) ! area surface integration measure 
            enddo    
            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay
        
            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+(log(fdisS(2))+qS(2)*psiSurfR(s)) 
            enddo    

            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurf(RIGHT) =(log(fdisS(2))+log(fdisS(5)))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) - &
                2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
        
            FEchemSurf(RIGHT)=(log(fdisTaR(2))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
        
            FEchemSurf(RIGHT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 

            FEchemSurf(LEFT)= log(fdisTaL(2))*sigmaSurf(LEFT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        
        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurf(LEFT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 

        !     .. total free energy per area of surface 

        FE = FEq  + FEpi + FErho + FEel + FEelSurf(RIGHT) + FEelSurf(LEFT)
        FE = FE + FEchemSurf(RIGHT)+FEchemSurf(LEFT) + FEVdW + FEbind
        
        do i=LEFT,RIGHT
            qsurf(i)=0.0_dp
            do s=1,nx*ny     
                qsurf(i) = qsurf(i)+sigmaqSurf(i,s)
            enddo
            qsurf(i)=(qsurf(i)/(4.0_dp*pi*lb*delta))*areacell  ! areacell=area size one surface element, 4*pi*lb*delta make correct dimensional full unit    
        enddo

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

        !  .. computation ofalternative computation free energy

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

            FEchemSurfalt(RIGHT) = (FEchemSurftmp)*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)

        elseif(bcflag(RIGHT)=="cl" ) then  ! clay        
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(1))+qS(1)*psiSurfR(s)
            enddo

            FEchemSurfalt(RIGHT) = (FEchemSurftmp)*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)

        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurfalt(RIGHT) =(log(fdisS(2))+log(fdisS(5)))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) - &
                2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
            FEchemSurftmp=0.0_dp
            do s=1,nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisTaR(1))+qTA(1)*psiSurfR(s)
            enddo

            FEchemSurfalt(RIGHT)=(FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
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
            FEchemSurfalt(LEFT)=FEchemSurftmp*sigmaSurf(LEFT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        

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
        if(systype/="brush_mul") then 
            FEchem = FEchem_react()
        else
            FEchem = FEchem_react_multi()
        endif    

        print*,"isVdW in fcnenergy_elect_alternative=",isVdW
        if(isVdW) FEalt = FEalt-FEVdW ! add Van der Waals
            
        FEalt = FEalt + FEconf + FEchem 
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

        volumelat = volcell*nsize   ! volume lattice 
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
    

    subroutine fcnenergy_electbrush_mul()

        use globals
        use volume
        use parameters
        use field
        use VdW
        use surface

        !  .. local arguments 
    
        real(dp) :: sigmaq0,psi0
        real(dp) :: qsurf(2)           ! total charge on surface 
        real(dp) :: qsurfg             ! total charge on grafting surface  
        integer  :: i,j,s,g,t          ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
        integer  :: nzadius
        real(dp) :: sigmaSurf(2),sigmaqSurf(2,ny*nx),sigmaq0Surf(2,nx*ny),psiSurf(2,nx*ny)
        real(dp) :: FEchemSurftmp
        integer  :: ier
        logical  :: alloc_fail
        
        if (.not. allocated(sumphi))  then 
            allocate(sumphi(nsegtypes),stat=ier)
            if( ier/=0 ) alloc_fail=.true.
        endif

        !  .. computation of free energy

        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        sumphi = 0.0_dp
 
        FEq    = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdW  = 0.0_dp
        qres   = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - (xsol(i) + xHplus(i) + xOHmin(i)+ xNa(i)/vNa + xCa(i)/vCa + xCl(i)/vCl+xK(i)/vK +&
                xNaCl(i)/vNaCl +xKCl(i)/vKCl)                 ! sum over  rho_i 
            FEel  = FEel  - rhoq(i) * psi(i)
            qres = qres + rhoq(i)
        enddo

        checkphi=nseg*ngr
        do t=1,nsegtypes
            sumphi(t)=0.0_dp
            do i=1,nsize    
                sumphi(t) = sumphi(t) + rhopol(i,t)
            enddo
            sumphi(t) = volcell*sumphi(t)
            checkphi = checkphi-sumphi(t)
            print*,"fcnenergy brush mul rank=",rank,"t=",t,"sum= ",sumphi(t),"check=",checkphi
        enddo

        FEel  = (volcell/vsol)*FEel/2.0_dp  ! carefully check this
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
        qres  = (volcell/vsol)*qres


        ! .. calcualtion of FEVdW
        if(isVdW) then 
            FEVdW=-VdW_energy(rhopol)
        else   
            FEVdW=0.0_dp
        endif

        !  .. calcium and magnesium binding contribution to minimized free energy  

        if(systype/="brush_mul") then 
            do i=1,nsize
                FEbind = FEbind + (fdisA(i,5)+fdisA(i,7))*rhopol(i,tA)
            enddo                
            FEbind = volcell*FEbind /2.0_dp                                 
        endif      

        
        ! .. calcualtion of FEq
        
        FEq=0.0_dp
        do g=1,ngr
            FEq=FEq-log(q(g))  
        enddo
            
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
            FEelsurf(i)=FEelsurf(i)*areacell
        enddo    

        if(bcflag(RIGHT)=='qu') then ! quartz

            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+log(fdisS(2)) ! area surface integration measure 
            enddo    
            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cl" ) then  ! clay
        
            FEchemSurftmp=0.0_dp
            do s=1, nx*ny
                FEchemSurftmp=FEchemSurftmp+(log(fdisS(2))+qS(2)*psiSurfR(s)) 
            enddo    

            FEchemSurf(RIGHT) = FEchemSurftmp*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ca" ) then ! calcite
        
            FEchemSurf(RIGHT) =(log(fdisS(2))+log(fdisS(5)))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb) - &
                2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="ta" ) then ! taurine 
        
            FEchemSurf(RIGHT)=(log(fdisTaR(2))*sigmaSurf(RIGHT)*areacell/(delta*4.0_dp*pi*lb)) -2.0_dp*FEelsurf(RIGHT)
        
        elseif(bcflag(RIGHT)=="cc") then  
        
            FEchemSurf(RIGHT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(RIGHT) : ",bcflag(RIGHT)
            stop
        endif 

        if(bcflag(LEFT)=="ta" ) then ! taurine 

            FEchemSurf(LEFT)= log(fdisTaL(2))*sigmaSurf(LEFT)*areacell/(delta*4.0_dp*pi*lb) -2.0_dp*FEelsurf(LEFT)
        
        elseif(bcflag(LEFT)=="cc") then  
        
            FEchemSurf(LEFT)=0.0_dp
        
        else
            print*,"Error in fcnenergy"
            print*,"Wrong value bcflag(LEFT) : ",bcflag(LEFT)
        endif 


        !     .. total free energy per area of surface 

        FE = FEq  + FEpi + FErho + FEel + FEelSurf(RIGHT) + FEelSurf(LEFT)
        FE = FE + FEchemSurf(RIGHT)+FEchemSurf(LEFT) + FEVdW + FEbind
        
!        print*,"FE = " ,FE
        
        do i=LEFT,RIGHT
            qsurf(i)=0.0_dp
            do s=1,nx*ny     
                qsurf(i) = qsurf(i)+sigmaqSurf(i,s)
            enddo
            qsurf(i)=(qsurf(i)/(4.0_dp*pi*lb*delta))*areacell  ! areacell=area size one surface element, 4*pi*lb*delta make correct dimensional full unit    
        enddo

!        print*,"qsurf(LEFT)=",qsurf(LEFT),"qsurf(RIGHT)=",qsurf(RIGHT),"qres=",qres    

        qres = qres + (qsurf(RIGHT)+qsurf(LEFT))  ! total residual charge 

        
        volumelat= volcell*nsize !nz*delta   ! volume lattice

        FEbulk   = log(xbulk%sol)-(xbulk%sol+xbulk%Hplus +xbulk%OHmin+ & 
            xbulk%Na/vNa +xbulk%Ca/vCa +xbulk%Cl/vCl+ xbulk%K/vK + xbulk%NaCl/vNaCl +xbulk%KCl/vKCl )
        FEbulk = volumelat*FEbulk/(vsol)

        deltaFE = FE - FEbulk
    
    end subroutine fcnenergy_electbrush_mul




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
        integer  :: i,j,g               ! dummy variables 
        real(dp) :: volumelat          ! volume lattice 
    !    integer :: nzadius

        !     .. computation of free energy 
    
        FEpi  = 0.0_dp
        FErho = 0.0_dp
        FEel  = 0.0_dp
        FEelsurf = 0.0_dp
        sumphiA = 0.0_dp
        sumphiB = 0.0_dp

        FEq = 0.0_dp
        FEbind = 0.0_dp
        FEchem = 0.0_dp
        FEVdW = 0.0_dp 
        qres = 0.0_dp

        do i=1,nsize
            FEpi = FEpi  + log(xsol(i))
            FErho = FErho - xsol(i) 
            sumphiA = sumphiA +  rhopolA(i)
            sumphiB = sumphiB +  rhopolB(i)
        enddo
    
        FEpi  = (volcell/vsol)*FEpi
        FErho = (volcell/vsol)*FErho
    
        sumphiA = volcell*sumphiA
        sumphiB = volcell*sumphiB
        FEVdW=-VdW_energy_homo(rhopolB)
        
        do g=1,ngr
            FEq=FEq-log(qAB(g))  
        enddo

        FEel = 0.0_dp    
        FEelsurf =0.0_dp
        FEchem =0.0_dp
        FEbind =0.0_dp
        qres =0.0_dp
        qsurf =0.0_dp 

        !  .. total free energy per area of surface 

        FE = FEq + FEpi + FErho + FEVdW 
        
       
        volumelat= volcell*nsize  ! nz*delta  ! volume lattice divide by area surface
        FEbulk = log(xbulk%sol)-xbulk%sol
        FEbulk = volumelat*FEbulk/(vsol)
        deltaFE  = FE -FEbulk
    
        ! altnative computation free energy
        call FEconf_entropy(FEconf,Econf)

        FEtrans%sol=FEtrans_entropy(xsol,xbulk%sol,vsol,"w")     
        
        FEalt= FEconf+FEtrans%sol-FEVdW
    
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

        use globals, only : systype, nsize
        use field
        use volume, only : volcell
        use parameters, only : vpolA,vsol,zpolA,vpolB,zpolB

        implicit none

        integer :: i, k
        real(dp) :: lambdaA, lambdaB, rhopolAq, rhopolBq, xpolA, xpolB
        real(dp) :: betapi
        integer, parameter :: A=1, B=2


        select case (systype)
        case ("elect","electVdWAB","electdouble") 

            FEchem_react = 0.0_dp

            do i=1,nsize

                betapi=-log(xsol(i))/vsol

                lambdaA=-log(fdisA(i,1)) -psi(i)*zpolA(1)-betapi*vpolA(1)*vsol
                lambdaB=-log(fdisB(i,1)) -psi(i)*zpolB(1)-betapi*vpolB(1)*vsol

                rhopolAq = 0.0_dp
                rhopolBq = 0.0_dp
                xpolA=0.0_dp
                xpolB=0.0_dp

                do k=1,4
                    rhopolAq=rhopolAq+ zpolA(k)*fdisA(i,k)*rhopol(i,A)
                    xpolA   =xpolA + rhopol(i,A)*fdisA(i,k)*vpolA(k)*vsol
                    rhopolBq=rhopolBq+ zpolB(k)*fdisB(i,k)*rhopol(i,B)
                    xpolB   =xpolB + rhopol(i,B)*fdisB(i,k)*vpolB(k)*vsol
                enddo   
                xpolA = xpolA +rhopol(i,A)*fdisA(i,5)*vpolA(5)*vsol/2.0_dp
                xpolB = xpolB +rhopol(i,B)*fdisB(i,5)*vpolB(5)*vsol/2.0_dp
                
                FEchem_react = FEchem_react + (- rhopol(i,A)*lambdaA -psi(i)*rhopolAq -betapi*xpolA &
                    +fdisA(i,5)*rhopol(i,A)/2.0_dp)

                FEchem_react = FEchem_react + (- rhopol(i,B)*lambdaB -psi(i)*rhopolBq -betapi*xpolB &
                    +fdisB(i,5)*rhopol(i,B)/2.0_dp)

            enddo

            FEchem_react=volcell*FEChem_react    

        case("electA") 

            FEchem_react = 0.0_dp

            do i=1,nsize

                betapi=-log(xsol(i))/vsol

                lambdaA=-log(fdisA(i,1)) -psi(i)*zpolA(1)-betapi*vpolA(1)*vsol

                rhopolAq = 0.0_dp
                xpolA=0.0_dp

                do k=1,4
                    rhopolAq=rhopolAq+ zpolA(k)*fdisA(i,k)*rhopolA(i)
                    xpolA   =xpolA + rhopolA(i)*fdisA(i,k)*vpolA(k)*vsol
                enddo   
                xpolA = xpolA +rhopolA(i)*fdisA(i,5)*vpolA(5)*vsol/2.0_dp
                
                FEchem_react = FEchem_react + (- rhopolA(i)*lambdaA -psi(i)*rhopolAq -betapi*xpolA &
                    +fdisA(i,5)*rhopolA(i)/2.0_dp)
            enddo

            FEchem_react=volcell*FEChem_react 

        case("neutral","electnopoly") 
            FEchem_react=0.0_dp
        case default
            print*,"systype in FEchem_react_wrong"  
        end select

    end function FEchem_react

    function FEchem_react_multi()result(FEchem_react)

        use globals, only : systype,nsize,nsegtypes
        use field, only : xsol,psi,rhopol,fdis,fdisA !, epsfcn,Depsfcn
        use field, only : xNa,xCl,xHplus,xOHmin,xCa,xMg,xRb
        use volume, only : volcell
        use parameters, only : vsol, vpol,zpol,vpolAA,zpolAA !, constqE, lb, bornrad
        use parameters, only : vNa,vCl,vRb,vCa,vMg,zNa,zCl,zRb,zCa,zMg, tA
        use chains, only : ismonomer_chargeable
       ! use dielectric_const, only : born

        real(dp) :: FEchem_react

        integer :: i, k, t
        real(dp) :: lambda,  rhopolq, betapi, xpol, Eself ,bornene, lbr, Ebornself

        FEchem_react = 0.0_dp

        select case (systype)
        case("brush_mul") 
        
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then
                    do i=1,nsize

                        betapi=-log(xsol(i))/vsol
                        lambda=-log(fdis(i,t)) -psi(i)*zpol(t,2)-betapi*vpol(t)*vsol
                        rhopolq=zpol(t,2)*fdis(i,t)*rhopol(i,t)
            
                        FEchem_react = FEchem_react + &
                            (- rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*rhopol(i,t)*vpol(t)*vsol )
                    enddo
                endif        
            enddo
        
        case("brush","brush_neq") 

            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then

                    do i=1,nsize

                        betapi=-log(xsol(i))/vsol
                        lambda=-log(fdisA(i,1)) -psi(i)*zpolAA(1)-betapi*vpolAA(1)*vsol
                       
            
                        rhopolq = 0.0_dp
                        xpol=0.0_dp
                        do k=1,4
                            rhopolq=rhopolq+ zpolAA(k)*fdisA(i,k)*rhopol(i,t)
                            xpol  =xpol + rhopol(i,t)*fdisA(i,k)*vpolAA(k)*vsol
                        enddo   
                        rhopolq=rhopolq+ zpolAA(6)*fdisA(i,6)*rhopol(i,t)
                        xpol = xpol +rhopol(i,t)*(fdisA(i,5)*vpolAA(5)*vsol/2.0_dp +&
                                                  fdisA(i,6)*vpolAA(6)*vsol +&
                                                  fdisA(i,7)*vpolAA(7)*vsol/2.0_dp)

            
                        FEchem_react = FEchem_react + &
                            ( -rhopol(i,t)*lambda -psi(i)*rhopolq -betapi*xpol+(fdisA(i,5)+fdisA(i,7))*rhopol(i,t)/2.0_dp)
                       
                    enddo

                endif        
            enddo

        case default
            print*,"systype in FEchem_react_multi wrong"  
        end select      

        FEchem_react=volcell*FEChem_react    
        

    end function FEchem_react_multi


end module energy
