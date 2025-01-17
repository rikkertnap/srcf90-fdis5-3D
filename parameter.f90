 module parameters

    use physconst
    use mathconst
    !use random
    use volume
    use molecules
    use loopvar

    implicit none

    !  .. list of parameters

    type(moleclist) :: xbulk,expmu
    type(dipolelist):: dipole

    !  .. volume 
  
    real(dp) :: vsol               ! volume of solvent  in nm^3       
    real(dp) :: vpolB(5)           ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolA(5)           ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolC              ! volume of one polymer segment hydrocarbon, vpol  in units of vsol
  
    real(dp) :: deltavA(4)
    real(dp) :: deltavB(4)
  
    real(dp) :: vNa                ! volume positive ion in units of vsol
    real(dp) :: vK                 ! volume positive ion in units of vsol
    real(dp) :: vCl                ! volume negative ion in units of vsol   
    real(dp) :: vCa                ! volume positive divalent ion in units of vsol
    real(dp) :: vNaCl
    real(dp) :: vKCl
  
    !  .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RCl
    real(dp) :: RCa
    real(dp) :: RNaCl
    real(dp) :: RKCl
    
    real(dp) :: lsegAB
    real(dp) :: lsegA              ! segment length of A polymer in nm
    real(dp) :: lsegB              ! segment length of B polymer in nm
    real(dp) :: lsegC              ! segment length of C polymer in nm
    real(dp) :: lsegCH2 
    real(dp) :: lsegPAA   
    real(dp) :: lsegPEG  
    real(dp) :: lsegPAMPS
    real(dp) :: lsegPS  
    
    integer :: period            ! peridociy of repeat of A or B block 
  
    type (looplist), target :: VdWeps            ! strenght VdW interaction in units of kT
    type (looplist), target :: VdWepsAA, VdWepsBB,VdWepsAB            ! strenght VdW interaction in units of kT
    logical :: isVdW  

    !real(dp) :: VdWepsB            ! strenght VdW interaction in units of kT
    !real(dp) :: VdWepsC            ! strenght VdW interaction in units of kT
   
    integer :: VdWcutoff         ! cutoff VdW interaction in units of lseg 	
    integer :: VdWcutoffdelta    ! cutoff VdW interaction in units of delta
    integer :: layeroffset
  

    integer :: zpolA(5)          ! valence charge polymer
    integer :: zpolB(5)          ! valence charge polymer
    integer :: zNa               ! valence charge positive ion 
    integer :: zK                ! valence charge positive ion 
    integer :: zCa               ! valence charge divalent positive ion 
    integer :: zCl               ! valence charge negative ion 
    integer :: zsurf             ! valence surface charge 
  
    real(dp) :: Tref             ! temperature in K
    real(dp) :: dielectW         ! dielectric constant of water 
    real(dp) :: lb               ! Bjerrum length	   
    real(dp) :: constqW          ! constant in Poisson eq dielectric constant of water 
    real(dp) :: lb0              ! Bjerrum length in vacuum       
    real(dp) :: constq0          ! constant in Poisson eq dielectric constant of vacuum 


    real(dp) :: sigmaAB            ! sigma AB polymer coated on planar surface
    real(dp) :: sigmaABL           ! sigma AB polymer coated on planar surface
    real(dp) :: sigmaABR           ! sigma AB polymer coated on planar surface
    
    real(dp) :: sigmaC             ! sigma C polymer coated on planar surface
  
    integer :: itmax             ! maximum number of iterations
    real(dp) :: error              ! error imposed accuaracy
    real(dp) :: fnorm              ! L2 norm of residual vector function fcn  
    integer :: infile            ! infile=1 read input files infile!=1 no input files 
    integer :: iter              ! counts number of iterations
  
    character(len=8) :: method           ! method="kinsol" or "zspow"  
    character(len=8) :: chainmethod      ! method of generating chains ="MC" or "FILE" 
    character(len=8) :: chaintype        ! type of chain: diblock,alt
    integer :: readinchains              ! nunmber of used/readin chains
    character(len=3) ::  verboseflag       ! select input falg 

    real(dp) :: heightAB           ! average height of layer
    real(dp) :: heightC            ! average height of layer 
    real(dp) :: qpolA              ! charge poly A of layer 
    real(dp) :: qpolB              ! charge poly B of layer 
    real(dp) :: qpol_tot           ! charge poly A+B of layer 
    real(dp) :: avfdisA(5)         ! average degree of dissociation 
    real(dp) :: avfdisB(5)         ! average degree of dissociation
  
  !     .. weak polyelectrolyte variables 
  !     .. equibrium constant
  
    real(dp) :: K0A(4)              ! intrinsic equilibruim constant
    real(dp) :: KA(4)               ! experimemtal equilibruim constant 
    real(dp) :: pKA(4)              ! experimental equilibruim constant pKa= -log[Ka]
    real(dp) :: K0B(4)              ! intrinsic equilibruim constant
    real(dp) :: KB(4)               ! experimemtal equilibruim constant 
    real(dp) :: pKB(4)              ! experimental equilibruim constant pKa= -log[Ka]
  
    real(dp) :: pKw                 ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-] 
  
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
  
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
  
    !     .. bulk volume fractions 
  
    real(dp) :: pibulk             ! -ln(xsolbulk)
    real(dp) :: xNaClsalt          ! volume fraction of salt in bulk
    real(dp) :: xKClsalt           ! volume fraction of salt in bulk
    real(dp) :: xCaCl2salt         ! volume fraction of divalent salt in bulk
  
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp) :: cNaCl              ! concentration of salt in bulk in mol/liter
    real(dp) :: cKCl               ! concentration of salt in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of divalent salt in bulk in mol/liter
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
  
    type (looplist), target :: pH
        
contains

    ! determine total number of non linear equations

    subroutine set_size_neq()

        use globals
        use volume, only : nx, ny, nz

        implicit none

        integer(8) :: neq_bc

        nsize= nx*ny*nz
        neq_bc=0 
        if(bcflag(LEFT)/="cc") neq_bc=neq_bc+nx*ny
        if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+nx*ny

        select case (systype)
            case ("dipolarweak") 
                neq = 2 * nsize + neq_bc
            case ("dipolarstrong")  
                neq = 2 * nsize + neq_bc 
            case ("dipolarweakA") 
                neq = 2 * nsize + neq_bc
            case ("elect") 
                neq = 4 * nsize + neq_bc
            case ("electA") 
                neq = 3 * nsize + neq_bc
            case("electVdWAB")             ! copolymer weak polyacid, VdW
                neq = 4 * nsize + neq_bc
            case ("electdouble")  
                neq = 4 * nsize 
            case ("electnopoly") 
                neq = 2 * nsize + neq_bc
            case ("dipolarnopoly") 
                neq = 2 * nsize + neq_bc
            case ("neutral") 
                neq = nsize
            case ("bulk water") 
                neq = 5 
            case default
                print*,"Wrong value systype:  ",systype
                stop
        end select  

        neqint =neq ! used for MPI func binding, MPI has no integer(8)
         
    end subroutine set_size_neq

    
    real(dp) function BjerrumLenght(T)

        use mathconst
        use physconst
        
        implicit none
        
        real(dp), intent(in) :: T     
        real(dp) :: lb

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*T) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp              ! bjerrum length in water in nm
        BjerrumLenght=lb

    end function BjerrumLenght
        

    ! unit of dipole moment need to be in e unit of unit charge X nm
    
    subroutine init_dipolesmoment(dipolemoment)

        type(dipolelist), intent(inout) :: dipolemoment

        real(dp) :: unit_Debye_to_enm=0.020819434_dp  ! 1D= 0.020819434 e·nm 

        dipolemoment%sol=4.8_dp*unit_Debye_to_enm ! 1.8546_dp*unit_Debye_to_enm !
        dipolemoment%pol=1.7_dp*unit_Debye_to_enm !
         
    end subroutine


    !     purpose: initialize all constants parameter 
    !     pre: first read_inputfile has to be called   
    
    subroutine init_constants()

        use globals
        use volume
        use random
        use physconst
        
        implicit none      
        
        real(dp) :: vA,vB, vAA, vAMPS, vPEG
        
        !  .. initializations of variables
 
        pi=acos(-1.0_dp)          ! pi = arccos(-1)
        itmax=2000                ! maximum number of iterations      
        !nz=nzmax

        !     .. charges
        
        zsurf =-1                 ! valence surface charge 
        zNa   = 1                 ! valence positive charged ion
        zK    = 1                 ! valence positive charged ion
        zCa   = 2                 ! valence divalent positive charged ion
        zCl   =-1                 ! valence negative charged ion
        
        zpolA(1)=-1 ! A-
        zpolA(2)= 0 ! AH
        zpolA(3)= 0 ! ANa
        zpolA(4)= 1 ! ACa+
        zpolA(5)= 0 ! A2Ca
        
        zpolB(1)=-1 ! B-
        zpolB(2)= 0 ! BH
        zpolB(3)= 0 ! BNa
        zpolB(4)= 1 ! BCa+
        zpolB(5)= 0 ! B2Ca
        
        !     .. radii
        
        RNa = 0.102_dp             ! radius of Na+ in nm
        RK  = 0.138_dp             ! radius of K+ in nm
        RCl = 0.181_dp             ! radius of Cl- in nm
        RCa = 0.106_dp             ! radius of Ca2+ in nm
        RNaCl = 0.26_dp            ! radius of ion pair: this value is strange  and not used !!
        RKCl = 0.26_dp             ! radius of ion pair
        
        !     .. volume
        if(systype/="neutral") then 
            vsol = 0.030_dp              ! volume water solvent molecule in (nm)^3
        elseif(systype=="neutral") then 
            vsol  = 0.030_dp
           ! vsol = 0.218_dp             ! volume hexane Mw=86.18 g/mol and rho=0.6548 g/ml  
        else 
            print*,"Error in call to init_constants subroutine"   
            print*,"Wrong system flag"
            stop
        endif   

        vNa  = ((4.0_dp/3.0_dp)*pi*(RNa)**3)/vsol 
        vK   = ((4.0_dp/3.0_dp)*pi*(RK)**3)/vsol 
        vCl  = ((4.0_dp/3.0_dp)*pi*(RCl)**3)/vsol 
        vCa  = ((4.0_dp/3.0_dp)*pi*(RCa)**3)/vsol 
        vNaCl= (vNa+vCl)          ! contact ion pair
        vKCl = (vK+vCl)           ! contact ion pair
        
        !     .. volume polymer segments
        !     .. all volume scaled by vsol
        
        vAA  =  0.07448_dp/vsol ! volume based on VdW radii 
        vAMPS = 0.2134_dp/vsol
        vPEG  = 0.065_dp/vsol

        vA = vAA
        vB = vAMPS

        vpolA(1) = vA              ! vA-
        vpolA(2) = vA              ! vAH
        vpolA(3) = vA+vNa          ! vANa
        vpolA(4) = vA+vCa          ! vACa
        vpolA(5) = 2.0_dp*vA+vCa   ! vA2Ca
        
        vpolB(1) = vB              ! vB-
        vpolB(2) = vB              ! vBH
        vpolB(3) = vB+vNa          ! vBNa
        vpolB(4) = vB+vCa          ! vBCa
        vpolB(5) = 2.0_dp*vB+vCa   ! vB2Ca
        
        deltavA(1) = vpolA(1)+1.0_dp-vpolA(2) ! vA-+vH+-vAH
        deltavA(2) = vpolA(1)+vNa-vpolA(3)    ! vA-+vNa+-vANa+
        deltavA(3) = vpolA(1)+vCa-vpolA(4)    ! vA- + vCa2+ -vACa+
        deltavA(4) = 2.0_dp*vpolA(1)+vCa-vpolA(5) ! 2vA- + vCa2+ -vA2Ca
        
        deltavB(1) = vpolB(1)+1.0_dp-vpolB(2) ! vB-+vH+-vBH
        deltavB(2) = vpolB(1)+vNa-vpolB(3)    ! vB-+vNa+-vBNa+
        deltavB(3) = vpolB(1)+vCa-vpolB(4)    ! vB- +vCa2+ -vBCa+
        deltavB(4) = 2.0_dp*vpolB(1)+vCa-vpolB(5) ! 2vB- + vCa2+ -vB2Ca+
        
        vpolC = vPEG  ! volume PEG   !0.0270_dp/vsol     ! volume CH2
       
        !     .. other physical varaibles
        lsegPAA  = 0.36287_dp       ! segment length in nm
        lsegPAMPS = 0.545_dp        ! segment length in nm
        lsegCH2 =  0.153_dp         ! segment length in nm of CH2 check  
        lsegPEG =  0.3_dp           ! segment length in nm 
        
        lsegAB=lsegPAMPS          
        lsegA=lsegPAA            
        lsegB=lsegPAMPS          
        lsegC=lsegPEG               

        ! see also subroutine set_chain_properties 

        pKw=14.0_dp                 ! water equilibruim constant
        Tref=298.0_dp               ! temperature in Kelvin
        dielectW=78.54_dp           ! dielectric constant water
        
        seed  = 435672              ! seed for random number generator

        call init_elect_constants(Tref)  

        max_conforAB=cuantasAB
        max_conforC=cuantasC


    end subroutine init_constants


    function dielectric_constant_water(Temp) result(eps)
        implicit none
        real(dp), intent(in) :: Temp ! temperature in Kelvin
        real(dp) :: eps, Tc

        Tc=Temp-273.15_dp
        eps=(3.70886e4_dp - 8.2168e1_dp*Tc)/( 4.21854e2_dp + Tc)        !   Meissner and  Wentz 
        !eps= 87.7410_dp - 0.4008_dp*Tc + 9.398e-4_dp*Tc**2 - 1.410e-6_dp*Tc**3 ! Malmber and Maryott 

    end function

    subroutine init_elect_constants(Temp)
        
        use globals
        use volume, only : delta
    !    use physconst

        real(dp), intent(in) :: Temp 

        dielectW=dielectric_constant_water(Temp)

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*Temp) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp                           ! bjerrum length in water in nm
        constqW = delta*delta*(4.0_dp*pi*lb)/vsol ! multiplicative constant Poisson Eq. 

        if(systype.eq."dipolarweak".or.systype.eq."dipolarstrong".or.&
            systype.eq."dipolarweakA".or.systype.eq."dipolarpoly") then 

            lb=(elemcharge**2)/(4.0_dp*pi*dielect0*kBoltzmann*Temp)  ! bjerrum length in vacum in m
            lb= lb/1.0e-9_dp                                         ! bjerrum length in vacum in nm
            constqW = delta*delta*(4.0_dp*pi*lb)/vsol                ! multiplicative constant Poisson Eq. 

            lb0=(elemcharge**2)/(4.0_dp*pi*dielect0*kBoltzmann*Temp) ! bjerrum length in vacum in m
            lb0= lb0/1.0e-9_dp                                       ! bjerrum length in vacum in nm
            constq0 = delta*delta*(4.0_dp*pi*lb0)/vsol               ! multiplicative constant Poisson Eq. 

        endif     

    end subroutine init_elect_constants
   
    ! compute surface coverge based on number of grafted point (ngr) 
    subroutine init_sigma()

        use globals, only : systype
        use volume, only : ngr, nsurf, delta

        select case (systype) 
        case("elect") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = 0.0_dp
            sigmaAB  = sigmaABL
        case("electA") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = 0.0_dp
            sigmaAB  = sigmaABL    
        case("electdouble")
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = sigmaABL
            sigmaAB  = sigmaABL
        case('electnopoly')    
            sigmaABL = 0.0_dp
            sigmaABR = 0.0_dp
            sigmaAB  = 0.0_dp
        case('dipolarnopoly')    
            sigmaABL = 0.0_dp
            sigmaABR = 0.0_dp
            sigmaAB  = 0.0_dp    
        case("dipolarstrong") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = sigmaABL
            sigmaAB  = sigmaABL
        case("dipolarweak") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = sigmaABL
            sigmaAB  = sigmaABL
        case("dipolarweakA") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = 0.0_dp
            sigmaAB  = sigmaABL   
        case("neutral") 
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = 0.0_dp
            sigmaAB  = sigmaABL    
        case default
            print*,"Error: init_sigma: systype wrong value"
            print*,"stopping program"
            stop
        end select

    end subroutine
         
   
    !     purpose: initialize expmu needed by fcn 
    !     pre: first read_inputfile has to be called

    subroutine init_expmu_elect()
 
        use globals
        use physconst
        
        implicit none 
        
        !     .. local variable
        
        real(dp),  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
        real(dp),  dimension(:), allocatable :: xguess  
        integer :: i
        character(len=15) :: systype_old
        logical :: issolution
        
        allocate(x(5))
        allocate(xguess(5))
        
        !     .. initializations of input dependent variables, electrostatic part 
        
        pHbulk=pH%val ! transfer pH value 

        cHplus = (10.0_dp)**(-pHbulk) ! concentration H+ in bulk
        pOHbulk = pKw -pHbulk       
        cOHmin  = (10.0_dp)**(-pOHbulk) ! concentration OH- in bulk
        
        xbulk%Hplus = (cHplus*Na/(1.0e24_dp))*(vsol) ! volume fraction H+ in bulk vH+=vsol
        xbulk%OHmin = (cOHmin*Na/(1.0e24_dp))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
        
        xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
        
        if(pHbulk.le.7) then      ! pH<= 7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
        else                      ! pH >7
            xbulk%Na=xNaClsalt*vNa/(vNa+vCl) +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH  
            xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)  
        endif
        

        xKClsalt = (cKCl*Na/(1.0e24_dp))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
        xbulk%K = xKClsalt*vK/(vK+vCl)  
        xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl)  
        
        xCaCl2salt = (cCaCl2*Na/(1.0e24_dp))*((vCa+2.0_dp*vCl)*vsol) ! volume fraction CaCl2 in mol/l
        xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0_dp*vCl)
        xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0_dp*vCl/(vCa+2.0_dp*vCl)
        
        xbulk%NaCl=0.0_dp    ! no ion pairing
        xbulk%KCl=0.0_dp     ! no ion pairing
        
        xbulk%sol=1.0_dp -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca   
        
        !     .. if Kion neq 0 ion pairing !
        
        !     .. intrinstic equilibruim constant acid        
        !     Kion  = 0.246_dp ! unit 1/M= liter per mol !!!
        K0ionK  = KionK /(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        K0ionNa = KionNa/(vsol*Na/1.0e24_dp) ! intrinstic equilibruim constant 
        
        if((KionNa.ne.0.0_dp).or.(KionK.ne.0.0_dp)) then  
            systype_old=systype 
            systype="bulk water"        ! set solver to fcnbulk
            call set_size_neq()         ! number of nonlinear equations
            
            x(1)=xbulk%Na
            x(2)=xbulk%Cl
            x(3)=xbulk%NaCl
            x(4)=xbulk%K
            x(5)=xbulk%KCl
            
            xguess(1)=x(1)
            xguess(2)=x(2)
            xguess(3)=x(3)
            xguess(4)=x(4)
            xguess(5)=x(5)
           
            call solver(x, xguess, error, fnorm, issolution) 
            
            !     .. return solution
            
            xbulk%Na  =x(1)
            xbulk%Cl  =x(2)
            xbulk%NaCl=x(3)
            xbulk%K   =x(4)
            xbulk%KCl =x(5)

            ! reset of flags
            iter=0
            systype=systype_old         ! switch solver back
            call set_size_neq()         ! set number of non-linear equation  
            !call set_fcn()              ! set fcnptr to correct fcn        
            
            xbulk%sol=1.0_dp-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca 
            
        endif
         
        !     .. intrinstic equilibruim constants      
        do i=1,4
             Ka(i)  = 10.0_dp**(-pKa(i)) ! experimental equilibruim constant acid 
             Kb(i)  = 10.0_dp**(-pKb(i)) ! experimental equilibruim constant acid
             K0a(i) = (Ka(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
             K0b(i) = (Kb(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        !     .. rescale for i=4 2A- Ca <=> A2Ca
          
        K0a(4) = (K0a(4)*vsol)*(Na/1.0e24_dp)
        K0b(4) = (K0b(4)*vsol)*(Na/1.0e24_dp)
         

        pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
        ! exp(beta mu_i) = (rhobulk_i v_i) / exp(- beta pibulk v_i) 
        expmu%Na    = xbulk%Na   /(xbulk%sol**vNa) 
        expmu%K     = xbulk%K    /(xbulk%sol**vK)
        expmu%Ca    = xbulk%Ca   /(xbulk%sol**vCa) 
        expmu%Cl    = xbulk%Cl   /(xbulk%sol**vCl)
        expmu%NaCl  = xbulk%NaCl /(xbulk%sol**vNaCl)
        expmu%KCl   = xbulk%KCl  /(xbulk%sol**vKCl)
        expmu%Hplus = xbulk%Hplus/xbulk%sol ! vsol = vHplus 
        expmu%OHmin = xbulk%OHmin/xbulk%sol ! vsol = vOHmin 
          
        !     .. end init electrostatic part 
        

        deallocate(x)
        deallocate(xguess)
        
    end subroutine init_expmu_elect


    subroutine init_expmu_neutral

        implicit none
          
        xbulk%sol=1.0_dp ! only solvent 

        !VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
        !VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled

    end subroutine init_expmu_neutral

    subroutine init_expmu

        use globals, only : systype
        implicit none

        if(systype=="elect".or.systype=="electA") then 
            call init_expmu_elect()
        elseif(systype=="electdouble") then 
            call init_expmu_elect()
        elseif(systype=="electnopoly") then 
            call init_expmu_elect()
        elseif(systype=="neutral") then
            call init_expmu_neutral()
        elseif(systype=="dipolarstrong") then 
            call init_expmu_elect()
        elseif(systype=="dipolarweak") then 
            call init_expmu_elect()
        else
            print*,"Error in call to init_expmu subroutine"    
            print*,"Wrong value systype : ", systype
            stop        
        endif   

    end subroutine init_expmu


    ! inits chem potential 

    subroutine init_vars_input()

        use globals, only : systype
     
        implicit none
        
        ! local variable
        integer :: i

        select case (systype)
        case ("elect")
            call init_expmu_elect() 
        !    call init_elect_constants(T%val) 
        !    this call is not neccesary subroutine allready called in init_constants
        case ("electA")
            call init_expmu_elect() 
        case ("dipolarstrong")
            call init_expmu() 
            call init_dipolesmoment(dipole)
        case ("dipolarweak")
            call init_expmu() 
            call init_dipolesmoment(dipole)
         case ("dipolarweaA")
            call init_expmu() 
            call init_dipolesmoment(dipole)    
        case ("electdouble") 
            call init_expmu_elect() 
        case ("electnopoly")
            call init_expmu_elect()
        case ("dipolarnopoly")
            call init_expmu_elect()
        case ("neutral")
            call init_expmu_neutral()   
        case default   
            print*,"Error: systype incorrect at init_vars_input" 
            print*,"Wrong value systype : ", systype
            print*,"stopping program"
            stop
         end select
             
     end subroutine init_vars_input

 end module parameters
