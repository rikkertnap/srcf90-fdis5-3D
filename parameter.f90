 module parameters

    !use physconst
    !use mathconst
    use volume
    use molecules
    use loopvar

    implicit none

    !  .. list of parameters

    type(moleclist) :: xbulk
    type(moleclist) :: expmu
    type(bornmoleclist) :: bornrad,bornbulk 

    !  .. volume 
  
    real(dp) :: vsol                 ! volume of solvent  in nm^3       
    real(dp) :: vpolA(5),deltavA(4)  ! volume of one polymer segment, vpol  in units of vsol
    real(dp) :: vpolB(5),deltavB(4)
    real(dp) :: vpolAA(7),deltavAA(6)
    real(dp), dimension(:), allocatable :: vpol  ! volume of polymer segment of given type, vpol in units of vsol
    
    real(dp) :: vNa                ! volume positive ion in units of vsol
    real(dp) :: vK                 ! volume positive ion in units of vsol
    real(dp) :: vRb                 ! volume positive ion in units of vsol
    real(dp) :: vCl                ! volume negative ion in units of vsol   
    real(dp) :: vCa                ! volume positive divalent ion in units of vsol
    real(dp) :: vMg                ! volume positive divalent ion in units of vsol
    real(dp) :: vNaCl
    real(dp) :: vKCl
  
    !  .. radii
  
    real(dp) :: RNa
    real(dp) :: RK
    real(dp) :: RRb
    real(dp) :: RCl
    real(dp) :: RCa
    real(dp) :: RMg
    
    real(dp), dimension(:,:), allocatable :: VdWeps    ! strenght VdW interaction in units of kT
    real(dp) :: VdWepsAA, VdWepsBB,VdWepsAB            ! strenght VdW interaction in units of kT
    logical :: isVdW  
    integer :: VdWcutoff         ! cutoff VdW interaction in units of lseg 	
    integer :: VdWcutoffdelta    ! cutoff VdW interaction in units of delta
    integer :: layeroffset
    
    !    .. charges 
    integer :: zpolAA(7)
    integer, dimension(:,:), allocatable :: zpol          ! valence charge polymer
    integer :: zpolA(5)          ! valence charge polymer
    integer :: zpolB(5)          ! valence charge polymer
    integer :: zNa               ! valence charge positive ion 
    integer :: zK                ! valence charge positive ion 
    integer :: zRb               ! valence charge positive ion 
    integer :: zCa               ! valence charge divalent positive ion 
    integer :: zMg               ! valence charge divalent positive ion 
    integer :: zCl               ! valence charge negative ion 
    !integer :: zsurf             ! valence surface charge 
  
     !  .. input filenames select if chainmethod==file
    integer, parameter :: lenfname=40
    character(len=lenfname) :: chainfname,vpolfname,pKafname,typesfname,lsegfname

    ! .. other physical quanties
  
    real(dp) :: Tref             ! temperature in K
    real(dp) :: dielectW         ! dielectric constant of water
    real(dp) :: dielectP            ! dielectric constant of hydrocarbons/PA
    character(len=15) :: dielect_env ! selects dielectric fun 
    real(dp) :: lb,lb0           ! Bjerrum lengtin water and vacuum   
    real(dp) :: constqW          ! constant in Poisson eq dielectric constant of water       
    real(dp) :: constq0          ! constant in Poisson eq dielectric constant of vacuum 
    real(dp) :: constqE          ! electrostatic pre-factor in pdf 

    real(dp) :: sigmaAB          ! sigma AB polymer coated on planar surface
    real(dp) :: sigmaABL         ! sigma AB polymer coated on planar surface
    real(dp) :: sigmaABR         ! sigma AB polymer coated on planar surface
    real(dp) :: sigma            ! sigma  polymer coated on planar surface
  
     !  .. solver variables

    integer  :: itmax            ! maximum number of iterations
    real(dp) :: tol_conv         ! error imposed accuaracy
    real(dp) :: fnorm            ! L2 norm of residual vector function fcn  
    integer  :: infile           ! infile=1 read input files infile!=1 no input files 
    integer  :: iter             ! counts number of iterations
    character(len=8) :: method   ! method="kinsol"  
    logical ::  precondition     ! controls use of precondtioner 

    !  .. output control
    character(len=3) ::  verboseflag    ! select input falg 

    ! .. chain variables 
    real(dp) :: lseg    ! segment length of A polymer in nm
    real(dp), dimension(:), allocatable :: lsegAA 
    real(dp) :: lsegA              ! segment length of A polymer in nm
    real(dp) :: lsegB              ! segment length of B polymer in nm
    real(dp) :: lsegPAA   
    real(dp) :: lsegPEG  
    real(dp) :: lsegPAMPS
    real(dp) :: lsegPS  

    character(len=15) :: chainmethod   ! method of generating chains ="MC" or "FILE" 
    character(len=8)  :: chaintype     ! type of chain: diblock,alt
    integer :: readinchains           ! nunmber of used/readin chains
    integer :: chainperiod            ! peridociy of repeat of A or B block 
    integer :: maxnchainsrotations    ! number of rotations read in from input.in, assigned to maxnchain in chaingenerator default 12  
    integer :: tA             ! segment number type of monomer type A  

    ! ..average structural properties of layer

    real(dp) :: height,heightAB    ! average height of layer 
    real(dp) :: qpolA              ! charge poly A of layer 
    real(dp) :: qpolB              ! charge poly B of layer 
    real(dp) :: qpol_tot           ! charge poly A+B of layer 
  
    real(dp), dimension(:), allocatable :: qpol                ! charge poly of layer 
    real(dp), dimension(:), allocatable :: avfdis              ! average degree of dissociation
    real(dp) :: avfdisA(7)         ! average degree of dissociation 
    real(dp) :: avfdisB(5)         ! average degree of dissociation

    !  .. weak polyelectrolyte variables 
    !  .. equibrium constant
    real(dp), dimension(:), allocatable :: K0a              ! intrinsic equilibruim constant
    real(dp), dimension(:), allocatable :: Ka               ! experimemtal equilibruim constant 
    real(dp), dimension(:), allocatable :: pKa              ! experimental equilibruim constant pKa= -log[Ka]
    
    real(dp) :: KaA(4),K0aA(4),pKaA(4)     !  ... constant for  acid 
    real(dp) :: KaB(4),K0aB(4),pKaB(4)   
    real(dp) :: KaAA(6),K0aAA(6),pKaAA(6)    
      
     ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-]   
    real(dp) :: pKw                 
  
    real(dp) :: K0ionNa             ! intrinsic equilibruim constant
    real(dp) :: KionNa              ! experimemtal equilibruim constant 
    real(dp) :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
  
    real(dp) :: K0ionK              ! intrinsic equilibruim constant
    real(dp) :: KionK               ! experimemtal equilibruim constant 
    real(dp) :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
  
    !     .. bulk volume fractions 
  
  !  real(dp) :: pibulk             ! -ln(xsolbulk)
  !  real(dp) :: xNaClsalt          ! volume fraction of salt in bulk
  !  real(dp) :: xKClsalt           ! volume fraction of salt in bulk
  !  real(dp) :: xCaCl2salt         ! volume fraction of divalent salt in bulk
  
    real(dp) :: cHplus             ! concentration of H+ in bulk in mol/liter
    real(dp) :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real(dp) :: cNaCl              ! concentration of salt in bulk in mol/liter
    real(dp) :: cKCl               ! concentration of salt in bulk in mol/liter
    real(dp) :: cRbCl              ! concentration of RbCl in bulk in mol/liter
    real(dp) :: cCaCl2             ! concentration of CaCl2 in bulk in mol/liter
    real(dp) :: cMgCl2             ! concentration of MgCl2 in bulk in mol/liter


    type (looplist), target :: pH
    real(dp) :: pHbulk             ! pH of bulk pH = -log([H+])
    real(dp) :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
  
  

    !  retrun error of subroutine read_pKds
    integer, parameter ::  err_pKdfile_noexist = 1
    integer, parameter ::  err_pKdfile         = 2 
    integer, parameter ::  err_pKderror        = 3

    private :: err_pKdfile_noexist,err_pKdfile,err_pKderror  

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
            case ("brush_mul") 
                neq = (2+nsegtypes) * nsize + neq_bc
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
      
        !     .. charges  
       ! zsurf =-1                 ! valence surface charge 
        zNa   = 1                 ! valence positive charged ion
        zK    = 1                 ! valence positive charged ion
        zRb   = 1                 ! valence positive charged ion
        zCa   = 2                 ! valence divalent positive charged ion
        zMg   = 2                 ! valence divalent positive charged ion
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
        
        zpolAA(1)=-1 ! A-
        zpolAA(2)= 0 ! AH
        zpolAA(3)= 0 ! ANa
        zpolAA(4)= 1 ! ACa+
        zpolAA(5)= 0 ! A2Ca
        zpolAA(6)= 1 ! AMg+
        zpolAA(7)= 0 ! A2Mg

        !  .. radii
        !  .. ionic radii
        !  .. https://www.chemguide.co.uk/atoms/properties/atradius.html and http://abulafia.mt.ic.ac.uk/shannon/ptable.php

        RNa = 0.102_dp             ! radius of Na+ in nm
        RK  = 0.138_dp             ! radius of K+ in nm
        RCl = 0.181_dp             ! radius of Cl- in nm
        RCa = 0.106_dp             ! radius of Ca2+ in nm
        RRb = 0.152_dp             ! radius of Rb+ in nm 
        RMg = 0.072_dp             ! radius of Mg2+ in nm 
        
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
        vRb  = ((4.0_dp/3.0_dp)*pi*(RRb)**3)/vsol 
        vCl  = ((4.0_dp/3.0_dp)*pi*(RCl)**3)/vsol 
        vCa  = ((4.0_dp/3.0_dp)*pi*(RCa)**3)/vsol 
        vMg  = ((4.0_dp/3.0_dp)*pi*(RMg)**3)/vsol 

        vNaCl= (vNa+vCl)          ! contact ion pair
        vKCl = (vK+vCl)           ! contact ion pair
        
        !     .. volume polymer segments
        !     .. all volume scaled by vsol
        
        vAA  =  0.07448_dp/vsol ! volume based on VdW radii 
        vAMPS = 0.2134_dp/vsol
        vPEG  = 0.065_dp/vsol

        !    .. volume AA and vAMPS

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

        ! dissociation constant AA and AMPS

        pKaA(1) =  5.0_dp
        pKaA(2) = -0.4_dp
        pKaA(3) =  1.0_dp
        pKaA(4) =  4.0_dp
        pKaB(1) = -2.0_dp
        pKaB(2) = -0.42_dp
        pKaB(3) = -0.72243_dp
        pKaB(4) = -10.0_dp

        !     .. other physical varaibles
        lsegPAA  = 0.36287_dp       ! segment length in nm
        lsegPAMPS = 0.545_dp        ! segment length in nm
        lsegPEG =  0.3_dp           ! segment length in nm 
              
        lsegA=lsegPAA            
        lsegB=lsegPAMPS 

        ! ! see also subroutine set_chain_properties 

        pKw=14.0_dp                 ! water equilibruim constant
        Tref=298.0_dp               ! temperature in Kelvin
        dielectW=78.54_dp           ! dielectric constant water
        
        seed  = 435672              ! seed for random number generator

        call init_elect_constants(Tref)  

        max_confor=cuantas

    end subroutine init_constants


    function volumesphere(radius)result(volume)

        use mathconst

        real(dp), intent(in) :: radius
        real(dp) :: volume

        volume=(4.0_dp/3.0_dp)*pi*(radius**3)
    
    end function
      

    function radiussphere(volume)result(radius)

        use mathconst

        real(dp), intent(in) :: volume
        real(dp) :: radius

        radius=(volume*3.0_dp/(4.0_dp*pi))**(1.0_dp/3.0_dp)
    
    end function
      

    ! init variables specific for systype=brush, brushborn etc 
    ! variable are  constants, deltaG and K for sytyep=brush,bruhborn etc.
    ! pre : nsegtype, vsol,vpol, vNa etc and ismonomer_chargable need to be set 
    ! post : equlibriuem constant and volume set for charge state of 
    !       carboxylic group of systype ==brushborn are initliazed 

    subroutine init_dna  

        use globals, only : nsegtypes,nseg,systype
        use chains, only : type_of_monomer_char,type_of_monomer,ismonomer_chargeable
        use physconst, only : Na

        real(dp) :: KAA(6)
        real(dp) :: vA    
        integer  :: tAA,i,tt,s,flag_one
        logical  :: isOandNpresent,  isApresent
        integer   :: info

        ! determine segment type number of phosphate constaining segments

        do s=1, nseg
            if(type_of_monomer_char(s)=="P")  tA =type_of_monomer(s)
        enddo

        isApresent=(tA/=0) ! check if phosphate acid monomer is defined in list of typesfname

        if(.not.isApresent) then
            print*,"Error in init_dna:"
            print*,"A momomer is not defined in typesfname"
            print*,"tA= ",tA
            stop
        endif    

        call read_pKds(pKaAA,info)
        if(info/=0) then 
            if(info==err_pKdfile_noexist) then 
                ! set equilbrium constant for acrylic acid 
                pKaAA(1)=5.0_dp
                pKaAA(2)=-0.4_dp
                pKaAA(3)=1.0_dp
                pKaAA(4)=4.0_dp
                pKaAA(5)=1.0_dp
                pKaAA(6)=4.0_dp
            else ! errro
                print*,"Error in init_dna:"
                print*,"Failure to init pKaAA: info=",info
                stop
            endif    


        endif
        
        do i=1,6
            KaAA(i)=10.0_dp**(-pKaAA(i))  
            K0aAA(i) = KaAA(i)*(vsol*Na/1.0e24_dp)
        enddo

        K0aAA(4) = K0aAA(4)*(vsol*Na/1.0e24_dp) ! A2Ca
        K0aAA(6) = K0aAA(6)*(vsol*Na/1.0e24_dp) ! A2Mg 

        ! set volumes 
         
        vA=vpol(tA) 
        vpolAA(1)= vA              ! vA-
        vpolAA(2)= vA              ! vAH
        vpolAA(3)= vA+vNa          ! vANa
        vpolAA(4)= vA+vCa          ! vACa
        vpolAA(5)= 2.0_dp*vA+vCa   ! vA2Ca 
        vpolAA(6)= vA+vMg          ! vAMg       
        vpolAA(7)= 2.0_dp*vA+vMg   ! vA2Mg        

        deltavAA(1)=vpolAA(1)+1.0_dp-vpolAA(2) ! vA- + vH+ - vAH
        deltavAA(2)=vpolAA(1)+vNa-vpolAA(3)    ! vA- + vNa+ - vANa+
        deltavAA(3)=vpolAA(1)+vCa-vpolAA(4)    ! vA- + vCa2+ - vACa+
        deltavAA(4)=2.0_dp*vpolAA(1)+vCa-vpolAA(5) ! 2vA- + vCa2+ -vA2Ca 
        deltavAA(5)=vpolAA(1)+vMg-vpolAA(6)    ! vA- + vMg2+ - vAMg+
        deltavAA(6)=2.0_dp*vpolAA(1)+vMg-vpolAA(7) ! 2vA- + vMg2+ -vA2Mg     


        ! determine if there is only one seg type is chargeable
        flag_one=0
        do tt=1,nsegtypes
            if(ismonomer_chargeable(tt)) flag_one=flag_one+1    
        enddo 

        if(flag_one==0) then
            print*,"Warning: init_dna: zero chargeable acid monomers types"
        else if(flag_one>1) then
            print*,"Warning: init_dna: more then one chargeable acid monomer types"
        endif    

        if(systype=="brushborn") then

            ! make zpolAA zero if monomor t=tA has zpol(1)=zpol(2)=0
            ! see routine read_pKas_and_zpol
            if(.not.ismonomer_chargeable(tA)) zpolAA=0
        
            ! assign born radius charged states of acrylic acid monomer  
            ! see also init_constants
       
            bornrad%pol   = radiussphere(vpolAA(1)*vsol)  ! A^-
            bornrad%polCa = radiussphere(vpolAA(4)*vsol)  ! ACa^+ 
            bornrad%polMg = radiussphere(vpolAA(6)*vsol)  ! AMg^+
            
        endif

    end subroutine init_dna
     


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
        use physconst

        real(dp), intent(in) :: Temp 

        dielectW=dielectric_constant_water(Temp)

        lb=(elemcharge**2)/(4.0_dp*pi*dielectW*dielect0*kBoltzmann*Temp) ! bjerrum length in water=solvent in m
        lb=lb/1.0e-9_dp                           ! bjerrum length in water in nm
        constqW = delta*delta*(4.0_dp*pi*lb)/vsol ! multiplicative constant Poisson Eq. 

    end subroutine init_elect_constants
   
    ! compute surface coverge based on number of grafted point (ngr) 
    subroutine init_sigma()

        use globals, only : systype
        use volume, only : ngr, nsurf, delta

        select case (systype) 
        case("elect","electA","neutral") 
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
        case('brush_mul')    
            sigmaABL = ngr/(nsurf*delta*delta)
            sigmaABR = 0.0_dp
            sigmaAB  = 0.0_dp    
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
        
        real(dp) :: xNaClsalt          ! volume fraction of salt in bulk
        real(dp) :: xKClsalt           ! volume fraction of salt in bulk
        real(dp) :: xCaCl2salt         ! volume fraction of divalent salt in bulk
  

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
           
            call solver(x, xguess, tol_conv, fnorm, issolution) 
            
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
             KaA(i)  = 10.0_dp**(-pKaA(i)) ! experimental equilibruim constant acid 
             KaB(i)  = 10.0_dp**(-pKaB(i)) ! experimental equilibruim constant acid
             K0aA(i) = (KaA(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
             K0aB(i) = (KaB(i)*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
        enddo
        !     .. rescale for i=4 2A- Ca <=> A2Ca
          
        K0aA(4) = (K0aA(4)*vsol)*(Na/1.0e24_dp)
        K0aB(4) = (K0aB(4)*vsol)*(Na/1.0e24_dp)
         
        ! pKa constant from file assogned from read_pKas_and_zpol

        Ka  = 10.0_dp**(-pKa) ! experimental equilibruim constant acid 
        K0a = (Ka*vsol)*(Na/1.0e24_dp) ! intrinstic equilibruim constant 
 

        !pibulk = -log(xbulk%sol)  ! pressure (pi) of bulk
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

    end subroutine init_expmu_neutral

    subroutine init_expmu

        use globals, only : systype
        implicit none

        select case (systype)
        case ( "elect","electA","electdouble","electnopoly")
            call init_expmu_elect()
        case("neutral") 
            call init_expmu_neutral()
        case default 
            print*,"Error in call to init_expmu subroutine"    
            print*,"Wrong value systype : ", systype
            stop        
        end select  

    end subroutine init_expmu


    ! inits chem potential 

    subroutine init_vars_input()

        use globals, only : systype
     
        implicit none
        
        ! local variable
        integer :: i

        select case (systype)
        case ("elect","electA","electdouble","electnopoly")
            call init_expmu_elect()
        case ("neutral")
            call init_expmu_neutral()   
        case("brush_mul") 
            call init_expmu_elect()      
        case("brushssdna") 
            call init_dna  
            !call init_expmu_elect()
        case default   
            print*,"Error: systype incorrect at init_vars_input" 
            print*,"Wrong value systype : ", systype
            print*,"stopping program"
            stop
         end select
             
    end subroutine init_vars_input



   
    subroutine allocate_chain_parameters
        
        use globals, only : nsegtypes
        
        !  allocate array depending on nsegtypes

        allocate(vpol(nsegtypes)) !  volume polymer segments, all volume scaled by vsol
        allocate(pKa(nsegtypes))  !  equilibrium constants 
        allocate(Ka(nsegtypes))   
        allocate(K0a(nsegtypes))  
        allocate(zpol(nsegtypes,2)) !  charge of segment of given type       
        allocate(qpol(nsegtypes)) ! total charge of polymer type 
        allocate(avfdis(nsegtypes)) !  fraction of charge of polymer type 
        allocate(lsegAA(nsegtypes))

    end subroutine allocate_chain_parameters

    !   call to init_chain_parameters 
    !   post/after  allocate_chain_parameters and init_constants    

    subroutine init_chain_parameters
                
        call init_volume_pol 
        call init_pKas_and_zpol
        call init_lseg  ! init segment length

    end subroutine init_chain_parameters


    subroutine init_volume_pol

        use globals, only : nsegtypes

        call read_volume_pol(vpol,vsol,vpolfname, nsegtypes)
        
    end subroutine init_volume_pol
        

    subroutine init_pKas_and_zpol

        use globals, only : nsegtypes

        pKa  =0.0_dp
        zpol =0.0_Dp
    
        call read_pKas_and_zpol(pKa,zpol,pKafname, nsegtypes) 
       
    end subroutine init_pKas_and_zpol


    ! Warning !!this routine is semi redundant because lseg not used in VdW yet !!!

    subroutine init_lseg 

       use globals, only :nsegtypes
       
       call read_lseg(lsegAA,lsegfname, nsegtypes)

    end subroutine init_lseg


    !  .. assign vpol from values in file named filename
    !  .. values vpol are normalized by vsol
    !  .. checked if file exists- content and length not checked     

    subroutine read_lseg(lsegAA,filename, ntypes)

        use  myutils

        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout)  :: lsegAA(:) 
        character(40), intent(in) :: filename  
        integer,  intent(in) :: ntypes 

        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

        !     .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error read_lseg: opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
            stop
        endif
    
        do t=1,ntypes   
            read(un,*)lsegAA(t)
        enddo    

        close(un)
   
    end subroutine read_lseg


    !  .. assign vpolfrom values in file named filename
    !  .. values vpol are normalized by vsol
    !  .. checked if file exists- content and length not checked     

    subroutine read_volume_pol(vpol,vsol,filename, ntypes)

        use  myutils

        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout)  :: vpol(:) 
        real(dp), intent(in) :: vsol 
        character(40), intent(in) :: filename  
        integer,  intent(in) :: ntypes 

        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

    
        !     .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error read_volume_pol: opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
        !      .. local variables
            stop
        endif
    
        do t=1,ntypes   
            read(un,*)vpol(t)
            vpol(t)=vpol(t)/vsol
        enddo    

        close(un)
   
    end subroutine read_volume_pol

    
    !  .. assign  pKa and zpol from values in file named filename
   
    subroutine read_pKas_and_zpol(pKa,zpol,filename, ntypes)
   
        use  myutils
        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout) :: pKa(:)
        integer, intent(inout) ::  zpol(:,:)
        character(lenfname), intent(in) :: filename
        integer,  intent(in) :: ntypes

        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: t 
        character(80) :: istr,str

        !     .. reading in of variables from file
        open(unit=newunit(un),file=filename,iostat=ios,status='old')
        if(ios/=0 ) then
            write(istr,'(I2)')ios
            str='Error opening file '//trim(adjustl(filename))//' : iostat = '//istr
            print*,str
            stop
        endif
        
        do t=1,ntypes   
            read(un,*)pKa(t),zpol(t,1),zpol(t,2)
        enddo    
        
        close(un)


    end subroutine read_pKas_and_zpol

    ! read pKd for acid group tA including acid-base equilbrium Na condensation etc
    ! four return value 
    ! info=0 , correct, 
    ! info= err_pKdfile,err_pKdfile_noexist or err_pKderror : failure

    subroutine read_pKds(pKd,info)
   
        use  myutils
        implicit none 
        
        !     .. arguments 
        real(dp), intent(inout) :: pKd(:)
        integer, intent(out) :: info 
       
        !      .. local variables
        integer :: ios, un  ! un = unit number
        integer :: line, maxline
        character(len=80) :: istr,str
        logical :: exist
        character(len=10) :: fname

        info = 0
        !     .. reading in of variables from file
        write(fname,'(A10)')'pKdacid.in'
        inquire(file=fname,exist=exist)
        if(exist) then
            open(unit=newunit(un),file=fname,status='old',iostat=ios)
            if(ios >0 ) then
                print*, 'Error opening file : iostat =', ios
                info = err_pKdfile
                return
            endif
        else
            info= err_pKdfile_noexist 
            return
        endif
        
        line=0
        ios=0
        !maxline=size(pKd)
        
        do while (line<maxline.and.ios==0)
            line=line+1
            read(un,*,iostat=ios)pKd(line)
        enddo
        
        close(un)

    
        if(line/=maxline.or.ios/=0) then 
            str="reached end of file before all elements read or ios error"
            print*,str
            str="read file "//trim(adjustl(fname))//" failed"
            print*,str
            info = err_pKderror
            return
        endif
        

    end subroutine read_pKds


 end module parameters
