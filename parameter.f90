module parameters

    use physconst
    use mathconst
    use random
    use volume
    use molecules

    implicit none

    !  .. list of parameters

    type(moleclist) :: xbulk,expmu

    !  .. volume 
  
    real*8 :: vsol               ! volume of solvent  in nm^3       
    real*8 :: vpolB(5)           ! volume of one polymer segment, vpol  in units of vsol
    real*8 :: vpolA(5)           ! volume of one polymer segment, vpol  in units of vsol
    real*8 :: vpolC              ! volume of one polymer segment hydrocarbon, vpol  in units of vsol
  
    real*8 :: deltavA(4)
    real*8 :: deltavB(4)
  
    real*8 :: vNa                ! volume positive ion in units of vsol
    real*8 :: vK                 ! volume positive ion in units of vsol
    real*8 :: vCl                ! volume negative ion in units of vsol   
    real*8 :: vCa                ! volume positive divalent ion in units of vsol
    real*8 :: vNaCl
    real*8 :: vKCl
  
    !  .. radii
  
    real*8 :: RNa
    real*8 :: RK
    real*8 :: RCl
    real*8 :: RCa
    real*8 :: RNaCl
    real*8 :: RKCl
    
    real*8 :: lsegAB
    real*8 :: lsegA              ! segment length of A polymer in nm
    real*8 :: lsegB              ! segment length of B polymer in nm
    real*8 :: lsegC              ! segment length of C polymer in nm
  
    integer :: period            ! peridociy of repeat of A or B block 
  
    real*8 :: VdWepsB            ! strenght VdW interaction in units of kT
    real*8 :: VdWepsC            ! strenght VdW interaction in units of kT
    real*8 :: chibulk            ! value of chibulk 
    integer :: numlayers
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
  
    real*8 :: T                  ! temperature in K
    real*8 :: dielectW           ! dielectric constant of water 
    real*8 :: lb                 ! Bjerrum length	   
    real*8 :: constqW            ! constant in Poisson eq dielectric constant of water 

    real*8 :: sigmaAB            ! sigma AB polymer coated on planar surface
    real*8 :: sigmaABL           ! sigma AB polymer coated on planar surface
    real*8 :: sigmaABR           ! sigma AB polymer coated on planar surface
    
    real*8 :: sigmaC             ! sigma C polymer coated on planar surface
  
    integer :: itmax             ! maximum number of iterations
    real*8 :: error              ! error imposed accuaracy
    real*8 :: fnorm              ! L2 norm of residual vector function fcn  
    integer :: infile            ! infile=1 read input files infile!=1 no input files 
    integer :: iter              ! counts number of iterations
  
    character(len=8) :: method           ! method="kinsol" or "zspow"  
    character(len=8) :: chainmethod      ! method of generating chains ="MC" or "FILE" 
    character(len=8) :: chaintype        ! type of chain: diblock,alt
    integer :: readinchains              ! nunmber of used/readin chains
    integer, parameter :: numsys= 6      ! number of method   
    character(len=15) :: sysvalues(numsys) ! different system methodes 
    character(len=2) :: bcvalues(2,5)      ! boundary condition bc="qu" quartz,"cl" clay or "ca" calcite
    character(len=3) ::  verboseflag       ! select input falg 

    real*8 :: heightAB           ! average height of layer
    real*8 :: heightC            ! average height of layer 
    real*8 :: qpolA              ! charge poly A of layer 
    real*8 :: qpolB              ! charge poly B of layer 
    real*8 :: qpol_tot           ! charge poly A+B of layer 
    real*8 :: avfdisA(5)         ! average degree of dissociation 
    real*8 :: avfdisB(5)         ! average degree of dissociation
  
  !     .. weak polyelectrolyte variables 
  !     .. equibrium constant
  
    real*8 :: K0A(4)              ! intrinsic equilibruim constant
    real*8 :: KA(4)               ! experimemtal equilibruim constant 
    real*8 :: pKA(4)              ! experimental equilibruim constant pKa= -log[Ka]
    real*8 :: K0B(4)              ! intrinsic equilibruim constant
    real*8 :: KB(4)               ! experimemtal equilibruim constant 
    real*8 :: pKB(4)              ! experimental equilibruim constant pKa= -log[Ka]
  
    real*8 :: pKw                 ! water equilibruim constant pKw= -log[Kw] ,Kw=[H+][OH-] 
  
    real*8 :: K0ionNa             ! intrinsic equilibruim constant
    real*8 :: KionNa              ! experimemtal equilibruim constant 
    real*8 :: pKionNa             ! experimental equilibruim constant pKion= -log[Kion]	 
  
    real*8 :: K0ionK              ! intrinsic equilibruim constant
    real*8 :: KionK               ! experimemtal equilibruim constant 
    real*8 :: pKionK              ! experimental equilibruim constant pKion= -log[Kion]	 
  
    !     .. bulk volume fractions 
  
    real*8 :: pibulk             ! -ln(xsolbulk)
    real*8 :: xNaClsalt          ! volume fraction of salt in bulk
    real*8 :: xKClsalt           ! volume fraction of salt in bulk
    real*8 :: xCaCl2salt         ! volume fraction of divalent salt in bulk
  
    real*8 :: cHplus             ! concentration of H+ in bulk in mol/liter
    real*8 :: cOHmin             ! concentration of OH- in bulk in mol/liter
    real*8 :: cNaCl              ! concentration of salt in bulk in mol/liter
    real*8 :: cKCl               ! concentration of salt in bulk in mol/liter
    real*8 :: cCaCl2             ! concentration of divalent salt in bulk in mol/liter
    real*8 :: pHbulk             ! pH of bulk pH = -log([H+])
    real*8 :: pOHbulk            ! p0H of bulk p0H = -log([0H-])
  
        
    contains

        ! determine total number of non linear equations

        subroutine set_size_neq()

            use globals
            use volume, only : nz

            implicit none

            integer*8 :: neq_bc

            neq_bc=0 
            if(bcflag(LEFT)/="cc") neq_bc=neq_bc+1
            if(bcflag(RIGHT)/="cc") neq_bc=neq_bc+1

            select case (sysflag)
                case ("elect") 
                    neq = 4 * nz + neq_bc
                case ("electdouble")  
                    neq = 4 * nz 
                case ("electnopoly") 
                    neq = 2 * nz + neq_bc
                case ("electHC") 
                    neq = 5 * nz +neq_bc
                case ("neutral") 
                    neq = 2 * nz
                case ("bulk water") 
                    neq = 5 
                case default
                    print*,"Wrong value sysflag:  ",sysflag
                    stop
            end select  
            
            print*,"neq = ",neq        
        
        end subroutine set_size_neq

        ! defines number of different systems

        subroutine init_allowed_flags()

            use globals, only : LEFT,RIGHT
            implicit none

            sysvalues(1)="elect"
            sysvalues(2)="electdouble"
            sysvalues(3)="electnopoly"
            sysvalues(4)="electHC" 
            sysvalues(5)="neutral"
            sysvalues(6)="bulk water"

            bcvalues(RIGHT,1)="qu"
            bcvalues(RIGHT,2)="cl"
            bcvalues(RIGHT,3)="ca"
            bcvalues(RIGHT,4)="ta"
            bcvalues(RIGHT,5)="cc"

            bcvalues(LEFT,1)="ta"
            bcvalues(LEFT,2)="cc"


        end subroutine init_allowed_flags


        subroutine check_value_sysflag(fcnname)

            use globals
            
            implicit none

            character(len=80), optional, intent(in) :: fcnname

            integer :: i
            logical :: flag


            flag=.FALSE.

            do i=1,numsys
                if(sysflag==sysvalues(i)) flag=.TRUE.
            enddo
            if (flag.eqv. .FALSE.) then 
                print*,"Error value of sysflag is not permissible"
                print*,"sysflag = ",sysflag
                if(present(fcnname)) print*,"Error in ",fcnname
                stop
            endif
    
        end subroutine 

        subroutine check_value_bcflag(fcnname)
            use globals

            implicit none

            integer :: i
            logical :: flag

            character(len=80), optional, intent(in) :: fcnname

            flag=.FALSE.

            do i=1,5
                if(bcflag(RIGHT)==bcvalues(RIGHT,i)) flag=.TRUE.
            enddo
            if (flag.eqv. .FALSE.) then
                print*,"Error value of bcflag is not permissible"
                print*,"bcflag(RIGHT) = ",bcflag(RIGHT)
                if(present(fcnname)) print*,"Error in ",fcnname
                stop
            endif

            flag=.FALSE.
            do i=1,2
                if(bcflag(LEFT)==bcvalues(LEFT,i)) flag=.TRUE.
            enddo
            if (flag.eqv. .FALSE.) then
                print*,"Error value of bcflag is not permissible"
                print*,"bcflag(LEFT) = ",bcflag(LEFT)
                if(present(fcnname)) print*,"Error in ",fcnname
                stop
            endif

        end subroutine  

        real*8 function BjerrumLenght(T)

            use mathconst
            use physconst
            
            implicit none
            
            real*8, intent(in) :: T     
            real*8 :: lb

            lb=(elemcharge**2)/(4.0d0*pi*dielectW*dielect0*kBoltzmann*T) ! bjerrum length in water=solvent in m
            lb=lb/1.0d-9              ! bjerrum length in water in nm
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
            
            integer :: i
            real*8 :: vA,vB
            
            !  .. initializations of variables
     
            pi=dacos(-1.0d0)          ! pi = arccos(-1)
            itmax=2000                ! maximum number of iterations
            nz=nsize                  ! size of lattice in z-direction 
            
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
            
            RNa = 0.102d0             ! radius of Na+ in nm
            RK  = 0.138d0             ! radius of K+ in nm
            RCl = 0.181d0             ! radius of Cl- in nm
            RCa = 0.106d0             ! radius of Ca2+ in nm
            RNaCl = 0.26d0            ! radius of ion pair: this value is strange  and not used !!
            RKCl = 0.26d0             ! radius of ion pair
            
            !     .. volume
            if(sysflag/="neutral") then 
                vsol = 0.030d0              ! volume water solvent molecule in (nm)^3
            elseif(sysflag=="neutral") then 
                vsol = 0.218d0              ! volume hexane Mw=86.18 g/mol and rho=0.6548 g/ml  
            else 
                print*,"Error in call to init_constants subroutine"   
                print*,"Wrong system flag"
                stop
            endif   

            vNa  = ((4.0d0/3.0d0)*pi*(RNa)**3)/vsol 
            vK   = ((4.0d0/3.0d0)*pi*(RK)**3)/vsol 
            vCl  = ((4.0d0/3.0d0)*pi*(RCl)**3)/vsol 
            vCa  = ((4.0d0/3.0d0)*pi*(RCa)**3)/vsol 
            vNaCl= (vNa+vCl)          ! contact ion pair
            vKCl = (vK+vCl)           ! contact ion pair
            
            !     .. volume polymer segments
            !     .. all volume scaled by vsol
            
            vA   = 0.07448d0/vsol ! volume based on VdW radii 
            vB   = 0.2134d0/vsol
            
            vpolA(1)= vA              ! vA-
            vpolA(2)= vA              ! vAH
            vpolA(3)= vA+vNa          ! vANa
            vpolA(4)= vA+vCa          ! vACa
            vpolA(5)= 2.0d0*vA+vCa    ! vA2Ca
            
            vpolB(1)= vB              ! vB-
            vpolB(2)= vB              ! vBH
            vpolB(3)= vB+vNa          ! vBNa
            vpolB(4)= vB+vCa          ! vBCa
            vpolB(5)= 2.0d0*vB+vCa    ! vB2Ca
            
            deltavA(1)=vpolA(1)+1.0d0-vpolA(2) ! vA-+vH+-vAH
            deltavA(2)=vpolA(1)+vNa-vpolA(3) ! vA-+vNa+-vANa+
            deltavA(3)=vpolA(1)+vCa-vpolA(4) ! vA- + vCa2+ -vACa+
            deltavA(4)=2.0d0*vpolA(1)+vCa-vpolA(5) ! 2vA- + vCa2+ -vA2Ca
            
            deltavB(1)=vpolB(1)+1.0d0-vpolB(2) ! vB-+vH+-vBH
            deltavB(2)=vpolB(1)+vNa-vpolB(3) ! vB-+vNa+-vBNa+
            deltavB(3)=vpolB(1)+vCa-vpolB(4) ! vB- +vCa2+ -vBCa+
            deltavB(4)=2.0d0*vpolB(1)+vCa-vpolB(5) ! 2vB- + vCa2+ -vB2Ca+
            
            vpolC  = 0.0270d0/vsol     ! volume CH2
           
            !     .. other physical varaibles
            
            lsegAB=0.545d0             ! segment length in nm
            lsegA=0.545d0              ! segment length in nm
            lsegB=0.545d0              ! segment length in nm 
            lsegC=0.153d0              ! segment length in nm od CH2 check 
            
            pKw=14.0d0                 ! water equilibruim constant
            
            T=298.0d0                  ! temperature in Kelvin
            dielectW=78.54d0           ! dielectric constant water
            
            lb=BjerrumLenght(T)        ! bjerrum length in water in nm

            delta = 0.30d0             ! size lattice spacing
            seed  = 435672             ! seed for random number generator
            
            constqW = delta*delta*4.0d0*pi*lb/vsol ! multiplicative constant Poisson Eq. 
            
            !  .. initializations of input dependent variables 
            
            sigmaABL = sigmaABL * (1.0d0/(delta)) ! dimensionless sigma no vpol*vsol !!!!!!!!!!!! 
            sigmaABR = sigmaABR * (1.0d0/(delta)) 
            sigmaAB = sigmaAB * (1.0d0/(delta))  
            sigmaC   = sigmaC  * (1.0d0/(delta)) ! dimensionless sigma no vpol*vsol !!!!!!!!!!!!
            
            ! VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
            ! VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 
            
            ! .. make radius integer multiply of delta
            ! .. needed because VdW-coefficeint computed on grid 
            ! Íradius=delta*int(radius/delta)

            max_conforAB=cuantasAB
            max_conforC=cuantasC



        end subroutine init_constants
   
   
        !     purpose: initialize expmu needed by fcn 
        !     pre: first read_inputfile has to be called

        subroutine init_expmu_elect()
     
            use globals
            use physconst
            
            implicit none 
            
            !     .. local variable
            
            real*8,  dimension(:), allocatable :: x         ! volume fraction solvent iteration vector 
            real*8,  dimension(:), allocatable :: xguess  
            
            integer :: i
            character(len=15) :: sysflag_old
            
            allocate(x(5))
            allocate(xguess(5))
            
            !     .. initializations of input dependent variables, electrostatic part 
            
            cHplus = (10d0)**(-pHbulk) ! concentration H+ in bulk
            pOHbulk = pKw -pHbulk       
            cOHmin  = (10d0)**(-pOHbulk) ! concentration OH- in bulk
            
            xbulk%Hplus = (cHplus*Na/(1.0d24))*(vsol) ! volume fraction H+ in bulk vH+=vsol
            xbulk%OHmin = (cOHmin*Na/(1.0d24))*(vsol) ! volume fraction OH- in bulk vOH-=vsol
            
            xNaClsalt = (cNaCl*Na/(1.0d24))*((vNa+vCl)*vsol) ! volume fraction NaCl salt in mol/l
            
            if(pHbulk.le.7) then      ! pH<= 7
                xbulk%Na=xNaClsalt*vNa/(vNa+vCl)  
                xbulk%Cl=xNaClsalt*vCl/(vNa+vCl) +(xbulk%Hplus -xbulk%OHmin)*vCl  ! NaCl+ HCl
            else                      ! pH >7
                xbulk%Na=xNaClsalt*vNa/(vNa+vCl) +(xbulk%OHmin -xbulk%Hplus)*vNa ! NaCl+ NaOH  
                xbulk%Cl=xNaClsalt*vCl/(vNa+vCl)  
            endif
            
            xKClsalt = (cKCl*Na/(1.0d24))*((vK+vCl)*vsol) ! volume fraction KCl salt in mol/l
            xbulk%K = xKClsalt*vK/(vK+vCl)  
            xbulk%Cl = xbulk%Cl+xKClsalt*vCl/(vK+vCl)  
            
            xCaCl2salt = (cCaCl2*Na/(1.0d24))*((vCa+2.0d0*vCl)*vsol) ! volume fraction CaCl2 in mol/l
            xbulk%Ca=xCaCl2salt*vCa/(vCa+2.0d0*vCl)
            xbulk%Cl=xbulk%Cl+ xCaCl2salt*2.0d0*vCl/(vCa+2.0d0*vCl)
            
            xbulk%NaCl=0.0d0    ! no ion pairing
            xbulk%KCl=0.0d0     ! no ion pairing
            
            xbulk%sol=1.0d0 -xbulk%Hplus -xbulk%OHmin -xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca   
            
            !     .. if Kion neq 0 ion pairing !
            
            !     .. intrinstic equilibruim constant acid        
            !     Kion  = 0.246d0 ! unit 1/M= liter per mol !!!
            K0ionK  = KionK /(vsol*Na/1.0d24) ! intrinstic equilibruim constant 
            K0ionNa = KionNa/(vsol*Na/1.0d24) ! intrinstic equilibruim constant 
            
            if((KionNa.ne.0.0d0).or.(KionK.ne.0.0d0)) then  
                sysflag_old=sysflag 
                sysflag="bulk water"        ! set solver to fcnbulk
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
               
                call solver(x, xguess, error, fnorm) 
                
                !     .. return solution
                
                xbulk%Na  =x(1)
                xbulk%Cl  =x(2)
                xbulk%NaCl=x(3)
                xbulk%K   =x(4)
                xbulk%KCl =x(5)

                ! reset of flags
                iter=0
                sysflag=sysflag_old         ! switch solver back
                call set_size_neq()         ! number of non-linear  equation        
                
                xbulk%sol=1.0d0-xbulk%Hplus-xbulk%OHmin - xbulk%Cl -xbulk%Na -xbulk%K-xbulk%NaCl-xbulk%KCl-xbulk%Ca 
                
            endif
             
            !     .. intrinstic equilibruim constants      
            do i=1,4
                 Ka(i)  = 10d0**(-pKa(i)) ! experimental equilibruim constant acid 
                 Kb(i)  = 10d0**(-pKb(i)) ! experimental equilibruim constant acid
                 K0a(i) = (Ka(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
                 K0b(i) = (Kb(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
            enddo
            !     .. rescale for i=4 2A- Ca <=> A2Ca
              
            K0a(4) = (Ka(4)*vsol)*(Na/1.0d24)
            K0b(4) = (Kb(4)*vsol)*(Na/1.0d24)
             
              
            pibulk = -dlog(xbulk%sol)  ! pressure (pi) of bulk
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
                
            VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
            VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled 

            deallocate(x)
            deallocate(xguess)
            
        end subroutine init_expmu_elect


        subroutine init_expmu_neutral

            implicit none
              
            xbulk%sol=1.0d0 ! only solvent 
    
            VdWepsC  = VdWepsC/(vpolC*vsol) ! VdW eps scaled 
            VdWepsB  = VdWepsB/(vpolB(3)*vsol) ! VdW eps scaled

        end subroutine init_expmu_neutral

        subroutine init_expmu

            use globals, only : sysflag
            implicit none


            if(sysflag=="elect") then 
                call init_expmu_elect()
            elseif(sysflag=="electdouble") then 
                call init_expmu_elect()
            elseif(sysflag=="electnopoly") then 
                call init_expmu_elect()
            elseif(sysflag=="neutral") then
                call init_expmu_neutral()
            else
                print*,"Error in call to init_expmu subroutine"    
                print*,"Wrong value sysflag : ", sysflag
                stop        
            endif   

        end subroutine init_expmu
   
 end module parameters
