
module myio

    use precision_definition
    implicit none

    ! return error values

    integer, parameter ::  myio_err_systype   = 1
 !   integer, parameter ::  myio_err_runtype   = 2
 !   integer, parameter ::  myio_err_geometry  = 3
    integer, parameter ::  myio_err_method    = 4
    integer, parameter ::  myio_err_chaintype = 5
    integer, parameter ::  myio_err_domain    = 6
    integer, parameter ::  myio_err_inputfile = 7
    integer, parameter ::  myio_err_input     = 8
    integer, parameter ::  myio_err_bcflag    = 9 

contains


subroutine read_inputfile(info)

    use globals
    use parameters
    use surface 
    use myutils, only : newunit

    implicit none

    integer, intent(out), optional :: info

    ! .. local arguments

    integer :: info_sys, info_bc, info_run, info_geo, info_meth, info_chaintype, info_combi
    character(len=8) :: fname
    integer :: ios,un_input  ! un = unit number    
    character(len=80) :: fcnname

    if (present(info)) info = 0
    
    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    read(un_input,*)method
    read(un_input,*)systype
    read(un_input,*)bcflag(LEFT)
    read(un_input,*)bcflag(RIGHT)
    read(un_input,*)chainmethod
    read(un_input,*)chaintype
    read(un_input,*)sigmaABL
    read(un_input,*)sigmaABR
    read(un_input,*)sigmaC
    read(un_input,*)error             
    read(un_input,*)infile              ! guess  1==yes
    read(un_input,*)pHbulk
    read(un_input,*)KionNa
    read(un_input,*)KionK
    read(un_input,*)sigmaSurfL
    read(un_input,*)sigmaSurfR
    read(un_input,*)cNaCl
    read(un_input,*)cKCl
    read(un_input,*)cCaCl2
    read(un_input,*)pKa(1)           !   AH   <=> A- + H+ 
    read(un_input,*)pKa(2)           !   ANa  <=> A- + Na+  
    read(un_input,*)pKa(3)           !   ACa+ <=> A- + Ca2+ 
    read(un_input,*)pKa(4)           !   A2Ca <=> 2A- + Ca2+
    read(un_input,*)pKb(1)           !   BH   <=> B- + H+ 
    read(un_input,*)pKb(2)           !   BNa  <=> B- + Na+ 
    read(un_input,*)pKb(3)           !   BCa+ <=> B- + Ca2+   
    read(un_input,*)pKb(4)           !   B2Ca <=> 2B- + Ca2+   
    read(un_input,*)period
    read(un_input,*)nsize
    read(un_input,*)nsegAB
    read(un_input,*)cuantasAB
    read(un_input,*)nsegC
    read(un_input,*)cuantasC
    read(un_input,*)VdWepsC
    read(un_input,*)VdWepsB    
    read(un_input,*)VdWcutoff
    read(un_input,*)nzmax            ! max distance
    read(un_input,*)nzmin            ! min distance
    read(un_input,*)nzstep           ! step distance  
    read(un_input,*)verboseflag     

    close(un_input)


    write(fcnname,'(A14)')'read_inputfile'
    
    ! override input bcflags 
    
    if(systype=="electdouble") then
        bcflag(LEFT)="cc"
        bcflag(RIGHT)="cc"  
    endif
    if(systype=="elect") then
        sigmaAB=sigmaABL
        sigmaABR=0.0_dp  
    endif
    
    
    ! .. check error flag

    call check_value_systype(systype,info_sys) 
    if (info_sys == myio_err_systype) then
        if (present(info)) info = info_sys
        return
    endif

    call check_value_bcflag(bcflag,info_bc) 
    if (info_bc == myio_err_bcflag) then
        if (present(info)) info = info_bc
        return
    endif

    call check_value_method(method,info_meth)
    if (info_meth == myio_err_method) then
        if (present(info)) info = info_meth
        return
    endif

    call check_value_chaintype(chaintype,info_chaintype)
    if (info_chaintype == myio_err_chaintype) then
        if (present(info)) info = info_chaintype
        return
    endif

end subroutine read_inputfile
 

subroutine check_value_systype(systype,info)

    implicit none

    character(len=15), intent(in) :: systype
    integer, intent(out),optional :: info

    character(len=15) :: systypestr(6)
    integer :: i
    logical :: flag

    ! permissible values of systype

    systypestr(1)="elect"
    systypestr(2)="bulk water"
    systypestr(3)="neutral"
    systypestr(4)="electdouble"
    systypestr(5)="electnopoly"
    systypestr(6)="electHC"

    flag=.FALSE.

    do i=1,6
        if(systype==systypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of systype is not permissible"
        print*,"systype = ",systype
        if (present(info)) info = myio_err_systype
        return
    end if

end subroutine check_value_systype


subroutine check_value_bcflag(bcflag,info)

    use globals, only : LEFT, RIGHT

    implicit none

    character(len=2), intent(in), dimension(2) :: bcflag
    integer, intent(out), optional :: info

    character(len=15) :: bcvalues(2,5)
    integer :: i
    logical :: flag

    ! permissible values of bcflag

    bcvalues(RIGHT,1)="qu"
    bcvalues(RIGHT,2)="cl"
    bcvalues(RIGHT,3)="ca"
    bcvalues(RIGHT,4)="ta"
    bcvalues(RIGHT,5)="cc"
    bcvalues(LEFT,1)="ta"
    bcvalues(LEFT,2)="cc"

    flag=.FALSE.
            
    do i=1,5
        if(bcflag(RIGHT)==bcvalues(RIGHT,i)) flag=.TRUE.
    enddo
    if (flag.eqv. .FALSE.) then
        print*,"Error value of bcflag is not permissible"
        print*,"bcflag(RIGHT) = ",bcflag(RIGHT)
        if (present(info)) info = myio_err_systype
        !if (present(fcnname)) print*,"Error in ",fcnname
        stop
    endif

    flag=.FALSE.
    do i=1,2
        if(bcflag(LEFT)==bcvalues(LEFT,i)) flag=.TRUE.
    enddo
    if (flag.eqv. .FALSE.) then
        print*,"Error value of bcflag is not permissible"
        print*,"bcflag(LEFT) = ",bcflag(LEFT)
        if (present(info)) info = myio_err_bcflag
        !if(present(fcnname)) print*,"Error in ",fcnname
        stop
    endif


end subroutine check_value_bcflag


subroutine check_value_chaintype(chaintype,info)

        implicit none

        character(len=8), intent(in) :: chaintype
        integer, intent(out),optional :: info

        logical :: flag
        character(len=8) :: chaintypestr(3)
        integer :: i

        ! permissible values of chaintype

        chaintypestr(1)="diblock"
        chaintypestr(2)="altA"
        chaintypestr(3)="altB"

        flag=.FALSE.

        do i=1,3
            if(chaintype==chaintypestr(i)) flag=.TRUE.
        enddo

        if (present(info)) info = 0

        if (flag.eqv. .FALSE.) then
            print*,"Error: value of chaintype is not permissible"
            print*,"chaintype = ",chaintype
            if (present(info)) info = myio_err_chaintype
            return
        endif

    end subroutine check_value_chaintype

    subroutine check_value_method(method,info)


        implicit none

        character(len=8), intent(in) :: method
        integer, intent(out),optional :: info

        character(len=8) :: methodstr
        integer :: i
        logical :: flag

        ! permissible values of runtype

        methodstr="kinsol"

        flag=.FALSE.

        if (method==methodstr) flag=.TRUE.

        if (present(info)) info = 0

        if (flag.eqv. .FALSE.) then
            print*,"Error: value of method is not permissible"
            print*,"method = ",method
            if (present(info)) info = myio_err_method
            return
        endif

    end subroutine check_value_method




    subroutine output_elect(countfile)
  
       !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    use myutils, only : newunit
  
    implicit none
      
    !     .. scalar arguments
    
    integer :: countfile
    
    !     .. local arguments

    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolCfilename 
    character(len=90) :: xpolendfilename 
    character(len=90) :: xNafilename
    character(len=90) :: xKfilename
    character(len=90) :: xCafilename
    character(len=90) :: xNaClfilename
    character(len=90) :: xKClfilename
    character(len=90) :: xClfilename
    character(len=90) :: potentialfilename
    character(len=90) :: chargefilename
    character(len=90) :: xHplusfilename
    character(len=90) :: xOHminfilename
    character(len=90) :: densfracAfilename
    character(len=90) :: densfracBfilename
    character(len=90) :: densfracionpairfilename

    integer :: un_sys,un_xpol,un_xpolC,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair  ! unit numbers

    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals,fmt   

    integer :: i,j,k          ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    ! .. executable statements 

    ! ..format specifiers 

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 

    !     .. make label filenames 

    write(rstr,'(F5.3)')sigmaAB*delta 
    fnamelabel="sg"//trim(adjustl(rstr)) 
    write(rstr,'(F5.3)')cNaCl
    fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
    write(rstr,'(F5.3)')cCaCl2
    fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
    write(rstr,'(F7.3)')pHbulk
    fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
    write(rstr,'(I5.5)')countfile
    fnamelabel=trim(fnamelabel)//"."//trim(adjustl(rstr))//".dat"


    sysfilename='system.'//trim(fnamelabel)
    xpolABfilename='xpolAB.'//trim(fnamelabel)
    xpolCfilename='xpolC.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    densfracAfilename='densityAfrac.'//trim(fnamelabel)
    densfracBfilename='densityBfrac.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)

    !     .. opening files        
    open(unit=newunit(un_sys),file=sysfilename)       
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)

    if(systype/="electnopoly") then          
        open(unit=newunit(un_xpol),file=xpolABfilename)
        open(unit=newunit(un_xpolC),file=xpolCfilename)
        open(unit=newunit(un_fdisA),file=densfracAfilename) 
        open(unit=newunit(un_fdisB),file=densfracBfilename) 
    endif       

    if(verboseflag=="yes") then    
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif

    ! .. writting files

    write(un_psi,*)0.0_dp,psiSurfL
    do i=1,nz
        write(un_xsol,*)zc(i),xsol(i)
        write(un_psi,*)zc(i),psi(i)
    enddo    
    write(un_psi,*)nz*delta,psiSurfR 

    if(systype/="electnopoly") then 
        do i=1,nz
            write(un_xpol,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
            write(un_xpolC,fmt2reals)zc(i),xpolC(i)
            write(un_fdisA,fmt6reals)zc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
            write(un_fdisB,fmt6reals)zc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)
        enddo
    endif   
    
    if(verboseflag=="yes") then 
        do i=1,nz
            write(un_xNa,*)zc(i),xNa(i)
            write(un_xK,*)zc(i),xK(i)
            write(un_xCa,*)zc(i),xCa(i)
            write(un_xNaCl,*)zc(i),xNaCl(i)
            write(un_xKCl,*)zc(i),xKCl(i)
            write(un_xpair,*)zc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)zc(i),xCl(i)
            write(un_charge,*)zc(i),rhoq(i)
            write(un_xHplus,*)zc(i),xHplus(i)
            write(un_xOHmin,*)zc(i),xOHmin(i)    
        enddo    
    endif

    ! .. writing system information 
    
    write(un_sys,*)'system      = planar weakpolyelectrolyte brush'
    write(un_sys,*)'version     = ',VERSION
    write(un_sys,*)'chainmethod = ',chainmethod
    write(un_sys,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
       write(un_sys,*)'readinchains = ',readinchains
    endif
    write(un_sys,*)'systype     = ',systype
    write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    write(un_sys,*)'sigmaAB     = ',sigmaAB*delta
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'sumphiC     = ',sumphiC
    write(un_sys,*)'check phi   = ',checkphi 
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(un_sys,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEVdW       = ',FEVdW 
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'qAB         = ',qAB
    write(un_sys,*)'qC          = ',qC
    write(un_sys,*)'muAB        = ',-log(qAB)
    write(un_sys,*)'muC         = ',-log(qC)
    write(un_sys,*)'nsegAB      = ',nsegAB
    write(un_sys,*)'lsegAB      = ',lsegAB
    write(un_sys,*)'nsegC       = ',nsegC
    write(un_sys,*)'lsegC       = ',lsegC
    write(un_sys,*)'period      = ',period
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'delta       = ',delta 
    write(un_sys,*)'distance    = ',nz*delta    
    write(un_sys,*)'vsol        = ',vsol
    write(un_sys,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_sys,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_sys,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_sys,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_sys,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_sys,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_sys,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_sys,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_sys,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_sys,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_sys,*)'vpolC       = ',vpolC*vsol
    write(un_sys,*)'vNa         = ',vNa*vsol
    write(un_sys,*)'vCl         = ',vCl*vsol
    write(un_sys,*)'vCa         = ',vCa*vsol
    write(un_sys,*)'vK          = ',vK*vsol
    write(un_sys,*)'vNaCl       = ',vNaCl*vsol
    write(un_sys,*)'vKCl        = ',vKCl*vsol
    write(un_sys,*)'cNaCl       = ',cNaCl
    write(un_sys,*)'cKCl        = ',cKCl
    write(un_sys,*)'cCaCl2      = ',cCaCl2
    write(un_sys,*)'pHbulk      = ',pHbulk
    write(un_sys,*)'pKa         = ',pKa(1)      
    write(un_sys,*)'pKaNa       = ',pKa(2)
    write(un_sys,*)'pKaACa      = ',pKa(3)
    write(un_sys,*)'pKaA2Ca     = ',pKa(4)
    write(un_sys,*)'pKb         = ',pKb(1)      
    write(un_sys,*)'pKbNa       = ',pKb(2)
    write(un_sys,*)'pKbBCa      = ',pKb(3)
    write(un_sys,*)'pKbB2Ca     = ',pKb(4)
    write(un_sys,*)'KionNa      = ',KionNa
    write(un_sys,*)'KionK       = ',KionK
    write(un_sys,*)'K0ionNa     = ',K0ionNa
    write(un_sys,*)'K0ionK      = ',K0ionK
    write(un_sys,*)'xNabulk     = ',xbulk%Na
    write(un_sys,*)'xClbulk     = ',xbulk%Cl
    write(un_sys,*)'xKbulk      = ',xbulk%K
    write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
    write(un_sys,*)'xKClbulk    = ',xbulk%KCl
    write(un_sys,*)'xCabulk     = ',xbulk%Ca
    write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
    write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
    write(un_sys,*)'sigmaAB     = ',sigmaAB*delta
    write(un_sys,*)'sigmaC      = ',sigmaC*delta
    write(un_sys,*)'heightAB    = ',heightAB
    write(un_sys,*)'heightC     = ',heightC
    write(un_sys,*)'qpolA       = ',qpolA
    write(un_sys,*)'qpolB       = ',qpolB
    write(un_sys,*)'qpoltot     = ',qpol_tot
    write(un_sys,*)'avfdisA(1)  = ',avfdisA(1)
    write(un_sys,*)'avfdisA(2)  = ',avfdisA(2)
    write(un_sys,*)'avfdisA(3)  = ',avfdisA(3)
    write(un_sys,*)'avfdisA(4)  = ',avfdisA(4)
    write(un_sys,*)'avfdisA(5)  = ',avfdisA(5)
    write(un_sys,*)'avfdisB(1)  = ',avfdisB(1)
    write(un_sys,*)'avfdisB(2)  = ',avfdisB(2)
    write(un_sys,*)'avfdisB(3)  = ',avfdisB(3)
    write(un_sys,*)'avfdisB(4)  = ',avfdisB(4)
    write(un_sys,*)'avfdisB(5)  = ',avfdisB(5)
    write(un_sys,*)'dielectW    = ',dielectW
    write(un_sys,*)'lb          = ',lb
    write(un_sys,*)'T           = ',T
    write(un_sys,*)'VdWepsC     = ',VdWepsC*vpolC*vsol 
    write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(un_sys,*)'zpolA(1)    = ',zpolA(1)
    write(un_sys,*)'zpolA(2)    = ',zpolA(2)
    write(un_sys,*)'zpolA(3)    = ',zpolA(3)
    write(un_sys,*)'zpolA(4)    = ',zpolA(4)
    write(un_sys,*)'zpolB(1)    = ',zpolB(1)
    write(un_sys,*)'zpolB(2)    = ',zpolB(2)
    write(un_sys,*)'zpolB(3)    = ',zpolB(3)
    write(un_sys,*)'zpolB(4)    = ',zpolB(4)
    write(un_sys,*)'zpolB(5)    = ',zpolB(5)
    write(un_sys,*)'zNa         = ',zNa
    write(un_sys,*)'zCa         = ',zCa
    write(un_sys,*)'zK          = ',zK
    write(un_sys,*)'zCl         = ',zCl
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    write(un_sys,*)'iterations  = ',iter

    write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(un_sys,*)'sigmaSurfL  = ',sigmaSurfL/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaSurfR  = ',sigmaSurfR/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaqSurfL = ',sigmaqSurfL/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaqSurfR = ',sigmaqSurfR/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'psiSurfL    = ',psiSurfL
    write(un_sys,*)'psiSurfR    = ',psiSurfR

    fmt = "(A9,I1,A5,ES25.16)"
    if(bcflag(LEFT)=='ta') then
      do i=1,4   
        write(10,fmt)'fdisTaL(',i,')  = ',fdisTaL(i)
      enddo  
    endif
    if(bcflag(RIGHT)=='ta') then   
      do i=1,4   
        write(10,fmt)' fdisTaR(',i,')  = ',fdisTaR(i)
      enddo  
    endif
    if((bcflag(RIGHT)/='ta').and.(bcflag(RIGHT)/='cc') ) then
      do i=1,6   
        write(10,fmt)' fdisSuR(',i,')  = ',fdisS(i)
      enddo  
    endif

    ! .. closing files

    close(un_sys)
    close(un_xsol)
    close(un_psi)

    if(systype/="electnopoly") then
        close(un_xpol)   
        close(un_xpolC)
        close(un_fdisA)
        close(un_fdisB)
    endif

    if(verboseflag=="yes") then 
        close(un_xNa)   
        close(un_xK)
        close(un_xCa)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
    endif

end subroutine output_elect

 
subroutine output_electdouble(countfile)
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface
    use myutils, only : newunit
    !     use endpoint
  
    implicit none
    
    !     .. scalar arguments

    integer :: countfile
    
    !     .. local arguments

    !     .. output file names       

    character(len=16) :: sysfilename     
    character(len=24) :: xsolfilename 
    character(len=26) :: xpolABfilename
    character(len=26) :: xpolendfilename 
    character(len=24) :: xNafilename
    character(len=24) :: xKfilename
    character(len=26) :: xCafilename
    character(len=27) :: xNaClfilename
    character(len=27) :: xKClfilename
    character(len=24) :: xClfilename
    character(len=19) :: potentialfilename
    character(len=16) :: chargefilename
    character(len=22) :: xHplusfilename
    character(len=22) :: xOHminfilename
    character(len=22) :: densfracAfilename
    character(len=22) :: densfracBfilename
    character(len=28) :: densfracionpairfilename

    integer :: un_sys,un_xpol,un_xpolC,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair  ! unit numbers

    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals,fmt   

    integer :: i,j,k        ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    !     .. executable statements 

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
  
    !     .. make label filenames 

    write(rstr,'(F5.3)')sigmaAB*delta 
    fnamelabel="sg"//trim(adjustl(rstr)) 
    write(rstr,'(F5.3)')cNaCl
    fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
    write(rstr,'(F5.3)')cCaCl2
    fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
    write(rstr,'(F7.3)')pHbulk
    fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
    write(rstr,'(I5.5)')countfile
    fnamelabel=trim(fnamelabel)//"."//trim(adjustl(rstr))//".dat"

    !     .. make filenames 
    sysfilename='system.'//trim(fnamelabel)
    xpolABfilename='xpolAB.'//trim(fnamelabel)
    xsolfilename='xsol.'//trim(fnamelabel)
    xNafilename='xNaions.'//trim(fnamelabel)
    xKfilename='xKions.'//trim(fnamelabel)
    xCafilename='xCaions.'//trim(fnamelabel)
    xNaClfilename='xNaClionpair.'//trim(fnamelabel)
    xKClfilename='xKClionpair.'//trim(fnamelabel)
    xClfilename='xClions.'//trim(fnamelabel)
    potentialfilename='potential.'//trim(fnamelabel)
    chargefilename='charge.'//trim(fnamelabel)
    xHplusfilename='xHplus.'//trim(fnamelabel)
    xOHminfilename='xOHmin.'//trim(fnamelabel)
    densfracAfilename='densityAfrac.'//trim(fnamelabel)
    densfracBfilename='densityBfrac.'//trim(fnamelabel)
    densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
    ! write(xpolendfilename,'(A8,BZ,I5.5,A4)')'xpolend.', countfile,'.dat'
   
    !     .. opening files        
    open(unit=newunit(un_sys),file=sysfilename)       
    open(unit=newunit(un_xsol),file=xsolfilename)
    open(unit=newunit(un_psi),file=potentialfilename)
    open(unit=newunit(un_xpol),file=xpolABfilename)
    open(unit=newunit(un_fdisA),file=densfracAfilename) 
    open(unit=newunit(un_fdisB),file=densfracBfilename) 
    
    if(verboseflag=="yes") then    
        open(unit=newunit(un_xNa),file=xNafilename)
        open(unit=newunit(un_xK),file=xKfilename)
        open(unit=newunit(un_xCa),file=xCafilename)
        open(unit=newunit(un_xNaCl),file=xNaClfilename)
        open(unit=newunit(un_xKCl),file=xKClfilename)
        open(unit=newunit(un_xpair),file=densfracionpairfilename)
        open(unit=newunit(un_xCl),file=xClfilename)
        open(unit=newunit(un_charge),file=chargefilename)
        open(unit=newunit(un_xHplus),file=xHplusfilename)
        open(unit=newunit(un_xOHmin),file=xOHminfilename)
    endif

    write(un_psi,*)0.0_dp,psiSurfL

    do i=1,nz
         
        write(un_xpol,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
        write(un_xsol,*)zc(i),xsol(i)
        write(un_fdisA,fmt6reals)zc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
        write(un_fdisA,fmt6reals)zc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)   
        write(un_psi,*)zc(i),psi(i)
    enddo
    write(70,*)nz*delta,psiSurfR


    if(verboseflag=="yes") then  
        do i=1,nz   
            write(un_xNa,*)zc(i),xNa(i)
            write(un_xK,*)zc(i),xK(i)
            write(un_xCa,*)zc(i),xCa(i)
            write(un_xNaCl,*)zc(i),xNaCl(i)
            write(un_xKCl,*)zc(i),xKCl(i)
            write(un_xpair,*)zc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)zc(i),xCl(i)
            write(un_psi,*)zc(i),psi(i)
            write(un_charge,*)zc(i),rhoq(i)
            write(un_xHplus,*)zc(i),xHplus(i)
            write(un_xOHmin,*)zc(i),xOHmin(i)
                  
        enddo
    endif  
  
    !     .. system information 


    write(un_xsol,*)'system      = planar weakpolyelectrolyte brush'
    write(un_xsol,*)'version     = ',VERSION
    write(un_xsol,*)'chainmethod = ',chainmethod
    write(un_xsol,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
        write(un_xsol,*)'readinchains = ',readinchains
    endif
    write(un_xsol,*)'systype     = ',systype


    write(un_xsol,*)'free energy = ',FE
    write(un_xsol,*)'energy bulk = ',FEbulk 
    write(un_xsol,*)'deltafenergy = ',deltaFE
    write(un_xsol,*)'fnorm       = ',fnorm
    write(un_xsol,*)'q residual  = ',qres
    write(un_xsol,*)'error       = ',error
    write(un_xsol,*)'sigmaABL    = ',sigmaABL*delta
    write(un_xsol,*)'sigmaABR    = ',sigmaABR*delta
    write(un_xsol,*)'sumphiA     = ',sumphiA
    write(un_xsol,*)'sumphiB     = ',sumphiB
    write(un_xsol,*)'check phi   = ',checkphi 
    write(un_xsol,*)'FEq         = ',FEq 
    write(un_xsol,*)'FEpi        = ',FEpi
    write(un_xsol,*)'FErho       = ',FErho
    write(un_xsol,*)'FEel        = ',FEel
    write(un_xsol,*)'FEelsurf(LEFT) = ',FEelsurf(LEFT)
    write(un_xsol,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_xsol,*)'FEbind      = ',FEbind
    write(un_xsol,*)'FEVdW       = ',FEVdW 
    write(un_xsol,*)'FEalt       = ',FEalt
    write(un_xsol,*)'qABL        = ',qABL
    write(un_xsol,*)'qABR        = ',qABR
    write(un_xsol,*)'muAB        = ',-dlog(qAB)
    write(un_xsol,*)'nsegAB      = ',nsegAB
    write(un_xsol,*)'lsegAB      = ',lsegAB
    write(un_xsol,*)'period      = ',period
    write(un_xsol,*)'nz          = ',nz
    write(un_xsol,*)'delta       = ',delta 
    write(un_xsol,*)'distance    = ',nz*delta   
    write(un_xsol,*)'vsol        = ',vsol
    write(un_xsol,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_xsol,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_xsol,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_xsol,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_xsol,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_xsol,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_xsol,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_xsol,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_xsol,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_xsol,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_xsol,*)'vNa         = ',vNa*vsol
    write(un_xsol,*)'vCl         = ',vCl*vsol
    write(un_xsol,*)'vCa         = ',vCa*vsol
    write(un_xsol,*)'vK          = ',vK*vsol
    write(un_xsol,*)'vNaCl       = ',vNaCl*vsol
    write(un_xsol,*)'vKCl        = ',vKCl*vsol
    write(un_xsol,*)'cNaCl       = ',cNaCl
    write(un_xsol,*)'cKCl        = ',cKCl
    write(un_xsol,*)'cCaCl2      = ',cCaCl2
    write(un_xsol,*)'pHbulk      = ',pHbulk
    write(un_xsol,*)'pKa         = ',pKa(1)      
    write(un_xsol,*)'pKaNa       = ',pKa(2)
    write(un_xsol,*)'pKaACa      = ',pKa(3)
    write(un_xsol,*)'pKaA2Ca     = ',pKa(4)
    write(un_xsol,*)'pKb         = ',pKb(1)      
    write(un_xsol,*)'pKbNa       = ',pKb(2)
    write(un_xsol,*)'pKbBCa      = ',pKb(3)
    write(un_xsol,*)'pKbB2Ca     = ',pKb(4)
    write(un_xsol,*)'KionNa      = ',KionNa
    write(un_xsol,*)'KionK       = ',KionK
    write(un_xsol,*)'K0ionNa     = ',K0ionNa
    write(un_xsol,*)'K0ionK      = ',K0ionK
    write(un_xsol,*)'xNabulk     = ',xbulk%Na
    write(un_xsol,*)'xClbulk     = ',xbulk%Cl
    write(un_xsol,*)'xKbulk      = ',xbulk%K
    write(un_xsol,*)'xNaClbulk   = ',xbulk%NaCl
    write(un_xsol,*)'xKClbulk    = ',xbulk%KCl
    write(un_xsol,*)'xCabulk     = ',xbulk%Ca
    write(un_xsol,*)'xHplusbulk  = ',xbulk%Hplus
    write(un_xsol,*)'xOHminbulk  = ',xbulk%OHmin
    write(un_xsol,*)'sigmaABL    = ',sigmaABL*delta
    write(un_xsol,*)'sigmaABR    = ',sigmaABR*delta
    write(un_xsol,*)'psiSurfL    = ',psiSurfL
    write(un_xsol,*)'psiSurfR    = ',psiSurfR
    write(un_xsol,*)'heightAB    = ',heightAB
    write(un_xsol,*)'qpolA       = ',qpolA
    write(un_xsol,*)'qpolB       = ',qpolB
    write(un_xsol,*)'qpoltot     = ',qpol_tot

    write(un_xsol,*)'avfdisA(1)  = ',avfdisA(1)
    write(un_xsol,*)'avfdisA(2)  = ',avfdisA(2)
    write(un_xsol,*)'avfdisA(3)  = ',avfdisA(3)
    write(un_xsol,*)'avfdisA(4)  = ',avfdisA(4)
    write(un_xsol,*)'avfdisA(5)  = ',avfdisA(5)
    write(un_xsol,*)'avfdisB(1)  = ',avfdisB(1)
    write(un_xsol,*)'avfdisB(2)  = ',avfdisB(2)
    write(un_xsol,*)'avfdisB(3)  = ',avfdisB(3)
    write(un_xsol,*)'avfdisB(4)  = ',avfdisB(4)
    write(un_xsol,*)'avfdisB(5)  = ',avfdisB(5)
    write(un_xsol,*)'dielectW    = ',dielectW
    write(un_xsol,*)'lb          = ',lb
    write(un_xsol,*)'T           = ',T 
    write(un_xsol,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(un_xsol,*)'zpolA(1)    = ',zpolA(1)
    write(un_xsol,*)'zpolA(2)    = ',zpolA(2)
    write(un_xsol,*)'zpolA(3)    = ',zpolA(3)
    write(un_xsol,*)'zpolA(4)    = ',zpolA(4)
    write(un_xsol,*)'zpolB(1)    = ',zpolB(1)
    write(un_xsol,*)'zpolB(2)    = ',zpolB(2)
    write(un_xsol,*)'zpolB(3)    = ',zpolB(3)
    write(un_xsol,*)'zpolB(4)    = ',zpolB(4)
    write(un_xsol,*)'zpolB(5)    = ',zpolB(5)
    write(un_xsol,*)'zNa         = ',zNa
    write(un_xsol,*)'zCa         = ',zCa
    write(un_xsol,*)'zK          = ',zK
    write(un_xsol,*)'zCl         = ',zCl
    write(un_xsol,*)'nsize       = ',nsize  
    write(un_xsol,*)'cuantasAB   = ',cuantasAB
    write(un_xsol,*)'iterations  = ',iter


    ! .. closing files

    close(un_sys)
    close(un_xsol)
    close(un_psi)
    close(un_xpol)   
    close(un_fdisA)
    close(un_fdisB)
    
    if(verboseflag=="yes") then 
        close(un_xNa)   
        close(un_xK)
        close(un_xCa)
        close(un_xNaCl)
        close(un_xKCl)
        close(un_xpair)
        close(un_xCl)
        close(un_charge)
        close(un_xHplus)
        close(un_xOHmin)
    endif
  
end subroutine output_electdouble

subroutine output_neutral(countfile)
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters    
    use field
    use energy
    use myutils, only : newunit
    !     use endpoint
  
    implicit none
    
    !     .. scalar arguments
  
    integer :: countfile
    !     .. output file names       
  
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolCfilename 
    character(len=90) :: xpolendfilename 
    
    integer :: un_sys,un_xpolAB,un_xpolC,un_xsol ! unit numbers


    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

    !     .. local arguments
    
    integer :: i,j,k,n           ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    !     .. executable statements 

    n=nz

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
    
    !     .. make label filenames 

    write(rstr,'(F5.3)')sigmaAB*delta 
    fnamelabel="sg"//trim(adjustl(rstr)) 
    write(rstr,'(F5.3)')VdWepsB
    fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))
    write(rstr,'(I5.5)')countfile
    fnamelabel=trim(fnamelabel)//"."//trim(adjustl(rstr))//".dat"

    !     .. make filenames 
    
    sysfilename='system.'//trim(fnamelabel)
    xpolABfilename='xpolAB.'//trim(fnamelabel)   
    xpolCfilename='xpolC.'//trim(fnamelabel)   
    xsolfilename='xsol.'//trim(fnamelabel)   
    xpolendfilename='xpolend.'//trim(fnamelabel)   
    
    !      .. opening files

    open(unit=newunit(un_sys),file=sysfilename)   
    open(unit=newunit(un_xpolAB),file=xpolABfilename)
    open(unit=newunit(un_xpolC),file=xpolCfilename)
    open(unit=newunit(un_xsol),file=xsolfilename)
        
    do i=1,n    
       write(un_xpolAB,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
       write(un_xpolC,fmt2reals)zc(i),xpolC(i)
       write(un_xsol,fmt2reals)zc(i),xsol(i)
    !     write(40,*)zc(i),endpol(i)
    enddo
    
    !     .. system information 

    write(un_xsol,*)'system      = planar  brush' 
    write(un_xsol,*)'version     = ',VERSION
    write(un_xsol,*)'free energy = ',FE  
    write(un_xsol,*)'energy bulk = ',FEbulk 
    write(un_xsol,*)'deltafenergy = ',deltaFE
    write(un_xsol,*)'FEalt       = ',FEalt
    write(un_xsol,*)'FEconfC     = ',FEconfC
    write(un_xsol,*)'FEconfAB    = ',FEconfAB
    write(un_xsol,*)'FEtrans%sol = ',FEtrans%sol  
    write(un_xsol,*)'fnorm       = ',fnorm
    write(un_xsol,*)'error       = ',error
    write(un_xsol,*)'sumphiA     = ',sumphiA
    write(un_xsol,*)'sumphiB     = ',sumphiB
    write(un_xsol,*)'sumphiC     = ',sumphiC
    write(un_xsol,*)'check phi   = ',checkphi 
    write(un_xsol,*)'FEq         = ',FEq 
    write(un_xsol,*)'FEpi        = ',FEpi
    write(un_xsol,*)'FErho       = ',FErho
    write(un_xsol,*)'FEVdW       = ',FEVdW
    write(un_xsol,*)'qAB         = ',qAB
    write(un_xsol,*)'qC          = ',qC
    write(un_xsol,*)'muAB        = ',-dlog(qAB)
    write(un_xsol,*)'muC         = ',-dlog(qC)
    write(un_xsol,*)'nsegAB      = ',nsegAB
    write(un_xsol,*)'lsegAB      = ',lsegAB
    write(un_xsol,*)'nsegC       = ',nsegC
    write(un_xsol,*)'lsegC       = ',lsegC
    write(un_xsol,*)'period      = ',period
    write(un_xsol,*)'nz          = ',nz
    write(un_xsol,*)'delta       = ',delta
    write(un_xsol,*)'distance    = ',nz*delta   
    write(un_xsol,*)'vsol        = ',vsol
    write(un_xsol,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(un_xsol,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(un_xsol,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(un_xsol,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(un_xsol,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(un_xsol,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(un_xsol,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(un_xsol,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(un_xsol,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(un_xsol,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(un_xsol,*)'vpolC       = ',vpolC*vsol
    write(un_xsol,*)'T           = ',T
    write(un_xsol,*)'VdWepsC     = ',VdWepsC*vpolC*vsol
    write(un_xsol,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(un_xsol,*)'heightAB    = ',heightAB
    write(un_xsol,*)'heightC     = ',heightC
    write(un_xsol,*)'nsize       = ',nsize  
    write(un_xsol,*)'cuantasAB   = ',cuantasAB
    write(un_xsol,*)'cuantasC    = ',cuantasC
    write(un_xsol,*)'iterations  = ',iter
    write(un_xsol,*)'chainmethod = ',chainmethod
    write(un_xsol,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
       write(un_xsol,*)'readinchains = ',readinchains
        endif
    write(un_xsol,*)'systype     = ',systype


    close(un_xsol)
    close(un_xpolAB)
    close(un_xpolC)
    close(un_sys)
  
end subroutine output_neutral

subroutine output(countfile)

    use globals, only : systype
    implicit none

    integer :: countfile 

    if(systype=="elect") then 
        call output_elect(countfile)
    elseif(systype=="electdouble") then
        call output_electdouble(countfile)
    elseif(systype=="neutral") then
        call output_neutral(countfile)
    elseif(systype=="electnopoly") then
        call output_elect(countfile)
        call output_individualcontr_fe(countfile)
    else
        print*,"Error in output subroutine"
        print*,"Wrong value systype : ", systype
    endif     

end subroutine output


subroutine output_individualcontr_fe(countfile)

    use globals, only : LEFT,RIGHT, systype
    use energy
    use myutils, only : newunit
    use parameters, only : sigmaAB,cNaCl,cCaCl2,pHbulk,VdWepsB
    use volume, only : delta

    implicit none 

    integer, intent(in) :: countfile 
    
    integer :: un_fe

    character(len=100) :: fenergyfilename   
    character(len=100) :: fnamelabel
    character(len=20) :: rstr


    !     .. make label filename
    
    if(systype=="elect".or.systype=="electdouble".or.systype=="electnopoly") then 
        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')cCaCl2
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
        write(rstr,'(I5.5)')countfile
        fnamelabel=trim(fnamelabel)//"."//trim(adjustl(rstr))//".dat"
    elseif(systype=="neutral") then 
        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')VdWepsB
        fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))
        write(rstr,'(I5.5)')countfile
        fnamelabel=trim(fnamelabel)//"."//trim(adjustl(rstr))//".dat"
    else
        print*,"Error in output_individualcontr_fe subroutine"
        print*,"Wrong value systype : ", systype
    endif    

    fenergyfilename='energy.'//trim(fnamelabel)   
        
    !     .. opening files        

    open(unit=newunit(un_fe),file=fenergyfilename) 

    write(un_fe,*)'FE              = ',FE  
    write(un_fe,*)'FEbulk          = ',FEbulk 
    write(un_fe,*)'deltaFE         = ',deltaFE
    write(un_fe,*)'FEalt           = ',FEalt  
    write(un_fe,*)'FEbulkalt       = ',FEbulkalt 
    write(un_fe,*)'deltaFEalt      = ',deltaFEalt
    
    write(un_fe,*)"FEtrans%sol     = ",FEtrans%sol   
    write(un_fe,*)"FEtrans%Na      = ",FEtrans%Na  
    write(un_fe,*)"FEtrans%Cl      = ",FEtrans%Cl  
    write(un_fe,*)"FEtrans%Ca      = ",FEtrans%Ca  
    write(un_fe,*)"FEtrans%K       = ",FEtrans%K
    write(un_fe,*)"FEtrans%KCl     = ",FEtrans%KCl
    write(un_fe,*)"FEtrans%NaCl    = ",FEtrans%NaCl  
    write(un_fe,*)"FEtrans%Hplus   = ",FEtrans%Hplus  
    write(un_fe,*)"FEtrans%OHmin   = ",FEtrans%OHmin  
    
    write(un_fe,*)"FEchempot%Na    = ",FEchempot%Na
    write(un_fe,*)"FEchempot%Cl    = ",FEchempot%Cl
    write(un_fe,*)"FEchempot%Ca    = ",FEchempot%Ca
    write(un_fe,*)"FEchempot%K     = ",FEchempot%K
    write(un_fe,*)"FEchempot%KCl   = ",FEchempot%KCl
    write(un_fe,*)"FEchempot%NaCl  = ",FEchempot%NaCl
    write(un_fe,*)"FEchempot%Hplus = ",FEchempot%Hplus
    write(un_fe,*)"FEchempot%OHmin = ",FEchempot%OHmin

    write(un_fe,*)"FEchemsurf(LEFT)= ",FEchemsurf(LEFT)
    write(un_fe,*)"FEchemsurf(RIGHT)= ",FEchemsurf(RIGHT)
    write(un_fe,*)"FEchemsurfalt(LEFT)= ",FEchemsurfalt(LEFT)
    write(un_fe,*)"FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)

    write(un_fe,*)"delta FEchemsurfalt(LEFT)= ",FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
    write(un_fe,*)"delta FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)


    close(un_fe)


end subroutine   output_individualcontr_fe

end module


