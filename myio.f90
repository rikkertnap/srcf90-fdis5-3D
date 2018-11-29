module myio

    use precision_definition
    implicit none

    ! return error values

    integer, parameter ::  myio_err_sysflag   = 1
    integer, parameter ::  myio_err_runflag   = 2
    integer, parameter ::  myio_err_geometry  = 3
    integer, parameter ::  myio_err_method    = 4
    integer, parameter ::  myio_err_chaintype = 5
    integer, parameter ::  myio_err_domain    = 6
    integer, parameter ::  myio_err_inputfile = 7
    integer, parameter ::  myio_err_input     = 8
    integer, parameter ::  myio_err_bcflag    = 9 
    integer, parameter ::  myio_err_label     = 10 

    ! unit number 
    integer :: un_sys,un_xpolAB,un_xpolC,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair, un_rhopolAB, un_fe, un_q
    integer :: un_dip ,un_dielec,un_xpolABz  
   
    ! format specifiers 
    character(len=80), parameter  :: fmt = "(A9,I1,A5,ES25.16)"
    character(len=80), parameter  :: fmt1reals = "(ES25.16E3)"   
    character(len=80), parameter  :: fmt2reals = "(2ES25.16E3)"   
    character(len=80), parameter  :: fmt3reals = "(3ES25.16E3)"  
    character(len=80), parameter  :: fmt4reals = "(4ES25.16E3)" 
    character(len=80), parameter  :: fmt5reals = "(5ES25.16E3)"
    character(len=80), parameter  :: fmt6reals = "(6ES25.16E3)" 
    
    private
    public :: read_inputfile, output_individualcontr_fe, output, compute_vars_and_output
    
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
    character(len=100) :: buffer, label
    integer :: pos
    integer :: line



    if (present(info)) info = 0
    
    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    ios=0 
    line = 0

    ! ios<0 : if an end of record condition is encountered or if an endfile condition was detected.  
    ! ios>0 : if an error occured 
    ! ios=0 : otherwise.

    do while (ios == 0)

        read(un_input, '(A)', iostat=ios) buffer

        if (ios == 0) then
        
            line = line + 1

            !  Split label and data based on first occurence of a whitespace
            pos = scan(buffer, '     ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            select case (label) !list-directed The CHARACTER variable is treated as an 'internal file'
            case ('method')
                read(buffer, *,iostat=ios) method
            case ('sysflag')
                read(buffer, *,iostat=ios) sysflag
            case ('runflag')
                read(buffer, *,iostat=ios) runflag
            case ('bcflag(LEFT)')
                read(buffer,*,iostat=ios) bcflag(LEFT)
            case ('bcflag(RIGHT)')
                read(buffer,*,iostat=ios) bcflag(RIGHT)
            case ('chainmethod')
                read(buffer,*,iostat=ios) chainmethod
            case ('chaintype')
                read(buffer,*,iostat=ios) chaintype
            case ('tolerance')
                read(buffer,*,iostat=ios) error             
            case ('infile')
                read(buffer,*,iostat=ios) infile              ! guess  1==yes
            case ('pH%val')
                read(buffer,*,iostat=ios) pH%val
            case ('pH%min')
                read(buffer,*,iostat=ios) pH%min
            case ('pH%max')
                read(buffer,*,iostat=ios) pH%max
            case ('pH%stepsize')
                read(buffer,*,iostat=ios) pH%stepsize
            case ('pH%delta')
                read(buffer,*,iostat=ios) pH%delta
            case ('KionNa')
                read(buffer,*,iostat=ios) KionNa
            case ('KionK')
                read(buffer,*,iostat=ios) KionK
            case ('sigmaSurfL')
                read(buffer,*,iostat=ios) sigmaSurfL
            case ('sigmaSurfR')
                read(buffer,*,iostat=ios) sigmaSurfR
            case ('cNaCl')
                read(buffer,*,iostat=ios) cNaCl
            case ('cKCl')
                read(buffer,*,iostat=ios) cKCl
            case ('cCaCl2')
                read(buffer,*,iostat=ios) cCaCl2
            case ('pKa(1)')
                read(buffer,*,iostat=ios) pKa(1)           !   AH   <=> A- + H+ 
            case ('pKa(2)')
                read(buffer,*,iostat=ios) pKa(2)           !   ANa  <=> A- + Na+  
            case ('pKa(3)')
                read(buffer,*,iostat=ios) pKa(3)           !   ACa+ <=> A- + Ca2+ 
            case ('pKa(4)')
                read(buffer,*,iostat=ios) pKa(4)           !   A2Ca <=> 2A- + Ca2+
            case ('pKb(1)')
                read(buffer,*,iostat=ios) pKb(1)           !   BH   <=> B- + H+ 
            case ('pKb(2)')
                read(buffer,*,iostat=ios) pKb(2)           !   BNa  <=> B- + Na+ 
            case ('pKb(3)')
                read(buffer,*,iostat=ios) pKb(3)           !   BCa+ <=> B- + Ca2+   
            case ('pKb(4)')
                read(buffer,*,iostat=ios) pKb(4)           !   B2Ca <=> 2B- + Ca2+   
            case ('chainperiod')
                read(buffer,*,iostat=ios) period
            case ('nsize')
                read(buffer,*,iostat=ios) nsize
            case ('nsegAB')
                read(buffer,*,iostat=ios) nsegAB
            case ('cuantasAB')
                read(buffer,*,iostat=ios) cuantasAB
            case ('nsegC')
                read(buffer,*,iostat=ios) nsegC
            case ('cuantasC')
                read(buffer,*,iostat=ios) cuantasC
            case ('VdWepsC')
                read(buffer,*,iostat=ios) VdWepsC
            case ('VdWepsB ')
                read(buffer,*,iostat=ios) VdWepsB    
            case ('VdWcutoff')                             
                read(buffer,*,iostat=ios) VdWcutoff  
            case ('nx')
                read(buffer,*,iostat=ios) nx
            case ('ny')
                read(buffer,*,iostat=ios) ny
            case ('nzmax')
                read(buffer,*,iostat=ios) nzmax
            case ('nzmin')
                read(buffer,*,iostat=ios) nzmin
            case ('nzstep')
                read(buffer,*,iostat=ios) nzstep
            case ('verboseflag  ')
                read(buffer,*,iostat=ios) verboseflag  
            case ('delta')
                read(buffer,*,iostat=ios) delta   
            case ('geometry')
                read(buffer,*,iostat=ios) geometry
            case ('ngr_freq')
                read(buffer,*,iostat=ios) ngr_freq 
            case ('gamma')
                read(buffer,*,iostat=ios) gamma  
            case ('isRandom_pos_graft')
                read(buffer,*,iostat=ios) isRandom_pos_graft  
            case default
                if(pos>1) then 
                    print *, 'Invalid label at line', line  ! empty lines are skipped
                endif
            end select
        end if
    end do


    if(ios >0 ) then
        print*, 'Error parsing file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    close(un_input)
    
    ! override input bcflags 
    
    if(sysflag=="electdouble") then
        bcflag(LEFT)="cc"
        bcflag(RIGHT)="cc"  
    endif
    if(sysflag=="elect") then
        sigmaAB=sigmaABL
        sigmaABR=0.0_dp  
    endif
    
    ! .. check error flag

    call check_value_sysflag(sysflag,info_sys) 
    if (info_sys == myio_err_sysflag) then
        if (present(info)) info = info_sys
        return
    endif

    call check_value_runflag(runflag,info_run) 
    if (info_sys == myio_err_runflag) then
        if (present(info)) info = info_run
        return
    endif

    call check_value_geometry(geometry,info_geo)
    if (info_geo == myio_err_geometry) then
        if (present(info)) info = info_geo
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

    !  potentially overide input values

    call set_value_nzmin(runflag,nzmin,nzmax)

end subroutine read_inputfile
 

subroutine check_value_sysflag(sysflag,info)

    implicit none

    character(len=15), intent(in) :: sysflag
    integer, intent(out),optional :: info

    character(len=15) :: sysflagstr(9)
    integer :: i
    logical :: flag

    ! permissible values of sysflag

    sysflagstr(1)="elect"
    sysflagstr(2)="bulk water"
    sysflagstr(3)="neutral"
    sysflagstr(4)="electdouble"
    sysflagstr(5)="electnopoly"
    sysflagstr(6)="electHC"
    sysflagstr(7)="dipolarweak"
    sysflagstr(8)="dipolarstrong"
    sysflagstr(9)="electA"

    flag=.FALSE.

    do i=1,9
        if(sysflag==sysflagstr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of sysflag is not permissible"
        print*,"sysflag = ",sysflag
        if (present(info)) info = myio_err_sysflag
        return
    end if

end subroutine check_value_sysflag


subroutine check_value_runflag(runflag,info)

    implicit none

    character(len=15), intent(in) :: runflag
    integer, intent(out),optional :: info

    character(len=15) :: runflagstr(2)
    integer :: i
    logical :: flag

    ! permissible values of runflag

    runflagstr(1)="rangepH"
    runflagstr(2)="rangedist"

    flag=.FALSE.

    do i=1,2
        if(runflag==runflagstr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of runflag is not permissible"
        print*,"runflag = ",runflag
        if (present(info)) info = myio_err_runflag
        return
    end if

end subroutine check_value_runflag

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
        if (present(info)) info = myio_err_sysflag
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

subroutine check_value_geometry(geometry,info)
        
    implicit none

    character(len=11), intent(in) :: geometry
    integer, intent(out),optional :: info

    logical :: flag
    character(len=11) :: geometrystr(2) 
    integer :: i

    ! permissible values of geometry

    geometrystr(1)="cubic"
    geometrystr(2)="prism"
   
    flag=.FALSE.

    do i=1,2
        if(geometry==geometrystr(i)) flag=.TRUE.
    enddo
        
    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then 
        print*,"Error: value of geometry is not permissible"
        print*,"geometry = ",geometry
        if (present(info)) info = myio_err_geometry 
        return
    endif
    
end subroutine


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

    ! permissible values of runflag

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


! override input value nzmin
! Sets nzmin=nzmax for which runflag that do not loop over nz 
! ensuring that the output files are properly close once

subroutine set_value_nzmin(runflag,nzmin,nzmax)

    character(len=15), intent(in) :: runflag
    integer, intent(inout) :: nzmin
    integer, intent(in) :: nzmax
   
    if (runflag/="rangedist") nzmin=nzmax 

end subroutine



subroutine output()

    use globals, only : sysflag
    implicit none

    if(sysflag=="elect") then 
        call output_elect
        call output_individualcontr_fe
    elseif(sysflag=="electA") then
        call output_elect
        call output_individualcontr_fe
    elseif(sysflag=="electdouble") then
        call output_electdouble
        call output_individualcontr_fe
    elseif(sysflag=="neutral") then
        call output_neutral
    elseif(sysflag=="electnopoly") then
        call output_elect
        call output_individualcontr_fe
    elseif(sysflag=="dipolarstrong") then
        call output_elect
    elseif(sysflag=="dipolarweak") then
        call output_elect
    else
        print*,"Error in output subroutine"
        print*,"Wrong value sysflag : ", sysflag
    endif     

end subroutine output


subroutine output_elect
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    use myutils, only : newunit
    use chains, only : isHomopolymer
  
    implicit none
    
    !     .. local arguments
   
    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolABzfilename 
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
    character(len=90) :: qfilename
    character(len=90) :: densfracionpairfilename
    character(len=90)  :: dipolefilename
    character(len=100) :: dielecfilename
    integer :: i,j,k          ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen

    ! .. executable statements 

    if(nz.eq.nzmax)  then 
        !     .. make label filenames 

        write(rstr,'(F5.3)')sigmaABL
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        if(cCaCl2>=0.001) then 
            write(rstr,'(F5.3)')cCaCl2
        elseif(cCaCl2>0.0) then  
            write(rstr,'(ES9.2E2)')cCaCl2
        else 
            write(rstr,'(F3.1)')cCaCl2
        endif 
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"

        sysfilename='system.'//trim(fnamelabel)
        xpolABfilename='xpolAB.'//trim(fnamelabel)
        xpolCfilename='xpolC.'//trim(fnamelabel)
        xpolABzfilename='xpolABz.'//trim(fnamelabel)
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
        qfilename='q.'//trim(fnamelabel)
        dipolefilename='dipole.'//trim(fnamelabel)//'.dat'
        dielecfilename='dielectrel.'//trim(fnamelabel)//'.dat'

        !     .. opening files        
        
        open(unit=newunit(un_sys),file=sysfilename)       
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_psi),file=potentialfilename)

        if(sysflag/="electnopoly") then          
            open(unit=newunit(un_xpolAB),file=xpolABfilename)
            open(unit=newunit(un_xpolC),file=xpolCfilename)
            open(unit=newunit(un_fdisA),file=densfracAfilename) 
            open(unit=newunit(un_fdisB),file=densfracBfilename) 
            open(unit=newunit(un_q),file=qfilename)
            open(unit=newunit(un_xpolABz),file=xpolABzfilename)
        endif     
        if(sysflag=="dipolarweak".or.sysflag=="dipolarstrong") then          
            open(unit=newunit(un_dielec),file=dielecfilename)
            open(unit=newunit(un_dip),file=dipolefilename)
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

    else ! check that files are open
        inquire(unit=un_sys, opened=isopen)
        inquire(unit=un_xpolAB, opened=isopen)
        inquire(unit=un_xpolC, opened=isopen)
        inquire(unit=un_xsol, opened=isopen)

        if(.not.isopen) write(*,*)"un_xsol is not open"

    endif       

    !   .. writting files
    !   .. this line seperates different distances 
      
    write(un_xsol,*)'#D    = ',nz*delta 
    write(un_psi,*)'#D    = ',nz*delta
    if(sysflag/="electnopoly") then         
        write(un_xpolAB,*)'#D    = ',nz*delta 
        write(un_xpolC,*)'#D    = ',nz*delta 
        write(un_fdisA,*)'#D    = ',nz*delta
        write(un_fdisB,*)'#D    = ',nz*delta
    endif      
    if(verboseflag=="yes") then    
        write(un_xNa,*)'#D    = ',nz*delta 
        write(un_xK,*)'#D    = ',nz*delta 
        write(un_xCa,*)'#D    = ',nz*delta 
        write(un_xNaCl,*)'#D    = ',nz*delta 
        write(un_xKCL,*)'#D    = ',nz*delta 
        write(un_xpair,*)'#D    = ',nz*delta 
        write(un_charge,*)'#D    = ',nz*delta 
        write(un_xCl,*)'#D    = ',nz*delta 
        write(un_xHplus,*)'#D    = ',nz*delta 
        write(un_xOHMin,*)'#D    = ',nz*delta 
    endif

    do i=1,nsurf
        write(un_psi,*)psiSurfL(i)
    enddo 
        
    do i=1,nsize
        write(un_xsol,*)xsol(i)
        write(un_psi,*)psi(i)
    enddo    

    do i=1,nsurf
        write(un_psi,*)psiSurfR(i)
    enddo         

    if(sysflag/="electnopoly") then 

        do i =1, ngr
            write(un_q,*)qABL(i),qABR(i)
        enddo 

        do i=1,nsize
            write(un_xpolAB,fmt3reals)xpolAB(i),rhopolA(i),rhopolB(i)
            write(un_xpolC,fmt1reals)xpolC(i)
            write(un_fdisA,fmt5reals)fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
            write(un_fdisB,fmt5reals)fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)
        enddo
        do i=1,nz
            write(un_xpolABz,fmt1reals)xpolABz(i)
        enddo    
    endif   
    
    if(verboseflag=="yes") then 
        do i=1,nsize
            write(un_xNa,*)xNa(i)
            write(un_xK,*)xK(i)
            write(un_xCa,*)xCa(i)
            write(un_xNaCl,*)xNaCl(i)
            write(un_xKCl,*)xKCl(i)
            write(un_xpair,*)(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)xCl(i)
            write(un_charge,*)rhoq(i)
            write(un_xHplus,*)xHplus(i)
            write(un_xOHmin,*)xOHmin(i)    
        enddo    
    endif

    if(sysflag=="dipolarweak".or.sysflag=="dipolarstrong") then          
        do i = 1, nsize
            write(un_dip,*)rhob(i),electPol(1,i),electPol(2,i),electPol(3,i)
        !    write(un_dielec,*)dielectrel(i)
        enddo   
    endif

    ! .. writing system information 
    if(nz.eq.nzmax) then 
        write(un_sys,*)'===begin distance independent settings=='  
        write(un_sys,*)'system      = planar weakpolyelectrolyte brush'
        write(un_sys,*)'version     = ',VERSION
        ! chain description 
        write(un_sys,*)'chainmethod = ',chainmethod
        write(un_sys,*)'chaintype   = ',chaintype
        write(un_sys,*)"isHomopolymer= ",isHomopolymer
        if(chainmethod.eq."FILE") then
           write(un_sys,*)'readinchains = ',readinchains
        endif
        write(un_sys,*)'nsegAB      = ',nsegAB
        write(un_sys,*)'lsegAB      = ',lsegAB
        write(un_sys,*)'nsegC       = ',nsegC
        write(un_sys,*)'lsegC       = ',lsegC
        write(un_sys,*)'period      = ',period
        write(un_sys,*)'cuantasAB   = ',cuantasAB
        write(un_sys,*)'sigmaAB     = ',sigmaAB
        write(un_sys,*)'sigmaC      = ',sigmaC

        ! system description
        write(un_sys,*)'sysflag     = ',sysflag
        write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
        write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
        write(un_sys,*)'delta       = ',delta
        write(un_sys,*)'ngr         = ',ngr  
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conf    = ',error
        ! concentration 
        write(un_sys,*)'cNaCl       = ',cNaCl
        write(un_sys,*)'cKCl        = ',cKCl
        write(un_sys,*)'cCaCl2      = ',cCaCl2
        write(un_sys,*)'xNabulk     = ',xbulk%Na
        write(un_sys,*)'xClbulk     = ',xbulk%Cl
        write(un_sys,*)'xKbulk      = ',xbulk%K
        write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
        write(un_sys,*)'xKClbulk    = ',xbulk%KCl
        write(un_sys,*)'xCabulk     = ',xbulk%Ca
        write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
        write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
        write(un_sys,*)'pHbulk      = ',pHbulk
        ! disociation constants 
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
        ! other physcial parameters
        if(sysflag=="dipolarweak".or.sysflag=="dipolarstrong") then 
            write(un_sys,*)'dipole pol  = ',dipole%pol
            write(un_sys,*)'dipole sol  = ',dipole%sol
        endif    

        write(un_sys,*)'dielectW    = ',dielectW
        write(un_sys,*)'lb          = ',lb
        write(un_sys,*)'T           = ',Tref
        write(un_sys,*)'VdWepsC     = ',VdWepsC*vpolC*vsol 
        write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
        ! charge components
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
         ! volume  
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
    
        write(un_sys,*)'===end distance independent settings=='
    endif
    write(un_sys,*)'D plates    = ',nz*delta 
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'error       = ',error
    write(un_sys,*)'sigmaAB     = ',sigmaAB
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
    write(un_sys,*)'sigmaSurfL  = ',sigmaSurfL/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaSurfR  = ',sigmaSurfR/((4.0_dp*pi*lb)*delta)
    !write(un_sys,*)'sigmaqSurfL = ',sigmaqSurfL/((4.0_dp*pi*lb)*delta)
    !write(un_sys,*)'sigmaqSurfR = ',sigmaqSurfR/((4.0_dp*pi*lb)*delta)
   
    if(bcflag(LEFT)=='ta') then
      do i=1,4   
        write(un_sys,fmt)'fdisTaL(',i,')  = ',fdisTaL(i)
      enddo  
    endif
    if(bcflag(RIGHT)=='ta') then   
      do i=1,4   
        write(un_sys,fmt)' fdisTaR(',i,')  = ',fdisTaR(i)
      enddo  
    endif
    if((bcflag(RIGHT)/='ta').and.(bcflag(RIGHT)/='cc') ) then
      do i=1,6   
        write(un_sys,fmt)' fdisSuR(',i,')  = ',fdisS(i)
      enddo  
    endif
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    write(un_sys,*)'iterations  = ',iter
   
    ! .. closing files

    if(nz==nzmin) then 
        close(un_sys)
        close(un_xsol)
        close(un_psi)
        if(sysflag/="electnopoly") then
            close(un_xpolAB)   
            close(un_xpolC)
            close(un_fdisA)
            close(un_fdisB)
            close(un_xpolABz)
            close(un_q)
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
        if(sysflag=="dipolarstrong".or.sysflag=="dipolarweak") then
            close(un_dip)    
         !   close(un_dielec)
        endif    

    endif
        

end subroutine output_elect

 
subroutine output_electdouble()
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface
    use myutils, only : newunit
    use chains, only : isHomopolymer
    !     use endpoint
  
    implicit none
    
    !     .. scalar arguments
    
    !     .. local arguments

    !     .. output file names       

    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename
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
    character(len=90) :: rhopolABfilename

    integer :: i,j,k        ! dummy indexes
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    !     .. executable statements 

    if(nz.eq.nzmax) then 

        !     .. make label filenames 
        write(rstr,'(F5.3)')sigmaABL
        fnamelabel="sgL"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')sigmaABR
        fnamelabel=trim(fnamelabel)//"sgR"//trim(adjustl(rstr))  
        write(rstr,'(F5.3)')cNaCl
        fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
        if(cCaCl2>=0.001) then 
            write(rstr,'(F5.3)')cCaCl2
        elseif(cCaCl2>0.0) then  
            write(rstr,'(ES9.2E2)')cCaCl2
        else 
            write(rstr,'(F3.1)')cCaCl2
        endif    
        fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"

        ! print*,fnamelabel

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
        rhopolABfilename='rhopolAB.'//trim(fnamelabel)
       
        !     .. opening files        
        open(unit=newunit(un_sys),file=sysfilename)       
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_psi),file=potentialfilename)
        open(unit=newunit(un_xpolAB),file=xpolABfilename)
        open(unit=newunit(un_fdisA),file=densfracAfilename) 
        open(unit=newunit(un_fdisB),file=densfracBfilename) 
        open(unit=newunit(un_rhopolAB),file=rhopolABfilename)    

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
    endif    
     
    !   .. this line seperates different distances 
      
    write(un_xsol,*)'#D    = ',nz*delta 
    write(un_psi,*)'#D    = ',nz*delta    
    write(un_xpolAB,*)'#D    = ',nz*delta 
    write(un_fdisA,*)'#D    = ',nz*delta    
    write(un_fdisB,*)'#D    = ',nz*delta 
    write(un_rhopolAB,*)'#D    = ',nz*delta 
     
   
    if(verboseflag=="yes") then    
        write(un_xNa,*)'#D    = ',nz*delta 
        write(un_xK,*)'#D    = ',nz*delta 
        write(un_xCa,*)'#D    = ',nz*delta 
        write(un_xNaCl,*)'#D    = ',nz*delta 
        write(un_xKCl,*)'#D    = ',nz*delta 
        write(un_xpair,*)'#D    = ',nz*delta 
        write(un_xCl,*)'#D    = ',nz*delta 
        write(un_charge,*)'#D    = ',nz*delta 
        write(un_xHplus,*)'#D    = ',nz*delta 
        write(un_xOHMin,*)'#D    = ',nz*delta 
    endif    

    do i=1,nsurf
        write(un_psi,*)psiSurfL(i)
    enddo     
    do i=1,nsize
        write(un_xpolAB,fmt3reals) xpolAB(i),rhopolA(i),rhopolB(i)
        write(un_xsol,fmt1reals) xsol(i)
        write(un_fdisA,fmt5reals) fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
        write(un_fdisB,fmt5reals) fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)   
        write(un_psi,fmt1reals) psi(i)
        write(un_rhopolAB,fmt4reals) rhopolAL(i),rhopolBL(i),rhopolAR(i),rhopolBR(i)
    enddo

    do i=1,nsurf
        write(un_psi,*)psiSurfR(i)
    enddo     


    if(verboseflag=="yes") then  
        do i=1,nsize   
            write(un_xNa,*)xNa(i)
            write(un_xK,*)xK(i)
            write(un_xCa,*)xCa(i)
            write(un_xNaCl,*)xNaCl(i)
            write(un_xKCl,*)xKCl(i)
            write(un_xpair,*)(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(un_xCl,*)xCl(i)
            write(un_charge,*)rhoq(i)
            write(un_xHplus,*)xHplus(i)
            write(un_xOHmin,*)xOHmin(i)           
        enddo
    endif  
  
    !     .. system information 
    if(nz.eq.nzmax) then 

        write(un_sys,*)'===begin distance independent settings=='  
        write(un_sys,*)'system      = planar weakpolyelectrolyte brush'
        write(un_sys,*)'version     = ',VERSION
        ! chain description 
        write(un_sys,*)'chainmethod = ',chainmethod
        write(un_sys,*)'chaintype   = ',chaintype
        write(un_sys,*)"isHomopolymer= ",isHomopolymer
        if(chainmethod.eq."FILE") then
           write(un_sys,*)'readinchains = ',readinchains
        endif
        write(un_sys,*)'nsegAB      = ',nsegAB
        write(un_sys,*)'lsegAB      = ',lsegAB
        write(un_sys,*)'period      = ',period
        write(un_sys,*)'cuantasAB   = ',cuantasAB
        write(un_sys,*)'sigmaAB     = ',sigmaAB
        
        
        ! system description
        write(un_sys,*)'sysflag     = ',sysflag
        write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
        write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
        write(un_sys,*)'delta       = ',delta  
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conf    = ',error
        ! concentration 
         write(un_sys,*)'cNaCl       = ',cNaCl
        write(un_sys,*)'cKCl        = ',cKCl
        write(un_sys,*)'cCaCl2      = ',cCaCl2
        write(un_sys,*)'pHbulk      = ',pHbulk
        write(un_sys,*)'xNabulk     = ',xbulk%Na
        write(un_sys,*)'xClbulk     = ',xbulk%Cl
        write(un_sys,*)'xKbulk      = ',xbulk%K
        write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
        write(un_sys,*)'xKClbulk    = ',xbulk%KCl
        write(un_sys,*)'xCabulk     = ',xbulk%Ca
        write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
        write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
        ! disociation constant 
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
        !
        write(un_sys,*)'sigmaABL    = ',sigmaABL
        write(un_sys,*)'sigmaABR    = ',sigmaABR
        write(un_sys,*)'dielectW    = ',dielectW
        write(un_sys,*)'lb          = ',lb
        write(un_sys,*)'T           = ',Tref 
        write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
        ! charge
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
        ! volume   
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
        write(un_sys,*)'vNa         = ',vNa*vsol
        write(un_sys,*)'vCl         = ',vCl*vsol
        write(un_sys,*)'vCa         = ',vCa*vsol
        write(un_sys,*)'vK          = ',vK*vsol
        write(un_sys,*)'vNaCl       = ',vNaCl*vsol
        write(un_sys,*)'vKCl        = ',vKCl*vsol


        write(un_sys,*)'===end distance independent settings=='
    endif
    write(un_sys,*)'D   plates  = ',nz*delta 
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'free energy = ',FE
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'q residual  = ',qres
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'check phi   = ',checkphi 
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEel        = ',FEel
    write(un_sys,*)'FEelsurf(LEFT) = ',FEelsurf(LEFT)
    write(un_sys,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(un_sys,*)'FEbind      = ',FEbind
    write(un_sys,*)'FEVdW       = ',FEVdW 
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'qABL        = ',qABL
    write(un_sys,*)'qABR        = ',qABR
    write(un_sys,*)'muAB        = ',-log(qAB)
    write(un_sys,*)'qpolA       = ',qpolA
    write(un_sys,*)'qpolB       = ',qpolB
    write(un_sys,*)'qpoltot     = ',qpol_tot
    write(un_sys,*)'psiSurfL    = ',psiSurfL
    write(un_sys,*)'psiSurfR    = ',psiSurfR
    write(un_sys,*)'heightAB    = ',heightAB
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
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'iterations  = ',iter

    ! .. closing files
    if(nz.eq.nzmin) then 
        close(un_sys)
        close(un_xsol)
        close(un_psi)
        close(un_xpolAB)   
        close(un_fdisA)
        close(un_fdisB)
        close(un_rhopolAB)
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
    endif    

end subroutine output_electdouble

subroutine output_neutral
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters    
    use field
    use energy
    use myutils, only : newunit
    use chains, only : isHomopolymer
   
  
    implicit none

    !     .. output file names         
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolABfilename 
    character(len=90) :: xpolCfilename 
    character(len=90) :: xpolendfilename 
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

    !     .. local arguments
    integer :: i
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    !     .. executable statements 

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
    
    if(nz==nzmax) then 

        !     .. make label filenames 
        write(rstr,'(F5.3)')sigmaAB*delta 
        fnamelabel="sg"//trim(adjustl(rstr)) 
        write(rstr,'(F5.3)')VdWepsB
        fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))//".dat"

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
        
    else ! check that files are open
        inquire(unit=un_sys, opened=isopen)
        if(.not.isopen) write(*,*)"un_sys is not open" 
        inquire(unit=un_xpolAB, opened=isopen)
        if(.not.isopen) write(*,*)"un_xpolAB is not open"
        inquire(unit=un_xpolC, opened=isopen)
        if(.not.isopen) write(*,*)"un_xpolC is not open"
        inquire(unit=un_xsol, opened=isopen)
        if(.not.isopen) write(*,*)"un_xsol is not open"
    endif    
        

    do i=1,nsize    
       write(un_xpolAB,fmt3reals)xpolAB(i),rhopolA(i),rhopolB(i)
       write(un_xpolC,fmt1reals)xpolC(i)
       write(un_xsol,fmt1reals)xsol(i)
    !     write(40,*)zc(i),endpol(i)
    enddo
        
    !     .. system information 

    if(nz.eq.nzmax) then 
        write(un_sys,*)'===begin distance independent settings=='  
        write(un_sys,*)'system      = planar  brush' 
        write(un_sys,*)'version     = ',VERSION
        ! chain description 
        write(un_sys,*)'chainmethod = ',chainmethod
        write(un_sys,*)'chaintype   = ',chaintype
        write(un_sys,*)"isHomopolymer= ",isHomopolymer
        if(chainmethod.eq."FILE") then
           write(un_sys,*)'readinchains = ',readinchains
        endif
        write(un_sys,*)'nsegAB      = ',nsegAB
        write(un_sys,*)'lsegAB      = ',lsegAB
        write(un_sys,*)'nsegC       = ',nsegC
        write(un_sys,*)'lsegC       = ',lsegC
        write(un_sys,*)'period      = ',period
        write(un_sys,*)'cuantasAB   = ',cuantasAB
        write(un_sys,*)'cuantasC    = ',cuantasC
        ! system description 
        write(un_sys,*)'sysflag     = ',sysflag
        write(un_sys,*)'delta       = ',delta  
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conf    = ',error
        ! other physcial parameters
        write(un_sys,*)'T           = ',Tref
        write(un_sys,*)'VdWepsC     = ',VdWepsC*vpolC*vsol
        write(un_sys,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
        ! volume 
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
        write(un_sys,*)'===end distance independent settings=='
    endif
    
    write(un_sys,*)'distance    = ',nz*delta 
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'free energy = ',FE  
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'FEconfC     = ',FEconfC
    write(un_sys,*)'FEconfAB    = ',FEconfAB
    write(un_sys,*)'FEtrans%sol = ',FEtrans%sol  
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
    write(un_sys,*)'sumphiC     = ',sumphiC
    write(un_sys,*)'check phi   = ',checkphi 
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'qAB         = ',qAB
    write(un_sys,*)'qC          = ',qC
    write(un_sys,*)'muAB        = ',-log(qAB)
    write(un_sys,*)'muC         = ',-log(qC)
    write(un_sys,*)'heightAB    = ',heightAB
    write(un_sys,*)'heightC     = ',heightC
    write(un_sys,*)'nsize       = ',nsize  
    write(un_sys,*)'cuantasAB   = ',cuantasAB
    write(un_sys,*)'cuantasC    = ',cuantasC
    write(un_sys,*)'iterations  = ',iter
    
    ! .. closing files
    if(nz.eq.nzmin) then 
        close(un_xsol)
        close(un_xpolAB)
        close(un_xpolC)
        close(un_sys)
    endif
  
end subroutine output_neutral


subroutine output_individualcontr_fe

    use globals, only : LEFT,RIGHT, sysflag
    use energy
    use myutils, only : newunit
    use parameters, only : sigmaAB,sigmaABL,sigmaABR,cNaCl,cCaCl2,pHbulk,VdWepsB
    use volume, only : delta,nz,nzmax,nzmin

    implicit none 

    ! local arguments

    character(len=100) :: fenergyfilename   
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    if(nz==nzmax) then 

        !     .. make label filename

        if(sysflag=="electdouble") then 

            write(rstr,'(F5.3)')sigmaABL
            fnamelabel="sgL"//trim(adjustl(rstr))
            write(rstr,'(F5.3)')sigmaABR
            fnamelabel=trim(fnamelabel)//"sgR"//trim(adjustl(rstr))  
            write(rstr,'(F5.3)')cNaCl
            fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
            if(cCaCl2>=0.001) then 
                write(rstr,'(F5.3)')cCaCl2
            elseif(cCaCl2>0.0) then  
                write(rstr,'(ES9.2E2)')cCaCl2
            else 
                write(rstr,'(F3.1)')cCaCl2
            endif    
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
            write(rstr,'(F7.3)')pHbulk
            fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"


        elseif(sysflag=="elect".or.sysflag=="electnopoly".or.sysflag=="electA") then 


            write(rstr,'(F5.3)')sigmaABL
            fnamelabel="sg"//trim(adjustl(rstr)) 
            write(rstr,'(F5.3)')cNaCl
            fnamelabel=trim(fnamelabel)//"cNaCl"//trim(adjustl(rstr))
            if(cCaCl2>=0.001) then 
                write(rstr,'(F5.3)')cCaCl2
            elseif(cCaCl2>0.0) then  
                write(rstr,'(ES9.2E2)')cCaCl2
            else 
                write(rstr,'(F3.1)')cCaCl2
            endif 
            fnamelabel=trim(fnamelabel)//"cCaCl2"//trim(adjustl(rstr))
            write(rstr,'(F7.3)')pHbulk
            fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))//".dat"

        elseif(sysflag=="neutral") then 
            
            write(rstr,'(F5.3)')sigmaAB*delta 
            fnamelabel="sg"//trim(adjustl(rstr)) 
            write(rstr,'(F5.3)')VdWepsB
            fnamelabel=trim(fnamelabel)//"VdWepsB"//trim(adjustl(rstr))//".dat"
        else
            
            print*,"Error in output_individualcontr_fe subroutine"
            print*,"Wrong value sysflag : ", sysflag
        endif    

        fenergyfilename='energy.'//trim(fnamelabel)   
            
        !     .. opening files        

        open(unit=newunit(un_fe),file=fenergyfilename) 

    endif      
    write(un_fe,*)'#D              = ',nz*delta 
    write(un_fe,*)'FE              = ',FE  
    write(un_fe,*)'FEbulk          = ',FEbulk 
    write(un_fe,*)'deltaFE         = ',deltaFE
    write(un_fe,*)'FEalt           = ',FEalt  
    write(un_fe,*)'FEbulkalt       = ',FEbulkalt 
    write(un_fe,*)'deltaFEalt      = ',deltaFEalt
    write(un_fe,*)'deltadeltaFE    = ',deltaFE-deltaFEalt
    
    write(un_fe,*)'FEq             = ',FEq  
    write(un_fe,*)'FEpi            = ',FEpi
    write(un_fe,*)'FErho           = ',FErho
    write(un_fe,*)'FEel            = ',FEel
    write(un_fe,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(un_fe,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT) 
    write(un_fe,*)'FEbind          = ',FEbind 
    write(un_fe,*)'FEchem          = ',FEchem
    write(un_fe,*)'FEconfAB        = ',FEconfAB
    write(un_fe,*)'FEconfC         = ',FEconfC
    
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

    write(un_fe,*)"FEchemsurf(LEFT)     = ",FEchemsurf(LEFT)
    write(un_fe,*)"FEchemsurf(RIGHT)    = ",FEchemsurf(RIGHT)
    write(un_fe,*)"FEchemsurfalt(LEFT)  = ",FEchemsurfalt(LEFT)
    write(un_fe,*)"FEchemsurfalt(RIGHT) = ",FEchemsurfalt(RIGHT)

    write(un_fe,*)"delta FEchemsurfalt(LEFT) = ",FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
    write(un_fe,*)"delta FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)


    if(nz==nzmin) close(un_fe)

end subroutine output_individualcontr_fe

     
subroutine copy_solution(x)

    use globals, only : sysflag, neq, nsize, bcflag, LEFT, RIGHT
    use volume, only  : nx,ny
    use surface, only : psiSurfL, psiSurfR
    use field
 
    implicit none

    real(dp), dimension(neq) :: x  ! expliciet size array 
    
    ! local variable
    integer :: i, neq_bc

    select case (sysflag)
    case ("elect")   

        do i=1,nsize                   
            xsol(i)= x(i)              
            psi(i) = x(i+nsize)        
            rhopolA(i)=x(i+2*nsize)
            rhopolB(i)=x(i+3*nsize)
        enddo
        
        neq_bc=0 ! surface potential  
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(4*nsize+i) 
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(4*nsize+neq_bc+i)         
            enddo
            neq_bc=neq_bc+nx*ny
        endif    
       
    case ("electdouble")

        do i=1,nsize                     
            xsol(i)= x(i)            
            psi(i) = x(i+nsize)         
            rhopolA(i)=x(i+2*nsize)
            rhopolB(i)=x(i+3*nsize)
        enddo

    case ("electnopoly")
    
        do i=1,nsize                   
            xsol(i)= x(i)          
            psi(i) = x(i+nsize)     
        enddo
            
        neq_bc=0
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(2*nsize+i)                 
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(2*nsize+neq_bc+i)          
            enddo
            neq_bc=neq_bc+nx*ny
        endif    

    case ("neutral")

        do i=1,nsize                   
            xsol(i)= x(i)             
            rhopolB(i)=x(i+nsize)       
        enddo
        print*,"copy_solution: case neutral not correct"
        print*,"stopping program"
        stop

    case ("dipolarstrong")
    
        do i=1,nsize                   
            xsol(i)= x(i)          
            psi(i) = x(i+nsize)     
        enddo
            
        neq_bc=0
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(2*nsize+i)                 
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(2*nsize+neq_bc+i)          
            enddo
            neq_bc=neq_bc+nx*ny
        endif    
    
    case ("dipolarweak")
    
        do i=1,nsize                   
            xsol(i)= x(i)          
            psi(i) = x(i+nsize)     
        enddo
            
        neq_bc=0
        if(bcflag(RIGHT)/="cc") then
            neq_bc=nx*ny
            do i=1,neq_bc
                psiSurfR(i) =x(2*nsize+i)                 
            enddo
        endif   
        if(bcflag(LEFT)/="cc") then 
            do i=1,nx*ny
                psiSurfL(i) =x(2*nsize+neq_bc+i)          
            enddo
            neq_bc=neq_bc+nx*ny
        endif 

    case default   

        print*,"Error: sysflag incorrect at copy_solution"
        print*,"stopping program"
        stop
    
    end select
         
end subroutine copy_solution

! output routine 

subroutine compute_vars_and_output()

    use globals, only : sysflag
    use energy
    use field
    use parameters, only : heightAB
    implicit none

    select case (sysflag)
    case ("elect")
    
        call fcnenergy()       
        call average_density_z(xpolAB,xpolABz,heightAB)  
        call output()

    case ("electA")
    
        call fcnenergy()       
        call average_density_z(xpolAB,xpolABz,heightAB)  
        call output()

    case ("electdouble")
    
        call fcnenergy()       
        call average_density_z(xpolAB,xpolABz,heightAB)  
        call output()
    
    case ("electnopoly")

        call fcnenergy()        
        call output()           ! writing of output
     
    case ("neutral")  

        call output()           ! writing of output
    
    case ("dipolarweak")

        call average_density_z(xpolAB,xpolABz,heightAB)  
        call output()           ! writing of output 

    case ("dipolarstrong")
        call average_density_z(xpolAB,xpolABz,heightAB)  
        call output()           ! writing of output 

    case default   
        print*,"Error: sysflag incorrect at compute_vars_and_output"
        print*,"stopping program"
        stop
    end select
         
end subroutine compute_vars_and_output



end module


