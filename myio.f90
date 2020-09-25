module myio

    use precision_definition
    implicit none

    ! return error values

    integer, parameter ::  myio_err_systype   = 1
    integer, parameter ::  myio_err_runtype   = 2
    integer, parameter ::  myio_err_geometry  = 3
    integer, parameter ::  myio_err_method    = 4
    integer, parameter ::  myio_err_chaintype = 5
    integer, parameter ::  myio_err_domain    = 6
    integer, parameter ::  myio_err_inputfile = 7
    integer, parameter ::  myio_err_input     = 8
    integer, parameter ::  myio_err_bcflag    = 9 
    integer, parameter ::  myio_err_label     = 10
    integer, parameter ::  myio_err_VdWeps    = 11 
    integer, parameter ::  myio_err_nsegtypes = 12
    integer, parameter ::  myio_err_chainmethod = 13  
    integer, parameter ::  myio_err_chainsfile = 14 
    integer, parameter ::  myio_err_energyfile = 15
    integer, parameter ::  myio_err_dielect   = 16
    integer, parameter ::  myio_err_graft     = 17
    integer, parameter ::  myio_err_index     = 18
    integer, parameter ::  myio_err_conf      = 19
    integer, parameter ::  myio_err_nseg      = 20


    ! unit number 
    integer :: un_sys,un_xpolAB,un_xsol,un_xNa,un_xCl,un_xK,un_xCa,un_xNaCl,un_xKCl
    integer :: un_xOHmin,un_xHplus,un_fdisA,un_fdisB,un_psi,un_charge, un_xpair, un_rhopolAB, un_fe, un_q
    integer :: un_dip ,un_dielec,un_xpolABz, un_xpolz, un_xpol, un_fdis, un_xpro, un_fdisP
   
    ! format specifiers 
    character(len=80), parameter  :: fmt = "(A9,I1,A5,ES25.16)"
    character(len=80), parameter  :: fmt1reals = "(ES25.16E3)"   
    character(len=80), parameter  :: fmt2reals = "(2ES25.16E3)"   
    character(len=80), parameter  :: fmt3reals = "(3ES25.16E3)"  
    character(len=80), parameter  :: fmt4reals = "(4ES25.16E3)" 
    character(len=80), parameter  :: fmt5reals = "(5ES25.16E3)"
    character(len=80), parameter  :: fmt6reals = "(6ES25.16E3)" 
    
    private
    public :: read_inputfile, output_individualcontr_fe, output, compute_vars_and_output,write_chain_config
    public :: myio_err_chainsfile, myio_err_energyfile, myio_err_chainmethod, myio_err_geometry
    public :: myio_err_graft, myio_err_index, myio_err_conf, myio_err_nseg
    
contains

subroutine read_inputfile(info)

    use globals
    use parameters
    use surface 
    use myutils, only : newunit

    integer, intent(out), optional :: info

    ! .. local arguments

    integer :: info_sys, info_bc, info_run, info_geo, info_meth, info_chaintype, info_combi, info_VdWeps
    integer :: info_chainmethod, info_dielect
    character(len=8) :: fname
    integer :: ios,un_input  ! un = unit number    
    character(len=100) :: buffer, label
    integer :: pos
    integer :: line
    logical :: isSet_maxnchains, isSet_maxnchainsxy, isSet_precondition, isSet_savePalpha

    if (present(info)) info = 0
    
    !     .. reading in of variables from file
    write(fname,'(A8)')'input.in'
    open(unit=newunit(un_input),file=fname,iostat=ios,status='old')
    if(ios >0 ) then
        print*, 'Error opening input.in file : iostat =', ios
        if (present(info)) info = myio_err_inputfile
        return
    endif

    ! defaults
    isSet_maxnchains  =.false.
    isSet_maxnchainsxy=.false.  
    isSet_precondition=.false.
    isSet_savePalpha  =.false.  
    write_mc_chains   =.false.

    ! defailt concentrations
    cKCl=0.0_dp
    cCaCl2=0.0_dp
    cMgCl2=0.0_dp
    cRbCl=0.0_dp
    

    ios=0 
    line = 0

    ! ios<0 : if an end of record condition is encountered or if an end of file condition was detected.  
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
            case ('systype')
                read(buffer, *,iostat=ios) systype
            case ('runtype')
                read(buffer, *,iostat=ios) runtype
            case ('bcflag(LEFT)')
                read(buffer,*,iostat=ios) bcflag(LEFT)
            case ('bcflag(RIGHT)')
                read(buffer,*,iostat=ios) bcflag(RIGHT)
            case ('chainmethod')
                read(buffer,*,iostat=ios) chainmethod
            case ('chaintype')
                read(buffer,*,iostat=ios) chaintype
            case ('tolerance')
                read(buffer,*,iostat=ios) tol_conv           
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
            case ('cNaCl')
                read(buffer,*,iostat=ios) cNaCl
            case ('cKCl')
                read(buffer,*,iostat=ios) cKCl
            case ('cRbCl')
                read(buffer,*,iostat=ios) cRbCl
            case ('cCaCl2')
                read(buffer,*,iostat=ios) cCaCl2
            case ('cMgCl2')
                read(buffer,*,iostat=ios) cMgCl2
            case ('cpro%val')
                read(buffer,*,iostat=ios) cpro%val
            case ('cpro%min')
                read(buffer,*,iostat=ios) cpro%min
            case ('cpro%max')
                read(buffer,*,iostat=ios) cpro%max
            case ('cpro%stepsize')
                read(buffer,*,iostat=ios) cpro%stepsize
            case ('cpro%delta')
                read(buffer,*,iostat=ios) cpro%delta
            case ('Rpro')
                read(buffer,*,iostat=ios) Rpro
            case ('nsize')
                read(buffer,*,iostat=ios) nsize
            case ('nseg')
                read(buffer,*,iostat=ios) nseg
            case ('nsegtypes')    
                read(buffer,*,iostat=ios) nsegtypes        ! carefully need to be overwriiten depending on value systype and or chainmethod    
            case ('maxnchains')    
                read(buffer,*,iostat=ios) maxnchainsrotations  
                isSet_maxnchains=.true.      
            case ('maxnchainsxy')    
                read(buffer,*,iostat=ios) maxnchainsrotationsXY
                isSet_maxnchainsxy=.true.        
            case ('cuantas')
                read(buffer,*,iostat=ios) max_confor
            case ('chainperiod')
                read(buffer,*,iostat=ios) chainperiod
            case('vpolfname')
                read(buffer,*,iostat=ios) vpolfname
            case('pKafname')
                read(buffer,*,iostat=ios) pKafname      
            case('typesfname') 
                read(buffer,*,iostat=ios) typesfname  
            case('lsegfname') 
                read(buffer,*,iostat=ios) lsegfname  
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
            case('unit_conv')
                read(buffer,*,iostat=ios) unit_conv
            case ('geometry')
                read(buffer,*,iostat=ios) geometry
            case ('ngr_freq')
                read(buffer,*,iostat=ios) ngr_freq 
            case ('sgraft')
                read(buffer,*,iostat=ios) sgraft
            case ('nset_per_graft')
                read(buffer,*,iostat=ios) nset_per_graft      
            case ('gamma')
                read(buffer,*,iostat=ios) gamma 
            case ('isRandom_pos_graft')
                read(buffer,*,iostat=ios) isRandom_pos_graft  
            case ('seed_graft')
                read(buffer,*,iostat=ios) seed_graft  
            case ('write_mc_chains')
                read(buffer,*,iostat=ios) write_mc_chains        
            case ('precondition')    
                read(buffer,*,iostat=ios) precondition
                isSet_precondition=.true.   
            case ('dielect_env')
                read(buffer,*,iostat=ios) dielect_env
            case ('VdWscale%val')
                read(buffer,*,iostat=ios) VdWscale%val
            case ('VdWscale%min')
                read(buffer,*,iostat=ios) VdWscale%min
            case ('VdWscale%max')
                read(buffer,*,iostat=ios) VdWscale%max
            case ('VdWscale%stepsize')
                read(buffer,*,iostat=ios) VdWscale%stepsize
            case ('VdWscale%delta')
                read(buffer,*,iostat=ios) VdWscale%delta     
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
    
    bcflag(LEFT)="cc"
    bcflag(RIGHT)="cc"  
    

    ! .. check error flag

    call check_value_systype(systype,info_sys) 
    if (info_sys == myio_err_systype) then
        if (present(info)) info = info_sys
        return
    endif

    call check_value_runtype(runtype,info_run) 
    if (info_sys == myio_err_runtype) then
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

    call check_value_chainmethod(chainmethod,info_chainmethod)
    if (info_chainmethod == myio_err_chainmethod) then            
        if (present(info)) info = info_chainmethod
        return
    endif

    call check_value_dielect_env(dielect_env,info_dielect)
    if (info_dielect == myio_err_dielect) then
        if (present(info)) info = info_dielect
        return
    endif

    !  .. set input values

    call set_value_nzmin(runtype,nzmin,nzmax)
    call set_value_isVdW(systype,isVdW)
    call set_value_isVdWintEne(systype, isVdWintEne)
    call set_value_nsegtypes(nsegtypes,chaintype,systype,info)
    call set_value_maxnchains(maxnchainsrotations,isSet_maxnchains)
    call set_value_maxnchainsxy(maxnchainsrotationsxy,isSet_maxnchainsxy)
    call set_value_precondition(precondition,isSet_precondition)


    ! after set_value_isVdW
    call check_value_VdWeps(systype,isVdW,info_VdWeps)
    if (info_VdWeps == myio_err_VdWeps) then
        if (present(info)) info = info_VdWeps
        return
    endif

    ! overide certain input values
    if(systype=="brushssdna".or.systype=="brushborn".or.systype=="brush_mul") then
        KionNa=0.0_dp 
        KionK=0.0_dp
        cpro%val =0.0_dp
        Rpro=0.1_dp
    endif    

end subroutine read_inputfile
 

subroutine check_value_systype(systype,info)

    character(len=15), intent(in) :: systype
    integer, intent(out),optional :: info

    character(len=15) :: systypestr(8)
    integer :: i
    logical :: flag

    ! permissible values of systype

    systypestr(1)="elect"
    systypestr(2)="neutral"
    systypestr(3)="brush_mul" 
    systypestr(4)="brush_mulnoVdW"
    systypestr(5)="brushssdna"
    systypestr(6)="brushborn"
    systypestr(7)="bulk water"
    systypestr(8)="neutralnoVdW"
   
    flag=.FALSE.

    do i=1,8
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


subroutine check_value_runtype(runtype,info)

    character(len=15), intent(in) :: runtype
    integer, intent(out),optional :: info

    character(len=15) :: runtypestr(4)
    integer :: i
    logical :: flag

    ! permissible values of runtype

    runtypestr(1)="rangepH"
    runtypestr(2)="rangeVdWeps"
    runtypestr(3)="rangedist"
    runtypestr(4)="rangecpro"

    flag=.FALSE.

    do i=1,4
        if(runtype==runtypestr(i)) flag=.TRUE.
    enddo

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error: value of runtype is not permissible"
        print*,"runtype = ",runtype
        if (present(info)) info = myio_err_runtype
        return
    end if

end subroutine check_value_runtype

subroutine check_value_bcflag(bcflag,info)

    use globals, only : LEFT, RIGHT

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

subroutine check_value_geometry(geometry,info)
        
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

    character(len=8), intent(in) :: chaintype
    integer, intent(out),optional :: info

    logical :: flag
    character(len=8) :: chaintypestr(6)
    integer :: i

    ! permissible values of chaintype

    chaintypestr(1)="diblockA"
    chaintypestr(2)="diblockB"
    chaintypestr(3)="altA"
    chaintypestr(4)="altB"
    chaintypestr(5)="copolyAB"
    chaintypestr(6)="multi"

    flag=.FALSE.

    do i=1,6
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

subroutine check_value_chainmethod(chainmethod,info)
      
    character(len=15), intent(in) :: chainmethod
    integer, intent(out),optional :: info

    logical :: flag
    character(len=15) :: chainmethodstr(3) 
    integer :: i

    ! permissible values of chainmethod

    chainmethodstr(1)="MC"
    chainmethodstr(2)="FILE_lammps_xyz"
    chainmethodstr(3)="FILE_lammps_trj"
       
    flag=.FALSE.

    do i=1,3
        if(chainmethod==chainmethodstr(i)) flag=.TRUE.
    enddo
        
    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then 
        print*,"Error: value of chainmethod is not permissible"
        print*,"chainmethod = ",chainmethod
        if (present(info)) info = myio_err_chainmethod
        return
    endif

end subroutine check_value_chainmethod


subroutine check_value_method(method,info)

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

subroutine check_value_dielect_env(dielect_env,info)
        
    character(len=15), intent(in) :: dielect_env
    integer, intent(out),optional :: info

    logical :: flag
    character(len=15) :: dielect_env_str(3) 
    integer :: i

    ! permissible values of dielect_env

    dielect_env_str(1)="constant"
    dielect_env_str(2)="linear"
    dielect_env_str(3)="MaxwellGarnett"
    
    flag=.FALSE.

    do i=1,3
        if(dielect_env==dielect_env_str(i)) flag=.TRUE.
    enddo
        
    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then 
        print*,"Error: value of dielect_env is not permissible"
        print*,"dielect_env = ",dielect_env
        if (present(info)) info = myio_err_dielect 
        return
    endif
    
end subroutine check_value_dielect_env



subroutine check_value_VdWeps(systype,isVdW,info)

    logical, intent(in) :: isVdW
    character(len=15), intent(in) :: systype
    integer, intent(out), optional :: info

    character(len=15) :: systypestr(4)
    integer :: i
    logical :: flag

    flag=.false.

    if(isVdW) then ! check for correct combination systype 
     
        systypestr(1)="neutral"
        systypestr(2)="brush_mul"
        systypestr(3)="brushssdna"
        systypestr(4)="brushborn"
       
        do i=1,4! sofar only electA works with VdW 
            if(systype==systypestr(i)) flag=.TRUE.
        enddo
    else  ! isVdW=.false. so oke 
        flag=.true.
    endif        

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then
        print*,"Error:  combination systype and isVdW is not permissible"
        print*,"systype = ",systype, " isVdW =",isVdW
        if (present(info)) info = myio_err_VdWeps
        return
    end if

end subroutine check_value_VdWeps



! override input value nzmin
! Sets nzmin=nzmax for which runtype that do not loop over nz 
! ensuring that the output files are properly close once

subroutine set_value_nzmin(runtype,nzmin,nzmax)

    character(len=15), intent(in) :: runtype
    integer, intent(inout) :: nzmin
    integer, intent(in) :: nzmax
   
    if (runtype/="rangedist") nzmin=nzmax 

end subroutine

! isVdW = .true. then uses position dependent Van der Waals interaction

subroutine set_value_isVdW(systype, isVdW)

    character(len=15), intent(in) :: systype
    logical, intent(inout)  :: isVdW

    character(len=15) :: systypestr(3) 
    integer :: i

    ! permissible values of systype that involve VdW interaction 
    
    systypestr(1)="elect"
    systypestr(2)="neutralnoVdW" 
    systypestr(3)="brush_mulnoVdW"  
    isVdW=.True.

    do i=1,3     
        if(systype==systypestr(i)) isVdW=.FALSE.
    enddo    

end subroutine

! isVdWintEne = .true. then  compute internal VdW energy chain 

subroutine set_value_isVdWintEne(systype, isVdWintEne)

    character(len=15), intent(in) :: systype
    logical, intent(inout)  :: isVdWintEne

    character(len=15) :: systypestr(2) 
    integer :: i

    ! permissible values of systype that involve internal VdW chain energy 
    
    systypestr(1)="neutralnoVdW"   
    systypestr(2)="brush_mulnoVdW"  

    isVdWintEne=.False.
   
    do i=1,2
        if(systype==systypestr(i)) isVdWintEne=.true.
    enddo

end subroutine

! changes value isVdW based on values of VdWeps parameter
! pre :VdWeps need to be initialized see module VdW

subroutine set_value_isVdW_on_values(nsegtypes, VdWeps, isVdW)

    integer, intent(in) :: nsegtypes
    real(dp), intent(in) :: VdWeps(:,:)
    logical, intent(inout)  :: isVdW

    integer :: s,t
  
    do s=1,nsegtypes
        do t=1,nsegtypes   
            if (abs(VdWeps(s,t))>1.0e-4_dp) isVdW=.true.
        enddo
    enddo    

end subroutine


subroutine set_value_nsegtypes(nsegtypes,chaintype,systype,info)
            

    integer, intent(inout) :: nsegtypes
    integer, intent(out),optional :: info
    character(len=8), intent(in) :: chaintype
    character(len=15), intent(in) :: systype

    logical :: flag
    character(len=8) :: chaintypestr(5) 
    integer :: i

    ! permissible values of chaintype

    chaintypestr(1)="diblockA"
    chaintypestr(2)="diblockB"
    chaintypestr(3)="altA"
    chaintypestr(4)="altB"
    chaintypestr(5)="copolyAB"
    

    ! two components
    flag=.FALSE.
    
    do i=1,5
        if(chaintype==chaintypestr(i)) flag=.TRUE.
    enddo
    if(flag) nsegtypes=2

    if(chaintype=="multi") flag=.true.

    if(systype=="neutral") flag=.true.    

    if (present(info)) info = 0

    if (flag.eqv. .FALSE.) then 
        print*,"Error: value of nsegtype could not be set for given "
        print*,"chaintype = ",chaintype
        if (present(info)) info = myio_err_nsegtypes
        return
    endif
    
end subroutine set_value_nsegtypes

    
subroutine set_value_maxnchains(maxnchainsrotations,isSet_maxnchains)

    integer, intent(inout) :: maxnchainsrotations
    logical, intent(in)  :: isSet_maxnchains

    if(.not.isSet_maxnchains) then 
        maxnchainsrotations=12 ! default value
    else
        maxnchainsrotations=abs(maxnchainsrotations) !  make positive   
    endif
        
end subroutine set_value_maxnchains


    
subroutine set_value_maxnchainsxy(maxnchainsrotationsxy,isSet_maxnchainsxy)

    integer, intent(inout) :: maxnchainsrotationsxy
    logical, intent(in)  :: isSet_maxnchainsxy

    if(.not.isSet_maxnchainsxy) then 
        maxnchainsrotationsxy=1 ! default value
    else
        maxnchainsrotationsxy=abs(maxnchainsrotationsxy) !  make positive   
    endif
        
end subroutine set_value_maxnchainsxy


subroutine set_value_precondition(precondition,isSet_precondition)

    logical, intent(inout) :: precondition
    logical, intent(in)  :: isSet_precondition

    if(.not.isSet_precondition) precondition=.false. ! default value
   
end subroutine set_value_precondition

subroutine output()

    use globals, only : systype

    select case (systype)
    case ("elect")
        call output_elect  
    case("neutral","neutralnoVdW") 
        call output_neutral
        call output_individualcontr_fe
    case("brush_mul","brush_mulnoVdW","brushssdna") 
        call output_brush_mul  
        call output_individualcontr_fe  
    case default
        print*,"Error in output subroutine"
        print*,"Wrong value systype : ", systype
    end select     

end subroutine output


subroutine output_brush_mul
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    use myutils, only : newunit
    use chains, only : isHomopolymer
    
    !     .. local arguments
   
    integer :: t

    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolfilename 
    character(len=90) :: xpolzfilename 
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
    character(len=90) :: densfracfilename
    character(len=90) :: densfracPfilename
    character(len=90) :: qfilename
    character(len=90) :: densfracionpairfilename
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    integer :: i,j,k          ! dummy indexes
    real(dp) :: denspol

    ! .. executable statements 

    denspol=init_denspol()

    if(nz.eq.nzmax)  then 
        !     .. make label filenames f
        call make_filename_label(fnamelabel)

        sysfilename='system.'//trim(fnamelabel)
        xpolfilename='xpol.'//trim(fnamelabel)
        xpolzfilename='xpolz.'//trim(fnamelabel)
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
        densfracfilename='densityfrac.'//trim(fnamelabel)
        densfracPfilename='densityfracP.'//trim(fnamelabel)
        densfracionpairfilename='densityfracionpair.'//trim(fnamelabel)
        qfilename='q.'//trim(fnamelabel)
       
        !     .. opening files        
        
        open(unit=newunit(un_sys),file=sysfilename)       
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_psi),file=potentialfilename)
         
        open(unit=newunit(un_xpol),file=xpolfilename)
        open(unit=newunit(un_fdis),file=densfracfilename)  
        if(systype=="brushssdna") open(unit=newunit(un_fdisP),file=densfracPfilename) 
        open(unit=newunit(un_q),file=qfilename)        
        open(unit=newunit(un_xpolz),file=xpolzfilename)
        
       
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
        inquire(unit=un_xpol, opened=isopen)
        inquire(unit=un_xsol, opened=isopen)

        if(.not.isopen) write(*,*)"un_xsol is not open"

    endif       

    !   .. writting files
    !   .. this line seperates different distances 
      
    write(un_xsol,*)'#D    = ',nz*delta 
    write(un_psi,*)'#D    = ',nz*delta
    write(un_xpol,*)'#D    = ',nz*delta 
    write(un_fdis,*)'#D    = ',nz*delta
    if(systype=="brushssdna") write(un_fdisP,*)'#D    = ',nz*delta
    
     
    if(verboseflag=="yes") then    
        write(un_xNa,*)'#D    = ',nz*delta 
        write(un_xK,*)'#D    = ',nz*delta 
        write(un_xCa,*)'#D    = ',nz*delta 
        write(un_xNaCl,*)'#D    = ',nz*delta 
        write(un_xKCl,*)'#D    = ',nz*delta 
        write(un_xpair,*)'#D    = ',nz*delta 
        write(un_charge,*)'#D    = ',nz*delta 
        write(un_xCl,*)'#D    = ',nz*delta 
        write(un_xHplus,*)'#D    = ',nz*delta 
        write(un_xOHMin,*)'#D    = ',nz*delta 
    endif

!    do i=1,nsurf
!        write(un_psi,*)psiSurfL(i)
!    enddo 
        
    do i=1,nsize
        write(un_xsol,*)xsol(i)
        write(un_psi,*)psi(i)
    enddo    

!~    do i=1,nsurf
!       write(un_psi,*)psiSurfR(i)
!    enddo         

    write(un_q,*)q

    do i=1,nsize
        write(un_xpol,*)xpol(i),(rhopol(i,t),t=1,nsegtypes)
        write(un_fdis,*)(fdis(i,t),t=1,nsegtypes)
    enddo

    do i=1,nz
        write(un_xpolz,fmt1reals)xpolz(i)
    enddo     
    
    if(systype=="brushssdna")then
        do i=1,nsize
            write(un_fdisP,*)(fdisA(i,k),k=1,7)
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

    ! .. writing system information 
    if(nz.eq.nzmax) then 
        write(un_sys,*)'===begin distance independent settings=='  
        write(un_sys,*)'system      = planar weakpolyelectrolyte brush'
        write(un_sys,*)'version     = ',VERSION
        ! chain description 
        write(un_sys,*)'chainmethod = ',chainmethod
        write(un_sys,*)'chaintype   = ',chaintype
        write(un_sys,*)"isHomopolymer= ",isHomopolymer
        write(un_sys,*)'nseg        = ',nseg
        write(un_sys,*)'nsegtypes   = ',nsegtypes        
        write(un_sys,*)'lseg        = ',(lsegAA(t),t=1,nsegtypes)
        write(un_sys,*)'cuantas     = ',cuantas
        write(un_sys,*)'denspol     = ',denspol
       
        ! system description
        write(un_sys,*)'systype     = ',systype
        write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
        write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
        write(un_sys,*)'delta       = ',delta
        write(un_sys,*)'ngr         = ',ngr  
        write(un_sys,*)'ngr_node    = ',ngr_node  
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conv    = ',tol_conv
        ! concentration 
        write(un_sys,*)'cNaCl       = ',cNaCl
        write(un_sys,*)'cKCl        = ',cKCl
        write(un_sys,*)'cCaCl2      = ',cCaCl2
        write(un_sys,*)'xsolbulk    = ',xbulk%sol
        write(un_sys,*)'xNabulk     = ',xbulk%Na
        write(un_sys,*)'xClbulk     = ',xbulk%Cl
        write(un_sys,*)'xKbulk      = ',xbulk%K
        write(un_sys,*)'xRbbulk     = ',xbulk%Rb
        write(un_sys,*)'xNaClbulk   = ',xbulk%NaCl
        write(un_sys,*)'xKClbulk    = ',xbulk%KCl
        write(un_sys,*)'xCabulk     = ',xbulk%Ca
        write(un_sys,*)'xMgbulk     = ',xbulk%Mg
        write(un_sys,*)'xHplusbulk  = ',xbulk%Hplus
        write(un_sys,*)'xOHminbulk  = ',xbulk%OHmin
        write(un_sys,*)'pHbulk      = ',pHbulk

        ! disociation constants 
        write(un_sys,*)'pKa         = ',(pKa(t),t=1,nsegtypes)      
      
        write(un_sys,*)'KionNa      = ',KionNa
        write(un_sys,*)'KionK       = ',KionK
        write(un_sys,*)'K0ionNa     = ',K0ionNa
        write(un_sys,*)'K0ionK      = ',K0ionK
        write(un_sys,*)'dielectW    = ',dielectW
        write(un_sys,*)'lb          = ',lb
        write(un_sys,*)'T           = ',Tref

        ! charge components
        do t=1,nsegtypes
            write(un_sys,*)'zpol(',t,')     = ',zpol(t,1),zpol(t,2) 
        enddo
        write(un_sys,*)'zNa         = ',zNa
        write(un_sys,*)'zCa         = ',zCa
        write(un_sys,*)'zK          = ',zK
        write(un_sys,*)'zCl         = ',zCl 
         ! volume  
        write(un_sys,*)'vsol        = ',vsol
        write(un_sys,*)'vpol        = ',(vpol(t)*vsol,t=1,nsegtypes)
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
    write(un_sys,*)'tol_conv    = ',tol_conv
    write(un_sys,*)'denspol     = ',denspol
    write(un_sys,*)'sumphi      = ',(sumphi(t),t=1,nsegtypes)
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
    write(un_sys,*)'height      = ',height
    write(un_sys,*)'qpol        = ',(qpol(t),t=1,nsegtypes)
    write(un_sys,*)'qpoltot     = ',qpol_tot
    write(un_sys,*)'avfdis      = ',(avfdis(t),t=1,nsegtypes)
    if(systype=="brushssdna")then
        write(un_sys,*)'avfdisA      = ',(avfdisA(k),k=1,7)
    endif    
    write(un_sys,*)'sigmaSurfL  = ',sigmaSurfL/((4.0_dp*pi*lb)*delta)
    write(un_sys,*)'sigmaSurfR  = ',sigmaSurfR/((4.0_dp*pi*lb)*delta)

   
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
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'iterations  = ',iter
   
    ! .. closing files

    if(nz==nzmin) then 
        close(un_sys)
        close(un_xsol)
        close(un_psi)
        close(un_xpol)
        close(un_fdis)
        if(systype=="brushssdna") close(un_fdisP)
        close(un_xpolz)            
        close(un_q)
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

end subroutine output_brush_mul


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
    
    !     .. local arguments
   
    !     .. output file names       
    
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolfilename 
    character(len=90) :: xpolzfilename 
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
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    integer :: i,j,k          ! dummy indexes
    real(dp) :: denspol

    ! .. executable statements 

    denspol=init_denspol()

    if(nz.eq.nzmax)  then 
        
        !     .. make label filename
        call make_filename_label(fnamelabel)

        fnamelabel=trim(fnamelabel)//".dat"

        sysfilename='system.'//trim(fnamelabel)
        xpolfilename='xpol.'//trim(fnamelabel)
        xpolzfilename='xpolz.'//trim(fnamelabel)
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
       
        !     .. opening files        
        
        open(unit=newunit(un_sys),file=sysfilename)       
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_psi),file=potentialfilename)

                 
        open(unit=newunit(un_xpol),file=xpolfilename)
        open(unit=newunit(un_fdisA),file=densfracAfilename) 
        open(unit=newunit(un_fdisB),file=densfracBfilename)             
        open(unit=newunit(un_q),file=qfilename)
        open(unit=newunit(un_xpolz),file=xpolzfilename)
        
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
        inquire(unit=un_xpol, opened=isopen)
        inquire(unit=un_xsol, opened=isopen)

        if(.not.isopen) write(*,*)"un_xsol is not open"

    endif       

    !   .. writting files
    !   .. this line seperates different distances 
      
    write(un_xsol,*)'#D    = ',nz*delta 
    write(un_psi,*)'#D    = ',nz*delta       
    write(un_xpolAB,*)'#D    = ',nz*delta 
    write(un_fdisA,*)'#D    = ',nz*delta
    write(un_fdisB,*)'#D    = ',nz*delta
    
    if(verboseflag=="yes") then    
        write(un_xNa,*)'#D    = ',nz*delta 
        write(un_xK,*)'#D    = ',nz*delta 
        write(un_xCa,*)'#D    = ',nz*delta 
        write(un_xNaCl,*)'#D    = ',nz*delta 
        write(un_xKCl,*)'#D    = ',nz*delta 
        write(un_xpair,*)'#D    = ',nz*delta 
        write(un_charge,*)'#D    = ',nz*delta 
        write(un_xCl,*)'#D    = ',nz*delta 
        write(un_xHplus,*)'#D    = ',nz*delta 
        write(un_xOHMin,*)'#D    = ',nz*delta 
    endif

    !do i=1,nsurf
    !    write(un_psi,*)psiSurfL(i)
    !enddo 
        
    do i=1,nsize
        write(un_xsol,*)xsol(i)
        write(un_psi,*)psi(i)
    enddo    

    !do i=1,nsurf
    !    write(un_psi,*)psiSurfR(i)
    !enddo         

    write(un_q,*)q    

    do i=1,nsize
        write(un_xpol,fmt3reals)xpol(i),rhopol(i,1),rhopol(i,2)
        write(un_fdisA,fmt5reals)(fdisA(i,k),k=1,5)        
        write(un_fdisB,fmt5reals)(fdisB(i,k),k=1,5)
    enddo
    do i=1,nz
        write(un_xpolz,fmt1reals)xpolz(i)
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
        write(un_sys,*)'nseg        = ',nseg
        write(un_sys,*)'lseg        = ',lseg
        write(un_sys,*)'chainperiod = ',chainperiod
        write(un_sys,*)'cuantas     = ',cuantas
        write(un_sys,*)'denspol     = ',denspol
       
        ! system description
        write(un_sys,*)'systype     = ',systype
        write(un_sys,*)'bcflag(LEFT)  = ',bcflag(LEFT)
        write(un_sys,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
        write(un_sys,*)'delta       = ',delta
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conv    = ',tol_conv
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
        write(un_sys,*)'pKa         = ',pKaA(1)      
        write(un_sys,*)'pKaNa       = ',pKaA(2)
        write(un_sys,*)'pKaACa      = ',pKaA(3)
        write(un_sys,*)'pKaA2Ca     = ',pKaA(4)
        write(un_sys,*)'pKb         = ',pKaB(1)      
        write(un_sys,*)'pKbNa       = ',pKaB(2)
        write(un_sys,*)'pKbBCa      = ',pKaB(3)
        write(un_sys,*)'pKbB2Ca     = ',pKaB(4)
        write(un_sys,*)'KionNa      = ',KionNa
        write(un_sys,*)'KionK       = ',KionK
        write(un_sys,*)'K0ionNa     = ',K0ionNa
        write(un_sys,*)'K0ionK      = ',K0ionK
        write(un_sys,*)'dielectW    = ',dielectW
        write(un_sys,*)'lb          = ',lb
        write(un_sys,*)'T           = ',Tref
        if(isVdW) then 
            write(un_sys,*)'VdWeps  = ',VdWeps
            write(un_sys,*)'VdWepsAA = ',VdWepsAA 
            write(un_sys,*)'VdWepsBB = ',VdWepsBB
            write(un_sys,*)'VdWepsAB = ',VdWepsAB
        endif    
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
    write(un_sys,*)'tol_conv    = ',tol_conv
    write(un_sys,*)'denspol     = ',denspol
    write(un_sys,*)'sumphiA     = ',sumphiA
    write(un_sys,*)'sumphiB     = ',sumphiB
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
    write(un_sys,*)'height      = ',height
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
    write(un_sys,*)'cuantas     = ',cuantas
    write(un_sys,*)'iterations  = ',iter
   
    ! .. closing files

    if(nz==nzmin) then 
        close(un_sys)
        close(un_xsol)
        close(un_psi)
        close(un_xpol) 
        close(un_fdisA)
        close(un_fdisB)
        close(un_xpolz)
        close(un_q)
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
        
end subroutine output_elect


subroutine output_neutral
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters    
    use field
    use energy
    use myutils, only : newunit
    use chains, only : isHomopolymer

    !     .. output file names         
    character(len=90) :: sysfilename     
    character(len=90) :: xsolfilename 
    character(len=90) :: xpolfilename 
    character(len=90) :: xpolzfilename 
    character(len=90) :: xprofilename 
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

    !     .. local arguments
    integer :: i, t
    character(len=100) :: fnamelabel
    character(len=20) :: rstr
    logical :: isopen
    real(dp) :: denspol

    !     .. executable statements 

    denspol=init_denspol()

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
    
    if(nz==nzmax) then 

        !     .. make label filename
        call make_filename_label(fnamelabel)
       
        !     .. make filenames 
        sysfilename='system.'//trim(fnamelabel)
        xpolfilename='xpol.'//trim(fnamelabel)
        xpolzfilename='xpolz.'//trim(fnamelabel)  
        xsolfilename='xsol.'//trim(fnamelabel)   
        xprofilename='xpro.'//trim(fnamelabel)   


        !      .. opening files
        open(unit=newunit(un_sys),file=sysfilename)   
        open(unit=newunit(un_xpol),file=xpolfilename)
        open(unit=newunit(un_xpolz),file=xpolzfilename)
        open(unit=newunit(un_xsol),file=xsolfilename)
        open(unit=newunit(un_xpro),file=xprofilename)
        
        
    else ! check that files are open
        inquire(unit=un_sys, opened=isopen)
        if(.not.isopen) write(*,*)"un_sys is not open" 
        inquire(unit=un_xpol, opened=isopen)
        if(.not.isopen) write(*,*)"un_xpol is not open"
        inquire(unit=un_xsol, opened=isopen)
        if(.not.isopen) write(*,*)"un_xsol is not open"
        inquire(unit=un_xpolz, opened=isopen)
        if(.not.isopen) write(*,*)"un_xpolz is not open"
        inquire(unit=un_xpro, opened=isopen)
        if(.not.isopen) write(*,*)"un_xpro is not open"
    endif    
        
     !   .. writting files
    !   .. this line seperates different distances 
      
    write(un_xsol,*)'#D    = ',nz*delta 
    write(un_xpol,*)'#D    = ',nz*delta       
    write(un_xpro,*)'#D    = ',nz*delta 
    write(un_xpolz,*)'#D    = ',nz*delta     

    do i=1,nsize    
       write(un_xpol,fmt3reals)xpol(i),(rhopol(i,t),t=1,nsegtypes)
       write(un_xsol,fmt1reals)xsol(i)
       write(un_xpro,fmt1reals)xpro(i)
       
    enddo
    do i=1,nz   
       write(un_xpolz,fmt1reals)xpolz(i)
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
        write(un_sys,*)'nseg        = ',nseg
        write(un_sys,*)'lseg        = ',lseg
        write(un_sys,*)'chainperiod = ',chainperiod
        write(un_sys,*)'cuantas     = ',cuantas
        write(un_sys,*)'maxnchainsrot      = ',maxnchainsrotations     
        write(un_sys,*)'maxnchainsrotxy    = ',maxnchainsrotationsxy  
        ! system description 
        write(un_sys,*)'systype     = ',systype
        write(un_sys,*)'delta       = ',delta  
        write(un_sys,*)'nsize       = ',nsize  
        write(un_sys,*)'nx          = ',nx
        write(un_sys,*)'ny          = ',ny
        write(un_sys,*)'ngr_freq    = ',ngr_freq
        write(un_sys,*)'ngr         = ',ngr
        write(un_sys,*)'nzmax       = ',nzmax
        write(un_sys,*)'nzmin       = ',nzmin
        write(un_sys,*)'nzstep      = ',nzstep
        write(un_sys,*)'tol_conv    = ',tol_conv
        ! other physcial parameters
        write(un_sys,*)'T           = ',Tref
        ! volume 
        write(un_sys,*)'vsol        = ',vsol
        write(un_sys,*)'vpol        = ',(vpol(t)*vsol,t=1,nsegtypes)
        write(un_sys,*)'vpro        = ',vpro
        write(un_sys,*)'cpro        = ',cpro%val
        write(un_sys,*)'isVdW       = ',isVdW
        
        write(un_sys,*)'===end distance independent settings=='
    endif
    
    write(un_sys,*)'distance    = ',nz*delta 
    write(un_sys,*)'nz          = ',nz
    write(un_sys,*)'nsize       = ',nsize
    write(un_sys,*)'free energy = ',FE  
    write(un_sys,*)'energy bulk = ',FEbulk 
    write(un_sys,*)'deltafenergy = ',deltaFE
    write(un_sys,*)'FEalt       = ',FEalt
    write(un_sys,*)'FEconf      = ',FEconf
    write(un_sys,*)'Econf       = ',Econf  
    write(un_sys,*)'FEtrans%sol = ',FEtrans%sol  
    write(un_sys,*)'FEtrans%pro = ',FEtrans%pro  
    write(un_sys,*)'fnorm       = ',fnorm
    write(un_sys,*)'sumphi      = ',(sumphi(t),t=1,nsegtypes)
    write(un_sys,*)'check phi   = ',checkphi
    write(un_sys,*)'denspol     = ',denspol
    write(un_sys,*)'FEq         = ',FEq 
    write(un_sys,*)'FEpi        = ',FEpi
    write(un_sys,*)'FErho       = ',FErho
    write(un_sys,*)'FEVdW       = ',FEVdW
    write(un_sys,*)'q           = ',q
    write(un_sys,*)'mu          = ',-log(q)
    write(un_sys,*)'height      = ',height
    write(un_sys,*)'iterations  = ',iter
    
    ! .. closing files
    if(nz.eq.nzmin) then 
        close(un_xsol)
        close(un_xpol)
        close(un_sys)
        close(un_xpolz)
        close(un_xpro)
    endif
  
end subroutine output_neutral


subroutine output_individualcontr_fe

    use globals, only : LEFT,RIGHT, systype
    use energy
    use myutils, only : newunit
    use volume, only : delta,nz,nzmax,nzmin

    ! local arguments

    character(len=100) :: fenergyfilename   
    character(len=100) :: fnamelabel
    character(len=20) :: rstr

    if(nz==nzmax) then 
        !     .. make label filename
        call make_filename_label(fnamelabel)
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
    write(un_fe,*)'FEVdW           = ',FEVdW
    write(un_fe,*)'FEconf          = ',FEconf
    write(un_fe,*)'Econf           = ',Econf
    write(un_fe,*)'Eshift          = ',Eshift
    
    
    write(un_fe,*)"FEtrans%sol     = ",FEtrans%sol   
    write(un_fe,*)"FEtrans%Na      = ",FEtrans%Na  
    write(un_fe,*)"FEtrans%Cl      = ",FEtrans%Cl  
    write(un_fe,*)"FEtrans%Ca      = ",FEtrans%Ca
    write(un_fe,*)"FEtrans%Mg      = ",FEtrans%Mg  
    write(un_fe,*)"FEtrans%K       = ",FEtrans%K
    write(un_fe,*)"FEtrans%KCl     = ",FEtrans%KCl
    write(un_fe,*)"FEtrans%NaCl    = ",FEtrans%NaCl  
    write(un_fe,*)"FEtrans%Hplus   = ",FEtrans%Hplus  
    write(un_fe,*)"FEtrans%OHmin   = ",FEtrans%OHmin  
    write(un_fe,*)"FEtrans%pro     = ",FEtrans%pro
    
    
    write(un_fe,*)"FEchempot%Na    = ",FEchempot%Na
    write(un_fe,*)"FEchempot%Cl    = ",FEchempot%Cl
    write(un_fe,*)"FEchempot%Ca    = ",FEchempot%Ca
    write(un_fe,*)"FEchempot%Mg    = ",FEchempot%Mg
    write(un_fe,*)"FEchempot%K     = ",FEchempot%K
    write(un_fe,*)"FEchempot%KCl   = ",FEchempot%KCl
    write(un_fe,*)"FEchempot%NaCl  = ",FEchempot%NaCl
    write(un_fe,*)"FEchempot%Hplus = ",FEchempot%Hplus
    write(un_fe,*)"FEchempot%OHmin = ",FEchempot%OHmin
    write(un_fe,*)"FEchempot%pro   = ",FEchempot%pro

    write(un_fe,*)"FEchemsurf(LEFT)     = ",FEchemsurf(LEFT)
    write(un_fe,*)"FEchemsurf(RIGHT)    = ",FEchemsurf(RIGHT)
    write(un_fe,*)"FEchemsurfalt(LEFT)  = ",FEchemsurfalt(LEFT)
    write(un_fe,*)"FEchemsurfalt(RIGHT) = ",FEchemsurfalt(RIGHT)

    write(un_fe,*)"delta FEchemsurfalt(LEFT) = ",FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
    write(un_fe,*)"delta FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)


    if(nz==nzmin) close(un_fe)

end subroutine output_individualcontr_fe


subroutine make_filename_label(fnamelabel)
 
    use globals, only : LEFT,RIGHT, systype, runtype
    use parameters, only : cNaCl,cCaCl2,cMgCl2,pHbulk,VdWepsBB,init_denspol,cpro,VdWscale

    character(len=*), intent(inout) :: fnamelabel    

    character(len=20) :: rstr
    real(dp) :: denspol


    denspol=init_denspol()

    !     .. make label filename

    select case(systype) 
    case("elect","electnopoly","electA")  


        write(rstr,'(F5.3)')denspol
        fnamelabel="phi"//trim(adjustl(rstr)) 
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

    case("neutral","neutralnoVdW")
        
        write(rstr,'(F5.3)')denspol
        fnamelabel="phi"//trim(adjustl(rstr)) 
        if(cpro%val>=0.001) then 
            write(rstr,'(F5.3)')cpro%val
        elseif(cpro%val>0.0) then  
            write(rstr,'(ES9.2E2)')cpro%val
        else 
            write(rstr,'(F3.1)')cpro%val
        endif 
        fnamelabel=trim(fnamelabel)//"cPro"//trim(adjustl(rstr))
        write(rstr,'(F5.3)')VdWscale%val
        fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))//".dat"

    case("brush_mul","brush_mulnoVdW","brushssdna","brushborn")  
        
        write(rstr,'(F5.3)')denspol
        fnamelabel="phi"//trim(adjustl(rstr)) 
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
    
        if(cMgCl2>=0.001) then 
            write(rstr,'(F5.3)')cMgCl2
        elseif(cMgCl2>0.0) then  
            write(rstr,'(ES9.2E2)')cMgCl2
        else 
            write(rstr,'(F3.1)')cMgCl2
        endif 
        fnamelabel=trim(fnamelabel)//"cMgCl2"//trim(adjustl(rstr))
    
        write(rstr,'(F7.3)')pHbulk
        fnamelabel=trim(fnamelabel)//"pH"//trim(adjustl(rstr))
        if(runtype=="rangeVdWeps")then
            write(rstr,'(F5.3)')VdWscale%val
            fnamelabel=trim(fnamelabel)//"VdWscale"//trim(adjustl(rstr))//".dat"
        else
            fnamelabel=trim(fnamelabel)//".dat"
        endif


    case default
        print*,"Error in output_individualcontr_fe subroutine"
        print*,"Wrong value systype : ", systype
    endselect

end subroutine

subroutine copy_solution(x)

    use globals, only : systype, neq, nsize, bcflag, LEFT, RIGHT
    use volume, only  : nx,ny
    use surface, only : psiSurfL, psiSurfR
    use field

    real(dp), dimension(neq) :: x  ! expliciet size array 
    
    ! local variable
    integer :: i, neq_bc
    integer, parameter :: A=1, B=2

    select case (systype)
    case ("elect")   

        do i=1,nsize                   
            xsol(i)= x(i)              
            psi(i) = x(i+nsize)        
            rhopol(i,A)=x(i+2*nsize)
            rhopol(i,B)=x(i+3*nsize)
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

    case ("neutral")

        do i=1,nsize                   
            xsol(i)= x(i)     
        enddo

    case ("neutralnoVdW")

        do i=1,nsize                   
            xsol(i)= x(i)     
        enddo 
           
    case ("brush_mulnoVdW")   

        do i=1,nsize                   
            xsol(i)= x(i)              
            psi(i) = x(i+nsize)
        enddo

    case default   

        print*,"Error: systype incorrect in copy_solution"
        print*,"stopping program"
        stop
    
    end select
         
end subroutine copy_solution

! output routine 

subroutine compute_vars_and_output()

    use globals, only : systype
    use energy
    use field
    use parameters, only : height

    select case (systype)
    case ("elect")
    
        call fcnenergy()   
        call charge_polymer()
        call average_charge_polymer()     
        call average_density_z(xpol,xpolz,height)  
        call output()
    
    case ("neutral","neutralnoVdW")  

        call fcnenergy() 
        call average_density_z(xpol,xpolz,height) 
        call output()           ! writing of output
    
    case ("brush_mul","brush_mulnoVdW","brushssdna","brushborn")

        call fcnenergy()   
        call charge_polymer() 
        call average_charge_polymer()     
        call average_density_z(xpol,xpolz,height)    
        call output()           ! writing of output

    case default   
        print*,"Error: systype incorrect in compute_vars_and_output"
        print*,"stopping program"
        stop
    end select
         
end subroutine compute_vars_and_output


subroutine write_chain_config()

    use mpivars, only : rank
    use globals, only: nseg, nsegtypes
    use chains, only : type_of_monomer,type_of_monomer_char,isAmonomer
    use parameters, only : lseg,lsegAA,vpol, vsol

    use myutils, only : newunit

    character(len=100) :: fname 
    integer :: i, un_cc
    character(len=10) ::istr
    
    if(rank==0) then
        write(istr,'(I4)')rank
        fname='chain_config.'//trim(adjustl(istr))//'.log'
        !     .. opening file        
        open(unit=newunit(un_cc),file=fname) 

        write(un_cc,*)"chain configuration summary"
        write(un_cc,*)"###"   
        write(un_cc,*)"lseg=",lseg
        write(un_cc,*)"#nsegtypes   lsegAA    vpol"
        do i=1,nsegtypes
            write(un_cc,*)i,lsegAA(i),vpol(i)*vsol
        enddo
        write(un_cc,*)"#segment type_num  type_char  isAmonomer"
        do i=1,nseg
            write(un_cc,*)i,type_of_monomer(i),type_of_monomer_char(i),isAmonomer(i)
        enddo 
        write(un_cc,*)"###"  
        close(un_cc)

    endif    


end subroutine write_chain_config


end module


