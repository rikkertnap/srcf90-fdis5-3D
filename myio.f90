! --------------------------------------------------------------|
! myio.f90:                                                     |
! --------------------------------------------------------------|

subroutine read_inputfile

    use globals  
    use parameters
    use surface 
  
    implicit none
  
    character(len=8) :: fname
    integer :: ios
    character(len=80) :: fcnname

    !     .. reading in of variables from file
 
  
    write(fname,'(A8)')'input.in'
    open(unit=1,file=fname,iostat=ios,status='old')
    if(ios > 0 ) then
        print*, 'Error opening file : iostat =', ios
        stop
    endif
  
    read(1,*)method
    read(1,*)sysflag
    read(1,*)bcflag(LEFT)
    read(1,*)bcflag(RIGHT)
    read(1,*)chainmethod
    read(1,*)chaintype
    read(1,*)sigmaABL
    read(1,*)sigmaABR
    read(1,*)sigmaC
    read(1,*)error             
    read(1,*)infile              ! guess  1==yes
    read(1,*)pHbulk
    read(1,*)KionNa
    read(1,*)KionK
    read(1,*)sigmaSurfL
    read(1,*)sigmaSurfR
    read(1,*)cNaCl
    read(1,*)cKCl
    read(1,*)cCaCl2
    read(1,*)pKa(1)           !   AH   <=> A- + H+ 
    read(1,*)pKa(2)           !   ANa  <=> A- + Na+  
    read(1,*)pKa(3)           !   ACa+ <=> A- + Ca2+ 
    read(1,*)pKa(4)           !   A2Ca <=> 2A- + Ca2+
    read(1,*)pKb(1)           !   BH   <=> B- + H+ 
    read(1,*)pKb(2)           !   BNa  <=> B- + Na+ 
    read(1,*)pKb(3)           !   BCa+ <=> B- + Ca2+   
    read(1,*)pKb(4)           !   B2Ca <=> 2B- + Ca2+   
    read(1,*)period
    read(1,*)nsize
    read(1,*)nsegAB
    read(1,*)cuantasAB
    read(1,*)nsegC
    read(1,*)cuantasC
    read(1,*)VdWepsC
    read(1,*)VdWepsB    
    read(1,*)VdWcutoff
    read(1,*)nzmax            ! max distance
    read(1,*)nzmin            ! min distance
    read(1,*)nzstep           ! step distance  
    read(1,*)verboseflag     

    call init_allowed_flags()
    write(fcnname,'(A14)')'read_inputfile'
   
    call check_value_sysflag(fcnname)
    call check_value_bcflag()
   
    ! override input bcflags 
    if(sysflag=="electdouble") then
        bcflag(LEFT)="cc"
        bcflag(RIGHT)="cc"  
    endif
    if(sysflag=="elect") then
        sigmaAB=sigmaABL
        sigmaABR=0.0d0  
    endif

end subroutine read_inputfile
 

subroutine output_elect(countfile)
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters
    use field
    use energy
    use surface 
    !     use endpoint
  
    implicit none
      
    !     .. scalar arguments
    
    integer :: countfile
    !     .. output file names       
    
    character(len=16) :: sysfilename     
    character(len=24) :: xsolfilename 
    character(len=26) :: xpolABfilename 
    character(len=25) :: xpolCfilename 
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
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals,fmt   

    !     .. local arguments
    
    integer :: i,j,k,n           ! dummy indexes
    
    !     .. executable statements 

    n=nz

    ! ..format specifiers 

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
   
    ! .. make filenames 
    
    write(sysfilename,'(A7,BZ,I5.5,A4)')'system.',countfile,'.dat'
    write(xpolABfilename,'(A7,BZ,I5.5,A4)')'xpolAB.',countfile,'.dat'
    write(xpolCfilename,'(A6,BZ,I5.5,A4)')'xpolC.',countfile,'.dat'
    write(xsolfilename,'(A5,BZ,I5.5,A4)')'xsol.', countfile,'.dat'
    write(xpolendfilename,'(A8,BZ,I5.5,A4)')'xpolend.', countfile,'.dat'
    write(xNafilename,'(A8,BZ,I5.5,A4)')'xNaions.', countfile,'.dat'
    write(xKfilename,'(A7,BZ,I5.5,A4)')'xKions.', countfile,'.dat'
    write(xCafilename,'(A8,BZ,I5.5,A4)')'xCaions.', countfile,'.dat'
    write(xNaClfilename,'(A13,BZ,I5.5,A4)') 'xNaClionpair.', countfile,'.dat'
    write(xKClfilename,'(A12,BZ,I5.5,A4)') 'xKClionpair.', countfile,'.dat'
    write(xClfilename,'(A8,BZ,I5.5,A4)') 'xClions.', countfile,'.dat'
    write(potentialfilename,'(A10,BZ,I5.5,A4)')  'potential.', countfile,'.dat'
    write(chargefilename,'(A7,BZ,I5.5,A4)') 'charge.', countfile,'.dat'
    write(xHplusfilename,'(A7,BZ,I5.5,A4)') 'xHplus.', countfile,'.dat'
    write(xOHminfilename,'(A7,BZ,I5.5,A4)')  'xOHmin.', countfile,'.dat'
    write(densfracAfilename,'(A13,BZ,I5.5,A4)')'densityAfrac.', countfile,'.dat'
    write(densfracBfilename,'(A13,BZ,I5.5,A4)')'densityBfrac.', countfile,'.dat'
    write(densfracionpairfilename,'(A19,BZ,I5.5,A4)')'densityfracionpair.', countfile,'.dat'
        
    !     .. opening files        
    open(unit=10,file=sysfilename)       
    open(unit=30,file=xsolfilename)
    open(unit=70,file=potentialfilename)

    if(sysflag/="electnopoly") then          
        open(unit=20,file=xpolABfilename)
        open(unit=21,file=xpolCfilename)
        open(unit=110,file=densfracAfilename) 
        open(unit=120,file=densfracBfilename) 
    endif       

    if(verboseflag=="yes") then    
        open(unit=50,file=xNafilename)
        open(unit=51,file=xKfilename)
        open(unit=55,file=xCafilename)
        open(unit=58,file=xNaClfilename)
        open(unit=57,file=xKClfilename)
        open(unit=59,file=densfracionpairfilename)
        open(unit=60,file=xClfilename)
        open(unit=70,file=potentialfilename)
        open(unit=80,file=chargefilename)
        open(unit=90,file=xHplusfilename)
        open(unit=100,file=xOHminfilename)
    endif

    ! .. writting files

    if(sysflag/="electnopoly") then 
        do i=1,n
            write(20,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
            write(21,fmt2reals)zc(i),xpolC(i)
            write(110,fmt6reals)zc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
            write(120,fmt6reals)zc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)
        enddo
    endif   
        
    write(70,*),0.0d0,psiSurfL
    do i=1,n
        write(30,*)zc(i),xsol(i)
        write(70,*)zc(i),psi(i)
        if(verboseflag=="yes") then 
            write(50,*)zc(i),xNa(i)
            write(51,*)zc(i),xK(i)
            write(55,*)zc(i),xCa(i)
            write(58,*)zc(i),xNaCl(i)
            write(57,*)zc(i),xKCl(i)
            write(59,*)zc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
            write(60,*)zc(i),xCl(i)
            write(80,*)zc(i),rhoq(i)
            write(90,*)zc(i),xHplus(i)
            write(100,*)zc(i),xOHmin(i)
        endif
    enddo    
    write(70,*),n*delta,psiSurfR  
   

    ! .. writing system information 
    
    write(10,*)'system      = spherical weakpolyelectrolyte brush'
    write(10,*)'chainmethod = ',chainmethod
    write(10,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
       write(10,*)'readinchains = ',readinchains
    endif
    write(10,*)'sysflag     = ',sysflag
    write(10,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(10,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(10,*)'free energy = ',FE
    write(10,*)'energy bulk = ',FEbulk 
    write(10,*)'deltafenergy = ',deltaFE
    write(10,*)'fnorm       = ',fnorm
    write(10,*)'q residual  = ',qres
    write(10,*)'error       = ',error
    write(10,*)'sigmaAB     = ',sigmaAB*delta
    write(10,*)'sumphiA     = ',sumphiA
    write(10,*)'sumphiB     = ',sumphiB
    write(10,*)'sumphiC     = ',sumphiC
    write(10,*)'check phi   = ',checkphi 
    write(10,*)'FEq         = ',FEq 
    write(10,*)'FEpi        = ',FEpi
    write(10,*)'FErho       = ',FErho
    write(10,*)'FEel        = ',FEel
    write(10,*)'FEelsurf(LEFT)  = ',FEelsurf(LEFT)
    write(10,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
    write(10,*)'FEbind      = ',FEbind
    write(10,*)'FEVdW       = ',FEVdW 
    write(10,*)'FEalt       = ',FEalt
    write(10,*)'qAB         = ',qAB
    write(10,*)'qC          = ',qC
    write(10,*)'muAB        = ',-dlog(qAB)
    write(10,*)'muC         = ',-dlog(qC)
    write(10,*)'nsegAB      = ',nsegAB
    write(10,*)'lsegAB      = ',lsegAB
    write(10,*)'nsegC       = ',nsegC
    write(10,*)'lsegC       = ',lsegC
    write(10,*)'period      = ',period
    write(10,*)'nz          = ',nz
    write(10,*)'delta       = ',delta 
    write(10,*)'distance    = ',nz*delta    
    write(10,*)'vsol        = ',vsol
    write(10,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(10,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(10,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(10,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(10,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(10,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(10,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(10,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(10,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(10,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(10,*)'vpolC       = ',vpolC*vsol
    write(10,*)'vNa         = ',vNa*vsol
    write(10,*)'vCl         = ',vCl*vsol
    write(10,*)'vCa         = ',vCa*vsol
    write(10,*)'vK          = ',vK*vsol
    write(10,*)'vNaCl       = ',vNaCl*vsol
    write(10,*)'vKCl        = ',vKCl*vsol
    write(10,*)'cNaCl       = ',cNaCl
    write(10,*)'cKCl        = ',cKCl
    write(10,*)'cCaCl2      = ',cCaCl2
    write(10,*)'pHbulk      = ',pHbulk
    write(10,*)'pKa         = ',pKa(1)      
    write(10,*)'pKaNa       = ',pKa(2)
    write(10,*)'pKaACa      = ',pKa(3)
    write(10,*)'pKaA2Ca     = ',pKa(4)
    write(10,*)'pKb         = ',pKb(1)      
    write(10,*)'pKbNa       = ',pKb(2)
    write(10,*)'pKbBCa      = ',pKb(3)
    write(10,*)'pKbB2Ca     = ',pKb(4)
    write(10,*)'KionNa      = ',KionNa
    write(10,*)'KionK       = ',KionK
    write(10,*)'K0ionNa     = ',K0ionNa
    write(10,*)'K0ionK      = ',K0ionK
    write(10,*)'xNabulk     = ',xbulk%Na
    write(10,*)'xClbulk     = ',xbulk%Cl
    write(10,*)'xKbulk      = ',xbulk%K
    write(10,*)'xNaClbulk   = ',xbulk%NaCl
    write(10,*)'xKClbulk    = ',xbulk%KCl
    write(10,*)'xCabulk     = ',xbulk%Ca
    write(10,*)'xHplusbulk  = ',xbulk%Hplus
    write(10,*)'xOHminbulk  = ',xbulk%OHmin
    write(10,*)'sigmaAB     = ',sigmaAB*delta
    write(10,*)'sigmaC      = ',sigmaC*delta
    write(10,*)'heightAB    = ',heightAB
    write(10,*)'heightC     = ',heightC
    write(10,*)'qpolA       = ',qpolA
    write(10,*)'qpolB       = ',qpolB
    write(10,*)'qpoltot     = ',qpol_tot
    write(10,*)'avfdisA(1)  = ',avfdisA(1)
    write(10,*)'avfdisA(2)  = ',avfdisA(2)
    write(10,*)'avfdisA(3)  = ',avfdisA(3)
    write(10,*)'avfdisA(4)  = ',avfdisA(4)
    write(10,*)'avfdisA(5)  = ',avfdisA(5)
    write(10,*)'avfdisB(1)  = ',avfdisB(1)
    write(10,*)'avfdisB(2)  = ',avfdisB(2)
    write(10,*)'avfdisB(3)  = ',avfdisB(3)
    write(10,*)'avfdisB(4)  = ',avfdisB(4)
    write(10,*)'avfdisB(5)  = ',avfdisB(5)
    write(10,*)'dielectW    = ',dielectW
    write(10,*)'lb          = ',lb
    write(10,*)'T           = ',T
    write(10,*)'VdWepsC     = ',VdWepsC*vpolC*vsol 
    write(10,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(10,*)'zpolA(1)    = ',zpolA(1)
    write(10,*)'zpolA(2)    = ',zpolA(2)
    write(10,*)'zpolA(3)    = ',zpolA(3)
    write(10,*)'zpolA(4)    = ',zpolA(4)
    write(10,*)'zpolB(1)    = ',zpolB(1)
    write(10,*)'zpolB(2)    = ',zpolB(2)
    write(10,*)'zpolB(3)    = ',zpolB(3)
    write(10,*)'zpolB(4)    = ',zpolB(4)
    write(10,*)'zpolB(5)    = ',zpolB(5)
    write(10,*)'zNa         = ',zNa
    write(10,*)'zCa         = ',zCa
    write(10,*)'zK          = ',zK
    write(10,*)'zCl         = ',zCl
    write(10,*)'nsize       = ',nsize  
    write(10,*)'cuantasAB   = ',cuantasAB
    write(10,*)'cuantasC    = ',cuantasC
    write(10,*)'iterations  = ',iter

    write(10,*)'bcflag(LEFT)  = ',bcflag(LEFT)
    write(10,*)'bcflag(RIGHT) = ',bcflag(RIGHT)
    write(10,*)'sigmaSurfL  = ',sigmaSurfL/((4.0d0*pi*lb)*delta)
    write(10,*)'sigmaSurfR  = ',sigmaSurfR/((4.0d0*pi*lb)*delta)
    write(10,*)'sigmaqSurfL = ',sigmaqSurfL/((4.0d0*pi*lb)*delta)
    write(10,*)'sigmaqSurfR = ',sigmaqSurfR/((4.0d0*pi*lb)*delta)
    write(10,*)'psiSurfL    = ',psiSurfL
    write(10,*)'psiSurfR    = ',psiSurfR

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

    close(10)
    close(30)
    close(70)

    if(sysflag/="electnopoly") then
        close(20)   
        close(21)
        close(100)
        close(120)
    endif   

    if(verboseflag=="yes") then 
        close(50)   
        close(55)
        close(58)
        close(59)
        close(60)
        close(70)
        close(80)
        close(90)
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
    !     use endpoint
  
    implicit none
    
  !     .. scalar arguments
  
  integer :: countfile
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
  
  character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

  !     .. local arguments
  
  integer :: i,j,k,n           ! dummy indexes
  
  !     .. executable statements 

  n=nz

  fmt2reals = "(2ES25.16)"  
  fmt3reals = "(3ES25.16)"  
  fmt4reals = "(4ES25.16)"  
  fmt5reals = "(5ES25.16)" 
  fmt6reals = "(6ES25.16)" 
  
  
  !     .. make filenames 
  
  write(sysfilename,'(A7,BZ,I5.5,A4)')'system.',countfile,'.dat'
  write(xpolABfilename,'(A7,BZ,I5.5,A4)')'xpolAB.',countfile,'.dat'
  write(xsolfilename,'(A5,BZ,I5.5,A4)')'xsol.', countfile,'.dat'
  write(xpolendfilename,'(A8,BZ,I5.5,A4)')'xpolend.', countfile,'.dat'
  write(xNafilename,'(A8,BZ,I5.5,A4)')'xNaions.', countfile,'.dat'
  write(xKfilename,'(A7,BZ,I5.5,A4)')'xKions.', countfile,'.dat'
  write(xCafilename,'(A8,BZ,I5.5,A4)')'xCaions.', countfile,'.dat'
  write(xNaClfilename,'(A13,BZ,I5.5,A4)') 'xNaClionpair.', countfile,'.dat'
  write(xKClfilename,'(A12,BZ,I5.5,A4)') 'xKClionpair.', countfile,'.dat'
  write(xClfilename,'(A8,BZ,I5.5,A4)') 'xClions.', countfile,'.dat'
  write(potentialfilename,'(A10,BZ,I5.5,A4)')  'potential.', countfile,'.dat'
  write(chargefilename,'(A7,BZ,I5.5,A4)') 'charge.', countfile,'.dat'
  write(xHplusfilename,'(A7,BZ,I5.5,A4)') 'xHplus.', countfile,'.dat'
  write(xOHminfilename,'(A7,BZ,I5.5,A4)')  'xOHmin.', countfile,'.dat'
  write(densfracAfilename,'(A13,BZ,I5.5,A4)')'densityAfrac.', countfile,'.dat'
  write(densfracBfilename,'(A13,BZ,I5.5,A4)')'densityBfrac.', countfile,'.dat'
  write(densfracionpairfilename,'(A19,BZ,I5.5,A4)')'densityfracionpair.', countfile,'.dat'
      
!     .. opening files        
           
  open(unit=10,file=sysfilename)   
  open(unit=20,file=xpolABfilename)
  open(unit=30,file=xsolfilename)
!  open(unit=40,file=xpolendfilename)
  open(unit=50,file=xNafilename)
  open(unit=51,file=xKfilename)
  open(unit=55,file=xCafilename)
  open(unit=58,file=xNaClfilename)
  open(unit=57,file=xKClfilename)
  open(unit=59,file=densfracionpairfilename)
  open(unit=60,file=xClfilename)
  open(unit=70,file=potentialfilename)
  open(unit=80,file=chargefilename)
  open(unit=90,file=xHplusfilename)
  open(unit=100,file=xOHminfilename)
  open(unit=110,file=densfracAfilename) 
  open(unit=120,file=densfracBfilename) 
  
  write(70,*),0.0d0,psiSurfL
  
  do i=1,n
     
     write(20,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
     write(30,*)zc(i),xsol(i)
!     write(40,*)zc(i),endpol(i)
     write(50,*)zc(i),xNa(i)
     write(51,*)zc(i),xK(i)
     write(55,*)zc(i),xCa(i)
     write(58,*)zc(i),xNaCl(i)
     write(57,*)zc(i),xKCl(i)
     write(59,*)zc(i),(xNaCl(i)/vNaCl)/(xNa(i)/vNa+xCl(i)/vCl+xNaCl(i)/vNaCl)
     write(60,*)zc(i),xCl(i)
     write(70,*)zc(i),psi(i)
     write(80,*)zc(i),rhoq(i)
     write(90,*)zc(i),xHplus(i)
     write(100,*)zc(i),xOHmin(i)
     write(110,fmt6reals)zc(i),fdisA(1,i),fdisA(2,i),fdisA(3,i),fdisA(4,i),fdisA(5,i)        
     write(120,fmt6reals)zc(i),fdisB(1,i),fdisB(2,i),fdisB(3,i),fdisB(4,i),fdisB(5,i)        
  enddo
  
  write(70,*),n*delta,psiSurfR
  
  !     .. system information 
  
  
  write(10,*)'system      = spherical weakpolyelectrolyte brush'
  write(10,*)'chainmethod = ',chainmethod
  write(10,*)'chaintype   = ',chaintype
  if(chainmethod.eq."FILE") then
     write(10,*)'readinchains = ',readinchains
  endif
  write(10,*)'sysflag     = ',sysflag




  write(10,*)'free energy = ',FE
  write(10,*)'energy bulk = ',FEbulk 
  write(10,*)'deltafenergy = ',deltaFE
  write(10,*)'fnorm       = ',fnorm
  write(10,*)'q residual  = ',qres
  write(10,*)'error       = ',error
  write(10,*)'sigmaABL    = ',sigmaABL*delta
  write(10,*)'sigmaABR    = ',sigmaABR*delta
  write(10,*)'sumphiA     = ',sumphiA
  write(10,*)'sumphiB     = ',sumphiB
  write(10,*)'check phi   = ',checkphi 
  write(10,*)'FEq         = ',FEq 
  write(10,*)'FEpi        = ',FEpi
  write(10,*)'FErho       = ',FErho
  write(10,*)'FEel        = ',FEel
  write(10,*)'FEelsurf(LEFT) = ',FEelsurf(LEFT)
  write(10,*)'FEelsurf(RIGHT) = ',FEelsurf(RIGHT)
  write(10,*)'FEbind      = ',FEbind
  write(10,*)'FEVdW       = ',FEVdW 
  write(10,*)'FEalt       = ',FEalt
  write(10,*)'qABL        = ',qABL
  write(10,*)'qABR        = ',qABR
  write(10,*)'muAB        = ',-dlog(qAB)
  write(10,*)'nsegAB      = ',nsegAB
  write(10,*)'lsegAB      = ',lsegAB
  write(10,*)'period      = ',period
  write(10,*)'nz          = ',nz
  write(10,*)'delta       = ',delta 
  write(10,*)'distance    = ',nz*delta   
  write(10,*)'vsol        = ',vsol
  write(10,*)'vpolA(1)    = ',vpolA(1)*vsol
  write(10,*)'vpolA(2)    = ',vpolA(2)*vsol
  write(10,*)'vpolA(3)    = ',vpolA(3)*vsol
  write(10,*)'vpolA(4)    = ',vpolA(4)*vsol
  write(10,*)'vpolA(5)    = ',vpolA(5)*vsol
  write(10,*)'vpolB(1)    = ',vpolB(1)*vsol
  write(10,*)'vpolB(2)    = ',vpolB(2)*vsol
  write(10,*)'vpolB(3)    = ',vpolB(3)*vsol
  write(10,*)'vpolB(4)    = ',vpolB(4)*vsol
  write(10,*)'vpolB(5)    = ',vpolB(5)*vsol
  write(10,*)'vNa         = ',vNa*vsol
  write(10,*)'vCl         = ',vCl*vsol
  write(10,*)'vCa         = ',vCa*vsol
  write(10,*)'vK          = ',vK*vsol
  write(10,*)'vNaCl       = ',vNaCl*vsol
  write(10,*)'vKCl        = ',vKCl*vsol
  write(10,*)'cNaCl       = ',cNaCl
  write(10,*)'cKCl        = ',cKCl
  write(10,*)'cCaCl2      = ',cCaCl2
  write(10,*)'pHbulk      = ',pHbulk
  write(10,*)'pKa         = ',pKa(1)      
  write(10,*)'pKaNa       = ',pKa(2)
  write(10,*)'pKaACa      = ',pKa(3)
  write(10,*)'pKaA2Ca     = ',pKa(4)
  write(10,*)'pKb         = ',pKb(1)      
  write(10,*)'pKbNa       = ',pKb(2)
  write(10,*)'pKbBCa      = ',pKb(3)
  write(10,*)'pKbB2Ca     = ',pKb(4)
  write(10,*)'KionNa      = ',KionNa
  write(10,*)'KionK       = ',KionK
  write(10,*)'K0ionNa     = ',K0ionNa
  write(10,*)'K0ionK      = ',K0ionK
  write(10,*)'xNabulk     = ',xbulk%Na
  write(10,*)'xClbulk     = ',xbulk%Cl
  write(10,*)'xKbulk      = ',xbulk%K
  write(10,*)'xNaClbulk   = ',xbulk%NaCl
  write(10,*)'xKClbulk    = ',xbulk%KCl
  write(10,*)'xCabulk     = ',xbulk%Ca
  write(10,*)'xHplusbulk  = ',xbulk%Hplus
  write(10,*)'xOHminbulk  = ',xbulk%OHmin
  write(10,*)'sigmaABL    = ',sigmaABL*delta
  write(10,*)'sigmaABR    = ',sigmaABR*delta
  write(10,*)'psiSurfL    = ',psiSurfL
  write(10,*)'psiSurfR    = ',psiSurfR
  write(10,*)'heightAB    = ',heightAB
  write(10,*)'qpolA       = ',qpolA
  write(10,*)'qpolB       = ',qpolB
  write(10,*)'qpoltot     = ',qpol_tot
  
  write(10,*)'avfdisA(1)  = ',avfdisA(1)
  write(10,*)'avfdisA(2)  = ',avfdisA(2)
  write(10,*)'avfdisA(3)  = ',avfdisA(3)
  write(10,*)'avfdisA(4)  = ',avfdisA(4)
  write(10,*)'avfdisA(5)  = ',avfdisA(5)
  write(10,*)'avfdisB(1)  = ',avfdisB(1)
  write(10,*)'avfdisB(2)  = ',avfdisB(2)
  write(10,*)'avfdisB(3)  = ',avfdisB(3)
  write(10,*)'avfdisB(4)  = ',avfdisB(4)
  write(10,*)'avfdisB(5)  = ',avfdisB(5)
  write(10,*)'dielectW    = ',dielectW
  write(10,*)'lb          = ',lb
  write(10,*)'T           = ',T 
  write(10,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
  write(10,*)'zpolA(1)    = ',zpolA(1)
  write(10,*)'zpolA(2)    = ',zpolA(2)
  write(10,*)'zpolA(3)    = ',zpolA(3)
  write(10,*)'zpolA(4)    = ',zpolA(4)
  write(10,*)'zpolB(1)    = ',zpolB(1)
  write(10,*)'zpolB(2)    = ',zpolB(2)
  write(10,*)'zpolB(3)    = ',zpolB(3)
  write(10,*)'zpolB(4)    = ',zpolB(4)
  write(10,*)'zpolB(5)    = ',zpolB(5)
  write(10,*)'zNa         = ',zNa
  write(10,*)'zCa         = ',zCa
  write(10,*)'zK          = ',zK
  write(10,*)'zCl         = ',zCl
  write(10,*)'nsize       = ',nsize  
  write(10,*)'cuantasAB   = ',cuantasAB
  write(10,*)'iterations  = ',iter
 


  close(10)
  close(20)
  close(21)
  close(30)
  close(50)
  close(55)
  close(58)
  close(59)
  close(60)
  close(70)
  close(80)
  close(90)
  close(100)
  close(120)
  
end subroutine output_electdouble

subroutine output_neutral(countfile)
  
    !     .. variables and constant declaractions
    use globals 
    use volume
    use parameters    
    use field
    use energy
    !     use endpoint
  
    implicit none
    
    !     .. scalar arguments
  
    integer :: countfile
    !     .. output file names       
  
    character(len=16) :: sysfilename     
    character(len=24) :: xsolfilename 
    character(len=26) :: xpolABfilename 
    character(len=25) :: xpolCfilename 
    character(len=26) :: xpolendfilename 
    
    character(len=80) :: fmt2reals,fmt3reals,fmt4reals,fmt5reals,fmt6reals   

    !     .. local arguments
    
    integer :: i,j,k,n           ! dummy indexes
    
    !     .. executable statements 

    n=nz

    fmt2reals = "(2ES25.16)"  
    fmt3reals = "(3ES25.16)"  
    fmt4reals = "(4ES25.16)"  
    fmt5reals = "(5ES25.16)" 
    fmt6reals = "(6ES25.16)" 
    
    
    !     .. make filenames 
    
    write(sysfilename,'(A7,BZ,I5.5,A4)')'system.',countfile,'.dat'
    write(xpolABfilename,'(A7,BZ,I5.5,A4)')'xpolAB.',countfile,'.dat'
    write(xpolCfilename,'(A6,BZ,I5.5,A4)')'xpolC.',countfile,'.dat'
    write(xsolfilename,'(A5,BZ,I5.5,A4)')'xsol.', countfile,'.dat'
    write(xpolendfilename,'(A8,BZ,I5.5,A4)')'xpolend.', countfile,'.dat'
        
  !     .. opening files        
             
    open(unit=10,file=sysfilename)   
    open(unit=20,file=xpolABfilename)
    open(unit=21,file=xpolCfilename)
    open(unit=30,file=xsolfilename)
  !  open(unit=40,file=xpolendfilename)
    
    do i=1,n
       
       write(20,fmt4reals)zc(i),xpolAB(i),rhopolA(i),rhopolB(i)
       write(21,fmt2reals)zc(i),xpolC(i)
       write(30,*)zc(i),xsol(i)
  !     write(40,*)zc(i),endpol(i)
    enddo
    
    !     .. system information 
    
    
    write(10,*)'system      = planar  brush'
    write(10,*)'free energy = ',FE  
    write(10,*)'energy bulk = ',FEbulk 
    write(10,*)'deltafenergy = ',deltaFE
    write(10,*)'FEalt       = ',FEalt
    write(10,*)'FEconfC     = ',FEconfC
    write(10,*)'FEconfAB    = ',FEconfAB
    write(10,*)'FEtrans%sol = ',FEtrans%sol  
    write(10,*)'fnorm       = ',fnorm
    write(10,*)'error       = ',error
    write(10,*)'sumphiA     = ',sumphiA
    write(10,*)'sumphiB     = ',sumphiB
    write(10,*)'sumphiC     = ',sumphiC
    write(10,*)'check phi   = ',checkphi 
    write(10,*)'FEq         = ',FEq 
    write(10,*)'FEpi        = ',FEpi
    write(10,*)'FErho       = ',FErho
    write(10,*)'FEVdW       = ',FEVdW
    write(10,*)'qAB         = ',qAB
    write(10,*)'qC          = ',qC
    write(10,*)'muAB        = ',-dlog(qAB)
    write(10,*)'muC         = ',-dlog(qC)
    write(10,*)'nsegAB      = ',nsegAB
    write(10,*)'lsegAB      = ',lsegAB
    write(10,*)'nsegC       = ',nsegC
    write(10,*)'lsegC       = ',lsegC
    write(10,*)'period      = ',period
    write(10,*)'nz          = ',nz
    write(10,*)'delta       = ',delta
    write(10,*)'distance    = ',nz*delta   
    write(10,*)'vsol        = ',vsol
    write(10,*)'vpolA(1)    = ',vpolA(1)*vsol
    write(10,*)'vpolA(2)    = ',vpolA(2)*vsol
    write(10,*)'vpolA(3)    = ',vpolA(3)*vsol
    write(10,*)'vpolA(4)    = ',vpolA(4)*vsol
    write(10,*)'vpolA(5)    = ',vpolA(5)*vsol
    write(10,*)'vpolB(1)    = ',vpolB(1)*vsol
    write(10,*)'vpolB(2)    = ',vpolB(2)*vsol
    write(10,*)'vpolB(3)    = ',vpolB(3)*vsol
    write(10,*)'vpolB(4)    = ',vpolB(4)*vsol
    write(10,*)'vpolB(5)    = ',vpolB(5)*vsol
    write(10,*)'vpolC       = ',vpolC*vsol
    write(10,*)'T           = ',T
    write(10,*)'VdWepsC     = ',VdWepsC*vpolC*vsol
    write(10,*)'VdWepsB     = ',VdWepsB*vpolB(3)*vsol
    write(10,*)'heightAB    = ',heightAB
    write(10,*)'heightC     = ',heightC
    write(10,*)'nsize       = ',nsize  
    write(10,*)'cuantasAB   = ',cuantasAB
    write(10,*)'cuantasC    = ',cuantasC
    write(10,*)'iterations  = ',iter
    write(10,*)'chainmethod = ',chainmethod
    write(10,*)'chaintype   = ',chaintype
    if(chainmethod.eq."FILE") then
       write(10,*)'readinchains = ',readinchains
        endif
    write(10,*)'sysflag     = ',sysflag


    close(10)
    close(20)
    close(21)
    close(30)
  
end subroutine output_neutral

subroutine output(countfile)

    use globals, only : sysflag
    implicit none

    integer :: countfile 

    if(sysflag=="elect") then 
        call output_elect(countfile)
    elseif(sysflag=="electdouble") then
        call output_electdouble(countfile)
    elseif(sysflag=="neutral") then
        call output_neutral(countfile)
    elseif(sysflag=="electnopoly") then
        call output_elect(countfile)
        call output_individualcontr_fe(countfile)
    else
        print*,"Error in output subroutine"
        print*,"Wrong value sysflag : ", sysflag
    endif     
end subroutine output


subroutine output_individualcontr_fe(countfile)

    use globals, only : LEFT,RIGHT
    use energy

    implicit none 

    integer, intent(in) :: countfile 

    character(len=16) :: fenergyfilename   

    write(fenergyfilename,'(A7,BZ,I5.5,A4)')'energy.', countfile,'.dat'
        
    !     .. opening files        

    open(unit=10,file=fenergyfilename) 

    write(10,*)'FE              = ',FE  
    write(10,*)'FEbulk          = ',FEbulk 
    write(10,*)'deltaFE         = ',deltaFE
    write(10,*)'FEalt           = ',FEalt  
    write(10,*)'FEbulkalt       = ',FEbulkalt 
    write(10,*)'deltaFEalt      = ',deltaFEalt
    
    write(10,*)"FEtrans%sol     = ",FEtrans%sol   
    write(10,*)"FEtrans%Na      = ",FEtrans%Na  
    write(10,*)"FEtrans%Cl      = ",FEtrans%Cl  
    write(10,*)"FEtrans%Ca      = ",FEtrans%Ca  
    write(10,*)"FEtrans%K       = ",FEtrans%K
    write(10,*)"FEtrans%KCl     = ",FEtrans%KCl
    write(10,*)"FEtrans%NaCl    = ",FEtrans%NaCl  
    write(10,*)"FEtrans%Hplus   = ",FEtrans%Hplus  
    write(10,*)"FEtrans%OHmin   = ",FEtrans%OHmin  
    
    write(10,*)"FEchempot%Na    = ",FEchempot%Na
    write(10,*)"FEchempot%Cl    = ",FEchempot%Cl
    write(10,*)"FEchempot%Ca    = ",FEchempot%Ca
    write(10,*)"FEchempot%K     = ",FEchempot%K
    write(10,*)"FEchempot%KCl   = ",FEchempot%KCl
    write(10,*)"FEchempot%NaCl  = ",FEchempot%NaCl
    write(10,*)"FEchempot%Hplus = ",FEchempot%Hplus
    write(10,*)"FEchempot%OHmin = ",FEchempot%OHmin

    write(10,*)"FEchemsurf(LEFT)= ",FEchemsurf(LEFT)
    write(10,*)"FEchemsurf(RIGHT)= ",FEchemsurf(RIGHT)
    write(10,*)"FEchemsurfalt(LEFT)= ",FEchemsurfalt(LEFT)
    write(10,*)"FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)

    write(10,*)"delta FEchemsurfalt(LEFT)= ",FEchemsurfalt(LEFT)-FEchemsurf(LEFT)-diffFEchemsurf(LEFT)
    write(10,*)"delta FEchemsurfalt(RIGHT)= ",FEchemsurfalt(RIGHT)-FEchemsurf(RIGHT)-diffFEchemsurf(RIGHT)



    close(10)


end subroutine   output_individualcontr_fe