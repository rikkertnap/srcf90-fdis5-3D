module surface 
   
     use globals, only : LEFT, RIGHT
     use mathconst
   
    implicit none
    
    !   different surface states

    real*8 :: fdisS(6)            ! fraction of different surface states
    real*8 :: KS(5)               ! experimemtal equilibruim constant 
    real*8 :: pKS(5)              ! experimental equilibruim constant pKS= -log[KS]   
    real*8 :: K0S(5)              ! intrinsic equilibruim constant
    real*8 :: qS(6)               ! charge 
    real*8 :: cap                 ! capacitance
 
    real*8 :: sigmaSurfL          ! surface density of acid on surface in nm^2
    real*8 :: sigmaSurfR          ! surface density of acid on surface in nm^2
    
    real*8 :: sigmaqSurfL         ! surface charge density on surface in nm^2
    real*8 :: sigmaqSurfR         ! surface charge density on surface in nm^2
    
    real*8 :: psiSurfL            ! surface potential     
    real*8 :: psiSurfR            ! surface potential 

!    real*8 :: sigmaSurf
    
   ! taurine
    real*8 :: fdisTaL(4),fdisTaR(4) ! fraction of different surface states
    real*8 :: KTa(3)               ! experimemtal equilibruim constant 
    real*8 :: pKTa(3)              ! experimental equilibruim constant pKS= -log[KS]   
    real*8 :: K0Ta(3)              ! intrinsic equilibruim constant
    real*8 :: qTa(4)               ! charge 

    contains
 
        subroutine init_surface(bc)
    
            implicit none

            character(len=2) :: bc(2)
        
            select case (bc(RIGHT))
                case ("qu")  
                    print*,"Hello: case quartz selected"
                    call  init_surface_quartz()
                case ("cl")
                    call  init_surface_clay()
                case ("ca")  
                    call init_surface_calcite()
                case ("ta") 
                    call init_surface_taurine(RIGHT)
                case ("cc")
                    call init_surface_constcharge(RIGHT)
                case default
                    print*,"bc(RIGHT) does not match qu, cl, ca, ta, or cc"
            end select 
                
            select case (bc(LEFT))
                case ("ta") 
                    call init_surface_taurine(LEFT)
                case ("cc")
                    call init_surface_constcharge(LEFT)
                case default
                    print*,"bc(LEFT) does not match ta, or cc"
            end select 
        end subroutine init_surface
   
        real*8 function surface_charge(bc,psiSurf,side)
        
            implicit none
            
            real*8 :: psiSurf
            character(len=2) :: bc
            integer :: side 

            real*8 :: sigmaqSurf
            
            if(side==RIGHT) then 
                select case (bc)
                    case ("qu") 
                        sigmaqSurf = surface_charge_quartz(psiSurf)
                    case ("cl") 
                        sigmaqSurf = surface_charge_clay(psiSurf)
                    case ("ca") 
                        sigmaqSurf = surface_charge_calcite(psiSurf)
                    case ("ta")
                        sigmaqSurf = surface_charge_taurine(psiSurf,RIGHT)
                    case ("cc")
                        sigmaqSurf = sigmaSurfR
                    case default
                        print*,"bc does not match qu, cl, ca, or ta"    
                        sigmaqSurf = 0.0d0
                end select
            elseif (side==LEFT) then
                select case (bc)
                    case ("ta")
                        sigmaqSurf = surface_charge_taurine(psiSurf,LEFT)
                    case ("cc")
                        sigmaqSurf = sigmaSurfL
                    case default
                        print*,"bc does not match qu, cl, ca, or ta"    
                        sigmaqSurf = 0.0d0
                end select
            else 
                print*,"Error : side value is not equal to  LEFT or RIGHT"
                stop
            endif
                    
            surface_charge=sigmaqSurf

            return

        end function surface_charge
     
        subroutine init_surface_quartz()
       
            use mathconst 
            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            implicit none

            integer :: i
        
            cap=1.0d-18 ! in F/nm^2 1F= 1C/V
     
            pKS(1)=   8.1d0  !  8.1d0 !  >SOH <=> >SO- + H+ ! see Luetzenkirchen book chapter 14
            pKS(2)=  -4.4d0  ! -4.4d0 !  >SOH2+ <=> >SOH + H+
            pKS(3)=  -1.5d0  ! >SONa <=> >SO- + Na+  
            pKS(4)=  -5.20d0 ! >SOCa+ <=> >SO- + Ca2+  Ber. Bunsenges. Phys. Chem., 98, pp 1062-1067, 1994 
            pKS(5)=  -1.5d0  ! >SOH2Cl <=>>SOH2+ Cl- 

            do i=1,5
                KS(i)  = 10d0**(-pKS(i))       ! experimental equilibruim constant surface acid
                K0S(i) = (KS(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
            enddo

            ! charges surface states
            qS(1)=-1.0d0 ! SO^-
            qS(2)=0.0d0  ! SOH
            qS(3)=1.0d0  ! SOH_2^+
            qS(4)=0.0d0  ! SONa
            qS(5)=1.0d0  ! SOCa
            qS(6)=0.0d0  ! SOH2Cl
        
            ! sites density
            sigmaSurfR=8.0d0
            sigmaSurfR = sigmaSurfR * (4.0d0*pi*lb)*delta ! dimensionless surface charge     

        end subroutine init_surface_quartz

        subroutine init_surface_calcite()
            
            use mathconst
            use physconst, only : Na
            use parameters,  only : vsol,vCa,delta,lb

            implicit none

            integer :: i

        !    cap=1.0d-18 ! in F/nm^2 1F= 1C/V                                                                                                                              
            !   1=CO3-, 2=CO3H, 3=CO3Ca+, 4=CaO^-, 5=CaOH, 6= CaOH2^+ 
            
            pKS(1)=   4.9d0   ! >CO3H         <=> >CO3^- + H+ 
            pKS(2)=   2.8d0   ! >CO3H + Ca^2+ <=> >CO3Ca^+ + H+ 
            pKS(3)=   12.2d0  ! >CaOH2+       <=> >CaOH + H+     
            pKS(4)=   17.0d0  ! >CaOH         <=> >CaO- + H+                                                                                                                     

            do i=1,4
               KS(i)  = 10d0**(-pKS(i))       ! experimental equilibruim constant surface acid                                                                           
               K0S(i) = (KS(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant                                                                                       
            enddo
            ! i=2 is different
            K0S(2)=KS(2)

            ! charges surface states 
            qS(1)=-1.0d0 ! >CO3^-
            qS(2)=0.0d0  ! >CO3H                                                
            qS(3)=1.0d0  ! >CO3Ca^+

            qS(4)=-1.0d0 ! >CaO^-                              
            qS(5)=0.0d0  ! >CaOH 
            qS(6)=1.0d0  ! >CaOH2^+

            ! site density
            sigmaSurfR=5.0d0
            sigmaSurfR = sigmaSurfR * (4.0d0*pi*lb)*delta ! dimensionless surface charge     

        end subroutine init_surface_calcite

       
        subroutine init_surface_clay()

            use mathconst 
            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            implicit none

            integer :: i
            
            ! see Luetzenkirchen book chapter 7 pages 204-206                                     

            cap=0.90d-18 ! in F/nm^2 1F= 1C/V                                                                                         
            pKS(1)=   10.0d0  ! 10.0  >SOH^0.5 <=> >SO^-0.5 + H+ 
            pKS(2)=   5.3d0   !  5.3  >SOH_2^0.5 <=> >SOH_1.5^0 + 1/2H+ 
            pKS(3)=   0.20d0  !   >SONa^0.5 <=> >SO^-0.5 + Na+ 
            pKS(4)=  -5.20d0  !    >SOCa^1.5 <=> >SO^-0.5 + Ca2+             
            pKS(5)=   -0.40d0 !   >SOH_2Cl^-0.5 <=> >SOH_2^0.5+Cl-         
           
            do i=1,5
               KS(i)  = 10d0**(-pKS(i))       ! experimental equilibruim constant surface acid  
               K0S(i) = (KS(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
            enddo
            ! i=2
            !K0S(2) = KS(2)*((vsol*Na/1.0d24)**0.5) ! intrinstic equilibruim constant                                                  

            ! charges surface states                                                                                                 
            qS(1)=-0.5d0 ! SOH^-0.5                                                                            
            qS(2)=0.5d0  ! SOH_2^0.5   
            qS(3)=0.0d0  ! SOH_1.5^0  
            qS(4)=0.5d0  ! SOHNa^0.5                                  
            qS(5)= 1.5d0 ! SOHCa^1.5
            qS(6)=-0.5d0 ! SOH_2Cl^-0.5                                                                                                        

            ! site density 
            sigmaSurfR = 8.5d0
            sigmaSurfR = sigmaSurfR * (4.0d0*pi*lb)*delta ! dimensionless surface charge     

        end subroutine init_surface_clay

        subroutine init_surface_taurine(side)

            use mathconst
            use physconst, only : Na
            use parameters,  only : vsol,delta,lb

            implicit none
                
            integer, intent(in) :: side        

            integer :: i

            pKTa(1)=   -2.0d0    !  >SOH <=> >SO^- + H+ 
            pKTa(2)=   -0.42d0   !  >SONa <=> >SO^- + Na+
            pKTa(3)=   -0.72243d0  !  >SOCa^+ <=> >SO^- +Ca2+
           
            do i=1,3
                KTa(i)  = 10d0**(-pKTa(i))       ! experimental equilibruim constant surface acid  
                K0Ta(i) = (KTa(i)*vsol)*(Na/1.0d24) ! intrinstic equilibruim constant 
            enddo

            ! charges surface states                                                                                                 
            qTa(1)=-1.0d0 ! SO^-                                                                            
            qTa(2)=0.0d0  ! SOH  
            qTa(3)=0.0d0  ! SONa 
            qTa(4)=1.0d0  ! SOHCa^+                                  
                                                                                                                

            ! site density
            if(side==RIGHT) sigmaSurfR = sigmaSurfR * (4.0d0*pi*lb)*delta ! dimensionless surface charge     
            if(side==LEFT)  sigmaSurfL = sigmaSurfL * (4.0d0*pi*lb)*delta ! dimensionless surface charge     

        end subroutine init_surface_taurine

        subroutine init_surface_constcharge(side)

            use parameters,  only : delta,lb

            implicit none

            integer, intent(in) :: side
            
            ! site density
            if(side==RIGHT) sigmaSurfR = sigmaSurfR * (4.0d0*pi*lb)*delta ! dimensionless surface charge     
            if(side==LEFT)  sigmaSurfL = sigmaSurfL * (4.0d0*pi*lb)*delta ! dimensionless surface charge     

        end subroutine init_surface_constcharge

        real*8 function surface_charge_quartz(psiS)
     
            use physconst
            use mathconst
            use parameters
        
            implicit none
        
            real*8 :: psiS
         
            ! .. local variables
        
            real*8 :: xS(6)
            real*8 :: A,avfdis
            integer :: i
        
            xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! SOH/SO-
            xS(2)= ((xbulk%Hplus/xbulk%sol)/K0S(2))*dexp(-psiS) ! SOH2+/SOH
            xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*dexp(-psiS) ! SONa/SO-
            xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*dexp(-2.0d0*psiS) ! SOCa+/SO-
            xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*dexp(psiS) ! SOH2Cl/SOH2+

            A = xS(1)*(1.0d0+xS(2))+xS(3)+xS(4)+xS(1)*xS(2)*xS(5)
            fdisS(1)  = 1.0d0/(1.0d0 +A) ! SO-
            fdisS(2)  = fdisS(1)*xS(1) ! SOH
            fdisS(3)  = fdisS(2)*xS(2) ! SOH2+
            fdisS(4)  = fdisS(1)*xS(3) ! SONa
            fdisS(5)  = fdisS(1)*xS(4) ! SOCa+
            fdisS(6)  = fdisS(1)*xS(5)*xS(2)*xS(1)  ! SOH2Cl  
        
            avfdis=0.0d0
            do i=1,6
                avfdis=avfdis +qS(i)*fdisS(i)
            enddo
        
            surface_charge_quartz=sigmaSurfR*avfdis
        
        end function surface_charge_quartz


        real*8 function surface_charge_calcite(psiS)

            use physconst
            use mathconst
            use parameters

            implicit none

            real*8 :: psiS

            ! .. local variables                                                                                                                                          
            real*8 :: xS(6)
            real*8 :: A,avfdis
            integer :: i

            xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! CO3H/CO3- =fdisS(2)/fdisS(1)

            xS(2)= (xbulk%sol/xbulk%Hplus)*((xbulk%Ca*vsol/vCa)/(xbulk%sol**(vCa/vsol)))*K0S(2)*dexp(-psiS)   ! CO3Ca+/CO3H = fdisS(3)/fdisS(2)             
            xS(3)= ((xbulk%Hplus/xbulk%sol)/K0S(4))*dexp(-psiS) ! CaOH/CaO- =fdisS(5)/fdisS(4)   
            xS(4)= ((xbulk%Hplus/xbulk%sol)/K0S(3))*dexp(-psiS) ! CaOH2^+/CaOH = fdisS(6)/fdisS(5)

            A = xS(1)*(1.0d0+xS(2))
            fdisS(1)  = 1.0d0/(1.0d0 +A) ! CO3^- 
            fdisS(2)  = fdisS(1)*xS(1)   ! CO3H 
            fdisS(3)  = fdisS(2)*xS(2)   ! CO3Ca+ 
        
            A=xS(3)*(1.0d0+xS(4))
            fdisS(4)  = 1.0d0/(1.0d0 +A) ! CaO^- 
            fdisS(5)  = fdisS(4)*xS(3)   ! CaOH  
            fdisS(6)  = fdisS(5)*xS(4)   ! CaOH2^+                                                                                                              
            avfdis=0.0d0
            do i=1,6
                avfdis=avfdis +qS(i)*fdisS(i)
            enddo

            surface_charge_calcite=sigmaSurfR*avfdis

        end function surface_charge_calcite

      
        real*8 function surface_charge_clay(psiS)

            use physconst
            use mathconst
            use parameters

            implicit none

            real*8 :: psiS

            ! .. local variables                                                                                                  

            real*8 :: xS(6)
            real*8 :: A,avfdis
            integer :: i

            xS(1)= ((xbulk%Hplus/xbulk%sol)/K0S(1))*dexp(-psiS) ! SOH_2^0.5/SOH^-0.5                       
            xS(2)= (((xbulk%Hplus/xbulk%sol)**0.5)/K0S(2))*dexp(-0.5d0*psiS) ! SOH_1.5^0/SOH_2^0.5   
            xS(3)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0S(3))*dexp(-psiS) ! SOHNa^0.5/SOH^-0.5
            xS(4)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0S(4))*dexp(-2.0d0*psiS) ! SOHCa^1.5/SOH^-0.5
            xS(5)= (((xbulk%Cl/vCl)/(xbulk%sol**(vCl)))/K0S(5))*dexp(psiS) ! SOH2Cl^-0.5/SOH_2^0.5

            xS(2)=0.0 ! make zero 

            A = xS(1)*(1.0d0+xS(2))+xS(3)+xS(4)+xS(1)*xS(5)
            fdisS(1)  = 1.0d0/(1.0d0 +A) ! SOH^-0.5                                                                                 
            fdisS(2)  = fdisS(1)*xS(1) ! SOH_2^0.5 
            fdisS(3)  = fdisS(2)*xS(2) ! SOH_1.5^0                                                                                  
            fdisS(4)  = fdisS(1)*xS(3) ! SOHNa^0.5 
            fdisS(5)  = fdisS(1)*xS(4) ! SOHCa^1.5            
            fdisS(6)  = fdisS(1)*xS(1)*xS(5)  ! SOH_2Cl^-0.5                                                

            avfdis=0.0d0
            do i=1,6
                avfdis=avfdis +qS(i)*fdisS(i)
            enddo

            surface_charge_clay=sigmaSurfR*avfdis

        end function surface_charge_clay

        real*8 function surface_charge_taurine(psiS,side)

            use physconst
            use mathconst
            use parameters

            implicit none

            real*8 :: psiS
            integer :: side
        
            ! .. local variables                                                                                                  

            real*8 :: xS(3)
            real*8 :: A,avfdis
            integer :: i

            xS(1)= ((xbulk%Hplus/xbulk%sol)/K0Ta(1))*dexp(-psiS) ! SOH/SO^-                       
            xS(2)= (((xbulk%Na/vNa)/(xbulk%sol**(vNa)))/K0Ta(2))*dexp(-psiS) ! SO^-/SONa
            xS(3)= (((xbulk%Ca/vCa)/(xbulk%sol**(vCa)))/K0Ta(3))*dexp(-2.0d0*psiS) ! SO^-/SOCa^+
         
            A = xS(1)+xS(2)+xS(3)
         
            if(side==LEFT) then 
                fdisTaL(1)  = 1.0d0/(1.0d0 +A) ! SO^-
                fdisTaL(2)  = fdisTal(1)*xS(1) ! SOH                                                                                 
                fdisTal(3)  = fdisTaL(1)*xS(2) ! SONa 
                fdisTaL(4)  = fdisTaL(1)*xS(3) ! SOCa^+                                                                                 

                avfdis=0.0d0
                do i=1,4
                    avfdis=avfdis +qTa(i)*fdisTaL(i)   
                enddo
            
                surface_charge_taurine=sigmaSurfL*avfdis
            
            else if(side==RIGHT) then
                fdisTaR(1)  = 1.0d0/(1.0d0 +A) ! SO^-
                fdisTaR(2)  = fdisTaR(1)*xS(1) ! SOH                                                                                 
                fdisTaR(3)  = fdisTaR(1)*xS(2) ! SONa 
                fdisTaR(4)  = fdisTaR(1)*xS(3) ! SOCa^+                                                                                 

                avfdis=0.0d0
                do i=1,4
                    avfdis=avfdis +qTa(i)*fdisTaR(i)
                enddo
                surface_charge_taurine=sigmaSurfR*avfdis
            endif
        
        end function surface_charge_taurine


end module surface
