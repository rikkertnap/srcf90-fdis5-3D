! --------------------------------------------------------------|
! ConfEntropy.f90:                                              |
! constructs the free energy Conformational Entropy             |
!  beta Fconf = sigma \sum_alpha P(\alpha)lnP{\alpha)           |
! --------------------------------------------------------------|


module conform_entropy

    implicit none

    private                    ! default all routines in this module private 
    public  ::  FEconf_entropy ! only this subroutine  public
   

contains

    subroutine FEconf_entropy(FEconfAB,FEconfC)

 
        use globals
        use myutils
        
        implicit none
    
        real(dp), intent(out) :: FEconfAB, FEconfC  

        character(len=lenText) :: text 

        select case (sysflag) 
        case ("elect")
            call FEconf_elect(FEconfAB)
            FEconfC=0.0_dp
        case ("electdouble")
            call FEconf_electdouble(FEconfAB)
            FEconfC=0.0_dp
        case ("neutral")
            call FEconf_neutral(FEconfAB,FeconfC)
        case ("electHC")     
            text="fcenergy: sysflag: "//sysflag//" not implemented yet"
            call print_to_log(LogUnit,text)
            print*,text
        case default
            text="FEconf_entropy: wrong sysflag: "//sysflag//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        end select   

    end subroutine FEconf_entropy

  ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral(FEconfAB,FeconfC)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW

        implicit none
        
        real(dp), intent(out) :: FEconfAB, FEconfC  

        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize),exppiC(nsize)    ! auxilairy variable for computing P(\alpha) 

        integer :: i,j,k,c,s         ! dummy indices
        real(dp) :: pro,tmp,expVdW 
        integer :: conf              ! counts number of conformations
        

        real(dp), parameter :: tolconst = 1.0e-9_dp  ! tolerance for constA and constB 


        !     .. executable statements 

        do i=1,nz    
            exppiA(i)=(xsol(i)**vpolA(3)) !*dexp(-zpolA(3)*psi(i))/fdisA(3,i) ! auxiliary variable
            exppiB(i)=(xsol(i)**vpolB(3)) !*dexp(-zpolB(3)*psi(i))/fdisB(3,i) ! auxiliary variable   
            exppiC(i)=(xsol(i)**vpolC)
       
            !     .. VdW interaction   
            tmp = 0.0_dp
            if((i+VdWcutoffdelta)<=nsize) then 
                do j=minrange(i),i+VdWcutoffdelta
                    tmp = tmp + chis(i,j)*rhopolB(j)*vpolB(3)*vsol
                enddo
            endif
            expVdW=dexp(-VdWepsB*tmp)
            exppiB(i)=exppiB(i)*expVdW ! auxiliary variable
        enddo

        
        FEconfAB=0.0_dp

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0_dp                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo
            FEconfAB=FEconfAB+pro*log(pro)
        enddo

        ! normalize
        FEconfAB=(FEconfAB/qAB-log(qAB))*(sigmaAB*delta)    

        FEconfC = 0.0_dp                   
        do c=1,cuantasC            ! loop over cuantas                                                      
            pro=1.0_dp               ! initial weight conformation                                                   
            do s=1,nsegC            ! loop over segments                
                k=indexchainC(c,s)
                pro = pro*exppiC(k)
            enddo
            FEconfC = FEconfC+pro*log(pro)
        enddo
        ! normalize
        FEconfC=(FEConfC/qC -log(qC))*(sigmaC*delta)    

    end subroutine FEconf_neutral



    subroutine FEconf_elect(FEconfAB)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconfAB
        
        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer :: i,k,c,s         ! dummy indices
        real(dp) :: pro


        !     .. executable statements 

       
        do i=1,nz
              exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
              exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable
        enddo
        

        FEconfAB=0.0_dp

        do c=1,cuantasAB            ! loop over cuantas
            pro=1.0_dp                ! initial weight conformation 
            do s=1,nsegAB            ! loop over segments 
                k=indexchainAB(c,s)
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo
            FEconfAB=FEconfAB+pro*log(pro)
        enddo

        ! normalize

        FEconfAB=(FEconfAB/qAB-log(qAB))*(sigmaAB*delta)     ! check sigma

        
    end subroutine FEconf_elect


    subroutine FEconf_electdouble(FEconfAB)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconfAB
        
        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer :: i,k,c,s, kL,kR         ! dummy indices
        real(dp) :: proL,proR
        real(dp) :: FEconfABL,FEconfABR


        !     .. executable statements 

       
        do i=1,nz
              exppiA(i)=(xsol(i)**vpolA(1))*dexp(-zpolA(1)*psi(i))/fdisA(1,i) ! auxiliary variable
              exppiB(i)=(xsol(i)**vpolB(1))*dexp(-zpolB(1)*psi(i))/fdisB(1,i) ! auxiliary variable
        enddo
       
       
        FEconfABL=0.0_dp
        FEconfABR=0.0_dp

        do c=1,cuantasAB            ! loop over cuantas
      
            proL=1.0_dp                ! initial weight conformation 
            proR=1.0_dp

            do s=1,nsegAB
                kL=indexchainAB(c,s)
                kR=nz+1-kL

                if(isAmonomer(s)) then ! A segment 
                    proL = proL*exppiA(kL)
                    proR = proR*exppiA(kR)
                else
                    proL = proL*exppiB(kL)
                    proR = proR*exppiB(kR)
                endif
            enddo

            FEconfABL=FEconfABL+proL*log(proL)
            FEconfABR=FEconfABR+proR*log(proR)

        enddo

        ! normalize
        
        FEconfABL=(FEconfABL/qABL-log(qABL))*(sigmaABL*delta)     ! check sigma
        FEconfABR=(FEconfABR/qABR-log(qABR))*(sigmaABR*delta)     ! check sigma


        FEconfAB=FEconfABL+FEconfABR

        
    end subroutine FEconf_electdouble
    

end module conform_entropy
