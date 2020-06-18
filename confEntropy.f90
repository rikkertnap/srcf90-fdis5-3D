! --------------------------------------------------------------|
! ConfEntropy.f90:                                              |
! constructs the free energy Conformational Entropy             |
!  beta Fconf = sigma \sum_alpha P(\alpha)lnP{\alpha)           |
! --------------------------------------------------------------|


module conform_entropy

    use mpivars
    use precision_definition
    implicit none

    private                    ! default all routines in this module private 
    public  ::  FEconf_entropy ! only this subroutine  public
   
contains

    subroutine FEconf_entropy(FEconf,Econf)

        use globals, only : systype
        use myutils, only : lenTEXT, print_to_log, LogUnit
    
        real(dp), intent(out) :: FEconf 
        real(dp), intent(out) :: Econf

        character(len=lenText) :: text 

        select case (systype) 
        case ("elect")
            call FEconf_elect(FEconf)
            Econf=0.0_dp
        case ("neutral","neutralnoVdW")
            call FEconf_neutral(FEconf,Econf)
        case ("brush_mul","brushssdna")
            call FEconf_brush_mul(FEconf,Econf)
        case ("brushborn")
            print*,"FEconfig_entropy not implemented yet"
        case default
            text="FEconf_entropy: wrong systype: "//systype//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        end select   

    end subroutine FEconf_entropy

   ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral(FEconf,Econf)
    
        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, energychain
        use field, only : xsol, rhopol, q
        use parameters, only : vpol, isVdW, VdWscale
        use VdW, only : VdW_contribution_exp

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: exppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            do t=1,nsegtypes
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            do i=1,nsize
                exppi(i,t) = xsol(i)**vpol(t)
            enddo    
        enddo      
       
        if(isVdW) then 
            do t=1,nsegtypes  
                call VdW_contribution_exp(rhopol,exppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
            
        do c=1,cuantas         ! loop over cuantas
            pro=exp(-VdWscale%val*energychain(c))     
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                pro = pro*exppi(k,t)
            enddo    
            FEconf_local=FEconf_local+pro*log(pro)
            Econf_local= Econf_local+pro*energychain(c)*VdWscale%val
        enddo
        
        ! communicate FEconf

        if(rank==0) then
            FEconf=FEconf_local
            Econf=Econf_local
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                FEconf=FEconf +FEconf_local
                Econf =Econf+Econf_local
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf = FEconf/q-log(q)  
            Econf = Econf/q  
        endif

    end subroutine FEconf_neutral


    subroutine FEconf_brush_mul(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, energychain
        use field, only : xsol,psi, fdis,rhopol,q
        use parameters
        use VdW, only : VdW_contribution_exp

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: exppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        real(dp) :: FEconf_local
        real(dp) :: Econf_local

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                do i=1,nsize                                              
                    exppi(i,t) = (xsol(i)**vpol(t))*exp(-zpol(t,2)*psi(i) )/fdis(i,t)   ! auxilary variable palpha
                enddo  
            else
                do i=1,nsize
                     exppi(i,t) = xsol(i)**vpol(t)
                enddo  
            endif   
        enddo      
       
        if(isVdW) then 
            do t=1,nsegtypes  
                call VdW_contribution_exp(rhopol,exppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
            
        do c=1,cuantas         ! loop over cuantas
            pro=exp(-energychain(c))     
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                pro = pro*exppi(k,t)
            enddo    
            FEconf_local=FEconf_local+pro*log(pro)
            Econf_local= Econf_local+pro*energychain(c)
        enddo
        
        ! communicate FEconf

        if(rank==0) then

            FEconf=FEconf_local
            Econf=Econf_local
         

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                
                FEconf=FEconf +FEconf_local
                Econf =Econf+Econf_local                
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf = FEconf/q-log(q)  
            Econf = Econf/q  
        endif

    end subroutine FEconf_brush_mul


    subroutine FEconf_elect(FEconf)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        real(dp), intent(out) :: FEconf
        
        !     .. declare local variables
        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s,g,gn         ! dummy indices
        real(dp) :: pro
        real(dp) :: FEconf_local
        real(dp) :: q_local
       
        ! .. executable statements 

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisA(:,1),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisB(:,1),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(fdisA(:,1), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)    
            call MPI_RECV(fdisB(:,1), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
        endif
            
        do i=1,nsize
              exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable
              exppiB(i)=(xsol(i)**vpolB(1))*exp(-zpolB(1)*psi(i))/fdisB(i,1) ! auxiliary variable
        enddo
       
        FEconf=0.0_dp

        FEconf_local= 0.0_dp

        do c=1,cuantas             ! loop over cuantas
        
            pro=1.0_dp                ! initial weight conformation 
            
            
            do s=1,nseg              ! loop over segments 
                k=indexchain(s,c)         
                if(isAmonomer(s)) then ! A segment 
                    pro = pro*exppiA(k)
                else
                    pro = pro*exppiB(k)
                endif
            enddo

           FEconf_local=FEconf_local+pro*log(pro)
                  
        enddo  

        ! communicate FEconfABL and FEcopnd_ABR

        if(rank==0) then

            FEconf=FEconf_local
            
            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, 1 , MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                FEconf=FEconf+FEconf_local
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf = FEconf/q-log(q)    
        endif

    end subroutine FEconf_elect




end module conform_entropy
