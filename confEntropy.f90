! --------------------------------------------------------------|
! ConfEntropy.f90:                                              |
! constructs the free energy Conformational Entropy             |
!  beta Fconf = sigma \sum_alpha P(\alpha)lnP{\alpha)           |
! --------------------------------------------------------------|


module conform_entropy

    use mpivars
    implicit none

    private                    ! default all routines in this module private 
    public  ::  FEconf_entropy ! only this subroutine  public
   

contains

    subroutine FEconf_entropy(FEconf,Econf)

        use globals
        use myutils
    
        real(dp), intent(out) :: FEconf 
        real(dp), intent(out) :: Econf


        character(len=lenText) :: text 


        select case (systype) 
        case ("elect")
            call FEconf_elect(FEconf)
            Econf=0.0_dp
        case ("electA")
            call FEconf_electA(FEconf) 
            Econf=0.0_dp
        case ("electdouble")
            call FEconf_electdouble(FEconf) 
            Econf=0.0_dp
        case ("neutral")
            call FEconf_neutral(FEconf)
            Econf=0.0_dp
        case ("electnopoly")
            FEconf=0.0_dp
            Econf=0.0_dp
        case ("brush_mul")
            call FEconf_brush_mul(FEconf,Econf)
        case default
            text="FEconf_entropy: wrong systype: "//systype//"stopping program"
            call print_to_log(LogUnit,text)
            stop
        end select   

    end subroutine FEconf_entropy

  ! computes conformational entropy in neutral state 

    subroutine FEconf_neutral(FEconf)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters
        use VdW

        implicit none
        
        real(dp), intent(out) :: FEconf  

        !     .. declare local variables
        ! still to be implemented 

        FEconf=0.0_dp

    end subroutine FEconf_neutral


    subroutine FEconf_brush_mul(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use volume, only : ngr, ngr_node
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable ,exp_energychain
        use field, only : xsol,psi, fdis,rhopol,q
        use parameters
        use VdW, only : VdW_contribution_exp

        implicit none

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: exppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        
        real(dp) :: FEconf_local(ngr_node)
        real(dp) :: FEconf_array(size*ngr_node)
        real(dp) :: Econf_local(ngr_node)
        real(dp) :: Econf_array(size*ngr_node)


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
       
        FEconf=0.0_dp!

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
            FEconf_local(gn)= 0.0_dp !init FEconf
            Econf_local=0.0_dp ! init FEconf
            
            g=gn+rank*ngr_node
            do c=1,cuantas         ! loop over cuantas
                pro=exp_energychain(gn,c)     
                do s=1,nseg        ! loop over segments 
                    k=indexchain(s,gn,c)
                    t=type_of_monomer(s)                
                    pro = pro*exppi(k,t)
                enddo    
                FEconf_local(gn)=FEconf_local(gn)+pro*log(pro)
                Econf_local= Econf_local-(pro)*log(exp_energychain(gn,c))
            enddo
        enddo   

        ! communicate FEconf

        if(rank==0) then

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                FEconf_array(g)=FEconf_local(gn)
                Econf_array(g)=Econf_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    FEconf_array(g)=FEconf_local(gn)
                    Econf_array(g)=Econf_local(gn)
                enddo
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g)  
            enddo    
        endif

    end subroutine FEconf_brush_mul

    subroutine FEconf_elect(FEconf)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconf
        
        !     .. declare local variables

        real(dp) :: exppiA(nsize),exppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s,g,gn         ! dummy indices
        real(dp) :: pro
        real(dp) :: FEconf_local(ngr_node)
        real(dp) :: FEconf_array(size*ngr_node)
        real(dp) :: q_local(ngr_node)
       
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

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            FEconf_local(gn)= 0.0_dp

            g=gn+rank*ngr_node
     
            do c=1,cuantas             ! loop over cuantas
            
                pro=1.0_dp                ! initial weight conformation 
                
                if(weightchain(gn,c)) then ! initial weight conformation 

                    pro=1.0_dp
                
                    do s=1,nseg              ! loop over segments 
                        k=indexchain(s,gn,c)         
                        if(isAmonomer(s)) then ! A segment 
                            pro = pro*exppiA(k)
                        else
                            pro = pro*exppiB(k)
                        endif
                    enddo

                   FEconf_local(gn)=FEconf_local(gn)+pro*log(pro)
                  
                endif
            
            enddo

        enddo   

        ! communicate FEconfABL and FEcopnd_ABR

        if(rank==0) then

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                FEconf_array(g)=FEconf_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(FEconf_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    FEconf_array(g)=FEconf_local(gn)
                enddo
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf+ (FEconf_array(g)/q(g)-log(q(g)))    
            enddo    
        endif

    end subroutine FEconf_elect



    subroutine FEconf_electA(FEconfAB)

        !  .. variables and constant declaractions 

        use globals
        use volume
        use chains
        use field
        use parameters

        implicit none

        real(dp), intent(out) :: FEconfAB
        
        !     .. declare local variables

        real(dp) :: exppiA(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s, kL, g, gn         ! dummy indices
        real(dp) :: proL
        real(dp) :: FEconfABL
        real(dp) :: FEconfABL_local(ngr_node)
        real(dp) :: FEconfABL_array(size*ngr_node)
        real(dp) :: qABL_local(ngr_node)
        real(dp) :: qABL_array(size*ngr_node)
        
        ! .. executable statements 
        ! .. communicate xsol,psi and fdsiA(:,2) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, size-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(fdisA(:,1),nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(fdisA(:,1), nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)          
        endif
            
        do i=1,nsize
              exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable
        enddo
       
        FEconfABL=0.0_dp

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            FEconfABL_local(gn)= 0.0_dp

            g=gn+rank*ngr_node
     
            do c=1,cuantas               ! loop over cuantas
            
                proL=1.0_dp                ! initial weight conformation 
            
                if(weightchain(gn,c)) then ! initial weight conformation 

                    proL=1.0_dp
                
                    do s=1,nseg             ! loop over segments 
                        kL=indexchain(s,gn,c)         
                        proL = proL*exppiA(kL)
                    enddo

                   FEconfABL_local(gn)=FEconfABL_local(gn)+proL*log(proL)
                  
                endif
            
            enddo

        enddo   

        ! communicate FEconfABL and FEcopnd_ABR

        if(rank==0) then

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                FEconfABL_array(g)=FEconfABL_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(FEconfABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    FEconfABL_array(g)=FEconfABL_local(gn)
                enddo
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconfABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconfABL=0.0_dp
            do g=1,ngr 
                FEconfABL = FEconfABL+ (FEconfABL_array(g)/qABL(g)-log(qABL(g)))    
            enddo    
            FEconfAB=FEconfABL
        endif

    end subroutine FEconf_electA


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
        integer  :: i,k,c,s, kL,kR, g,gn         ! dummy indices
        real(dp) :: proL,proR
        real(dp) :: FEconfABL,FEconfABR
        real(dp) :: FEconfABL_local(ngr_node),FEconfABR_local(ngr_node)
        real(dp) :: FEconfABL_array(size*ngr_node),FEconfABR_array(size*ngr_node)
        real(dp) :: qABL_local(ngr_node),qABR_local(ngr_node)
        real(dp) :: qABL_array(size*ngr_node),qABR_array(size*ngr_node)

        ! .. executable statements 

        ! .. communicate xsol,psi and fdsiA(1:) and fdisB(1,:) to other nodes 

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
            
        do i=1,nz
              exppiA(i)=(xsol(i)**vpolA(1))*exp(-zpolA(1)*psi(i))/fdisA(i,1) ! auxiliary variable
              exppiB(i)=(xsol(i)**vpolB(1))*exp(-zpolB(1)*psi(i))/fdisB(i,1) ! auxiliary variable
        enddo
       
        FEconfABL=0.0_dp
        FEconfABR=0.0_dp

        do gn=1,ngr_node              ! loop over grafted points <=>  grafted area on different nodes 
 
            FEconfABL_local(gn)= 0.0_dp
            FEconfABR_local(gn)= 0.0_dp

            g=gn+rank*ngr_node
     
            do c=1,cuantas               ! loop over cuantas
            
                proL=1.0_dp                ! initial weight conformation 
                proR=1.0_dp
            
                if(weightchain(gn,c)) then ! initial weight conformation 

                    proL=1.0_dp
                    proR=1.0_dp

                    do s=1,nseg             ! loop over segments 
                        kL=indexchain(s,gn,c)      
                        kR=mirror_index(kL,nz)    
                        if(isAmonomer(s)) then ! A segment 
                            proL = proL*exppiA(kL)
                            proR = proR*exppiA(kR)
                        else
                            proL = proL*exppiB(kL)
                            proR = proR*exppiB(kR)
                        endif
                    enddo

                   FEconfABL_local(gn)=FEconfABL_local(gn)+proL*log(proL)
                   FEconfABR_local(gn)=FEconfABR_local(gn)+proR*log(proR)
                
                endif
            
            enddo

        enddo   

        ! communicate FEconfABL and FEcopnd_ABR

        if(rank==0) then

            do gn=1,ngr_node
                g = (0)*ngr_node+gn
                FEconfABL_array(g)=FEconfABL_local(gn)
                FEconfABR_array(g)=FEconfABR_local(gn)
            enddo

            do i=1, size-1
                source = i
                call MPI_RECV(FEconfABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(FEconfABR_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(qABL_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(qABR_local, ngr_node, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                
                do gn=1,ngr_node
                    g = (i)*ngr_node+gn
                    FEconfABL_array(g)=FEconfABL_local(gn)
                    FEconfABR_array(g)=FEconfABR_local(gn)
                    qABL_array(g)=qABL_local(gn)
                    qABR_array(g)=qABR_local(gn)
                enddo
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconfABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(FEconfABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABL_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(qABR_local, ngr_node , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            do g=1,ngr
                FEconfABL = FEconfABL+ (FEconfABL_array(g)/qABL_array(g)-log(qABL_array(g)))     ! check sigma
                FEconfABR = FEconfABR+ (FEconfABR_array(g)/qABR_array(g)-log(qABR_array(g)))    
            enddo    

            FEconfAB=FEconfABL+FEconfABR
        endif


    end subroutine FEconf_electdouble
    

end module conform_entropy
