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
            call FEconf_elect(FEconf,Econf)
        case ("neutral")
            call FEconf_neutral(FEconf,Econf)
        case ("neutralnoVdW")
            call FEconf_neutral_noVdW(FEconf,Econf)
        case ("brush_mul","brushdna")
            call FEconf_brush_mul(FEconf,Econf)
        case ("brush_mulnoVdW")
            call FEconf_brush_mulnoVdW(FEconf,Econf)
        case ("brushborn")
            call FEconf_brush_born(FEconf,Econf) 
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
        use chains, only : indexchain, type_of_monomer, logweightchain
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field, only : xsol, rhopol, q, lnproshift
        use parameters, only : vpol, isVdW, VdWscale
        use VdW, only : VdW_contribution_lnexp
        use volume, only : ngr, nset_per_graft

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
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
                lnexppi(i,t) = log(xsol(i))*vpol(t)
            enddo    
        enddo      
       
        if(isVdW) then 
            do t=1,nsegtypes  
                call VdW_contribution_lnexp(rhopol,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init  Econf
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do     
    
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)     
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro=lnpro+lnexppi(k,t)
            enddo  
            pro=exp(lnpro-lnproshift)   
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local=Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local=Rendsqr_local+Rendsqr(c)*pro
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do
        enddo
        
        ! communicate FEconf

        if(rank==0) then

            ! normalize
            FEconf_array=0.0_dp
            Econf_array=0.0_dp  
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp
            As_mtrx_array=0.0_dp

            FEconf_array(1)=FEconf_local
            Econf_array(1)=Econf_local
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do
            
            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)

                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local
                Econf_array(g) =Econf_array(g) +Econf_local
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do
                
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
             ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g)
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)  
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do
            enddo      
        endif

    end subroutine FEconf_neutral


    ! computes conformational entropy in neutral state 
     
    subroutine FEconf_neutral_noVdW(FEconf,Econf) 
    
        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, logweightchain
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field, only : xsol, rhopol, q, lnproshift
        use parameters, only : vpol, isVdW, VdWscale
        use VdW, only : VdW_contribution_exp
        use volume, only : ngr, nset_per_graft

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                !do t=1,nsegtypes
                !    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                !enddo
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            !do t=1,nsegtypes
            !    call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            !enddo
        endif    

        !     .. executable statements 

        do t=1,nsegtypes
            do i=1,nsize
                lnexppi(i,t) = log(xsol(i))*vpol(t)
            enddo    
        enddo      

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do
            
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)        ! internal energy  
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro+ lnexppi(k,t)
            enddo    
            pro=exp(lnpro-lnproshift)
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do
        enddo
        
        ! communicate FEconf

        if(rank==0) then
            ! normalize
            FEconf_array=0.0_dp
            Econf_array=0.0_dp
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp
            As_mtrx_array=0.0_dp  

            FEconf_array(1)=FEconf_local
            Econf_array(1)=Econf_local
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do

            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)            

                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local
                Econf_array(g) =Econf_array(g) +Econf_local
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
             ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g) 
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do 
            enddo    
        endif

    end subroutine FEconf_neutral_noVdW


    subroutine FEconf_brush_mul(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field, only : xsol,psi, fdis,rhopol,q, lnproshift
        use parameters
        use VdW, only : VdW_contribution_lnexp

        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
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
                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i)-log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            else
                do i=1,nsize
                     lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
       
        if(isVdW) then 
            do t=1,nsegtypes  
                call VdW_contribution_lnexp(rhopol,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do
         
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)     
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro+lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)      
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do
        enddo
        
        ! communicate FEconf

        if(rank==0) then

             ! normalize
            FEconf_array=0.0_dp
            Econf_array=0.0_dp  
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp
            As_mtrx_array=0.0_dp

            FEconf_array(1)=FEconf_local
            Econf_array(1)=Econf_local
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do
            
            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                
                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local
                Econf_array(g) =Econf_array(g) +Econf_local 
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local 
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do        
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g)
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do 
            enddo    
        endif

    end subroutine FEconf_brush_mul


    subroutine FEconf_brush_mulnoVdW(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field, only : xsol, psi, fdis, rhopol, q ,lnproshift
        use parameters
        use volume, only : ngr, nset_per_graft
        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k       ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
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
                    lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i)-log(fdis(i,t))   ! auxilary variable palpha
                enddo  
            else
                do i=1,nsize
                     lnexppi(i,t) = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      
       
    
        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do
            
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)       ! internal energy  
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro+ lnexppi(k,t)
            enddo    
            pro=exp(lnpro-lnproshift)
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do
        enddo
        
        ! communicate FEconf

        if(rank==0) then
            ! normalize
            FEconf_array=0.0_dp
            Econf_array=0.0_dp  
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp
            As_mtrx_array=0.0_dp

            FEconf_array(1)=FEconf_local
            Econf_array(1)=Econf_local
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do
            
            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)

                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local
                Econf_array(g) =Econf_array(g) +Econf_local
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local  
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do          
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g)
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do  
            enddo    
           
        endif

    end subroutine FEconf_brush_mulnoVdW


    subroutine FEconf_elect(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use volume, only : ngr, nset_per_graft
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain, isAmonomer
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field,  only : xsol, psi, fdisA,fdisB, rhopol, q ,lnproshift
        use parameters

        real(dp), intent(out) :: FEconf
        real(dp), intent(out) :: Econf
        
        !     .. declare local variables
        real(dp) :: lnexppiA(nsize),lnexppiB(nsize)    ! auxilairy variable for computing P(\alpha) 
        integer  :: i,k,c,s,g,gn         ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: pro,lnpro
        real(dp) :: FEconf_local, Econf_local
        real(dp) :: q_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)
        ! .. executable statements 

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        if(rank==0) then
            do i = 1, numproc-1
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
              lnexppiA(i)=log(xsol(i))*vpolA(1)-zpolA(1)*psi(i)-log(fdisA(i,1)) ! auxiliary variable
              lnexppiB(i)=log(xsol(i))*vpolB(1)-zpolB(1)*psi(i)-log(fdisB(i,1)) ! auxiliary variable
        enddo
       
    
        FEconf_local= 0.0_dp
        Econf_local=0.0_dp 
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do

        do c=1,cuantas             ! loop over cuantas
        
            lnpro=logweightchain(c)               ! initial weight conformation 
            do s=1,nseg              ! loop over segments 
                k=indexchain(s,c)         
                if(isAmonomer(s)) then ! A segment 
                    lnpro = lnpro+lnexppiA(k)
                else
                    lnpro = lnpro+lnexppiB(k)
                endif
            enddo
            pro=exp(lnpro-lnproshift)
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro 
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do     
        enddo  


        ! communicate FEconf

        if(rank==0) then
            ! normalize
            FEconf_array=0.0_dp
            FEconf_array(1)=FEconf_local
            Econf_array=0.0_dp
            Econf_array(1)=FEconf_local
            
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local

            As_mtrx_array=0.0_dp
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do

            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr) 

                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local 
                Econf_array(g)= Econf_array(g)+ Econf_local 
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local 
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do   
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g))) 
                Econf = Econf + Econf_array(g)/q(g)
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do     
            enddo     
        endif

        Econf=0.0_dp

    end subroutine FEconf_elect



    subroutine FEconf_brush_born(FEconf,Econf)

        !  .. variables and constant declaractions 

        use globals, only : nseg, nsegtypes, nsize, cuantas
        use chains, only : indexchain, type_of_monomer, ismonomer_chargeable, logweightchain
        use chains, only : Rgsqr, Rendsqr, As_mtrx, avRgsqr, avRendsqr, avAs_mtrx
        use field, only : xsol,psi, fdis,rhopol,q, lnproshift, fdisA, epsfcn, Depsfcn
        use field, only : xOHmin,xHplus,xNa,xCl,xMg,xCa,xRb
        use parameters, only : bornrad, lb, VdWscale, tA, isrhoselfconsistent, isVdW
        use parameters, only : vpolAA, vsol, vNa, vCl, vRb, vMg, vCa ,vpol
        use parameters, only : zNa, zCl, zRb, zMg, zCa, zpolAA
        use volume, only : ngr, nset_per_graft
        use VdW, only : VdW_contribution_lnexp
        use Poisson, only : Poisson_Equation_Eps, Poisson_Equation_Surface_Eps, grad_pot_sqr_eps_cubic
        use dielectric_const, only : dielectfcn, born
        
        real(dp), intent(out) :: FEconf,Econf
        
        ! .. declare local variables
        real(dp) :: lnexppi(nsize,nsegtypes)          ! auxilairy variable for computing P(\alpha)  
        real(dp) :: pro,lnpro
        integer  :: i,t,g,gn,c,s,k,tc       ! dummy indices
        integer  :: r1,r2,r3,r4,r5,c1,c2,c3,c4,c5  ! dummy indices for indexing As_mtrx
        real(dp) :: FEconf_local
        real(dp) :: Econf_local
        real(dp) :: FEconf_array(ngr)
        real(dp) :: Econf_array(ngr)
        real(dp) :: Rgsqr_local
        real(dp) :: Rendsqr_local
        real(dp) :: Rgsqr_array(ngr)
        real(dp) :: Rendsqr_array(ngr)
        real(dp) :: As_mtrx_local(3,3)
        real(dp) :: As_mtrx_array(ngr,3,3)
        integer  :: tcfdis(3)
        real(dp) :: rhopolAA(nsize),rhopolACa(nsize), rhopolAMg(nsize)
        real(dp) :: lbr,expborn,Etotself,expsqrgrad, Eself
        real(dp) :: expsqrgradpsi(nsize),expEtotself(nsize)

        ! .. communicate xsol,psi and fdsiA(:,1) and fdisB(:,1) to other nodes 

        tcfdis(1)=1
        tcfdis(2)=4
        tcfdis(3)=6
        
        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsize+1 , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(epsfcn, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(Depsfcn, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,3
                    tc=tcfdis(t)
                    call MPI_SEND(fdisA(:,tc) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo    
                do t=1,nsegtypes
                    call MPI_SEND(rhopol(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                enddo
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(psi , nsize+1, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
            call MPI_RECV(epsfcn, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            call MPI_RECV(Depsfcn, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            do t=1,3
                tc=tcfdis(t)
                call MPI_RECV(fdisA(:,tc) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)    
            enddo
            do t=1,nsegtypes
                call MPI_RECV(rhopol(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)  
            enddo
        endif    

        !     .. executable statements 

        ! gradient potential contribution to PDF
        call grad_pot_sqr_eps_cubic(psi,epsfcn, Depsfcn,expsqrgradpsi)

       
        do i=1,nsize
           
            rhopolAA(i) =fdisA(i,1)*rhopol(i,tA)
            rhopolACa(i)=fdisA(i,4)*rhopol(i,tA) 
            rhopolAMg(i)=fdisA(i,6)*rhopol(i,tA) 


            lbr = lb/epsfcn(i)     ! local Bjerrum length

            !xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-born(lbr,bornrad%Na,zNa)-psi(i)*zNa) ! Na+ volume fraction 
            !xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-born(lbr,bornrad%Cl,zCl)-psi(i)*zCl) ! Cl- volume fraction
            !xHplus(i)  = expmu%Hplus*(xsol(i))  *exp(-born(lbr,bornrad%Hplus,1)-psi(i))    ! H+  volume fraction
            !xOHmin(i)  = expmu%OHmin*(xsol(i))  *exp(-born(lbr,bornrad%OHmin,-1)+psi(i))   ! OH- volume fraction
            !xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-born(lbr,bornrad%Rb,zRb)-psi(i)*zRb) ! Rb+ volume fraction
            !xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-born(lbr,bornrad%Ca,zCa)-psi(i)*zCa) ! Ca++ volume fraction 
            !xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-born(lbr,bornrad%Mg,zMg)-psi(i)*zMg) ! Mg++ volume fraction 



            Etotself = &           ! total self energy    
                born(lbr,bornrad%pol  ,zpolAA(1))*rhopolAA(i)  + & ! rhpolAA(i)  = fdisA(i,1)*rhopolin(i,tA)
                born(lbr,bornrad%polCa,zpolAA(4))*rhopolACa(i) + & ! rhopolACa(i)= fdisA(i,4)*rhopolin(i,tA)
                born(lbr,bornrad%polMg,zpolAA(6))*rhopolAMg(i) + & ! rhopolAMg(i)= fdisA(i,6)*rhopolin(i,tA)
                born(lbr,bornrad%Na,zNa)*xNa(i)/(vNa*vsol)     + & 
                born(lbr,bornrad%Cl,zCl)*xCl(i)/(vCl*vsol)     + &
                born(lbr,bornrad%Rb,zRb)*xRb(i)/(vRb*vsol)     + & 
                born(lbr,bornrad%Ca,zCa)*xCa(i)/(vCa*vsol)     + &
                born(lbr,bornrad%Mg,zMg)*xMg(i)/(vMg*vsol)     + &
                born(lbr,bornrad%Hplus,1 )*xHplus(i)/vsol      + &
                born(lbr,bornrad%OHmin,-1)*xOHmin(i)/vsol

            expEtotself(i) = Etotself*(Depsfcn(i)/epsfcn(i))  
        
        enddo 




        do t=1,nsegtypes
      
            if(ismonomer_chargeable(t)) then    
                do i=1,nsize  
                    lbr=lb/epsfcn(i)
                    expborn    = -born(lbr,bornrad%pol,-1)+ expEtotself(i)*vpolAA(1)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpolAA(1)      ! no mutipilcation with vsol because defintion constqE
                    !exppi(i,t) = (xsol(i)**vpolAA(1))*exp(psi(i)+expsqrgrad+expborn) /fdisA(i,1)   ! auxilary variable palpha
                    lnexppi(i,t) = log(xsol(i))*vpolAA(1)+psi(i)+expsqrgrad+expborn -log(fdisA(i,1))   ! auxilary variable palpha  
    
                enddo  
            else

                do i=1,nsize
                    expborn    = expEtotself(i)*vpol(t)*vsol 
                    expsqrgrad = expsqrgradpsi(i)*vpol(t)
                    !exppi(i,t) = (xsol(i)**vpol(t))*exp(expsqrgrad+expborn)
                    lnexppi(i,t) = log(xsol(i))*vpol(t)+expsqrgrad+expborn
                    
                enddo  
            endif   
        enddo      
      
        ! Van der Waals   
        if(isVdW) then 
            do t=1,nsegtypes  
                if(isrhoselfconsistent(t)) call VdW_contribution_lnexp(rhopol,lnexppi(:,t),t)
            enddo
        endif 


        !  .. computation polymer volume fraction      
       
        FEconf_local= 0.0_dp !init FEconf
        Econf_local=0.0_dp ! init FEconf
        Rgsqr_local=0.0_dp ! init Rgsqr
        Rendsqr_local=0.0_dp ! init Rendsqr
        do r1=1,3
            do c1=1,3
                As_mtrx_local(r1,c1)=0.0_dp ! init Asphericity matrix (gyration tensor)
            end do
        end do
         
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c)     
            do s=1,nseg        ! loop over segments                     
                k=indexchain(s,c)
                t=type_of_monomer(s)                
                lnpro = lnpro+lnexppi(k,t)
            enddo 
            pro=exp(lnpro-lnproshift)      
            FEconf_local=FEconf_local+pro*(log(pro)-logweightchain(c))
            Rgsqr_local = Rgsqr_local+Rgsqr(c)*pro
            Rendsqr_local = Rendsqr_local+Rendsqr(c)*pro
            do r2=1,3
                do c2=1,3
                    As_mtrx_local(r2,c2)=As_mtrx_local(r2,c2)+As_mtrx(c,r2,c2)*pro
                end do
            end do
        enddo

        ! communicate FEconf

        if(rank==0) then

             ! normalize
            FEconf_array=0.0_dp
            Econf_array=0.0_dp
            Rgsqr_array=0.0_dp
            Rendsqr_array=0.0_dp 
            As_mtrx_array=0.0_dp 

            FEconf_array(1)=FEconf_local
            Econf_array(1)=Econf_local
            Rgsqr_array(1)=Rgsqr_local
            Rendsqr_array(1)=Rendsqr_local
            do r3=1,3
                do c3=1,3
                    As_mtrx_array(1,r3,c3)=As_mtrx_local(r3,c3)
                end do
            end do
           
            do i=1, numproc-1
                source = i
                call MPI_RECV(FEconf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Econf_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rgsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(Rendsqr_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)
                call MPI_RECV(As_mtrx_local, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)

                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                FEconf_array(g)=FEconf_array(g)+FEconf_local
                Econf_array(g) =Econf_array(g) +Econf_local  
                Rgsqr_array(g) =Rgsqr_array(g) +Rgsqr_local
                Rendsqr_array(g) =Rendsqr_array(g) +Rendsqr_local
                do r4=1,3
                    do c4=1,3
                        As_mtrx_array(g,r4,c4)=As_mtrx_array(g,r4,c4)+As_mtrx_local(r4,c4)
                    end do
                end do           
            enddo 
        else     ! Export results
            dest = 0
            call MPI_SEND(FEconf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Econf_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rgsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(Rendsqr_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(As_mtrx_local, 1 , MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
        endif


        if(rank==0) then
            ! normalize
            FEconf=0.0_dp
            Econf=0.0_dp
            do g=1,ngr 
                FEconf = FEconf + (FEconf_array(g)/q(g)-log(q(g)))  
                Econf = Econf + Econf_array(g)/q(g)
                avRgsqr(g) = Rgsqr_array(g)/q(g)
                avRendsqr(g) = Rendsqr_array(g)/q(g)
                do r5=1,3
                    do c5=1,3
                        avAs_mtrx(g,r5,c5) = As_mtrx_array(g,r5,c5)/q(g)
                    end do
                end do  
            enddo    
        endif

    end subroutine FEconf_brush_born


end module conform_entropy
