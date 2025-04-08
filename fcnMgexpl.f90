! --------------------------------------------------------------|
! fcnMg-excpl.f90:                                              |
! Constructs the vector function  needed by the                 |
! routine solver, which solves the SCMFT eqs for                |
! polyelctolyte  with ionbinding and  Mg binding                |
! using phosphate pairs.                                        |
! --------------------------------------------------------------|

module modfcnMgexpl
   
    use mpivars
    use precision_definition

    implicit none


contains
   
    subroutine compute_fdisPP(fdisPP,fdisP2Mg,position1,position2)

        use field, only : xHplus, xNa, xK, xMg, xsol
        use parameters, only : vNa, vK, vMg, deltavAA
        use parameters, only : K0aAA
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 

        real(dp), intent(inout), dimension(:,:) :: fdisPP
        real(dp), intent(inout) :: fdisP2Mg
        integer , intent(in) :: position1, position2

        real(dp) :: xP(5,2), xP2Mg, fPP, sumxP         ! disociation variables  
        integer :: i, j     
        integer  :: JJ, KK

        ! .. executable statements 
                                                     
        xP(Phos,1)   = 1.0_dp 
        xP(Phos,2)   = 1.0_dp    
      
        i = position1 ! position in lattice numbers
        j = position2

        !print*,"i=",i,' j=',j

        xP(PhosH,1)  = xHplus(i)/(K0aAA(1)*(xsol(i)**deltavAA(1)))      !  (PH)/P-    : f(PH)P(i,j)/fPP(i,j)
        xP(PhosH,2)  = xHplus(j)/(K0aAA(1)*(xsol(j)**deltavAA(1)))      !  (PH)/P-    : fP(PH)(i,j)/fPP(i,j)
      
        xP(PhosNa,1) = (xNa(i)/vNa)/(K0aAA(2)*(xsol(i)**deltavAA(2)))   !  PNa/P-     : f(PNa)P(i,j)/fPP(i,j) 
        xP(PhosNa,2) = (xNa(j)/vNa)/(K0aAA(2)*(xsol(j)**deltavAA(2)))   !  PNa/P-     : fPPNa(i,j)/fPP(i,j) 
      
        xP(PhosK,1)  = (xK(i)/vK)/(K0aAA(7)*(xsol(i)**deltavAA(7)))     !  PK/P-      : f(PK)P(i,j)/fPP(i,j) 
        xP(PhosK,2)  = (xK(j)/vK)/(K0aAA(7)*(xsol(j)**deltavAA(7)))     !  PK/P-      : fPPK(i,j)/fPP(i,j) 
     
        xP(PhosMg,1) = (xMg(i)/vMg)/(K0aAA(5)*(xsol(i)**deltavAA(5)))   ! PMgP+/PP2-
        xP(PhosMg,2) = (xMg(j)/vMg)/(K0aAA(5)*(xsol(j)**deltavAA(5)))   ! PPAMg+/PP2-

        xP2Mg  = sqrt( (xMg(i)/vMg)*(xMg(j)/vMg)/ ((K0aAA(6)**2) *(xsol(i)**deltavAA(6))*(xsol(j)**deltavAA(6)))) ! P2Mg/PP2- 
            ! deltavAA(6) = 2.0_dp*vpolAA(1)+vMg-vpolAA(7) ! 2vA- + vMg2+ -vA2Mg == "total volume vP2Mg thus divison by 2 aka sqrt !!
       
        sumxP = 0.0_dp
        do JJ=1,5
            do KK=1,5
                sumxP = sumxP + xP(JJ,1) * xP(KK,2)
            enddo
        enddo
        sumxP=sumxP+xP2Mg

        fPP = 1.0_dp/sumxP    ! fraction of phophate pairs that are both charged
          
        do JJ=1,5             ! fraction of phophate pairs that form a bind with H^+,Na^+,K^+
             do KK=1,5
                 fdisPP(JJ,KK) = fPP * xP(JJ,1) * xP(KK,2)
             enddo
        enddo
            
        fdisP2Mg = fPP * xP2Mg  ! fraction of phophate pairs that form a Mg-bridge
                        
    end subroutine compute_fdisPP


    ! polyelectrolyte with phosphate acid groups
    ! with ion chargeable group being on one acid (tA) with counterion binding etc 
    ! distribute volume of neighboring cells

    subroutine fcn_Mg_expl(x,f,nn)

        use mpivars
        use globals, only    : nsize, nsegtypes, nseg, neq, neqint, cuantas, bcflag
        use parameters, only : expmu 
        use parameters, only : vsol, vpol, vNa, vK, vCl, vRb, vCa, vMg, vpro, vPP ! , vpolAA ! ,deltavAA,vnucl,vPP
                                                                          ! check vpolAA and etc 
        use parameters, only : zpol,zNa,zK,zCl,zRb,zCa,zMg,qPP, K0a  ! ,K0aAA,K0a !, K0aion
        use parameters, only : ta, isVdW, iter ! isrhoselfconsistent
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 
        use volume, only     : volcell, nset_per_graft 
        use chains, only     : indexchain, type_of_monomer, logweightchain, ismonomer_chargeable
        use chains, only     : indexconfpair, nneigh
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xRb,xMg,xCa,rhopol,rhoqpol,rhoq
        use field, only      : xpro
        use field, only      : psi, fdis, rhopol_charge
        use field, only      : fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap, rhoqphos
        use field, only      : q, lnproshift
        use field, only      : xpol=>xpol_t, xpol_tot=>xpol
        use vectornorm, only : L2norm, L2norm_f90
        use VdW, only        : VdW_contribution_lnexp
        use surface, only    : LEFT, sigmaqSurfL, psiSurfL, RIGHT, sigmaqSurfR, psiSurfR, surface_charge
        use Poisson, only    : Poisson_Equation, Poisson_Equation_Surface

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)                     ! local density nucleosome
        real(dp) :: local_xpolphos(nsize)                             ! local volume  fraction of phophates  
        real(dp) :: local_rhoqphos(nsize)                             ! local charge density of phosphates     
        real(dp) :: local_q                                           ! local normalization q     
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize) 
        real(dp) :: pro,lnpro
        integer  :: n, i, j, k, c, s, t , m, g        ! dummy indices
        integer  :: JJ, KK 
        real(dp) :: norm, normvol, normPE
        real(dp) :: rhopol0 
        real(dp) :: locallnproshift(2), globallnproshift(2)
        real(dp) :: sum_rhoqphos,sum_xphos

        !     .. executable statements 
        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        n=nsize
        ! read out x 
        k=n
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            psi(i)  = x(i+k)      ! potential
        enddo  
             
        !  .. assign global and local polymer density 
        do t=1,nsegtypes
            do i=1,n
                xpol(i,t)  = 0.0_dp 
                rhopol(i,t) = 0.0_dp 
                local_rhopol(i,t) = 0.0_dp
                rhopol_charge(i,t) = 0.0_dp
            enddo    
        enddo    
       
        do i=1,n
            local_xpolphos(i)  = 0.0_dp 
            local_rhoqphos(i) = 0.0_dp 
        enddo

        do i=1,n     ! init volume fractions
            xpol_tot(i) = 0.0_dp                                  ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xpro(i)    = expmu%pro*(xsol(i)**vpro)                ! crowder volume fraction  

            lnexppivw(i) = log(xsol(i))/vsol                      ! auxilary variable  divide by vsol  !!
            local_rhoqphos(i) = 0.0_dp 
        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    ! charged not phosphate 
                    do i=1,n                                         
                        fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                        lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                    enddo  
                else
                    ! t=ta : phosphate 
                    do i=1,nsize  
                        lnexppi(i,t) = psi(i)         ! auxilary variable palpha
                    enddo
                endif
            else    
                ! neutral monomomer
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      

        if(isVdW) then 
           ! Van der Waals
            print*,"isVdW true for fcn_Mg_expl, stop!!"
            stop 
        endif 

        !  .. computation polymer volume fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                t=type_of_monomer(s)
                if(t/=ta) then 
                    k=indexchain(s,c)                
                    lnpro = lnpro +lnexppi(k,t)
                else 
                    ! phosphates 
                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                        call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        lnpro =lnpro + (lnexppi(k,ta) +lnexppi(m,ta)+(lnexppivw(k)+lnexppivw(m))*(vpol(tA)*vsol) &
                                    -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                    enddo
                endif           
            enddo 

        enddo
 
        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                t=type_of_monomer(s)
                if(t/=ta) then 
                    k=indexchain(s,c)                
                    lnpro = lnpro +lnexppi(k,t)
                else 
                    ! phosphates 
                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                        call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) &
                                          -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))                     
                    enddo
                endif           
            enddo 

            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro
            do s=1,nseg
                t = type_of_monomer(s)
                if(t/=ta) then 
                    k = indexchain(s,c) 
                    local_rhopol(k,t) = local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
                else
                 
                    ! pair density of phosphates 
                    k = indexchain(s,c)

                    do j=1,nneigh(s,c)

                        m = indexconfpair(s,c)%elem(j)
                             
                        call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)
                        call  compute_fdisPP(fdisPP_loc_swap, fdisP2Mg_loc_swap,  m , k) 
                        
                        ! first part integral   

                        sum_rhoqphos = 0.0_dp
                        sum_xphos = 0.0_dp 
                        
                        do JJ=1,5
                            do KK=1,5
                                sum_rhoqphos = sum_rhoqphos+&
                                    (fdisPP_loc(JJ,KK)*qPP(JJ)+fdisPP_loc_swap(JJ,KK)*qPP(KK))/2.0_dp
                                sum_xphos = sum_xphos   +&
                                    (fdisPP_loc(JJ,KK)*vPP(JJ)+fdisPP_loc_swap(JJ,KK)*vPP(KK))/2.0_dp
                            enddo
                        enddo
        
                        sum_xphos = sum_xphos+(fdisP2Mg_loc+fdisP2Mg_loc_swap)*vPP(Phos2Mg)/4.0_dp 
                            ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                      
                        local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                        local_xpolphos(k) = local_xpolphos(k) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                        local_rhopol(k,ta)=local_rhopol(k,ta)+pro/(2.0_dp*nneigh(s,c))
                            
                        ! second integral contributes to location m of rhoqpos and xphol  xpol  
                         
                        sum_rhoqphos=0.0_dp
                        sum_xphos=0.0_dp 
                        
                        ! contributes to location m of rhoqpos and xol
                               
                        do JJ=1,5
                            do KK=1,5   
                                sum_rhoqphos = sum_rhoqphos+&
                                    (fdisPP_loc_swap(JJ,KK)*qPP(JJ)+fdisPP_loc(JJ,KK)*qPP(KK))/2.0_dp
                                sum_xphos = sum_xphos   +&
                                    (fdisPP_loc_swap(JJ,KK)*vPP(JJ)+fdisPP_loc(JJ,KK)*vPP(KK))/2.0_dp
                            enddo
                        enddo
        
                        sum_xphos=sum_xphos+(fdisP2Mg_loc_swap +fdisP2Mg_loc)*vPP(Phos2Mg)/4.0_dp

                        ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                      
                        local_rhoqphos(m) = local_rhoqphos(m) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                        local_xpolphos(m) = local_xpolphos(m) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                        local_rhopol(m,ta)=local_rhopol(m,ta)+pro/(2.0_dp*nneigh(s,c))

                    enddo 

                endif 
            
            enddo
        enddo

        !  .. import results 

        ! print*,"rank=",rank, " local_q", local_q, " lnproshift=", lnproshift, "locallnproshift=", locallnproshift

        if (rank==0) then 

            q=0.0_dp
            g=1  
            q(g)=local_q
            
             do i=1, numproc-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                q(g)=q(g)+local_q
            enddo

            ! first graft point/processor 
            do t=1,nsegtypes
                do k=1,n
                    rhopol(k,t) = local_rhopol(k,t)/q(1) ! polymer density 
                enddo
            enddo

            do k=1,nsize
                xpol(k,tA) = local_xpolphos(k)/q(1)
                rhoqphos(k) = local_rhoqphos(k)/q(1) 
            enddo

            ! other graft points/processors  
            do i=1, numproc-1
                source = i
                g =int(source/nset_per_graft)+1 
               
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t)/q(g) ! polymer density 
                    enddo    
                enddo

                call MPI_RECV(local_xpolphos, nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)     
                call MPI_RECV(local_rhoqphos, nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                
                do k=1,nsize
                    xpol(k,tA) = xpol(k,tA) + local_xpolphos(k)/q(g)      ! phosphate volume fraction
                    rhoqphos(k) = rhoqphos(k) + local_rhoqphos(k)/q(g)
                enddo 
            enddo  
        
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)! volume polymer segment per volume cell

            do t=1, nsegtypes

                if(ismonomer_chargeable(t)) then 

                    if(t/=ta) then
                                                              
                        do i=1,n
                            rhopol(i,t) = rhopol0 * rhopol(i,t)         ! density polymer of type t  
                            xpol(i,t)   = rhopol(i,t) * vpol(t) * vsol  ! volume fraction polymer
                            rhoqpol(i)  = rhoqpol(i) + &
                                    (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol ! total  charge density in units of vsol   
                        enddo  

                    else
                        ! phophate t=tA
                         do i=1,n
                            rhopol(i,ta) = rhopol0 * rhopol(i,ta) 
                            rhoqphos(i)  = rhopol0 * rhoqphos(i) 
                            xpol(i,ta)   = rhopol0 * xpol(i,tA)
                            rhoqpol(i)   = rhoqpol(i) + rhoqphos(i) * vsol ! total  charge density in units of vsol 
                        enddo           
                    endif    
                else  
                     ! neutral monomeer
                    do i=1,n
                        rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                        xpol(i,t)    = xpol(i,t) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                    enddo
                endif   

                do i=1,n
                    xpol_tot(i) = xpol_tot(i)+xpol(i,t)  
                enddo

            enddo    

            do i=1,n
                f(i) = xpol_tot(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i)+xpro(i) -1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK ! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 

            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)
            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
    
            ! .. boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)    

           

            norm=l2norm_f90(f)
            iter=iter+1
                        
            normvol = L2norm_f90(f(1:nsize))
            normPE  = L2norm_f90(f(nsize+1:2*nsize))
                        
            print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE

        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo
           
            call MPI_SEND(local_xpolphos, nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)      
            call MPI_SEND(local_rhoqphos, nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
    
        endif

    end subroutine fcn_Mg_expl



    ! compute the average fraction of charged state of the phosphate pairs 

    subroutine compute_average_charge_PP_expl(avfdisP2Mg,avfdisPP)

    
        !     .. local variables
        use mpivars
        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, cuantas, DEBUG
        use parameters, only : vsol, vpol
        use parameters, only : zpol, Phos
        use parameters, only : ta !, isVdW! isrhoselfconsistent 
        use volume, only     : nx, ny, ngr
        use volume, only     : nset_per_graft
        use chains, only     : indexchain, type_of_monomer, logweightchain, ismonomer_chargeable
        use chains, only     : indexconfpair, nneigh
        use field, only      : xsol,psi,fdis, rhopol_charge, fdisPP_loc,  fdisP2Mg_loc ,fdisPP_loc_swap , fdisP2Mg_loc_swap
        use field, only      : q, lnproshift
        use myutils, only    : error_handler

        real(dp), intent(inout) :: avfdisP2Mg
        real(dp), intent(inout) :: avfdisPP(5,5)

        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,c,s,m,t,g                 ! dummy indices
        integer  :: JJ, KK
        real(dp) :: local_avfdisP2Mg,local_avfdisPP(5,5)
        real(dp) :: sumrhopairs
        integer  :: nsizepsi

        ! .. executable statements 

        ! .. communication between processors 

       
        nsizepsi = nsize + 2 * nx * ny
        local_avfdisPP = 0.0_dp
        local_avfdisP2Mg = 0.0_dp

        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsizepsi , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(rhopol_charge(:,ta) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then 
                        call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    endif
                enddo
                call MPI_SEND(q , ngr , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(psi , nsizepsi, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(rhopol_charge(:,ta) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then 
                    call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                endif
            enddo

            call MPI_RECV(q , ngr, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    

        n=nsize

        do i=1,nsize
            lnexppivw(i)=log(xsol(i))/vsol    ! auxilary variable for t=tA  
        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    ! charged not phosphate 
                    do i=1,n                                         
                        lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                    enddo   
                else
                    ! t=ta : phosphate
                    do i=1,n  
                        lnexppi(i,t) =  psi(i)!!   ! auxilary variable palpha
                    enddo

                endif
            else  
                !lnexppi(:,t) = 0.0_dp
                lnexppi(i,t) = log(xsol(i))*vpol(t)
            endif   
        enddo   


        !  .. computation of probability 

        lnpro = 0.0_dp
              
       
         do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                t=type_of_monomer(s)
                if(t/=ta) then ! not phosphate either charged or neutral
                    k=indexchain(s,c)                
                    lnpro = lnpro +lnexppi(k,t)
                else 
                    ! phosphates 
                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                        call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) &
                            -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                        ! lnpro =lnpro +  (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) /(2.0_dp*nneigh(s,c))    

                    enddo
                endif           
            enddo 

            pro=exp(lnpro-lnproshift)   
        
            do s=1,nseg
                    
                t=type_of_monomer(s)

                if(t==ta) then 
                                    
                    ! pair density of phosphates 
                    k = indexchain(s,c)
       
                    do j=1,nneigh(s,c)

                        m = indexconfpair(s,c)%elem(j)

                        call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k ,m)
                        call compute_fdisPP(fdisPP_loc_swap,fdisP2Mg_loc_swap, m ,k)

                        do JJ=1,5
                            do KK=1,5
                                !local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                !    fdisPP_loc(JJ,KK)*pro/nneigh(s,c)

                                local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                    (fdisPP_loc(JJ,KK)+fdisPP_loc_swap(JJ,KK))*pro/(2.0_dp*nneigh(s,c)) 
                                
                                !local_avfdisPP(JJ,KK) = local_avfdisPP(JJ,KK)+&
                                !    (fdisPP_loc(JJ,KK))*pro/(2.0_dp*nneigh(s,c))
                            
                            enddo
                        enddo
                            
                        local_avfdisP2Mg=local_avfdisP2Mg+fdisP2Mg_loc*pro/nneigh(s,c)
        
                    enddo
                endif
            enddo   
        enddo

       
        !   .. import results 

        if (rank==0) then 
        
            avfdisP2Mg = local_avfdisP2Mg/q(1)
            do i=1, numproc-1
                source = i
                g =int(source/nset_per_graft)+1
                call MPI_RECV(local_avfdisP2Mg, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                  avfdisP2Mg = avfdisP2Mg + local_avfdisP2Mg/q(g)
            enddo

            avfdisPP = local_avfdisPP/q(1)
            do i=1, numproc-1
                source = i
                g =int(source/nset_per_graft)+1
                call MPI_RECV(local_avfdisPP, 25, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                avfdisPP = avfdisPP + local_avfdisPP/q(g)
            enddo

            ! .. construction of avfdisP2Mg and avfdisPP 
            ! .. normalized avfdisPP with number of average number pairs = integral of rhopol_charge(:,ta) in Nucleosome prog . 

            ! sumrhopairs=sum(rhopol_charge(:,tA)) 
            ! sumrhopairs=sumrhopairs*volcell

            ! .. here rho_charge not computed 
            ! .. alternative computations
            sumrhopairs=0.0_dp
            do s=1,nseg
                t=type_of_monomer(s)
                if(t==ta) then 
                    sumrhopairs = sumrhopairs + 1
                endif 
            enddo  
            sumrhopairs = sumrhopairs * ngr  ! /2.0_dp   

            !sumrhopairs=sum(rhopol_charge(:,tA)) 
            !sumrhopairs=sumrhopairs*volcell

            avfdisPP=avfdisPP/(sumrhopairs)  !*q) ! also norm with q
           ! avfdisP2Mg=avfdisP2Mg/(sumrhopairs*q)
            avfdisP2Mg=avfdisP2Mg/(sumrhopairs)
    
        
        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_avfdisP2Mg, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            call MPI_SEND(local_avfdisPP,25, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

        endif

    end subroutine compute_average_charge_PP_expl



    ! compute FEchem contribution in case of binding of Mg tio pairs of phosphate   
    !
    ! Note FEchem does implicity involve volcell= delta^3 ( volume integration )
    ! instead we evalue delta function 
    ! becuase \int \int dV dv" <\rho(r,r')> f_IJ(r,r') lne  = \sum_a Pro(a) delta (r-r(a,s)) delta(r'-r(a,t)) 
    ! Compare density instead pair density contribution \int dV \rho(x) f(x) ln(f(x) + ..) 
    ! this invlove a explicet sum of  hte volume  \int dV f(v) = delta^3 \sum(i,j,k) f(i,j,k)

    subroutine compute_FEchem_react_PP_expl(FEchemPP)

        use precision_definition
        use globals, only    : nsize, nsegtypes, nseg, cuantas, DEBUG
        use parameters, only : vsol, vpol, zpol
        use parameters, only : vPP, qPP,  Phos, Phos2Mg, ta 
        use volume, only     : nx, ny, ngr
        use volume, only     : nset_per_graft
        use chains, only     : indexchain, type_of_monomer, logweightchain, ismonomer_chargeable
        use chains, only     : indexconfpair, nneigh
        use field, only      : xsol, psi, fdis
        use field, only      : fdisPP_loc, fdisP2Mg_loc! ,fdisPP_loc_swap, fdisP2Mg_loc_swap
        use field, only      : q, lnproshift
        use myutils, only    : error_handler

        real(dp), intent(inout) :: FEchemPP

        !     .. local variables
        
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize)
        real(dp) :: pro,lnpro
        integer  :: n,i,j,k,c,s,m,t,g                ! dummy indices
        integer  :: JJ, KK
        real(dp) :: local_FEchempair, FEchempair
        !real(dp) :: K0aPP   ! Kdis of P2Mg pair temporarily define 
        integer  :: nsizepsi
        real(dp) :: betapi_k, betapi_m, psi_k, psi_m, lambda, sum_pi, sum_psi


        ! .. executable statements 

        ! .. communication between processors 

        !K0aPP = K0aAA(6) ! P2Mg
        nsizepsi = nsize + 2 * nx * ny
        local_FEchempair = 0.0_dp
       
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 

        if(rank==0) then
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(xsol, nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(psi , nsizepsi , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                do t=1,nsegtypes
                    if(ismonomer_chargeable(t)) then 
                        call MPI_SEND(fdis(:,t) , nsize , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
                    endif
                enddo
                call MPI_SEND(q , ngr , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        else
            source = 0 
            call MPI_RECV(xsol, nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            call MPI_RECV(psi , nsizepsi, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr)   
            do t=1,nsegtypes
                if(ismonomer_chargeable(t)) then 
                    call MPI_RECV(fdis(:,t) , nsize, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
                endif
            enddo

            call MPI_RECV(q , ngr, MPI_DOUBLE_PRECISION, source,tag, MPI_COMM_WORLD,stat, ierr) 
        endif    
        
        n=nsize

        do i=1,nsize
            lnexppivw(i)=log(xsol(i))/vsol
        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    ! charged not phosphate 
                    do i=1,n                                         
                        lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                    enddo   
                else
                    ! t=ta : phosphate
                    do i=1,n  
                        lnexppi(i,t) = psi(i)!!   ! auxilary variable palpha
                    enddo

                endif
            else  
                lnexppi(:,t) =  log(xsol(i))*vpol(t)
            endif   
        enddo   


        !  .. computation of probability 

        lnpro = 0.0_dp
              
        do c=1,cuantas         ! loop over cuantas

            lnpro = logweightchain(c) 
        
            do s=1,nseg        ! loop over segments 
        
                t = type_of_monomer(s)

                if(t/=ta) then 
                    k = indexchain(s,c)                
                    lnpro = lnpro + lnexppi(k,t)
                else 
                    ! phosphates 

                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                        call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        lnpro = lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) &
                             -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                            
                    enddo
                endif           
            enddo 

            pro=exp(lnpro-lnproshift)   
        
            do s=1,nseg
                    
                t=type_of_monomer(s)

                if(t==ta) then 
                                
                    ! pair density of phosphates 
                    k = indexchain(s,c)
                
                    betapi_k=-log(xsol(k))/vsol
                    psi_k = psi(k)
    
                    do j=1,nneigh(s,c)

                        m = indexconfpair(s,c)%elem(j)

                        betapi_m= -log(xsol(m))/vsol
                        psi_m = psi(m)
                    
                        call compute_fdisPP(fdisPP_loc,fdisP2Mg_loc, k, m)

                        ! Lagrange multiplier lambda(r,r') 

                        lambda = -(betapi_k +betapi_m)*vPP(Phos) -(psi_k+psi_m)*qPP(Phos) &
                            -log(fdisPP_loc(Phos,Phos))

    
                        lambda = lambda*pro/nneigh(s,c)        

                        sum_pi  = 0.0_dp
                        sum_psi = 0.0_dp

                        do JJ=1,5
                            do KK=1,5
                                sum_pi=sum_pi-(vPP(JJ)*betapi_k+vPP(KK)*betapi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                                sum_psi=sum_psi-(qPP(JJ)*psi_k+qPP(KK)*psi_m)*fdisPP_loc(JJ,KK)*pro/nneigh(s,c)
                            enddo
                        enddo
                       ! print*,"sum_psi=",sum_psi, " sum_pi=",sum_pi

                        sum_pi=sum_pi-(vPP(Phos2Mg)/2.0_dp)*(betapi_k+betapi_m)*fdisP2Mg_loc*pro/nneigh(s,c)

                        ! division 2.0_dp  because  vPP(Phos2Mg)/2 is volume change per phosphate 
                        
                        local_FEchempair = local_FEchempair+(-lambda +sum_pi+sum_psi)/2.0_dp             
                    
                    enddo 

                endif

            enddo

        enddo

        !   .. import results 

        if (rank==0) then 
            g=1
            FEchempair = local_FEchempair/q(g)

            do i=1, numproc-1
                source = i
                g = int(source/nset_per_graft)+1
                call MPI_RECV(local_FEchempair, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                FEchempair = FEchempair + local_FEchempair/q(g)
            enddo

            ! .. normalized FEchempair  with q 
            !  FEchempair = FEchempair /q 

            FEchemPP=FEchempair 

        else                      ! Export results 
            
            dest = 0 

            call MPI_SEND(local_FEchempair, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            
        endif

    end subroutine compute_FEchem_react_PP_expl






    ! Neutral polymer  using pairs of neutral  phosphate groups
    ! neutral version of fcn_Mg_exp

    subroutine fcn_neutral_expl(x,f,nn)

        use mpivars
        use globals, only    : nsize, nsegtypes, nseg, neq, neqint, cuantas, bcflag
        use parameters, only : expmu 
        use parameters, only : vsol, vpol, vNa, vK, vCl, vRb, vCa, vMg, vpro, vPP ! , vpolAA ! ,deltavAA,vnucl,vPP
                                                                          ! check vpolAA and etc 
        use parameters, only : zpol,zNa,zK,zCl,zRb,zCa,zMg,qPP, K0a  ! ,K0aAA,K0a !, K0aion
        use parameters, only : ta, isVdW, iter ! isrhoselfconsistent
        use parameters, only : Phos, PhosH, PhosK, PhosNa, PhosMg, Phos2Mg 
        use volume, only     : volcell, nset_per_graft 
        use chains, only     : indexchain, type_of_monomer, logweightchain, ismonomer_chargeable
        use chains, only     : indexconfpair, nneigh
        use field, only      : xsol,xNa,xCl,xK,xHplus,xOHmin,xRb,xMg,xCa,rhopol,rhopolin,rhoqpol,rhoq
        use field, only      : xpro
        use field, only      : psi, fdis, rhopol_charge
        use field, only      : fdisPP_loc, fdisPP_loc_swap, fdisP2Mg_loc, fdisP2Mg_loc_swap, rhoqphos
        use field, only      : q, lnproshift
        use field, only      : xpol=>xpol_t, xpol_tot=>xpol
        use vectornorm, only : L2norm, L2norm_f90
        use VdW, only        : VdW_contribution_lnexp
        use surface, only    : LEFT, sigmaqSurfL, psiSurfL, RIGHT, sigmaqSurfR, psiSurfR, surface_charge
        use Poisson, only    : Poisson_Equation, Poisson_Equation_Surface

        !     .. scalar arguments

        integer(8), intent(in) :: nn

        !     .. array arguments

        real(dp), intent(in) :: x(neq)
        real(dp), intent(out) :: f(neq)

        !     .. local variables
        
        real(dp) :: local_rhopol(nsize,nsegtypes)                     ! local density nucleosome
        ! real(dp) :: local_xpol(nsize,nsegtypes)                       ! local volume  fraction polymer
        real(dp) :: local_xpolphos(nsize)                             ! local volume  fraction of phophates  
        real(dp) :: local_rhoqphos(nsize)                             ! local charge density of phosphates     
        real(dp) :: local_q                                           ! local normalization q     
        real(dp) :: lnexppi(nsize,nsegtypes)                          ! auxilairy variable for computing P(\alpha) 
        real(dp) :: lnexppivw(nsize) 
        real(dp) :: pro,lnpro
        integer  :: n, i, j, k, c, s, t , m, g        ! dummy indices
        integer  :: JJ, KK
        
        real(dp) :: norm, normPE, normvol
        real(dp) :: rhopol0 
        real(dp) :: locallnproshift(2), globallnproshift(2)
        real(dp) :: sum_rhoqphos,sum_xphos

        !     .. executable statements 
        !     .. communication between processors 

        if (rank.eq.0) then 
            flag_solver = 1      !  continue program  
            do i = 1, numproc-1
                dest = i
                call MPI_SEND(flag_solver, 1, MPI_INTEGER,dest, tag,MPI_COMM_WORLD,ierr)
                call MPI_SEND(x, neqint , MPI_DOUBLE_PRECISION, dest, tag,MPI_COMM_WORLD,ierr)
            enddo
        endif

        n=nsize
        ! read out x 
        k=n
        do i=1,n                     
            xsol(i) = x(i)        ! volume fraction solvent
            psi(i)  = x(i+k)      ! potential
        enddo  
     
       ! do t=1,nsegtypes
       !     k=(t+1)*n
       !     do i=1,n 
       !         rhopolin(i,t) = x(i+k) ! density 
       !     enddo    
       ! enddo
             
        !  .. assign global and local polymer density 
        do t=1,nsegtypes
            do i=1,n
                xpol(i,t)  = 0.0_dp 
                rhopol(i,t) = 0.0_dp 
                local_rhopol(i,t) = 0.0_dp
                rhopol_charge(i,t) = 0.0_dp
            enddo    
        enddo    
       
        do i=1,n
            local_xpolphos(i)  = 0.0_dp 
            local_rhoqphos(i) = 0.0_dp 
        enddo

        ! set fdisPP to neutral only 
        ! used to ease of evaluation of sumxphos see below
        do JJ=1,5
            do KK=1,5
                fdisPP_loc(JJ,KK)=0.0_dp
                fdisPP_loc_swap(JJ,KK)=0.0_dp
            enddo
        enddo 
           
        fdisP2Mg_loc =0.0_dp 
        fdisP2Mg_loc_swap  =0.0_dp 
        fdisPP_loc(PhosH,PhosH)=1.0_dp
        fdisPP_loc_swap(PhosH,PhosH)=1.0_dp

        do i=1,n                  ! init volume fractions
            xpol_tot(i) = 0.0_dp    ! volume fraction polymer
            rhoqpol(i) = 0.0_dp                                   ! charge density AA monomoer
            xNa(i)     = expmu%Na*(xsol(i)**vNa)*exp(-psi(i)*zNa) ! Na+ volume fraction
            xK(i)      = expmu%K*(xsol(i)**vK)*exp(-psi(i)*zK)    ! K+ volume fraction
            xCl(i)     = expmu%Cl*(xsol(i)**vCl)*exp(-psi(i)*zCl) ! Cl- volume fraction
            xHplus(i)  = expmu%Hplus*(xsol(i))*exp(-psi(i))       ! H+  volume fraction
            xOHmin(i)  = expmu%OHmin*(xsol(i))*exp(+psi(i))       ! OH- volume fraction
            xRb(i)     = expmu%Rb*(xsol(i)**vRb)*exp(-psi(i)*zRb) ! Rb+ volume fraction
            xCa(i)     = expmu%Ca*(xsol(i)**vCa)*exp(-psi(i)*zCa) ! Ca++ volume fraction
            xMg(i)     = expmu%Mg*(xsol(i)**vMg)*exp(-psi(i)*zMg) ! Mg++ volume fraction
            xpro(i)    = expmu%pro*(xsol(i)**vpro)                ! crowder volume fraction  

            lnexppivw(i) = log(xsol(i))/vsol                       ! auxilary variable  divide by vsol  !!
            local_rhoqphos(i) = 0.0_dp 

        enddo

        do t=1,nsegtypes
            if(ismonomer_chargeable(t)) then
                if(t/=ta) then
                    ! charged not phosphate 
                    do i=1,n                                         
                        fdis(i,t)  = 1.0_dp/(1.0_dp+xHplus(i)/(K0a(t)*xsol(i)))      
                        lnexppi(i,t) = log(xsol(i))*vpol(t) -zpol(t,2)*psi(i) -log(fdis(i,t))   ! auxilary variable palpha
                    enddo  
                else
                   ! t=ta : phosphate neutral !!
                   
                endif
            else    
                ! neutral  monomomer
                do i=1,n
                    fdis(i,t)  = 0.0_dp
                    lnexppi(i,t)  = log(xsol(i))*vpol(t)
                enddo  
            endif   
        enddo      

        ! Van der Waals   
        if(isVdW) then 
            print*,"isVdW true for fcn_neutral_expl, stop!!"
            stop 
            do t=1,nsegtypes  
                call VdW_contribution_lnexp(rhopolin,lnexppi(:,t),t)
            enddo
        endif 

        !  .. computation polymer volume fraction      
 
        local_q = 0.0_dp    ! init q
        lnpro = 0.0_dp
        
        do c=1,cuantas         ! loop over cuantas

            lnpro=lnpro+logweightchain(c)        ! internal weight

            do s=1,nseg        ! loop over segments 
                t=type_of_monomer(s)
                if(t/=ta) then 
                    k=indexchain(s,c)                
                    lnpro = lnpro +lnexppi(k,t)
                else 
                    ! phosphates 
                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                       ! call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        !lnpro =lnpro + (lnexppi(k,ta) +lnexppi(m,ta)+(lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) &
                        !                  -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))    
                        lnpro =lnpro +  (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) /(2.0_dp*nneigh(s,c))    

                    enddo
                endif           
            enddo 

        enddo
 
        locallnproshift(1)=lnpro/cuantas
        locallnproshift(2)=rank  
    
        call MPI_Barrier(  MPI_COMM_WORLD, ierr) ! synchronize 
        call MPI_ALLREDUCE(locallnproshift, globallnproshift, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_WORLD,ierr)
       
        lnproshift=globallnproshift(1)
             
        do c=1,cuantas         ! loop over cuantas
            lnpro=logweightchain(c) 
            do s=1,nseg        ! loop over segments 
                t=type_of_monomer(s)
                if(t/=ta) then 
                    k=indexchain(s,c)                
                    lnpro = lnpro +lnexppi(k,t)
                else 
                    ! phosphates 
                    k = indexchain(s,c)

                    do jj=1,nneigh(s,c)           ! loop neighbors 

                        m = indexconfpair(s,c)%elem(jj)

                        !call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)

                        !lnpro =lnpro + (lnexppi(k,ta) + lnexppi(m,ta)+ (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) &
                        !                  -log(fdisPP_loc(Phos,Phos))  )/(2.0_dp*nneigh(s,c))   

                        lnpro =lnpro +  (lnexppivw(k) + lnexppivw(m))*(vpol(tA)*vsol) /(2.0_dp*nneigh(s,c))                     
                    enddo
                endif           
            enddo 

            pro=exp(lnpro-lnproshift)   
            local_q = local_q+pro

            do s=1,nseg
                t = type_of_monomer(s)
                if(t/=ta) then 
                    k = indexchain(s,c) 
                    local_rhopol(k,t) = local_rhopol(k,t)+pro ! unnormed polymer density at k given that the 'beginning'of chain is at l
                else
                 
                    ! pair density of phosphates 
                    k = indexchain(s,c)

                    ! k_ind = inverse_index_phos(k) 

                    do j=1,nneigh(s,c)

                        m = indexconfpair(s,c)%elem(j)
                             
                        ! call  compute_fdisPP(fdisPP_loc, fdisP2Mg_loc, k , m)
                        ! call  compute_fdisPP(fdisPP_loc_swap, fdisP2Mg_loc_swap,  m , k) 
                        ! first part integral   

                        sum_rhoqphos = 0.0_dp
                        sum_xphos = 0.0_dp 
                        
                        do JJ=1,5
                            do KK=1,5
                                sum_rhoqphos = sum_rhoqphos+&
                                    (fdisPP_loc(JJ,KK)*qPP(JJ)+fdisPP_loc_swap(JJ,KK)*qPP(KK))/2.0_dp
                                sum_xphos = sum_xphos   +&
                                    (fdisPP_loc(JJ,KK)*vPP(JJ)+fdisPP_loc_swap(JJ,KK)*vPP(KK))/2.0_dp
                            enddo
                        enddo
        
                        sum_xphos = sum_xphos+(fdisP2Mg_loc+fdisP2Mg_loc_swap)*vPP(Phos2Mg)/4.0_dp 
                            ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                      
                        local_rhoqphos(k) = local_rhoqphos(k) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                        local_xpolphos(k) = local_xpolphos(k) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                        local_rhopol(k,ta)=local_rhopol(k,ta)+pro/(2.0_dp*nneigh(s,c))
                            
                        ! second integral contributes to location m of rhoqpos and xphol  xpol  
                         
                        sum_rhoqphos=0.0_dp
                        sum_xphos=0.0_dp 
                        
                        ! contributes to location m of rhoqpos and xol
                               
                        do JJ=1,5
                            do KK=1,5   
                                sum_rhoqphos = sum_rhoqphos+&
                                    (fdisPP_loc_swap(JJ,KK)*qPP(JJ)+fdisPP_loc(JJ,KK)*qPP(KK))/2.0_dp
                                sum_xphos = sum_xphos   +&
                                    (fdisPP_loc_swap(JJ,KK)*vPP(JJ)+fdisPP_loc(JJ,KK)*vPP(KK))/2.0_dp
                            enddo
                        enddo
        
                        sum_xphos=sum_xphos+(fdisP2Mg_loc_swap +fdisP2Mg_loc)*vPP(Phos2Mg)/4.0_dp

                        ! division 4.0_dp  because symmetry and  vPP(Phos2Mg)/2 is volume change per phosphate 
                      
                        local_rhoqphos(m) = local_rhoqphos(m) + pro * sum_rhoqphos /(2.0_dp*nneigh(s,c)) ! nneigh could be zero  hence with in loop 
                        local_xpolphos(m) = local_xpolphos(m) + pro * sum_xphos /(2.0_dp*nneigh(s,c))

                        local_rhopol(m,ta)=local_rhopol(m,ta)+pro/(2.0_dp*nneigh(s,c))

                    enddo 

                endif 
            
            enddo
        enddo

        !   .. import results 

       

        if (rank==0) then 

            q=0.0_dp
            g=1  
            q(g)=local_q
            
             do i=1, numproc-1
                source = i
                call MPI_RECV(local_q, 1, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat, ierr)             
                g =int(source/nset_per_graft)+1  ! nset_per_graft = int(size/ngr)
                q(g)=q(g)+local_q
            enddo

               

            ! first graft point 
            do t=1,nsegtypes
                do i=1,n
                    rhopol(i,t)=local_rhopol(i,t)/q(1) ! polymer density 
                enddo
            enddo

            do i=1,nsize
                xpol(i,tA) = local_xpolphos(i)/q(1)
                rhoqphos(i) = local_rhoqphos(i)/q(1) 
            enddo

            ! other graft points 
            do i=1, numproc-1
                source = i
                g =int(source/nset_per_graft)+1 
               
                do t=1,nsegtypes
                    call MPI_RECV(local_rhopol(:,t), nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                    do k=1,nsize
                        rhopol(k,t)=rhopol(k,t)+local_rhopol(k,t)/q(g) ! polymer density 
                    enddo    
                enddo

                call MPI_RECV(local_xpolphos, nsize,   MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)     
                call MPI_RECV(local_rhoqphos, nsize, MPI_DOUBLE_PRECISION,source,tag,MPI_COMM_WORLD,stat,ierr)
                
                do k=1,nsize
                    xpol(k,tA) = xpol(k,tA) + local_xpolphos(k)/q(g)      ! phosphate volume fraction
                    rhoqphos(i) = rhoqphos(i) + local_rhoqphos(k)/q(g)
                enddo 
            enddo  
        
            !     .. construction of fcn and volume fraction polymer             
            rhopol0=(1.0_dp/volcell)! volume polymer segment per volume cell

            do t=1, nsegtypes

                if(ismonomer_chargeable(t)) then 

                    if(t/=ta) then
                                                              
                        do i=1,n
                            rhopol(i,t) = rhopol0 * rhopol(i,t)     ! density polymer of type t  
                            xpol(i,t)   = rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer

                            rhoqpol(i)  =  rhoqpol(i) + &
                                    (zpol(t,2)*fdis(i,t)+zpol(t,1)*(1.0_dp-fdis(i,t)))*rhopol(i,t)*vsol ! total  charge density in units of vsol 
                            
                           ! f(i+(t+1)*n)    = rhopol(i,t) - rhopolin(i,t)         ! scf eq for density
                        enddo  

                    else
                        ! phophate t=tA
                         do i=1,n
                            rhopol(i,ta) = rhopol0 * rhopol(i,ta) 
                            rhoqphos(i)  = rhopol0 * rhoqphos(i) 
                            xpol(i,ta)   = rhopol0 * xpol(i,tA)

                            rhoqpol(i)   = rhoqpol(i) + rhoqphos(i) * vsol ! total  charge density in units of vsol 

                        enddo           
                    endif    
                else  
                     ! neutral monomeer
                    do i=1,n
                        rhopol(i,t)  = rhopol0 * rhopol(i,t)               ! density polymer of type t  
                        xpol(i,t)    = xpol(i,t) + rhopol(i,t)*vpol(t)*vsol  ! volume fraction polymer
                       ! f(i+(t+1)*n) = rhopol(i,t) - rhopolin(i,t)         ! scf eq for density
                    enddo
                endif   

                do i=1,n
                    xpol_tot(i) = xpol_tot(i)+xpol(i,t)  
                enddo

            enddo    

            do i=1,n
                f(i) = xpol_tot(i)+xsol(i)+xNa(i)+xCl(i)+xHplus(i)+xOHmin(i)+xRb(i)+xCa(i)+xMg(i)+xK(i)+xpro(i) -1.0_dp
                rhoq(i) = rhoqpol(i)+zNa*xNa(i)/vNa +zCl*xCl(i)/vCl +xHplus(i)-xOHmin(i)+ &
                    zCa*xCa(i)/vCa +zMg*xMg(i)/vMg+zRb*xRb(i)/vRb +zK*xK(i)/vK ! total charge density in units of vsol  
            enddo
          
            ! .. end computation polymer density and charge density  

            ! .. electrostatics 

            sigmaqSurfR=surface_charge(bcflag(RIGHT),psiSurfR,RIGHT)
            sigmaqSurfL=surface_charge(bcflag(LEFT),psiSurfL,LEFT)
            
            ! .. Poisson Eq 
            call Poisson_Equation(f,psi,rhoq,sigmaqSurfR,sigmaqSurfL)
    
            ! .. boundary conditions
            call Poisson_Equation_Surface(f,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)    

            norm=l2norm_f90(f)
            iter=iter+1
                        
            normvol = L2norm_f90(f(1:nsize))
            normPE  = L2norm_f90(f(nsize+1:2*nsize))
                        
            print*,'iter=', iter ,'norm=',norm, "normvol=",normvol,"normPE=",normPE
           
            
        else                      ! Export results 
            
            dest = 0 
           
            call MPI_SEND(local_q, 1 , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)

            do t=1,nsegtypes
                call MPI_SEND(local_rhopol(:,t),nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
            enddo
           
            call MPI_SEND(local_xpolphos, nsize, MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)      
            call MPI_SEND(local_rhoqphos, nsize , MPI_DOUBLE_PRECISION, dest,tag, MPI_COMM_WORLD, ierr)
    
        endif

    end subroutine fcn_neutral_expl



end module modfcnMgexpl

   

