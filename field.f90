module field
  
  !     .. variables
    use precision_definition

    implicit none
    
    real(dp), dimension(:), allocatable :: xpol     ! volume fraction of polymer   
    real(dp), dimension(:), allocatable :: xpolz ! total volume fraction of polymer in z-direction  
    real(dp), dimension(:,:), allocatable :: rhopol ! density  monomer of polymer in layer i of type t
    real(dp), dimension(:,:), allocatable :: rhopolin 
    real(dp), dimension(:), allocatable :: rhoqpol  ! charge density  monomer of polymer in layer i 

    real(dp), dimension(:), allocatable :: xpolAB  ! total volume fraction of polymer 
    !real(dp), dimension(:), allocatable :: xpolABz ! total volume fraction of polymer in z-direction
    real(dp), dimension(:), allocatable :: rhopolA ! density A monomer of polymer 
    real(dp), dimension(:), allocatable :: rhopolB ! density B monomer of polymer 

    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: psi     ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xK      ! volume fraction of positive K+ ion
    real(dp), dimension(:), allocatable :: xRb      ! volume fraction of positive Rb+ ion
    real(dp), dimension(:), allocatable :: xCa     ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xMg     ! volume fraction of positive Mg2+ ion    
    real(dp), dimension(:), allocatable :: xNaCl   ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl    ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl     ! volume fraction of negative ion
    real(dp), dimension(:), allocatable :: xHplus  ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin  ! volume fraction of OHmin 

    real(dp), dimension(:), allocatable :: rhoq    ! total free charge density in units of vsol  
    !real(dp), dimension(:), allocatable :: qpol    ! charge density of polymer
    real(dp),dimension(:), allocatable :: epsfcn    ! dielectric constant 
    real(dp),dimension(:), allocatable :: Depsfcn   ! derivative dielectric constant

    real(dp), dimension(:,:), allocatable :: fdis   ! degree of dissociation of acid monomer
    real(dp), dimension(:,:), allocatable :: fdisA   ! degree of dissociation 
    real(dp), dimension(:,:), allocatable :: fdisB   ! degree of dissociation
     
    real(dp), dimension(:), allocatable :: q      ! normalization partion fnc polymer 
    
    real(dp), dimension(:), allocatable :: qAB      ! normalization partion fnc polymer 
    real(dp), dimension(:), allocatable :: qABL,qABR

    real(dp), dimension(:), allocatable :: rhopolAL ! density A monomer of polymer z=0 surface
    real(dp), dimension(:), allocatable :: rhopolBL ! density B monomer of polymer z=0 surface
    real(dp), dimension(:), allocatable :: rhopolAR ! density A monomer of polymer z=nz delta surface 
    real(dp), dimension(:), allocatable :: rhopolBR ! density B monomer of polymer z=nz delta surface

  
contains

    subroutine allocate_field(Nx,Ny,Nz,nsegtypes)
 
        integer, intent(in) :: Nx,Ny,Nz,nsegtypes
        
        integer :: N
        integer :: ier

        N=Nx*Ny*Nz

        allocate(xpol(N),stat=ier)
        allocate(xpolz(Nz),stat=ier)
        allocate(rhopol(N,nsegtypes),stat=ier) 
        allocate(rhopolin(N,nsegtypes),stat=ier) 
        allocate(rhoqpol(N),stat=ier) 
        allocate(xpolAB(N),stat=ier)
        allocate(rhopolA(N),stat=ier)
        allocate(rhopolB(N),stat=ier)
        allocate(xsol(N),stat=ier)
        allocate(psi(N+2*Nx*Ny))
        allocate(xNa(N),stat=ier)
        allocate(xK(N),stat=ier)
        allocate(xRb(N),stat=ier)
        allocate(xCa(N),stat=ier)
        allocate(xMg(N),stat=ier)
        allocate(xNaCl(N),stat=ier) 
        allocate(xKCl(N),stat=ier) 
        allocate(xCl(N),stat=ier) 
        allocate(xHplus(N),stat=ier)
        allocate(xOHmin(N),stat=ier)
        allocate(rhoq(N),stat=ier)
        allocate(fdis(N,nsegtypes),stat=ier)
        allocate(fdisA(N,7),stat=ier)
        allocate(fdisB(N,5),stat=ier)
        allocate(rhopolAL(N),stat=ier)
        allocate(rhopolAR(N),stat=ier)
        allocate(rhopolBL(N),stat=ier)
        allocate(rhopolBR(N),stat=ier)
        allocate(epsfcn(N),stat=ier)    ! relative dielectric constant
        allocate(Depsfcn(N),stat=ier)   ! derivate relative dielectric constant


        if( ier/=0 ) then
            print*, 'Allocation error : stat =', ier
            stop
        endif
        
    end subroutine allocate_field


    subroutine deallocate_field()
        
        deallocate(xpol)
        deallocate(xpolz)
        deallocate(rhopol)
        deallocate(rhoqpol)
        deallocate(xpolAB)
        deallocate(rhopolA)
        deallocate(rhopolB)
        deallocate(xsol)
        deallocate(psi)
        deallocate(xNa)
        deallocate(xK)
        deallocate(xCa)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl) 
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(rhoq)
        deallocate(fdisA)
        deallocate(fdisB)
        deallocate(rhopolAL)
        deallocate(rhopolAR)
        deallocate(rhopolBL)
        deallocate(rhopolBR)
        deallocate(epsfcn)
        deallocate(Depsfcn)
        
    end subroutine deallocate_field


    subroutine allocate_part_fnc(N)
       
        integer, intent(in) :: N

        allocate(qAB(N))
        allocate(qABL(N))
        allocate(qABR(N))
        allocate(q(N))

    end subroutine allocate_part_fnc

    ! set all densities to zero
    
    subroutine init_field()

        xpol=0.0_dp
        xpolAB=0.0_dp
        rhopol=0.0_dp
        rhoqpol=0.0_dp
        rhopolA=0.0_dp
        rhopolB=0.0_dp
        xsol=0.0_dp
        xNa=0.0_dp
        xK=0.0_dp
        xCa=0.0_dp
        xNaCl=0.0_dp 
        xKCl =0.0_dp
        xCl=0.0_dp
        xHplus=0.0_dp
        xOHmin=0.0_dp
        rhoq=0.0_dp
        rhopolAL=0.0_dp
        rhopolAR=0.0_dp
        rhopolBL=0.0_dp
        rhopolBR=0.0_dp
        fdisA=0.0_dp
        fdisB=0.0_dp
        xpolz=0.0_dp
        psi=0.0_dp
           
    end subroutine init_field


    !  debug routine

    subroutine check_integral_rholpol_multi(sumrhopol, checkintegral)

        use volume, only : volcell, ngr
        use globals, only : nsize, systype, nseg, nsegtypes

        real(dp), intent(inout) :: sumrhopol,checkintegral 
        integer :: t,i
        real(dp) :: intrhopol

        sumrhopol=0.0_dp
        do t=1,nsegtypes
            do i=1,nsize
                sumrhopol=sumrhopol+rhopol(i,t)
            enddo  
        enddo      
        sumrhopol=sumrhopol*volcell

        intrhopol=nseg*ngr
        if(systype=="electdouble") intrhopol=intrhopol*2.0_dp  

        checkintegral=sumrhopol-intrhopol

    end subroutine

    subroutine check_integral_rholpolAB(sumrhopol, checkintegral)

        use volume, only : volcell, ngr
        use globals, only : nsize, systype, nseg

        real(dp), intent(inout) :: sumrhopol,checkintegral 
        integer :: i
        real(dp) :: intrhopol

        sumrhopol=0.0_dp
        do i=1,nsize
            sumrhopol=sumrhopol+(rhopolA(i)+rhopolB(i))
        enddo    
        sumrhopol=sumrhopol*volcell

        intrhopol=nseg*ngr
        if(systype=="electdouble") intrhopol=intrhopol*2.0_dp  

        checkintegral=sumrhopol-intrhopol

    end subroutine
       

    subroutine charge_polymer()

        use globals, only : systype
        
        select case (systype) 
        case ("brush_mul")
            call charge_polymer_multi()
        case ("brushssdna","brushborn")
            call charge_polymer_dna()
        case ("elect","electA","electVdWAB","electdouble")   
            call charge_polymer_binary()
        case default
            print*,"Error in average_charge_polymer subroutine"    
            print*,"Wrong value systype : ", systype
            stop
        end select  
       

    end subroutine charge_polymer


    subroutine charge_polymer_dna()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, qpol, qpol_tot, tA

        integer :: i, t


        qpol_tot=0.0_dp
        do t=1,nsegtypes
            qpol(t)=0.0_dp
            if(t/=tA) then    
                do i=1,nsize
                    qpol(t)=qpol(t)+(fdis(i,t)*zpol(t,1)+(1.0_dp-fdis(i,t))*zpol(t,2))*rhopol(i,t)
                enddo
            else
                do i=1,nsize
                    qpol(t)=qpol(t)+ (-fdisA(i,1)*rhopol(i,1)+fdisA(i,4)*rhopol(i,t))
                enddo
            endif    

            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_polymer_dna





    subroutine charge_polymer_multi()

        use globals, only : nsize, nsegtypes
        use volume, only : volcell
        use parameters, only : zpol, qpol, qpol_tot

        integer :: i, t

        qpol_tot=0.0_dp
        do t=1,nsegtypes
            qpol(t)=0.0_dp
            do i=1,nsize
                qpol(t)=qpol(t)+(fdis(i,t)*zpol(t,1)+(1.0_dp-fdis(i,t))*zpol(t,2))*rhopol(i,t)
            enddo
            qpol(t)=qpol(t)*volcell
            qpol_tot=qpol_tot+qpol(t)
        enddo

    end subroutine charge_polymer_multi

    subroutine charge_polymer_binary()

        use globals, only : nsize
        use volume, only : volcell
        use parameters, only : zpolA, zpolB, qpolA,qpolB, qpol_tot

        integer :: i

        qpolA=0.0_dp
        qpolB=0.0_dp

        do i=1,nsize
            qpolA=qpolA+(zpolA(1)*fdisA(i,1)*rhopol(i,1)+zpolA(4)*fdisA(i,4)*rhopol(i,1))
            qpolB=qpolB+(zpolB(1)*fdisB(i,1)*rhopol(i,2)+zpolB(4)*fdisB(i,4)*rhopol(i,2))
        enddo

        qpolA=qpolA*volcell
        qpolB=qpolB*volcell
        qpol_tot=qpolA+qpolB

    end subroutine charge_polymer_binary

    ! .. post : return average charge of state of polymer

    subroutine average_charge_polymer()

        use globals, only : systype
        
        select case (systype) 
        case ("brush_mul")
            call average_charge_polymer_multi()
        case ("brushssdna","brushborn")
            call average_charge_polymer_dna()
        case ("elect","electA","electVdWAB","electdouble")   
            call average_charge_polymer_binary()
        case default
            print*,"Error in average_charge_polymer subroutine"    
            print*,"Wrong value systype : ", systype
            stop
        end select  

    end subroutine average_charge_polymer
        

     subroutine average_charge_polymer_dna()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell,ngr
        use parameters, only : zpol, avfdis, avfdisA, tA
        use chains, only: type_of_monomer

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t,k

        allocate(npol(nsegtypes))
        npol=0.0_dp

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   

        do t=1,nsegtypes
            avfdis(t)=0.0_dp
            if(npol(t)/=0) then
                if(t/=tA) then    
                    do i=1,nsize
                        avfdis(t)=avfdis(t)+(fdis(i,t)*zpol(t,1)+(1.0_dp-fdis(i,t))*zpol(t,2))*rhopol(i,t)
                    enddo
                    avfdis(t)=avfdis(t)*volcell/(npol(t)*ngr)        
                else
                    do k=1,7
                        avfdisA(k)=0.0_dp
                        do i=1,nsize
                            avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol(i,t)
                        enddo
                        avfdisA(k)=avfdisA(k)*volcell/(npol(t)*ngr)
                    enddo
                    avfdis(t)=avfdisA(1)
                endif       
            endif
        enddo         

        deallocate(npol)    

    end subroutine average_charge_polymer_dna

    subroutine average_charge_polymer_multi()

        use globals, only : nseg,nsize,nsegtypes
        use volume, only : volcell,ngr
        use parameters, only : zpol, avfdis
        use chains, only: type_of_monomer

        integer, dimension(:), allocatable   :: npol
        integer :: i,s,t

        allocate(npol(nsegtypes))
        npol=0.0_dp

        do s=1,nseg
            t=type_of_monomer(s)
            npol(t)=npol(t)+1
        enddo   

        do t=1,nsegtypes
            if(npol(t)/=0) then
                avfdis(t)=0.0_dp
                do i=1,nsize
                    avfdis(t)=avfdis(t)+(fdis(i,t)*zpol(t,1)+(1.0_dp-fdis(i,t))*zpol(t,2))*rhopol(i,t)
                enddo
                avfdis(t)=avfdis(t)*volcell/(npol(t)*ngr)        
            else
                avfdis(t)=0.0_dp
            endif
        enddo         

        deallocate(npol)    

    end subroutine average_charge_polymer_multi
        

    subroutine average_charge_polymer_binary()
        
        use globals, only : nseg,nsize
        use volume, only : volcell
        use parameters
        use chains, only : isAmonomer

        integer :: i,s,k
        integer   :: npolA,npolB
        integer, parameter :: A=1, B=2

        ! .. number of A and B monomors 
        npolA=0
        do s=1,nseg
           if(isAmonomer(s).eqv..true.) then
              npolA=npolA+1
           endif
        enddo
        npolB=nseg-npolA
          

        if(npolA/=0) then
           do k=1,5
              avfdisA(k)=0.0_dp
              do i=1,nsize
                 avfdisA(k)=avfdisA(k)+fdisA(i,k)*rhopol(i,A)
              enddo
              avfdisA(k)=avfdisA(k)*volcell/(npolA*ngr)
           enddo
        else
           avfdisA=0.0_dp
        endif

        if(npolB/=0) then
           do k=1,5
              avfdisB(k)=0.0_dp
              do i=1,nsize
                 avfdisB(k)=avfdisB(k)+fdisB(i,k)*rhopol(i,B)
              enddo
              avfdisB(k)=avfdisB(k)*volcell/(npolB*ngr)
           enddo
        else
           do k=1,5
              avfdisB(k)=0.0_dp
           enddo
        endif

    end subroutine average_charge_polymer_binary

    !     .. compute average height of denisty provile
    !     .. first moment of density profile 
  
    function average_height_z(rho) result(meanz)

        use volume, only : nz,nx,ny,delta,linearIndexFromCoordinate

        real(dp), intent(in) :: rho(:)
        real(dp) :: meanz    

        integer :: ix, iy, iz, id
        real(dp) :: sumrhoz, sumrho

        sumrhoz = 0.0_dp
        meanz= 0.0_dp

        do iz = 1, nz

            sumrho = 0.0_dp       
            do ix=1, nx
                do iy=1, ny
                    call linearIndexFromCoordinate(ix,iy,iz ,id)
                    sumrho=sumrho+ rho(id)
                enddo
            enddo

            meanz=meanz+sumrho*(iz-0.5_dp)*delta
            sumrhoz=sumrhoz+sumrho
        enddo

        if(sumrhoz>0.0_dp) then 
            meanz=meanz/sumrhoz
        else
            meanz=0.0_dp
        endif

    end function average_height_z

    !     .. compute average of density or volume fraction profile in z-direction
  
    subroutine average_density_z(xvol,xvolz,meanz)

        use volume, only : nz,nx,ny,delta,linearIndexFromCoordinate

        real(dp), intent(in) :: xvol(:)
        real(dp), intent(out) :: xvolz(:)
        real(dp), intent(out), optional :: meanz

        integer :: ix, iy, iz, id
        real(dp) :: sumrhoz, sumxvol

        if(present(meanz)) then 

            sumrhoz = 0.0_dp
            do iz = 1, nz
                sumxvol = 0.0_dp       
                do ix=1, nx
                    do iy=1, ny
                        call linearIndexFromCoordinate(ix,iy,iz ,id)
                        sumxvol=sumxvol+ xvol(id)
                    enddo
                enddo
                xvolz(iz)=sumxvol/(1.0_dp*nx*ny)

                meanz=meanz+sumxvol*(iz-0.5_dp)*delta
                sumrhoz=sumrhoz+sumxvol
            enddo

            if(sumrhoz>0.0_dp) then 
                meanz=meanz/sumrhoz
            else
                meanz=0.0_dp
            endif
        else

            sumrhoz = 0.0_dp
            do iz = 1, nz
                sumxvol = 0.0_dp       
                do ix=1, nx
                    do iy=1, ny
                        call linearIndexFromCoordinate(ix,iy,iz ,id)
                        sumxvol=sumxvol+ xvol(id)
                    enddo
                enddo
                xvolz(iz)=sumxvol/(1.0_dp*nx*ny)
            enddo

        endif
            
    end subroutine average_density_z
  
end module field

