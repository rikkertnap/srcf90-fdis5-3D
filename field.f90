module field
  
  !     .. variables
    use precision_definition

    implicit none
    
    real(dp), dimension(:), allocatable :: xpolAB  ! total volume fraction of polymer 
    real(dp), dimension(:), allocatable :: xpolABz ! total volume fraction of polymer in z-direction
    real(dp), dimension(:), allocatable :: xpolC   ! total volume fraction of polymer 
    real(dp), dimension(:), allocatable :: rhopolA ! density A monomer of polymer 
    real(dp), dimension(:), allocatable :: rhopolB ! density B monomer of polymer 
    real(dp), dimension(:), allocatable :: rhopolC ! density C monomer of polymer
    real(dp), dimension(:), allocatable :: xsol    ! volume fraction solvent
    real(dp), dimension(:), allocatable :: psi     ! electrostatic potential 
    real(dp), dimension(:), allocatable :: xNa     ! volume fraction of positive Na+ ion
    real(dp), dimension(:), allocatable :: xK      ! volume fraction of positive K+ ion
    real(dp), dimension(:), allocatable :: xCa     ! volume fraction of positive Ca2+ ion
    real(dp), dimension(:), allocatable :: xNaCl   ! volume fraction of NaCl ion pair
    real(dp), dimension(:), allocatable :: xKCl    ! volume fraction of KCl  ion pair
    real(dp), dimension(:), allocatable :: xCl     ! volume fraction of negative ion
    real(dp), dimension(:), allocatable :: xHplus  ! volume fraction of Hplus
    real(dp), dimension(:), allocatable :: xOHmin  ! volume fraction of OHmin 
    real(dp), dimension(:), allocatable :: rhoq    ! total free charge density in units of vsol  
    real(dp), dimension(:), allocatable :: rhob    ! total bond charge density in units of vsol
    real(dp), dimension(:,:), allocatable :: electPol ! elect polarization in 3d-direction 
    real(dp), dimension(:), allocatable :: qpol    ! charge density of polymer
    real(dp), dimension(:,:), allocatable :: fdisA   ! degree of dissociation 
    real(dp), dimension(:,:), allocatable :: fdisB   ! degree of dissociation
     
    real(dp), dimension(:), allocatable :: qAB      ! normalization partion fnc polymer 
    real(dp), dimension(:), allocatable :: qC      ! normalization partion fnc polymer 

    real(dp), dimension(:), allocatable :: rhopolAL ! density A monomer of polymer z=0 surface
    real(dp), dimension(:), allocatable :: rhopolBL ! density B monomer of polymer z=0 surface
    real(dp), dimension(:), allocatable :: rhopolAR ! density A monomer of polymer z=nz delta surface 
    real(dp), dimension(:), allocatable :: rhopolBR ! density B monomer of polymer z=nz delta surface

    real(dp), dimension(:), allocatable :: qABL,qABR

  
contains

    subroutine allocate_field(Nx,Ny,Nz)
 
        integer, intent(in) :: Nx,Ny,Nz
        integer :: N

        N=Nx*Ny*Nz

        allocate(xpolAB(N))
        allocate(xpolABz(Nz))
        allocate(xpolC(N))
        allocate(rhopolA(N))
        allocate(rhopolB(N))
        allocate(rhopolC(N))
        allocate(xsol(N))
        allocate(psi(N+2*Nx*Ny))
        allocate(electPol(3,N+2*Nx*Ny))
        allocate(xNa(N))
        allocate(xK(N))
        allocate(xCa(N))
        allocate(xNaCl(N)) 
        allocate(xKCl(N)) 
        allocate(xCl(N)) 
        allocate(xHplus(N))
        allocate(xOHmin(N))
        allocate(rhoq(N))
        allocate(rhob(N))
        allocate(qpol(N))
        allocate(fdisA(5,N))
        allocate(fdisB(5,N))
        allocate(rhopolAL(N))
        allocate(rhopolAR(N))
        allocate(rhopolBL(N))
        allocate(rhopolBR(N))
        
    end subroutine allocate_field


    subroutine deallocate_field()
        
        deallocate(xpolAB)
        deallocate(xpolABz)
        deallocate(xpolC)
        deallocate(rhopolA)
        deallocate(rhopolB)
        deallocate(rhopolC)
        deallocate(xsol)
        deallocate(psi)
        deallocate(electPol)
        deallocate(xNa)
        deallocate(xK)
        deallocate(xCa)
        deallocate(xNaCl) 
        deallocate(xKCl) 
        deallocate(xCl) 
        deallocate(xHplus)
        deallocate(xOHmin)
        deallocate(rhoq)
        deallocate(rhob)
        deallocate(qpol)
        deallocate(fdisA)
        deallocate(fdisB)
        deallocate(rhopolAL)
        deallocate(rhopolAR)
        deallocate(rhopolBL)
        deallocate(rhopolBR)
        
    end subroutine deallocate_field


    subroutine allocate_part_fnc(N)
       
        integer, intent(in) :: N

        allocate(qAB(N))
        allocate(qABL(N))
        allocate(qABR(N))
        allocate(qC(N))

    end subroutine allocate_part_fnc


    !  debug routine

    subroutine check_integral_rholpolAB(sumrhopolAB, checkintegralAB)

        use volume, only : volcell, ngr
        use globals, only : nsize, nsegAB, sysflag

        real(dp), intent(inout) :: sumrhopolAB,checkintegralAB 
        integer :: i
        real(dp) :: intrhopolAB

        sumrhopolAB=0.0_dp
        do i=1,nsize
            sumrhopolAB=sumrhopolAB+(rhopolA(i)+rhopolB(i))
        enddo    
        sumrhopolAB=sumrhopolAB*volcell

        intrhopolAB=nsegAB*ngr
        if(sysflag=="electdouble") intrhopolAB=intrhopolAB*2.0_dp  

        checkintegralAB=sumrhopolAB-intrhopolAB

    end subroutine
        

    subroutine charge_polymer()

        use globals
        use volume
        use parameters

        implicit none

        integer :: i

        qpolA=0.0_dp
        qpolB=0.0_dp

        do i=1,nz
            qpolA=qpolA+(zpolA(1)*fdisA(1,i)*rhopolA(i)+&
                zpolA(4)*fdisA(4,i)*rhopolA(i))*deltaG(i) 
            qpolB=qpolB+(zpolB(1)*fdisB(1,i)*rhopolB(i)+&
                zpolB(4)*fdisB(4,i)*rhopolB(i))*deltaG(i)
        enddo

        qpolA=qpolA*delta
        qpolB=qpolB*delta
        qpol_tot=qpolA+qpolB

    end subroutine charge_polymer

    ! .. post : return average charge of state of 
    !   of polymers

    subroutine average_charge_polymer()
        
        use globals
        use volume
        use parameters
        use chains

        implicit none 

        integer :: i,s,k
        integer   :: npolA,npolB
        real(dp)  :: sigmaLR  

        ! .. number of A and B monomors 
        npolA=0
        do s=1,nsegAB
           if(isAmonomer(s).eqv..true.) then
              npolA=npolA+1
           endif
        enddo
        npolB=nsegAB-npolA
          
        sigmaLR=sigmaABL+sigmaABR

        if(npolA/=0 .and. sigmaLR/=0 ) then
           do k=1,5
              avfdisA(k)=0.0_dp
              do i=1,nz
                 avfdisA(k)=avfdisA(k)+fdisA(k,i)*rhopolA(i)*deltaG(i) 
              enddo
              avfdisA(k)=avfdisA(k)*delta/(sigmaLR*delta*npolA)
           enddo
        else
           do k=1,5
              avfdisA(k)=0.0_dp
           enddo
        endif

        if(npolB/=0 .and. sigmaLR/=0) then
           do k=1,5
              avfdisB(k)=0.0_dp
              do i=1,nz
                 avfdisB(k)=avfdisB(k)+fdisB(k,i)*rhopolB(i)*deltaG(i) 
              enddo
              avfdisB(k)=avfdisB(k)*delta/(sigmaLR*delta*npolB)
           enddo
        else
           do k=1,5
              avfdisB(k)=0.0_dp
           enddo
        endif

    end subroutine average_charge_polymer
  

  !     .. compute average height of tethered layer 
  !     .. first moment of density profile 
  
    subroutine average_polymer(xpol,xpolz,meanz)

        use volume, only : nz,nx,ny,delta,linearIndexFromCoordinate

        real(dp), intent(in) :: xpol(:)
        real(dp), intent(out) :: xpolz(:)
        real(dp), intent(out) :: meanz    

        integer :: ix, iy, iz, id
        real(dp) :: sumrhoz, sumxpol

        sumrhoz = 0.0_dp
        meanz= 0.0_dp

        do iz = 1, nz

            sumxpol = 0.0_dp       
            do ix=1, nx
                do iy=1, ny
                    call linearIndexFromCoordinate(ix,iy,iz ,id)
                    sumxpol=sumxpol+ xpol(id)
                enddo
            enddo

            xpolz(iz)=sumxpol/(1.0_dp*nx*ny)

            meanz=meanz+sumxpol*(iz-0.5_dp)*delta
            sumrhoz=sumrhoz+sumxpol
        enddo

        if(sumrhoz>0.0_dp) then 
            meanz=meanz/sumrhoz
        else
            meanz=0.0_dp
        endif

    end subroutine average_polymer

   
  
end module field

