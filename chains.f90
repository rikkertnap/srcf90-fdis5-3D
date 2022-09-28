!     .. module file of chains variables
module chains
  
    use globals
    use volume, only: ngr

    implicit none
  
    integer, dimension(:,:), allocatable :: indexchain               ! index(alpha,s)= layer number of conf alpha and segment number s
    integer, dimension(:,:), allocatable :: indexchain_init 
    logical, dimension(:), allocatable      :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable      :: type_of_monomer          ! type of monomer represented as a number
    character(len=2), dimension(:), allocatable :: type_of_monomer_char ! type of monomer represented as two letters
    logical, dimension(:,:), allocatable    :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable      :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    real(dp), dimension(:), allocatable     :: energychain              ! energy chain   
    real(dp), dimension(:), allocatable     :: energychain_init         ! energy chain   
    real(dp) :: energychain_min              ! mimimum energy chain
    real(dp), dimension(:),   allocatable    :: logweightchain
    logical :: isHomopolymer
    double precision, dimension(:),allocatable :: lsegseq  ! only needed for copolymer

    ! chain structural quantities
    real(dp), dimension(:), allocatable     :: Rgsqr                    ! radius of gyration (for all conformations) 
    real(dp), dimension(:), allocatable     :: Rendsqr                  ! end-to-end distance (for all conformations)
    real(dp), dimension(:,:,:), allocatable   :: As_mtrx                   ! "asphericity matrix"/gyration tensor (for all conformations)

    real(dp), dimension(:), allocatable     :: avRgsqr                    ! average radius of gyration (for each graft point)
    real(dp), dimension(:), allocatable     :: avRendsqr                  ! average end-to-end distance (for each graft point)
    real(dp), dimension(:,:,:), allocatable   :: avAs_mtrx                  ! average "asphericity matrix"/gyration tensor (for each graft point)

contains


    subroutine allocate_chains(cuantas,nseg,nsegtypes,maxnchains,maxnchainsxy)

        integer, intent(in) :: cuantas,nseg,nsegtypes
        integer, intent(in) :: maxnchains,maxnchainsxy

        integer :: maxcuantas
    
        maxcuantas=cuantas+maxnchains*maxnchainsxy     ! .. extra  because of  nchain rotations
        allocate(indexchain(nseg,maxcuantas))
        allocate(indexchain_init(nseg,maxcuantas))
        allocate(energychain(maxcuantas))
        allocate(energychain_init(maxcuantas))
        allocate(logweightchain(maxcuantas))
        allocate(isAmonomer(nseg)) 
        allocate(type_of_monomer(nseg)) 
        allocate(type_of_monomer_char(nseg))
        allocate(ismonomer_of_type(nseg,nsegtypes)) 
        allocate(ismonomer_chargeable(nsegtypes))

        ! chain structural quantities

        allocate(Rgsqr(maxcuantas))
        allocate(Rendsqr(maxcuantas))
        allocate(As_mtrx(maxcuantas,3,3))
        allocate(avRgsqr(ngr))
        allocate(avRendsqr(ngr))
        allocate(avAs_mtrx(ngr,3,3))

    end subroutine allocate_chains
  
end module chains
