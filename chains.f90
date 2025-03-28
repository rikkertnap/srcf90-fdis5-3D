!     .. module file of chains variables
module chains
  
    use globals
    use volume, only: ngr

    implicit none
  
    type var_iarray
        integer, allocatable :: elem(:)
    end type var_iarray

    type(var_iarray), allocatable               :: indexconfpair(:,:)       ! indexconfpair(s,alpha)%elem(j) = layer number of conf alpha and 
                                                                            ! segment number s and neighbor j used for distributed volume 
    integer, dimension(:,:), allocatable        :: nneigh                   ! number of neigbors or pairs of segment s in conf alpha used only phosphates  
    integer, dimension(:), allocatable          :: max_nneigh_phos          ! maximum number of neigbors or pairs of phosphate segments in conf alpha and seggment s

    integer, dimension(:,:), allocatable    :: indexchain               ! index(alpha,s)= layer number of conf alpha and segment number s
    integer, dimension(:,:), allocatable    :: indexchain_init 
    logical, dimension(:), allocatable      :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable      :: type_of_monomer          ! type of monomer represented as a number
    character(len=2), dimension(:), allocatable :: type_of_monomer_char ! type of monomer represented as two letters
    logical, dimension(:,:), allocatable    :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable      :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    real(dp), dimension(:), allocatable     :: energychain              ! energy chain   
    real(dp), dimension(:), allocatable     :: energychain_init         ! energy chain   
    real(dp) :: energychain_min                                         ! mimimum energy chain
    real(dp), dimension(:),   allocatable    :: logweightchain
    logical :: isHomopolymer
    double precision, dimension(:),allocatable :: lsegseq               ! only needed for copolymer

    ! chain structural quantities

    real(dp), dimension(:), allocatable       :: Rgsqr                  ! radius of gyration (for all conformations) 
    real(dp), dimension(:), allocatable       :: Rendsqr                ! end-to-end distance (for all conformations)
    real(dp), dimension(:), allocatable       :: avRgsqr                ! average radius of gyration (for each graft point)
    real(dp), dimension(:), allocatable       :: avRendsqr              ! average end-to-end distance (for each graft point)
    real(dp), dimension(:), allocatable       :: Asphparam              ! Asphericity parameter invariant of gyration tensor 
    real(dp), dimension(:), allocatable       :: avAsphparam            ! Average asphericity parameter for each graft point 
   
    ! .. pairing parameters 

    real(dp) :: distphoscutoff ! distance allow between two phosphate to be a pair
    integer  :: maxneigh       ! maximum of neigbors 
    integer  :: len_index_phos ! length of array index_phos

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
        allocate(Asphparam(maxcuantas))
        allocate(avRgsqr(ngr))
        allocate(avRendsqr(ngr))
        allocate(avAsphparam(ngr))
    
    end subroutine allocate_chains

    ! Allocates indexconfpair(s,alpha)%elem(j) = layer number of conf alpha and segment number s and neigbor pair j
    ! used for distributed volume
    ! When used indexchain is not needed and can be deallocated
    ! inputs: dimension of indexconfpair: cuantas, nseg and  nelem(:) 

    subroutine allocate_indexconfpair(cuantas,nseg)

        integer, intent(in) :: cuantas,nseg

        allocate(indexconfpair(nseg,cuantas))  
        
    end subroutine allocate_indexconfpair

   
    ! Allocates neigh : neigbors that segment number s in conf alpha has 
    ! used only for segment s that is a phosphate

    subroutine allocate_nneighbor(cuantas,nseg)

        integer, intent(in) :: cuantas,nseg

        allocate(nneigh(nseg,cuantas))
                  
    end subroutine allocate_nneighbor

     ! Allocates addition maxneigh_phos

    subroutine allocate_max_nneighbor_phos(cuantas)

        integer, intent(in) :: cuantas

        allocate(max_nneigh_phos(cuantas))

    end subroutine allocate_max_nneighbor_phos


end module chains
