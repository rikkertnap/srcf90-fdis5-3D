!     .. module file of chains variables
module chains
  
    use globals
    implicit none
  
    integer(2), dimension(:,:), allocatable :: indexchain               ! index(alpha,s)= layer number of conf alpha and segment number s
    integer(2), dimension(:,:), allocatable :: indexchain_init 
    logical, dimension(:), allocatable      :: isAmonomer               ! isAmonomer(s) =.true. if s is a "A" monomoer  
    integer, dimension(:), allocatable      :: type_of_monomer          ! type of monomer represented as a number
    character(len=2), dimension(:), allocatable :: type_of_monomer_char ! type of monomer represented as two letters
    logical, dimension(:,:), allocatable    :: ismonomer_of_type        ! ismomomer_of_type(s,t)= true if segment number "s" is of type "t" otherwise false 
    logical, dimension(:), allocatable      :: ismonomer_chargeable     ! ismonomer_chargeabl(s)=true if segment number type "t" is acid or base  
    real(dp), dimension(:), allocatable     :: exp_energychain          ! energy chain   
    logical, dimension(:),   allocatable    :: weightchain

    logical :: isHomopolymer
    double precision, dimension(:),allocatable :: lsegseq  ! only needed for copolymer

contains
  
    subroutine allocate_chains(cuantas,nseg,nsegtypes,ngr_node)

        implicit none
    
        integer, intent(in) :: cuantas,nseg,nsegtypes,ngr_node
    
        !     .. extra 100 in index because of  nchain rotations
        allocate(indexchain(nseg,cuantas+100))
        allocate(indexchain_init(nseg,cuantas+100))
        allocate(exp_energychain(cuantas+100))
        allocate(weightchain(cuantas+100))
        allocate(isAmonomer(nseg)) 
        allocate(type_of_monomer(nseg)) 
        allocate(type_of_monomer_char(nseg))
        allocate(ismonomer_of_type(nseg,nsegtypes)) 
        allocate(ismonomer_chargeable(nsegtypes))

    end subroutine allocate_chains
  
end module chains
