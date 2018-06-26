!     .. module file of chains variables
module chains
  
    use globals
    implicit none
  
    integer, dimension(:,:,:), allocatable :: indexchainAB 
    integer, dimension(:,:,:), allocatable :: indexchainC 
    integer, dimension(:,:,:), allocatable :: indexchainAB_init 
    integer, dimension(:,:,:), allocatable :: indexchainC_init 
    logical, dimension(:), allocatable :: isAmonomer
    logical, dimension(:,:), allocatable :: weightchainAB

    logical :: isHomopolymer
    double precision, dimension(:),allocatable :: lsegseq  ! only needed for copolymer

contains
  
    subroutine allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC,ngr_node)
        implicit none
    
        integer, intent(in) :: cuantasAB,nsegAB,cuantasC,nsegC,ngr_node
    
        !     .. extra 100 in index because of  nchain rotations

        allocate(indexchainAB(nsegAB,ngr_node,cuantasAB+100))
        allocate(indexchainAB_init(nsegAB,ngr_node,cuantasAB+100))
        allocate(isAmonomer(nsegAB)) 
        allocate(weightchainAB(ngr_node,cuantasAB+100))

        allocate(indexchainC(nsegC,ngr_node,cuantasC+100))
        allocate(indexchainC_init(nsegC,ngr_node,cuantasC+100))
    
    end subroutine allocate_chains
  
end module chains
