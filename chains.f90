!     .. module file of chains variables
module chains
  
    use globals
    implicit none
  
    integer*2, dimension(:,:), allocatable :: indexchainAB 
    integer*2, dimension(:,:), allocatable :: indexchainC 
    integer*2, dimension(:,:), allocatable :: indexchainAB_init 
    integer*2, dimension(:,:), allocatable :: indexchainC_init 
    logical, dimension(:), allocatable :: isAmonomer

contains
  
    subroutine allocate_chains(cuantasAB,nsegAB,cuantasC,nsegC)
        implicit none
    
        integer, intent(in) :: cuantasAB,nsegAB,cuantasC,nsegC
    
        !     .. extra 100 in index because of  nchain rotations

        allocate(indexchainAB(cuantasAB+100,nsegAB))
        allocate(indexchainC(cuantasC+100,nsegC))
        allocate(indexchainAB_init(cuantasAB+100,nsegAB))
        allocate(indexchainC_init(cuantasC+100,nsegC))

        allocate(isAmonomer(nsegAB))
    
    end subroutine allocate_chains
  
end module chains
