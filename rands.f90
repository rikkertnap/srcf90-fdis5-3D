!     .. modules for random number generator
 
module random
      
    implicit none
      
    integer :: seed           ! seed of random generator 
  
    contains 
  
    double precision function rands (SEED)
    
        integer :: SEED
        integer :: I1,I2
    
        I2 = MOD (SEED, 32768) * 30485
        I1 = MOD(SEED/ 32768 * 30485,65536) + MOD(MOD(SEED,32768)*48603,65536) + I2/32768 + MOD(seed/32768,2) * 32768 +13849   
        I2 = MOD (I2, 32768) + 12659
        I1 = I1 + I2 / 32768
        I2 = MOD (I2, 32768)
        I1 = MOD (I1, 65536)
        seed = I1 * 32768 + I2
        rands = REAL(I1 * 32768 + I2) * 4.65661287308E-10
    
        return
    end function rands
  
end module random

  

