module molecules

    implicit none
	
	type moleclist
        real*8 :: sol
        real*8 :: Na
        real*8 :: Cl
        real*8 :: K
        real*8 :: Ca
        real*8 :: NaCl
        real*8 :: KCl
        real*8 :: Hplus
        real*8 :: OHmin
    end type moleclist

end module 