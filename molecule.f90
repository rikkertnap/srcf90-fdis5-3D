module molecules

  use precision_definition
  implicit none
	
  type moleclist
     real(dp) :: sol
     real(dp) :: Na
     real(dp) :: Cl
     real(dp) :: K
     real(dp) :: Ca
     real(dp) :: NaCl
     real(dp) :: KCl
     real(dp) :: Hplus
     real(dp) :: OHmin
  end type moleclist
  
end module 
