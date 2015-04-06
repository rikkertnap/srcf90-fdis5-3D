 PROGRAM test_trim
   implicit none
   real :: num 
   integer :: i
   character(len=3) :: istr
   character(len=8) :: rstr
   character(len=100):: outfile
 
   





   num=13.456
   write(rstr,"(F8.2)")num
   write(*,*) trim(adjustl(rstr)), rstr
   
  
END PROGRAM
