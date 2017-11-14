! real precision  dp equal to machine defined double precision

module precision_definition
  
   implicit none

   integer,  parameter :: sp=kind(1.0)
   integer,  parameter :: dp=kind(1.0D0)

end module precision_definition


! real precision  sp,dp and qp defined excplicitly 

module precision_definition_def

  implicit none

  integer , parameter :: sp=selected_real_kind(6,37)                                     
  integer , parameter :: dp=selected_real_kind(15,307)       
  integer , parameter :: qp=selected_real_kind(33,4931)        

end module precision_definition_def


! real precision defined from fortran iso 2008 definitions

module precision_definition_f2008 
 
  use, intrinsic :: iso_fortran_env
  
  implicit none
  integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
   
 
end module precision_definition_f2008
