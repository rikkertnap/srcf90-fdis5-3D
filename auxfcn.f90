! auxilari functions used in fcndipolar

module auxfcn

    use precision_definition
    implicit none
 
    real(dp), parameter :: xeps       = 1.0e-8_dp  ! tolerance for function argement to be one for sinhc
    real(dp), parameter :: xepslan    = 1.0e-4_dp  ! tolerance for Langevin function evalution 

    private  :: xeps ,xepslan

contains  
   

    ! sinhc(x)=sinh(x)/x
    ! sinhc(x)=1.0+(x**2)/6.0 +(x**4)/120.0+(x**6)/5040 ! Taylor truncation error x^8
    ! for x<= xeps  x**2 < 1.0e-16 hence more then 16 signification number are need to represent value sinhc
    function sinhc(x) result(fval)

       implicit none

        real(dp), intent(in) :: x
        real(dp) :: fval
    
        if(x>xeps) then
            fval=sinh(x)/x
        else if(x<=xeps) then
            fval=1.0_dp
        endif        

    end function sinhc

    ! Evaluate L(x) by backward recursion of the continued fraction
    ! L(x) = x / ( 3 + x^2 / ( 5 + x^2 / ( 7 + x^2 / ( 9 + ... ) ) ) )  
    ! see https://en.wikipedia.org/wiki/Brillouin_and_Langevin_functions
    ! see https://www.mathworks.com/matlabcentral/fileexchange/38405-langevin-function-accurate-evaluation
    ! L(x)=x/3.0-(x**3)/45.0 +2.0*(x**5)/945.0 -(x**7)/4725.0 +2.0_dp*(x**9)/(93555.0) ! Taylor truncation error x^11
    ! Taylor expansion osccilation cause significant  number loss

    function Langevin_continuedfrac(x) result(fval)

        real(dp), intent(in) :: x
        real(dp) :: fval 

        ! local variables
        real(dp) :: xsqr, frac, bk 
        integer :: Nlevel, k

        Nlevel = 10     ! default truncation level for continued fraction
        xsqr   = x * x  ! square of x
        
        frac=0.0_dp

        do k = Nlevel, -1, 2          ! backward recursion  from k=N to k=2
            bk = 2.0_dp * k + 1.0_dp 
            frac = xsqr/ ( bk + frac ) 
        enddo

        fval = x / ( 3.0_dp + frac )   ! last step k=1

     end function



    ! langevin function L(x) =coth(x) -1/x=1/tanh(x) -1/x
   
    function Langevin(x) result(fval)

        real(dp), intent(in) :: x
        real(dp) :: fval

        if(x>xepslan) then

            fval=(1.0_dp/tanh(x)) -1.0_dp/x
        
        else if(x>0.0_dp.and.x<=xepslan) then
        
            fval=Langevin_continuedfrac(x)
        
        else  ! x=0
        
            fval=0.0_dp   
        
        endif        

    end function Langevin

    ! returns 1 if x>0 and -1 otherwise
    ! gives sign of direction in x-direction x being the vector element 
     
    function sign_fnc(x) result(ival)

        real(dp), intent(in) :: x
        integer :: ival
        if(x>0) then 
            ival=1
        else
            ival=-1     
        endif
    end function sign_fnc

    ! tests evaluation of sinhc for selected values of x 
    subroutine unit_test_sinhc

        use myutils, only : newunit

        real(dp):: xval,deltax
        integer :: i
        character(len=9) :: fname
        integer :: ios, un
        
        write(fname,'(A9)')'sinhc.dat'
        open(newunit=un,file=fname,iostat=ios)   

        xval=1.0e-18_dp
        do i=0,20
            xval=xval*10
            write(un,*)xval,sinhc(xval),sinh(xval)/xval 
        enddo    
        deltax=1.0e-9_dp
        do i=1,30
            xval=deltax*(i)
            write(un,*)xval,sinhc(xval),sinh(xval)/xval
        enddo   
        
        close(un)

    end subroutine


    ! tests evaluation of Langevin fnc for selected values of x 
    subroutine unit_test_Langevin

        use myutils, only : newunit

        real(dp):: fval,xval,deltax
        integer :: i
        character(len=12) :: fname        
        integer :: ios, un
        
        write(fname,'(A12)')'Langevin.dat'
        open(newunit=un,file=fname,iostat=ios)

        xval=1.0e-18_dp
        do i=0,20
            xval=xval*10
            write(un,*)xval,Langevin(xval),((1.0_dp/tanh(xval)) -1.0_dp/xval),(1.0_dp/tanh(xval)), -1.0_dp/xval 
        enddo    
        deltax=1.0e-9_dp
        do i=1,30
            xval=deltax*(i)
            write(un,*)xval,Langevin(xval),(1.0_dp/tanh(xval)) -1.0_dp/xval,1.0_dp/tanh(xval),-1.0_dp/xval 
        enddo   
        
        close(un)

    end subroutine


end module auxfcn
