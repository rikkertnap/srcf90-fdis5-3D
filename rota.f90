! ------------------------------------------------------------|
! rotates a given chains conformation                         |  
! pre: xend = input chain                                     |
! post: xendr= rotated chain                                  |
!       return= is_postive                                    |
! ------------------------------------------------------------|

module chain_rotation

  implicit none

contains  

!Euler rotation Z1X2Z3

function rotation(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac=rands(seed)
    fac1=rands(seed)
    fac2=rands(seed)
    alfa=fac*2.0_dp*pi
    cbe=fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)
        xendr(1,i)=a*(-cbe*sal*sga+cal*cga) -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
        xendr(2,i)=a*(cbe*cal*sga+sal*cga)+b*(cbe*cal*cga-sal*sga)-c*sbe*cal
        xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe 
    enddo
  
    is_positive_z=.true.
    i=2 

    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(1,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotation


! rotation only around x-axis
! Euler rotatation X1Z2X3 

function rotationXaxis(xend,xendr,nseg)result(is_positive_z)
  
    use random
    use mathconst
    
    implicit none

    integer,  intent(in)    :: nseg  
    real(dp), intent(in)    :: xend(:,:)
    real(dp), intent(inout) :: xendr(:,:)  
    logical                 :: is_positive_z  ! output 
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac = rands(seed)
    fac1= rands(seed)
    fac2= rands(seed)

    alfa=fac*2.0_dp*pi
    cbe=1.0_dp ! fac1*2.0_dp-1.0_dp
    gama=fac2*2.0_dp*pi
  
    sbe=(1.0_dp-cbe**2)**0.5_dp
    cal=cos(alfa)
    sal=sin(alfa)
    cga=cos(gama)
    sga=sin(gama)
  
    do i=1,nseg+1     ! rotation of xend result stored in xendr
        a=xend(1,i)
        b=xend(2,i)
        c=xend(3,i)

        xendr(1,i)=a ! a*cbe + b*(-cga*sbe)+c*sbe*sga
        xendr(2,i)=a*cal*sbe +b*(cal*cbe*cga-sal*sga) +c*(-cga*sal-cal*cbe*sga)
        xendr(3,i)=a*sbe*sga +b*(cal*sga+cbe*cga*sal) +c*( cal*cga-cbe*sal*sga)
    enddo
  
    is_positive_z=.true.
    i=2 
    do while((i.le.(nseg+1)).and.(is_positive_z)) 
        if(xendr(3,i)<=0.0_dp) is_positive_z=.false.
        i=i+1
    enddo

end function rotationXaxis

end module chain_rotation
    
