! ------------------------------------------------------------|
! rotates a given chains conformation                         |  
! pre: xend = input chain                                     |
! post: xendr= rotated chain                                  |
!       return= is_postive                                    |
! ------------------------------------------------------------|


subroutine rotation(xend,xendr,nseg,is_positive_z,lseg)
  
    use random
    use mathconst
    
    implicit none

    integer :: nseg  
    logical :: is_positive_z
    real(dp)  :: lseg
    real(dp)  :: xend(3,nseg+5),xendr(3,nseg+5)  
  
    ! ..local argument 

    integer :: i
    real(dp)  :: fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
    real(dp)  :: alfa,gama,cga,a,b,c
    
    fac=rands(seed)
    fac1=rands(seed)
    fac2=rands(seed)
    alfa=fac*2*pi
    cbe=fac1*2.0-1.0
    gama=fac2*2*pi
  
    sbe=(1-cbe**2)**0.5
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
         if(xendr(1,i)<=0) is_positive_z=.false.
         i=i+1
    enddo

end subroutine rotation


    
