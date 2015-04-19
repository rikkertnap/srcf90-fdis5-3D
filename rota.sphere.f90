!------------------------------------------------------------|
! rotates a given chains conformation                        |  
! pre: xend = input chain                                    |
! post: xendr= rotated chain                                 |
!       return= test                                         |
!  (02-04-2014)                                              |
! TODO: remove goto statement and make into function         |
!------------------------------------------------------------|


subroutine rota(xend,xendr,nseg,test,lseg)
  
  use random
  use mathconst
  use volume, only : delta,radius
  
  implicit none
  
  integer :: nseg
  real(dp) :: xend(3,nseg+5),xendr(3,nseg+5)
  real(dp) :: lseg
  
  ! ..local argument 

  real(dp) :: theta,theta1,x(120),y(120),z(120),xp(120),yp(120)
  
  character(len=1) :: test
  integer :: k,m1,i
  real(dp) :: fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
  real(dp) :: alfa,gama,cga,a,b,c
   

  theta1=pi/6.0
  
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
  
  do m1=1,1    ! further rotation and translation to 
     ! grafted onto sphere 
     
     theta=m1*theta1
     
     do k=2,nseg+1
        
        x(k)=xendr(2,k)*cos(theta)+xendr(3,k)*sin(theta)+radius*sin(theta)
        y(k)=-xendr(2,k)*sin(theta)+xendr(3,k)*cos(theta)+radius*cos(theta)
        z(k)=xendr(1,k)
     enddo
     
     do i=2,nseg+1
        
        if((x(i)**2+y(i)**2 + z(i)**2)**0.5.lt.radius) then
           test='N'  
           !               print*,test,i
           go to 2
        endif
        
     enddo
     
  enddo
  
2 return
end subroutine rota


    
