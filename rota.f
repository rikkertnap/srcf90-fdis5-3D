      subroutine rota(xend,xendr,n,test)

      implicit none
     
      integer n

      real*8 xend(3,n+5),rands,xendr(3,n+5)
      real*8 delta,radius

      real*8 theta,theta1,x(120),y(120),z(120),xp(120),yp(120),pi
   
      character*1 test
      integer seed,k,m1,i
      real*8 fac,fac1,fac2,sbe,cbe,sal,cal,sga,dist
      real*8 alfa,gama,cga,a,b,c

      common /seedrandom/ seed
      common /layer/ delta,radius  ! see  volume.h

      pi=acos(-1.0)

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

      do i=1,n+1     ! rotation of xend result stored in xendr

         a=xend(1,i)
         b=xend(2,i)
         c=xend(3,i)

         xendr(1,i)=a*(-cbe*sal*sga+cal*cga)
     &        -b*(cbe*sal*cga+cal*sga)+c*sbe*sal
         xendr(2,i)=a*(cbe*cal*sga+sal*cga)+
     &        b*(cbe*cal*cga-sal*sga)-c*sbe*cal
         xendr(3,i)=a*sbe*sga+b*sbe*cga+c*cbe
 
      enddo


      do i=2,n+1

         if(xendr(1,i).le.0) then
            test='N'  
            go to 2
         endif

      enddo

 2    return
      end
      
