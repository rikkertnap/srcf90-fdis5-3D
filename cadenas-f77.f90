!-----------------------------------------------------------|
! Generates a chains conformations bases                    |
! on a three state RIS-model see Flory book                 |
!                                                           | 
!-----------------------------------------------------------|

      subroutine cadenas(chains,nseg,lseg,ncha)
 
      use random
      use mathconst
      use matrices

      implicit none
	
      integer nseg              ! number of segments on sphere
      real*8 lseg               ! lenght segment in nm
      real*8 chains(3,nseg,200) ! 
      integer ncha              ! number of rotations

      integer i,state,ii,j,k1,k2
      real*8 rn,dista,angle
      real*8 m(3,3),mm(3,3) 
      real*8 x(3),xend(3,nseg+5),xendr(3,nseg+5)
      character*1 test
      logical is_positive_z


 223  xend(1,1)=0.0
      xend(2,1)=0.0
      xend(3,1)=0.0

      rn=rands(seed)

      angle=0.0
      
      m(1,1)=cotheta
      m(1,2)=sitheta
      m(1,3)=0.0
      m(2,1)=cos(angle)*sitheta
      m(2,2)=-cos(angle)*cotheta
      m(2,3)=sin(angle)
      m(3,1)=sin(angle)*sitheta
      m(3,2)=-sin(angle)*cotheta
      m(3,3)=-cos(angle)
      
      x(1)=m(1,1)*lseg
      x(2)=m(2,1)*lseg
      x(3)=m(3,1)*lseg
      
      xend(1,2)=xend(1,1)+x(1)
      xend(2,2)=xend(2,1)+x(2)
      xend(3,2)=xend(3,1)+x(3)

      do 10 i=3,nseg+1
         
 123     rn=rands(seed)
         state=int(rn*3)

         if (state.eq.3) then 
            state=2
         endif
         if (state.eq.0) then
!*********************************** TRANS     

            call mrrrr(m,tt,mm)
            do 30 ii=1,3
               do 40 j=1,3
                  m(ii,j)=mm(ii,j)
 40            continue
 30         continue


         elseif (state.eq.1) then

!********************************** GAUCHE +

            call mrrrr(m,tp,mm)
            do 31 ii=1,3
               do 41 j=1,3
                  m(ii,j)=mm(ii,j)
 41            continue
 31         continue
            


         elseif (state.eq.2) then

!********************************** GAUCHE -

            call mrrrr(m,tm,mm)
            do 32 ii=1,3
               do 42 j=1,3
                  m(ii,j)=mm(ii,j)
 42            continue
 32         continue
            
         endif
         
         x(1)=m(1,1)*lseg
         x(2)=m(2,1)*lseg
         x(3)=m(3,1)*lseg
         
         xend(1,i)=xend(1,i-1)+x(1)
         xend(2,i)=xend(2,i-1)+x(2)
         xend(3,i)=xend(3,i-1)+x(3)
         
 10   continue

      dista=0.0
      do k1=4,nseg+1
         do k2=1,k1-3
            dista=(xend(1,k2)-xend(1,k1))**(2.0)
            dista=dista+(xend(2,k2)-xend(2,k1))**(2.0)
            dista=dista+(xend(3,k2)-xend(3,k1))**(2.0)
            dista=sqrt(dista)
            if (dista.lt.lseg) then
!c               print*,'no self-a.',k1,k2
               goto 223
            endif
         enddo
      enddo

     
      ncha=0
      do 400 i=1,72
!         test='S'
         call rotation(xend,xendr,nseg,is_positive_z,lseg)
!        if (test.eq.'N') goto 400
         if (is_positive_z.eqv..false.) goto 400
         ncha=ncha+1
         do 401 j=1,nseg
            chains(1,j,ncha)=xendr(1,j+1)
            chains(2,j,ncha)=xendr(2,j+1)
            chains(3,j,ncha)=xendr(3,j+1)
 401     continue
         if (ncha.eq.12) goto 402
      
 400  continue
 402  if (ncha.eq.0) goto 223

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        subroutine mrrrr(a,b,c)

        real*8 a(3,3),b(3,3),c(3,3)

        do 1 i=1,3

           do 1 j=1,3

              c(i,j)=0

 1      continue


        do 2 i=1,3

           do 2 j=1,3

              do 2 k=1,3

                 c(i,j)=c(i,j)+a(i,k)*b(k,j)

 2      continue

        return
        end 

