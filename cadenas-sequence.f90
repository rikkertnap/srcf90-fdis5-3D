! ----------------------------------------------------------|
! Generates a chains conformations bases                    |
! on a three state RIS-model see Flory book                 |
!                                                           | 
! ----------------------------------------------------------|

module cadenas_sequence

  implicit none

  private

  double precision :: lsegAA ,lsegBB, lsegAB
  
  public   :: make_linear_seq_chains, make_lsegseq

contains  

subroutine make_linear_seq_chains(chain,nchains,maxnchains,nseg)
  
  use mathconst
  use random
  use matrices
  use chain_rotation, only : rotation
  use chains, only : lsegseq

  implicit none
  
  !     .. scalar arguments

 
  integer :: nseg               ! number of segments on sphere
  integer :: nchains            ! number of rotations
  integer :: maxnchains

  !  .. array arguments
  
  double precision  :: chain(3,nseg,200)  

  !  .. local variables
  
  
  
  integer :: i,j,state
  double precision :: rn, angle
  double precision, dimension(3,3) ::  m, mm
  double precision, dimension(3) ::  x
  double precision, dimension(3,nseg+5) :: xend, xendr
  integer :: maxattempts
  logical :: is_selfavoid,is_positive_rot,is_positive_z
  character(len=1) :: test

  
  !     .. executable statements 

  maxattempts =72 ! max number of attempts to rotate a chain
 
  do while(nchains==0) ! make at least one chain 
     
     is_selfavoid=.FALSE.

     do while(is_selfavoid.eqv..FALSE.)

        xend(1,1)=0.0
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
        
        x(1)=m(1,1)*lsegseq(1)
        x(2)=m(2,1)*lsegseq(1)
        x(3)=m(3,1)*lsegseq(1)
        
        xend(1,2)=xend(1,1)+x(1) ! second position 
        xend(2,2)=xend(2,1)+x(2)
        xend(3,2)=xend(3,1)+x(3)
        
        do i=3,nseg+1
     
           rn=rands(seed)
           state=int(rn*3)
     
           if (state.eq.3) then 
              state=2
           endif
           if (state.eq.0) then   ! trans 
              call mrrrr(m,tt,mm)
           elseif (state.eq.1) then ! gauche plus 
              call mrrrr(m,tp,mm)
           elseif (state.eq.2) then ! gauche minus
              call mrrrr(m,tm,mm)
           endif
               
           x(1)=m(1,1)*lsegseq(i-1)
           x(2)=m(2,1)*lsegseq(i-1)
           x(3)=m(3,1)*lsegseq(i-1)
           
           xend(1,i)=xend(1,i-1)+x(1)
           xend(2,i)=xend(2,i-1)+x(2)
           xend(3,i)=xend(3,i-1)+x(3)
         
        enddo
        
        is_selfavoid=selfavoidance(xend,nseg)
    
     enddo
     
     i=1

     do while((i.le.maxattempts).and.(nchains.lt.maxnchains)) 
      
        is_positive_rot=rotation(xend,xendr,nseg) 
        if (is_positive_rot) then
           nchains=nchains+1
           do j=1,nseg
              chain(1,j,nchains)=xendr(1,j+1)
              chain(2,j,nchains)=xendr(2,j+1)
              chain(3,j,nchains)=xendr(3,j+1)
           enddo
        endif

        i=i+1

     enddo
    
  enddo

end subroutine make_linear_seq_chains


! check self avoidance 
! pre  xend 
! post return selfavoidance=true or false(=intersection)

logical function selfavoidance(xend,nseg) 
 
  use chains, only  : isAmonomer

  implicit none

  integer :: nseg
  double precision :: xend(3,nseg+5)
   
  ! .. local vars 
  double precision :: dista
  integer :: k,l
  logical :: selfavoid, is_k, is_l
  
  dista=0.0
  selfavoid=.TRUE.
  

  do k=4,nseg+1
     is_k=isAmonomer(k-1) 
     do l=1,k-3
        is_l=isAmonomer(l)
        dista=(xend(1,l)-xend(1,k))**(2.0)
        dista=dista+(xend(2,l)-xend(2,k))**(2.0)
        dista=dista+(xend(3,l)-xend(3,k))**(2.0)
        if (dista.lt.segment_lengthsqr(is_l,is_k))then
           selfavoid=.FALSE.
           !print*,'no self-avoidance.',k,l,"dista=",dista
        endif
     enddo
  enddo

  selfavoidance = selfavoid
  return

end function selfavoidance


subroutine mrrrr(a,b,c)
  
  implicit none
  
  double precision ::  a(3,3),b(3,3),c(3,3)
  integer :: i,j,k
  
  
  do i=1,3
     do  j=1,3
        c(i,j)=0.0
     enddo
  enddo
  
  do i=1,3
     do j=1,3
        do k=1,3
           c(i,j)=c(i,j)+a(i,k)*b(k,j)
        enddo
     enddo
  enddo
  
  do i=1,3
     do  j=1,3
        a(i,j)=c(i,j)
     enddo
  enddo

end subroutine mrrrr


function segment_lengthsqr(s,t)


  implicit none

  logical :: s, t
  double precision :: segment_lengthsqr,lseg2


  if (s) then
     if (t) then
        lseg2  = lsegAA
     else
        lseg2 = lsegAB
     endif
  else
     if (t) then
        lseg2= lsegAB
     else
        lseg2= lsegBB
     endif
  endif  
  segment_lengthsqr=lseg2
  return

end function   


subroutine make_lsegseq(lsegseq,nseg)

    use parameters, only : lsegPAA, lsegPAMPS
    use chains, only  : isAmonomer

    implicit none

    integer :: nseg
    double precision, dimension(nseg) :: lsegseq

    integer :: s

    lsegAA=lsegPAA*lsegPAA
    lsegBB=lsegPAMPS*lsegPAMPS
    lsegAB=((lsegPAA+lsegPAMPS)/2.0)**2


    do s=1,nseg
        if(isAmonomer(s)) then
            lsegseq(s)=lsegPAA
        else
            lsegseq(s)=lsegPAMPS
        endif
    enddo

end subroutine 
  
  


end module    
   