module matrices 

    implicit none

    real*8 :: sitheta,cotheta,siphip,cophip
    real*8 :: tt(3,3),tp(3,3),tm(3,3)      
    
    contains

        subroutine init_matrices

            use mathconst

            !      common /mathconst/ pi
            !      common /matrices/ tt,tp,tm
            !      common /senos/ sitheta,cotheta,siphip,cophip

            !      pi=acos(-1.0000000e0)
            sitheta=sin(68.0*pi/180.0)
            cotheta=cos(68.0*pi/180.0)
            siphip=sin(120.0*pi/180.0)
            cophip=cos(120.0*pi/180.0)

            tt(1,1)=cotheta
            tt(1,2)=sitheta
            tt(1,3)=0.0
            tt(2,1)=sitheta
            tt(2,2)=-cotheta
            tt(2,3)=0.0
            tt(3,1)=0.0
            tt(3,2)=0.0
            tt(3,3)=-1.0

            tp(1,1)=cotheta
            tp(1,2)=sitheta
            tp(1,3)=0.0
            tp(2,1)=sitheta*cophip
            tp(2,2)=-cotheta*cophip
            tp(2,3)=siphip
            tp(3,1)=sitheta*siphip
            tp(3,2)=-cotheta*siphip
            tp(3,3)=-cophip

            tm(1,1)=cotheta
            tm(1,2)=sitheta
            tm(1,3)=0.0
            tm(2,1)=sitheta*cophip
            tm(2,2)=-cotheta*cophip
            tm(2,3)=-siphip
            tm(3,1)=-sitheta*siphip
            tm(3,2)=cotheta*siphip
            tm(3,3)=-cophip

        end subroutine init_matrices

end module matrices
