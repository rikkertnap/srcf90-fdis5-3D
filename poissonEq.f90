module Poisson  

    use precision_definition

    implicit none

    private :: ipbc 

contains

    subroutine Poisson_Equation(fvec,psi,rhoq,sigmaqSurfR,sigmaqSurfL)

        use globals, only : nsize, neq
        use parameters, only : constqW
        use volume, only : nx,ny,nz, linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout) :: fvec(:)
        real(dp), intent(in) :: psi(:)
        real(dp), intent(in) :: rhoq(:)
        real(dp), intent(in) :: sigmaqSurfR(:)
        real(dp), intent(in) :: sigmaqSurfL(:)
        
        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: idxR, idxL, id2D  
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin

        ! .. electrostatics 
       
        noffset=nsize

        ! do iz=1,nz
        !     do iy=1,ny
        !         do ix=1,nx
        !             call linearIndexFromCoordinate(ix,iy,iz  ,id)
        !             print*,"(",ix,iy,iz,") index ",id 
        !         enddo    
        !     enddo
        ! enddo        



        do ix=1,nx
            do iy=1,ny
                do iz=2,nz-1
                    call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                    call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                    call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                    call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                    call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                    call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                    call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

                    fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls)+psi(idzmin) &
                        -6.0_dp*psi(id) +rhoq(id)*constqW)
                enddo
            enddo
        enddo    

        ! boundary iz=1 

        do ix=1,nx
            do iy=1,ny
                iz=1
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz+1,idzpls)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)

            !    print*,ix,iy,iz,id,idypls,ipbc(iy+1,ny)

                fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzpls) +sigmaqSurfL(id) &
                    - 5.0_dp*psi(id) +rhoq(id)*constqW)
            
            enddo
        enddo    

        ! boundary iz=nz 

        do ix=1,nx
            do iy=1,ny
                iz=nz
                call linearIndexFromCoordinate(ix,           iy,iz  ,id)
                call linearIndexFromCoordinate(ipbc(ix+1,nx),iy,iz  ,idxpls)
                call linearIndexFromCoordinate(ipbc(ix-1,nx),iy,iz  ,idxmin)
                call linearIndexFromCoordinate(ix,           iy,iz-1,idzmin)
                call linearIndexFromCoordinate(ix,ipbc(iy+1,ny),iz  ,idypls)
                call linearIndexFromCoordinate(ix,ipbc(iy-1,ny),iz  ,idymin)
                
                id2D=id-(nsize-nx*ny)
                
                fvec(noffset+id)= -0.5_dp*( psi(idxpls)+psi(idxmin) +psi(idypls)+psi(idymin)+psi(idzmin) +sigmaqSurfR(id2D) &
                    -5.0_dp*psi(id) +rhoq(id)*constqW)
            enddo
        enddo    


    end subroutine


    
    subroutine Poisson_Equation_Surface(fvec,psi,rhoq,psisurfR,psisurfL,sigmaqSurfR,sigmaqSurfL,bcflag)

        use globals, only : nsize, neq, LEFT, RIGHT, sysflag
        use parameters, only : constqW
        use volume, only : nx,ny,nz, linearIndexFromCoordinate

        implicit none

        ! input arguments 
        real(dp), intent(inout)  :: fvec(:)
        real(dp), intent(in)  :: psi(:)
        real(dp), intent(in)  :: rhoq(:)
        real(dp), intent(inout) :: psisurfR(:)
        real(dp), intent(inout) :: psisurfL(:)
        real(dp), intent(in)  :: sigmaqSurfR(:)
        real(dp), intent(in)  :: sigmaqSurfL(:)
        character(len=2), intent(in)  :: bcflag(2)

        ! local variables
        integer :: ix, iy, iz, noffset
        integer :: idxR, idxL, idxR2D
        integer :: id, idxpls, idxmin, idypls, idymin, idzpls, idzmin
        integer :: neq_bc

        ! .. electrostatics: self consistent boundary conditions
        
        if(sysflag=="elect") then 
            noffset=4*nsize
        else if(sysflag=="electdouble") then 
            noffset=4*nsize
        else if(sysflag=="electnopoly") then 
            noffset=2*nsize
        else 
            print*,"error: sysflag wrong value for Poisson_equation_surface "    
        endif     

        neq_bc=0

        if(bcflag(RIGHT)/='cc') then
            neq_bc=nx*ny
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny) ! check this 
                    fvec(noffset+idxR2D)=psisurfR(idxR2D)-psi(idxR)-sigmaqSurfR(idxR2D)/2.0_dp
                enddo
            enddo
        else
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,nz,idxR)
                    idxR2D=idxR-(nsize-nx*ny)
                    psisurfR(idxR2D) = psi(idxR)+sigmaqSurfR(idxR2D)/2.0_dp 
                enddo
            enddo
        endif    

        noffset=noffset +neq_bc

        if(bcflag(LEFT)/='cc') then 
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    fvec(noffset+idxL)=psi(idxL)-psisurfL(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        else    
            do ix=1,nx
                do iy=1,ny
                    call linearIndexFromCoordinate(ix,iy,1,idxL)
                    psisurfL(idxL) = psi(idxL)+sigmaqSurfL(idxL)/2.0_dp
                enddo
            enddo
        endif   

    end subroutine

        


        

    function ipbc(ival,imax) result(intpbc)
        implicit none 
        integer, intent(in) :: ival
        integer, intent(in) :: imax
        integer :: intpbc

        if(ival>0) then
            intpbc=ival-int((ival-1)/imax)*imax
        else
            intpbc=ival-(int((ival-1)/imax)-1)*imax
        endif

    end function


end module