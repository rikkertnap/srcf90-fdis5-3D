module lateral_Rgsqr
        
        use precision_definition
        implicit none

contains 

! Computes the lateral radius of gyration as a function of the distance from the tethering surface. 
! Equation from: Szleifer and Carignano, Macromolecules(1995) 'On the structure and presure of tethered polymer layers in good solvent

!Note2Self: Created a separate module rather than incorporate the function into chaingenerator.f90 because of the order in which the files are compiled. 

function calc_lateral_Rgsqr(conf) result (Rgsqr_lateral)
     
    use globals, only: nseg, cuantas
    use volume, only: delta, nz
    use chains, only: indexchain
    use myutils, only : error_handler
    use volume, only: CoordinateFromLinearIndex

    ! Input
    integer, intent(in) :: conf

    ! Output
    real(dp) :: Rgsqr_lateral(nz)

    ! Local variables
    integer :: c, s, idx, x, y, z, nseg_z(nz), z_max
    real(dp):: x_cm(nz), y_cm(nz)
    character(20) :: text

    ! Initialize arrays
    x_cm = 0.0_dp
    y_cm = 0.0_dp
    nseg_z = 0.0_dp
    Rgsqr_lateral = 0.0_dp
    z_max = 0

    ! Step 1: Compute center of mass for each z-layer
    do s = 1, nseg
        idx = indexchain(s, conf)
        call CoordinateFromLinearIndex(idx, x, y, z)

        x_cm(z) = x_cm(z) + x
        y_cm(z) = y_cm(z) + y
        nseg_z(z) = nseg_z(z) + 1
        
        z_max = max(z_max, z)  ! Track the highest z value    
    enddo

    ! Normalize center of mass
    do z = 1, z_max
        if (nseg_z(z) > 0) then
            x_cm(z) = x_cm(z) / nseg_z(z)
            y_cm(z) = y_cm(z) / nseg_z(z)
        endif
    enddo

    ! Step 2: Compute lateral radius of gyration for each z-layer
    do s = 1, nseg
        idx = indexchain(s, conf)
        call CoordinateFromLinearIndex(idx, x, y, z)

        Rgsqr_lateral(z) = Rgsqr_lateral(z) + (x - x_cm(z))**2 + (y - y_cm(z))**2
    enddo

    ! Final normalization
    do z = 1, z_max
        if (nseg_z(z) > 0) then
            Rgsqr_lateral(z) = (Rgsqr_lateral(z) / nseg_z(z)) * (delta**2)
        endif
    enddo

 ! check z_max < nz
      if(z_max > nz) then
         text="Error: z_max != nz"
         call error_handler(1,text)
      endif
       end function calc_lateral_Rgsqr

end module lateral_Rgsqr
