module gizmo_hdf5_input
    use hdf5
    implicit none
    
    contains
    
    subroutine open_hdf5(filename)
        implicit none
        character*512 :: filename
        integer :: hdf5_err
        
        call HFOPEN_F(hdf5_err)
    end subroutine open_hdf5
    
end module gizmo_hdf5_input
