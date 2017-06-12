! read gadget ICs in Fortran
! because C is a load of nonsense

module gadget_ics_IO
    type header_t
        integer(kind=4), dimension(6) :: np ! number of particles of each type in this file
        real(kind=8), dimension(6) :: mp_g ! mass of each group of particles
        real(kind=8) :: time, redshift
        integer(kind=4) :: flagSFR, flagFeed
        integer(kind=4), dimension(6) :: Nall ! total number of particles of each type in all files
        integer(kind=4) :: flagCooling, numFiles
        real(kind=8) :: boxSize, omega0, omegaLambda, hubbleParam
        integer(kind=4) :: flagAge, flagMetals
        integer(kind=4) :: hashtabsize
        
        character*84 :: fill
        
    end type header_t
    type p_data_t
        integer :: np_tot,ng
        integer :: Nm ! number of particles with masses
    
        real(kind=4), allocatable, dimension(:,:) :: r_p, v_p
        real(kind=4), allocatable, dimension(:) :: m_p, u_p
        integer(kind=4), allocatable, dimension(:) :: id_p
        
        contains
        
        procedure :: doAllocations => p_doAllocations
        procedure :: sizesFromInput => p_sizesFromInput

        
    end type p_data_t

    
    type(header_t) :: header
    type(p_data_t) :: p_data
    
    
    contains
    
    subroutine p_sizesFromInput(this)
        implicit none
        class(p_data_t) :: this
        integer :: ig ! group loop variable
    
        this%Nm = 0
        do ig=1,6
            if ( header%mp_g(ig)==0 ) then
                this%Nm = this%Nm + header%np(ig)
            endif
        end do

        this%np_tot = sum(header%np)
        this%ng = header%np(1)
    end subroutine p_sizesFromInput

    subroutine p_doAllocations(this)
        implicit none
        class(p_data_t) :: this

        allocate(this%r_p(3,this%np_tot))
        allocate(this%v_p(3,this%np_tot))

        allocate(this%m_p(this%Nm))
    
        allocate(this%u_p(this%ng))

        allocate(this%id_p(this%np_tot))
    end subroutine p_doAllocations
    
    subroutine read_data(filename)
        implicit none
        character(len=*) :: filename

        open(unit=10,file=filename,form="unformatted")
        read(10) header
    
        call p_data%doAllocations
    
        read(10) p_data%r_p
        read(10) p_data%v_p
        read(10) p_data%id_p

        if ( p_data%Nm>0 ) then
            read(10) p_data%m_p
        endif
        
        read(10) p_data%u_p
    
        close(10)
    end subroutine read_data
    
    subroutine write_ics(filename)
        implicit none
        character(len=*) :: filename
        
        open(unit=10,file=filename,form="unformatted")
        write(10) header
        write(10) p_data%r_p
        write(10) p_data%v_p
        write(10) p_data%id_p
        write(10) p_data%u_p
        close(10)

    
    end subroutine write_ics

end module gadget_ics_IO
