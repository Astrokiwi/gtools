program blobfile_ICs
    use gadget_ics_IO
    implicit none
    
    character(len=256) :: pramfile,intext,invariable,invalue,outfile,trimtext
    integer :: ios
    integer :: whitey

    integer :: np
    
    logical :: valid_command
    
    integer :: ib, n_so_far
    
    real(kind=8) :: mtot,T0,blob_rad
    character(len=256) :: blobloc_file

    real(kind=8) :: blob_m
    integer :: nblobs
    integer :: blob_n
    real(kind=8), dimension(3) :: blob_r,blob_v
    
    if ( command_argument_count()<1 ) then
        print *,"Please give the parameter file in the command line"
        print *,"e.g. ./disc_ICs prams/fatDiscPrams.dat"
        stop
    endif
    
    call get_command_argument(1,pramfile)
    
    
    ! set default values
    
    mtot = 1.d4 ! units of Msun
    
    np = 10**4
    
    T0 = 1.e4 ! in K

    nblobs = 0

    outfile = "data/testout.dat"
    
    ! load in parameters to overload
    open(15,file=pramfile)
    ios = 0
    do while ( ios==0 )
        read(15,"(A)",iostat=ios) intext
        valid_command = .true.
        if ( len_trim(intext)<1 ) then
            valid_command=.false.
        else
            trimtext = trim(intext)
            if ( scan(trimtext,"#")>0 ) then
                valid_command=.false.
            endif
        endif
        if ( ios/=0 ) then
            valid_command=.false.
        endif
        if ( valid_command ) then
            read(intext,*,iostat=ios) invariable,invalue
            if ( ios/=0 ) then
                print *,"INCORRECT FORMAT:"
                print *,intext
                print *,"QUITTING"
                stop
            endif
            select case(invariable)
                case("blobfile")
                    whitey = scan(intext,' ')
                    blobloc_file = adjustl(trim(intext(whitey:)))
                case("outfile")
                    whitey = scan(intext,' ')
                    outfile = adjustl(trim(intext(whitey:)))
                case("T0")
                    read(invalue,*) T0
                case("nblobs")
                    read(invalue,*) nblobs
                case("N")
                    read(invalue,*) np
                case("mtot")
                    read(invalue,*) mtot
                case("blob_rad")
                    read(invalue,*) blob_rad
                    blob_rad=blob_rad*1.e-3 ! kpc to pc
                case default
                    print *,invariable," not understood"
                    
                    
            end select
            if ( ios/=0 ) then
                print *,"INCORRECT FORMAT:"
                print *,intext
                print *,"QUITTING"
                stop
            endif

        endif
    end do
    close(15)
    
    if ( nblobs==0 ) then
        print *,"Must specify at least one blob!"
        stop
    endif
    
    ! units of 1e10 Msun, kpc
    mtot = mtot/1.e10

    blob_n = np/nblobs
    np=blob_n*nblobs
    
    call set_N_onlygas(np)

    ! only mp_g(1) is actually read, but we need to give a value for 
    ! every entry, otherwise gadget/gizmo will try to read an array of 
    ! mass values per particle from the IC file
    header%mp_g(1:6) = mtot/header%np(1)

    blob_m=header%mp_g(1)*blob_n
    
    
    call set_cosmo_noncosmo
    
    call set_flags

    header%boxSize = 0. ! ignored with periodic boundaries turned off

    call p_data%doAllocations
    
    n_so_far = 0
    open(unit=15,file=blobloc_file)
    do ib=1,nblobs
        read(15,*) blob_r,blob_v
        call blobfile_locs(n_so_far,blob_n,blob_r,-blob_v,T0,blob_rad)
        n_so_far = n_so_far + blob_n
    end do
    close(15)
    
    call write_ICs(outfile)
    
end program blobfile_ICs

! set N etc arrays for 1 file IC, with only gas, and all particles
! having the same mass
subroutine set_N_onlygas(ng)
    use gadget_ics_IO
    implicit none
    integer :: ng
    
    ! set N
    p_data%ng = ng
    p_data%np_tot = p_data%ng
    
    header%Nall(1) = p_data%ng
    header%Nall(2:6) = 0
    
    p_data%Nm = 0
    
    ! I/O stuff
    header%numfiles = 1
    header%np(1:6) = header%Nall(1:6)

end subroutine set_N_onlygas

! set cosmological values & time to sensible values
! most of these are ignored because this is a non-cosmo sim
subroutine set_cosmo_noncosmo
    use gadget_ics_IO
    implicit none
    
    ! cosmological stuff that gets ignored
    header%omega0 = 0
    header%omegaLambda = 0
    header%hubbleParam = 0
    header%redshift = 0
    header%time = 0

end subroutine set_cosmo_noncosmo

subroutine set_flags
    use gadget_ics_IO
    implicit none
    
    ! physics stuff that is probably ignored at the moment
    header%flagAge = 0
    header%flagCooling = 0
    header%flagMetals = 0
    header%flagFeed = 0
    header%flagSFR = 0

    header%hashTabSize = 0


end subroutine 

subroutine blobfile_locs(offset, blobn, blob_cent, blob_v, T0,blob_rad)
    use gadget_ics_IO
    implicit none
    
    integer :: offset, blobn
    
    integer :: ip
    real(kind=8) :: T0,blob_rad
    real(kind=8), dimension(3) :: blob_cent, blob_v
    
    real(kind=8), dimension(:,:), allocatable :: rans
    
    real(kind=8), parameter :: pi = 4.d0*atan(-1.d0)
    
    ! gadget/gizmo internal units
    real(kind=8), parameter :: munit_cgs = 1.989d43 ! 10^10 solar masses
    real(kind=8), parameter :: tunit_cgs = 3.08568d16 ! just under a Gyr
    real(kind=8), parameter :: runit_cgs = 3.0857d21 ! 1 kpc
    real(kind=8), parameter :: uunit_cgs = 1.d10 ! 1e10 erg/g
    !real(kind=8), parameter :: G = 6.67408d-8 * tunit_cgs**2 * munit_cgs / runit_cgs**3 ! G=43007.1
    real(kind=8), parameter :: G = 43007.1 ! from Gizmo output, consist with above units
    
    real(kind=8), parameter :: molecular_mass = 4./(1.+3.*.76), proton_mass_cgs = 1.6726d-24
    real(kind=8), parameter :: gamma_minus_one = 5./3.-1., boltzmann_cgs = 1.38066d-16
    real(kind=8), parameter :: u_to_TK = gamma_minus_one/boltzmann_cgs*(molecular_mass*proton_mass_cgs)
    
    ! unit conversions
     ! n.b. tunit/runit is basically 1 km/s anyway
    blob_v = blob_v*1.e5*(tunit_cgs/runit_cgs)
    blob_cent=blob_cent/1.e3 ! pc to kpc
    
    allocate(rans(blobn,3))

    call random_number(rans) ! generate random numbers

    rans(:,1) = rans(:,1)**(1./3.) * blob_rad ! radius
    rans(:,2) = 2.d0*pi*rans(:,2) ! theta (0,2pi)
    rans(:,3) = acos(2.d0*rans(:,3)-1.d0) ! phi (0,pi)

    p_data%r_p(1,offset+1:offset+blobn) = real( rans(:,1)*cos(rans(:,2))*sin(rans(:,3))+blob_cent(1) )
    p_data%r_p(2,offset+1:offset+blobn) = real( rans(:,1)*sin(rans(:,2))*sin(rans(:,3))+blob_cent(2) )
    p_data%r_p(3,offset+1:offset+blobn) = real( rans(:,1)*cos(rans(:,3))               +blob_cent(3) )

    
    deallocate(rans)
    
    ! convert from temperature in K to internal units (1e10 erg/g)
    p_data%u_p(offset+1:offset+blobn-1) = real(T0/u_to_TK/uunit_cgs)

    p_data%v_p(1,offset+1:offset+blobn) = real(blob_v(1))
    p_data%v_p(2,offset+1:offset+blobn) = real(blob_v(2))
    p_data%v_p(3,offset+1:offset+blobn) = real(blob_v(3))

    ! set initial IDs
    do ip=offset,offset+blobn-1
        p_data%id_p(ip) = ip
    end do
    

end subroutine blobfile_locs













