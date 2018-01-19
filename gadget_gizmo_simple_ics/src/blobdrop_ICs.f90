! create a random blob of just gas

module blobdata
    integer, parameter :: n_maxblobs = 100

    real(kind=8),dimension(n_maxblobs) :: mtot, blob_dist, blob_rad, blob_vcirc, blob_vrad
    real(kind=8),dimension(n_maxblobs) :: T0
    integer,dimension(n_maxblobs)      :: blobn
    
    integer :: nblobs

    contains
    
    subroutine parse_blobline(intext)
        implicit none
        character(len=256) :: intext,textleft
        character(len=256) :: word
        character(len=8) :: format_str
        character(len=256) :: command
        
        nblobs = nblobs + 1
        
        read(intext,"(A4,A)") word,textleft ! remove "blob" part
        command = "0"
        
        do while (len(trim(textleft))>0 )
            textleft = adjustl(textleft)
            read(textleft,*) word
            write(format_str,"('(A',I2.2,',A)')") len(trim(word))
            read(textleft,format_str) word,textleft
            if ( command=="0" ) then
                command = word
                !print *,trim(word),"command"
            else
                select case(trim(command))
                    case("mtot")
                        read(word,*) mtot(nblobs)
                    case("temp")
                        read(word,*) T0(nblobs)
    
                    case("dist")
                        read(word,*) blob_dist(nblobs)
                    case("rad")
                        read(word,*) blob_rad(nblobs)
    
                    case("vcirc")
                        read(word,*) blob_vcirc(nblobs)
                    case("vrad")
                        read(word,*) blob_vrad(nblobs)
    
                    case("N")
                        read(word,*) blobn(nblobs)
                    case default
                        print *,trim(command),trim(word)," not understood"
                end select
                !print *,trim(word),"value"
                command = "0"
            endif
        enddo
        
    end subroutine
end module blobdata

program blobdrop_ICs
    use gadget_ics_IO
    use blobdata
    implicit none
    
    character(len=256) :: pramfile,intext,invariable,invalue,outfile,trimtext
    integer :: ios
    integer :: whitey

    integer :: np
    
    logical :: valid_command
    
    integer :: ib, n_so_far
    
    if ( command_argument_count()<1 ) then
        print *,"Please give the parameter file in the command line"
        print *,"e.g. ./disc_ICs prams/fatDiscPrams.dat"
        stop
    endif
    
    call get_command_argument(1,pramfile)
    
    
    ! set default values
    
    mtot = 1.d4 ! units of Msun
    
    np = 10**4
    
    blob_dist = 5.d0 ! in pc
    blob_rad = 1.d0 ! in pc
    
    blob_vcirc = 0. ! in km/s
    blob_vrad = 0. ! in km/s

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
                case("blob")
                    call parse_blobline(intext)
                case("outfile")
                    whitey = scan(intext,' ')
                    outfile = adjustl(trim(intext(whitey:)))
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
    blob_dist = blob_dist/1.e3
    blob_rad = blob_rad/1.e3
    
    np = sum(blobn(1:nblobs))
    
    call set_N_onlygas(np)
    
    ! only mp_g(1) is actually read, but we need to give a value for 
    ! every entry, otherwise gadget/gizmo will try to read an array of 
    ! mass values per particle from the IC file
    header%mp_g(1:6) = sum(mtot(1:nblobs))/header%np(1)
    
    call set_cosmo_noncosmo
    
    call set_flags

    header%boxSize = 0. ! ignored with periodic boundaries turned off

    call p_data%doAllocations
    
    n_so_far = 0
    do ib=1,nblobs
        call blobdrop_locs(n_so_far,blobn(ib),blob_dist(ib), blob_rad(ib), blob_vcirc(ib), blob_vrad(ib), T0(ib))
        n_so_far = n_so_far + blobn(ib)
    end do
    
    call write_ICs(outfile)
    
end program blobdrop_ICs

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

subroutine blobdrop_locs(offset, blobn, blob_dist, blob_rad, blob_vcirc, blob_vrad, T0)
    use gadget_ics_IO
    implicit none
    
    integer :: offset, blobn
    
    integer :: ip
    real(kind=8) :: blob_dist, blob_rad, blob_vcirc, blob_vrad, T0
    
    real(kind=8), dimension(:,:), allocatable :: rans
    real(kind=8), dimension(:),allocatable :: rads
    
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
    
    real(kind=8), dimension(3) :: blob_cent
    
    ! unit conversions
     ! n.b. tunit/runit is basically 1 km/s anyway
    blob_vrad = blob_vrad*1.e5*(tunit_cgs/runit_cgs)
    blob_vcirc = blob_vcirc*1.e5*(tunit_cgs/runit_cgs)
    
    call random_number(blob_cent)
    blob_cent(3) = 2.d0*pi*blob_cent(3)
    blob_cent(1) = sin(blob_cent(3)) * blob_dist
    blob_cent(2) = cos(blob_cent(3)) * blob_dist
    blob_cent(3) = 0.d0
    
    allocate(rans(blobn,3))

    call random_number(rans) ! generate random numbers

    rans(:,1) = rans(:,1)**(1./3.) * blob_rad ! radius
    rans(:,2) = 2.d0*pi*rans(:,2) ! theta (0,2pi)
    rans(:,3) = acos(2.d0*rans(:,3)-1.d0) ! phi (0,pi)

    p_data%r_p(1,offset+1:offset+blobn) = real( rans(:,1)*cos(rans(:,2))*sin(rans(:,3))+blob_cent(1) )
    p_data%r_p(2,offset+1:offset+blobn) = real( rans(:,1)*sin(rans(:,2))*sin(rans(:,3))+blob_cent(2) )
    p_data%r_p(3,offset+1:offset+blobn) = real( rans(:,1)*cos(rans(:,3))               +blob_cent(3) )

    
    deallocate(rans)
    
    ! Velocities - use spherical approximation, which I think is what Widrow Pym Dubinski 2008 use
    ! so a_r=GM(r)/r^2 and v_circ=sqrt(r*a_r)
    ! vertical eqm is just from temperature, and we don't need velocity dispersion because we have temperature
    !
    ! density is constant, so M(r)=G * M_tot * r^2 / r_max^2
    ! add external potential to this
    
    allocate(rads(blobn))
    print *,shape(rads),blobn
    
    rads = sqrt(p_data%r_p(1,offset+1:offset+blobn)**2+&
                p_data%r_p(2,offset+1:offset+blobn)**2+&
                p_data%r_p(3,offset+1:offset+blobn)**2) ! sqrt is slow, but this is O(N) so okay

    ! convert from temperature in K to internal units (1e10 erg/g)
    p_data%u_p(offset+1:offset+blobn-1) = real(T0/u_to_TK/uunit_cgs)

    ! convert to cartesian
    p_data%v_p(1,offset+1:offset+blobn) =  real(blob_vcirc * p_data%r_p(2,offset+1:offset+blobn) / rads)
    p_data%v_p(2,offset+1:offset+blobn) = real(-blob_vcirc * p_data%r_p(1,offset+1:offset+blobn) / rads)
    p_data%v_p(3,offset+1:offset+blobn) = 0.

    ! include radial impulse
    p_data%v_p(:,offset+1:offset+blobn) =  real(p_data%v_p(:,offset+1:offset+blobn) &
                + p_data%r_p(:,offset+1:offset+blobn)  * spread(blob_vrad/rads,1,3))
    deallocate(rads)

    ! set initial IDs
    do ip=offset,offset+blobn-1
        p_data%id_p(ip) = ip
    end do
    

end subroutine blobdrop_locs













