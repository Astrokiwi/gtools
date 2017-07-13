! create a random blob of just gas

program disc_ICs
    use gadget_ics_IO
    implicit none
    
    character(len=256) :: pramfile,intext,invariable,invalue,outfile
    integer :: ios
    integer :: whitey
    
    real(kind=8) :: mtot, disc_rad, disc_inner_rad, disc_thick ! disc parameters
    real(kind=8) :: m_smbh,v_large,a_scale,c_scale ! potential parameters
    
    integer :: np
    
    logical :: selfgrav
    
    if ( command_argument_count()<1 ) then
        print *,"Please give the parameter file in the command line"
        print *,"e.g. ./disc_ICs prams/fatDiscPrams.dat"
        stop
    endif
    
    call get_command_argument(1,pramfile)
    
    
    ! set default values
    
    ! Softened Keplerian in centre, flat rotation curve at large radii
    m_smbh = 1.d6! in solar masses *SOLAR_MASS/ALL.UnitMass_in_g;
    !m_smbh = 1.d7! in solar masses *SOLAR_MASS/ALL.UnitMass_in_g;
    v_large=146. ! in km/s/UnitVelocity_in_cm_per_s;
    a_scale=0.01d0 ! in kpc/CM_PER_MPC/UnitLength_in_cm;
    c_scale=0.0001d0 ! in kpc/CM_PER_MPC/UnitLength_in_cm;

    !mtot = 2.d6/1.d10 ! units of 1.d10 msun
    mtot = 8.426453d6 ! units of Msun
    !mtot = 1.d2/1.d10 ! units of 1.d10 msun
    
    np = 276118

    disc_rad = 20.d0 ! in pc
    disc_inner_rad = 0.5d0 ! pc
    !disc_thick = 0.001d0 ! in kpc
    disc_thick = 0.1d0 ! in pc
    !disc_thick = 0.01d0 ! in kpc/kpc - this is a slope
    
    selfgrav = .true.
    
    outfile = "data/testout.dat"

    ! load in parameters to overload
    open(15,file=pramfile)
    ios = 0
    do while ( ios==0 )
        read(15,"(A)",iostat=ios) intext
        if ( ios==0 ) then
            read(intext,*,iostat=ios) invariable,invalue
            if ( ios/=0 ) then
                print *,"INCORRECT FORMAT:"
                print *,intext
                print *,"QUITTING"
                stop
            endif
            select case(invariable)
                case("m_smbh")
                    read(invalue,*) m_smbh
                case("v_large")
                    read(invalue,*) v_large
                case("a_scale")
                    read(invalue,*) a_scale
                case("c_scale")
                    read(invalue,*) c_scale

                case("mtot")
                    read(invalue,*) mtot

                case("disc_rad")
                    read(invalue,*) disc_rad
                case("disc_inner_rad")
                    read(invalue,*) disc_inner_rad
                case("disc_thick")
                    read(invalue,*) disc_thick

                case("selfgrav")
                    read(invalue,*) selfgrav

                case("N")
                    read(invalue,*) np

                case("outfile")
                    whitey = scan(intext,' ')
                    outfile = adjustl(trim(intext(whitey:)))
                case default
                    print *,invariable," not understood"
            end select
        endif
    end do
    close(15)
    
    ! units of 1e10 Msun, kpc
    mtot = mtot/1.e10
    disc_rad = disc_rad/1.e3
    disc_inner_rad = disc_inner_rad/1.e3
    disc_thick = disc_thick/1.e3
    
    
    !call set_N_onlygas(8192)
    !call set_N_onlygas(65536)
    call set_N_onlygas(np)
    
    ! only mp_g(1) is actually read, but we need to give a value for 
    ! every entry, otherwise gadget/gizmo will try to read an array of 
    ! mass values per particle from the IC file
    header%mp_g(1:6) = mtot/header%np(1)
    
    call set_cosmo_noncosmo
    
    call set_flags

    header%boxSize = 0. ! ignored with periodic boundaries turned off

    call p_data%doAllocations
    
    call disc_locs(disc_rad, disc_inner_rad, disc_thick, mtot,&
                   m_smbh,v_large,a_scale,c_scale,selfgrav)
    
    !call write_ICs("data/discICs_midres.dat")
    !call write_ICs("data/holy_razorthin_discICs_midres.dat")
    call write_ICs(outfile)
    !call write_ICs("data/holy_discICs_noself.dat")
    
end program disc_ICs

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
!     header%omega0 = .27
!     header%omegaLambda = .73
!     header%hubbleParam = .7
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

subroutine disc_locs(disc_rad, disc_inner_rad, disc_thick, mtot,&
                   m_smbh,v_large,a_scale,c_scale,selfgrav)
    use gadget_ics_IO
    implicit none
    
    integer :: ip
    real(kind=8) :: disc_rad, disc_inner_rad, disc_thick, mtot
    real(kind=8) :: m_smbh,v_large,a_scale,c_scale ! potential parameters
    real(kind=8), dimension(:,:), allocatable :: rans
    real(kind=8), dimension(:),allocatable :: rads, vcircs, mrads
    logical :: selfgrav
    
    real(kind=8), parameter :: pi = 4.d0*atan(-1.d0)
    
    ! gadget/gizmo internal units
    real(kind=8), parameter :: munit_cgs = 1.989d43 ! 10^10 solar masses
    real(kind=8), parameter :: tunit_cgs = 3.08568d16 ! just under a Gyr
    real(kind=8), parameter :: runit_cgs = 3.0857d21 ! 1 kpc
    !real(kind=8), parameter :: G = 6.67408d-8 * tunit_cgs**2 * munit_cgs / runit_cgs**3 ! G=43007.1
    real(kind=8), parameter :: G = 43007.1 ! from Gizmo output, consist with above units
    
    ! unit conversions
    m_smbh = m_smbh/1.d10
    v_large = v_large*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    ! a_scale and c_scale are already in kpc
    
    allocate(rans(p_data%ng,3))

    call random_number(rans) ! generate random numbers
    
    ! radius - for constant density
    !rans(:,1) = dsqrt(rans(:,1))*disc_rad ! no inner radius
    rans(:,1) = dsqrt(rans(:,1)*(disc_rad**2-disc_inner_rad**2)+disc_inner_rad**2)  
    !rans(:,1) = (rans(:,1)*(disc_rad**3-disc_inner_rad**3)+disc_inner_rad**3)**(1.d0/3.d0) ! with flaring disc
    ! theta coordinate
    rans(:,2) = rans(:,2)*2.*pi
    ! height above/below plane
    !rans(:,3) = (rans(:,3)-.5d0)*2.d0*(disc_thick*rans(:,1)) ! thickness is disk_thick*rad
    rans(:,3) = (rans(:,3)-.5d0)*2.d0*disc_thick ! thickness is disk_thick
    
    p_data%r_p(1,:) = rans(:,1)*dsin(rans(:,2))
    p_data%r_p(2,:) = rans(:,1)*dcos(rans(:,2))
    p_data%r_p(3,:) = rans(:,3)
    
    deallocate(rans)
    
    ! TODO: inner cutoff radius
    
    ! Velocities - use spherical approximation, which I think is what Widrow Pym Dubinski 2008 use
    ! so a_r=GM(r)/r^2 and v_circ=sqrt(r*a_r)
    ! vertical eqm is just from temperature, and we don't need velocity dispersion because we have temperature
    !
    ! density is constant, so M(r)=G * M_tot * r^2 / r_max^2
    ! add external potential to this
    
    allocate(rads(p_data%ng))
    allocate(vcircs(p_data%ng))
    allocate(mrads(p_data%ng))
    
    rads = sqrt(p_data%r_p(1,:)**2+p_data%r_p(2,:)**2+p_data%r_p(3,:)**2) ! sqrt is slow, but this is O(N) so okay
    !
    !mrads = mtot*rads**2/disc_rad**2 ! disc mass for self gravity
    
    if ( selfgrav ) then
        mrads = mtot*(rads-disc_inner_rad)**2/(disc_rad-disc_inner_rad)**2 ! disc mass for self gravity with inner cutoff
    else
        mrads = 0.d0 ! only external potential
    endif
    
    mrads = mrads + m_smbh*rads**2/(rads**2+c_scale**2) + v_large**2*rads**2/G/dsqrt(rads**2+a_scale**2)
    vcircs = dsqrt(G*mrads/rads)
    
    !vcircs = dsqrt(G*mtot*rads/disc_rad**2)
    !print *,vcircs*runit_cgs/tunit_cgs*1.e-5
    
    ! convert to cartesian
    p_data%v_p(1,:) =  vcircs * p_data%r_p(2,:) / rads
    p_data%v_p(2,:) = -vcircs * p_data%r_p(1,:) / rads
    p_data%v_p(3,:) = 0.d0
    
    deallocate(mrads)
    deallocate(rads)
    deallocate(vcircs)
    
    !
    do ip=1,p_data%ng
        p_data%id_p(ip) = ip
    end do
    

end subroutine disc_locs


