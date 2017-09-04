! create a random blob of just gas

program disc_ICs
    use gadget_ics_IO
    implicit none
    
    character(len=256) :: pramfile,intext,invariable,invalue,outfile,trimtext
    character(len=256) :: tproffile
    integer :: ios
    integer :: whitey
    
    real(kind=8) :: mtot, disc_rad, disc_inner_rad, disc_thick ! disc parameters
    real(kind=8) :: m_smbh,v_large,a_scale,c_scale ! potential parameters
    real(kind=8) :: inwards_v
    real(kind=8) :: v0,v0rad,v_index
    real(kind=8) :: T0
    
    real(kind=8) :: dense_index
    
    integer :: np
    
    logical :: selfgrav
    
    logical :: valid_command
    
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
    inwards_v=0.

    !mtot = 2.d6/1.d10 ! units of 1.d10 msun
    mtot = 8.426453d6 ! units of Msun
    !mtot = 1.d2/1.d10 ! units of 1.d10 msun
    
    np = 276118

    disc_rad = 20.d0 ! in pc
    disc_inner_rad = 0.5d0 ! pc
    !disc_thick = 0.001d0 ! in kpc
    disc_thick = 0.1d0 ! in pc
    !disc_thick = 0.01d0 ! in kpc/kpc - this is a slope
    
    v0rad = -1. ! in km/s
    
    selfgrav = .true.
    
    outfile = "data/testout.dat"
    
    tproffile = 'none'
    T0 = -1 ! i.e. unset - if still unset, then set to default below
    
    dense_index = 0.d0 ! flat surface density

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
                case("m_smbh")
                    read(invalue,*) m_smbh
                case("v_large")
                    read(invalue,*) v_large
                case("a_scale")
                    read(invalue,*) a_scale
                case("c_scale")
                    read(invalue,*) c_scale
                case("inv")
                    read(invalue,*) inwards_v

                case("vlaw")
                    !read(invalue,*) v_flat
                    read(intext,*,iostat=ios) invariable,v0,v0rad,v_index

                case("mtot")
                    read(invalue,*) mtot

                case("disc_rad")
                    read(invalue,*) disc_rad
                case("disc_inner_rad")
                    read(invalue,*) disc_inner_rad
                case("disc_thick")
                    read(invalue,*) disc_thick
                case("index")
                    read(invalue,*) dense_index

                case("selfgrav")
                    read(invalue,*) selfgrav

                case("tproffile")
                    read(intext,"(A10,A)",iostat=ios) invariable,tproffile
                    if ( T0>=0 ) then
                        print *,"Warning - tproffile overrides flat temperature"
                    endif
                case("temp")
                    read(invalue,*) T0

                case("N")
                    read(invalue,*) np

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
    
    if ( T0==-1 ) then
        T0 = 30. ! default, arbitrary
    endif
    
    ! units of 1e10 Msun, kpc
    mtot = mtot/1.e10
    disc_rad = disc_rad/1.e3
    disc_inner_rad = disc_inner_rad/1.e3
    disc_thick = disc_thick/1.e3
    v0rad = v0rad/1.e3
    
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
                   m_smbh,v_large,a_scale,c_scale,selfgrav,&
                   inwards_v,v0,v0rad,v_index,tproffile,T0,&
                   dense_index)
    
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
                   m_smbh,v_large,a_scale,c_scale,selfgrav,inwards_v,&
                   v0,v0rad,v_index,tproffile,T0,dense_index)
    use gadget_ics_IO
    implicit none
    
    integer :: ip
    real(kind=8) :: disc_rad, disc_inner_rad, disc_thick, mtot
    real(kind=8) :: m_smbh,v_large,a_scale,c_scale ! potential parameters
    real(kind=8) :: inwards_v
    real(kind=8) :: v0,v0rad,v_index
    real(kind=8) :: T0
    real(kind=8) :: dense_index
    character(len=256) :: tproffile
    
    real(kind=8), dimension(:,:), allocatable :: rans
    real(kind=8), dimension(:),allocatable :: rads, vcircs, mrads
    logical :: selfgrav
    
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
    m_smbh = m_smbh/1.d10
    v_large = v_large*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    inwards_v = inwards_v*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    ! a_scale and c_scale are already in kpc
    
    allocate(rans(p_data%ng,3))

    call random_number(rans) ! generate random numbers
    
    ! radius - for constant density
    !rans(:,1) = dsqrt(rans(:,1))*disc_rad ! no inner radius
    !rans(:,1) = dsqrt(rans(:,1)*(disc_rad**2-disc_inner_rad**2)+disc_inner_rad**2)  
    
    ! radius - for power-law density
    rans(:,1) = (  rans(:,1)*(disc_rad**(dense_index+2.d0)-disc_inner_rad**(dense_index+2.d0)) &
                    +disc_inner_rad**(dense_index+2.d0)   )**(1./(dense_index+2.d0))  

    !rans(:,1) = (rans(:,1)*(disc_rad**3-disc_inner_rad**3)+disc_inner_rad**3)**(1.d0/3.d0) ! with flaring disc
    ! theta coordinate
    rans(:,2) = rans(:,2)*2.*pi
    ! height above/below plane
    !rans(:,3) = (rans(:,3)-.5d0)*2.d0*(disc_thick*rans(:,1)) ! thickness is disk_thick*rad
    rans(:,3) = (rans(:,3)-.5d0)*2.d0*disc_thick ! thickness is disk_thick
    
    p_data%r_p(1,:) = real(rans(:,1)*dsin(rans(:,2)))
    p_data%r_p(2,:) = real(rans(:,1)*dcos(rans(:,2)))
    p_data%r_p(3,:) = real(rans(:,3))
    
    deallocate(rans)
    
    ! Velocities - use spherical approximation, which I think is what Widrow Pym Dubinski 2008 use
    ! so a_r=GM(r)/r^2 and v_circ=sqrt(r*a_r)
    ! vertical eqm is just from temperature, and we don't need velocity dispersion because we have temperature
    !
    ! density is constant, so M(r)=G * M_tot * r^2 / r_max^2
    ! add external potential to this
    
    allocate(rads(p_data%ng))
    allocate(vcircs(p_data%ng))
    
    rads = sqrt(p_data%r_p(1,:)**2+p_data%r_p(2,:)**2+p_data%r_p(3,:)**2) ! sqrt is slow, but this is O(N) so okay
    !
    !mrads = mtot*rads**2/disc_rad**2 ! disc mass for self gravity

    
    ! Temperatures - use constant temperature if tproffile is not assigned
    ! otherwise, interpolate from tproffile
    if ( trim(tproffile)=='none' ) then
        p_data%u_p(:)=real(T0)

    
        ! Calculate rotation curve from potential, if no vlaw is given (i.e. input v0rad<0), and no explicit acceleration profile is given
        ! Otherwise, set vcircs(:)=v0*(R/v0rad)**v_index
        if ( v0rad<0. ) then
            allocate(mrads(p_data%ng))
    
            if ( selfgrav ) then
                !mrads = mtot*(rads-disc_inner_rad)**2/(disc_rad-disc_inner_rad)**2 ! disc mass for self gravity with inner cutoff
                 ! disc mass for self gravity with inner cutoff and power-law surface density
                mrads = mtot*(rads-disc_inner_rad)**(2.d0+dense_index)/(disc_rad-disc_inner_rad)**(2.d0+dense_index)
            else
                mrads = 0.d0 ! only external potential
            endif
    
            mrads = mrads + m_smbh*rads**2/(rads**2+c_scale**2) + v_large**2*rads**2/G/dsqrt(rads**2+a_scale**2)
            vcircs = dsqrt(G*mrads/rads)
        
            deallocate(mrads)
        else
            vcircs = v0*(rads/v0rad)**v_index
        endif
    
    else
        call temp_profile(tproffile,rads,vcircs)
    endif 

    ! convert from temperature in K to internal units (1e10 erg/g)
    p_data%u_p(:) = real(p_data%u_p(:)/u_to_TK/uunit_cgs)

    !vcircs = dsqrt(G*mtot*rads/disc_rad**2)
    !print *,vcircs*runit_cgs/tunit_cgs*1.e-5
    
    ! convert to cartesian
    p_data%v_p(1,:) =  real(vcircs * p_data%r_p(2,:) / rads)
    p_data%v_p(2,:) = real(-vcircs * p_data%r_p(1,:) / rads)
    p_data%v_p(3,:) = 0.
    
    ! include radial impulse
    p_data%v_p(:,:) = real(p_data%v_p(:,:) - p_data%r_p(:,:) * inwards_v / spread(rads(:),1,3))
    
    deallocate(rads)
    deallocate(vcircs)
    
    ! set initial IDs
    do ip=1,p_data%ng
        p_data%id_p(ip) = ip
    end do
    

end subroutine disc_locs


! set temperature profile by interpolating from points in a file
subroutine temp_profile(tproffile, rads, vcircs)
    use gadget_ics_IO
    implicit none
    
    character(len=256) :: tproffile
    real(kind=8), dimension(p_data%ng), intent(in) :: rads
    real(kind=8), dimension(p_data%ng), intent(out) :: vcircs
    
    real(kind=8), dimension(p_data%ng) :: accels

    integer :: ip, ibin, jbin
    integer :: nbins
    integer :: ios
    real(kind=8) :: fbin
    
    real(kind=8), dimension(:), allocatable :: tab_rad, tab_temp, tab_accel
    
    real(kind=8) :: dummy
    logical :: finished
    
    ios = 0
    
    
    nbins = 0
    open(15,file=tproffile)
    do while ( ios==0 )
        read(15,*,iostat=ios)
        if ( ios==0 ) then
            nbins = nbins + 1
        endif
    end do
    close(15)
    
    allocate(tab_rad(nbins))
    allocate(tab_temp(nbins))
    allocate(tab_accel(nbins))
    
    open(15,file=tproffile)
    do ibin=1,nbins
        read(15,*) tab_rad(ibin),dummy,tab_temp(ibin),tab_accel(ibin)
    end do
    close(15)
    
    tab_rad = tab_rad/1.d3 ! convert to kpc
    
    do ip=1,p_data%ng
        if ( rads(ip)>tab_rad(nbins) ) then
            ! assume flat T profile at large radius
            p_data%u_p(ip) = real(tab_temp(nbins))
            accels(ip) = tab_accel(nbins)*tab_rad(nbins)/rads(ip)
            !print *,rads(ip),"far",p_data%u_p(ip)
        elseif ( rads(ip)<tab_rad(1) ) then ! assume power law basically
            fbin = (rads(ip)-tab_rad(1))/(tab_rad(2)-tab_rad(1))
            p_data%u_p(ip) = real(tab_temp(1)**(1-fbin)*tab_temp(2)**(fbin))
            accels(ip) = tab_accel(1)*(1-fbin)+tab_accel(2)*(fbin)
        else
            jbin = -1
            ibin = 1
            finished = .false.
            do while ( .not.finished )
                if ( tab_rad(ibin)>rads(ip) ) then
                    jbin = ibin
                    finished = .true.
                else
                    ibin = ibin + 1
                endif
                if ( ibin==nbins+1 ) then
                    finished = .true.
                endif
            end do
            !print *,rads(ip),jbin,tab_rad(jbin)
            if ( jbin==-1 ) then
                stop 'BIN ERROR'
            endif
            if ( jbin==1 ) then
                stop 'OUT OF BOUNDS ERROR'
            endif
            fbin = (rads(ip)-tab_rad(jbin-1))/(tab_rad(jbin)-tab_rad(jbin-1))
            if ( fbin<0. .or. fbin>1. ) then
                stop 'fbin BOUND ERROR'
            endif
            p_data%u_p(ip) = real(tab_temp(jbin-1)**(1-fbin)*tab_temp(jbin)**(fbin))
            accels(ip) = tab_accel(jbin-1)*(1-fbin)+tab_accel(jbin)*(fbin)
        endif
    end do
    
    vcircs = sqrt(-accels*rads)*555488.7d0 ! sqrt(kpc*cm/s/s) to km/s
    
end subroutine temp_profile
















