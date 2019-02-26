! create a disc of gas from an IC file

module IC_parameters
    
    real(kind=8) :: mtot, disc_rad, disc_inner_rad, disc_thick ! disc parameters
    real(kind=8) :: m_smbh,v_large,a_scale,c_scale ! potential parameters
    real(kind=8) :: inwards_v
    real(kind=8) :: v0,v0rad,v_index
    real(kind=8) :: T0
    real(kind=8) :: vloop0,vlooprad
    real(kind=8) :: hotcore
    integer :: nhotcore,ndisc
    
    real(kind=8) :: dense_index
    
    real(kind=8) :: hernquist_mass,hernquist_scale
    
    real(kind=8) :: Q_target, rho_target, surf_target
    
    integer :: pot_type
    integer, parameter :: FLAT_POT=0, HERNQUIST_POT=1
    
    logical :: selfgrav
    logical :: constant_surf

    character(len=256) :: position_file


    real(kind=8), parameter :: cooled_sound_speed = .67d0


end module IC_parameters

module Units
    real(kind=8), parameter :: pi = 4.d0*atan(1.d0)
    
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

    real(kind=8), parameter :: nunit = munit_cgs/runit_cgs**3/molecular_mass/proton_mass_cgs
end module Units

program disc_ICs
    use gadget_ics_IO
    use IC_parameters
    implicit none
    
    character(len=256) :: pramfile,intext,invariable,invalue,outfile,trimtext
    character(len=256) :: tproffile
    integer :: ios
    integer :: whitey
    
    integer :: ip,np
    
    logical :: valid_command
    
    if ( command_argument_count()<1 ) then
        print *,"Please give the parameter file in the command line"
        print *,"e.g. ./disc_ICs prams/fatDiscPrams.dat"
        stop
    endif
    
    call get_command_argument(1,pramfile)
    
    
    ! set default values
    
    ! Softened Keplerian in centre, flat rotation curve at large radii
    pot_type = 0
    m_smbh = 1.d6! in solar masses *SOLAR_MASS/ALL.UnitMass_in_g;
    !m_smbh = 1.d7! in solar masses *SOLAR_MASS/ALL.UnitMass_in_g;
    v_large=146. ! in km/s/UnitVelocity_in_cm_per_s;
    a_scale=0.01d0 ! in kpc/CM_PER_MPC/UnitLength_in_cm;
    c_scale=0.0001d0 ! in kpc/CM_PER_MPC/UnitLength_in_cm;
    inwards_v=0.
    
    Q_target = -1 ! i.e. don't fit Q
    rho_target = -1 ! i.e. don't fit density
    surf_target = -1 ! i.e. don't fit surface density
    constant_surf = .false. ! if true, Q_target is *minimum* Q, and try to make surface density & 3d density constant 
    
    position_file = "" ! generate positions randomly
    
    ! hernquist type bulge potential
    hernquist_mass = 1.d9
    hernquist_scale = 250.d0
    
    ! negative = not initial hot core, otherwise = temp of hot core
    hotcore = -1.
    nhotcore = 0

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
    
    vlooprad = -1. ! in pc (i.e. don't do vloop)
    
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
                    if ( pot_type/=FLAT_POT ) then
                        print *,"WARNING - potential parameters set for wrong profile"
                    endif
                case("a_scale")
                    read(invalue,*) a_scale
                    if ( pot_type/=FLAT_POT ) then
                        print *,"WARNING - potential parameters set for wrong profile"
                    endif
                case("c_scale")
                    read(invalue,*) c_scale
                case("inv")
                    read(invalue,*) inwards_v
                    
                case("Q")
                    read(invalue,*) Q_target
                case("CONST_SURF")
                    constant_surf = .true.
                case("nH")
                    read(invalue,*) rho_target
                case("surf")
                    read(invalue,*) surf_target
                
                case("flat")
                    pot_type = FLAT_POT
                case("hernquist")
                    pot_type = HERNQUIST_POT
                
                case("hernquist_scale")
                    read(invalue,*) hernquist_scale
                    if ( pot_type/=HERNQUIST_POT ) then
                        print *,"WARNING - potential parameters set for wrong profile"
                    endif
                case("hernquist_mass")
                    read(invalue,*) hernquist_mass
                    if ( pot_type/=HERNQUIST_POT ) then
                        print *,"WARNING - potential parameters set for wrong profile"
                    endif

                case("vlaw")
                    !read(invalue,*) v_flat
                    read(intext,*,iostat=ios) invariable,v0,v0rad,v_index

                case("vloop")
                    !read(invalue,*) v_flat
                    read(intext,*,iostat=ios) invariable,vloop0,vlooprad

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
                case("hotcore")
                    read(intext,*,iostat=ios) invariable,hotcore,nhotcore

                case("pfile")
                    read(intext,"(A6,A)",iostat=ios) invariable,position_file
                    ! read in N too

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
    vlooprad = vlooprad/1.e3
    
    ndisc = np
    if ( hotcore>0. .and. nhotcore>0 ) then
        print *,"Increasing np from ",np," to ",np+nhotcore," to produce hot core"
        print *,"Increasing mtot from ",mtot," to ",mtot/np*(np+nhotcore)," to maintain resolution"
        mtot = mtot/np*(np+nhotcore)
        np = np+nhotcore
    endif
    
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
    
    call disc_locs(tproffile)
    
    if ( hotcore>0. .and. nhotcore>0 ) then
        call hotcore_locs()
    endif
    
    ! set initial IDs
    do ip=1,p_data%ng
        p_data%id_p(ip) = ip
    end do

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

function epicyclic_parameter(r) result(kappa)
    use IC_parameters
    use Units
    implicit none
    
    real(kind=8) :: r
    
    real(kind=8) :: kappa

    if ( pot_type/=HERNQUIST_POT ) then
        print*,"CAN ONLY USE CONSTANT Q SURFACE DENSITY PROFILE WITH HERNQUIST POTENTIAL"
        stop
    endif
    
    kappa = 0.

    ! smbh component
    kappa = m_smbh/(r**2+c_scale**2) * (2*c_scale**2/(r**2+c_scale**2)+1)

    ! hernquist component
    kappa = kappa + hernquist_mass*(hernquist_scale+r)**2*(2*hernquist_scale/(hernquist_scale+r)+1)
    
    kappa = kappa * G/r
    
    kappa = sqrt(kappa)
    
    return
end function
    

function constant_q_surf(r,cs,Q) result(surf)
    use IC_parameters
    use Units
    implicit none

    real(kind=8) :: r,cs,Q
    real(kind=8) :: surf
    real(kind=8) :: epicyclic_parameter
    
    surf = cs*epicyclic_parameter(r)/pi/G/Q
    
    return
end function

! assume surface density decreases monotonically with distance at constant Q
subroutine calc_max_q_surf_profile
    use IC_parameters
    use Units
    implicit none

    real(kind=8) :: constant_q_surf

    
    real(kind=8) :: r_guess, r_lowerbound,r_upperbound
    real(kind=8) :: mass_guess
    real(kind=8) :: sensitivity, error
    
    sensitivity = 1.d-4
    
    error = sensitivity+1.d0
    
    r_lowerbound = disc_inner_rad
    r_upperbound = disc_inner_rad*1.d3
    
    do while ( error>=sensitivity )
        r_guess = (r_lowerbound+r_upperbound)*.5d0
        mass_guess = constant_q_surf(r_guess,cooled_sound_speed,Q_target)
        mass_guess = mass_guess * pi * (r_guess**2 - disc_inner_rad**2)
        error = abs((mass_guess-mtot)/mtot)

        !print *,r_lowerbound,r_guess,r_upperbound,mass_guess,mtot,error

        if ( mass_guess>mtot ) then
            r_upperbound = r_guess
        else
            r_lowerbound = r_guess
        endif
    end do
    
    dense_index = 0.d0
    disc_rad = r_guess
    
    disc_thick = constant_q_surf(r_guess,cooled_sound_speed,Q_target)/rho_target/2.d0
end subroutine

subroutine calc_const_q_surf_profile(rans)
    use IC_parameters
    use gadget_ics_IO
    use Units
    implicit none
    
    integer :: ibin
    integer, parameter :: nbin=100
    
    real(kind=8) :: constant_q_surf


    real(kind=8) :: r_inner,r_outer,rbin
    real(kind=8) :: m_so_far,mbin,bin_surf
    
    
    real(kind=8) :: stepsize
    
    real(kind=8), dimension(:), allocatable :: fbin
    real(kind=8), dimension(p_data%ng,3) :: rans
    
    integer :: ip
    
    stepsize = (disc_rad-disc_inner_rad)/nbin
    r_inner = disc_inner_rad
    m_so_far = 0.d0
    
    
    ! recalculate disc outer radius, given surface density profile from Q=constant criterion
    do while (m_so_far<mtot)
        r_outer = r_inner+stepsize
        rbin = (r_inner+r_outer)/2.d0
        bin_surf = constant_q_surf(rbin,cooled_sound_speed,Q_target)
        mbin = bin_surf*pi*(r_outer**2-r_inner**2)
        m_so_far = m_so_far + mbin
!         print *,rbin*1d3,bin_surf*(1d4),mbin*1d10,m_so_far*1d10
        r_inner = r_outer
    end do

    if ( m_so_far==mtot ) then
        disc_rad = r_outer
    else
        ! interpolate to guess final position
        disc_rad = sqrt((mtot-m_so_far+mbin)/pi/bin_surf+(r_outer-stepsize)**2)
    endif
    

    ! redo integration, tracking cumulative mass to make probability frequency
    allocate(fbin(nbin))
    
    stepsize = (disc_rad-disc_inner_rad)/nbin
    r_inner = disc_inner_rad
    m_so_far = 0.d0
    do ibin=1,nbin
        r_outer = r_inner+stepsize
        rbin = (r_inner+r_outer)/2.d0
        bin_surf = constant_q_surf(rbin,cooled_sound_speed,Q_target)
        mbin = bin_surf*pi*(r_outer**2-r_inner**2)
        m_so_far = m_so_far + mbin
        fbin(ibin) = m_so_far
        r_inner = r_outer
    end do
    
    fbin = fbin/fbin(nbin)
    
    do ip=1,ndisc
        ibin = 1
        do while (fbin(ibin)<rans(ip,1) .and. ibin<=nbin)
            ibin = ibin + 1
        end do
        r_inner = disc_inner_rad+(ibin-1)*stepsize
        if ( ibin==1 ) then
            rans(ip,1) = r_inner+rans(ip,1)/fbin(1)*stepsize
        else
            rans(ip,1) = r_inner+(rans(ip,1)-fbin(ibin-1))/(fbin(ibin)-fbin(ibin-1))*stepsize
        endif
        if ( rho_target>=1 ) then
            ! set thickness to fit rho
            bin_surf = constant_q_surf(rans(ip,1),cooled_sound_speed,Q_target)
            rans(ip,3) = (rans(ip,3)-.5d0)*(bin_surf/rho_target)
        endif
    end do
    
end subroutine

subroutine disc_locs(tproffile)
    use gadget_ics_IO
    use IC_parameters
    use Units
    implicit none
    
    integer :: ip
    character(len=256) :: tproffile
    
    real(kind=8), dimension(:,:), allocatable :: rans
    real(kind=8), dimension(:),allocatable :: rads, vcircs, mrads, vrads, vphis
    real(kind=8), dimension(:),allocatable :: rad2d, z_over_r, rad_transformed
    real(kind=8), dimension(:,:),allocatable :: r_transformed
    
    ! unit conversions
    m_smbh = m_smbh/1.d10
    hernquist_mass = hernquist_mass/1.d10
    v_large = v_large*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    inwards_v = inwards_v*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    vloop0 = vloop0*1.d5*(tunit_cgs/runit_cgs) ! tunit/runit is basically 1 km/s anyway
    ! a_scale and c_scale are already in kpc

    if ( rho_target>0. ) then
        rho_target = rho_target / nunit
    endif
    
    allocate(rans(ndisc,3))
    call random_number(rans) ! generate random numbers
    if ( position_file=="" ) then

    
        ! radius - for constant density
        !rans(:,1) = dsqrt(rans(:,1))*disc_rad ! no inner radius
        !rans(:,1) = dsqrt(rans(:,1)*(disc_rad**2-disc_inner_rad**2)+disc_inner_rad**2)  
    
        if ( Q_target>0. .and. .not. constant_surf) then
            call calc_const_q_surf_profile(rans)
            !call calc_max_q_surf_profile
        else
            if ( constant_surf ) then
                if ( Q_target<=0. ) then
                    print *,"Need to set Q_target to have constant_surf."
                    print *," To set surface manually, turn off constant_surf and use 'surf [x]' instead"
                    stop
                endif
                
                call calc_max_q_surf_profile
            endif
            
            if ( surf_target>0. ) then
                disc_rad=sqrt(disc_inner_rad**2+mtot/pi/surf_target)
                print *,"Outer radius set to ",disc_rad
            endif
            
            ! radius - for power-law density
            rans(:,1) = (  rans(:,1)*(disc_rad**(dense_index+2.d0)-disc_inner_rad**(dense_index+2.d0)) &
                        +disc_inner_rad**(dense_index+2.d0)   )**(1./(dense_index+2.d0))  
            if ( dense_index==0.d0 ) then
                print *,"Surface density:",mtot/pi/(disc_rad**2-disc_inner_rad**2)
            endif
        endif


        !rans(:,1) = (rans(:,1)*(disc_rad**3-disc_inner_rad**3)+disc_inner_rad**3)**(1.d0/3.d0) ! with flaring disc
        ! theta coordinate
        rans(:,2) = rans(:,2)*2.*pi
        ! height above/below plane
        !rans(:,3) = (rans(:,3)-.5d0)*2.d0*(disc_thick*rans(:,1)) ! thickness is disk_thick*rad
        if ( rho_target<=0 .or. constant_surf ) then
            rans(:,3) = (rans(:,3)-.5d0)*2.d0*disc_thick ! thickness is disk_thick
        endif
    
        p_data%r_p(1,1:ndisc) = real(rans(:,1)*dsin(rans(:,2)))
        p_data%r_p(2,1:ndisc) = real(rans(:,1)*dcos(rans(:,2)))
        p_data%r_p(3,1:ndisc) = real(rans(:,3))
    
    else
        open(unit=16,file=position_file)
        do ip=1,ndisc
            read(16,*) p_data%r_p(2,ip), p_data%r_p(1,ip)
            ! centre arbitrarily
            p_data%r_p(1,ip) = (p_data%r_p(1,ip)-325)*1.e-5
            p_data%r_p(2,ip) = (-p_data%r_p(2,ip)+180)*1.e-5
            p_data%r_p(3,ip) = real((rans(ip,3)-.5e0)*2.e-5+.5e-3)
        end do
        
        close(16)
    endif
    deallocate(rans)
    
    ! Velocities - use spherical approximation, which I think is what Widrow Pym Dubinski 2008 use
    ! so a_r=GM(r)/r^2 and v_circ=sqrt(r*a_r)
    ! vertical eqm is just from temperature, and we don't need velocity dispersion because we have temperature
    !
    ! density is constant, so M(r)=G * M_tot * r^2 / r_max^2
    ! add external potential to this
    
    allocate(rads(ndisc))
    allocate(vcircs(ndisc))
    allocate(vrads(ndisc))
    
    rads = sqrt(p_data%r_p(1,1:ndisc)**2+p_data%r_p(2,1:ndisc)**2+p_data%r_p(3,1:ndisc)**2) ! sqrt is slow, but this is O(N) so okay
    !
    !mrads = mtot*rads**2/disc_rad**2 ! disc mass for self gravity

    
    ! Temperatures - use constant temperature if tproffile is not assigned
    ! otherwise, interpolate from tproffile
    if ( trim(tproffile)=='none' ) then
        p_data%u_p(1:ndisc)=real(T0)

    
        ! Calculate rotation curve from potential, if no vlaw is given (i.e. input v0rad<0), and no explicit acceleration profile is given
        ! Otherwise, set vcircs(:)=v0*(R/v0rad)**v_index
        if ( v0rad<0. ) then
            allocate(mrads(ndisc))
    
            if ( selfgrav ) then
                !mrads = mtot*(rads-disc_inner_rad)**2/(disc_rad-disc_inner_rad)**2 ! disc mass for self gravity with inner cutoff
                 ! disc mass for self gravity with inner cutoff and power-law surface density
                mrads = mtot*(rads-disc_inner_rad)**(2.d0+dense_index)/(disc_rad-disc_inner_rad)**(2.d0+dense_index)
            else
                mrads = 0.d0 ! only external potential
            endif
    
            ! smbh potential
            mrads = mrads + m_smbh*rads**2/(rads**2+c_scale**2)
            
            select case(pot_type)
                case(FLAT_POT)
                    mrads = mrads + v_large**2*rads**2/G/dsqrt(rads**2+a_scale**2)
                case(HERNQUIST_POT)
                    mrads = mrads + hernquist_mass*(rads/(hernquist_scale+rads))**2
                case default
                    print*,"POT NOT FOUND"
            end select
            vcircs = dsqrt(G*mrads/rads)
        
            deallocate(mrads)
        else
            vcircs = v0*(rads/v0rad)**v_index
        endif
        print *,"vcirc - min,max",minval(vcircs),maxval(vcircs)
    else
        call temp_profile(tproffile,rads,vcircs)
    endif 
    
    if ( vlooprad<0 ) then
        vrads = - inwards_v
    else
        allocate(rad2d(ndisc))
        allocate(z_over_r(ndisc))
        allocate(vphis(ndisc))
        allocate(r_transformed(3,ndisc))
        allocate(rad_transformed(ndisc))
    

        rad2d = sqrt(p_data%r_p(1,1:ndisc)**2+p_data%r_p(2,1:ndisc)**2) - disc_inner_rad
        z_over_r = rad2d/p_data%r_p(3,1:ndisc)
        
        vrads = (1.d0-z_over_r**2)/(z_over_r**2+1.d0)
        vphis = 2.d0*z_over_r/(z_over_r**2+1.d0)
        
        r_transformed(3,:) = p_data%r_p(3,1:ndisc)
        r_transformed(1:2,:) = p_data%r_p(1:2,1:ndisc)*spread(1.d0-disc_inner_rad/rad2d,1,2)
        
        rad_transformed=sqrt(r_transformed(1,:)**2+r_transformed(2,:)**2+r_transformed(3,:)**2)
        
        vrads = vloop0*vrads*erfc(rad_transformed/vlooprad)/2.
        vphis = vloop0*vphis*erfc(rad_transformed/vlooprad)/2.
    endif

    ! convert from temperature in K to internal units (1e10 erg/g)
    p_data%u_p(1:ndisc) = real(p_data%u_p(1:ndisc)/u_to_TK/uunit_cgs)

    if ( position_file/="" ) then
        !vcircs = 0.d0
        !print *,"velocities=0 for great justice"
    endif

    !vcircs = dsqrt(G*mtot*rads/disc_rad**2)
    !print *,vcircs*runit_cgs/tunit_cgs*1.e-5
!     print *,p_data%v_p(:,1),p_data%r_p(:,1)

    
    ! convert to cartesian
    p_data%v_p(1,1:ndisc) =  real(vcircs * p_data%r_p(2,1:ndisc) / rads)
    p_data%v_p(2,1:ndisc) = real(-vcircs * p_data%r_p(1,1:ndisc) / rads)
    p_data%v_p(3,1:ndisc) = 0.

    print *,p_data%v_p(:,1),p_data%r_p(:,1)
    
    if ( vlooprad>=0 ) then
    ! TODO: include vphi
        p_data%v_p(1:2,1:ndisc) = real(p_data%v_p(1:2,1:ndisc) + &
                p_data%r_p(1:2,1:ndisc)*spread(vphis*z_over_r/sqrt(z_over_r**2+1.)/rad2d,1,2))
        p_data%v_p(3,1:ndisc) = real(p_data%v_p(3,1:ndisc) - vphis/sqrt(z_over_r**2+1.))
        
        ! include radial impulse
        p_data%v_p(:,1:ndisc) = real(p_data%v_p(:,1:ndisc) + r_transformed(:,:) * spread(vrads(:)/rad_transformed(:),1,3))
    else
        ! include radial impulse
        p_data%v_p(:,1:ndisc) = real(p_data%v_p(:,1:ndisc) + p_data%r_p(:,1:ndisc) * spread(vrads(:)/rads,1,3))
    endif
    
    print *,p_data%v_p(:,1),p_data%r_p(:,1)

    if ( vlooprad>0 ) then
        deallocate(z_over_r)
        deallocate(rad2d)
        deallocate(vphis)
        deallocate(r_transformed)
    endif

    deallocate(rads)
    deallocate(vrads)
    deallocate(vcircs)
    

end subroutine disc_locs

subroutine hotcore_locs
    use gadget_ics_IO
    use IC_parameters
    use Units
    implicit none
    
    real(kind=8), dimension(:,:), allocatable :: rans
    real(kind=8), dimension(:), allocatable :: rad,phi,theta
    
    allocate(rans(nhotcore,3))
    allocate(rad(nhotcore))
    allocate(phi(nhotcore))
    allocate(theta(nhotcore))
    call random_number(rans) ! generate random numbers

    phi = rans(:,1)*2.d0*pi
    theta = acos(rans(:,2)*2.d0-1.d0)
    rad = disc_inner_rad*rans(:,3)**(1.d0/3.d0)

    p_data%r_p(1,ndisc+1:ndisc+nhotcore) = real(rad*sin(theta)*cos(phi))
    p_data%r_p(2,ndisc+1:ndisc+nhotcore) = real(rad*sin(theta)*sin(phi))
    p_data%r_p(3,ndisc+1:ndisc+nhotcore) = real(rad*cos(theta))

    deallocate(rans)
    deallocate(rad)
    deallocate(phi)
    deallocate(theta)
    
    p_data%v_p(:,ndisc+1:ndisc+nhotcore)=0.d0

    p_data%u_p(ndisc+nhotcore)=real(T0/u_to_TK/uunit_cgs)

end subroutine hotcore_locs

! set temperature profile by interpolating from points in a file
subroutine temp_profile(tproffile, rads, vcircs)
    use gadget_ics_IO
    use IC_parameters, only: ndisc 
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
    
    do ip=1,ndisc
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
















