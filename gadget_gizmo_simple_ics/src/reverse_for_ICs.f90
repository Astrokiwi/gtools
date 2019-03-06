! This code reads in some parameters and uses them to integrate some blobs outwards
! The blobs are "frozen" after some time, in turn
! Use blobfile_ICs to convert this into GIZMO/GADGET ICs
! With the given potential, the blobs will then infall with the given inflow rate (in blobs/yr)

module blob_reverse_prams
    real(kind=8), parameter ::  G_cgs = 6.67259e-8,&
                                solarmass_in_g = 1.989e33,&
                                pc_in_cm = 3.086e18,&
                                yr_in_s = 31556926,&
                                km_in_cm = 1.e5,&
                                G_internal = G_cgs/pc_in_cm**2*solarmass_in_g/km_in_cm,&
                                M_PI = 4.0_8*atan(1.0_8)

    real(kind=8) ::  a_hernquist=250.,&
                                m_hernquist=1.d9,&
                                m_smbh=1.d6,&
                                c_scale=1.d-2

    integer :: N=100
    real(kind=8) ::  rmin=0.005,&
                                inflow_rate=0.01,&
                                mean_rot_vel=600.!100.
    character(len=256) :: outfile
    
    contains
        subroutine get_prams(pramfile)
            use ic_formatter
            implicit none
            character(len=*) :: pramfile
            
            call parse_prams(pramfile)
            
            N = get_int('N',100)
            a_hernquist = get_real('hernquist_scale',250.d0)
            m_hernquist = get_real('hernquist_mass',1.d9)
            m_smbh = get_real('m_smbh',1.d6)
            c_scale = get_real('c_scale',1.d-2)
            rmin = get_real('rmin',0.005d0)
            inflow_rate = get_real('inflow_rate',0.01d0)
            mean_rot_vel = get_real('vrot',600.d0)
            outfile = get_string('outfile',"blobdata/generic_out.dat")
            
            print *,N,a_hernquist,m_hernquist,m_smbh,N,rmin,inflow_rate,mean_rot_vel,trim(outfile)

            call tidy_up
        end subroutine get_prams

    

end module

program reverse_for_ICs
    use blob_reverse_prams
    implicit none
    character(len=256) :: pramfile
    
    ! Read input from parameter file
    if ( command_argument_count()<1 ) then
        print *,"Please give the parameter file in the command line"
        print *,"e.g. ./reverse_for_ICs blobprams/fastRotPrams.dat"
        stop
    endif
    
    call get_command_argument(1,pramfile)
    
    call get_prams(pramfile)

    call integrate_and_save
end program reverse_for_ICs
    

subroutine integrate_and_save
    use blob_reverse_prams
    implicit none
    real(kind=8) :: v_esc
    real(kind=8) :: r2d,r3d,x,m,vcirc
    real(kind=8) :: vcirc_avg,dvcirc
    real(kind=8), dimension(N,3) :: xyz,vel
    real(kind=8), dimension(3) :: accel
    integer :: i,tick,frame,frameperiod
    integer :: n_active
    
    real(kind=8) :: phi,theta
    
    real(kind=8) :: t
    real(kind=8),parameter :: dt=1.d-2
    character(len=128) :: fname

    ! random initial positions, moving outwards at escape velocity, then add spin
    do i=1,N
        call random_number(phi)
        call random_number(theta)
        phi=2*M_PI*phi
        theta=acos(2.*theta-1.)
        xyz(i,1) = rmin*sin(theta)*sin(phi)
        xyz(i,2) = rmin*sin(theta)*cos(phi)
        xyz(i,3) = rmin*cos(theta)
    end do

    v_esc=sqrt(G_internal*m_smbh/rmin*(pc_in_cm/km_in_cm))
    do i=1,N
        r3d=sqrt(sum(xyz(i,:)**2))
        vel(i,:)=xyz(i,:)*v_esc/r3d*1.1
    end do

    ! distribute angular momentum randomly
    vcirc_avg=0.d0
    do i=1,N
        call random_number(vcirc)
        vcirc = (vcirc*3.-1.)*mean_rot_vel
        r2d=sqrt(xyz(i,1)**2+xyz(i,2)**2)
        vel(i,1)=vel(i,1)-vcirc*xyz(i,2)/r2d
        vel(i,2)=vel(i,2)+vcirc*xyz(i,1)/r2d
        vcirc_avg=vcirc_avg+ (-vel(i,1)*xyz(i,2) + vel(i,2)*xyz(i,1))/r2d
    end do
    vcirc_avg = vcirc_avg/N
    ! correct to make sure mean specific angular momentum is constant between runs
    dvcirc=mean_rot_vel-vcirc_avg

    do i=1,N
        r2d=sqrt(xyz(i,1)**2+xyz(i,2)**2)
        vel(i,1)=vel(i,1)-dvcirc*xyz(i,2)/r2d
        vel(i,2)=vel(i,2)+dvcirc*xyz(i,1)/r2d
    end do
    
    tick=0
    frame=0
    t=0.d0
    frameperiod=1./(10.*inflow_rate*dt)

    n_active=N
    do while (n_active>0)
        t=t+dt
        do i=1,n_active
            r3d=sqrt(sum(xyz(i,:)**2))
            x=r3d/a_hernquist
            m=m_smbh*r3d**2/(r3d**2+c_scale**2) + m_hernquist*(x/(1.+x))**2
    
            accel=(-G_internal*m/r3d**3) * xyz(i,:)
            
            xyz(i,:) = xyz(i,:) + vel(i,:)*dt/977813.106
            vel(i,:) = vel(i,:) + accel*dt/3.16887646e-8
        end do
        if ( mod(tick,frameperiod)==0 ) then
            print *,tick,t,n_active
            write(fname,"('blob_reverse_frames/frame',I6.6,'.dat')") frame
            open(unit=17,file=fname)
            do i=1,N
                write(17,"(6E13.5)") xyz(i,:),vel(i,:)
            end do
            close(17)
            frame=frame+1
        endif
        tick=tick+1
        n_active=N-inflow_rate*t
    end do
    print *,"nframes=",frame

    open(unit=17,file=outfile)
    do i=1,N
        write(17,"(6E13.5)") xyz(i,:),vel(i,:)
    end do
    close(17)

end subroutine integrate_and_save
