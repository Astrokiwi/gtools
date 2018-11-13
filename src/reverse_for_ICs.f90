! 
! function accel(xyz,N)
!     use global_constants
!     implicit none
!     
!     integer :: N
!     real(kind=8),dimension(N,3) :: xyz,accel
!     real(kind=8) :: r,x,m
!     integer :: i
! 
!     do i=1,N
!         r=sqrt(sum(xyz(i,:)**2))
!         x=r/a_hernquist
!         m=m_smbh*r**2/(r**2+c_scale**2) + m_hernquist*(x/(1.+x))**2
!     
!         accel(i,:)=(-G_internal*m/r**3) * xyz(i,:)
!     end do
! 
!     return
! end function

program reverse_for_ICs
!     use global_constants
    implicit none


    real(kind=8), parameter ::  G_cgs = 6.67259e-8,&
                                solarmass_in_g = 1.989e33,&
                                pc_in_cm = 3.086e18,&
                                yr_in_s = 31556926,&
                                km_in_cm = 1.e5,&
                                G_internal = G_cgs/pc_in_cm**2*solarmass_in_g/km_in_cm,&
                                M_PI = 4.0_8*atan(1.0_8)

    real(kind=8), parameter ::  a_hernquist=250.,&
                                m_hernquist=1.d9,&
                                m_smbh=1.d6,&
                                c_scale=1.d-2
    
    integer, parameter :: N=100
    real(kind=8), parameter ::  rmin=0.005,&
                                inflow_rate=0.01,&
                                mean_rot_vel=600.!100.

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

    open(unit=17,file="data/stronger_rot.dat")
    do i=1,N
        write(17,"(6E13.5)") xyz(i,:),vel(i,:)
    end do
    close(17)

end program reverse_for_ICs
