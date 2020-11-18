module sph_plotter
    implicit none
    
    integer, parameter :: nkern = 1000
    real(kind=8),dimension(nkern) :: kern_tab,kern_tab2
    real(kind=8), save :: kern_norm
    
    logical, save :: kernel_initialized = .false.
    logical, save :: parallel = .false.

    
    integer(kind=1), parameter :: DENSE_MODE = B'00000000' ! default
    integer(kind=1), parameter :: WEIGHT_MODE = B'00000001'
    integer, parameter :: DENSE_WEIGHT_POS = 0
    integer(kind=1), parameter :: ZSLICE_MODE = B'00000010'
    integer, parameter :: ZSLICE_POS = 1
    integer(kind=1), parameter :: MIN_MODE = B'00000100'
    integer, parameter :: MIN_POS = 2
    integer(kind=1), parameter :: VORINOI_MODE = B'00001000'
    integer, parameter :: VORINOI_POS = 3
    integer(kind=1), parameter :: MAX_MODE = B'00010000'
    integer, parameter :: MAX_POS = 4
    integer(kind=1), parameter :: SDEV_MODE = B'00100000'
    integer, parameter :: SDEV_POS = 5

    logical, save :: verbose = .false.

    real(kind=8), parameter :: M_PI = 2.d0*acos(0.d0)

    contains
    
    subroutine set_parallel
        implicit none
        
        parallel = .true.

        return
    end subroutine

    subroutine set_serial
        implicit none
        
        parallel = .false.

        return
    end subroutine
    
    function sph_vel_absorb(x,y,z,m,h,u,opac,vx,vy,vz,gx,gy,theta,phi,dv,nv,n) result(prof)
        implicit none
        integer :: nv,n
        real(kind=8), dimension(n) :: m,h ! mass (Msun), smoothing (kpc)
        real(kind=8), dimension(n) :: vx,vy,vz ! particle velocities (km/s)
        real(kind=8), dimension(n) :: x,y,z ! particle positions (kpc)
        real(kind=8), dimension(n) :: u ! internal energy (erg/g)
        real(kind=8), dimension(n) :: opac ! cross section per unit mass 
        real(kind=8) :: gx,gy ! how far this ray is off-centre
        real(kind=8) :: theta, phi ! range of velocities to cover
        real(kind=8) :: dv ! range of velocities to cover

        real(kind=8) :: rx,ry,rz ! rotated and transposed particle positions
        real(kind=8) :: dr ! distance from ray to particle centre
        real(kind=8) :: vlos ! line of site velocity
        real(kind=8) :: vwidth2 ! gaussian width squared
        real(kind=8) :: mdepth,vdepth ! depth for this particular particle
        real(kind=8), dimension(n) :: zlos ! line of sight position
        integer, dimension(n) :: zarg ! line of sight order
        real(kind=8), dimension(3) :: ray ! ray direction vector
        real(kind=8), dimension(3) :: cent ! transpose particles by this much so that ray goes through (0,0)
        real(kind=8), dimension(3) :: b1,b2 ! basic vectors for calculating cent
        
        real(kind=8), dimension(nv) :: prof ! output profile
        
        integer :: i,ip ! this particle
        integer :: iv ! position in velocity profile array
        real(kind=8) :: v_here
        
        real(kind=8) :: cum_depth

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif

        ray = (/cos(theta)*cos(phi),sin(theta)*cos(phi),sin(phi)/)
        b1 = (/cos(theta)*sin(phi),sin(theta)*sin(phi),-cos(phi)/)
        b2 = (/ ray(2)*b1(3)-ray(3)*b1(2),&
                ray(3)*b1(1)-ray(1)*b1(3),&
                ray(1)*b1(2)-ray(2)*b1(1)   /)
        
        cent = b1*gx + b2*gy

        ! for extinction
        zlos = x*ray(1)+y*ray(2)+z*ray(3)
        
!         print *,"sorting"
        
        zarg = rargsort(zlos) ! SLOW! presort maybe? store values?
!         print *,"sorted"

        prof = 0.d0
        cum_depth = 0.d0
        i = 1
        
        do ip=1,n
!         do while (i<=n .and. cum_depth<10.d0)
!             ip = zarg(i)
        
            ! `ray` has magnitude 1
            ! So: distance from ray to particle is just the norm of the cross product
            rx = x(ip)-cent(1)
            ry = y(ip)-cent(2)
            rz = z(ip)-cent(3)

            dr =      (ry*ray(3)-rz*ray(2))**2
            dr = dr + (rz*ray(1)-rx*ray(3))**2
            dr = dr + (rx*ray(2)-ry*ray(1))**2
            !dr = sqrt(dr) ! could use dr**2 vs h**2 and be faster
!             print *,rx,ry,rz
!             print *,ray
!             print *,ip,dr,h(ip)
        
            if ( dr<=h(ip)**2 ) then
            

                dr = sqrt(dr)
                mdepth = m(ip)*fkern(dr/h(ip)) ! UNITS???
                vdepth = mdepth * exp(-cum_depth)
                cum_depth = cum_depth + mdepth * opac(ip)

                !vdepth = mdepth
            
                vlos = vx(ip)*ray(1)+vy(ip)*ray(2)+vz(ip)*ray(3)
                vwidth2 = 10.d0/9.d0*u(ip) * 1.d-10 ! sound speed squared in km**2/s**2
                do iv=1,nv ! SLOW
                    v_here = (nv/2-iv)*(dv*2.d0/nv)
                    !prof(iv) = prof(iv) + vdepth
                    prof(iv) = prof(iv) + vdepth * exp(-(vlos-v_here)**2/2.d0/vwidth2)
                    !print *,iv,v_here,vlos-v_here,prof(iv)
                end do
            endif
!             i = i + 1
!             stop
        end do
        
    end function sph_vel_absorb

    function sph_dense_spherical(x,y,z,m,h,L1,L2,c,w,n) result(g)
        use omp_lib
        implicit none
        integer :: n
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        integer :: L1,L2 ! ncells in each direction, as parameters for implicit array declaration
        integer, dimension(2) :: L ! ncells in each direction, in array for tidiness
        real(kind=8), dimension(2) :: c ! top-left corner in radians
        real(kind=8), dimension(2) :: w ! size in radians in each direction
        
        real(kind=8), dimension(L1,L2) :: g ! mass/density grid, to return
        
        integer :: ip ! this particle
        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer, dimension(2) :: ih ! integer h, 2D because cells may not be square
        
        real(kind=8), dimension(2) :: r_cell, l_cell
        
        real(kind=8) :: rdist, weight

        real(kind=8) :: rad,phi,theta
        
        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = 0.d0
        
        L(1)=L1
        L(2)=L2
        
        l_cell = w/L ! cell size in radians

        
        do ip=1,n
            rad = sqrt(x(ip)**2+y(ip)**2+z(ip)**2)
            theta = atan2(y(ip),x(ip))
            phi = acos(z(ip)/rad) - M_PI/2.
            
            r_cell = rad*l_cell ! cell size in cartesian (2D)

            ix = int((phi-c(1))/l_cell(1))+1
            iy = int((theta-c(2))/l_cell(2))+1

            ih = ceiling(h(ip)/r_cell)
    
            ix0 = max(1,ix-ih(1))
            iy0 = max(1,iy-ih(2))

            ix1 = min(L(1),ix+ih(1))
            iy1 = min(L(2),iy+ih(2))
    
            if ( ix0<=L(1) .and. iy0<=L(2) .and. ix1>=1 .and. iy1>=1 ) then
                do hix=ix0,ix1
                    do hiy=iy0,iy1
                        rdist = sqrt(( (hix-ix)*r_cell(1))**2 + ((hiy-iy)*r_cell(2))**2 )
                        weight = fkern(rdist/h(ip))/h(ip)**2
                    
                        g(hix,hiy) = g(hix,hiy) + weight * m(ip)
                    end do
                end do

            endif

        end do
        return
    end function sph_dense_spherical
    
    
    function sph_general(x,y,m,h,v,L,c,w,f,z,zslice,mode,n) result(g)
        use omp_lib
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8) :: zslice ! slice to measure value along

        integer(kind=1) :: mode
        
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(:,:), allocatable :: mg ! mass or distance grid
        
        real(kind=8), dimension(:,:), allocatable :: g2 ! secondary grid, if needed
        
        integer :: ip ! this particle
        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, weight
        real(kind=8) :: dz
        
        if ( parallel ) then
            if ( verbose ) then
                print *,"parallel mode"!", nprocs=",omp_get_num_threads()
            endif
        else
            if ( verbose ) then
                print *,"serial mode"
            endif
            call    OMP_SET_NUM_THREADS(1)
        endif
        
        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        
        if ( btest(mode,MIN_POS) ) then
            g = huge(g(1,1))
        else if ( btest(mode,MAX_POS) ) then
            g = -huge(g(1,1))
        else
            g = 0.d0
        endif
        
        allocate(mg(L,L))
        if ( btest(mode,VORINOI_POS) ) then
            mg = huge(mg(1,1))
        else
            mg = 0.d0
        endif
        
        allocate(g2(L,L))
        if ( btest(mode,SDEV_POS) ) then
            g2 = 0.d0
        endif
        r_cell = w/L
        area_cell = r_cell**2
        
!         if (any(isnan(v)) .and. btest(mode,DENSE_WEIGHT_POS)) then
!             print *,"v",n,count(isnan(v))
!         endif
!         if (any(isnan(m))) then
!             print *,"m",n,count(isnan(v))
!         endif

! turn off parallel completely - gives segmentation faults??
!$OMP PARALLEL DO private(ip,ix,iy,hix,hiy,ix0,iy0,ix1,iy1,ih)&
!$OMP& private(rdist,weight,dz)&
!$OMP& shared(m,h,v,x,y,z,w,c,zslice,mode,f,L,r_cell,area_cell)&
!$OMP& reduction(+:g,g2,mg)&
!$OMP& default(none) schedule(dynamic,1)
        do ip=1,n
!             print *,"in loop: parallel mode, nprocs=",omp_get_num_threads()
            if ( f(ip) .and. &
             (.not.btest(mode,DENSE_WEIGHT_POS) .or. &
                .not. ( isnan(v(ip)) .or. v(ip)>HUGE(v(ip)) .or. v(ip)<-HUGE(v(ip)) ) )&
                   ) then


                ix = int((x(ip)-c(1))/r_cell)+1
                iy = int((y(ip)-c(2))/r_cell)+1
                if ( btest(mode,MIN_POS) ) then
                    if ( ix<=L .and. iy<=L .and. ix>=1 .and. iy>=1 ) then
                        g(ix,iy) = min(v(ip),g(ix,iy))
                    endif
                else
                    if ( btest(mode,ZSLICE_POS) ) then
                        dz = abs(z(ip)-zslice)
                    else
                        dz = -1
                    endif
                    
                    if ( dz<h(ip) ) then
                
                        ih = ceiling(h(ip)/r_cell)
                
                        ix0 = max(1,ix-ih)
                        iy0 = max(1,iy-ih)

                        ix1 = min(L,ix+ih)
                        iy1 = min(L,iy+ih)
                
                        if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
                            do hix=ix0,ix1
                                do hiy=iy0,iy1
                                    if ( btest(mode,ZSLICE_POS) ) then
                                        rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell+dz**2)
                                        weight = kern(rdist/h(ip))/h(ip)**3
                                    else
                                        rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                                        weight = fkern(rdist/h(ip))/h(ip)**2
                                    endif

                                    if ( btest(mode,MAX_POS) ) then
                                        if ( rdist<h(ip) ) then
                                            g(hix,hiy) = max(v(ip),g(hix,hiy))
                                        endif
                                    elseif ( btest(mode,VORINOI_POS) ) then
                                        if ( mg(hix,hiy)>rdist .and. ( rdist<=h(ip) .or. ih==1 ) ) then
                                            mg(hix,hiy)=rdist
                                            g(hix,hiy)=v(ip)
                                        endif
                                    else
                                        if ( btest(mode,DENSE_WEIGHT_POS) ) then
                                            g(hix,hiy) = g(hix,hiy) + weight * m(ip) * v(ip)
                                            if ( btest(mode,SDEV_POS) ) then
                                                g2(hix,hiy) = g2(hix,hiy) + weight * m(ip) * v(ip)**2
                                            endif
                                        endif
                                        mg(hix,hiy) = mg(hix,hiy) + weight * m(ip)
                                    endif
                                end do
                            end do

                        endif
                    endif

                endif
                
            endif
        end do
!$OMP END PARALLEL DO
        
        if ( .not. btest(mode,MAX_POS) ) then
            if ( btest(mode,MIN_POS) .or. btest(mode,VORINOI_POS) ) then
                where (g==huge(g(1,1))) g=-huge(g(1,1))
            else
                if ( btest(mode,DENSE_WEIGHT_POS) ) then
                    if ( btest(mode,SDEV_POS) ) then
                        g = sqrt(g2/mg-(g/mg)**2)
                    else
                        g = g/mg
                    endif
                else
                    g = mg
                endif
            endif
        endif
        
        deallocate(mg)
!         if ( btest(mode,SDEV_POS) ) then
            deallocate(g2)
!         endif

!         g(1,1) = count(isnan(v))

        return
    end function sph_general
    
    function sph_weight(x,y,m,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        
        g = sph_general(x,y,m,h,v,L,c,w,f,nonsense_array,0.d0,WEIGHT_MODE,n)

    end function sph_weight

    function sph_sdev(x,y,m,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        
        g = sph_general(x,y,m,h,v,L,c,w,f,nonsense_array,0.d0,ior(SDEV_MODE,WEIGHT_MODE),n)

    end function sph_sdev
    
    function sph_vorinoi(x,y,m,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        
        g = sph_general(x,y,m,h,v,L,c,w,f,nonsense_array,0.d0,VORINOI_MODE,n)

    end function sph_vorinoi
    
    function sph_max(x,y,m,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        
        g = sph_general(x,y,m,h,v,L,c,w,f,nonsense_array,0.d0,MAX_MODE,n)

    end function sph_max
    
    
    
    function sph_min(x,y,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        g = sph_general(x,y,nonsense_array,h,v,L,c,w,f,nonsense_array,0.d0,MIN_MODE,n)

    end function sph_min

    function sph_minslice(x,y,h,v,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid

        real(kind=8), dimension(n) :: nonsense_array

        
        g = sph_general(x,y,nonsense_array,h,v,L,c,w,f,nonsense_array,0.d0,ior(MIN_MODE,ZSLICE_MODE),n)

    end function sph_minslice
    
    function sph_dense(x,y,m,h,L,c,w,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        real(kind=8), dimension(n) :: nonsense_array

        g = sph_general(x,y,m,h,nonsense_array,L,c,w,f,nonsense_array,0.d0,DENSE_MODE,n)

    end function sph_dense
    
    
    function sph_dense_slice(x,y,m,h,L,c,w,z,zslice,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        real(kind=8) :: zslice
        
        real(kind=8), dimension(n) :: nonsense_array

        g = sph_general(x,y,m,h,nonsense_array,L,c,w,f,z,zslice,ior(DENSE_MODE,ZSLICE_MODE),n)

    end function sph_dense_slice

    function sph_weight_slice(x,y,m,h,v,L,c,w,z,zslice,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        real(kind=8) :: zslice

        g = sph_general(x,y,m,h,v,L,c,w,f,z,zslice,ior(WEIGHT_MODE,ZSLICE_MODE),n)

    end function sph_weight_slice

    function sph_vorinoi_slice(x,y,m,h,v,L,c,w,z,zslice,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        real(kind=8) :: zslice

        g = sph_general(x,y,m,h,v,L,c,w,f,z,zslice,ior(VORINOI_MODE,ZSLICE_MODE),n)

    end function sph_vorinoi_slice


    function sph_max_slice(x,y,m,h,v,L,c,w,z,zslice,f,n) result(g)
        implicit none
        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! values to smooth
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        real(kind=8) :: zslice

        g = sph_general(x,y,m,h,v,L,c,w,f,z,zslice,ior(MAX_MODE,ZSLICE_MODE),n)

    end function sph_max_slice


    function crossp(a,b) result(c)
        implicit none
        
        real(kind=8), dimension(3) :: a,b
        real(kind=8), dimension(3) :: c
        
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
        
    end function crossp

    function norm2crossp(a,b) result(norm2)
        implicit none
        
        real(kind=8), dimension(3) :: a,b
        real(kind=8) :: norm2
        
        norm2 =    (a(2)*b(3) - a(3)*b(2))**2 &
                 + (a(3)*b(1) - a(1)*b(3))**2 &
                 + (a(1)*b(2) - a(2)*b(1))**2
        
    end function norm2crossp


    ! histogram along some radial ray, taking opacity into account
    function sph_ray_histogram_opacity(xyz,m,h,vel,vmin,vmax,op,xyzray_in,rayoffset_in,f,broaden,nbins,nray,n) result(rayhist)
!     function sph_ray_histogram(xyz,m,h,v,vmin,vmax,xyzray_in,f,nbins,nray,n) result(rayhist)
        !$ use omp_lib
        implicit none
        
        integer :: n,nray,nbins
        
        real(kind=8), dimension(n,3), intent(in) :: xyz ! particle positions
        real(kind=8), dimension(n), intent(in) :: m,h ! "mass" = emission weight, smoothing
        real(kind=8), dimension(n,3), intent(in) :: vel ! 3D velocity
!         real(kind=8), dimension(n) :: broaden ! if <0, no broadening. If >0, width of gaussian broadening
        real(kind=8), dimension(n), intent(in) :: broaden ! width of gaussian broadening

        real(kind=8), dimension(n), intent(in) :: op ! particle cross section
        logical, dimension(n), intent(in) :: f ! mask
        
        
        real(kind=8), dimension(nray,3), intent(in) :: xyzray_in,rayoffset_in ! Ray directions and offset
        real(kind=8), dimension(nray,3) :: xyzray ! Ray from 0,0 in these directions - normalised
        
        real(kind=8), intent(in) :: vmin,vmax
        
        real(kind=8), dimension(nbins,nray) :: rayhist ! ray value along los (surface density)

        real(kind=8) :: ray_norm
        real(kind=8) :: h2 ! h squared to avoid square roots later
        real(kind=8) :: weight,broad_weight
        real(kind=8) :: dv2_normed, gauss_weight
!         real(kind=8) :: opac_weight,optical_depth

        real(kind=8), dimension(n) :: gauss_norm
        integer, dimension(:),allocatable :: zarg ! order along ray
        real(kind=8), dimension(:),allocatable :: z_ray ! position along ray
        real(kind=8) :: impact_pram

        real(kind=8) :: cum_depth,this_depth,vlos
        real(kind=8) :: min_vlos,max_vlos
        
        real(kind=8) :: bin_width
        
!         character(len=256) :: filename
!         integer :: ib
        


        
        integer :: iray,i,ip,ibin
        
        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
                
        ! 
        do iray=1,nray
            ray_norm = sqrt(sum(xyzray_in(iray,:)**2))
            xyzray(iray,:) = xyzray_in(iray,:)/ray_norm
        end do
        
        gauss_norm = broaden*sqrt(4.*atan(1.d0))
        
!         rayhist = 0.d0
        bin_width = (vmax-vmin)/nbins

!         print *,vmin,vmax,maxval(vel(:,2)),minval(vel(:,2))
!         min_vlos = 1.e9
!         max_vlos = -1.e9


!$OMP PARALLEL DO private(iray,ibin,i,ip,h2)&
!$OMP& private(z_ray,zarg)&
!$OMP& private(cum_depth,this_depth,vlos,dv2_normed)&
!$OMP& private(weight,gauss_weight,broad_weight,impact_pram)&
!$OMP& shared(nray,n,nbins,bin_width,vmin,vmax)&
!$OMP& shared(xyz,xyzray,rayhist,m,h,rayoffset_in,f,op,broaden,vel,gauss_norm)&
!$OMP& default(none) schedule(dynamic,1)
        do iray=1,nray
            rayhist(:,iray)=0.d0
            allocate(z_ray(n))
            allocate(zarg(n))
            cum_depth = 0.
            do ip=1,n
                z_ray(ip) = sum(xyz(ip,:)*xyzray(iray,:))
            end do
            call merge_argsort(z_ray,zarg)
            deallocate(z_ray)
!             write(filename,"('visualisation/test_dumps/ray',I3.3,'.dat')") iray
!             open(unit=15+iray,file=filename)
!             ib = 0
            !particleloop:  
            do i=1,n
                ip = zarg(i)
                if ( .not. f(ip) ) then
                    cycle! particleloop
                endif

!                 vlos = sum(vel(ip,:)*xyzray(iray,:))
!                 if ( vlos>50. .or. vlos<-50. ) then
!                     print *,vlos," tau ",cum_depth
!                 endif

!                 cum_depth = 0. ! hack to include all particles
                if ( cum_depth>25. ) then
                    exit! particleloop
                endif
                
                h2 = h(ip)**2
                impact_pram = norm2crossp(xyz(ip,:)-rayoffset_in(iray,:),xyzray(iray,:))
!                 if ( vlos>50. .or. vlos<-50. ) then
!                     print *,vlos," tau ",cum_depth,h(ip),sqrt(impact_pram)
!                     if (impact_pram>=h2) then
!                         write(*,"(16F10.6)") xyz(ip,:),rayoffset_in(iray,:),xyz(ip,:)-rayoffset_in(iray,:),&
!                                              xyzray(iray,:),vlos,vel(ip,:)
!                         stop
!                     endif
!                 endif
!                 impact_pram = 0. ! hack to include all particles
                
                if ( impact_pram>=h2 ) then
                    cycle! particleloop
                endif
!                 impact_pram = sqrt(impact_pram)
!                 weight = m(ip)*fkern(impact_pram/h(ip))/h2
                weight = m(ip)*fkern2(impact_pram/h2)/h2
                this_depth = cum_depth
                cum_depth = cum_depth+op(ip)*fkern(impact_pram/h(ip))/h2
                vlos = sum(vel(ip,:)*xyzray(iray,:))
!                 min_vlos = min(vlos,min_vlos)
!                 max_vlos = max(vlos,max_vlos)
                
                if ( broaden(ip)>bin_width ) then
                    do ibin=1,nbins
                        dv2_normed = ((vlos-vmin-(bin_width*(ibin-1)))/broaden(ip))**2
        !                 if ( dv2_normed>50. ) then
                        if ( dv2_normed>25. ) then
                            cycle! bin loop!
                        endif
                        gauss_weight = exp(-dv2_normed)/gauss_norm(ip)
                        broad_weight = weight * gauss_weight/bin_width
                        rayhist(ibin,iray) = rayhist(ibin,iray) + broad_weight*exp(-this_depth)
                    end do
!                     write(15+iray,"(9E15.7)") z_ray(ip),vel(ip,:),vlos,this_depth,cum_depth,impact_pram,h(ip)
                else
!                     print *,bin_width/broaden(ip)
                    ibin = int((vlos-vmin)/bin_width)+1
                    if ( ibin<1 .or. ibin>nbins ) then
!                         ib=ib+1
                        cycle ! particleloop
                    endif
                    rayhist(ibin,iray) = rayhist(ibin,iray) + weight/bin_width*exp(-this_depth)
!                     write(15+iray,"(9E15.7)") z_ray(ip),vel(ip,:),vlos,this_depth,cum_depth,impact_pram,h(ip)
                endif
                
            end do! particleloop
!             close(15+iray)
            deallocate(zarg)
!             print *,iray,min_vlos,max_vlos
!             print *,iray,"final optical depth=",cum_depth,"outside vrange=",ib,"/",n
        end do

    end function

    ! histogram along some radial ray
    function sph_ray_histogram(xyz,m,h,vel,vmin,vmax,xyzray_in,rayoffset_in,f,broaden,nbins,fullgal,nray,n) result(rayhist)
!     function sph_ray_histogram(xyz,m,h,v,vmin,vmax,xyzray_in,f,nbins,nray,n) result(rayhist)
        !$ use omp_lib
        implicit none
        
        integer :: n,nray,nbins
        
        real(kind=8), dimension(n,3) :: xyz ! particle positions
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
!         real(kind=8), dimension(n) :: v ! value to sum along ray
        real(kind=8), dimension(n,3), intent(in) :: vel ! 3D velocity
!         real(kind=8), dimension(n) :: broaden ! if <0, no broadening. If >0, width of gaussian broadening
        real(kind=8), dimension(n) :: broaden ! width of gaussian broadening
        logical :: fullgal ! true = (optically thin) emission lines from whole galaxy, false = absorption line from ray from offset
        
        real(kind=8), dimension(nray,3) :: xyzray_in,rayoffset_in ! Ray directions and offset
        real(kind=8), dimension(nray,3) :: xyzray ! Ray from 0,0 in these directions - normalised

        real(kind=8), allocatable, dimension(:,:) :: impact_pram_stored
        real(kind=8), dimension(nray) :: local_rayhist
        
        real(kind=8) :: vmin,vmax
        
        real(kind=8), dimension(nbins,nray) :: rayhist ! ray value along los (surface density)

        real(kind=8) :: ray_norm, dotprod
        real(kind=8) :: h2 ! h squared to avoid square roots later
        real(kind=8) :: weight,broad_weight
        real(kind=8) :: dv2_normed, gauss_weight
        real(kind=8), dimension(n) :: gauss_norm
        
        real(kind=8) :: local_broaden,local_impact_pram,local_v,local_m
        
        real(kind=8) :: bin_width
        

        logical, dimension(n) :: f ! mask

        
        integer :: iray,ip,ibin
        
        integer :: kk ! vector loop variable

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
                
        ! 
        do iray=1,nray
            ray_norm = sqrt(sum(xyzray_in(iray,:)**2))
            xyzray(iray,:) = xyzray_in(iray,:)/ray_norm
        end do
        
        gauss_norm = broaden*sqrt(4.*atan(1.d0))
        
        rayhist = 0.d0
        bin_width = (vmax-vmin)/nbins
        
        if ( .not. fullgal ) then
            allocate(impact_pram_stored(n,nray))
!             dotprod_stored = 0.

            do ip=1,n
                do iray=1,nray
!                     dotprod = 0.
!                     do kk=1,3
!                         dotprod = dotprod + xyzray(iray,kk)*(xyz(ip,kk)-rayoffset_in(iray,kk))
!                     end do
        
                    impact_pram_stored(ip,iray) = norm2crossp(xyz(ip,:)-rayoffset_in(iray,:),xyzray(iray,:))
                end do
            end do
        endif
        
!         print *,"integrating: MAX THREADS=",omp_get_max_threads(),nbins,nray,nbins*nray,nbins*nray/omp_get_max_threads()

!$OMP PARALLEL DO private(ibin,iray,dv2_normed,broad_weight)&
!$OMP& private(ip,weight,h2,gauss_weight,kk)&
!$OMP& private(local_rayhist,local_broaden,local_impact_pram,local_v,local_m)&
!$OMP& firstprivate(fullgal)&
!$OMP& shared(vel,vmin,bin_width,broaden,gauss_norm)&
!$OMP& shared(rayoffset_in,rayhist,nray,nbins,xyzray,xyz,h,m,n,f)&
!$OMP& shared(impact_pram_stored)&
!$OMP& default(none) schedule(dynamic,1)
        do ibin=1,nbins
!             print *,ibin,omp_get_thread_num()
            local_rayhist = 0.
            do ip=1,n
!             do ip=1,20000
!                 if ( mod(ip,1000)==0 .and. ibin==1) then
!                     print *,ip,"/",n
!                 endif
                if ( .not. f(ip) ) then
                    cycle
                endif
                local_broaden = broaden(ip)
!                 local_v = v(ip)

!                 else
!                     if ( local_v<vmin+bin_width*(ibin-1) ) then
!                         cycle
!                     endif
!                     if ( local_v>vmin+bin_width*(ibin) ) then
!                         cycle
!                     endif
!                     gauss_weight = 1.
!                 endif
                
                local_m = m(ip)
                
                if ( .not. fullgal ) then
                    h2 = h(ip)**2
                    do iray=1,nray
                        local_v = sum(vel(ip,:)*xyzray(iray,:))
    !                 if ( local_broaden>bin_width ) then
                        dv2_normed = ((local_v-vmin-(bin_width*(ibin-1)))/local_broaden)**2
        !                 if ( dv2_normed>50. ) then
                        if ( dv2_normed>25. ) then
                            cycle
                        endif
                        gauss_weight = exp(-dv2_normed)/gauss_norm(ip)

                        local_impact_pram = impact_pram_stored(ip,iray)
                        if ( local_impact_pram>=h2 ) then
                            cycle
                        endif
!                         weight = local_m*fkern2(local_impact_pram/h2)/h2
                        weight = local_m*fkern2(local_impact_pram/h2)/h2
                        broad_weight = weight * gauss_weight/bin_width
                        local_rayhist(iray) = local_rayhist(iray) + broad_weight
                    end do
                else
                    do iray=1,nray
                        local_v = sum(vel(ip,:)*xyzray(iray,:))
                        dv2_normed = ((local_v-vmin-(bin_width*(ibin-1)))/local_broaden)**2
                        if ( dv2_normed>25. ) then
                            cycle
                        endif
                        gauss_weight = exp(-dv2_normed)/gauss_norm(ip)

                        broad_weight = local_m * gauss_weight/bin_width
                        local_rayhist(iray) = local_rayhist(iray) + broad_weight
                    end do
                endif


!                 do iray=1,nray
!                     if ( .not. fullgal ) then
! !                         dotprod = 0.
! !                         do kk=1,3
! ! !                         dotprod = sum(xyzray(iray,:)*(xyz(ip,:)-rayoffset_in(iray,:)))
! !                             dotprod = dotprod + xyzray(iray,kk)*(xyz(ip,kk)-rayoffset_in(iray,kk))
! !                         end do
! !                         if ( dotprod<=0. ) then
! !                             cycle
! !                         endif
!                 
! !                         impact_pram = norm2crossp(xyz(ip,:)-rayoffset_in(iray,:),xyzray(iray,:))
!                         local_impact_pram = impact_pram_stored(ip,iray)
!                         if ( local_impact_pram>=h2 ) then
!                             cycle
!                         endif
!                 
! !                         impact_pram = sqrt(impact_pram_stored(ip,iray))
!                         weight = m(ip)*fkern2(local_impact_pram/h2)/h2
!                     else
!                         weight = m(ip) ! m is actually luminosity
!                     endif
!                             ! should we also weight by bin width?
!                     broad_weight = weight * gauss_weight/bin_width
! !                     rayhist(ibin,iray) = rayhist(ibin,iray) + broad_weight
!                     local_rayhist(iray) = local_rayhist(iray) + broad_weight
!                 end do
            end do
            rayhist(ibin,:) = local_rayhist
        end do
!$OMP END PARALLEL DO

        if ( .not. fullgal ) then
!             deallocate(dotprod_stored)
            deallocate(impact_pram_stored)
        endif
    end function


    function sph_ray_histogram_2(xyz,m,h,v,vmin,vmax,xyzray_in,rayoffset_in,f,broaden,nbins,fullgal,nray,n) result(rayhist)
!     function sph_ray_histogram(xyz,m,h,v,vmin,vmax,xyzray_in,f,nbins,nray,n) result(rayhist)
        !$ use omp_lib
        implicit none
        
        integer :: n,nray,nbins
        
        real(kind=8), dimension(n,3) :: xyz ! particle positions
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! value to sum along ray
!         real(kind=8), dimension(n) :: broaden ! if <0, no broadening. If >0, width of gaussian broadening
        real(kind=8), dimension(n) :: broaden ! width of gaussian broadening
        logical :: fullgal ! true = (optically thin) emission lines from whole galaxy, false = absorption line from ray from offset

        
        real(kind=8), dimension(nray,3) :: xyzray_in,rayoffset_in ! Ray directions and offset
        real(kind=8), dimension(nray,3) :: xyzray ! Ray from 0,0 in these directions - normalised
        
        real(kind=8) :: vmin,vmax
        
        real(kind=8), dimension(nbins,nray) :: rayhist ! ray value along los (surface density)

        real(kind=8) :: ray_norm, impact_pram, dotprod
        real(kind=8) :: h2 ! h squared to avoid square roots later
        real(kind=8) :: weight,broad_weight
        real(kind=8) :: dv2_normed, gauss_weight
        real(kind=8), dimension(n) :: gauss_norm
        
        real(kind=8) :: bin_width
        
        integer :: index
        

        logical, dimension(n) :: f ! mask

        
        integer :: iray,ip,ibin

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        ! 
        do iray=1,nray
            ray_norm = sqrt(sum(xyzray_in(iray,:)**2))
            xyzray(iray,:) = xyzray_in(iray,:)/ray_norm
        end do
        
        gauss_norm = broaden*sqrt(4.*atan(1.d0))
        
        rayhist = 0.d0
        bin_width = (vmax-vmin)/nbins
        
!         print *,"MAX THREADS=",omp_get_max_threads()

!$OMP PARALLEL DO private(ibin,iray,dv2_normed,broad_weight,dotprod)&
!$OMP& private(impact_pram,ip,weight,h2,gauss_weight,index)&
!$OMP& shared(v,vmin,bin_width,broaden,gauss_norm,fullgal)&
!$OMP& shared(rayoffset_in,rayhist,nray,nbins,xyzray,xyz,h,m,n,f)&
!$OMP& default(none) schedule(dynamic,1)
        do index=0,nbins*nray-1
            ibin = mod(index,nbins)+1
            iray = index/nbins+1
            do ip=1,n
                if ( .not. f(ip) ) then
                    cycle
                endif
                if ( .not. fullgal ) then
                    h2 = h(ip)**2
                endif

                if ( .not. fullgal ) then
                    dotprod = sum(xyzray(iray,:)*(xyz(ip,:)-rayoffset_in(iray,:)))
                    if ( dotprod<=0. ) then
                        cycle
                    endif
            
                    impact_pram = norm2crossp(xyz(ip,:)-rayoffset_in(iray,:),xyzray(iray,:))
                    if ( impact_pram>=h2 ) then
                        cycle
                    endif
            
                    impact_pram = sqrt(impact_pram)
                    weight = m(ip)*fkern(impact_pram/h(ip))/h2
                else
                    weight = m(ip) ! m is actually luminosity
                endif
                
                if ( broaden(ip)>bin_width ) then
                    dv2_normed = ((v(ip)-vmin-(bin_width*(ibin-1)))/broaden(ip))**2
    !                 if ( dv2_normed>50. ) then
                    if ( dv2_normed>25. ) then
                        cycle
                    endif
                    gauss_weight = exp(-dv2_normed)/gauss_norm(ip)
                else
                    if ( v(ip)<vmin+bin_width*(ibin-1) ) then
                        cycle
                    endif
                    if ( v(ip)>vmin+bin_width*(ibin) ) then
                        cycle
                    endif
                    gauss_weight = 1.
                endif
                        ! should we also weight by bin width?
                broad_weight = weight * gauss_weight/bin_width
                rayhist(ibin,iray) = rayhist(ibin,iray) + broad_weight

            end do
        end do
!$OMP END PARALLEL DO
    end function

    ! surface density along some radial ray
    function sph_ray_integrate(xyz,m,h,xyzray_in,nray,n) result(rw)
        !$ use omp_lib
        implicit none
        
        integer :: n,nray
        
        real(kind=8), dimension(n,3) :: xyz ! particle positions
        real(kind=8), dimension(n) :: m,h ! mass, smoothing

        
        real(kind=8), dimension(nray,3) :: xyzray_in ! Ray from 0,0 in these directions
        real(kind=8), dimension(nray,3) :: xyzray ! Ray from 0,0 in these directions - normalised
        
        real(kind=8), dimension(nray) :: rw ! ray value along los (surface density)

        real(kind=8) :: ray_norm, impact_pram, dotprod
        real(kind=8) :: h2 ! h squared to avoid square roots later
        real(kind=8) :: weight
        
        integer :: iray,ip

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        do iray=1,nray
            ray_norm = sqrt(sum(xyzray_in(iray,:)**2))
            xyzray(iray,:) = xyzray_in(iray,:)/ray_norm
        end do
        
        print *,"Integrating:",n," particles ",nray," rays"
        print *,"MAX THREADS=",omp_get_max_threads()
        
        rw = 0.d0
        !$OMP PARALLEL DO private(iray,dotprod,impact_pram,weight,h2,ip)&
        !$OMP& shared(nray,xyzray,xyz,rw,h,m,n) default(none)
        do iray=1,nray
            do ip=1,n
!                 if ( omp_get_thread_num()==0 .and. mod(ip,10000)==1 ) then
!                     print *,ip,"/",n
!                 endif
                h2 = h(ip)**2
    !                 print *,omp_get_thread_num(),ip,iray
                dotprod = sum(xyzray(iray,:)*xyz(ip,:))
                if ( dotprod>0. ) then
                    impact_pram = norm2crossp(xyz(ip,:),xyzray(iray,:))
                    if ( impact_pram<h2 ) then
                        impact_pram = sqrt(impact_pram)
                        weight = m(ip)*fkern(impact_pram/h(ip))/h2
                    
                        rw(iray) = rw(iray) + weight
                    endif
                endif
            end do
        end do
        !$OMP END PARALLEL DO 

! 
!         do ip=1,n
! !             if ( mod(ip,10000)==1 ) then
!                 print *,ip,"/",n
! !             endif
!             h2 = h(ip)**2
!             !$OMP PARALLEL DO private(iray,dotprod,impact_pram,weight)&
!             !$OMP& shared(nray,xyzray,xyz,rw,h2,h,m,ip) default(none)
!             do iray=1,nray
! !                 print *,omp_get_thread_num(),ip,iray
!                 dotprod = sum(xyzray(iray,:)*xyz(ip,:))
!                 if ( dotprod>0. ) then
!                     impact_pram = norm2crossp(xyz(ip,:),xyzray(iray,:))
!                     if ( impact_pram<h2 ) then
!                         impact_pram = sqrt(impact_pram)
!                         weight = m(ip)*fkern(impact_pram/h(ip))/h2
!                         
!                         rw(iray) = rw(iray) + weight
!                     endif
!                 endif
!             end do
!             !$OMP END PARALLEL DO 
!         end do
        print *,"Integration complete"
        
    end function

    ! quick dotplot
    function sph_dot(x,y,v,L,c,w,overmode,f,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: v ! particle value
        real(kind=8), dimension(n) :: x,y ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        integer :: overmode ! 0 = max, 1 = min

        integer :: ip ! particle index

        integer :: ix,iy ! grid position of particle
        
        real(kind=8) :: r_cell, area_cell

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = -HUGE(v(ip))
! 
!         if ( overmode==0 ) then
!             g = -HUGE(v(ip))
!         else
!             g = HUGE(v(ip))
!         endif
        
        r_cell = w/L
        area_cell = r_cell**2

        do ip=1,n
            if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
                ix = nint((x(ip)-c(1))/r_cell)
                iy = nint((y(ip)-c(2))/r_cell)
                
                if ( ix<=L .and. iy<=L .and. ix>=1 .and. iy>=1 ) then
                    if ( overmode==0 ) then
                        g(ix,iy) = max(v(ip),g(ix,iy))
                    else
                        if ( g(ix,iy)==-HUGE(v(ip)) ) then
                            g(ix,iy)=v(ip)
                        else
                            g(ix,iy) = min(v(ip),g(ix,iy))
                        endif
                    endif

                endif
                
            endif
        end do

    end function
    
    ! Imaging
    function sph_optical_depth_los(x,y,m,h,v,op,L,c,w,z,inzarg,f,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! brightness of particle
        real(kind=8), dimension(n) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        integer, dimension(n) :: inzarg,zarg ! sorted positions of particles along line of sight
!         real(kind=8), dimension(L,L) :: opg ! grid of optical depths

        integer :: ip ! loop variable, particle index

        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
!         integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, planedist2, weight!, this_opac
        
        integer, parameter :: n_zdir_max=100000
        integer, parameter :: persamples=16
        integer :: iz,nz
        real(kind=8), dimension(n_zdir_max) :: binedges,opacgrid,emitgrid!,massgrid
        logical, dimension(n) :: inray
        !real(kind=8),allocatable,dimension(:,:,:) :: opac_samples,emit_samples! FIXXXXX
        real(kind=8) :: stepsize
        
        
        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = 0.d0
!         opg = 0.d0
        r_cell = w/L
        area_cell = r_cell**2

        !print *,"sorting"
        !zarg = rargsort(z)
        !print *,"sorted"
        zarg = inzarg+1 ! python to fortran numbering

        ! flip it, for backwards compatability
        z = maxval(v)-z ! lololol
        
!         print *,binedges(1:nz)
        ! parallelise this loop
!$OMP PARALLEL DO &
!$OMP& private(ix,iy,inray,ip,ih,hix,hiy,nz,binedges,opacgrid,emitgrid,iz,stepsize,planedist2,weight,rdist)&
!$OMP& shared(g,v,op,m,n,c,r_cell,f,L,h,x,y,z) default(none)
        do ix=1,L
            print *,ix,"/",L,"started"
            do iy=1,L
!                 print *,ix,"/",L,iy,"/",L," ray"
                inray = .false.
                do ip=1,n
                    if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
                        ih = ceiling(h(ip)/r_cell)
                        hix = nint((x(ip)-c(1))/r_cell)
                        hiy = nint((y(ip)-c(2))/r_cell)
                
!                         ix0 = hix-ih)
!                         iy0 = max(1,hiy-ih)
! 
!                         ix1 = min(L,hix+ih)
!                         iy1 = min(L,hiy+ih)
                
                        if (  hix+ih>=ix .and. hix-ih<=ix.and. hiy+ih>=iy .and. hiy-ih<=iy  ) then
                            inray(ip) = .true.
!                             print *,z(ip),h(ip)/persamples
                        endif
                    endif
                end do
!                 print *,"n=",count(inray)
!                 print *,pack(z,inray)
                nz = persamples
                binedges(1) = minval(z-h,inray)
                binedges(nz) = maxval(z+h,inray)
                opacgrid(1:nz)=0. 
!                 massgrid(1:nz)=0.
                emitgrid(1:nz)=0.
                do iz=2,nz-1
                    binedges(iz) = (binedges(nz)-binedges(1))*(iz-1.)/(nz-1)+binedges(1)
                end do
!                 print *,minval(h,inray)
!                 print *,binedges(nz)-binedges(1)
!                 print *,(binedges(nz)-binedges(1))/minval(h,inray)*persamples

                do ip=1,n
!                     if ( mod(ip,10000)==1 ) then
!                         print *,ip,"/",n,nz
!                     endif
                    if ( .not.inray(ip) ) then
                        cycle
                    endif
                    stepsize = h(ip)/persamples
        !             print *,stepsize
                    iz = 1
                    planedist2 = (x(ip)-(ix)*r_cell-c(1))**2+ &
                                 (y(ip)-(iy)*r_cell-c(1))**2
                    do while (iz<nz .and. binedges(iz)<=z(ip)+h(ip) )
                        if ( binedges(iz+1)>=z(ip)-h(ip) ) then
                            if ( binedges(iz+1)-binedges(iz)>stepsize ) then
                                ! halve size - refine
                                if ( nz+1 > n_zdir_max ) then
                                    print *,"over-refined",ip,nz
                                    stop
                                endif
!                                 print *,ip,inray(ip)
!                                 print *,binedges(1:nz)
                                ! shuffle up edges by one, insert one half-way in-between
                                binedges(iz+2:nz+1)=binedges(iz+1:nz)
                                binedges(iz+1) = (binedges(iz)+binedges(iz+1))/2.
                                ! for SPH values, these are all volumetric values, just copy them over
                                ! values iz+1 and iz+2 are now identical - i.e. below the res we cared about
                                emitgrid(iz+2:nz+1)=emitgrid(iz+1:nz) 
                                opacgrid(iz+2:nz+1)=opacgrid(iz+1:nz) 
!                                 massgrid(iz+2:nz+1)=massgrid(iz+1:nz) 
                                nz = nz + 1
!                                 print *,binedges(1:nz)
!                                 stop
                            else
                                ! we are within range, add our SPH weights to this point
                                rdist=sqrt(planedist2+(z(ip)-binedges(iz))**2)
                                if ( rdist<h(ip) ) then
                                    weight = kern(rdist/h(ip))/h(ip)**3
                                    emitgrid(iz)=emitgrid(iz)+v(ip)*op(ip)*weight*m(ip)
                                    opacgrid(iz)=opacgrid(iz)+op(ip)*weight*m(ip)
!                                     massgrid(iz)=massgrid(iz)+weight*m(ip)
                                endif
                                iz = iz + 1
                            endif
                        else
                            ! out of range, skip
                            iz = iz + 1
                        endif
                    end do
                end do

                ! splat it down flat
                weight = 0.
                do iz=1,nz-1
        !             write(12,"(I8,5E9.3)") iz,binedges(iz),binedges(iz+1)-binedges(iz),opacgrid(iz),emitgrid(iz),massgrid(iz)
                    weight = weight + (binedges(iz+1)-binedges(iz))*(emitgrid(iz))
                    weight = weight*exp(-(binedges(iz+1)-binedges(iz))*opacgrid(iz))
!                     write(12,*) iz,binedges(iz),binedges(iz+1)-binedges(iz),opacgrid(iz),emitgrid(iz),massgrid(iz),weight
                end do
                g(ix,iy) = weight

            end do
            print *,ix,"/",L,"done"
        end do
!         opacgrid(1:nz)=opacgrid(1:nz)!/massgrid(1:nz) ! 
!         emitgrid(1:nz)=emitgrid(1:nz)/massgrid(1:nz)
!         open(unit=12,file="test.out")
!         weight = 0.
!         do iz=1,nz-1
! !             write(12,"(I8,5E9.3)") iz,binedges(iz),binedges(iz+1)-binedges(iz),opacgrid(iz),emitgrid(iz),massgrid(iz)
!             weight = weight + (binedges(iz+1)-binedges(iz))*(emitgrid(iz)-opacgrid(iz)*weight)
!             write(12,*) iz,binedges(iz),binedges(iz+1)-binedges(iz),opacgrid(iz),emitgrid(iz),massgrid(iz),weight
!         end do
! !         print *,binedges(1:nz)
!         close(12)
!         stop
!         do ip=1,n
!             if ( mod(ip,1000)==1 ) then
!                 print *,ip,"/",n,nz
!             endif
!             if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
!                 ih = ceiling(h(ip)/r_cell)
!                 ix = nint((x(ip)-c(1))/r_cell)
!                 iy = nint((y(ip)-c(2))/r_cell)
!                 
!                 ix0 = max(1,ix-ih)
!                 iy0 = max(1,iy-ih)
! 
!                 ix1 = min(L,ix+ih)
!                 iy1 = min(L,iy+ih)
!                 
!                 if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
!                 
!                     do hix=ix0,ix1
!                         do hiy=iy0,iy1
!                             iz = 1
!                             do while ( iz<nz .and. binedges(iz)<=z(ip)+h(ip))
!                                 if ( binedges(iz)>=z(ip)+h(ip) ) then
!                                     rdist=sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell+(z(ip)-h(ip)**2))
!                                     weight = kern(rdist/h(ip))/kern_norm/h(ip)**3
!                                     
!                                 endif
!                                 iz = iz + 1
!                             end do
!                         end do
!                     end do
!                 endif
!             endif
! 
!         end do
        
        
!         open(unit=12,file="test.out")
!         do iz=1,nz-1
!             write(12,*) iz,binedges(iz+1)-binedges(iz)
!         end do
! !         print *,binedges(1:nz)
!         close(12)
!         stop
!         i_zedges = 0
!         do i=1,n
!             ip = zarg(i)
!             if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
!                 if ( i_zedges==0 ) then
!                     do j=0,persamples-1
!                         z_edges(j) = (j-persamples/2)*h(ip)/persamples+z_p(j)
!                     end do
!                 else
!                     
!                 endif
!             endif
!         end do
!         stop



! old version
!         do i=1,n
!             ip = zarg(i)
!             if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
!                 ih = ceiling(h(ip)/r_cell)
!                 ix = nint((x(ip)-c(1))/r_cell)
!                 iy = nint((y(ip)-c(2))/r_cell)
!                 
!                 ix0 = max(1,ix-ih)
!                 iy0 = max(1,iy-ih)
! 
!                 ix1 = min(L,ix+ih)
!                 iy1 = min(L,iy+ih)
!                 
!                 if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
!                 
!                     do hix=ix0,ix1
!                         do hiy=iy0,iy1
!                             !if ( opg(hix,hiy)<1.d0 ) then
!                                 rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
!                                 weight = fkern(rdist/h(ip))/h(ip)**2
!                             
! 
!                                 this_opac = weight * m(ip) * op(ip)
!                                 g(hix,hiy) = g(hix,hiy) + min(1.,this_opac)*v(ip)*exp(-opg(hix,hiy))
!                                 !g(hix,hiy) = g(hix,hiy) + v(ip)*exp(-opg(hix,hiy)) ! is this maybe correct?
!                                 opg(hix,hiy) = opg(hix,hiy) + this_opac
! !                                 if ( opg(hix,hiy)>=1.d0 ) then
! !                                     g(hix,hiy)=v(ip)
! !                                 endif
!                             !endif
!                             
!                         end do
!                     end do
! 
!                 endif
!                 
!             endif
!         end do

    end function



! basically, assume that dust is both producing and extincting
! v is emission in erg/s/cm^2 where /cm^2 is per surface area of dust
! op is opacity in cm^2/g where g is mass of total gas, and cm^2 is cross section for both emission and absorption
! op*(line of sight column density) = fraction of area covered by dust. Take max of 1 for 100% coverage
! multiply this by v to get erg/s/cm^2 column contribution
    function sph_optical_depth_los_area(x,y,m,h,v,op,L,c,w,z,inzarg,f,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! brightness of particle
        real(kind=8), dimension(n) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        integer, dimension(n) :: inzarg,zarg ! sorted positions of particles along line of sight
        real(kind=8), dimension(L,L) :: opg ! grid of optical depths

        integer :: i,ip ! loop variable, particle index

        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, weight, this_opac

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = 0.d0
        opg = 0.d0
        r_cell = w/L
        area_cell = r_cell**2

        !print *,"sorting"
        !zarg = rargsort(z)
        !print *,"sorted"
        zarg = inzarg+1 ! python to fortran numbering

        do i=1,n
            ip = zarg(i)
            if ( .not. f(ip) .or. isnan(v(ip)) .or. v(ip)>HUGE(v(ip)) .or. v(ip)<-HUGE(v(ip)) ) cycle
            ih = ceiling(h(ip)/r_cell)
            ix = nint((x(ip)-c(1))/r_cell)
            iy = nint((y(ip)-c(2))/r_cell)
            
            ix0 = max(1,ix-ih)
            iy0 = max(1,iy-ih)

            ix1 = min(L,ix+ih)
            iy1 = min(L,iy+ih)
            
            if ( ix0>L .or. iy0>L .or. ix1<1 .or. iy1<1 ) cycle
            
            do hix=ix0,ix1
                do hiy=iy0,iy1
                    !if ( opg(hix,hiy)<1.d0 ) then
                        rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                        weight = fkern(rdist/h(ip))/h(ip)**2 ! column density (normalised) of this ray through particle
                        ! m(ip) * weight = mass column density of this ray through particle

                        this_opac = weight * m(ip) * op(ip) ! optical depth of this ray through particle
                        g(hix,hiy) = g(hix,hiy) + min(1.,this_opac)*v(ip)*exp(-opg(hix,hiy))
                        !g(hix,hiy) = g(hix,hiy) + v(ip)*exp(-opg(hix,hiy)) ! is this maybe correct?
                        opg(hix,hiy) = opg(hix,hiy) + this_opac
!                                 if ( opg(hix,hiy)>=1.d0 ) then
!                                     g(hix,hiy)=v(ip)
!                                 endif
                    !endif
                    
                end do
            end do

        end do
    end function

! `emit` is total luminosity of particle
! `emit` weighted by 2D kernel gives erg/s/cm^2 flux
! weight by `v` to calculate "line" output, divide at the end
! opacity is completely independent, but still mass weighted
    function sph_optical_depth_los_weight(x,y,m,h,v,emit,op,L,c,w,z,inzarg,f,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! value to smooth over
        real(kind=8), dimension(n) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n) :: emit ! particle emission per unit mass
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        integer, dimension(n) :: inzarg ! sorted positions of particles along line of sight

        
        g = sph_optical_depth_los_weight_thresh(x,y,m,h,v,emit,op,L,c,w,z,inzarg,f,0.d0,n)
        return
    end function

! `emit` is total luminosity of particle
! `emit` weighted by 2D kernel gives erg/s/cm^2 flux
! weight by `v` to calculate "line" output, divide at the end
! opacity is completely independent, but still mass weighted
    function sph_optical_depth_los_weight_thresh(x,y,m,h,v,emit,op,L,c,w,z,inzarg,f,threshold,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! value to smooth over
        real(kind=8), dimension(n) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n) :: emit ! particle emission per unit mass
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        integer, dimension(n) :: inzarg,zarg ! sorted positions of particles along line of sight
        real(kind=8), dimension(L,L) :: opg, emitg ! grid of optical depths, grid of emissions
        
        real(kind=8) :: threshold

        integer :: i,ip ! loop variable, particle index

        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, weight, this_opac, this_emission

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = 0.d0
        opg = 0.d0
        emitg = 0.d0
        r_cell = w/L
        area_cell = r_cell**2

        !print *,"sorting"
        !zarg = rargsort(z)
        !print *,"sorted"
        zarg = inzarg+1 ! python to fortran numbering

        do i=1,n
            ip = zarg(i)
            if ( .not. f(ip) .or. isnan(v(ip)) .or. v(ip)>HUGE(v(ip)) .or. v(ip)<-HUGE(v(ip)) ) cycle
            ih = ceiling(h(ip)/r_cell)
            ix = nint((x(ip)-c(1))/r_cell)
            iy = nint((y(ip)-c(2))/r_cell)
            
            ix0 = max(1,ix-ih)
            iy0 = max(1,iy-ih)

            ix1 = min(L,ix+ih)
            iy1 = min(L,iy+ih)
            
            if ( ix0>L .or. iy0>L .or. ix1<1 .or. iy1<1 ) cycle
            
            do hix=ix0,ix1
                do hiy=iy0,iy1
                    !if ( opg(hix,hiy)<1.d0 ) then
                        rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                        weight = fkern(rdist/h(ip))/h(ip)**2
                    

                        this_opac = weight * m(ip) * op(ip)
                        this_emission = weight * emit(ip)
                        g(hix,hiy) = g(hix,hiy) + this_emission*v(ip)*exp(-opg(hix,hiy))
                        !g(hix,hiy) = g(hix,hiy) + v(ip)*exp(-opg(hix,hiy)) ! is this maybe correct?
                        emitg(hix,hiy) = emitg(hix,hiy) + this_emission*exp(-opg(hix,hiy))
                        opg(hix,hiy) = opg(hix,hiy) + this_opac
!                                 if ( opg(hix,hiy)>=1.d0 ) then
!                                     g(hix,hiy)=v(ip)
!                                 endif
                    !endif
                    
                end do
            end do

        end do
        
        if ( threshold>0 ) then
        
            do ix=1,L
                do iy=1,L
                    if ( emitg(ix,iy)<threshold ) then
                        g(ix,iy)=0.
                        emitg(ix,iy)=0.
                    endif
                end do
            end do

        endif
        
        g = g / emitg
        return
    end function

    function sph_mean_optical_depth(x,y,z,m,h,op,L,c,w,f,tau,n) result(opg)
        implicit none

        integer, intent(in) :: n,L
        real(kind=8), dimension(n), intent(in) :: m,h ! mass, smoothing
        real(kind=8), dimension(n), intent(in) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n),intent(in) :: x,y,z ! particle positions - z is los
        logical, dimension(n), intent(in) :: f ! mask
        real(kind=8), dimension(2), intent(in) :: c ! top-left corner
        real(kind=8), intent(in) :: w ! width

        real(kind=8), dimension(n) :: tau ! output optical depths



        integer, dimension(n) :: zarg ! sorted positions of particles along line of sight
        real(kind=8), dimension(L,L) :: opg ! grid of optical depths

        integer :: i,ip ! loop variable, particle index

        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, weight, this_opac
        real(kind=8) :: weighted_tau, weight_sum

        if ( .not. kernel_initialized ) then
            call kernel_init
        endif

        opg = 0.d0
        r_cell = w/L
        area_cell = r_cell**2
        
        call merge_argsort(z,zarg)

        tau=0.
        do i=1,n
            ip = zarg(i)
            if ( .not. f(ip) ) cycle
            ih = ceiling(h(ip)/r_cell)
            ix = int((x(ip)-c(1))/r_cell)+1
            iy = int((y(ip)-c(2))/r_cell)+1
            
            ix0 = max(1,ix-ih)
            iy0 = max(1,iy-ih)

            ix1 = min(L,ix+ih)
            iy1 = min(L,iy+ih)
            
            if ( ix0>L .or. iy0>L .or. ix1<1 .or. iy1<1 ) cycle
            
            weighted_tau=0.
            weight_sum=0.
            do hix=ix0,ix1
                do hiy=iy0,iy1
                    rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                    weight = fkern(rdist/h(ip))/h(ip)**2
                    weight_sum = weight_sum + weight*area_cell
            
                    this_opac = weight * m(ip) * op(ip)
                    weighted_tau=weighted_tau+weight*area_cell*opg(hix,hiy)

                    opg(hix,hiy) = opg(hix,hiy) + this_opac
!                     print *,ix,iy,hix,hiy,x(ip),y(ip),rdist,h(ip)
                    
                end do
            end do
            tau(ip)=weighted_tau/weight_sum
!             print *,i,z(ip),tau(ip)
        end do
    end function



! I think this is a selection sort, so it's slow
! taken from:
! Module for sorting arrays.
! Based on code written by John E. Pask, LLNL.
! https://github.com/certik/fortran-utils/blob/master/src/sorting.f90


    function rargsort(a) result(b)
        ! Returns the indices that would sort an array.
        !
        ! Arguments
        ! ---------
        !
        real(kind=8), intent(in):: a(:)   ! array of numbers
        integer :: b(size(a))         ! indices into the array 'a' that sort it
        !
        ! Example
        ! -------
        !
        ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

        integer :: N                           ! number of numbers/vectors
        integer :: i,imin                      ! indices: i, i of smallest
        integer :: temp1                       ! temporary
        real(kind=8) :: temp2
        real(kind=8) :: a2(size(a))
        a2 = a
        N=size(a)
        do i = 1, N
            b(i) = i
        end do
        do i = 1, N-1
            ! find ith smallest in 'a'
            imin = minloc(a2(i:),1) + i - 1
            ! swap to position i in 'a' and 'b', if not already there
            if (imin /= i) then
                temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
                temp1 = b(i); b(i) = b(imin); b(imin) = temp1
            end if
        end do
    end function

! my own sorting subroutine - should be faster
    subroutine merge_argsort(r,d)
        real(kind=8), intent(in), dimension(:) :: r
        integer, intent(out), dimension(size(r)) :: d
    
        integer, dimension(size(r)) :: il

        integer :: stepsize
        integer :: i,j,n,left,k,ksize
    
        n = size(r)
    
        do i=1,n
            d(i)=i
        end do
    
        if ( n==1 ) return
    
        stepsize = 1
        do while (stepsize<n)
            do left=1,n-stepsize,stepsize*2
                i = left
                j = left+stepsize
                ksize = min(stepsize*2,n-left+1)
                k=1
        
                do while ( i<left+stepsize .and. j<left+ksize )
                    if ( r(d(i))>r(d(j)) ) then
                        il(k)=d(i)
                        i=i+1
                        k=k+1
                    else
                        il(k)=d(j)
                        j=j+1
                        k=k+1
                    endif
                enddo
        
                if ( i<left+stepsize ) then
                    ! fill up remaining from left
                    il(k:ksize) = d(i:left+stepsize-1)
                else
                    ! fill up remaining from right
                    il(k:ksize) = d(j:left+ksize-1)
                endif
                d(left:left+ksize-1) = il(1:ksize)
            end do
            stepsize=stepsize*2
        end do

        return
      end subroutine    

    
        ! 
!         integer :: ip ! this particle
!         integer :: ix,iy ! grid position of particle
!         integer :: hix,hiy ! position in h circle
!         integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
!         
!         integer :: ih ! integer h
!         
!         real(kind=8) :: r_cell, area_cell
!         
!         real(kind=8) :: rdist, weight
!         
!         if ( .not. kernel_initialized ) then
!             call kernel_init
!         endif
!         
!         g = 0.
!         r_cell = w/L
!         area_cell = r_cell**2
!         
!         do ip=1,n
!             if ( f(ip) ) then
!                 ih = ceiling(h(ip)/r_cell)
!                 ix = nint((x(ip)-c(1))/r_cell)
!                 iy = nint((y(ip)-c(2))/r_cell)
!                 
!                 ix0 = max(1,ix-ih)
!                 iy0 = max(1,iy-ih)
! 
!                 ix1 = min(L,ix+ih)
!                 iy1 = min(L,iy+ih)
!                 
!                 if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
!                 
!                     do hix=ix0,ix1
!                         do hiy=iy0,iy1
!                             rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
!                             weight = fkern(rdist/h(ip))*area_cell/h(ip)**2
!                             g(hix,hiy) = g(hix,hiy) + weight * m(ip)
!                         end do
!                     end do
! 
!                 endif
!                 
!             endif
!         end do
!         
!         return
!     end function sph_dense
    
    subroutine kernel_init
        implicit none
        
        integer :: ir,iz
        real(kind=8) :: r,z,x
        integer, parameter :: zsteps = 1000
        
        real(kind=8) :: norm
        
!         print *,"initializing kernel"
        
        kern_tab = 0.
!        sum = 0.
        
        do ir=1,nkern
            r = ((ir-1.)/(nkern-1.))
            do iz=-zsteps,zsteps-1
                z = (iz+.5)/zsteps
                x = sqrt(r**2+z**2)
                kern_tab(ir) = kern_tab(ir) + kern(x)/zsteps
            end do
        end do

        ! normalise flattened kernel
        norm = 0.d0
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            !norm = norm+kern_tab(ir)/(nkern-1.)*3.14159265359d0*r
            norm = norm+kern_tab(ir)/(nkern-1.)*3.14159265359d0*r
!            write(*,"(I5,3E15.5)") ir,r,kern_tab(ir),sum
        end do
        kern_tab(:) = kern_tab(:)/norm ! Normalise so it all sums to 1.

        
        ! normalise 3d kernel
        kern_norm = 0.
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            kern_norm = kern_norm+fkern(r)/(nkern-1.)*3.14*r**2
        end do

        do ir=1,nkern
            r = ((ir-1.)/(nkern-1.))
            kern_tab2(ir) = fkern(sqrt(r))
        end do
        
        kernel_initialized = .true.
    end subroutine kernel_init
    
    function kern(x)
        implicit none
        real(kind=8) :: kern
        real(kind=8) :: x

        if ( x<0.5 ) then
            kern = 1-6*x**2+6*x**3
        else if ( x<=1 ) then
            kern = 2*(1-x)**3
        else
            kern = 0.
        endif
        return
    end function kern
    
    function fkern(x)
        implicit none
        
        real(kind=8) :: x

        integer :: ix
        real(kind=8) :: fw,fs
        real(kind=8) :: fkern
        
        ix = floor(x*(nkern-1)+1)
        
        if ( ix>=nkern ) then
            fkern = 0.
            return
        else
        
            fw = x-(ix-1.)/(nkern-1.)
            fs = 1./(nkern-1)
        
            fkern = kern_tab(ix)+fw/fs*(kern_tab(ix+1)-kern_tab(ix))
            if ( fkern<0 ) then
                stop '<0'
            endif
            return
        endif
        
    end function fkern
    
    function fkern2(x)
        implicit none
        
        real(kind=8) :: x

        integer :: ix
        real(kind=8) :: fw,fs
        real(kind=8) :: fkern2
        
        if ( x>=1. ) then
            fkern2 = 0.
            return
        endif
        
        ix = floor(x*(nkern-1)+1)
!         if ( ix<0 ) then
!             print *,"bad ix",ix,x,nkern
!             stop 'bad ix'
!         endif
        
        if ( ix>=nkern ) then
            fkern2 = 0.
            return
        else
        
            fw = x-(ix-1.)/(nkern-1.)
            fs = 1./(nkern-1)
        
            fkern2 = kern_tab2(ix)+fw/fs*(kern_tab2(ix+1)-kern_tab2(ix))
            if ( fkern2<0 ) then
                stop '<0'
            endif
            return
        endif
        
    end function fkern2
    

end module sph_plotter

program test
    use sph_plotter
    implicit none
    
    integer,parameter :: n=1328413
    real(kind=8),dimension(:),allocatable :: m,h,broaden,op,tau
    real(kind=8),dimension(:,:),allocatable :: r,v
    real(kind=8),dimension(:,:),allocatable :: opg
    logical,dimension(n) :: flag
    
    integer,parameter :: nray=10
    real(kind=8), dimension(nray,3) :: ray_dir,ray_offset
    
    integer,parameter :: nbins=200

    real(kind=8), dimension(2) :: c
    real(kind=8) :: w
    integer :: L
    
    integer :: i
    
    allocate(m(n))
    allocate(h(n))
    allocate(broaden(n))
    allocate(op(n))
    allocate(r(n,3))
    allocate(v(n,3))
    allocate(opg(L,L))
    allocate(tau(n))
    
    do i=1,nray
        ray_dir(i,:) = [1.,0.,1.*(i-1.)]
        ray_offset(i,:) = [0.,0.,0.]
    end do
    
    call random_number(r)
    r=r*150.d0-75.d0
    v=r ! outwards explosion
!     call random_number(v)
!     v=v*2.d0-1.d0
    

    m=1.d0
    h=10.
    broaden = 0.05
    op = 1.d-2
    flag = .true.
    
    w = 200.
    c=[-w/2,-w/2]
    L=1024
    
    call set_parallel
    opg = sph_dense(r(:,1),r(:,2),m,h,L,c,w,flag,n)
    
!     opg=sph_mean_optical_depth(r(:,1),r(:,2),r(:,3),m,h,op,L,c,w,flag,tau,n) 
    
!     rayhist =  sph_ray_histogram_opacity(r,m,h,v,-2.d0,2.d0,op,ray_dir,ray_offset,flag,broaden,nbins,nray,n)

!     open(unit=15,file="test_dumps/opacz.txt")
!     do i=1,n
!         write(15,*) r(i,:),tau(i)
!     end do
!     close(15)

    
!     open(unit=15,file="test_dumps/rayhist.txt")
!     do i=1,nbins
!         write(15,*) rayhist(i,:)
!     end do
!     close(15)

end program
