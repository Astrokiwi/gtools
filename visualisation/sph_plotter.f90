module sph_plotter
    implicit none
    
    integer, parameter :: nkern = 1000
    real(kind=8),dimension(nkern) :: kern_tab
    real(kind=8), save :: kern_norm
    
    logical :: kernel_initialized = .false.
    
    integer(kind=1), parameter :: DENSE_MODE = B'00000000' ! default
    integer(kind=1), parameter :: WEIGHT_MODE = B'00000001'
    integer, parameter :: DENSE_WEIGHT_POS = 0
    integer(kind=1), parameter :: ZSLICE_MODE = B'00000010'
    integer, parameter :: ZSLICE_POS = 1



    contains
    
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
        
!         zarg = rargsort(zlos) ! SLOW! presort maybe? store values?

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
                mdepth = m(ip)*fkern(dr/h(ip))
!                 vdepth = mdepth * exp(-cum_depth)
!                 cum_depth = cum_depth + mdepth * opac(ip)

                vdepth = mdepth
            
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
    
    function sph_general(x,y,m,h,v,L,c,w,f,z,zslice,mode,n) result(g)
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

        real(kind=8), dimension(L,L) :: mg ! mass grid
        
        integer :: ip ! this particle
        integer :: ix,iy ! grid position of particle
        integer :: hix,hiy ! position in h circle
        integer :: ix0,iy0,ix1,iy1 ! bounds of h circle
        
        integer :: ih ! integer h
        
        real(kind=8) :: r_cell, area_cell
        
        real(kind=8) :: rdist, weight
        real(kind=8) :: dz
        
        if ( .not. kernel_initialized ) then
            call kernel_init
        endif
        
        g = 0.d0
        mg = 0.d0
        r_cell = w/L
        area_cell = r_cell**2
        
        do ip=1,n
            if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
                ih = ceiling(h(ip)/r_cell)
                ix = nint((x(ip)-c(1))/r_cell)
                iy = nint((y(ip)-c(2))/r_cell)
                
                ix0 = max(1,ix-ih)
                iy0 = max(1,iy-ih)

                ix1 = min(L,ix+ih)
                iy1 = min(L,iy+ih)
                
                if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
                
                    do hix=ix0,ix1
                        do hiy=iy0,iy1
                            if ( btest(mode,ZSLICE_POS) ) then
                                dz = abs(z(ip)-zslice)
                                rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell+dz**2)
                                weight = kern(rdist/h(ip))/h(ip)**3
                            else
                                rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                                weight = fkern(rdist/h(ip))/h(ip)**2
                            endif
                            
                            if ( btest(mode,DENSE_WEIGHT_POS) ) then
                                g(hix,hiy) = g(hix,hiy) + weight * m(ip) * v(ip)
                            endif
                            mg(hix,hiy) = mg(hix,hiy) + weight * m(ip)
                            
                        end do
                    end do

                endif
                
            endif
        end do
        
        if ( btest(mode,DENSE_WEIGHT_POS) ) then
            g = g/mg
        else
            g = mg
        endif
        
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


    ! Imaging in optically thick limit
    !
    !  Emissivity in each pixel is e.g. sigma*T^4 where T is the dust temperature
    !  of the tau = 1 surface. (sigma is the stefan-boltzmann constant)
    ! 
    ! Algorithm:
    !  Sort particles in z order
    !  In z-order (from near to far), each particle uses its smoothing kernel
    !   to add its contribution to optical depth (tau) to each pixel.
    !   (Optical depth is particle opacity op * column density for that pixel)
    !  Once a pixel reaches tau=1, set the corresponding output pixel to the
    !   value of v of the particle the caused tau to cross 1.0,
    !   and don't continue to add optical depth to that particle any more.
    !   (typically v is something like sigma*T_d**4)
    function sph_optical_depth_los(x,y,m,h,v,op,L,c,w,z,f,n) result(g)
        implicit none

        integer :: n,L
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! value of particle that is visible if it's at the tau=1 surface
        real(kind=8), dimension(n) :: op ! particle opacity per unit mass
        real(kind=8), dimension(n) :: x,y,z ! particle positions
        logical, dimension(n) :: f ! mask
        real(kind=8), dimension(2) :: c ! top-left corner
        real(kind=8) :: w ! width
        real(kind=8), dimension(L,L) :: g ! output grid
        
        integer, dimension(n) :: zarg ! sorted positions of particles along line of sight
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

        zarg = rargsort(z)
        do i=1,n
            ip = zarg(i)
            if ( f(ip) .and. .not. isnan(v(ip)) .and. .not. v(ip)>HUGE(v(ip)) .and. .not. v(ip)<-HUGE(v(ip)) ) then
                ih = ceiling(h(ip)/r_cell)
                ix = nint((x(ip)-c(1))/r_cell)
                iy = nint((y(ip)-c(2))/r_cell)
                
                ix0 = max(1,ix-ih)
                iy0 = max(1,iy-ih)

                ix1 = min(L,ix+ih)
                iy1 = min(L,iy+ih)
                
                if ( ix0<=L .and. iy0<=L .and. ix1>=1 .and. iy1>=1 ) then
                
                    do hix=ix0,ix1
                        do hiy=iy0,iy1
                            !if ( opg(hix,hiy)<1.d0 ) then
                                rdist = sqrt(((hix-ix)**2+(hiy-iy)**2)*area_cell)
                                weight = fkern(rdist/h(ip))/h(ip)**2
                            

                                this_opac = weight * m(ip) * op(ip)
                                g(hix,hiy) = g(hix,hiy) + min(1.,this_opac)*v(ip)*exp(-opg(hix,hiy))
                                opg(hix,hiy) = opg(hix,hiy) + this_opac
!                                 if ( opg(hix,hiy)>=1.d0 ) then
!                                     g(hix,hiy)=v(ip)
!                                 endif
                            !endif
                            
                        end do
                    end do

                endif
                
            endif
        end do

    end function

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
        
        print *,"initializing kernel"
        
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
        kern_tab(ir) = kern_tab(ir)/norm ! Normalise so it all sums to 1.
        
        ! normalise 3d kernel
        kern_norm = 0.
        do ir=1,nkern
            r = (ir-1.)/(nkern-1.)
            kern_norm = kern_norm+fkern(r)/(nkern-1.)*3.14*r**2
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
    

end module sph_plotter
