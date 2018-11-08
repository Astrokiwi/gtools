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
    integer(kind=1), parameter :: MIN_MODE = B'00000100'
    integer, parameter :: MIN_POS = 2
    integer(kind=1), parameter :: VORINOI_MODE = B'00001000'
    integer, parameter :: VORINOI_POS = 3
    integer(kind=1), parameter :: MAX_MODE = B'00010000'
    integer, parameter :: MAX_POS = 4



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

        real(kind=8), dimension(L,L) :: mg ! mass or distance grid
        
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
        
        
        if ( btest(mode,MIN_POS) ) then
            g = huge(g(1,1))
        else if ( btest(mode,MAX_POS) ) then
            g = -huge(g(1,1))
        else
            g = 0.d0
        endif
        if ( btest(mode,VORINOI_POS) ) then
            mg = huge(mg(1,1))
        else
            mg = 0.d0
        endif
        r_cell = w/L
        area_cell = r_cell**2
        
!         if (any(isnan(v)) .and. btest(mode,DENSE_WEIGHT_POS)) then
!             print *,"v",n,count(isnan(v))
!         endif
!         if (any(isnan(m))) then
!             print *,"m",n,count(isnan(v))
!         endif

        do ip=1,n
            if ( f(ip) .and. &
             (.not.btest(mode,DENSE_WEIGHT_POS) .or. &
                .not. ( isnan(v(ip)) .or. v(ip)>HUGE(v(ip)) .or. v(ip)<-HUGE(v(ip)) ) )&
                   ) then


                ix = nint((x(ip)-c(1))/r_cell)
                iy = nint((y(ip)-c(2))/r_cell)
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
        
        if ( .not. btest(mode,MAX_POS) ) then
            if ( btest(mode,MIN_POS) .or. btest(mode,VORINOI_POS) ) then
                where (g==huge(g(1,1))) g=-huge(g(1,1))
            else
                if ( btest(mode,DENSE_WEIGHT_POS) ) then
                    g = g/mg
                else
                    g = mg
                endif
            endif
        endif
        
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


    ! histogram along some radial ray
    function sph_ray_histogram(xyz,m,h,v,vmin,vmax,xyzray_in,f,broaden,nbins,nray,n) result(rayhist)
!     function sph_ray_histogram(xyz,m,h,v,vmin,vmax,xyzray_in,f,nbins,nray,n) result(rayhist)
        !$ use omp_lib
        implicit none
        
        integer :: n,nray,nbins
        
        real(kind=8), dimension(n,3) :: xyz ! particle positions
        real(kind=8), dimension(n) :: m,h ! mass, smoothing
        real(kind=8), dimension(n) :: v ! value to sum along ray
!         real(kind=8), dimension(n) :: broaden ! if <0, no broadening. If >0, width of gaussian broadening
        real(kind=8), dimension(n) :: broaden ! width of gaussian broadening

        
        real(kind=8), dimension(nray,3) :: xyzray_in ! Ray from 0,0 in these directions
        real(kind=8), dimension(nray,3) :: xyzray ! Ray from 0,0 in these directions - normalised
        
        real(kind=8) :: vmin,vmax
        
        real(kind=8), dimension(nbins,nray) :: rayhist ! ray value along los (surface density)

        real(kind=8) :: ray_norm, impact_pram, dotprod
        real(kind=8) :: h2 ! h squared to avoid square roots later
        real(kind=8) :: weight,broad_weight
        real(kind=8) :: dv2_normed, gauss_weight
        real(kind=8), dimension(n) :: gauss_norm
        
        real(kind=8) :: bin_width

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
!$OMP& private(impact_pram,ip,weight,h2,gauss_weight)&
!$OMP& shared(v,vmin,bin_width,broaden,gauss_norm)&
!$OMP& shared(rayhist,nray,nbins,xyzray,xyz,h,m,n,f)&
!$OMP& default(none) schedule(static,1)
        do ibin=1,nbins
!             print *,ibin,omp_get_thread_num()
            do ip=1,n
!             do ip=1,20000
!                 if ( mod(ip,1000)==0 .and. ibin==1) then
!                     print *,ip,"/",n
!                 endif
                if ( .not. f(ip) ) then
                    continue
                endif
                dv2_normed = ((v(ip)-vmin-(bin_width*(ibin-1)))/broaden(ip))**2
!                 if ( dv2_normed>50. ) then
                if ( dv2_normed>25. ) then
                    continue
                endif
                gauss_weight = exp(-dv2_normed)/gauss_norm(ip)
                h2 = h(ip)**2
                do iray=1,nray
                    dotprod = sum(xyzray(iray,:)*xyz(ip,:))
                    if ( dotprod<=0. ) then
                        continue
                    endif
                
                    impact_pram = norm2crossp(xyz(ip,:),xyzray(iray,:))
                    if ( impact_pram>=h2 ) then
                        continue
                    endif
                
                    impact_pram = sqrt(impact_pram)
                    weight = m(ip)*fkern(impact_pram/h(ip))/h2
                            ! should we also weight by bin width?
                    broad_weight = weight * gauss_weight
                    rayhist(ibin,iray) = rayhist(ibin,iray) + broad_weight
                end do
            end do
        end do
!$OMP END PARALLEL DO

        
!         do ip=1,n
!             if ( mod(ip,1000)==0 ) then
!                 print *,ip,"/",n
!             endif
!             if ( .not. f(ip) ) then
!                 continue
!             endif
!             h2 = h(ip)**2
! !        !$OMP PARALLEL DO private(iray,dv2_normed,broad_weight,dotprod,impact_pram)&
! !        !$OMP& shared(v,vmin,bin_width,broaden,weight,gauss_norm,rayhist,nray,nbins,ip,xyzray,xyz,h,m,h2)&
! !        !$OMP& default(none)
!             do iray=1,nray
!                 dotprod = sum(xyzray(iray,:)*xyz(ip,:))
!                 if ( dotprod<=0. ) then
!                     continue
!                 endif
!                 
!                 impact_pram = norm2crossp(xyz(ip,:),xyzray(iray,:))
!                 if ( impact_pram>=h2 ) then
!                     continue
!                 endif
!                 
! !                 if ( broaden(ip)<0. ) then ! unbroadened
!                     ibin = floor((v(ip)-vmin)/bin_width+1)
!                     if ( ibin>=1 .and. ibin<=nbins ) then
!     
!                         impact_pram = sqrt(impact_pram)
!                         weight = m(ip)*fkern(impact_pram/h(ip))/h2
!         
!                         rayhist(ibin,iray) = rayhist(ibin,iray) + weight
!                     endif
! !                 else
! !                     impact_pram = sqrt(impact_pram)
! !                     weight = m(ip)*fkern(impact_pram/h(ip))/h2
! ! !        !$OMP PARALLEL DO private(ibin,dv2_normed,broad_weight)&
! ! !        !$OMP& shared(v,vmin,bin_width,broaden,weight,gauss_norm,rayhist,iray,nbins,ip) default(none)
! !                     do ibin=1,nbins
! !                         ! should we also weight by bin width?
! !                         dv2_normed = ((v(ip)-vmin-(bin_width*(ibin-1)))/broaden(ip))**2
! !                         broad_weight = weight * exp(-dv2_normed)/gauss_norm(ip)
! !                         rayhist(ibin,iray) = rayhist(ibin,iray) + broad_weight
! !                     end do
! ! !        !$OMP END PARALLEL DO
! !                 endif
!             end do
! ! !        !$OMP END PARALLEL DO
!         end do
        
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


    function sph_optical_depth_los_old(x,y,m,h,v,op,L,c,w,z,inzarg,f,n) result(g)
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
                                !g(hix,hiy) = g(hix,hiy) + v(ip)*exp(-opg(hix,hiy)) ! is this maybe correct?
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
