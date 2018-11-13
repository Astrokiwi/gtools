module fort_calc_neighbours
    contains
        function get_neigh(r,soft,N) result(neigh)
            implicit none
            
            integer :: N
            real(kind=8), dimension(N,3) :: r
            real(kind=8), dimension(N) :: soft
            integer, dimension(N) :: neigh
            
            integer :: i,j
            real(kind=8) :: dist2,soft2
            
            integer, parameter :: L=256
            integer, allocatable, dimension(:,:,:) :: grid_ll!,n_ll
            integer, dimension(N) :: p_ll
            integer, dimension(3) :: grid_r, loopleft, loopright
            integer :: ix,iy,iz,ih
            real(kind=8), dimension(3) :: grid_corner,grid_max
            real(kind=8) :: grid_size
            
            integer :: out_ll
            
            integer :: n_in,n_out
            
            
            
!             grid_corner = r(1,:)
!             grid_max = r(1,:)
!             do i=2,N
!                 do j=1,3
!                     if ( r(i,j)<grid_corner(j) ) then
!                         grid_corner(j) = r(i,j)
!                     endif
!                     if ( r(i,j)>grid_max(j) ) then
!                         grid_max(j) = r(i,j)
!                     endif
!                 end do
!             end do
!             grid_size = maxval(grid_max-grid_corner)
            grid_size = 1.e-2
            grid_corner = (/-grid_size/2.,-grid_size/2.,-grid_size/2./)
            
            allocate(grid_ll(L,L,L))
!             allocate(n_ll(L,L,L))
!             n_ll = 0
            grid_ll = -1
            p_ll = -1
            out_ll = -2
            n_in = 0
            n_out = 0
            do i=1,N
                grid_r = floor(L*(r(i,:)-grid_corner)/grid_size)+1
                if ( any(grid_r<1) .or. any(grid_r>L) ) then
                    if ( out_ll==-2 ) then
                        p_ll(i) = -2
                        out_ll = i
                    else
                        p_ll(i) = p_ll(out_ll)
                        out_ll = i
                    endif
!                     if ( i<100 ) then
!                         print *,"out",out_ll
!                     endif
                    n_out = n_out + 1
                else
                    if ( grid_ll(grid_r(1),grid_r(2),grid_r(3))==-1 ) then
                        grid_ll(grid_r(1),grid_r(2),grid_r(3)) = i
                    else
                        p_ll(i) = grid_ll(grid_r(1),grid_r(2),grid_r(3))
                        grid_ll(grid_r(1),grid_r(2),grid_r(3)) = i
                    endif
!                     if ( i<200 ) then
!                         print *,"in",grid_r,r(i,:),grid_corner,grid_size
!                     endif
                    n_in = n_in + 1
!                     n_ll(grid_r(1),grid_r(2),grid_r(3))=n_ll(grid_r(1),grid_r(2),grid_r(3))+1
                endif
            end do
!             print *,n_in,n_out,N
!             print *,minval(n_ll),maxval(n_ll)
!             print *,n_ll(126:130,126:130,126:130)
!             print *,count(n_ll==0)*L**(-3.)
!             deallocate(n_ll)
!             stop
            
            do i=1,N
                if ( mod(i-1,1000)==0 ) then
                    print *,i
                endif
                neigh(i) = 0
                soft2 = soft(i)**2
                grid_r = floor(L*(r(i,:)-grid_corner)/grid_size)+1
                if ( any(grid_r<1) .or. any(grid_r>L) ) then
                    ! out of grid, check all
                    do j=1,N
                        if ( i/=j ) then
                            dist2 = sum((r(i,:)-r(j,:))**2)
                            if ( dist2<=soft2 ) then
                                neigh(i)=neigh(i)+1
                            endif
                        endif
                    end do
                else
                    ih = ceiling(soft(i)*L/grid_size)
                    loopleft = grid_r-ih
                    loopright = grid_r+ih
                    do j=1,3
                        if ( loopleft(j)<1 ) then
                            loopleft(j)=1
                        endif
                        if ( loopright(j)>L ) then
                            loopright(j)=L
                        endif
                    end do
                    do ix=loopleft(1),loopright(1)
                        do iy=loopleft(2),loopright(2)
                            do iz=loopleft(3),loopright(3)
                                j = grid_ll(ix,iy,iz)
                                do while (j/=-1)
                                    if ( i/=j ) then
                                        dist2 = sum((r(i,:)-r(j,:))**2)
                                        if ( dist2<=soft2 ) then
                                            neigh(i)=neigh(i)+1
                                        endif
                                    endif
                                    j = p_ll(j)
                                end do
                            end do
                        end do
                    end do
                    ! check outside of grid p's too
                    j = out_ll
                    do while (j/=-2)
                        if ( i/=j ) then
                            dist2 = sum((r(i,:)-r(j,:))**2)
                            if ( dist2<=soft2 ) then
                                neigh(i)=neigh(i)+1
                            endif
                        endif
                        j = p_ll(j)
                    end do
                    
                endif
                
!                 print *,i,loopleft,loopright,neigh(i)
            end do
            
            deallocate(grid_ll)
            return
        end function
        
end module
