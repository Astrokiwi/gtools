module fof
    contains
    
!     real(kind=8) function norm_sq(r1,r2)
!         implicit none
!         real(kind=8), dimension(3), intent(in) :: r1,r2
!         !real(kind=8), intent(out) :: norm_sq
!         
!         norm_sq = sum((r1-r2)**2)
!     end function

    function calc_fof(r,rlink,n) result(grp)
        implicit none
        integer :: n
        real(kind=8), dimension(n,3), intent(in) :: r ! particle locations
        real(kind=8), intent(in) :: rlink ! linking distance
        integer, dimension(n) :: grp ! output group numbers
        
        real(kind=8) :: rlink_sq

        
        real(kind=8) :: rdist_sq
        !real(kind=8), dimension(3) :: dr
        integer :: i,j,k
        integer :: ngrp
        
        integer, dimension(:,:,:), allocatable :: grid_ll
        integer, dimension(n) :: p_ll
        real(kind=8), dimension(3) :: grid_corner
        integer, dimension(3) :: grid_L
        real(kind=8),dimension(3) :: rmin,rmax
        integer, dimension(3) :: igrid
        integer :: ix,iy,iz

        print *,"building grid"
        
        do i=1,3
            rmin(i)=minval(r(:,i))
            rmax(i)=maxval(r(:,i))
        end do
        
        grid_corner = rmin
        grid_L = ceiling((rmax-rmin)/rlink)
        grid_L(3) = 1 ! force to be 2D, for RAM
        
        print *,grid_L
        
        allocate(grid_ll(grid_L(1),grid_L(2),grid_L(3)))
        
        grid_ll = -1
        p_ll = -1
        do i=1,n
            if ( mod(i,10000)==0 ) then
                print *,i,"/",n
            endif
            igrid = floor((r(i,:)-grid_corner)/rlink)+1
            igrid(3) = 1 ! force 2D grid
            if ( any(igrid<1) .or. any(igrid>grid_L) ) then
                print *,"out of grid"
                return
            endif
            if ( grid_ll(igrid(1),igrid(2),igrid(3))==-1 ) then
                grid_ll(igrid(1),igrid(2),igrid(3))=i
            else
                p_ll(i) = grid_ll(igrid(1),igrid(2),igrid(3))
                grid_ll(igrid(1),igrid(2),igrid(3))=i
            endif
        end do
        
        grp = -1
        
!         do ix=1,grid_L(1)
!             do iy=1,grid_L(2)
!                 if ( grid_ll(ix,iy,1)/=-1 ) then
!                     j = grid_ll(ix,iy,1)
!                     do while ( j/=-1 )
!                         grp(j) = j
!                         j=p_ll(j)
!                     end do
!                 endif
!             end do
!         end do
!         
!         do j=1,n
!             if ( grp(j)/=j ) then
!                 print *,"j wrong:",j,grp(j)
!             endif
!         end do
        
        print *,"fof steps"
        
        rlink_sq = rlink**2 ! we can compare distance**2, not distance, for efficiency, without any inaccuracy
        
        ngrp = 0
        ! naive N^2 loop
        grp = -1
        !do i=2,n
        do i=1,n
            if ( mod(i,10000)==0 ) then
                print *,i,"/",n
            endif
            igrid = floor((r(i,:)-grid_corner)/rlink)+1
            igrid(3) = 1 ! force 2D grid
!             do ix=1,grid_L(1)
!                 do iy=1,grid_L(2)
            !do ix=max(1,igrid(1)-1),min(grid_L(1),igrid(1)+1)
                !do iy=max(1,igrid(2)-1),min(grid_L(2),igrid(2)+1)
                    !do iz=max(1,igrid(3)),min(grid_L(3),igrid(3))
!                     do iz=1,1
!                     
!                         if ( grid_ll(ix,iy,iz)/=-1 ) then
!                             j = grid_ll(ix,iy,iz)
!                             do while (j/=-1)
                               do j=1,n
                                if ( j/=i ) then
                                    !rdist_sq = norm_sq(r(i,:),r(j,:))
                                    rdist_sq = sum((r(i,:)-r(j,:))**2)
                                    if ( rdist_sq<=rlink_sq ) then
                                        if ( grp(i)==-1 ) then
                                            if ( grp(j)==-1 ) then
                                                ! two linked but groupless -> create a new group
                                                grp(i) = ngrp
                                                grp(j) = ngrp
                                                ngrp = ngrp + 1
                                                !print *,"ngroups=",ngrp
                                            else
                                                ! linked to only one group -> join that group
                                                grp(i) = grp(j)
                                            endif
                                        else
                                            if ( grp(j)==-1 ) then
                                                ! linked to a group already - you can join my group
                                                grp(j) = grp(i) ! TODO something about the ordering here is wrong
                                                ! do nothing?
                                            else if ( grp(i)/=grp(j) ) then
                                                !print *,"merged",rdist_sq,rlink_sq,i,j,grp(i),grp(j)
                                                ! linked to >1 group
                                                ! merge the groups
                                                do k=1,n
                                                    if ( grp(k)==grp(j) ) then
                                                        grp(k)=grp(i)
                                                    endif
                                                end do
                                            endif
                                        endif
                                    endif
                                endif
!                                 j = p_ll(j)
                                if ( j>n ) then
                                    print *,"bad j",j
                                    return
                                endif
                            end do
!                         endif
                        
                        
!                     end do
!                 end do
!             end do
        end do
        
        deallocate(grid_ll)
        
        return
    end function

end module fof