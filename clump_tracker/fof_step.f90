module fof_step

    integer :: n

    real(kind=8), dimension(:,:), allocatable :: r
    real(kind=8) :: rlink ! linking distance
    integer, dimension(:), allocatable :: grp ! output group numbers
    
    real(kind=8) :: rlink_sq

    
    real(kind=8) :: rdist_sq

    integer :: ngrp
    
    integer, dimension(:,:,:), allocatable :: grid_ll
    integer, dimension(:), allocatable :: p_ll
    
    real(kind=8), dimension(3) :: grid_corner
    integer, dimension(3) :: grid_L
    
    integer :: ip

    integer, parameter :: max_groups = 100000
    integer, dimension(max_groups) :: grp_ll
    integer, dimension(max_groups) :: grp_n
    integer, dimension(:), allocatable :: p_grp_ll

    contains


    logical function setup(r_in,rlink_in,n_in)
        implicit none

        integer, intent(in) :: n_in
        real(kind=8), dimension(n_in,3), intent(in) :: r_in ! particle locations
        real(kind=8), intent(in) :: rlink_in ! linking distance
        
        integer :: i
        real(kind=8), dimension(3) :: rmin,rmax
        integer, dimension(3) :: igrid

        setup = .false.

        n = n_in
        
        allocate(r(n,3))
        r = r_in
        rlink = rlink_in
        
        allocate(grp(n))
        allocate(p_ll(n))
        
        do i=1,3
            rmin(i)=minval(r(:,i))
            rmax(i)=maxval(r(:,i))
        end do
        
        grid_corner = rmin
        grid_L = ceiling((rmax-rmin)/rlink)
        grid_L(3) = 1 ! force to be 2D, for RAM
        
!         print *,grid_L
        
        allocate(grid_ll(grid_L(1),grid_L(2),grid_L(3)))
        
        grid_ll = -1
        p_ll = -1
        do i=1,n
!             if ( mod(i,10000)==0 ) then
!                 print *,i,"/",n
!             endif
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
        
!         print *,"grid assignment complete"
        
        grp = -1
        rlink_sq = rlink**2 ! we can compare distance**2, not distance, for efficiency, without any inaccuracy
        
        ngrp = 0
        ip = 0

        grp_ll = -1
        grp_n = 0
        allocate(p_grp_ll(n))
        p_grp_ll = -1

        
        setup = .true.
        return
    end function setup
    
    logical function create_group(i,j)
        implicit none
        integer, intent(in) :: i,j
        
        create_group = .false.
        
        if ( ngrp==max_groups ) then
            return
        endif
        
        if ( grp(i)/=-1 .or. grp(j)/=-1 ) then
            return
        endif
        
        grp(i) = ngrp
        grp(j) = ngrp
        p_grp_ll(i) = j
        grp_n(ngrp+1) = 2
        grp_ll(ngrp+1) = i
        
        ngrp = ngrp + 1
        
        create_group = .true.
        return
    end function

    subroutine add_to_group(i,igrp)
        implicit none
        integer, intent(in) :: i,igrp
        
        grp(i) = igrp
        if ( grp_ll(igrp+1)==-1 ) then
            ! create new group
            grp_ll(igrp+1) = i
        else
            p_grp_ll(i) = grp_ll(igrp+1)
            grp_ll(igrp+1) = i
        endif
        
        grp_n(igrp+1) = grp_n(igrp+1) + 1
    end subroutine add_to_group
    
    logical function merge_groups(igrp,jgrp)
        implicit none
        integer, intent(in) :: igrp,jgrp
        
        integer :: b_grp, k_grp ! big new group, killed group
        
        integer :: j, lastj
        
        merge_groups = .false.
        
        if (    grp_n(igrp+1)<=0 .or. grp_n(jgrp+1)<=0 .or. &
                grp_ll(igrp+1)==-1 .or. grp_ll(jgrp+1)==-1 ) then
                return
        endif
        
        if ( grp_n(igrp+1)>grp_n(jgrp+1) ) then
            b_grp = igrp+1
            k_grp = jgrp+1
        else
            b_grp = jgrp+1
            k_grp = igrp+1
        endif
        
        ! convert the particles to the new group
        j = grp_ll(k_grp)
        do while ( j/=-1 )
            grp(j) = b_grp-1
            lastj = j
            j = p_grp_ll(j)
        end do
        
        ! tie the group linked lists together
        p_grp_ll(lastj) = grp_ll(b_grp)
        grp_ll(b_grp) = grp_ll(k_grp)
        
        grp_n(b_grp) = grp_n(b_grp)+grp_n(k_grp)
        
        ! kill the other group
        grp_n(k_grp) = 0
        grp_ll(k_grp) = -1
        
        merge_groups = .true.
    end function merge_groups
        

    logical function step(nsteps)
        implicit none
        
        integer,intent(in) :: nsteps

        integer :: j!,k
        integer :: istep
        integer, dimension(3) :: igrid
        
        integer :: ix,iy,iz
        
        !integer :: oldgroup
        
        step = .false.

        istep = 0
        
        if ( ip>=n ) then
            return
        endif
        
        do while ( (istep<nsteps .or. nsteps<0) .and. ip<n )
            ip = ip + 1
            istep = istep+1
            
            !print *,ip,"/",n

!             if ( mod(ip,10000)==0 ) then
!                 print *,ip,"/",n
!             endif
            igrid = floor((r(ip,:)-grid_corner)/rlink)+1
            igrid(3) = 1 ! force 2D grid
!             do ix=1,grid_L(1)
!                 do iy=1,grid_L(2)
            do ix=max(1,igrid(1)-1),min(grid_L(1),igrid(1)+1)
                do iy=max(1,igrid(2)-1),min(grid_L(2),igrid(2)+1)
!                     do iz=max(1,igrid(3)),min(grid_L(3),igrid(3))
                     do iz=1,1
!                     
                        if ( grid_ll(ix,iy,iz)/=-1 ) then
                            j = grid_ll(ix,iy,iz)
                            do while (j/=-1)
!                                do j=1,n
                                if ( j/=ip ) then
                                    !print *,j
                                    !rdist_sq = norm_sq(r(i,:),r(j,:))
                                    rdist_sq = sum((r(ip,:)-r(j,:))**2)
                                    if ( rdist_sq<=rlink_sq ) then
                                        if ( grp(ip)==-1 ) then
                                            if ( grp(j)==-1 ) then
                                                ! two linked but groupless -> create a new group
                                                
                                                if ( .not. create_group(ip,j) ) then
                                                    return
                                                endif
                                                !grp(ip) = ngrp
                                                !grp(j) = ngrp
                                                !ngrp = ngrp + 1
                                                !print *,"ngroups=",ngrp
                                            else
                                                ! linked to only one group -> join that group
                                                !grp(ip) = grp(j)
                                                call add_to_group(ip,grp(j))
                                            endif
                                        else
                                            if ( grp(j)==-1 ) then
                                                ! linked to a group already - you can join my group
                                                !grp(j) = grp(ip)
                                                call add_to_group(j,grp(ip))
                                                ! do nothing?
                                            else if ( grp(ip)/=grp(j) ) then
                                                !print *,"merged",rdist_sq,rlink_sq,i,j,grp(i),grp(j)
                                                ! linked to >1 group
                                                ! merge the groups
                                                if ( .not. merge_groups(grp(j),grp(ip)) ) then
                                                    return
                                                endif
!                                                 oldgroup = grp(j)
!                                                 do k=1,n
!                                                     if ( grp(k)==oldgroup ) then
!                                                         grp(k)=grp(ip)
!                                                     endif
!                                                 end do
                                            endif
                                        endif
                                    endif
                                endif
                                j = p_ll(j)
                                if ( j>n ) then
                                    print *,"bad j",j
                                    return
                                endif
                            end do
                        endif
                        
                        
                    end do
                end do
            end do
        end do
        
        step = .true.
        return
    
    end function step

    logical function delete_group(igrp)
        implicit none
        integer, intent(in) :: igrp
        integer :: j,oldj
        
        delete_group = .false.
        
        if ( grp_n(igrp+1)>0 ) then
            j = grp_ll(igrp+1)
            do while (j/=-1)
                grp(j) = -1
                oldj = j
                j = p_grp_ll(j)
                p_grp_ll(oldj) = -1
            end do
            grp_ll(igrp+1) = -1
            grp_n(igrp+1) = 0
        endif
        
        delete_group = .true.
    end function delete_group

    logical function group_cut(ncut)
        implicit none
        integer, intent(in) :: ncut
        integer :: igrp
        
        group_cut = .false.
        
        do igrp=1,ngrp
            if ( grp_n(igrp)<ncut .and. grp_n(igrp)>0 ) then
                if ( .not. delete_group(igrp-1) ) then
                    return
                endif
            endif
        end do
        
        group_cut = .true.
        return
    end function group_cut


    logical function group_relabel()
        implicit none
        integer :: igrp
        integer :: new_ngrp
        
        integer :: j
        
        new_ngrp = 1
        
        group_relabel = .false.
        
        do igrp=1,ngrp
            if ( grp_n(igrp)>0 ) then
                if ( igrp/=new_ngrp ) then
                    j = grp_ll(igrp)
                    do while (j/=-1)
                        grp(j) = new_ngrp-1
                        j = p_grp_ll(j)
                    end do
                    grp_ll(new_ngrp) = grp_ll(igrp)
                    grp_n(new_ngrp) = grp_n(igrp)
                
                    new_ngrp = new_ngrp + 1
                endif
            endif
        end do
        
        ngrp = new_ngrp-1
        
        group_relabel = .true.
        return
    end function group_relabel

    logical function finish()
        implicit none
        
        deallocate(r)
        deallocate(grp)
        deallocate(grid_ll)
        deallocate(p_ll)
        deallocate(p_grp_ll)
        
        finish = .true.
    end function finish


        
    function fof(r_in,rlink_in,ncut,n_in) result(grp_out)
        implicit none
        
        integer :: n_in
        real(kind=8), dimension(n_in,3) :: r_in ! particle locations
        real(kind=8) :: rlink_in ! linking distance
        integer :: ncut
        
        logical :: dummy
        
        integer, dimension(n_in) :: grp_out
        
        dummy=setup(r_in,rlink_in,n_in)
        
        dummy=step(-1)
        
        dummy=group_cut(ncut)
        
        dummy=group_relabel()
        
        grp_out = grp

        dummy=finish()
        return
    end function fof
    
    
    
end module fof_step