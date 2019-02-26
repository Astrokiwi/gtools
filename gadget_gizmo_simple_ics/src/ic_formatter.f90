module ic_formatter
    integer, dimension(:), save, allocatable              :: int_prams
    real(kind=8), dimension(:), save, allocatable         :: real_prams
    character(len=256), dimension(:), save, allocatable   :: pram_names
    character(len=256), dimension(:), save, allocatable   :: string_prams
    logical, dimension(:), save, allocatable              :: bad_int
    
    integer, save :: n_prams
    
    contains
    
    subroutine parse_prams(infile)
        implicit none
        character(len=*) :: infile
        
        call allocate_n(infile)
        call parse_file(infile)
    end subroutine parse_prams
    
    subroutine allocate_n(infile)
        implicit none
        character(len=*) :: infile
        character(len=256) :: inline
        integer :: ierr
        n_prams = 0
        
        open(unit=12,file=infile)
        do
            read(12,*,iostat=ierr) inline
            if ( ierr/=0 ) then
                exit
            endif
            if ( len_trim(inline)<1 ) then
                exit
            endif
            n_prams = n_prams + 1
        end do
        close(12)
        
        allocate(int_prams(n_prams))
        allocate(real_prams(n_prams))
        allocate(pram_names(n_prams))
        allocate(bad_int(n_prams))
        allocate(string_prams(n_prams))
        
        print *,"nprams=",n_prams
        
    end subroutine allocate_n 
    
    subroutine parse_file(infile)
        implicit none
        character(len=*) :: infile        
        character(len=512) :: inline
        integer :: iline
        integer :: ierr
        
        open(unit=12,file=infile)
        do iline=1,n_prams
            read(12,"(A)") inline
            read(inline,*) pram_names(iline)
            string_prams(iline)=adjustl(trim(inline(len_trim(pram_names(iline))+1:)))
            read(string_prams(iline),*,iostat=ierr) int_prams(iline)
            if ( ierr/=0 ) then
                bad_int(iline)=.true.
            else
                bad_int(iline)=.false.
            endif
            read(string_prams(iline),*,iostat=ierr) real_prams(iline)
!             print *,trim(pram_names(iline)),int_prams(iline),real_prams(iline),bad_int(iline),trim(string_prams(iline))
!             print *,"X"//trim(string_prams(iline))//"X"
        end do
        close(12)

    end subroutine parse_file
    
    function get_int(pram,default_value)
        implicit none
        integer :: get_int,default_value
        character(len=*) :: pram
        
        integer :: ipram
        get_int = default_value
        
        do ipram=1,n_prams
            if ( trim(pram_names(ipram))==trim(pram) ) then
                if ( bad_int(ipram) ) then
                    print *,"BAD INT:",trim(pram_names(ipram)),real_prams(ipram)
                    return
                endif
                get_int = int_prams(ipram)
                return
            endif
        end do
!         print *,pram," not found in input "
        
        return
    end function get_int

    function get_real(pram,default_value)
        implicit none
        real(kind=8) :: get_real,default_value
        character(len=*) :: pram
        
        integer :: ipram
        get_real = default_value
        
        do ipram=1,n_prams
            if ( trim(pram_names(ipram))==trim(pram) ) then
                get_real = real_prams(ipram)
                return
            endif
        end do
!         print *,trim(pram)," not found in input "
        
        return
    end function get_real
    
    function get_string(pram,default_value)
        implicit none
        character(len=256) :: get_string,default_value
        character(len=*) :: pram
        
        integer :: ipram
        get_string = default_value
        
        do ipram=1,n_prams
            if ( trim(pram_names(ipram))==trim(pram) ) then
                get_string = string_prams(ipram)
                return
            endif
        end do
!         print *,trim(pram)," not found in input "

        return
    end function get_string
    
    subroutine tidy_up
        implicit none
        deallocate(int_prams)
        deallocate(real_prams)
        deallocate(pram_names)
        deallocate(bad_int)
        deallocate(string_prams)
    end subroutine tidy_up

end module ic_formatter
! 
! program formatter_test
!     use ic_formatter
!     implicit none
!     
!     call parse_prams("testprams.dat")
!     print *,get_real("a"),get_int("a")
!     print *,get_real("beta"),get_int("beta")
!     print *,get_real("gamma"),get_int("gamma")
!     print *,get_real("delta"),get_int("delta")
!     print *,get_real("zeta"),get_int("zeta")
!     print *,get_real("poos"),get_int("poos")
! 
! end program formatter_test
