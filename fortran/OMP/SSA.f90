subroutine SSA(Y0,Tfin,T,Y,info)
    !
    !   Copyright 2016, Viktor Reshiak, All rights reserved.    
    !    
    !   Purpose
    !   =======
    !   Find solution of a stochastic chemical system.
    !
    !
    !   Method
    !   ======
    !   Stochastic Simulation Algorithm:
    !
    !  / y1[k+1] \   / y1[k] \   / nu_1r \
    !  | y2[k+1] |   | y2[k] |   | nu_2r |
    !  |    .    | = |   .   | + |   .   |
    !  |    .    |   |   .   |   |   .   |
    !  \ yn[k+1] /   \ yn[k] /   \ nu_nr /
    !   
    !   where r = Multinomial( 1, prop(y[k]) / Sum(prop_j) )
    !
    !
    !   IN
    !   ==
    !   1) num_species      - dimension of the state vector
    !   2) num_reactions      - number of reactions
    !   3) prop   - function to propensity vector
    !   4) nu     - num_species-by-num_reactions stoichiometric matrix
    !   5) Y0     - num_species-dimensional vector with initial data
    !   6) Tfin   - final time
    !
    !
    !   OUT
    !   ===
    !   7) T      - K-dimensional vector of time points
    !   8) Y      - num_species-by-K solution array. Each row in Y is the solution of the corresponding equation
    !   9) info   - structure with information about solvers
    !        time - average time of solver per time step
    !

    implicit none

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! dummy arguments

       real(kind=8), intent(in),    dimension(num_species)      :: Y0
       real(kind=8), intent(in)                                 :: Tfin
       real(kind=8), intent(inout), dimension(:),   allocatable :: T
       real(kind=8), intent(out),   dimension(:,:), allocatable :: Y

    type(solver_info), intent(out), optional :: info


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    real(kind=8), dimension(num_reactions)   :: a, a0

    type(solver_info) :: inf

    integer(kind=8)   :: count1, count2, count_rate
    integer(kind=4)   :: i, r, K
       real(kind=8)   :: time_old, time_new

    real(kind=8), dimension(:),   allocatable  :: tmp_arr_T
    real(kind=8), dimension(:,:), allocatable  :: tmp_arr_Y
    real(kind=8), dimension(num_species)       :: Y_new, Y_old

    logical :: all_times

    integer,      dimension(8)              :: date_time
    integer,      dimension(:), allocatable :: seed
    real(kind=8), dimension(2)              :: rnd

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! all_times -  track state of the system after each reaction or at fixed time instances
    if ( allocated(T) ) then
        K = size(T)
        if ( .not. allocated(Y) ) then
            allocate(Y(num_species,K))
        endif
        all_times = .false.
    else
        K = 10000
        allocate(T(K), Y(num_species,K))
        all_times = .true.
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!     ! random seed
!     call date_and_time(values=date_time)
!     call random_seed(size=i)
!     allocate(seed(1:i))
!     seed(:) = date_time(8) + omp_get_thread_num()
!     call random_seed(put=seed)


    time_new = 0.d0
    Y_new    = Y0
    if ( all_times ) then
        Y(:,1) = Y_new
        T(1)   = time_new
    else
        if ( T(1) .eq. 0.d0 ) then
            Y(:,1) = Y_new
        endif
    endif

    call system_clock(count1)
    i = 2
	do while ( time_new <= Tfin )
        ! dynamically change sizes of arrays
        if (all_times) then
            if (i > K) then
                call change_size_of_arrays()
            endif
        endif

        Y_old    = Y_new
        time_old = time_new

        ! calculate propensities
        a  = propensities(Y_old)
        a0 = cumsum(num_reactions,a)

        ! choose reaction
        call rand(2,rnd)
        !call random_number(rnd)
        time_new = time_old - log(rnd(1)) / a0(num_reactions)
        r        = minloc( a0, 1, a0 > (a0(num_reactions)*rnd(2)) )

        ! update state
        Y_new = Y_old + nu(:,r)
        if ( all_times ) then
            Y(:,i) = Y_new
            T(i)   = time_new
            i      = i + 1
        else
            do while ( i <= K .and. time_new > T(i) ) 
                Y(:,i) = Y_old
                i      = i + 1
            enddo
            !if ( time_new > T(i) ) then
            !    Y(:,i) = Y_old
            !    i      = i + 1
            !endif
        endif
	enddo
    !print*, i-1, K, Y(:,K)
    !print*,Tfin,T(size(T)), Y(:,size(T)), Y(:,size(T)-1), i, shape(T)
    call system_clock(count2, count_rate)
    inf%time = (count2-count1) / real(count_rate,8)


    ! remove trailing zeros
    if ( i <= K .and. all_times ) then
        K = i-1
        allocate(tmp_arr_T(K))
        allocate(tmp_arr_Y(num_species,K))
        tmp_arr_T = T(1:K)
        tmp_arr_Y = Y(:,1:K)
        deallocate(Y, T)
        call move_alloc(tmp_arr_T,T)
        call move_alloc(tmp_arr_Y,Y)
    endif

    ! interpolate solution at last time point
    if (all_times) then
        T(K)   = Tfin
    endif
    if ( T(K) > Tfin ) then
        Y(:,K) = Y(:,K-1)
    endif

    if ( present(info) ) then
        info = inf
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains

    subroutine change_size_of_arrays()
        allocate(tmp_arr_T(2*K))
        allocate(tmp_arr_Y(num_species,2*K))
        tmp_arr_T(1:K)   = T
        tmp_arr_Y(:,1:K) = Y
        deallocate(Y, T)
        call move_alloc(tmp_arr_T,T)
        call move_alloc(tmp_arr_Y,Y)
        K = 2*K
    end subroutine change_size_of_arrays


    function cumsum(dim,vector)
        implicit none
        integer,      intent(in)                 :: dim
        real(kind=8), intent(in), dimension(dim) :: vector
        real(kind=8),             dimension(dim) :: cumsum

        integer :: i

        cumsum(1) = vector(1)
        do i = 2,dim
            cumsum(i) = cumsum(i-1) + vector(i)
        enddo

    end function cumsum

end subroutine SSA