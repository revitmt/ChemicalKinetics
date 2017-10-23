!   Copyright 2016, Viktor Reshiak, All rights reserved.    

subroutine rnd_init(num_streams)
! 
!   Purpose
!   =======
!   Initialize random number generator of the MKL VSL library
!
    implicit none

    integer, intent(in), optional :: num_streams

    integer :: date_time(8)
    integer :: brng, seed
    integer :: i

    if ( present(num_streams) ) then 
        num_rnd_streams = num_streams
    else
        num_rnd_streams = OMP_GET_MAX_THREADS()
    endif

    allocate( stream(num_rnd_streams) )

    ! random seed
    call date_and_time(values=date_time)
    seed = date_time(8)

    ! A set of 273 Wichmann-Hill combined multiplicative congruential generators
    brng = VSL_BRNG_WH

    ! A set of 6024 Mersenne Twister pseudorandom number generators
    !brng = VSL_BRNG_MT2203

    ! intialize random number generator
    do i = 1, num_rnd_streams
        if ( vslnewstream( stream(i), brng + i - 1,  seed ) /= VSL_STATUS_OK ) then
            print*, 'Error in initialization of random streams'
            call MKL_FREE_BUFFERS
            stop 1
        endif
    enddo

end subroutine rnd_init


subroutine rnd_delete()
! 
!   Purpose
!   =======
!   Delete random streams
!
    implicit none

    integer :: i

    ! delete random streams
    ! $OMP BARRIER
    do i = 1, num_rnd_streams
        if ( vsldeletestream( stream(i) ) /= VSL_STATUS_OK ) then
            print*, 'Error in deletion of random streams'
            call MKL_FREE_BUFFERS
            stop 1
        endif
    enddo

    ! cleanup
    ! $OMP BARRIER
    if ( omp_get_thread_num() == 0 ) then
        deallocate( stream )    
    endif

end subroutine rnd_delete


subroutine rand(N,rnd_vector)
!
!   Purpose
!   =======
!   Generate vector of Uniformly(0,1) distributed random values
!
!
!   Method
!   ======
!   Intel MKL VSL library.
!
!
!   IN
!   ==
!   1) N          - dimension of the random vector
!
!
!   OUT
!   ===
!   2) rnd_vector - N-by-1 random vector 
!
!
!   Notes
!   =====
!   Random streams are assumed to be initialized. 
!   The number of streams is the number of parallel threads.
!   The current stream has the number of the calling thread
!
    implicit none 

    ! dummy arguments

    integer(kind=4), intent(in)             :: N
    real(kind=8), intent(out), dimension(N) :: rnd_vector

    integer :: method, str_num
    integer :: i

    if ( .not. allocated(stream) ) then
        if ( omp_get_thread_num() == 0 ) then
            call rnd_init()
            print*, 'Init rnd'
        endif
        !$OMP BARRIER
    endif

    ! Acceptance/rejection method
    method  = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE

    ! number of the stream
    str_num = omp_get_thread_num() + 1

    ! generate random vector
    if ( vdrnguniform( method, stream(str_num), N, rnd_vector, 0.d0, 1.d0 ) /= VSL_STATUS_OK ) then 
        print*, 'Error in generation of random numbers'
        call MKL_FREE_BUFFERS
        stop 1
    endif 
    rnd_vector = 1.d0 - rnd_vector
end subroutine rand


subroutine poissrnd(N,lambda,rnd_vector)
!
!   Purpose
!   =======
!   Generate vector of Poisson distributed random values
!
!
!   Method
!   ======
!   Intel MKL VSL library.
!
!
!   IN
!   ==
!   1) N          - dimension of the random vector
!   2) lambda     - N-by-1 vector with intensities of Poisson random variables
!
!
!   OUT
!   ===
!   3) rnd_vector - N-by-1 random vector 
!
!
!   Notes
!   =====
!   Random streams are assumed to be initialized. 
!   The number of streams is the number of parallel threads.
!   The current stream has the number of the calling thread
!
    implicit none 

    ! dummy arguments

    integer(kind=4), intent(in)                :: N
       real(kind=8), intent(in),  dimension(N) :: lambda
    integer(kind=4), intent(out), dimension(N) :: rnd_vector

    integer :: method, str_num
    integer :: i


    if ( .not. allocated(stream) ) then
        if ( omp_get_thread_num() == 0 ) then
            call rnd_init()
        endif
        !$OMP BARRIER
    endif

    ! Acceptance/rejection method
    method  = VSL_RNG_METHOD_POISSON_PTPE

    ! number of the stream
    str_num = omp_get_thread_num() + 1

    ! generate random vector
    if ( virngpoissonv( method, stream(str_num), N, rnd_vector, lambda ) /= VSL_STATUS_OK ) then 
        print*, 'Error in generation of random numbers'
        call MKL_FREE_BUFFERS
        stop 1
    endif 

    !rnd_vector = 0
end subroutine poissrnd