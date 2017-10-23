subroutine inverse(N,A,inv_A)
!
!   Copyright 2017, Viktor Reshiak, All rights reserved.    
!    
!
!   Purpose
!   =======
!   Find inverse of a matrix
!
!
!   Method
!   ======
!   wrapper to MKL LAPACK
!
!
!   IN
!   ==
!   1) N  - dimension of the matrix
!   2) A  - matrix
!
!
	implicit none 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! dummy arguments

    integer(kind=4), intent(in)                            :: N
       real(kind=8), intent(in),  dimension(N,N)           :: A
       real(kind=8), intent(out), dimension(N,N), optional :: inv_A

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the sumsl/smsno solver

    integer, dimension(N)   :: ipiv
    integer, dimension(2*N) :: work
    integer                 :: info

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( present(inv_A) ) then
        inv_A = A
        call dgetrf( N, N, inv_A, N, ipiv, info )
        call dgetri( N, inv_A, N, ipiv, work, 2*N, info )
    else
        call dgetrf( N, N, A, N, ipiv, info )
        call dgetri( N, A, N, ipiv, work, 2*N, info )
    endif

end subroutine inverse