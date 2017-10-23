subroutine linsolve(n,nrhs,A,b)
!
!   Copyright 2017, Viktor Reshiak, All rights reserved.    
!    
!
!   Purpose
!   =======
!   Solve system of linear equations with Gaussian elimination
!
!
!   IN
!   ==
!   1) n    - number of linear equations
!   2) nrhs - number of right hand sides
!   2) A    - matrix
!   3) b    - right hand side
!
!   
!   INOUT
!   =====
!   3) b - solution
!
!
	implicit none 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! dummy arguments

    integer(kind=4), intent(in)                       :: n, nrhs
       real(kind=8), intent(in), dimension(n,n)       :: A
       real(kind=8), intent(inout), dimension(n,nrhs) :: b

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(kind=4), dimension(n)  :: ipiv
    integer(kind=4)                :: status

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call dgesv( n, nrhs, A, n, ipiv, b, n, status )

end subroutine linsolve