function kron(A,B) result(C)
!
!   Copyright 2017, Viktor Reshiak, All rights reserved.
!
!
!   Purpose
!   =======
!   Kronecker product of two matrices
!
!
!   IN
!   ==
!   1) na1 - first  dimension of array A
!   2) na2 - second dimension of array A
!   3) A   - array
!   4) nb1 - first  dimension of array B
!   5) nb2 - second dimension of array B
!   6) B   - array
!
!
!   INOUT
!   =====
!   1) C  - resulting matrix
!
!
	implicit none


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! dummy arguments
    real(kind=8), intent(in), dimension(:,:) :: A
    real(kind=8), intent(in), dimension(:,:) :: B

    real(kind=8), dimension(size(A,1)*size(B,1),size(A,2)*size(B,2))    :: C

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer(kind=4) :: i, j
    integer(kind=4), dimension(2) :: sh_A, sh_B

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    sh_A = shape(A)
    sh_B = shape(B)

    do i = 1,sh_A(1)
        do j = 1,sh_A(2)
            C((i-1)*sh_B(1)+1:i*sh_B(1),(j-1)*sh_B(2)+1:j*sh_B(2)) = A(i,j) * B
        enddo
    enddo

end function kron




! function kron(na1,na2,A,nb1,nb2,B) result(C)
! !
! !   Copyright 2017, Viktor Reshiak, All rights reserved.
! !
! !
! !   Purpose
! !   =======
! !   Kronecker product of two matrices
! !
! !
! !   IN
! !   ==
! !   1) na1 - first  dimension of array A
! !   2) na2 - second dimension of array A
! !   3) A   - array
! !   4) nb1 - first  dimension of array B
! !   5) nb2 - second dimension of array B
! !   6) B   - array
! !
! !
! !   INOUT
! !   =====
! !   1) C  - resulting matrix
! !
! !
!     implicit none


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     ! dummy arguments

!     integer(kind=4), intent(in)                             :: na1, na2, nb1, nb2
!        real(kind=8), intent(in), dimension(na1,na2)         :: A
!        real(kind=8), intent(in), dimension(nb1,nb2)         :: B

!        real(kind=8), dimension(na1*nb1,na2*nb2)             :: C

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     integer(kind=4) :: i, j

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     do i = 1,na1
!         do j = 1,na2
!             C((i-1)*nb1+1:i*nb1,(j-1)*nb2+1:j*nb2) = A(i,j) * B
!         enddo
!     enddo

! end function kron
