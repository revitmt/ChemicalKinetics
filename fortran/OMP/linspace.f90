subroutine linspace(x1, x2, dx, n, y)
	implicit none

	real(kind=8), intent(in)                :: x1, x2
	real(kind=8), intent(in)                :: dx
	real(kind=8), intent(out), dimension(:), allocatable :: y
	integer :: n

	integer :: i
	

	n = ceiling( ( x2 - x1 ) / dx ) + 1

	allocate(y(n))
		
	y    = (/ ( (x1+i*dx), i=0,n-1 ) /)
	y(n) = min(y(n),x2)

end subroutine linspace


! function dtime()
! 	implicit none

! 	real(kind=8) :: dtime

! 	integer(kind=8) :: clocks

! 	call mkl_get_cpu_clocks( clocks )

! ! 	dtime = mkl_get_clocks_frequency() / clocks
! 	dtime = clocks / mkl_get_clocks_frequency()

! end function dtime