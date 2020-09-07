include "mkl_rci.f90"
include "mkl_service.f90"
include "mkl_vsl.f90"
module MKL_wrappers
	use omp_lib
	use mkl_rci
	use mkl_service		! timing 
	use mkl_vsl_type	! random numbers
	use mkl_vsl			! random numbers
	implicit none


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! nonlinear solver

	! abstract interfaces to objective functions of nonlinear solver
	abstract interface
		subroutine fsolve_objective_fun(dim_f,dim_X,X,f)
        	integer,      intent(in)                          :: dim_f, dim_X
        	real(kind=8), intent(in),  dimension(dim_X )      :: X
        	real(kind=8), intent(out), dimension(dim_f)       :: f
		end subroutine fsolve_objective_fun
		subroutine fsolve_objective_Jac(dim_f,dim_X,X,J)
			integer,      intent(in)                          :: dim_f, dim_X
			real(kind=8), intent(in),  dimension(dim_X)       :: X
			real(kind=8), intent(out), dimension(dim_f,dim_X) :: J
		end subroutine fsolve_objective_Jac
	end interface

	type fsolve_info
		integer      :: iter           = -1
		integer      :: st_cr          = -1
		real(kind=8) :: init_residual  = -1.d0
		real(kind=8) :: final_residual = -1.d0
	end type fsolve_info

	type fsolve_options
		integer(kind=4)                 :: MaxIter      = 5000
		integer(kind=4)                 :: MaxIterTrial = 100
		   real(kind=8), dimension(6)   :: eps = (/ 1.d-6, 1.d-6, 1.d-6, 1.d-6, 1.d-6, 1.d-6 /)
		   real(kind=8)                 :: rs  = 1.d0
	end type fsolve_options


	type(fsolve_options) :: default_fsolve_options


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! optimization routines

	! abstract interfaces to objective functions of optimization routine
	abstract interface
		function fminunc_objective_fun(dim_X,X) result(f)
        	integer,      intent(in)                   :: dim_X
        	real(kind=8), intent(in), dimension(dim_X) :: X
        	real(kind=8)                               :: f
		end function fminunc_objective_fun
		function fminunc_objective_grad(dim_X,X) result(g)
			integer,      intent(in)                   :: dim_X
			real(kind=8), intent(in), dimension(dim_X) :: X
			real(kind=8),             dimension(dim_X) :: g
		end function fminunc_objective_grad
	end interface

	type fminunc_info
		integer      :: iter           = -1
		integer      :: st_cr          = -1
		real(kind=8) :: init_residual  = -1.d0
		real(kind=8) :: final_residual = -1.d0
	end type fminunc_info

	type fminunc_options
		integer(kind=4)                 :: MaxIter      = 1000 !5000
		integer(kind=4)                 :: MaxIterTrial = 100
		   real(kind=8), dimension(6)   :: eps = (/ 1.d-6, 1.d-6, 1.d-6, 1.d-6, 1.d-6, 1.d-6 /)
		   real(kind=8)                 :: rs  = 1.d0
	end type fminunc_options


	type(fminunc_options) :: default_fminunc_options


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! random number generator

	integer :: num_rnd_streams = 0
    type (vsl_stream_state), dimension(:), allocatable :: stream


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains
	include "MKL_fsolve.f90"
	include "MKL_poissrnd.f90"
	!include "fsolve.f90"
	include "linspace.f90"
! 	include "praxis.f90"
! 	include "toms611.f90"
	include "fminunc.f90"
	include "fmincon.f90"
	include "gradsearch.f90"
	include "inverse.f90"
	include "linsolve.f90"
	include "kron.f90"


end module MKL_wrappers