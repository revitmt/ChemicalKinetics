module methods
	use ChemicalSystem
	use MKL_wrappers
	implicit none


	type solver_options
		real(kind=4) :: theta    = -1.0
	    real(kind=4) :: eta1     =  1.0
	    real(kind=4) :: eta2     =  1.0
	    logical      :: Jacobian = .true.
	    real(kind=8), dimension(:,:),   allocatable :: tht, et1, et2
	    real(kind=8), dimension(:,:,:), allocatable :: tht_bi, et1_bi, et2_bi
	    real(kind=8), dimension(:),     allocatable :: stab_loc
	    real(kind=8), dimension(:),     allocatable :: stab_width
	end type solver_options

	type solver_info 
		real(kind=8)    :: time		 = -1.0
		real(kind=8)    :: residual	 = -1.0
		integer(kind=4) :: iter 	 = -1
		integer(kind=4) :: neg_steps = -1
	end type solver_info


	type (solver_options) :: default_solver_options


	! arrays with values local to each thread
    real(kind=8), dimension(:,:), allocatable :: theta_tau
    real(kind=8), dimension(:,:), allocatable :: Y_buf

contains
	include "ssTauLeap.f90"
	include "ssTauLeap_old.f90"
	include "ssTauLeapBistable.f90"
	include "ThetaTauLeap.f90"
	include "SSA.f90"
	include "estimate_parameters.f90"
	include "estimate_stationary_parameters.f90"
	include "estimate_parameters_bistable.f90"
	!include "estimate_parameters_old.f90"

end module methods
