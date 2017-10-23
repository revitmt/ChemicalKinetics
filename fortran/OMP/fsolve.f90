subroutine fsolve(N,M,fun,X,X0,Jac,options,info)
!
!   Copyright 2016, Viktor Reshiak, All rights reserved.    
!    
!   Purpose
!   =======
!   Find solution of a system of nonlinear equations.
!
!
!   Method
!   ======
!   Intel MKL Trust-Region algorithm.
!
!
!   IN
!   ==
!   1) N       - dimension of the input vector
!   2) M       - dimension of the objective function
!   3) fun     - objective function
!   5) X0      - (optional) initial guess as a separate vector
!   6) Jac     - (optional) N-by-N Jacobian of the objective function
!   7) options - (optional) struct with options of the solver
!
!   
!   INOUT
!   =====
!   4) X       - initial guess / final solution
!
!
!   OUT
!   ===
!   8) info    - (optional) structure with information about solution process
!
	implicit none 


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! dummy arguments
	
    procedure(objective_fun)                               :: fun
    procedure(objective_Jac)                               :: Jac

    integer(kind=4), intent(in)                            :: N, M
       real(kind=8), intent(inout), dimension(N)           :: X
       real(kind=8), intent(in),    dimension(N), optional :: X0

    type (fsolve_options), intent(in),  optional           :: options
    type (fsolve_info),    intent(out), optional           :: info


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the MKL solver

    ! dtrnlsp_init
    type(fsolve_options), automatic :: op

    ! dtrnlsp_solve
    real(kind=8), dimension(N), automatic   :: fvec    
    real(kind=8), dimension(N,N), automatic :: fjac
    real(kind=8), dimension(N), automatic   :: delta1, delta2

    ! dgesv
    integer(kind=4), automatic               :: status, iter
    integer(kind=4), dimension(N), automatic :: ipiv


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	logical, automatic:: solved


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( present(options) ) then
        op = options
    else
        op = default_fsolve_options
    endif

    if ( present(X0) ) then
        X = X0
    endif

    call fun(M, N, X, fvec)
    call Jac(M, N, X, fjac)
    delta1 = -fvec

    if ( present(info) ) then
        info%init_residual = norm(fvec)
    endif

    ! run solver
    solved = .false.
    iter   = 0
	do while (.not. solved)
        iter = iter + 1

        call dgesv( N, 1, fjac, N, ipiv, delta1, N, status )
        
        X = X + delta1

        if ( (status/=0) .or. (iter>op%MaxIter) ) then
            exit
        else
            call fun(M, N, X, fvec)
            !if ( norm(delta1) <= op%eps(1)*norm(X) .and. norm(fvec) <= op%eps(2) ) then 
            if ( norm(fvec) <= op%eps(2) ) then 
                solved = .true.
            else
                delta1 = -fvec
                call Jac(M, N, X, fjac)
            endif
        endif
	enddo

    ! retrieve information
    if ( present(info) ) then
        info%st_cr = status
        info%iter  = iter
        info%final_residual = norm(fvec)
    endif

    contains

    function norm(X) result(res)
        real(kind=8), dimension(N), intent(in) :: X
        real(kind=8) :: res

        res = sqrt(sum(X**2))
    end function norm

end subroutine fsolve