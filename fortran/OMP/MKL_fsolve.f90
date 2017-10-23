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
	
    procedure(fsolve_objective_fun)                        :: fun
    procedure(fsolve_objective_Jac), optional              :: Jac

    integer(kind=4), intent(in)                            :: N, M
       real(kind=8), intent(inout), dimension(N)           :: X
       real(kind=8), intent(in),    dimension(N), optional :: X0

    type (fsolve_options), intent(in),  optional           :: options
    type (fsolve_info),    intent(out), optional           :: info


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the MKL solver

    ! dtrnlsp_init
    type(fsolve_options) :: op
    type(HANDLE_TR)      :: handle

    ! dtrnlsp_solve
       real(kind=8), dimension(N)   :: fvec    
       real(kind=8), dimension(N,N) :: fjac
    integer(kind=4)                 :: RCI_Request


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	logical :: solved


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if ( present(options) ) then
        op = options
    else
        op = default_fsolve_options
    endif

    RCI_Request = 0
    solved      = .false.


    if ( present(X0) ) then
        X = X0
    endif

    call fun(M, N, X, fvec)
    if ( present(Jac) ) then
        call Jac(M, N, X, fjac)
    else
        if ( djacobi(fun,N,M,fjac,X,op%eps(1)) /= TR_SUCCESS ) then
            print*, '| error in djacobix'
            call MKL_FREE_BUFFERS
            stop 1;
        end if
    endif


    ! initialize solver
	if ( dtrnlsp_init(handle, N, M, X, op%eps, op%MaxIter, op%MaxIterTrial, op%rs) /= TR_SUCCESS ) then 
        print*, '| error in dtrnlsp_init'
        call MKL_FREE_BUFFERS
        stop 1;
    endif
    !tmp = dtrnlsp_check(handle, n, m, fjac, fvec, op%eps, inf1)
    

    ! run solver
	do while (.not. solved)
        if (dtrnlsp_solve (handle, fvec, fjac, RCI_Request) /= TR_SUCCESS) then
            print*, '| error in dtrnlsp_solve'
            call MKL_FREE_BUFFERS
            stop 1;
        end if

        if (RCI_Request < 0) then
            solved = .true.
        end if

        ! recalculate function value
        if (RCI_Request == 1) then
            call fun(M, N, X, fvec)
        end if

        ! recalculate Jacobian
        if (RCI_Request == 2) then
            if ( present(Jac) ) then
        	   call Jac(M, N, X, fjac)
            else
                if ( djacobi(fun,N,M,fjac,X,op%eps(1)) /= TR_SUCCESS ) then
                    print*, '| error in djacobix'
                    call MKL_FREE_BUFFERS
                    stop 1;
                end if
            endif
        end if
	enddo

    ! retrieve information
    if ( present(info) ) then
        if (dtrnlsp_get(handle, info%iter, info%st_cr, info%init_residual, info%final_residual) /= TR_SUCCESS) then
            print*, '| error in dtrnlsp_get'
            call MKL_FREE_BUFFERS
            stop 1;
        end if
        !if ( info%iter > 10 ) then
        !    print*, info%st_cr, info%iter, info%init_residual, info%final_residual
        !endif
    endif

    if (dtrnlsp_delete(handle) /= TR_SUCCESS) then
        print*, '| error in dtrnlsp_delete'
        call MKL_FREE_BUFFERS
        stop 1;
    end if

end subroutine fsolve