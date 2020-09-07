function gradsearch(N,fun,X,X0,grad,options,info,N_steps,MaxIter) result(fval)
!
!   Copyright 2017, Viktor Reshiak, All rights reserved.    
!    
!
!   Purpose
!   =======
!   Find minimum of the scalar function of several variables.
!
!
!   Method
!   ======
!   John Burkardt's f90 implementation of the TOMS 611 algorithm
!   ( model/trust-region approach and the double-dogleg technique )
!
!
!   IN
!   ==
!   1) N       - dimension of the input vector
!   3) fun     - objective function
!   5) X0      - (optional) initial guess as a separate vector
!   6) grad    - (optional) gradient of the objective function
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
	
    procedure(fminunc_objective_fun)                       :: fun
    procedure(fminunc_objective_grad), optional            :: grad

    integer(kind=4), intent(in)                            :: N
       real(kind=8), intent(inout), dimension(N)           :: X
       real(kind=8), intent(in),    dimension(N), optional :: X0

    integer(kind=4), intent(in), optional                  :: N_steps, MaxIter

    type (fsolve_options), intent(in),  optional           :: options
    type (fsolve_info),    intent(out), optional           :: info

        real(kind=8) :: fval


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the solver
    real(kind=8)               :: f, gamma, grad_norm
    real(kind=8), dimension(N) :: df


    ! dtrnlsp_init
    type(fsolve_options) :: op


    integer :: i, nn, max_nn


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

    if ( present(N_steps) ) then
        nn = N_steps
        max_nn = MaxIter
    else
        nn = 10000
        max_nn = 10*nn
    endif


    f  = fun(N,X)
    grad_norm = norma(N,gradient(N,X,f))
    if ( grad_norm >= 1.d-5 ) then
            gamma =  1.d-5 / grad_norm
        else
            gamma = 1.d0
    endif
    do i = 1,nn
        f  = fun(N,X)
        df = gradient(N,X,f)
        grad_norm = norma(N,df)
        do while ( fun(N,X-gamma*df) > fun(N,X) )
            gamma = gamma / 2
            nn = 2*nn - i
        enddo
        X = X - gamma * df
        if ( i > max_nn ) then
            exit
        endif
    enddo
    fval = fun(N,X)


    contains

    function gradient(dim_x,x,f_x) result(df)
        integer, intent(in)                        :: dim_x
        real(kind=8), dimension(dim_x), intent(in) :: x
        real(kind=8), intent(in)                   :: f_x
        real(kind=8), dimension(dim_x)             :: df

        real(kind=8), dimension(dim_x) :: x1
        integer      :: i
        real(kind=8) :: eps, dx

        eps = sqrt(epsilon(1.d0))

        x1 = x
        df = f_x
        do i = 1, dim_x
            dx    = max(eps*x(i),eps)
            !dx    = 1.d-10
            x1(i) = x(i) + dx
            df(i) = ( fun(dim_x,x1) - df(i) ) / dx
            x1(i) = x(i)
        enddo        
    end function gradient

    function norma(dim_x,x) result(len)
        integer, intent(in)            :: dim_x
        real(kind=8), dimension(dim_x) :: x
        real(kind=8)                   :: len

        len = sqrt(sum(x**2))
    end function norma

end function gradsearch