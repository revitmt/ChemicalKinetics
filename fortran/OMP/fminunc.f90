function fminunc(N,fun,X,X0,grad,options,info,typical_val) result(residual)
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
       real(kind=8), intent(in),    dimension(N), optional :: typical_val

    type (fsolve_options), intent(in),  optional           :: options
    type (fsolve_info),    intent(out), optional           :: info

        real(kind=8) :: residual


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the sumsl/smsno solver

       real(kind=8), dimension(N)   :: d
    integer(kind=4)                 :: lv
       real(kind=8), dimension(:), allocatable  :: v
    integer(kind=4), parameter      :: liv = 60
!     integer(kind=4), parameter      :: lv  = 77 + N*(N+17)/2
    integer(kind=4), dimension(liv) :: iv
!        real(kind=8), dimension(lv)  :: v

    integer,      dimension(1) :: uip
    real(kind=8), dimension(1) :: urp



    external smsno

    ! dtrnlsp_init
    type(fsolve_options) :: op

    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    lv  = 77 + N*(N+17)/2

    allocate( v(lv) )


    if ( present(options) ) then
        op = options
    else
        op = default_fsolve_options
    endif

    if ( present(X0) ) then
        X = X0
    endif


    ! scaling vector
    if ( present(typical_val) ) then
        ! TODO: is this correct?
        d = 1.0 / typical_val
    else
        d = 1.0d0
        !d(1:N/3) = 1.d3
    endif

    ! default options
    call deflt(2, iv, liv, lv, v)
    iv(17) = 10000
    iv(18) = 10000
    iv(21) = 0
    !v(31)  = 1.0e-10
    !v(32)  = 1.0e-10
    v(35) = 1.d-6

    call smsno ( N, d, X, objective_fun, iv, liv, lv, v, uip, urp, fun )

    !print*, iv(1), v(10)

    residual = v(10)

    contains

    subroutine objective_fun(dim_X,X,nf,f,uip,urp,ufp)
        integer,      intent(in)                     :: dim_X, nf
        real(kind=8), intent(in),  dimension(dim_X ) :: X
        real(kind=8), intent(out)                    :: f
        integer,      dimension(:)                   :: uip
        real(kind=8), dimension(:)                   :: urp
        external                                     :: ufp

        f = fun(dim_X,X)
    end

end function fminunc