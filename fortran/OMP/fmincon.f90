function fmincon(N,fun,X,X0,lb,ub,grad,options,info) result(fval)
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
       real(kind=8), intent(in),    dimension(N), optional :: lb, ub

    type (fsolve_options), intent(in),  optional           :: options
    type (fsolve_info),    intent(out), optional           :: info

        real(kind=8) :: fval


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! parameters of the solver
    integer(kind=4), parameter    :: M = 20, iprint = -1
    integer(kind=4)               :: iwa(3*N)
       real(kind=8)               :: wa(2*M*N+5*N+11*M*M+8*M)
    integer(kind=4), dimension(N) :: bd_type
       real(kind=8), parameter    :: factr = 0.0d0 !1.0d1 ! extremely high accuracy
       real(kind=8), parameter    :: pgtol = 0.0d0
       real(kind=8)               :: f
       real(kind=8), dimension(N) :: df
    character(len=60)             :: task, csave
    logical                       :: lsave(4)
    integer                       :: isave(44)
    real(kind=8)                  :: dsave(29)


    ! dtrnlsp_init
    type(fsolve_options) :: op



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ( present(lb) .and. present(ub) ) then
        bd_type = 2
    elseif ( present(lb) ) then
        bd_type = 1
    elseif ( present(ub) ) then
        bd_type = 3
    else
        bd_type = 0
    endif


    if ( present(options) ) then
        op = options
    else
        op = default_fsolve_options
    endif

    if ( present(X0) ) then
        X = X0
    endif

!     f  = fun(N,X)
!     df = gradient(N,X,f)
!     print*, f
!     print*, X

    task = 'START'
    do while( task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START') 
        call setulb(N,M,X,lb,ub,bd_type,f,df,factr,pgtol,wa,iwa,task,iprint,csave,lsave,isave,dsave)

        if (task(1:2).eq.'FG') then
            f  = fun(N,X)
            df = gradient(N,X,f)
        elseif (task(1:5).eq.'NEW_X') then   
            !print*, f
            !print*, X
            if ( isave(34) .ge. op%MaxIter ) then
                task = 'STOP'
            endif
            if ( dsave(13) .le. 1.d-10*(1.0d0 + abs(f)) ) then
                task = 'STOP'
            endif
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

end function fmincon