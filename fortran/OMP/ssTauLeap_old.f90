subroutine ssTauLeap_old(K,T,Y0,Y,options,info)
    !
    !   Copyright 2016, Viktor Reshiak, All rights reserved.    
    !    
    !   Purpose
    !   =======
    !   Find solution of a_tau_s stochastic chemical system.
    !
    !
    !   Method
    !   ======
    !   Two-stage split-step Tau-leaping:
    !                                                          /                                                   \
    !  / y1'[k] \    / y1[k] \    / nu_11  nu_12 ... nu_1m \   |        / prop_1'[k] \               / prop_1[k] \ | 
    !  | y2'[k] |    | y2[k] |    | nu_21  nu_22 ... nu_2m |   |        | prop_2'[k] |               | prop_2[k] | |
    !  |    .   |  = |   .   |  + |   .     .          .   | * | eta1 * |     .      |  + (1-eta1) * |   .       | | * (1-theta) * tau
    !  |    .   |    |   .   |    |   .           .    .   |   |        |     .      |               |   .       | |
    !  \ yn'[k] /    \ yn[k] /    \ nu_n1  nu_n2 ... nu_nm /   |        \ prop_n'[k] /               \ prop_n[k] / |  
    !                                                          \                                                   /
    !
    !                                                          /                                 \
    !  / y1''[k] \   / y1'[k] \   / nu_11  nu_12 ... nu_1m \   | / P_1' \   / prop_1'[k] \       |
    !  | y2''[k] |   | y2'[k] |   | nu_21  nu_22 ... nu_2m |   | | P_2' |   | prop_2'[k] |       |
    !  |    .    | = |   .    | + |   .     .          .   | * | |  .   | - |     .      | * tau |
    !  |    .    |   |   .    |   |   .           .    .   |   | |  .   |   |     .      |       |
    !  \ yn''[k] /   \ yn'[k] /   \ nu_n1  nu_n2 ... nu_nm /   | \ P_m' /   \ prop_n'[k] /       |
    !                                                          \                                 /
    !   where P_j' = Poisson( prop_j(y'[k]) * tau )
    !                                                           /                                                      \
    !  / y1[k+1] \   / y1''[k] \   / nu_11  nu_12 ... nu_1m \   |        / prop_1[k+1] \               / prop_1''[k] \ | 
    !  | y2[k+1] |   | y2''[k] |   | nu_21  nu_22 ... nu_2m |   |        | prop_2[k+1] |               | prop_2''[k] | |
    !  |   .     | = |   .     | + |   .     .          .   | * | eta2 * |     .       |  + (1-eta2) * |     .       | | * theta * tau
    !  |   .     |   |   .     |   |   .           .    .   |   |        |     .       |               |     .       | |
    !  \ yn[k+1] /   \ yn''[k] /   \ nu_n1  nu_n2 ... nu_nm /   |        \ prop_n[k+1] /               \ prop_n''[k] / |  
    !                                                           \                                                      /
    !
    !
    !   IN
    !   ==
    !   K         - number of time points
    !   T         - K-dimensional row vector of time points
    !   Y0        - num_species-dimensional column vector with initial data
    !   options   - (optional) struct with options of the solver
    !    theta    - implicitness parameter 
    !    eta1     - implicitness parameter of the ODE solver for the deterministic predictor 
    !    eta2     - implicitness parameter of the ODE solver for the deterministic corrector 
    !    Jacobian - true if Jacobain is supplied with 'prop' function, false otherwise
    !
    !
    !   OUT
    !   ===
    !   Y            - num_species-by-K solution array. Each row in Y is the solution of the corresponding equation
    !   info         - structure with information about solvers
    !     time       - average time of solver per time step
    !     iter       - average number of nonlinear iterations per time step
    !   dN           - num_reactions-by-K array of Poisson increments (reaction counts)
    !   poissmsr_tau - 2-by-K2 array of driving compound Poisson measure { (t_i,a_i), i = 1,...,K2 }.
    !
    !
    !   Notes
    !   =====
    !   It is assumed that time grid is uniform
    
    implicit none

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! dummy arguments

    integer(kind=4), intent(in)                               :: K
       real(kind=8), intent(in),  dimension(num_species)      :: Y0
       real(kind=8), intent(in),  dimension(K)                :: T
       real(kind=8), intent(out), dimension(:,:), allocatable :: Y

    type(solver_options), intent(in),  optional :: options
    type(solver_info),    intent(out), optional :: info


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    type(solver_options) :: op
    type(solver_info)    :: inf
    type(fsolve_info)    :: fsolve_inf

    real(kind=8)                                       :: tau
    real(kind=8), dimension(num_reactions)             :: theta_c, theta_p
    real(kind=8), dimension(num_reactions)             :: a_tau_p, a_tau_c, a_tau_s, a_buf
    integer,      dimension(num_reactions)             :: dN
    real(kind=8), dimension(num_reactions,num_species) :: Ja
    real(kind=8), dimension(num_species,num_species)   :: eye_N

    integer(kind=8)  :: count1, count2, count_rate
    integer(kind=4)  :: i, skipped
    integer(kind=4)  :: my_thread


    real(kind=8), dimension(num_reactions)             :: et1, et2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    if ( present(options) ) then
        op = options
    else
        op = default_solver_options
    endif


    if ( .not. allocated(Y) ) then
        allocate(Y(num_species,K))
    endif


    tau     = T(2) - T(1)
    theta_c = (/ (op%theta,i=1,num_reactions) /)
    theta_p = 1 - theta_c
    
    
    eye_N = 0.d0
    forall(i=1:num_species) eye_N(i,i) = 1.d0

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    my_thread = omp_get_thread_num()+1

    inf%iter     = 0
    inf%residual = 0.d0
    call system_clock(count1)

    i = 2
    skipped = 0
    Y(:,1)  = Y0
    do while (i <= K)
        call choose_theta()

        call deterministic_predictor()        

        ! stochastic step
        call poissrnd( num_reactions, a_tau_s, dN )
        Y(:,i) = Y(:,i) + matmul( nu, dN - a_tau_s )

        call deterministic_corrector()

        ! update solution (preserve integer states)
        Y(:,i) = Y(:,i-1) + matmul( nu, nint( dN + a_tau_p + a_tau_c - a_tau_s ) )

        call preserve_positivity()

        i = i + 1
        if ( inf%iter == default_fsolve_options%MaxIter ) then
            i = i - 1;
            skipped = skipped + 1
            if ( skipped > K * 0.5 ) then
                print*, 'Solution did not converge!', dN, Y(:,i)
                ! need this for script to keep the execution window active
                print*, 'Press any key'
                read(*,*) 
                call exit(0)
            endif
        endif
	enddo
    
    call system_clock(count2, count_rate)
    inf%time     = (count2-count1) / real(count_rate,8)
    inf%residual = inf%residual / (K-1);
    !inf%iter     = inf%iter     / (K-1);


    if ( present(info) ) then
        info = inf
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains    
    
    subroutine choose_theta()
        if ( op.theta == -1.0 ) then
            a_buf   = relaxRates(Y(:,i-1)) * tau
            theta_c = sqrt(2.d0/a_buf) - 1.d0/a_buf
            where ( a_buf < 2.0 ) theta_c = 0.633975 - 0.0566243 * a_buf
            theta_p = 1.d0 - theta_c
        endif

        et1 = 1.0
        et2 = 1.0

        ! choose theta based on a particular reaction
        if ( op.theta > 1 ) then
            a_buf   = relaxRates(Y(:,i-1))
            a_buf   = a_buf(int(op.theta)) * tau
            theta_c = sqrt(2.d0/a_buf) - 1.d0/a_buf
            where ( a_buf < 2.0 ) theta_c = 0.633975 - 0.0566243 * a_buf
            theta_c = maxval(theta_c)
            theta_p = 1.d0 - theta_c
        endif
    end subroutine choose_theta        

    subroutine deterministic_predictor()
        if ( op.eta1 /= 0 ) then
            theta_tau(:,my_thread) =                          theta_p *    et1  * tau
            a_tau_p                = propensities(Y(:,i-1)) * theta_p * (1-et1) * tau
            Y_buf(:,my_thread)     = Y(:,i-1) + matmul( nu,  a_tau_p )

            call fsolve( num_species, num_species, fun, Y(:,i), X0=Y(:,i-1), info=fsolve_inf, Jac=Jac )

            inf%iter     = inf%iter     + fsolve_inf%iter
            inf%residual = inf%residual + fsolve_inf%final_residual
            
            a_tau_s = propensities(Y(:,i)) * tau
            a_tau_p = a_tau_p + a_tau_s * theta_p * et1
        else
            a_tau_p = propensities(Y(:,i-1)) * theta_p * tau
            Y(:,i)  = Y(:,i-1) + matmul( nu, a_tau_p )
            a_tau_s = propensities(Y(:,i)) * tau
        endif
    end subroutine deterministic_predictor

    subroutine deterministic_corrector()
        if ( op.eta2 /= 0 ) then
            theta_tau(:,my_thread) =                        theta_c *    et2  * tau
            a_tau_c                = propensities(Y(:,i)) * theta_c * (1-et2) * tau
            Y_buf(:,my_thread)     = Y(:,i) + matmul( nu,  a_tau_c )

            call fsolve( num_species, num_species, fun, Y(:,i), X0=Y(:,i-1), info=fsolve_inf, Jac=Jac )

            inf%iter     = inf%iter     + fsolve_inf%iter
            inf%residual = inf%residual + fsolve_inf%final_residual

            a_tau_c = a_tau_c + propensities(Y(:,i)) * theta_c * et2 * tau
        else
            a_tau_c = propensities(Y(:,i)) * theta_c * tau
            Y(:,i)  = Y(:,i) + matmul( nu,  a_tau_c )
        endif
    end subroutine deterministic_corrector

    subroutine preserve_positivity()
        integer(kind=4) :: j, n

        if ( any( Y(:,i) < 0 ) ) then
            do j = 1,num_species
                do while ( Y(j,i)<0 )
                    do n = 1,num_reactions
                        if ( nu(j,n)<0 ) then
                            Y(:,i) = Y(:,i) - nu(:,n)
                        endif
                    enddo
                enddo
            enddo
        endif
    end subroutine preserve_positivity

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! objective function for deterministic predictor & corrector
    subroutine fun(dim_f,dim_X,X,f)
        integer,      intent(in)                     :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X ) :: X
        real(kind=8), intent(out), dimension(dim_f)  :: f
        integer :: my_thread

        my_thread = omp_get_thread_num()+1

        f = X - Y_buf(:,my_thread) - matmul( nu, propensities(X) * theta_tau(:,my_thread) )
    end subroutine fun

    ! objective Jacobian for deterministic predictor & corrector
    subroutine Jac(dim_f,dim_X,X,J)
        integer,      intent(in)                          :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X)       :: X
        real(kind=8), intent(out), dimension(dim_f,dim_X) :: J
        integer :: my_thread

        my_thread = omp_get_thread_num()+1

        Ja = propJacobian(X)

        forall(i=1:num_species) Ja(:,i) = Ja(:,i) * theta_tau(:,my_thread)

        J = eye_N - matmul( nu, Ja )
    end subroutine Jac


end subroutine ssTauLeap_old