subroutine ThetaTauLeap(K,T,Y0,Y,options,info)
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
    !   Theta Tau-leaping:
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
    !   1) num_species         - dimension of the state vector
    !   2) num_reactions         - number of reactions
    !   3) K         - number of time points
    !   4) prop      - function to propensity vector, its Jacobian and rates of reactions subsystems
    !   5) nu        - num_species-by-num_reactions stoichiometric matrix
    !   6) T         - K-dimensional vector of time points
    !   7) Y0        - num_species-dimensional vector with initial data
    !   8) options   - (optional) struct with options of the solver
    !       theta    - implicitness parameter 
    !       Jacobian - true if Jacobain is supplied with 'prop' function, false otherwise
    !
    !
    !   OUT
    !   ===
    !   Y            - num_species-by-K solution array. Each row in Y is the solution of the corresponding equation
    !   info         - structure with information about solvers
    !     time       - average time of solver per time step
    !     iter       - average number of nonlinear iterations per time step
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

    type(solver_options), intent(in),  optional    :: options
    type(solver_info),    intent(out), optional    :: info


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    type(solver_options)         :: op
    type(solver_info)            :: inf
    type(fsolve_info)            :: fsolve_info

    real(kind=8)                 :: tau, theta
    real(kind=8), dimension(num_reactions) :: a_tau_s
    integer,      dimension(num_reactions) :: dN
    real(kind=8), dimension(num_reactions,num_species) :: Ja
    real(kind=8), dimension(num_species,num_species) :: eye_N

    integer(kind=8)  :: count1, count2, count_rate
    integer(kind=4)  :: i, skipped
    integer(kind=4)  :: my_thread

    
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


    tau   = T(2) - T(1)
    theta = op%theta
    
    
    eye_N = 0.d0
    forall(i=1:num_species) eye_N(i,i) = 1.d0

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    my_thread = omp_get_thread_num()+1


    inf%iter      = 0
    inf%neg_steps = 0
    inf%residual  = 0.d0
    call system_clock(count1)
    Y(:,1) = Y0
	!do i = 2,K
    i = 2
    skipped = 0
    do while (i <= K)

        ! stochastic step
        a_tau_s = propensities(Y(:,i-1))
        a_tau_s = a_tau_s * tau
        call poissrnd( num_reactions, a_tau_s, dN )
        Y(:,i) = Y(:,i-1) + matmul( nu, dN - a_tau_s * theta )

        call deterministic_corrector()
        
        ! update solution (preserve integer states)
        Y(:,i) = Y(:,i-1) + matmul( nu, nint( dN + ( propensities(Y(:,i))*tau - a_tau_s ) * theta ) );

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
            print*, 'Skip iteration', inf%iter, skipped, i, K
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

    subroutine deterministic_corrector()
        theta_tau(:,my_thread) = theta * tau
        Y_buf(:,my_thread)     = Y(:,i)

        call fsolve( num_species, num_species, fun, Y(:,i), X0=Y(:,i-1), info=fsolve_info, Jac=Jac )
        
        inf%iter     = inf%iter     + fsolve_info%iter
        inf%residual = inf%residual + fsolve_info%final_residual
    end subroutine deterministic_corrector

    subroutine preserve_positivity()
        integer(kind=4) :: j, n

        if ( any( Y(:,i) < 0 ) ) then
            inf%neg_steps = inf%neg_steps + 1
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


end subroutine ThetaTauLeap