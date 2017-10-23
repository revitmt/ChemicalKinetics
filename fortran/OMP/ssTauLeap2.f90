subroutine ssTauLeap2(K,T,Y0,Y,options,info)
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

    real(kind=8), dimension(num_species,num_species)   :: covariance, covariance_old, true_covariance
    real(kind=8), dimension(num_species)               :: mean, mean_old
    real(kind=8), dimension(num_species,num_species)   :: covariance1, covariance2, covariance3
    real(kind=8), dimension(num_species)               :: mean1, mean2, mean3
    real(kind=8), dimension(num_reactions)             :: theta_c1, theta_c2, theta_c3

    integer(kind=8)  :: count1, count2, count_rate
    integer(kind=4)  :: i, skipped
    integer(kind=4)  :: my_thread


    real(kind=8), dimension(num_species,num_species)   :: R1
    real(kind=8), dimension(num_species)               :: pt, true_mean
    integer(kind=4), dimension(num_species)            :: ipiv
    integer(kind=4)                                    :: status, j, kk
    real(kind=8) :: xT


    real ( kind = 8 ), external :: praxis

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

    mean = Y0
    true_mean = Y0
    covariance = 0.0
    pt = 0.0
    pt(1) = 1.0
    xT = sum(X0)

    mean1 = Y0
    mean2 = Y0
    mean3 = Y0
    covariance1 = 0.0
    covariance2 = 0.0
    covariance3 = 0.0

    i = 2
    skipped = 0
    Y(:,1)  = Y0
    do while (i <= K)
        Ja = propJacobian(Y(:,i-1))
        R1 = eye_N - tau * matmul( nu, Ja )

        call dgesv( num_species, 1, R1, num_species, ipiv, pt, num_species, status )
        do j = 1,num_species
            do kk = 1,num_species
                true_covariance(j,kk) = -xT * pt(j) * pt(kk)
            enddo
            true_covariance(j,j) = xT * pt(j) * ( 1 - pt(j) )
        enddo
        true_mean = xT * pt

        mean_old       = mean
        covariance_old = covariance
        call choose_theta()

!         print*, true_mean
!         print*, mean

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

    print '(a3, 3f15.2)', '1) ', mean1
    print '(a3, 3f15.2)', '2) ', mean2
    print '(a3, 3f15.2)', '3) ', mean3
    print*, ' '
    print '(a3, 4f8.4, f15.2)', '1) ', theta_c1, sqrt(residual(covariance1))
    print '(a3, 4f8.4, f15.2)', '2) ', theta_c2, sqrt(residual(covariance2))
    print '(a3, 4f8.4, f15.2)', '3) ', theta_c3, sqrt(residual(covariance3))
    
    call system_clock(count2, count_rate)
    inf%time     = (count2-count1) / real(count_rate,8)
    inf%residual = inf%residual / (K-1);


    if ( present(info) ) then
        info = inf
    endif


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains    
    
    subroutine deterministic_predictor()
        if ( op.eta1 /= 0 ) then
            theta_tau(:,my_thread) =                          theta_p *    op.eta1  * tau
            a_tau_p                = propensities(Y(:,i-1)) * theta_p * (1-op.eta1) * tau
            Y_buf(:,my_thread)     = Y(:,i-1) + matmul( nu,  a_tau_p )

            call fsolve( num_species, num_species, fun, Y(:,i), X0=Y(:,i-1), info=fsolve_inf, Jac=Jac )

            inf%iter     = inf%iter     + fsolve_inf%iter
            inf%residual = inf%residual + fsolve_inf%final_residual
            
            a_tau_s = propensities(Y(:,i)) * tau
            a_tau_p = a_tau_p + a_tau_s * theta_p * op.eta1
        else
            a_tau_p = propensities(Y(:,i-1)) * theta_p * tau
            Y(:,i)  = Y(:,i-1) + matmul( nu, a_tau_p )
            a_tau_s = propensities(Y(:,i)) * tau
        endif
    end subroutine deterministic_predictor

    subroutine deterministic_corrector()
        if ( op.eta2 /= 0 ) then
            theta_tau(:,my_thread) =                        theta_c *    op.eta2  * tau
            a_tau_c                = propensities(Y(:,i)) * theta_c * (1-op.eta2) * tau
            Y_buf(:,my_thread)     = Y(:,i) + matmul( nu,  a_tau_c )

            call fsolve( num_species, num_species, fun, Y(:,i), X0=Y(:,i-1), info=fsolve_inf, Jac=Jac )

            inf%iter     = inf%iter     + fsolve_inf%iter
            inf%residual = inf%residual + fsolve_inf%final_residual

            a_tau_c = a_tau_c + propensities(Y(:,i)) * theta_c * op.eta2 * tau
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


    subroutine choose_theta()
        real(kind=8) :: residual, temp(2)


        mean_old       = mean1
        covariance_old = covariance1
        if ( op.theta == -1.0 ) then
            a_buf   = relaxRates(Y(:,i-1)) * tau
            theta_c = sqrt(2.d0/a_buf) - 1.d0/a_buf
            where ( a_buf < 2.0 ) theta_c = 0.633975 - 0.0566243 * a_buf
            theta_p = 1.d0 - theta_c
        endif
        call update_mean_covariance(theta_c)
        mean1       = mean
        covariance1 = covariance
        theta_c1    = theta_c

! print '(a,4e15.4)', '1)', theta_c
        mean_old       = mean2
        covariance_old = covariance2
        temp = 0.05
        call fsolve( 1, 1, covariance_norm, temp, info=fsolve_inf)
        theta_c = temp(1)
        mean2       = mean
        covariance2 = covariance
        theta_c2    = theta_c
! print '(a,4e15.4)', '2)', temp(1)
!         if (i==2) then
!             temp = 0.05
!             theta_c = temp(1)
!             temp(1) = theta_c(1)
!             temp(2) = theta_c(3)
!         endif
        mean_old       = mean3
        covariance_old = covariance3
        temp(1) = theta_c1(1)
        temp(2) = theta_c1(3)
        residual = fminunc( 2, covariance_norm_toms, temp )
        theta_c(1) = temp(1)
        theta_c(2) = temp(1)
        theta_c(3) = temp(2)
        theta_c(4) = temp(2)
        mean3       = mean
        covariance3 = covariance
        theta_c3    = theta_c
!         theta_c = (/ 0.036, 0.095, 0.036, 0.095/)
!         where(theta_c>1) theta_c = 0.99
!         where(theta_c<0) theta_c = 0.01
        theta_p  = 1.d0 - theta_c
! print '(a,4e15.4,e15.4)', '3)', theta_c, residual
!         print '(4e10.2)', theta_c

!         if (residual > 1.0 .or. residual < 1.e-9) then
!         print '(i4,e10.2,4e10.2)', i, residual, theta_c
!         print '(i4,e10.2,4e10.2)', i, residual, theta_p
!         endif

    end subroutine choose_theta  

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function covariance_norm_toms(dim_X,X) result(f)
    ! objective function for theta
        integer,      intent(in)                    :: dim_X
        real(kind=8), intent(in), dimension(dim_X ) :: X
        real(kind=8)                                :: f
        real(kind=8)                                :: f2

        real(kind=8), dimension(num_reactions)           :: theta
        real(kind=8), dimension(num_species,num_species) :: dcov
        real(kind=8), dimension(num_species)             :: dmean
        integer(kind=4)                                  :: j, k


        theta(1) = X(1)
        theta(2) = X(1)
        theta(3) = X(2)
        theta(4) = X(2)
        call update_mean_covariance(theta)
        dcov  = covariance - true_covariance
        dmean = mean       - true_mean

        f  = 0.0
        f2 = 0.0
        do j = 1,num_species
            do k = 1,num_species
                f = f + dcov(j,k)**2
            enddo
            f2 = f2 + dmean(j)**2
        enddo
        !f = f + 1000*f2
        !f = f + f2
    end


    subroutine covariance_norm(dim_f,dim_X,X,f)
    ! objective function for theta
        integer,      intent(in)                     :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X ) :: X
        real(kind=8), intent(out), dimension(dim_f)  :: f

        real(kind=8), dimension(num_reactions)           :: theta
        real(kind=8), dimension(num_species,num_species) :: dcov
        integer(kind=4)                                  :: j, k

        theta = X(1)
        call update_mean_covariance(theta)

        dcov = covariance - true_covariance

        f(1) = 0.0
        do j = 1,num_species
            do k = 1,num_species
                f(1) = f(1) + dcov(j,k)**2
            enddo
        enddo
    end


    subroutine update_mean_covariance(theta)
        real(kind=8), intent(in), dimension(num_reactions) :: theta

        real(kind=8), dimension(num_species)               :: temp1
        real(kind=8), dimension(num_reactions)             :: temp2
        real(kind=8), dimension(num_species,num_reactions) :: temp3, temp4, nu2, nu3
        real(kind=8), dimension(num_species,num_species)   :: R1, R3, A
        integer                                            :: j


        forall(j=1:num_reactions) nu2(:,j) = nu(:,j) * (1.0-theta(j))
        forall(j=1:num_reactions) nu3(:,j) = nu(:,j) * theta(j)

        R1 = eye_N - op.eta1 * tau * matmul( nu2, Ja )
        R3 = eye_N - op.eta2 * tau * matmul( nu3, Ja )
        call inverse(num_species,R1)
        call inverse(num_species,R3)

        ! propagation matrix
        A = matmul(R3,R1)

        temp1 = matmul(R1,mean_old)
        temp2 = matmul(Ja,temp1)
        temp3 = matmul(R3,nu)
        forall(j=1:num_reactions) temp4(:,j) = temp3(:,j) * temp2(j)
        R1 = tau * matmul(temp4,transpose(temp3))

        ! update mean
        mean = matmul(A,mean_old)

        ! update covariance
        covariance = matmul(A,covariance_old)
        covariance = matmul(covariance,transpose(A))
        covariance = covariance + R1
    end


    function residual(cov) result(f)
        real(kind=8), intent(in), dimension(num_species,num_species) :: cov
        real(kind=8) :: f

        integer :: j,k
        real(kind=8), dimension(num_species,num_species) :: dcov

        dcov = cov - true_covariance

        f  = 0.0
        do j = 1,num_species
            do k = 1,num_species
                f = f + dcov(j,k)**2
            enddo
        enddo
    end


    subroutine fun(dim_f,dim_X,X,f)
    ! objective function for deterministic predictor & corrector
        integer,      intent(in)                     :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X ) :: X
        real(kind=8), intent(out), dimension(dim_f)  :: f
        integer :: my_thread

        my_thread = omp_get_thread_num()+1

        f = X - Y_buf(:,my_thread) - matmul( nu, propensities(X) * theta_tau(:,my_thread) )
    end subroutine fun

    subroutine Jac(dim_f,dim_X,X,J)
    ! objective Jacobian for deterministic predictor & corrector
        integer,      intent(in)                          :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X)       :: X
        real(kind=8), intent(out), dimension(dim_f,dim_X) :: J
        integer :: my_thread

        my_thread = omp_get_thread_num()+1

        Ja = propJacobian(X)

        forall(i=1:num_species) Ja(:,i) = Ja(:,i) * theta_tau(:,my_thread)

        J = eye_N - matmul( nu, Ja )
    end subroutine Jac


end subroutine ssTauLeap2