subroutine estimate_parameters_old(K,T,Y0)
    !
    !   Copyright 2017, Viktor Reshiak, All rights reserved.    
    !    
    !   Purpose
    !   =======
    !   estimate parameters of the split-step method.
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
    !
    !
    implicit none

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! dummy arguments

    integer(kind=4), intent(in)                               :: K
       real(kind=8), intent(in),  dimension(num_species)      :: Y0
       real(kind=8), intent(in),  dimension(K)                :: T

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    real(kind=8)                                         :: tau
    real(kind=8), dimension(num_reactions,num_species)   :: C
    real(kind=8), dimension(num_reactions)               :: d, a_buf
    real(kind=8), dimension(num_reactions,num_reactions) :: diag_prop
    real(kind=8), dimension(num_species,num_species)     :: eye_N
    real(kind=8), dimension(3*num_reactions)             :: params

    real(kind=8), dimension(num_species,K)             :: mu
    real(kind=8), dimension(num_species*num_species,K) :: cov
    real(kind=8), dimension(num_species,num_species)   :: P1, P3, P4
    real(kind=8), dimension(num_species)               :: p2
    real(kind=8), dimension(num_species,num_species)   :: nu_C
    real(kind=8), dimension(num_species)               :: nu_d
    real(kind=8), dimension(num_species,num_reactions) :: nu_real

    real(kind=8), dimension(num_species**2,num_species**2) :: inv_kron_P3

    type(fsolve_info) :: fsolve_inf
    integer(kind=4)   :: i, j
       real(kind=8)   :: residual, res


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tau = T(2) - T(1)
    
    eye_N = 0.d0
    forall(i=1:num_species) eye_N(i,i) = 1.d0

    nu_real = nu
    
    allocate( default_solver_options%tht(num_reactions,K-1) )
    allocate( default_solver_options%et1(num_reactions,K-1) )
    allocate( default_solver_options%et2(num_reactions,K-1) )

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! intial guess for parameters
    a_buf   = relaxRates(Y0) * tau
    params(1:num_reactions) = sqrt(2.d0/a_buf) - 1.d0/a_buf
    where ( a_buf < 2.0 ) params(1:num_reactions) = 0.633975 - 0.0566243 * a_buf
    params(num_reactions+1:3*num_reactions) = 1.0d0

    !print*, ' '
    !print*, params
    !print*, ' '

    mu(:,1)   = Y0    
    cov(:,1)  = 0.d0 
    diag_prop = 0.d0;
    do i = 2,K
        ! estimate mean
        call fsolve( num_species, num_species, fun, mu(:,i), X0=mu(:,i-1), info=fsolve_inf, Jac=Jac )

        C = propJacobian(mu(:,i))
        d = propensities(mu(:,i)) - matmul(C,mu(:,i))
        nu_C = matmul(nu,C)
        nu_d = matmul(nu,d)

        call inverse( num_species, eye_N-tau*nu_C, P1 )
        p2 = matmul( P1, nu_d )
        P3 = 0.5d0 * eye_N - tau * nu_C
        P4 = 0.5d0 * eye_N

        inv_kron_P3 = kron(num_species,num_species,eye_N,num_species,num_species,P3) + kron(num_species,num_species,P3,num_species,num_species,eye_N)
        call inverse(num_species**2,inv_kron_P3)

        ! estimate covariance
        a_buf = propensities(mu(:,i))
        forall(j=1:num_reactions) diag_prop(j,j) = a_buf(j)
        cov(:,i) = matmul( kron( num_species,num_species,eye_N,num_species,num_species,P4 ) + kron( num_species,num_species,P4,num_species,num_species,eye_N ), cov(:,i-1) ) &
                 + tau * matmul( kron(num_species,num_reactions,nu_real,num_species,num_reactions,nu_real), reshape(diag_prop,(/num_reactions**2/)) )
        call linsolve( num_species*num_species, 1, kron( num_species,num_species,eye_N,num_species,num_species,P3 ) + kron( num_species,num_species,P3,num_species,num_species,eye_N ), cov(:,i) )


        ! find optimal parameter values
        residual = fminunc( 3*num_reactions, covariance_norm, params )

        default_solver_options%tht(:,i-1) = params(1:num_reactions)
        default_solver_options%et1(:,i-1) = params(num_reactions+1:2*num_reactions)
        default_solver_options%et2(:,i-1) = params(2*num_reactions+1:3*num_reactions)

        print*, residual
    enddo
    print*, 'Final mean:'
    print*, mu(:,K)
    print*, ' '
    print*, 'Final covariance:'
    print*, cov(:,K)


    print*, ' '
    !print*, mu
    print*, default_solver_options%tht(:,1)
    print*, default_solver_options%et1(:,1)
    print*, default_solver_options%et2(:,1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    contains    
    subroutine fun(dim_f,dim_X,X,f)
        integer,      intent(in)                     :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X ) :: X
        real(kind=8), intent(out), dimension(dim_f)  :: f

        f = X - mu(:,i-1) - matmul( nu, propensities(X) * tau )
    end subroutine fun
    subroutine Jac(dim_f,dim_X,X,J)
        integer,      intent(in)                          :: dim_f, dim_X
        real(kind=8), intent(in),  dimension(dim_X)       :: X
        real(kind=8), intent(out), dimension(dim_f,dim_X) :: J

        J = eye_N - tau * matmul( nu, propJacobian(X) )
    end subroutine Jac

    function covariance_norm(dim_X,X) result(f)
        integer,      intent(in)                    :: dim_X
        real(kind=8), intent(in), dimension(dim_X ) :: X
        real(kind=8)                                :: f

        real(kind=8), dimension(num_reactions)             :: tht, et1, et2
        real(kind=8), dimension(num_species,num_reactions) :: nu2, nu3, nu4, nu5, nu6, nu7

        real(kind=8), dimension(num_species,num_species) :: R1, R3, R3_R1
        real(kind=8), dimension(num_species)             :: r2, r4

        real(kind=8), dimension(num_species**2,num_reactions**2) :: kron_nu, kron_R3_nu
        real(kind=8), dimension(num_reactions)                   :: C_R1, C_P1, C_r2, C_p2
        real(kind=8), dimension(num_reactions,num_reactions)     :: diag_R1, diag_P1, diag_r2, diag_d

        real(kind=8), dimension(num_species)                  :: dmu
        real(kind=8), dimension(num_species*num_species)  :: dcov
        real(kind=8), dimension(num_species**2,num_species**2) :: e3
        real(kind=8), dimension(num_species**2)                :: e4, e5

        integer :: j

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        tht = X(1:dim_X/3)
        et1 = X(dim_X/3+1:2*dim_X/3)
        et2 = X(2*dim_X/3+1:dim_X)
        
        forall(j=1:num_reactions) nu2(:,j) = nu(:,j) * et1(j)     * (1-tht(j))
        forall(j=1:num_reactions) nu3(:,j) = nu(:,j) * (1-et1(j)) * (1-tht(j))
        forall(j=1:num_reactions) nu4(:,j) = nu(:,j) * (1-tht(j))
        forall(j=1:num_reactions) nu5(:,j) = nu(:,j) * et2(j) * tht(j)
        forall(j=1:num_reactions) nu6(:,j) = nu(:,j) * (1-et2(j)) * tht(j)
        forall(j=1:num_reactions) nu7(:,j) = nu(:,j) * tht(j)

        R1 = eye_N + tau * matmul(nu3,C)
        r2 = matmul(nu4,d)
        R3 = eye_N + tau * matmul(nu6,C)
        r4 = matmul(nu7,d)
        call linsolve( num_species, num_species, eye_N-tau*matmul(nu2,C), R1 )
        call linsolve( num_species, 1,           eye_N-tau*matmul(nu2,C), r2 )
        call linsolve( num_species, num_species, eye_N-tau*matmul(nu5,C), R3 )
        call linsolve( num_species, 1,           eye_N-tau*matmul(nu5,C), r4 )


        nu3 = matmul(R3,nu)
        R3_R1      = matmul(R3,R1)
        kron_nu    = kron(num_species,num_reactions,nu_real,num_species,num_reactions,nu_real)
        kron_R3_nu = kron(num_species,num_reactions,nu3,num_species,num_reactions,nu3)


        C_r2 = d + matmul(C,mu(:,i))
        diag_r2 = 0.d0
        forall(j=1:num_reactions) diag_r2(j,j) = C_r2(j)


        dmu  = matmul( R3_R1 - P1, mu(:,i) ) + tau * ( matmul(R3,r2) + r4 - p2 )
        dcov = matmul( kron(num_species,num_species,R3_R1,num_species,num_species,R3_R1) - matmul(inv_kron_P3,kron(num_species,num_species,eye_N,num_species,num_species,P4) + kron(num_species,num_species,P4,num_species,num_species,eye_N)), cov(:,i))  &
             + tau * matmul( kron_R3_nu - matmul(inv_kron_P3,kron_nu), reshape(diag_r2,(/num_reactions**2/)) )

        f = sum(dmu**2) + sum(dcov**2)
    end

!     function covariance_norm(dim_X,X) result(f)
!     ! objective function for theta
!         integer,      intent(in)                    :: dim_X
!         real(kind=8), intent(in), dimension(dim_X ) :: X
!         real(kind=8)                                :: f

!         real(kind=8), dimension(num_reactions)             :: tht, et1, et2
!         real(kind=8), dimension(num_species,num_reactions) :: nu2, nu3, nu4, nu5, nu6, nu7

!         real(kind=8), dimension(num_species,num_species) :: R1, R3, R3_R1
!         real(kind=8), dimension(num_species)             :: r2, r4

!         real(kind=8), dimension(num_species**2,num_reactions**2) :: kron_nu, kron_R3_nu
!         real(kind=8), dimension(num_reactions)                   :: C_R1, C_P1, C_r2, C_p2
!         real(kind=8), dimension(num_reactions,num_reactions)     :: diag_R1, diag_P1, diag_r2, diag_d

!         real(kind=8), dimension(num_species,num_species)       :: e1
!         real(kind=8), dimension(num_species)                   :: e2
!         real(kind=8), dimension(num_species**2,num_species**2) :: e3
!         real(kind=8), dimension(num_species**2)                :: e4, e5

!         integer :: j

!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!         tht = X(1:dim_X/3)
!         et1 = X(dim_X/3+1:2*dim_X/3)
!         et2 = X(2*dim_X/3+1:dim_X)

!         !et1 = 1.0d0
!         !et2 = 1.0d0
        
!         forall(j=1:num_reactions) nu2(:,j) = nu(:,j) * et1(j)     * (1-tht(j))
!         forall(j=1:num_reactions) nu3(:,j) = nu(:,j) * (1-et1(j)) * (1-tht(j))
!         forall(j=1:num_reactions) nu4(:,j) = nu(:,j) * (1-tht(j))
!         forall(j=1:num_reactions) nu5(:,j) = nu(:,j) * et2(j) * tht(j)
!         forall(j=1:num_reactions) nu6(:,j) = nu(:,j) * (1-et2(j)) * tht(j)
!         forall(j=1:num_reactions) nu7(:,j) = nu(:,j) * tht(j)

!         R1 = eye_N + tau * matmul(nu3,C)
!         r2 = matmul(nu4,d)
!         R3 = eye_N + tau * matmul(nu6,C)
!         r4 = matmul(nu7,d)
!         call linsolve( num_species, num_species, eye_N-tau*matmul(nu2,C), R1 )
!         call linsolve( num_species, 1,           eye_N-tau*matmul(nu2,C), r2 )
!         call linsolve( num_species, num_species, eye_N-tau*matmul(nu5,C), R3 )
!         call linsolve( num_species, 1,           eye_N-tau*matmul(nu5,C), r4 )


!         nu2 = nu
!         nu3 = matmul(R3,nu)
!         R3_R1      = matmul(R3,R1)
!         kron_nu    = kron(num_species,num_reactions,nu2,num_species,num_reactions,nu2)
!         kron_R3_nu = kron(num_species,num_reactions,nu3,num_species,num_reactions,nu3)


!         C_R1 = sum(matmul(C,R1),2)
!         C_P1 = sum(matmul(C,P1),2)
!         C_r2 = d + tau*matmul(C,r2)
!         C_p2 = d + tau*matmul(C,p2)
!         diag_R1 = 0.d0
!         diag_P1 = 0.d0
!         diag_r2 = 0.d0
!         diag_d  = 0.d0
!         forall(j=1:num_reactions) diag_R1(j,j) = C_R1(j)
!         forall(j=1:num_reactions) diag_P1(j,j) = C_P1(j)
!         forall(j=1:num_reactions) diag_r2(j,j) = C_r2(j)
!         forall(j=1:num_reactions) diag_d(j,j)  = C_p2(j)


!         e1 = R3_R1 - P1
!         e2 = matmul(R3,r2) + r4 - p2
!         e3 = kron(num_species,num_species,R3_R1,num_species,num_species,R3_R1) - matmul( inv_kron_P3, ( kron(num_species,num_species,eye_N,num_species,num_species,P4) + kron(num_species,num_species,P4,num_species,num_species,eye_N) ) )       
!         e4 = matmul( kron_R3_nu, reshape(diag_R1,(/num_reactions**2/)) ) - matmul( inv_kron_P3, matmul(kron_nu,reshape(diag_P1,(/num_reactions**2/))) )
!         e5 = matmul( kron_R3_nu, reshape(diag_r2,(/num_reactions**2/)) ) - matmul( inv_kron_P3, matmul(kron_nu,reshape(diag_d,(/num_reactions**2/))) )

!         f = sum(e1**2) + sum(e2**2) + sum(e3**2) + sum(e4**2) + sum(e5**2)
!     end

end subroutine estimate_parameters_old