program statistics
	use MKL_wrappers
	use ChemicalSystem
	use methods
	implicit none


	! command line arguments
	character(len=200) :: cmd_ch, cmd_name, cmd_body
	integer            :: pos1, pos2

	! parameters of the simulation
	character(len=200) :: system
	character(len=200) :: method
	character(len=200) :: filename, filepath
	character(len=200) :: prefix
	character(len=200) :: optim_mode = 'unc'
	character(len=200) :: lin_mode   = 'all'
	character(len=200) :: stat_mode  = 'all'
	real(kind=8)       :: theta
	real(kind=8)       :: dt
	integer            :: NPaths, Nthreads
	real(kind=8)       :: Tfin
	logical            :: writeFinalTime, writeHistogram

	! results of the simulation
	real(kind=8), dimension(:),     allocatable :: T
	real(kind=8), dimension(:,:),   allocatable :: Y
	real(kind=8), dimension(:,:),   allocatable :: mu_a, mu
    real(kind=8), dimension(:,:),   allocatable :: cov_a, cov
	real(kind=8), dimension(:,:,:), allocatable :: Y_stat
	real(kind=8), dimension(:,:),   allocatable :: Y_stat_final
	type (solver_info)                          :: info
	integer,      dimension(:), allocatable     :: iters, neg_steps
	type (solver_options)                       :: options
	integer                                     :: K

	integer(kind=4), dimension(:,:), allocatable :: hist

	integer(kind=8) :: count1, count2, count_rate
	   real(kind=8) :: seconds

	integer :: i, j, nn, kk

	character(len=200) :: buf
	character(len=200) :: X0_buf, c_buf



	! default values
	system   = 'Chemical System'
	dt       = -1.d0
	theta    = -1.0
	NPaths   = 1
	Nthreads = OMP_GET_MAX_THREADS()
	prefix   = ' '
	filepath = './results/'
	writeFinalTime = .false.
	writeHistogram = .false.

	! parse command line arguments
	do i = 1,command_argument_count()
		call get_command_argument(i,cmd_ch)
		
		pos1     = scan(cmd_ch,'=')
		cmd_name = cmd_ch(1:pos1-1)
		cmd_body = cmd_ch(pos1+1:)

		select case ( cmd_name )
			case ( 'system' )
				system = trim(cmd_body)
			case ( 'method' )
				method = trim(cmd_body)
			case ( 'optim_mode' )
				optim_mode = trim(cmd_body)
			case ( 'lin_mode' )
				lin_mode = trim(cmd_body)
			case ( 'stat_mode' )
				stat_mode = trim(cmd_body)
			case ( 'file' )
				filename = trim(cmd_body)
			case ( 'prefix' )
				prefix = trim(cmd_body)
			case ( 'X0' )
				X0_buf = trim(cmd_body)
				read(cmd_body,*) X0
			case ( 'rates' )
				c_buf = trim(cmd_body)
				read(cmd_body,*) c
			case ( 'Tfin' )
				read(cmd_body,*) Tfin
			case ( 'theta' )
				read(cmd_body,*) theta
			case ( 'Nthreads' )
				read(cmd_body,*) Nthreads
			case ( 'NPaths' )
				read(cmd_body,*) NPaths
			case ( 'dt' )
				read(cmd_body,*) dt
			case ( 'writeFinalTime' )
				read(cmd_body,*) writeFinalTime
			case ( 'writeHistogram' )
				read(cmd_body,*) writeHistogram
		end select
	enddo

	if ( theta /= -1 ) then
		options%theta = theta
	endif
	if ( dt > 0 ) then 
		call linspace(0.d0, Tfin, dt, K, T)
	endif


	! compose the filename	
	if ( len(trim(filename)) == 200 ) then
		! method
		select case ( method )
			case ( 'ThetaTauLeap' )
				write(filename, '(a)'), trim(prefix)
			case ( 'ssTauLeap' )
				write(filename, '(a,a)'), trim(prefix),'_ss'
			case ( 'ssTauLeapBistable' )
				write(filename, '(a,a)'), trim(prefix),'_ss_bi'
			case ( 'ssTauLeap_old' )
				write(filename, '(a,a)'), trim(prefix),'_ss'
			case ( 'ssTauLeap2' )
				write(filename, '(a,a)'), trim(prefix),'_ss2'
			case ( 'SSA' )
				write(filename, '(a,a)'), trim(prefix),'_SSA'
		end select
		filename = trim(filename) // '_X0=' // trim(adjustl(X0_buf))
		filename = trim(filename) // '_c=' // trim(adjustl(c_buf))
		! theta
		if ( theta /= -1 ) then
			write(buf,'(f4.2)'), theta
			filename = trim(filename) // '_th=' // trim(adjustl(buf))
		endif
		! number of samlpes
		write(buf,'(es8.2)'), real(NPaths)
		filename = trim(filename) // '_N=' // trim(adjustl(buf))
		! final time
		write(buf,'(es8.2)'), Tfin
		filename = trim(filename) // '_T=' // trim(adjustl(buf))
		! time step
		if ( dt .ge. Tfin ) then
			write(buf,'(es8.2)'), Tfin
			filename = trim(filename) // '_dt=' // trim(adjustl(buf))
		else
			write(buf,'(es8.2)'), dt
			filename = trim(filename) // '_dt=' // trim(adjustl(buf))
		endif
		filename = trim(filename) // '.txt'
	endif


	! print input parameters
	print*
	print '(a)',         '                 Parameters of the simulation                         '
	print '(a)',         '========================================================================'
	print '(a)',         'Chemical system :   '//trim(system)
	write(buf,'(a3,i2,a2,i2,a1,i1,a1)'), '(a,',num_species,'f',12,'.',1,')'
	print  buf,          'Initial state   :   ',X0
	write(buf,'(a3,i2,a2,i2,a1,i1,a1)'), '(a,',num_reactions,'f',12,'.',1,')'
	print  buf,          'Reaction rates  :   ',c
	print '(a)',         'Method          :   '//trim(method)
	if ( method /= 'SSA' ) then
	print '(a, f6.3)',   'Theta           :   ', theta
	endif
	print '(a, es10.3)', 'Final time      :   ', Tfin
	if ( dt > 0 .and. method /= 'SSA' ) then
	print '(a, es10.3)', 'Time step       :   ', dt
	endif
	print '(a, i7)',     'NPaths          :   ', NPaths
	print '(a, i7)',     'OMP threads     :   ', Nthreads
	print '(a)',         'Output file     :   '//trim(filename)
	print '(a)',         '========================================================================'


	
	if ( writeFinalTime ) then
		allocate( Y_stat_final(num_species,NPaths) )
	endif
	if ( (.not. writeFinalTime) .or. writeHistogram ) then
		allocate( Y_stat(num_species,K,NPaths) )
	endif
	allocate(mu_a(num_species,K),mu(num_species,K))
	allocate(cov_a(num_species*num_species,K),cov(num_species*num_species,K))

	allocate(theta_tau(num_reactions,Nthreads))
	allocate(Y_buf(num_species,Nthreads))
	allocate(iters(NPaths), neg_steps(NPaths))

	! run solvers
	call rnd_init(Nthreads)

	! estimate parameters of the split-step method
	select case ( method )
		case ( 'ssTauLeap' )
			select case ( stat_mode )
				case ( 'all' )
					call estimate_parameters(K,T,X0,method,optim_mode,lin_mode,options,mu_a,cov_a,mu,cov)
				case ( 'end' )
					call estimate_stationary_parameters(K,T,X0,method,optim_mode,lin_mode,options,mu_a,cov_a,mu,cov)
			end select
		case ( 'ssTauLeapBistable' )
			call estimate_parameters_bistable(K,T,X0,method,optim_mode,lin_mode,options,mu_a,cov_a,mu,cov)
		case ( 'ThetaTauLeap' )
			call estimate_parameters(K,T,X0,method,optim_mode,lin_mode,options,mu_a,cov_a,mu,cov)
	end select

	! estimate execution time
	if ( NPaths > 200 ) then
		call system_clock(count1)
		!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(Y) SHARED(NPaths, stream, num_rnd_streams, nu, K, T, X0, theta_tau, Y_buf) NUM_THREADS( Nthreads )
		do i = 1,100
			select case ( method )
		 		case ( 'SSA' )
					call SSA(X0,T(K),T,Y)
	 			case ( 'ssTauLeap' )
		 			call ssTauLeap(K,T,X0,Y,options=options)
				case ( 'ssTauLeapBistable' )
		 			call ssTauLeapBistable(K,T,X0,Y,options=options)
		 		case ( 'ssTauLeap_old' )
		 			call ssTauLeap_old(K,T,X0,Y,options=options)
	 			case ( 'ThetaTauLeap' )
	 				call ThetaTauLeap(K,T,X0,Y,options=options)
			end select
		enddo
		!$OMP END PARALLEL DO
    	call system_clock(count2, count_rate)
    	seconds = (count2-count1) / real(count_rate,8) / 100 * NPaths
		print*
    	print '(a,i3,a,i2,a,i2,a)', 'Estimated execution time : ', floor(seconds/3600), ' hours ', floor(seconds/60-60*floor(seconds/3600)), ' minutes ', floor(seconds-60*floor(seconds/60)), ' seconds'
	endif


	! run solver
	call system_clock(count1)
	!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) PRIVATE(Y,info)  SHARED(Y_stat, Y_stat_final, NPaths, stream, num_rnd_streams, nu, K, T, X0, theta_tau, Y_buf, iters, neg_steps) NUM_THREADS( Nthreads )
	do i = 1,NPaths
		select case ( method )
		 	case ( 'SSA' )
				call SSA(X0,Tfin,T,Y)
	 		case ( 'ssTauLeap' )
	 			call ssTauLeap(K,T,X0,Y,options=options,info=info)
	 		case ( 'ssTauLeapBistable' )
	 			call ssTauLeapBistable(K,T,X0,Y,options=options,info=info)
	 		case ( 'ssTauLeap_old' )
	 			call ssTauLeap_old(K,T,X0,Y,options=options,info=info)
	 		case ( 'ThetaTauLeap' )
 				call ThetaTauLeap(K,T,X0,Y,options=options,info=info)
		end select
		if ( writeFinalTime ) then
			Y_stat_final(:,i) = Y(:,K)
		endif
		if ( (.not. writeFinalTime) .or. writeHistogram ) then
			Y_stat(:,:,i) = Y
		endif
		iters(i)     = info%iter
		neg_steps(i) = info%neg_steps
	enddo
	!$OMP END PARALLEL DO
    call system_clock(count2, count_rate)
    seconds  = (count2-count1) / real(count_rate,8)
	call rnd_delete()
	print '(a,i1,a,i2,a,i2,a)', 'Execution time           : ', floor(seconds/3600), ' hours ', floor(seconds/60-60*floor(seconds/3600)), ' minutes ', floor(seconds-60*floor(seconds/60)), ' seconds'


	! write data to file
	print*
	print '(a,i7,a,i8,a)', 'Writing ',NPaths,' path(s) of solutions evaluated at ',K,' time points to file "'//trim(filename)//'" ...'
	filename = trim(filepath)//trim(filename)
	open(4,file=filename)
		if ( writeFinalTime ) then
			write(4,*) num_species, num_reactions, 1, NPaths, seconds, sum(iters), sum(neg_steps)
			write(4,*) X0
			write(4,*) c
			write(4,*) nu
			write(4,*) propJacobian(X0)
			!select case ( method )
	 		!	case ( 'ssTauLeap' )
		 			write(4,*) options%tht
	 				write(4,*) options%et1
	 				write(4,*) options%et2
	 		!	case ( 'ThetaTauLeap' )
	 		!		write(4,*) options%theta
	 		!		write(4,*) options%eta1
	 		!		write(4,*) options%eta2
			!end select
			write(4,*) T(K)
			write(4,*) mu_a(:,K)
			write(4,*) mu(:,K)
			write(4,*) cov_a(:,K)
			write(4,*) cov(:,K)
			write(4,*) Y_stat_final
		else
			write(4,*) num_species, num_reactions, K, NPaths, seconds, sum(iters), sum(neg_steps)
			write(4,*) X0
			write(4,*) c
			write(4,*) nu
			write(4,*) propJacobian(X0)
			!select case ( method )
	 		!	case ( 'ssTauLeap' )
		 			write(4,*) options%tht
	 				write(4,*) options%et1
	 				write(4,*) options%et2
	 		!	case ( 'ThetaTauLeap' )
	 		!		write(4,*) options%theta
	 		!		write(4,*) options%eta1
	 		!		write(4,*) options%eta2
			!end select
			write(4,*) T
			write(4,*) mu_a
			write(4,*) mu
			write(4,*) cov_a
			write(4,*) cov
			write(4,*) Y_stat
		endif
	close(4)

	if ( writeHistogram ) then
		filename = trim(filepath)//'hist_'//trim(filename)
		print '(a)', 'Writing histogram to file "'//trim(filename)//'" ...'
		open(4,file=filename)
			write(4,*) num_species, K, NPaths
			write(4,*) T
			do nn = 1,num_species
				do kk = 1,K
					allocate( hist( 2, int(minval(Y_stat(nn,kk,:))):int(maxval(Y_stat(nn,kk,:))) ) )
					hist(1,:) = -123456789
					hist(2,:) = 0
					do i = 1,NPaths
						hist( 1, int(Y_stat(nn,kk,i)) ) = int(Y_stat(nn,kk,i))
						hist( 2, int(Y_stat(nn,kk,i)) ) = hist( 2, int(Y_stat(nn,kk,i)) ) + 1
					enddo
					!write(4,*) size(hist,2)
					!write(4,*) hist
					write(4,*) count(hist(1,:) /= (-123456789))
					write(4,*) pack(hist(1,:),hist(1,:) /= (-123456789))
					write(4,*) pack(hist(2,:),hist(1,:) /= (-123456789))
					deallocate( hist )
				enddo
			enddo
		close(4)		
	endif


	! need this for script to keep the execution window active
	!print*, 'Press any key'
	!read(*,*) 


end program statistics
