function [h,err,order] = StrongConvergenceTester(method,input_data,FTime,NPaths,dt_min,dt_max)
%
%   Copyright 2016, Viktor Reshiak, All rights reserved.
%
%   Purpose
%   =======
%   Test strong convergence of the numerical method
%
%
%   Method
%   ======
%   Strong error is obtained using the exact solution of the system obtained with the Gillespie algorithm (SSA)
%
%
%   IN
%   ==
%   1) method     - cell array of function handles that calls corresponding numerical methods
%   2) input_data - input data of the chemical system
%   3) FTime      - final time of the simulation
%   4) NPaths     - (optional) number of samples for Monte Carlo estimator
%   5) dt_min     - (optional) finest time discretization step
%   6) dt_max     - (optional) coarsest time discretization step
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated orders of convergence

   
    % parameters of the stochastic chemical system
    [nu,prop,Y0] = input_data();

    
    % estimate maximal total propensity and minimal time between reactions
    num_of_react = 1e4;
    A = 0;
    dt_min_r = 1;
    for i = 1:10
        [Y, poissprc, ~, num_of_react] = SSA(prop,nu,FTime,Y0,num_of_react);
        for j = 1:num_of_react
            A = max(A,sum(prop(Y(:,j))));
            dt_min_r = min(dt_min_r,min(diff(poissprc(1,:))));
        end
        fprintf('%i ', i);
    end
    dt_min = max(dt_min,dt_min_r);
    A = 1.3 * A;


    switch nargin 
        case 3
            NPaths = 100;
            dt_min = 2^(-15);
            dt_max = 2^(-3);
        case 4
            dt_min = 2^(-15);
            dt_max = 2^(-3);
    end
    dt_max = 2^floor(log2(dt_max/dt_min)) * dt_min;

    
    Nm    = length(method);                 % number of methods to test
    Nlev  = log2(dt_max/dt_min) + 1;        % number of testing levels
    
    fprintf('\n\nmethods: %s', func2str(method{1}));
    for i = 2:Nm
        fprintf(', %s', func2str(method{i}));
    end
    fprintf('\nmax total propensity A = %6.4e\n', A);
    fprintf('Recommended dt_min = %6.4e\n', dt_min_r);
    fprintf('dt_min = %6.4e\n', dt_min);
    fprintf('dt_max = %6.4e\n', dt_max);
    fprintf('Nlev   = %i\n\n', Nlev);
       
    
    % allocate various arrays
    err      = zeros(Nm,Nlev,NPaths);       % strong errors
    h        = zeros(1,Nlev);               % time steps 
    tspan    = cell(1,Nlev);                % time points
	for p = 1:Nlev
        h(p)     = 2^(p-1) * dt_min;
        tspan{p} = 0:h(p):FTime;
	end
    
    parfor path = 1:NPaths  
        % generate reference compound Poisson process
        poissprc = CompoundPoissonProcess(FTime,A,num_of_react+1);
        
        % generate reference solution
        Y_ref     = SSAd(prop,nu,FTime,Y0,i-1,poissprc);
        Y_ref_end = Y_ref(:,end);
        
        % run test levels
        for p = 1:Nlev
            % test different methods
            for n = 1:Nm
                Y = method{n}(prop,nu,tspan{p},Y0,poissprc);
                err(n, p, path) = norm( Y_ref_end - Y(:,end) );
            end
            %fprintf('Path: %i of %i    Level: %i\n',path,NPaths,p);
        end
        fprintf('Path: %i of %i \n',path,NPaths);
    end
	
    % root-mean-square error
    err = sqrt( sum(err.^2,3) ./ NPaths );
%     err = sum(err,3) ./ NPaths;

    
    order = zeros(Nm, 1);
    for n = 1:Nm
        pp = polyfit( log(h), log(err(n,:)), 1 );
        order(n) = pp(1);
    end
    
    
    
    fprintf('\n\nmethods: %s', func2str(method{1}));
    for i = 2:Nm
        fprintf(', %s', func2str(method{i}));
    end
    fprintf('\nmax total propensity A = %6.4e\n', A);
    fprintf('Recommended dt_min = %6.4e\n', dt_min_r);
    fprintf('dt_min = %6.4e\n', dt_min);
    fprintf('dt_max = %6.4e\n', dt_max);
    fprintf('Nlev   = %i\n', Nlev);
    

end