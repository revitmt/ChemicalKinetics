function [h,err,order,tspan,Y_mean,Y_var,Y,SSA_tspan,SSA_mean,SSA_var,SSA_Y] = WeakConvergenceTester(methods,input_data,FTime,h,options)
%
%   Copyright 2016, Viktor Reshiak, All rights reserved.
%
%   Purpose
%   =======
%   Test weak convergence of the numerical method
%
%
%   Method
%   ======
%   Weak error is obtained using the exact solution of the system obtained with the Gillespie algorithm (SSA)
%
%
%   IN
%   ==
%   1) methods    - cell array of function handles that calls corresponding numerical methods
%   2) input_data - input data of the chemical system
%   3) FTime      - final time of the simulation
%   4) h          - vector with discretization step sizes
%   6) options    - (optional) struct with options 
%       test_SSA  - (optional) is SSA simulation required?
%       NPaths    - (optional) number of MC paths
%       solver_options - (optional) cell array with options for different methods
%
%
%   OUT
%   ===
%   h     - vector with time steps
%   err   - vector with corresponding strong errors
%   order - estimated orders of convergence

   
    % default options
    op = struct('test_SSA',1,'test_methods',1,'NPaths',100,'solver_options',[]);

    % overwrite default options
    if nargin == 5 && isstruct(options)
        for fn = fieldnames(options)'
            op.(fn{1}) = options.(fn{1});
        end
    end
    
    Nm   = length(methods);  % number of methods to test
    Nlev = length(h);        % number of testing levels

    % check if options for solvers are supplied
    if length(op.solver_options) == 1 && Nm > 1
        if isstruct(op.solver_options)
            temp = op.solver_options;
        elseif iscell(op.solver_options)
            temp = op.solver_options{1};
        else 
            error('options must be a struct or cell array of structs\n')
        end
        op.solver_options = cell(1,Nm);
        for n = 1:Nm
            op.solver_options{n} = temp;
        end
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % parameters of the stochastic chemical system
    [nu,prop,Y0] = input_data();
    
    % dimension of the system
    N = length(Y0);
    
    % allocate various arrays
    err      = zeros(Nm,Nlev);       % strong errors
    Y        = cell(Nm,Nlev);        % solutions
    Y_mean   = cell(Nm,Nlev);        % mean of solutions
    Y_var    = cell(Nm,Nlev);        % variance of solutions
    tspan    = cell(1,Nlev);         % time points
	for p = 1:Nlev
        tspan{p} = 0:h(p):FTime;
    end
    
    
    fprintf('\nmethods: %s', func2str(methods{1}));
    for i = 2:Nm
        fprintf(', %s', func2str(methods{i}));
    end
    fprintf('\ndt_min = %6.4e\n', h(1));
    fprintf('dt_max = %6.4e\n', h(end));
    fprintf('Nlev   = %i\n\n', Nlev);
    
    
    poolobj = gcp('nocreate');
	if isempty(poolobj)
        parpool('local');
        poolobj = gcp('nocreate');
    end
    
    fprintf('Start time is %s\n\n',datetime('now'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % if SSA simulation is required
    if op.test_SSA
        % generate reference time grid
        tic;
        [~,SSA_tspan,SSA_K] = SSA(prop,nu,FTime,Y0);
        seconds = toc*1.2*op.NPaths/poolobj.NumWorkers;
        fprintf('SSA:\n\tEstimated execution time is %i hours %i minutes %i seconds\n',floor(seconds/3600),floor(seconds/60-60*floor(seconds/3600)),floor(seconds-60*floor(seconds/60)));
        fprintf('\tEstimated final time is %s\n\n',datetime(datevec(datetime('now'))+[0,0,0,0,0,seconds]));
    
        % approximate number of reactions
        SSA_K = round(1.2*SSA_K);
        
        % time points of the true solution
        SSA_tspan = SSA_tspan(1,round(linspace(1,size(SSA_tspan,2),1000)));
        % reference SSA solution
        SSA_Y     = zeros(N,length(SSA_tspan),op.NPaths);
    
        % generate reference solution
    	parfor path = 1:op.NPaths 
            [Y_ref_temp,poissprc] = SSA(prop,nu,FTime,Y0);
            for i = 1:N
                SSA_Y(i,:,path) = interp1(poissprc(1,:),Y_ref_temp(i,:),SSA_tspan,'previous','extrap');
            end
            %task = getCurrentTask();
            %if task.ID == 1
            %    fprintf('Path: %i of %i    SSA\n',path,op.NPaths);
            %end
        end
        SSA_mean = mean(SSA_Y,3);
        SSA_var  = var(SSA_Y,0,3);
    else
        SSA_tspan = 0;
        SSA_mean  = 0;
        SSA_var   = 0;
        SSA_Y     = 0;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if op.test_methods
        tic;
        for p = 1:Nlev
            for n = 1:Nm
                methods{n}(prop,nu,tspan{p},Y0,op.solver_options{n});
            end
        end
        seconds = toc*1.2*op.NPaths/poolobj.NumWorkers;
        fprintf('Solvers:\n\tEstimated execution time is %i hours %i minutes %i seconds\n',floor(seconds/3600),floor(seconds/60-60*floor(seconds/3600)),floor(seconds-60*floor(seconds/60)));
        fprintf('\tEstimated final time is %s\n\n',datetime(datevec(datetime('now'))+[0,0,0,0,0,seconds]));
        
        % test different levels of discretization
        for p = 1:Nlev
            % test different methods
            for n = 1:Nm
                Y_temp = zeros(N,length(tspan{p}),op.NPaths);
                % simulate paths
                parfor path = 1:op.NPaths
                    Y_temp(:,:,path) = methods{n}(prop,nu,tspan{p},Y0,op.solver_options{n});
                    %task = getCurrentTask();
                    %if task.ID == 1
                    %    fprintf('Path: %i of %i    Level: %i\n',path,op.NPaths,p);
                    %end
                end
                Y{n, p} = Y_temp;
                if op.test_SSA
                    err(n,p) = norm( mean(Y_temp(:,end,:),3) - SSA_mean(:,end) );
                end
                Y_mean{n,p} = mean(Y_temp,3);
                Y_var{n,p}  = var(Y_temp,0,3);
            end
        end
    else
        for p = 1:Nlev
            tspan{p} = [0 1];
            for n = 1:Nm
                Y_mean{n,p} = ones(N,2);
                Y_var{n,p}  = ones(N,2);
                Y{n, p}     = ones(N,2);
            end
        end
    end

    
    order = zeros(Nm, 1);
    for n = 1:Nm
        pp = polyfit( log(h), log(err(n,:)), 1 );
        order(n) = pp(1);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('Final time is %s\n\n',datetime('now'));
    
    fprintf('\n\nmethods: %s', func2str(methods{1}));
    for i = 2:Nm
        fprintf(', %s', func2str(methods{i}));
    end
    fprintf('\ndt_min = %6.4e\n', h(1));
    fprintf('dt_max = %6.4e\n', h(end));
    fprintf('Nlev   = %i\n\n', Nlev);
    

end