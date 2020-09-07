classdef cMLMC < handle
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%
%   Purpose
%   =======
%   Simulate the chemical system
%
%
%   Method
%   ======
%   Multilevel Monte Carlo method
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        ChemSys;    % chemical system
        Y0;         % initial condition
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=false )
        T0;
        Tfin;
        Nlevels;    % number of levels
        Q;          % mesh refinement
        solvers;    % array of solvers at each level
        MC;         % array of Monte Carlo solvers at each level
        Nsamples;   % array of number of MC samples at each level
        fun;        % functional of the solution to be simulated
        stat;       % struct with generated statistic
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=true )
        info;       % struct with various parameters of the method
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cMLMC(varargin)
            % default parameter values
            obj.Q  = 2;
            obj.T0 = 0;
            
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'system'
                        obj.ChemSys = varargin{i+1};
                    case 'T0'
                        obj.T0 = varargin{i+1};
                    case 'Tfin'
                        obj.Tfin = varargin{i+1};
                    case 'X0'
                        obj.Y0 = varargin{i+1};
                    case 'functional'
                        obj.fun = varargin{i+1};
                    case 'Nlevels'
                        obj.Nlevels = varargin{i+1};
                    case 'tau'
                        obj.Nlevels = ceil( log2( varargin{i+1} / ( obj.Tfin - obj.T0 ) ) );
                    case 'refinement'
                        obj.Q = varargin{i+1};
                end
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                run simulation                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function run(obj)

        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                estimate cost                                            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function choose_solvers(obj)
            fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf('|  lev    tau      q     K_c  |     C_l       sig_l       mu_l    |     C_0       sig_0       mu_0   |     f_l        f_0     |    Cost      dCost    |\n');
            fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            
            int_states = false;
            dCost = 0;
            for lev = obj.Nlevels:-1:1
                obj.info.Nsamples(lev) = 100;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % explicit part
            
            for lev = obj.Nlevels:-1:1
                N_steps1 = obj.Q^(lev-1);
                N_steps2 = obj.Q^lev;
                
                obj.solvers{lev} = cCoupledPaths( obj.ChemSys, cThetaTauLeap, cThetaTauLeap );
                obj.solvers{lev}.set_parameters( 'X0',         obj.Y0,      ...
                                                 'Tfin',       obj.Tfin,    ...
                                                 'theta',      0,           ...
                                                 'int_states', int_states );
                obj.solvers{lev}.path1.set_parameters( 'N_steps', N_steps1 );
                obj.solvers{lev}.path2.set_parameters( 'N_steps', N_steps2 );
     
                tau = obj.solvers{lev}.path1.tau;
                run_MC();

                if dCost > ( 0.05 * obj.info.tot_cost_L0(lev) )
                    obj.info.interface_level = lev;
                    fprintf('  <-- remove this level\n');
                    break;
                else
                    fprintf('\n');
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interface part
            
            % how much implicit is more expensive than explicit?
            Solver1 = cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',obj.Q^(lev-1),'theta',1,'int_states',int_states);
            Solver2 = cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',obj.Q^(lev-1),'theta',0,'int_states',int_states);
            Solver1.generate();
            Solver2.generate();
            cost_ratio = Solver1.info.time / Solver2.info.time;
            
            %lev_jump = round( log(Solver1.info.time / Solver2.info.time) / log(obj.Q) ) ;
                
            fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------------\n');
                % relaxation time
                time_relax    = 0.01 * obj.Tfin;
                N_steps_relax = ceil( time_relax / obj.solvers{lev}.path2.tau );
                time_relax    = N_steps_relax * obj.solvers{lev}.path2.tau;
                
                %N_steps1 = obj.Q^(lev-1-lev_jump);
                %N_steps1 = N_steps1 * ( 1 - time_relax / obj.Tfin );
                %N_steps1 = obj.Q^floor( log( N_steps1 ) / log(obj.Q) );
                
                N_steps2 = obj.Q^lev - N_steps_relax;
                N_steps1 = obj.Q^( floor( log( N_steps2/cost_ratio ) / log(obj.Q) ) - 2 );
                lev_jump = lev - 1 - floor(log(N_steps1)/log(obj.Q)) + 1;

                fprintf('Interface level: lev = %d, T_relax = %9.2e, K_relax = %d, K_c = %d, K_f = %d \n',lev - 1 - lev_jump, time_relax, N_steps_relax, N_steps1, N_steps2);
                
                Solver1 = cThetaTauLeap('system',obj.ChemSys,'T0',obj.Tfin-time_relax,'Tfin',obj.Tfin,'N_steps',N_steps_relax,'theta',0,'int_states',int_states);
                Solver2 = cThetaTauLeap('system',obj.ChemSys,'T0',obj.Tfin-time_relax,'Tfin',obj.Tfin,'N_steps',N_steps_relax,'theta',0,'int_states',int_states);
                
                obj.solvers{lev} = cCoupledPaths( ...
                    cCompositeIntegrator( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin-time_relax,'N_steps',N_steps1,'theta',1,'int_states',int_states), Solver1 ), ...
                    cCompositeIntegrator( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin-time_relax,'N_steps',N_steps2,'theta',0,'int_states',int_states), Solver2 ) );
  
%                 obj.solvers{lev} = cCoupledPaths( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',N_steps1,'theta',1,'int_states',int_states), ...
%                                                   cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',N_steps2,'theta',0,'int_states',int_states) );

%                 tau = obj.solvers{lev}.path1.tau;
                tau = obj.solvers{lev}.path1.integrators{1}.tau;
                
                run_MC();
                fprintf('\n');
            fprintf('-------------------------------------------------------------------------------------------------------------------------------------------------------\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % implicit part
            
            for lev = obj.info.interface_level-1:-1:1+lev_jump
                N_steps1 = obj.Q^(lev-1-lev_jump);
                N_steps2 = obj.Q^(lev-lev_jump);
%                 obj.solvers{lev} = cCoupledPaths( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',N_steps1,'theta',1,'int_states',int_states), ...
%                                                   cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin,'N_steps',N_steps2,'theta',1,'int_states',int_states) );
%                 tau = obj.solvers{lev}.path1.tau;
                obj.solvers{lev} = cCoupledPaths( ...
                    cCompositeIntegrator( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin-time_relax,'N_steps',N_steps1,'theta',1,'int_states',int_states), Solver1 ), ...
                    cCompositeIntegrator( cThetaTauLeap('system',obj.ChemSys,'X0',obj.Y0,'Tfin',obj.Tfin-time_relax,'N_steps',N_steps2,'theta',1,'int_states',int_states), Solver2 ) );
                
                tau = obj.solvers{lev}.path1.integrators{1}.tau;
                run_MC();
                fprintf('\n');
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % common MC part
            
            function run_MC()
                obj.MC{lev} = cMonteCarlo('solver',obj.solvers{lev},'functional',obj.fun);
                obj.MC{lev}.run( obj.info.Nsamples(lev) );
                
                obj.info.mu_l(lev)    = obj.MC{lev}.stat.mu;
                obj.info.sigma_l(lev) = obj.MC{lev}.stat.sigma;
                obj.info.kurt_l(lev)  = obj.MC{lev}.stat.kurt;
                obj.info.cost_l(lev)  = obj.MC{lev}.stat.cost;
                
                obj.info.mu_c(lev)    = obj.MC{lev}.stat.mu1;
                obj.info.sigma_c(lev) = obj.MC{lev}.stat.sigma1;
                obj.info.kurt_c(lev)  = obj.MC{lev}.stat.kurt1;
                obj.info.cost_c(lev)  = obj.MC{lev}.stat.cost1;
                
                obj.info.lev_cost_l(lev) = sqrt( obj.info.sigma_l(lev) * obj.info.cost_l(lev) );
                obj.info.lev_cost_c(lev) = sqrt( obj.info.sigma_c(lev) * obj.info.cost_c(lev) );
                
                obj.info.tot_cost_L0(lev) = ( obj.info.lev_cost_c(lev) + sum( obj.info.lev_cost_l(lev:end) ) )^2;
                
                if lev==obj.Nlevels
                    dCost = 0;
                else
                    dCost = obj.info.tot_cost_L0(lev) - obj.info.tot_cost_L0(lev+1);
                end

                fprintf('|  %3d %9.2e %4d %6d  | %10.3e %10.3e %10.3e  | %10.3e %10.3e %10.3e | %10.3e %10.3e  | %10.3e %10.3e |', ...
                         lev, tau, obj.Q, N_steps1, ...
                         obj.info.cost_l(lev), obj.info.sigma_l(lev), obj.info.mu_l(lev),...
                         obj.info.cost_c(lev), obj.info.sigma_c(lev), obj.info.mu_c(lev),...
                         obj.info.lev_cost_l(lev), obj.info.lev_cost_c(lev), ...
                         obj.info.tot_cost_L0(lev), dCost );
            end
        end
        
        
        function plot_params(obj)
            L_0 = find(obj.info.sigma_c,1);
            l_i = obj.info.interface_level;
            
            fig_num = 1;
            
            figure(fig_num); fig_num = fig_num + 1;
            %ss0(L_0:obj.Nlevels) = smooth(obj.info.sigma_c(L_0:obj.Nlevels));
            y_c(L_0:obj.Nlevels) = (obj.info.sigma_c(L_0:obj.Nlevels));
            y_l(L_0:obj.Nlevels) = (obj.info.sigma_l(L_0:obj.Nlevels));
            plot_common_semilog();
            xlabel('level');
            xlim([ L_0 obj.Nlevels ]);
            legend('V_l^{''}','V_l','Location','east');
            title('Variance');
            set(gca,'FontSize',15);
            
            figure(fig_num); fig_num = fig_num + 1;
            y_c(L_0:obj.Nlevels) = (obj.info.mu_c(L_0:obj.Nlevels));
            y_l(L_0:obj.Nlevels) = (obj.info.mu_l(L_0:obj.Nlevels));
            plot_common();
            xlabel('level');
            xlim([ L_0 obj.Nlevels ]);
            legend('\mu_l^{''}','\mu_l','Location','east');
            title('Mean');
            set(gca,'FontSize',15);
            
            figure(fig_num); fig_num = fig_num + 1;
            y_c(L_0:obj.Nlevels) = (obj.info.cost_c(L_0:obj.Nlevels));
            y_l(L_0:obj.Nlevels) = (obj.info.cost_l(L_0:obj.Nlevels));
            plot_common_semilog();
            xlabel('level');
            xlim([ L_0 obj.Nlevels ]);
            legend('C_l^{''}','C_l','Location','northwest');
            title('Cost per level');
            set(gca,'FontSize',15);

            figure(fig_num); fig_num = fig_num + 1;
            y_c(L_0:obj.Nlevels) = (obj.info.kurt_c(L_0:obj.Nlevels));
            y_l(L_0:obj.Nlevels) = (obj.info.kurt_l(L_0:obj.Nlevels));
            plot_common_semilog();
            xlabel('level');
            xlim([ L_0 obj.Nlevels ]);
            legend('K_l^{''}','K_l','Location','east');
            title('Kurtosis');
            set(gca,'FontSize',15);

            figure(fig_num); fig_num = fig_num + 1;
            y_c(L_0:obj.Nlevels) = (obj.info.tot_cost_L0(L_0:obj.Nlevels));
            semilogy( (l_i+1:obj.Nlevels), y_c(l_i+1:obj.Nlevels),'--b', 'MarkerEdgeColor','b','LineWidth',2); hold on
            semilogy( (L_0:l_i), y_c(L_0:l_i),'--r','MarkerEdgeColor','r','LineWidth',2); hold on
            semilogy( (l_i:l_i+1), y_c(l_i:l_i+1),'--g','MarkerEdgeColor','g','LineWidth',2); hold on
            xlabel('Coarsest level');
            xlim([ L_0 obj.Nlevels ]);
            title('Cost as a function of coarsest level');
            set(gca,'FontSize',15);
            
            function plot_common_semilog()
                % plot explicit part
                semilogy( (l_i+1:obj.Nlevels), y_c(l_i+1:obj.Nlevels),'--b', ...
                          (l_i+1:obj.Nlevels), y_l(l_i+1:obj.Nlevels),'-b', 'MarkerEdgeColor','b','LineWidth',2); hold on
                % plot implicit part
                semilogy( (L_0:l_i), y_c(L_0:l_i),'--r', ...
                          (L_0:l_i), y_l(L_0:l_i),'-r', 'MarkerEdgeColor','r','LineWidth',2); hold on
                % plot interface part
                semilogy( (l_i:l_i+1), y_c(l_i:l_i+1),'--g', ...
                          (l_i:l_i+1), y_l(l_i:l_i+1),'-g', 'MarkerEdgeColor','g','LineWidth',2); hold on
            end
            function plot_common()
                % plot explicit part
                plot( (l_i+1:obj.Nlevels), y_c(l_i+1:obj.Nlevels),'--b', ...
                      (l_i+1:obj.Nlevels), y_l(l_i+1:obj.Nlevels),'-b', 'MarkerEdgeColor','b','LineWidth',2); hold on
                % plot implicit part
                plot( (L_0:l_i), y_c(L_0:l_i),'--r', ...
                      (L_0:l_i), y_l(L_0:l_i),'-r', 'MarkerEdgeColor','r','LineWidth',2); hold on
                % plot interface part
                plot( (l_i:l_i+1), y_c(l_i:l_i+1),'--g', ...
                      (l_i:l_i+1), y_l(l_i:l_i+1),'-g', 'MarkerEdgeColor','g','LineWidth',2); hold on
            end
        end
    end
    
end

