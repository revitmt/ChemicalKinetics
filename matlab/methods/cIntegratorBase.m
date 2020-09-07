classdef (Abstract) cIntegratorBase < handle
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%
%   Purpose
%   =======
%   Find solution of a stochastic chemical system.
%
%
%   Method
%   ======
%   Base class for integrators
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        ChemSys;            % chemical system
        Y0;                 % initial condition
    end
    properties ( GetAccess='public', SetAccess='protected', Hidden=false )
        tau;                % time step
        T;                  % vector of time points
        Y;                  % solution
        info;               % struct with info about solver
    end
    properties ( GetAccess='public', SetAccess='protected', Hidden=true )
        num_intervals=1;    % number of subintervals (always one in this case)
        K;                  % number of time points
        T0;                 % initial time
        Tfin;               % final time
    end
    properties ( GetAccess='protected', SetAccess='protected', Hidden=true )
        fslv_opt;           % options of the nonlinear solver
        curr_ind;           % index of the current time instance
        curr_time;          % current time instance
        a;                  % vector with current propensity values
        dN;                 % random increment
        generate_random_increment;  % pointer to the function generating random increments
    end
    properties ( Dependent, Hidden=true )
        N_steps;            % number of time steps
        Yfin;               % solution at final time
    end
    
    % properties for coupled paths
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        Q;                  % mesh refinement parameter
    end
	properties ( GetAccess='public', SetAccess='public', Hidden=true )
        coupledPath;        % pointer to the coupled integrator
        generate_fine;      % pointer to the coupled fine path gerator
        tau_pre;            % part of the fine step before the point on the coarse grid
        tau_post;           % part of the fine step after  the point on the coarse grid
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods (Abstract)
        single_step(obj)
    end
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cIntegratorBase(varargin)
            % default parameters
            obj.T0 = 0;

            % set parameters of the integrator
            if nargin == 1
                obj.ChemSys = varargin{1};
            else
                obj.set_parameters(varargin{:});
            end
            
            % default function for single path
            obj.generate_random_increment = @obj.default_random_increment;
            
            % options of nonlinear solver
            obj.fslv_opt = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','FunValCheck','off','MaxIter',100,'SpecifyObjectiveGradient',true);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         change parameters of the integrator                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_parameters(obj,varargin)
            reinit_T = false;

            for i = 1:2:nargin-1
                option = varargin{i};
                switch option
                    case 'system'
                        obj.ChemSys = varargin{i+1};
                    case 'X0'
                        obj.Y0 = varargin{i+1};
                    case 'T0'
                        obj.T0 = varargin{i+1};
                    case 'Tfin'
                        obj.Tfin = varargin{i+1};
                    case 'N_steps'
                        obj.K = varargin{i+1} + 1;
                        reinit_T = true;
                    case 'tau'
                        obj.K = ceil( ( obj.Tfin - obj.T0 ) / varargin{i+1} ) + 1;
                        reinit_T = true;
                end
            end
            
            if reinit_T
                obj.tau = ( obj.Tfin - obj.T0 ) / ( obj.K - 1 );
                obj.T   = obj.T0:obj.tau:obj.Tfin;
                obj.Y   = zeros(obj.ChemSys.num_species,obj.K);
            end
        end
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     setters and getters for dependent properties                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.N_steps(obj,value)
            obj.K   = value + 1;
            obj.tau = ( obj.Tfin - obj.T0 ) / value;
            obj.T   = obj.T0:obj.tau:obj.Tfin;
            obj.Y   = zeros(obj.ChemSys.num_species,obj.K);
        end
        function value = get.N_steps(obj)
            value = obj.K - 1;
        end
        function value = get.Yfin(obj)
            value = obj.Y(:,end);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                            generate random increment                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function default_random_increment(obj,Y)
            obj.a  = obj.ChemSys.propensities(Y);
            obj.dN = poissrnd( obj.a * obj.tau );
        end
        function coarse_path_random_increment(obj,Y)
            obj.a = obj.ChemSys.propensities(Y);
            obj.generate_fine();
        end
        function fine_matching_path_random_increment(obj,Y)
            obj.a              = obj.ChemSys.propensities(Y);
            a_common           = min(obj.coupledPath.a, obj.a );
            dN_common          = poissrnd( a_common * obj.tau );
            obj.dN             =                      dN_common + poissrnd( ( obj.a             - a_common ) * obj.tau );
            obj.coupledPath.dN = obj.coupledPath.dN + dN_common + poissrnd( ( obj.coupledPath.a - a_common ) * obj.tau );
        end  
        function fine_nonmatching_path_random_increment(obj,Y)
            obj.a    = obj.ChemSys.propensities(Y);
            a_common = min( obj.coupledPath.a, obj.a );
            
            if obj.tau_pre
                dN_common          = poissrnd( a_common * obj.tau_pre );
                obj.dN             =                      dN_common + poissrnd( ( obj.a             - a_common ) * obj.tau_pre );
                obj.coupledPath.dN = obj.coupledPath.dN + dN_common + poissrnd( ( obj.coupledPath.a - a_common ) * obj.tau_pre );
            elseif obj.tau_post
                dN_common          = poissrnd( a_common * obj.tau_post );
                obj.dN             = obj.dN + dN_common + poissrnd( ( obj.a             - a_common ) * obj.tau_post );
                obj.coupledPath.dN =          dN_common + poissrnd( ( obj.coupledPath.a - a_common ) * obj.tau_post );
            else
                dN_common          = poissrnd( a_common * obj.tau );
                obj.dN             =                      dN_common + poissrnd( ( obj.a             - a_common ) * obj.tau );
                obj.coupledPath.dN = obj.coupledPath.dN + dN_common + poissrnd( ( obj.coupledPath.a - a_common ) * obj.tau );
            end
        end  
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  generate path                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate(obj)
            tt = tic();
            
            obj.info.iter    = 0;
            obj.info.funeval = 0;
            
            obj.Y(:,1) = obj.Y0;
            for i = 2:obj.K
                % current index and time instance
                obj.curr_ind  = i;
                obj.curr_time = obj.T(i);
                
                obj.single_step();
            end
            
            time = toc(tt);
            obj.info.time = time;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % functions for coupled paths
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       entry point to the generate function                              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fhandle = choose_coarse_path(obj)
            % if 'obj' has coarser grid than 'obj.coupledPath'
            if (obj.N_steps <= obj.coupledPath.N_steps)
                path_c = obj;
                path_f = obj.coupledPath;
            % if 'obj.coupledPath' has coarser grid than 'obj'
            else
                path_f = obj;
                path_c = obj.coupledPath;
            end

            fhandle = @path_c.generate_coarse;
            path_c.generate_random_increment = @path_c.coarse_path_random_increment;
            path_c.Q = path_c.N_steps / path_f.N_steps;
            path_f.Q = path_f.N_steps / path_c.N_steps;
            if ~mod(path_f.N_steps,path_c.N_steps)
                path_c.generate_fine = @path_f.generate_matching_fine;
                path_f.generate_random_increment = @path_f.fine_matching_path_random_increment;
            else
                path_c.generate_fine = @path_f.generate_nonmatching_fine;
                path_f.generate_random_increment = @path_f.fine_nonmatching_path_random_increment;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                             outer 'coarse' integrator                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generate_coarse(obj)
            tt = tic();
            
            obj.coupledPath.info.iter    = 0;
            obj.coupledPath.info.funeval = 0;
            obj.coupledPath.info.time    = 0;
            obj.info.iter                = 0;
            obj.info.funeval             = 0;
            
            obj.coupledPath.curr_ind     = 1;
            obj.coupledPath.curr_time    = obj.coupledPath.T0;
            obj.curr_ind                 = 1;
            obj.curr_time                = obj.T0;

            obj.coupledPath.tau_post = 0;
            
            obj.coupledPath.Y(:,1) = obj.coupledPath.Y0;
            obj.Y(:,1)             = obj.Y0;
            for i = 2:obj.K
                % current index and time instance
                obj.curr_ind  = i;
                obj.curr_time = obj.T(i);
                
                obj.single_step();
            end
            
            obj.info.time = toc(tt) - obj.coupledPath.info.time;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                             inner 'fine' integrator                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generate_matching_fine(obj)
            tt = tic();
            
            obj.coupledPath.dN = zeros(size(obj.coupledPath.a));
            
            for i = 1:obj.Q
                % current index and time instance
                obj.curr_ind  = obj.curr_ind + 1;
                obj.curr_time = obj.T(i);
                
                obj.single_step();
            end
            
            obj.info.time = obj.info.time + toc(tt);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                             inner 'fine' integrator                                     %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generate_nonmatching_fine(obj)
            tt = tic();
            
            % when the point on coarse grid is inside the interval on fine grid
            if obj.tau_post
                obj.tau_pre  = 0;
                obj.single_step();
                obj.tau_post = 0;
            else
                obj.coupledPath.dN = zeros(size(obj.coupledPath.a));
            end
            
            while 1
                % current index and time instance
                obj.curr_ind  = obj.curr_ind + 1;
                obj.curr_time = obj.T(obj.curr_ind);
                
                if ( obj.curr_time < obj.coupledPath.curr_time )
                    obj.single_step();
                else
                    obj.tau_pre  = obj.coupledPath.curr_time - obj.T(obj.curr_ind-1);
                    obj.single_step();
                    obj.tau_post = obj.tau - obj.tau_pre;
                    break;
                end
            end
            
            obj.info.time = obj.info.time + toc(tt);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                 plot solution                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot(obj,var)
            if nargin == 1
                for s = 1:obj.ChemSys.num_species
                    plot(obj.T,obj.Y(s,:),'LineWidth',2); hold on;
                end
                set(gca,'FontSize',15);
                xlim([0 obj.Tfin]);
                xlabel('time');
                legend(obj.ChemSys.species);
                hold off;
            else
                for s = var
                    plot(obj.T,obj.Y(s,:),'LineWidth',2); hold on;
                end
                set(gca,'FontSize',15);
                xlim([0 obj.Tfin]);
                xlabel('time');
                legend(obj.ChemSys.species(var));
                hold off;
            end
        end
    end
    
end

