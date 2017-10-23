classdef cCompositeIntegrator < handle
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
%   Composite integrator
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        ChemSys;            % chemical system
        integrators;        % array of pointers to integrators at each subinterval
    end
    properties ( GetAccess='public', SetAccess='protected', Hidden=false )
        num_intervals;      % number of subintervals (always one in this case)
        info;               % struct with info about solver
    end
    properties ( Dependent, GetAccess='public', SetAccess='protected', Hidden=false )
        Y0;                 % initial condition
        Yfin;               % solution at final time
        T0;                 % array of initial times at each subinterval
        Tfin;               % array of final times at each subinterval
    end
    
    % properties for coupled paths
    properties ( GetAccess='public', SetAccess='public', Hidden=true )
        coupledPath;        % pointer to the coupled integrator
        generate_coarse_no; % array of pointers to coarse integrators at each subinterval
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cCompositeIntegrator(varargin)
            obj.num_intervals = nargin;
            obj.integrators   = cell(1,obj.num_intervals);
            
            for i = 1:nargin
                obj.integrators{i} = varargin{i};
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     setters and getters for dependent properties                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %function value = get.num_intervals(obj)
        %    value = numel(obj.integrators);
        %end
        function value = get.T0(obj)
            value = zeros(1,obj.num_intervals);
            for i = 1:numel(value)
                value(i) = obj.integrators{i}.T0;
            end
        end
        function value = get.Tfin(obj)
            value = zeros(1,obj.num_intervals);
            for i = 1:numel(value)
                value(i) = obj.integrators{i}.Tfin;
            end
        end
        function value = get.Y0(obj)
            value = obj.integrators{1}.Y0;
        end
        function value = get.Yfin(obj)
            value = obj.integrators{obj.num_intervals}.Yfin;
        end
        function value = get.ChemSys(obj)
            value = obj.integrators{1}.ChemSys;
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                 add integrator                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function add_integrator(obj, integrator)
            obj.num_intervals = obj.num_intervals + 1;
            obj.integrators{end+1} = integrator;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  generate path                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate(obj)
            tt = tic();
            
            obj.integrators{1}.generate();
            obj.info.iter    = obj.integrators{1}.info.iter;
            obj.info.funeval = obj.integrators{1}.info.funeval;
            for i = 2:obj.num_intervals
                obj.integrators{i}.Y0 = obj.integrators{i-1}.Y(:,end);
                obj.integrators{i}.generate();
                obj.info.iter    = obj.info.iter    + obj.integrators{i}.info.iter;
                obj.info.funeval = obj.info.funeval + obj.integrators{i}.info.funeval;
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
            for i = 1:obj.num_intervals
                obj.integrators{i}.coupledPath             = obj.coupledPath.integrators{i};
                obj.coupledPath.integrators{i}.coupledPath = obj.integrators{i};
                
                obj.generate_coarse_no{i} = obj.integrators{i}.choose_coarse_path();
            end
            fhandle = @obj.generate_coarse;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  generate path                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generate_coarse(obj)
            obj.generate_coarse_no{1}();
            
            obj.info             = obj.integrators{1}.info;
            obj.coupledPath.info = obj.integrators{1}.coupledPath.info;

            for i = 2:obj.num_intervals
                obj.integrators{i}.Y0 = obj.integrators{i-1}.Y(:,end);
                obj.coupledPath.integrators{i}.Y0 = obj.coupledPath.integrators{i-1}.Y(:,end);
                obj.generate_coarse_no{i}();

                obj.info.iter    = obj.info.iter    + obj.integrators{i}.info.iter;
                obj.info.funeval = obj.info.funeval + obj.integrators{i}.info.funeval;
                obj.info.time    = obj.info.time    + obj.integrators{i}.info.time;

                obj.coupledPath.info.iter    = obj.coupledPath.info.iter    + obj.integrators{i}.coupledPath.info.iter;
                obj.coupledPath.info.funeval = obj.coupledPath.info.funeval + obj.integrators{i}.coupledPath.info.funeval;
                obj.coupledPath.info.time    = obj.coupledPath.info.time    + obj.integrators{i}.coupledPath.info.time;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                 plot solution                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot(obj,var)
            T1 = [];
            Y1 = [];
            for i = 1:obj.num_intervals
                T1 = [ T1 obj.integrators{i}.T ];
                Y1 = [ Y1 obj.integrators{i}.Y ];
            end
            
            if nargin == 1
                for s = 1:obj.integrators{i}.ChemSys.num_species
                    plot(T1,Y1(s,:),'LineWidth',2); hold on;
                end
                legend(obj.integrators{1}.ChemSys.species);
            else
                for s = var
                    plot(T1,Y1(s,:),'LineWidth',2); hold on;
                end
                legend(obj.integrators{1}.ChemSys.species(var));
            end
            set(gca,'FontSize',15);
            xlim([obj.integrators{1}.T0 obj.integrators{obj.num_intervals}.Tfin]);
            xlabel('time');
            hold off;
        end
    end
    
end

