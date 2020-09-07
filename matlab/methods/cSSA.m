classdef cSSA < handle
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
%   Direct Stochastic Simulation Algorithm:
%
%  / y1[k+1] \   / y1[k] \   / nu_1r \
%  | y2[k+1] |   | y2[k] |   | nu_2r |
%  |    .    | = |   .   | + |   .   |,
%  |    .    |   |   .   |   |   .   |
%  \ yn[k+1] /   \ yn[k] /   \ nu_nr /
%   
%   where r = Multinomial( 1, prop(y[k]) / sum(propensities) )


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        ChemSys;    % chemical system
        Y0;         % initial condition
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=false )
        K;          % number of time points ( number of reactions fired )
        T;          % vector of jump times
        Y;          % solution
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=true )
        Tfin;       % final time
        generate;   % function handle to the function which generates solution
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cSSA(varargin)
            obj.generate = @obj.first_call_generate_all_times;
            obj.K = 0;
            
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'system'
                        obj.ChemSys = varargin{i+1};
                    case 'X0'
                        obj.Y0 = varargin{i+1};
                    case 'Tfin'
                        obj.Tfin = varargin{i+1};
                    case 'T'
                        switch numel(varargin{i+1})
                            case 1
                                obj.T    = varargin{i+1};
                                obj.Tfin = obj.T;
                                obj.generate = @obj.generate_final_time;
                            otherwise
                                obj.T    = varargin{i+1};
                                obj.Tfin = obj.T(end);
                                obj.generate = @obj.generate_given_times;
                        end
                end
            end
            
            obj.K = numel(obj.T);
            obj.Y = zeros(obj.ChemSys.num_species,obj.K);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    first call to the function 'generate_all_times'                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = first_call_generate_all_times(obj)
            obj.generate = @obj.generate_all_times;
            
            % estimate number of reactions
            if obj.K == 0
                Y1 = obj.Y0;
                dt = 0;
                for i = 1:10
                    a0 = cumsum( obj.ChemSys.propensities(Y1) );
                    dt = dt - log(rand()) / a0(end);
                    r  = find( a0 >= (a0(end)*rand()), 1 );
                    Y1 = Y1 + obj.ChemSys.nu(:,r);
                end
                dt     = dt / 10;
                obj.K  = ceil( 1.5 * obj.Tfin / dt );
            end
            
            obj.T = zeros(1,obj.K);
            obj.Y = zeros(obj.ChemSys.num_species,obj.K);
            
            time = obj.generate();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                       generate path with all jump times                                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate_all_times(obj)
            tic();
            
            obj.Y(:,1) = obj.Y0;
            i = 2;
            while obj.T(i-1) <= obj.Tfin
                % choose next reaction
                a0  = cumsum( obj.ChemSys.propensities( obj.Y(:,i-1) ) );
                rnd = rand(2,1);
                % time of the next reaction
                obj.T(i) = obj.T(i-1) - log(rnd(1)) / a0(end);
                % index of the next reaction
                r = find( a0 >= a0(end)*rnd(2), 1 );
                
                % update solution
                obj.Y(:,i) = obj.Y(:,i-1) + obj.ChemSys.nu(:,r);
          
                i = i + 1;
        
                % dynamically reallocate arrays
                if i > obj.K
                    obj.T = [obj.T zeros(2,ceil(0.5*obj.K))];
                    obj.Y = [obj.Y zeros(N,ceil(0.5*obj.K))];
                    obj.K = obj.K + ceil(0.5*obj.K);
                end
            end
            
            % remove trailing zeros
            if i <= obj.K
                obj.T(:,i:obj.K) = [];
                obj.Y(:,i:obj.K) = [];
            end
            obj.K = i-1;
            
            % interpolate solution at the endpoint
            if obj.T(1,obj.K) > obj.Tfin
                obj.T(1,obj.K) = obj.Tfin;
                obj.Y(:,obj.K) = obj.Y(:,obj.K-1);
            end
            
            time = toc();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                         generate path at given points                                   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate_given_times(obj)
            tic();
            
            obj.Y(:,1) = obj.Y0;
            Y_new      = obj.Y0;
            time_new   = 0;
            i = 2;
            while time_new <= obj.Tfin
                Y_old = Y_new;
                
                % choose next reaction
                a0  = cumsum( obj.ChemSys.propensities( Y_old ) );
                rnd = rand(2,1);
                % time of the next reaction
                time_new = time_new - log(rnd(1)) / a0(end);
                % index of the next reaction
                r = find( a0 >= a0(end)*rnd(2), 1 );
                
                Y_new = Y_old + obj.ChemSys.nu(:,r);
                
                % update solution
                while ( i <= obj.K && time_new > obj.T(i) )
                    obj.Y(:,i) = Y_old;
                    i = i + 1;
                end
            end
            
            time = toc();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     generate solution at final time only                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate_final_time(obj)
            tic();
            
            time_new = 0;
            obj.Y = obj.Y0;
            while true
                % choose next reaction
                a0  = cumsum( obj.ChemSys.propensities( obj.Y ) );
                rnd = rand(2,1);
                % time of the next reaction
                time_new = time_new - log(rnd(1)) / a0(end);
                % index of the next reaction
                r = find( a0 >= a0(end)*rnd(2), 1 );
                
                if time_new >= obj.Tfin
                    break;
                end
                
                % update solution
                obj.Y = obj.Y + obj.ChemSys.nu(:,r);
            end
            obj.T = obj.Tfin;
            
            time = toc();
        end
        
        
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

