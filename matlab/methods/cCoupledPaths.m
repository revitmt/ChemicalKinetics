classdef cCoupledPaths < handle
%	v1
%   Copyright 2017, Viktor Reshiak, All rights reserved.    
%
%
%   Purpose
%   =======
%   Generate coupled paths of the stochastic chemical system.
%
%
%   Method
%   ======
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='private', Hidden=false )
        path1;
        path2;
        Q;                  % mesh spacing parameter
        info;               % struct with info about solver
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=true )
        generate;           % function handle to generate function
    end
    properties ( Dependent, GetAccess='public', SetAccess='private', Hidden=false )
        Yfin;               % solution at final time
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cCoupledPaths(varargin) %system,path1,path2)
            if nargin==2
                obj.path1 = varargin{1};
                obj.path2 = varargin{2};
                obj.path1.coupledPath = obj.path2;
                obj.path2.coupledPath = obj.path1;
            elseif nargin==3
                obj.path1 = varargin{2};
                obj.path2 = varargin{3};
                obj.path1.ChemSys = varargin{1};
                obj.path2.ChemSys = varargin{1};
                obj.path1.coupledPath = obj.path2;
                obj.path2.coupledPath = obj.path1;
            end
            
            obj.generate = @obj.generate_first_call;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                          change parameters of both paths                                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_parameters(obj,varargin)
            obj.path1.set_parameters(varargin{:});
            obj.path2.set_parameters(varargin{:});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           change parameters of path1                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_parameters_for_path1(obj,varargin)
            obj.path1.set_parameters(varargin{:});
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           change parameters of path2                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set_parameters_for_path2(obj,varargin)
            obj.path2.set_parameters(varargin{:});
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                     setters and getters for dependent properties                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function value = get.Yfin(obj)
            value{1} = obj.path1.Yfin;
            value{2} = obj.path2.Yfin;
        end
        function value = get.info(obj)
            value.iter    = obj.path1.info.iter    + obj.path2.info.iter;
            value.funeval = obj.path1.info.funeval + obj.path2.info.funeval;
            value.time    = obj.path1.info.time    + obj.path2.info.time;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                 generate paths                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function generate_first_call(obj)
            obj.generate = obj.path1.choose_coarse_path();
            %obj.Q = max( obj.path1.Q, obj.path2.Q );
            obj.generate();
        end
                
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                 plot solution                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot(obj,var)
            obj.path1.plot(var); hold on;
            obj.path2.plot(var); hold off;
        end
    end
    
end

