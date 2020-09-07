classdef cThetaLangevin < cSinglePathBase
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
%   Theta tau-leap.
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        theta;      % implicitness parameter
    end
    properties ( GetAccess='private', SetAccess='private', Hidden=true )
        J_X;        % identity matrix used in Jacobian for nonlinear solver
        tau_theta;  % tau * theta
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cThetaLangevin(varargin)
            obj@cSinglePathBase(varargin{:});
            
            % default value (explicit tau-leap)
            obj.theta = 0;
            
            for i = 1:2:nargin-1
                option = varargin{i};
                switch option
                    case 'theta'
                        obj.theta = varargin{i+1};
                end
            end
            
            % identity matrix used in Jacobian for nonlinear solver
            if obj.theta ~= 0
                obj.J_X = eye(obj.ChemSys.num_species);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                  generate path                                          %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = generate(obj)
            tt = tic();

            obj.tau_theta = obj.tau * obj.theta;
            
            obj.info.iter    = 0;
            obj.info.funeval = 0;
            
            obj.Y(:,1) = obj.Y0;
            for i = 2:obj.K
                a = obj.ChemSys.propensities( max( obj.Y(:,i-1), 0 ) );
                
                % generate random Poisson increments
                dW = sqrt(obj.tau) * sqrt(a) .* randn( obj.ChemSys.num_reactions, 1 );
        
                % if implicit
                if obj.theta ~= 0
                    % stochastic step
                    Ys = obj.Y(:,i-1) + obj.ChemSys.nu * dW;
      
                    % deterministic corrector
                    [obj.Y(:,i),~,~,output] = fsolve( @fun, obj.Y(:,i-1), obj.fslv_opt);
                    obj.info.iter    = obj.info.iter    + output.iterations;
                    obj.info.funeval = obj.info.funeval + output.funcCount;
                    
                % if explicit
                else
                    obj.Y(:,i) = obj.Y(:,i-1) + obj.ChemSys.nu * ( a * obj.tau + dW );
                end
            end
            
            function [f,Jf] = fun(X)
                f  = obj.ChemSys.propensities(X);
                f  = X - Ys - obj.ChemSys.nu * ( f * obj.tau_theta + a * (1-obj.theta) * obj.tau );
                Jf = obj.ChemSys.propJacobian(X);
                Jf = obj.J_X - obj.ChemSys.nu * Jf * obj.tau_theta;
            end
            
            time = toc(tt);
            obj.info.time = time;
        end
    end
    
end

