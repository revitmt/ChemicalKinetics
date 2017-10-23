classdef cThetaTauLeap < cIntegratorBase
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
        int_states; % is integer states?
    end
    properties ( GetAccess='protected', SetAccess='protected', Hidden=true )
        J_X;        % identity matrix used in Jacobian for nonlinear solver
        tau_theta;  % tau * theta
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cThetaTauLeap(varargin)
            obj@cIntegratorBase(varargin{:});
            
            % default value (explicit tau-leap)
            obj.theta      = 0;
            obj.int_states = true;
            
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'theta'
                        obj.theta = varargin{i+1};
                    case 'int_states'
                        obj.int_states = varargin{i+1};
                end
            end
            
            obj.tau_theta = obj.tau * obj.theta;
            
            % identity matrix used in Jacobian for nonlinear solver
            if obj.theta ~= 0
                obj.J_X = eye(obj.ChemSys.num_species);
            end
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                            single step of the method                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function single_step(obj)
            obj.generate_random_increment(obj.Y(:,obj.curr_ind-1));
            
            % if implicit
            if obj.theta ~= 0
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % stochastic step
                    
                Ys = obj.Y(:,obj.curr_ind-1) + obj.ChemSys.nu * ( obj.dN - obj.a*obj.tau );
      
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % deterministic corrector

                [obj.Y(:,obj.curr_ind),~,flag,output] = fsolve( @fun, obj.Y(:,obj.curr_ind-1), obj.fslv_opt);
                if flag<=0
                    obj.Y(:,obj.curr_ind-1)
                    flag
                    output.message
                end
                obj.info.iter    = obj.info.iter    + output.iterations;
                obj.info.funeval = obj.info.funeval + output.funcCount;
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % preserve integer states
                    
                if obj.int_states
                    obj.Y(:,obj.curr_ind) = obj.Y(:,obj.curr_ind-1) + obj.ChemSys.nu * round( obj.dN + ( obj.ChemSys.propensities(obj.Y(:,obj.curr_ind)) - obj.a ) * obj.tau_theta );
                end
            % if explicit
            else
                obj.Y(:,obj.curr_ind) = obj.Y(:,obj.curr_ind-1) + obj.ChemSys.nu * obj.dN;
            end
               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % preserve nonnegative states
            
            obj.Y(:,obj.curr_ind) = max( obj.Y(:,obj.curr_ind), 0 );
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % nonlinear equation for deterministic step
            
            function [f,Jf] = fun(X)
                f  = obj.ChemSys.propensities(X);
                f  = X - Ys - obj.ChemSys.nu * ( f * obj.tau_theta );
                Jf = obj.ChemSys.propJacobian(X);
                Jf = obj.J_X - obj.ChemSys.nu * Jf * obj.tau_theta;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    update mean and covariance of the linear system                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ sigma, mu, A ] = update_mean_covariance(obj,mean_old,covariance_old)
            Ja = obj.ChemSys.propJacobian(obj.Y(:,obj.curr_ind-1));
            nu = obj.ChemSys.nu;

            R1 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
            R3 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
            A  = R3 * R1;

            temp  = R3 * nu;
            mu    = A * mean_old;
            sigma = A * covariance_old * A' + obj.tau * temp * diag(Ja*R1*mean_old) * temp';
        end
    end
    
end

