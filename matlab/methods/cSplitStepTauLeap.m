classdef cSplitStepTauLeap < cIntegratorBase
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
%   Split-step tau-leap.
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='public', Hidden=false )
        theta;      % splitting between predcitor and corrector steps
        eta1;       % implicitness parameter of predictor
        eta2;       % implicitness parameter of corrector
        int_states; % is integer states?
    end
    properties ( GetAccess='protected', SetAccess='protected', Hidden=true )
        J_X;        % identity matrix used in Jacobian for nonlinear solver
        SS_theta;   % adaptive split-step theta parameter
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cSplitStepTauLeap(varargin)
            obj@cIntegratorBase(varargin{:});
            
            % default value (explicit tau-leap)
            obj.theta      = -1;
            obj.eta1       = 1;
            obj.eta2       = 1;
            obj.int_states = true;
            
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'theta'
                        obj.theta = varargin{i+1};
                    case 'eta1'
                        obj.eta1 = varargin{i+1};
                    case 'eta2'
                        obj.eta2 = varargin{i+1};
                    case 'int_states'
                        obj.int_states = varargin{i+1};
                end
            end
            
            % adaptive split-step theta parameter
            obj.SS_theta = obj.theta * ones(obj.ChemSys.num_reactions,1);
            
            % identity matrix used in Jacobian for nonlinear solver
            if obj.theta ~= 0
                obj.J_X = eye(obj.ChemSys.num_species);
            end
            
            obj.a = zeros(obj.ChemSys.num_reactions,1);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                            single step of the method                                    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function single_step(obj)
            % choose theta
            if obj.theta < 0
                z   = obj.ChemSys.relaxRates(obj.Y(:,obj.curr_ind-1)) * obj.tau;
                ind = find(z<2);
                obj.SS_theta      = sqrt(2./z) - 1./z;
                obj.SS_theta(ind) = 0.633975 - 0.0566243 * z(ind);
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % deterministic predictor
                
            a_prev = obj.ChemSys.propensities(obj.Y(:,obj.curr_ind-1));
            if ( obj.eta1 ~= 0 && obj.theta ~= 1 )
                [Yd,~,~,output]  = fsolve( @fun_pred, obj.Y(:,obj.curr_ind-1), obj.fslv_opt); 
                obj.info.iter    = obj.info.iter    + output.iterations;
                obj.info.funeval = obj.info.funeval + output.funcCount;
            else
                Yd = obj.Y(:,obj.curr_ind-1) + obj.ChemSys.nu * ( a_prev .* (1-obj.SS_theta) * obj.tau );
            end
            Yd = max( Yd, 0 );
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % stochastic step

            obj.generate_random_increment(Yd);
            Ys = Yd + obj.ChemSys.nu * ( obj.dN - obj.a*obj.tau );
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % deterministic corrector
                
            a_s = obj.ChemSys.propensities(Ys);
            if obj.eta2 ~= 0
                [obj.Y(:,obj.curr_ind),~,~,output] = fsolve( @fun_corr, obj.Y(:,obj.curr_ind-1), obj.fslv_opt);
                obj.info.iter    = obj.info.iter    + output.iterations;
                obj.info.funeval = obj.info.funeval + output.funcCount;
            else
                obj.Y(:,obj.curr_ind) = Ys - obj.ChemSys.nu * ( a_s .* obj.SS_theta * obj.tau );
            end
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % preserve integer states
                
            if obj.int_states
                if obj.theta ~= 1
                    prop_d_s = (1-obj.SS_theta) .* ( obj.eta1*obj.a + (1-obj.eta1)*a_prev  );
                else
                    prop_d_s = zeros(obj.ChemSys.num_reactions,1);
                end
                prop_d_s = prop_d_s + obj.SS_theta  .* ( obj.eta2*obj.ChemSys.propensities(obj.Y(:,obj.curr_ind)) + (1-obj.eta2)*a_s );
                obj.Y(:,obj.curr_ind) = obj.Y(:,obj.curr_ind-1) + obj.ChemSys.nu * round( obj.dN + (prop_d_s - obj.a) * obj.tau );
            end
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % preserve nonnegative states
            
            obj.Y(:,obj.curr_ind) = max( obj.Y(:,obj.curr_ind), 0 );    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % nonlinear equations for deterministic predictor and corrector
            
            function [f,Jf] = fun_pred(X)
                f  = obj.ChemSys.propensities(X);
                f  = X - obj.Y(:,obj.curr_ind-1) - obj.ChemSys.nu * ( ( obj.eta1 * f + (1-obj.eta1) * a_prev ) .* (1-obj.SS_theta) ) * obj.tau;
                Jf = obj.ChemSys.propJacobian(X);
                Jf = obj.J_X - obj.ChemSys.nu * bsxfun(@times,Jf,1-obj.SS_theta) * obj.tau;
            end
            function [f,Jf] = fun_corr(X)
                f  = obj.ChemSys.propensities(X);
                f  = X - Ys - obj.ChemSys.nu * ( ( obj.eta2 * f + (1-obj.eta2) * a_s ) .* obj.SS_theta ) * obj.tau;
                Jf = obj.ChemSys.propJacobian(X);
                Jf = obj.J_X - obj.ChemSys.nu * bsxfun(@times,Jf,obj.SS_theta) * obj.tau;
            end
        end
        
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                    update mean and covariance of the linear system                      %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ mu, sigma ] = update_mean_covariance(obj,theta,eta1,eta2,R2theta,mean_old,covariance_old,X)
            if nargin < 8
                X = mean_old;
            end
            nu = obj.ChemSys.nu;
            Ja = obj.ChemSys.propJacobian(X);
            d  = obj.ChemSys.propensities(X) - Ja * X;

            tht = zeros(1,obj.ChemSys.num_reactions);
            et1 = zeros(1,obj.ChemSys.num_reactions);
            et2 = zeros(1,obj.ChemSys.num_reactions);
            for i = 1:numel(R2theta)
                for j = R2theta{i}
                    tht(j) = theta(i);
                    et1(j) = eta1(i);
                    et2(j) = eta2(i);
                end
            end

            R1    = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
            r2    = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( nu * diag(1-tht) * d );
            R3    = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
            r4    = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( nu * diag(tht) * d );
            R3_R1 = R3 * R1;
            R3_nu = R3 * nu;
            
            mu    = R3_R1 * mean_old + obj.tau * ( R3*r2 + r4 );
            sigma = R3_R1 * covariance_old * R3_R1' + obj.tau * R3_nu * diag(Ja*R1*mean_old + d + obj.tau*Ja*r2) * R3_nu';
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                frobenius norm of the covariance error and its gradient                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ f, df ] = covariance_norm(obj,theta,eta1,eta2,R2theta,mean_old,covariance_old,mu_ref,covariance_ref,mask)
            Ja = obj.ChemSys.propJacobian(obj.Y(:,obj.curr_ind-1));
            nu = obj.ChemSys.nu;

            tht = zeros(1,obj.ChemSys.num_reactions);
            et1 = zeros(1,obj.ChemSys.num_reactions);
            et2 = zeros(1,obj.ChemSys.num_reactions);
            for i = 1:numel(R2theta)
                for j = R2theta{i}
                    tht(j) = theta(i);
                    et1(j) = eta1(i);
                    et2(j) = eta2(i);
                end
            end

            % update mean and covariance
            R01 = inv( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja );
            R03 = inv( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja );
            R1  = R01 * ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
            R3  = R03 * ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
            A   = R3 * R1;
            A_t = A';

            R3_nu   = R3 * nu;
            R3_nu_t = R3_nu';
            mu      = A * mean_old;
            sigma   = A * covariance_old * A_t + obj.tau * R3_nu * diag(Ja*R1*mean_old) * R3_nu_t;
            
            dsigma = sigma - covariance_ref;
            % logical indexing
            if nargin == 10
                dsigma(~mask) = 0;
            end
            
            f  = norm( dsigma, 'fro' )^2 + norm( mu - mu_ref, 'fro' )^2;
            df = zeros(3*numel(R2theta),1);
            
            R01_nu = R01 * nu;
            R03_nu = R03 * nu;
            for j = 1:numel(R2theta)
                i = R2theta{j}(1);
                
                % theta derivative
                R01_nu_Ja = R01_nu(:,i) * Ja(i,:);
                R03_nu_Ja = R03_nu(:,i) * Ja(i,:);
                
                dR1   = -obj.tau * R01_nu_Ja * ( et1(i)*R1 + (1-et1(i))*obj.J_X );
                dR3   =  obj.tau * R03_nu_Ja * ( et2(i)*R3 + (1-et2(i))*obj.J_X );
                dR3R1 =  dR3*R1 + R3*dR1;
                
                omega = dR3R1 * covariance_old * A_t + obj.tau * dR3 * nu * diag(Ja*R1*mean_old) * R3_nu_t;
                dcov  = omega + omega' + obj.tau * R3_nu * diag(Ja*dR1*mean_old) * R3_nu_t;
                if nargin == 10
                    dcov(~mask) = 0;
                end
                df(j) = 2 * trace(dsigma*dcov);
                
                % eta1 derivative
                dR1   = -obj.tau * (1-tht(i)) * R01_nu_Ja * ( R1 - obj.J_X );
                dR3R1 =  R3*dR1;
                
                omega = dR3R1 * covariance_old * A_t;
                dcov  = omega + omega' + obj.tau * R3_nu * diag(Ja*dR1*mean_old) * R3_nu_t;
                if nargin == 10
                    dcov(~mask) = 0;
                end
                df(numel(R2theta)+j) = 2 * trace(dsigma*dcov);
                
                % eta2 derivative
                dR3   = obj.tau * tht(i) * R03_nu_Ja * ( R3 - obj.J_X );
                dR3R1 = dR3*R1;
                
                omega = dR3R1 * covariance_old * A_t + obj.tau * dR3 * nu * diag(Ja*R1*mean_old) * R3_nu_t;
                dcov  = omega + omega';
                if nargin == 10
                    dcov(~mask) = 0;
                end
                df(2*numel(R2theta)+j) = 2 * trace(dsigma*dcov);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                frobenius norm of the covariance error and its gradient                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [ f ] = Lyapunov_error_norm(obj,theta,eta1,eta2,R2theta,mu_X,cov_X)
            if nargin < 7
                mu_X  = zeros(obj.ChemSys.num_species,1);
                cov_X = eye(obj.ChemSys.num_species);
            end
            nu = obj.ChemSys.nu;
            Ja = obj.ChemSys.propJacobian(mu_X);
            d  = obj.ChemSys.propensities(mu_X) - Ja * mu_X; % zeros(obj.ChemSys.num_reactions,1);
            theta0 = 1;
            
            tht = zeros(1,obj.ChemSys.num_reactions);
            et1 = zeros(1,obj.ChemSys.num_reactions);
            et2 = zeros(1,obj.ChemSys.num_reactions);
            for i = 1:numel(R2theta)
                for j = R2theta{i}
                    tht(j) = theta(i);
                    et1(j) = eta1(i);
                    et2(j) = eta2(i);
                end
            end

            R1 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
            r2 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( nu * diag(1-tht) * d );
            R3 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
            r4 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( nu * diag(tht) * d );
            
            P1 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ ( obj.J_X + obj.tau * (1-theta0) * nu * Ja );
            p2 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ (nu*d);
            P3 = 0.5*obj.J_X - obj.tau *    theta0  * nu * Ja;
            P4 = 0.5*obj.J_X + obj.tau * (1-theta0) * nu * Ja;
            
            R3_R1   = R3 * R1;
            kron_P3 = kron(obj.J_X,P3) + kron(P3,obj.J_X);
            kron_P4 = kron(obj.J_X,P4) + kron(P4,obj.J_X);

            dd = diag(Ja*mu_X+d);
            
            dmu  = ( R3_R1 - P1 ) * mu_X + obj.tau * ( R3*r2 + r4 - p2 );
            dcov = ( kron(R3_R1,R3_R1) - kron_P3\kron_P4 ) * cov_X(:) + obj.tau * ( kron(R3,R3)*kron(nu,nu) - kron_P3\kron(nu,nu) ) * dd(:);

            f = norm( dmu )^2 + norm( dcov )^2;
        end
        
        function [ f ] = Lyapunov_local_error_norm(obj,theta,eta1,eta2,R2theta,mu_X,cov_X)
            if nargin < 7
                mu_X  = zeros(obj.ChemSys.num_species,1);
                cov_X = eye(obj.ChemSys.num_species);
            end
            nu = obj.ChemSys.nu;
            Ja = obj.ChemSys.propJacobian(mu_X);
            d  = obj.ChemSys.propensities(mu_X) - Ja * mu_X; % zeros(obj.ChemSys.num_reactions,1);
            theta0 = 1;
            
            tht = zeros(1,obj.ChemSys.num_reactions);
            et1 = zeros(1,obj.ChemSys.num_reactions);
            et2 = zeros(1,obj.ChemSys.num_reactions);
            for i = 1:numel(R2theta)
                for j = R2theta{i}
                    tht(j) = theta(i);
                    et1(j) = eta1(i);
                    et2(j) = eta2(i);
                end
            end

            R1 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
            r2 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( nu * diag(1-tht) * d );
            R3 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
            r4 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( nu * diag(tht) * d );
            
            P1 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ ( obj.J_X + obj.tau * (1-theta0) * nu * Ja );
            p2 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ (nu*d);
            P3 = 0.5*obj.J_X - obj.tau *    theta0  * nu * Ja;
            P4 = 0.5*obj.J_X + obj.tau * (1-theta0) * nu * Ja;
            
            R3_R1       = R3 * R1;
            kron_nu     = kron(nu,nu);
            kron_R3_nu  = kron(R3*nu,R3*nu);
            kron_P3     = kron(obj.J_X,P3) + kron(P3,obj.J_X);
            kron_P4 = kron(obj.J_X,P4) + kron(P4,obj.J_X);
            
            %diag_R1 = diag(sum(Ja*R1,2));
            %diag_P1 = diag(sum(Ja*P1,2));
            %diag_r2 = diag(d+obj.tau*Ja*r2);
            %diag_p2 = diag(d+obj.tau*Ja*p2);

            %e1 = R3_R1 - P1;
            %e2 = R3*r2 + r4 - p2;
            %e3 = kron(R3_R1,R3_R1)       - kron_P3 \ ( kron(obj.J_X,P4) + kron(P4,obj.J_X) );
            %e4 = kron_R3_nu * diag_R1(:) - kron_P3 \ ( kron_nu * diag_P1(:) );
            %e5 = kron_R3_nu * diag_r2(:) - kron_P3 \ ( kron_nu * diag_p2(:) );

            %f = norm( e1, 'fro')^2 + norm( e3, 'fro')^2 + obj.tau^2 * ( norm( e2, 'fro')^2 + norm( e4, 'fro')^2 + norm( e5, 'fro')^2 );
            
            e1 = R3_R1 - P1;
            e2 = R3*r2 + r4 - p2;
            e3 = kron(R3_R1,R3_R1) - kron_P3\kron_P4;
            e4 = kron(R3,R3) - inv(kron_P3);

            f = norm( e1, 'fro')^2 + norm( e2, 'fro')^2 + norm( e3, 'fro')^2 + norm( e4, 'fro')^2;
        end
%         function [ f ] = Lyapunov_error_norm(obj,theta,eta1,eta2,R2theta,mean_old,covariance_old,mu_ref,covariance_ref)
%             Ja = obj.ChemSys.propJacobian(obj.Y(:,obj.curr_ind-1));
%             nu = obj.ChemSys.nu;
% 
%             d = zeros(obj.ChemSys.num_reactions,1);
%             theta0 = 1;
%             
%             tht = zeros(1,obj.ChemSys.num_reactions);
%             et1 = zeros(1,obj.ChemSys.num_reactions);
%             et2 = zeros(1,obj.ChemSys.num_reactions);
%             for i = 1:numel(R2theta)
%                 for j = R2theta{i}
%                     tht(j) = theta(i);
%                     et1(j) = eta1(i);
%                     et2(j) = eta2(i);
%                 end
%             end
% 
%             % update mean and covariance
%             R1 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
%             r2 = ( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja ) \ ( nu * diag(1-tht) * d );
%             R3 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
%             r4 = ( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja ) \ ( nu * diag(tht) * d );
%             
%             P1 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ ( obj.J_X + obj.tau * (1-theta0) * nu * Ja );
%             p2 = ( obj.J_X - obj.tau * theta0 * nu * Ja ) \ (nu*d);
%             P3 = 0.5*obj.J_X - obj.tau * theta0 * nu * Ja ;
%             P4 = 0.5*obj.J_X + obj.tau * (1-theta0) * nu * Ja ;
%             
%             R3_R1 = R3 * R1;
%             R3_nu = R3 * nu;
% 
%             e1 = R3_R1 - P1;
%             e2 = R3*r2 + r4 - p2;
%             e3 = P3 * ( R3_R1 * R3_R1' ) - P4;
%             e3 = e3 + e3';
%             e4 = P3 * R3_nu * diag(sum(Ja*R1,2)) * R3_nu';
%             e4 = e4 + e4' - nu * diag(sum(Ja*P1,2)) * nu';
%             %e4 = P3 * R3_nu * diag(Ja*R1*mean_old) * R3_nu';
%             %e4 = e4 + e4' - nu * diag(Ja*P1*mean_old) * nu';
%             e5 = P3 * R3_nu * diag(d+obj.tau*Ja*r2) * R3_nu';
%             e5 = e5 + e5' - nu * diag(d) * nu';
% 
%             f = norm( e1, 'fro')^2 + norm( e2, 'fro')^2 + norm( e3, 'fro')^2 + norm( e4, 'fro')^2 + norm( e5, 'fro')^2;
%         end
%         function [ f ] = Lyapunov_error_norm(obj,theta,eta1,eta2,R2theta,mean_old,covariance_old,mu_ref,covariance_ref)
%             Ja = obj.ChemSys.propJacobian(obj.Y(:,obj.curr_ind-1));
%             nu = obj.ChemSys.nu;
% 
%             tht = zeros(1,obj.ChemSys.num_reactions);
%             et1 = zeros(1,obj.ChemSys.num_reactions);
%             et2 = zeros(1,obj.ChemSys.num_reactions);
%             for i = 1:numel(R2theta)
%                 for j = R2theta{i}
%                     tht(j) = theta(i);
%                     et1(j) = eta1(i);
%                     et2(j) = eta2(i);
%                 end
%             end
% 
%             % update mean and covariance
%             R01 = inv( obj.J_X - obj.tau * nu * diag(et1) * diag(1-tht) * Ja );
%             R03 = inv( obj.J_X - obj.tau * nu * diag(et2) * diag(tht)   * Ja );
%             R1  = R01 * ( obj.J_X + obj.tau * nu * diag(1-et1) * diag(1-tht) * Ja );
%             R3  = R03 * ( obj.J_X + obj.tau * nu * diag(1-et2) * diag(tht)   * Ja );
%             P1  = inv( obj.J_X - obj.tau * nu * Ja );
%             A   = R3 * R1;
% 
%             R1_nu   = R1 * nu;
%             R1_nu_t = R1_nu';
%             
%             dmu    = A - P1;
%             dsigma = nu * diag(Ja*R1*mean_old) * nu' - R1_nu * diag(Ja*P1*mean_old) * R1_nu_t;
% 
% %             f  = norm( dsigma, 'fro' )^2 + norm( dmu, 'fro' )^2;
%                         
%             P2 = (0.5*obj.J_X - nu*Ja);
%             
%             asd1 = diag( Ja*P1*mean_old );
% %             inv( kron(obj.J_X,P2) + kron(P2,obj.J_X) ) * covariance_old(:) 
% %             inv( kron(obj.J_X,P2) + kron(P2,obj.J_X) ) * kron(nu,nu) * asd1(:)
%             
%             asd2 = diag( Ja*R1*mean_old );
% %             kron(R3,R3)*kron(R1,R1) * covariance_old(:) 
% %             kron(R3,R3)*kron(nu,nu) * asd2(:)
% 
% 
% %             asf = kron(nu,nu)
% %             asg = kron(R3,R3)
% %             ash = inv(kron(obj.J_X,P2) + kron(P2,obj.J_X))
% %             asg*asf
% %             ash*asf
% 
% ss1 = covariance_old;
% ss2 = covariance_ref;
% % w = 0.5;
% %             f  = (1-w)*norm( kron(A,A)*ss1(:) - inv( kron(obj.J_X,P2) + kron(P2,obj.J_X) )*ss2(:), 'fro') + ...
% %                  (1-w)*norm( kron(R3*nu,R3*nu)*asd2(:) - inv( kron(obj.J_X,P2) + kron(P2,obj.J_X))*kron(nu,nu)*asd1(:), 'fro') + ...
% %                  w*norm( A*mean_old - P1*mu_ref, 'fro' )^2;
%              
%              f  = norm( kron(A,A)*ss1(:) - inv( kron(obj.J_X,P2) + kron(P2,obj.J_X) )*ss2(:), 'fro') + ...
%                   norm( kron(R3*nu,R3*nu)*asd2(:) - inv( kron(obj.J_X,P2) + kron(P2,obj.J_X))*kron(nu,nu)*asd1(:), 'fro') + ...
%                   norm( A*mean_old - P1*mu_ref, 'fro' )^2;
%         end
    end
    
end

