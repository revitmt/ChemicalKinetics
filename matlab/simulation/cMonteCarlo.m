classdef cMonteCarlo < handle
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
%   Monte Carlo method
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    properties ( GetAccess='public', SetAccess='private', Hidden=false )
        solver;     % solver for the model to be simulated
        Nsamples;   % number of MC samples
        fun;        % functional of the solution to be simulated
        stat;       % struct with generated statistic
        stat1;      % struct with generated statistic for path 1
        stat2;      % struct with generated statistic for path 2
    end
    properties ( GetAccess='public', SetAccess='private', Hidden=false )
        stat_data;
        stat_data1;
        stat_data2;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              object constructor                                         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = cMonteCarlo(varargin)
            for i = 1:2:nargin
                option = varargin{i};
                switch option
                    case 'solver'
                        obj.solver = varargin{i+1};
                    case 'functional'
                        obj.fun = varargin{i+1};
                end
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                run simulation                                           %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function run(obj,Nsamples)
            obj.Nsamples      = Nsamples;
            obj.stat.Nsamples = Nsamples;
            
            f1  = zeros(Nsamples,1);
            t1  = zeros(Nsamples,1);

            if isa(obj.solver,'cCoupledPaths')
                f2 = zeros(Nsamples,1);    
                t2 = zeros(Nsamples,1);
                parfor i = 1:obj.Nsamples
                    obj.solver.generate();
                    f1(i) = obj.fun( obj.solver.Yfin{1} );
                    f2(i) = obj.fun( obj.solver.Yfin{2} );
                    t1(i) = obj.solver.path1.info.time;
                    t2(i) = obj.solver.path2.info.time;
                end
                cc = corrcoef(f1,f2);
                obj.stat.corr = cc(1,2);       
                %obj.stat1.mu = mean(f1);     obj.stat1.sigma = var(f1);     obj.stat.cost1 = mean(t1);
                %obj.stat2.mu = mean(f2);     obj.stat2.sigma = var(f2);     obj.stat.cost2  = mean(t2);
                %obj.stat.mu  = mean(f2-f1);  obj.stat.sigma  = var(f2-f1);  obj.stat.cost   = mean(t1+t2);   
                obj.stat1.cost = mean(t1);
                obj.stat2.cost = mean(t2);
                obj.stat.cost  = mean(t1+t2);   
            else
                parfor i = 1:obj.Nsamples
                    obj.solver.generate();
                    f1(i) = obj.fun( obj.solver.Yfin );
                    t1(i) = obj.solver.info.time;
                end
                %obj.stat.mu    = mean(f1);
                %obj.stat.sigma = var(f1);
                obj.stat.cost  = mean(t1);
            end
            
            
            obj.stat_data.M1 = 0;  obj.stat_data1.M1 = 0;  obj.stat_data2.M1 = 0;
            obj.stat_data.M2 = 0;  obj.stat_data1.M2 = 0;  obj.stat_data2.M2 = 0;
            obj.stat_data.M3 = 0;  obj.stat_data1.M3 = 0;  obj.stat_data2.M3 = 0;
            obj.stat_data.M4 = 0;  obj.stat_data1.M4 = 0;  obj.stat_data2.M4 = 0;
            if isa(obj.solver,'cCoupledPaths')
                for i = 1:obj.Nsamples
                    temp1 = (i-1)*(i-2);
                    temp2 = (i-1)*(i^2-3*i+3);
                    
                    f = f2(i) - f1(i);
                    delta   = f - obj.stat_data.M1;
                    delta_i = delta / i;
                    obj.stat_data.M1 = obj.stat_data.M1 + delta_i;
                    M2_new = obj.stat_data.M2 + delta * delta_i * (i-1);
                    M3_new = obj.stat_data.M3 + delta * delta_i^2 * temp1 - 3*delta_i*obj.stat_data.M2;
                    obj.stat_data.M4 = obj.stat_data.M4 + delta*delta_i^3*temp2 + 6*delta_i^2*obj.stat_data.M2 - 4*delta_i*obj.stat_data.M3;
                    obj.stat_data.M3 = M3_new;
                    obj.stat_data.M2 = M2_new;
                    
                    delta   = f1(i) - obj.stat_data1.M1;
                    delta_i = delta / i;
                    obj.stat_data1.M1 = obj.stat_data1.M1 + delta_i;
                    M2_new = obj.stat_data1.M2 + delta * delta_i * (i-1);
                    M3_new = obj.stat_data1.M3 + delta * delta_i^2 * temp1 - 3*delta_i*obj.stat_data1.M2;
                    obj.stat_data1.M4 = obj.stat_data1.M4 + delta*delta_i^3*temp2 + 3*delta_i^2*obj.stat_data1.M2 - 4*delta_i*obj.stat_data1.M3;
                    obj.stat_data1.M3 = M3_new;
                    obj.stat_data1.M2 = M2_new;
                    
                    delta   = f2(i) - obj.stat_data2.M1;
                    delta_i = delta / i;
                    obj.stat_data2.M1 = obj.stat_data2.M1 + delta_i;
                    M2_new = obj.stat_data2.M2 + delta * delta_i * (i-1);
                    M3_new = obj.stat_data2.M3 + delta * delta_i^2 * temp1 - 3*delta_i*obj.stat_data2.M2;
                    obj.stat_data2.M4 = obj.stat_data2.M4 + delta*delta_i^3*temp2 + 3*delta_i^2*obj.stat_data2.M2 - 4*delta_i*obj.stat_data2.M3;
                    obj.stat_data2.M3 = M3_new;
                    obj.stat_data2.M2 = M2_new;
                end
            else
                for i = 1:obj.Nsamples
                    delta   = f1(i) - obj.stat_data.M1;
                    delta_i = delta / i;
                    obj.stat_data.M1 = obj.stat_data.M1 + delta_i;
                    M2_new = obj.stat_data.M2 + delta * delta_i * (i-1);
                    M3_new = obj.stat_data.M3 + delta * delta_i^2 * (i-1)*(i-2) - 3*delta_i*obj.stat_data.M2;
                    obj.stat_data.M4 = obj.stat_data.M4 + delta*delta_i^3*(i-1)*(i^2-3*i+3) + 3*delta_i^2*obj.stat_data.M2 - 4*delta_i*obj.stat_data.M3;
                    obj.stat_data.M3 = M3_new;
                    obj.stat_data.M2 = M2_new;
                end
            end

            
            obj.stat.mu    = obj.stat_data.M1;
            obj.stat.sigma = obj.stat_data.M2 / ( obj.Nsamples - 1 );
            obj.stat.skew  = obj.stat_data.M3 * sqrt( obj.Nsamples / obj.stat_data.M2^3 );
            obj.stat.kurt  = obj.stat_data.M4 * obj.Nsamples / obj.stat_data.M2^2;
            if isa(obj.solver,'cCoupledPaths')
                obj.stat1.mu    = obj.stat_data1.M1;
                obj.stat1.sigma = obj.stat_data1.M2 / ( obj.Nsamples - 1 );
                obj.stat1.skew  = obj.stat_data1.M3 * sqrt( obj.Nsamples / obj.stat_data1.M2^3 );
                obj.stat1.kurt  = obj.stat_data1.M4 * obj.Nsamples / obj.stat_data1.M2^2;
                
                obj.stat2.mu    = obj.stat_data2.M1;
                obj.stat2.sigma = obj.stat_data2.M2 / ( obj.Nsamples - 1 );
                obj.stat2.skew  = obj.stat_data2.M3 * sqrt( obj.Nsamples / obj.stat_data2.M2^3 );
                obj.stat2.kurt  = obj.stat_data2.M4 * obj.Nsamples / obj.stat_data2.M2^2;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                           correlation of coupled paths                                  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function coupled_correlation(obj,Nsamples)
            K = obj.solver.path1.integrators{1}.K + obj.solver.path1.integrators{2}.K - 1;
            
            f1 = zeros(Nsamples,K);
            f2 = zeros(Nsamples,K);
            corr  = zeros(1,K);
            sigma = zeros(1,K);
            mu    = zeros(1,K);

            parfor i = 1:Nsamples
                obj.solver.generate();
                temp1 = zeros(1,obj.solver.path1.integrators{1}.K);
                temp2 = zeros(1,obj.solver.path1.integrators{2}.K-1);
                
                for j = 1:obj.solver.path1.integrators{1}.K
                    temp1(j) = obj.fun( obj.solver.path1.integrators{1}.Y(:,j) );
                end
                for j = 2:obj.solver.path1.integrators{2}.K
                    temp2(j-1) = obj.fun( obj.solver.path1.integrators{2}.Y(:,j) );
                end
                f1(i,:) = [temp1 temp2];
                
                for j = 1:obj.solver.path1.integrators{1}.K
                    temp1(j) = obj.fun( obj.solver.path2.integrators{1}.Y(:,obj.solver.path2.integrators{2}.Q*(j-1)+1) ); 
                end
                for j = 2:obj.solver.path1.integrators{2}.K
                    temp2(j-1) = obj.fun( obj.solver.path2.integrators{2}.Y(:,obj.solver.path2.integrators{2}.Q*(j-1)+1) );
                end
                f2(i,:) = [temp1 temp2];
            end
            for j = 1:K
                cc = corrcoef(f1(:,j),f2(:,j));
                corr(j)  = cc(1,2);
                sigma(j) = var(f1(:,j)-f2(:,j));
                mu(j)    = mean(f1(:,j)-f2(:,j));
            end
            %f2
            %corr
            %hold on
            
            figure(1)
            plot( [ obj.solver.path1.integrators{1}.T obj.solver.path1.integrators{2}.T(2:end) ], corr,'LineWidth',2 );
            xlim( [obj.solver.path1.integrators{1}.T(1) obj.solver.path1.integrators{2}.T(end) ]);
            ylim([-1.2 1.2]);
            title('Correlation');
            
            figure(2)
            plot( [ obj.solver.path1.integrators{1}.T obj.solver.path1.integrators{2}.T(2:end) ], mu,'LineWidth',2 );
            xlim( [obj.solver.path1.integrators{1}.T(1) obj.solver.path1.integrators{2}.T(end) ]);
            title('Mean');
            
            figure(3)
            plot( [ obj.solver.path1.integrators{1}.T obj.solver.path1.integrators{2}.T(2:end) ], sigma,'LineWidth',2 );
            xlim( [obj.solver.path1.integrators{1}.T(1) obj.solver.path1.integrators{2}.T(end) ]);
            title('Variance');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                estimate cost                                            %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function time = estimate_cost(obj,Nsamples)
            tt = tic;
            for i = 1:5
                obj.solver.generate();
            end
            time = toc(tt);
            time = time * Nsamples / 5;
        end
    end
    
end

