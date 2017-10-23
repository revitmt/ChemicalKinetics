% generate ensembles of solutions for a given set of reactions
% test weak convergence of numerical schemes
% test evolution of distributions, mean values and variances

clc;
clearvars;

% add all folders to the path variable
run('../pathlist');


set(groot,'defaultAxesFontSize',15);
set(groot,'defaultLineLineWidth',2,'defaultLineMarkerSize',7,'defaultLineMarkerEdgeColor','black');


% input_data = @lotka_volterra_1_data;  Tfin = 5.0;  dt_min = 1e-2;  dt_max = 1e0;  Npaths = 1000;
% input_data = @virus_kinetics_data;  Tfin = 15.0;  dt_min = 1e-3;  dt_max = 1e0;  Npaths = 6;
input_data = @data_apnum2015_ex3;  Tfin = 1e-5;  h = [ 1e-6 5e-6 ]; Npaths = 1000;
% input_data = @isomerization_data;  Tfin = 100;  dt_min = 2^(-2);  dt_max = 2^(-1);  Npaths = 1000;
% input_data = @dimerization_data;  Tfin = 1e-1; h = 1e-2;
% input_data = @reversible_isomerization_data;  Tfin = 1e-1;  h = [ 1e-3 1e-2 5e-2 ];
% input_data = @fast_slow_1_data;  Tfin = 5e-2; h = [ 1e-3 5e-3 1e-2 ];
% input_data = @data_fast_slow_2;  Tfin = 5e-3; h = [ 5e-4 ];
% input_data = @Schlogl_data;  Tfin = 4; h = 1e0;

% methods = { @ssTauLeap };
methods = { @ThetaTauLeap };

solver_options = {{ struct('Jacobian',true) }};

options = struct('test_SSA',1,'test_methods',1,'NPaths',Npaths,'solver_options',solver_options);

leg = cell(1,length(methods));
for i = 1:length(leg)
	leg{i} = func2str(methods{i});
    %leg{i} = sprintf('theta = %0.1f',solver_options{1}{i}.theta);
end

Nmeth = length(methods);
Nlev  = length(h);


%%
[h,err,order,tspan,Y_mean,Y_var,Y,SSA_tspan,SSA_mean,SSA_var,SSA_Y] = WeakConvergenceTester(methods,input_data,Tfin,h,options);

N = size(Y_mean{1,1},1);

figure(2)
loglog(h,err,'-o');
legend(leg);
order


% plot means of solutions
figure(1);
for lev = 1:Nlev
    subplot(Nlev,1,lev);
    for method = 1:Nmeth
        for i = 1:N
            plot(tspan{lev},Y_mean{method,lev}(i,:),'-o'); hold on;
        end
    end
    if lev == 1
        title('Means of solutions');
    end
    legend(leg);
    xlabel('time');
    ylabel('Y');
end
hold off;


% plot variances of solutions
figure(2);
for lev = 1:Nlev
    subplot(Nlev,1,lev);
    for method = 1:Nmeth
        for i = 1:N
            plot(tspan{lev},Y_var{method,lev}(i,:),'-o'); hold on;
        end
    end
    if lev == 1
        title('Variances of solutions');
    end
    legend(leg);
    xlabel('time');
    ylabel('Y');
end
hold off;


% dist_supp = 1:500;
% 
% % plot distributions for different step sizes
% [XX,YY] = meshgrid(dist_supp,h);
% f = zeros(Nlev,length(dist_supp));
% figure(3);
% for method = 1:Nmeth
%     for i = 1:1%N
%         for lev = 1:Nlev
%             distribution = squeeze(Y{method,lev}(i,end,:));
%             distribution(distribution<=0) = 1e-5;
%             f(lev,:) = ksdensity(distribution,dist_supp,'kernel','triangle','Support','positive');
%             %plot(dist_supp,f(lev,:),'LineWidth',2); hold on;
%         end
%         mesh(XX,YY,f,'MeshStyle','row','FaceAlpha',0);hold on;
%         colormap([0,0,0]);
%     end
% end
% 
% % plot distributions for different values of theta
% [XX,YY] = meshgrid(dist_supp,0.0:0.1:1.0);
% f = zeros(11,length(dist_supp));
% figure(4);
% for lev = 1:Nlev
%     for i = 1:1%N
%         for method = 1:Nmeth
%             distribution = squeeze(Y{method,lev}(i,end,:));
%             distribution(distribution<=0) = 1e-5;
%             f(method,:) = ksdensity(distribution,dist_supp,'kernel','triangle','Support','positive');
%             %plot(dist_supp,f(lev,:),'LineWidth',2); hold on;
%         end
%         mesh(XX,YY,f,zeros(size(f))+lev,'MeshStyle','row','FaceAlpha',0,'LineWidth',1.0);hold on;
%         %colormap([0,0,0]);
%     end
%     ylabel('theta');
% end
% leg_h = cell(1,Nlev);
% for i = 1:length(leg_h)
%     leg_h{i} = sprintf('h = %0.4f',h(i));
% end
% legend(leg_h);
% 
% 
% 
% histogram(squeeze(Y{1,1}(1,end,:)),'BinMethod','integers')
% 
% 
% for i = 1:size(SSA_mean,1)
%     plot(SSA_tspan,SSA_mean(i,:));hold on;
% end
% hold off;
% 
% histogram(squeeze(SSA_Y(3,end,:)),'normalization','pdf','BinMethod','integers'); hold on;
% histogram(squeeze(Y{1,1}(3,end,:)),'normalization','pdf','BinMethod','integers'); hold on;
% histogram(squeeze(Y{2,1}(1,end,:)),'normalization','pdf','BinMethod','integers'); hold on;
% histogram(squeeze(Y{3,1}(1,end,:)),'normalization','pdf','BinMethod','integers'); hold off;
% % histogram(squeeze(SSA_Y(2,end,:)),'normalization','pdf','BinMethod','integers'); hold off;

