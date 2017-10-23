clc;
clearvars;

% add all folders to the path variable
run('../pathlist');


% input_data = @lotka_volterra_1_data;  Tfin = 5.0;  dt_min = 1e-2;  dt_max = 1e0;  Npaths = 1000;
% input_data = @virus_kinetics_data;  Tfin = 15.0;  dt_min = 1e-3;  dt_max = 1e0;  Npaths = 6;
input_data = @data_apnum2015_ex3;  Tfin = 1e-3;  dt_min = 1e-8;  dt_max = 1e-5;  Npaths = 6;
% input_data = @isomerization_data;  Tfin = 1;  dt_min = 2^(-12);  dt_max = 2^(-4);  Npaths = 10000;


options = struct('theta',1.0,'plot',0,'plot_ss',0,'Jacobian',true);

% methods = { @ssTauLeap };
methods = { @ThetaTauLeap };

leg = cell(1,length(methods));
for i = 1:length(leg)
	leg{i} = func2str(methods{i});
end




%% sample paths
[nu,prop,Y0] = input_data();
tspan = 0:1e-4:Tfin;

% figure(1)
% [Y_ssa,poissprc] = SSA(prop,nu,Tfin,Y0);
% stairs(poissprc(1,:),Y_ssa(1,:),'-x'); hold on;
% for i = 1:length(methods)
%     Y = methods{i}(prop,nu,tspan,Y0,options);
%     stairs(tspan,Y(1,:),'-','LineWidth',2); hold on;
% end
% hold off;
% xlim([0 1*Tfin]);
% % legend(['SSA' leg]);
% xlabel('time');
% ylabel('X1');
% legend('SSA', '\theta=0.95', '\theta=1' );
% set(gca,'FontSize',15);



%%
[h,err,order] = StrongConvergenceTester(methods,input_data,Tfin,Npaths,dt_min,dt_max);

figure(2)
loglog(h,err,'-o');
legend(leg);
order
