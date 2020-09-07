clearvars
clear all
clc

% add all folders to the path variable
run('../pathlist');

%% system setup

% Chemical system
cs = cChemicalSystem('kinetic_law','combinatorial');
cs.add_reaction('A + B <-> C');
cs.add_reaction('A + C <-> B');
cs.add_reaction('B + C <-> A');

% reaction rates
% cs.rates = 1e5*[ 1e0; 1e-2; 1e2; 1e-2; 1e0; 1e-1; 1e0; 1e0 ];
cs.rates = 1e0*[ 1e3; 1e3; 1e-5; 1e1; 1e0; 1e6 ];
% cs.rates = 1e0*[ 1e1; 1e3; 1e1; 1e3; 1e1; 1e3 ];
% cs.rates = 1e-1*[ 1e0; 1e0; 1e0; 1e0; 1e0; 1e0 ];
% cs.rates = 1e-1*[ 1e1; 1e0; 1e1; 1e1; 1e0; 0.5e1 ];
% cs.rates = 1e5*[ 1e0; 1e1; 1e0; 1e0; 2e0; 1e1 ];
% cs.rates = 1e4*[ 1e0; 1e0 ];

cs.summary;

%% solver setup

% initial state
X0 = 1e0*[ 1e3; 1e3; 1e6 ];
% X0 = [ 990; 1010; 1e6 ];
% X0 = [ 1e3; 0; 0 ];

% total number of species
xT = sum(X0);

Tfin      = 1e-2;
N_steps   = 10;
tau       = Tfin / N_steps;
solver_th = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',1,'eta1',1,'eta2',1);
solver_ss = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',-1);

%% Lyapunov mean and covariance

theta_l = 1;
I = eye(cs.num_species);

opt_fsl = optimoptions('fsolve','Display','off');

% mean
mu_l(:,1) = X0;
for i = 2:N_steps+1
    fun = @(X) X - tau * theta_l * cs.nu * cs.propensities(X) - mu_l(:,i-1) - tau * (1-theta_l) * cs.nu * cs.propensities(mu_l(:,i-1));
    mu_l(:,i) = fsolve(fun, mu_l(:,i-1), opt_fsl);
end

% covariance
cov_l(:,:,1) = zeros(cs.num_species);
for i = 2:N_steps+1
    Ja = cs.propJacobian(mu_l(:,i));
    P3 = 0.5*I -     theta_l * tau * cs.nu * Ja;
    P4 = 0.5*I + (1-theta_l) * tau * cs.nu * Ja;

    Q = P4 * cov_l(:,:,i-1) + cov_l(:,:,i-1) * P4' + tau * cs.nu * diag(cs.propensities(mu_l(:,i))) * cs.nu';
    cov_l(:,:,i) = lyap( P3, -Q );
end


%% theta = 1

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

mu1a    = X0;
sigma1a = zeros(cs.num_species);
tht_eta = ones(cs.num_reactions,1);
for i = 2:N_steps+1
    [ mu1a, sigma1a ] = solver_th.update_mean_covariance(tht_eta,tht_eta,tht_eta,R2theta,mu1a,sigma1a,mu_l(:,i));
end
residual1a = norm( sigma1a - cov_l(:,:,end), 'fro' );


%% theta for reversible pairs

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

z   = cs.relaxRates(mu_l(:,end))';
ind = find(z<2);
theta1      = sqrt(2./z) - 1./z;
theta1(ind) = 0.633975 - 0.0566243 * z(ind);

mu1    = X0;
sigma1 = zeros(cs.num_species);
eta12  = ones(cs.num_reactions,1);
for i = 2:N_steps+1
    [ mu1, sigma1 ] = solver_ss.update_mean_covariance(theta1,eta12,eta12,R2theta,mu1,sigma1,mu_l(:,i));
end
residual1 = norm( sigma1 - cov_l(:,:,i), 'fro' );


%% optimal parameters

N_theta = cs.num_reactions;
R2theta = cell(1,N_theta);
for i = 1:N_theta
    R2theta{i} = i;
end
temp0(1:N_theta)             = theta1;
temp0(N_theta+1:2*N_theta)   = 1;
temp0(2*N_theta+1:3*N_theta) = 1;

opt_unc = optimoptions(@fminunc,'Display','off','SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10); %,'TypicalX',temp0);
opt_con = optimoptions(@fmincon,'Algorithm','sqp','SpecifyObjectiveGradient',false,'Display','off','DiffMaxChange',1e-2,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'TypicalX',temp0);
lb  = horzcat( 0.0*ones(1,N_theta), 0.7*ones(1,N_theta), 0.7*ones(1,N_theta) );
ub  = horzcat( 0.1*ones(1,N_theta), 1.2*ones(1,N_theta), 1.2*ones(1,N_theta) );

fun = @(X) solver_ss.Lyapunov_error_norm(X(1:N_theta),X(N_theta+1:2*N_theta),X(2*N_theta+1:3*N_theta),R2theta,mu_l(:,end),cov_l(:,:,end));
% [ temp, residual, flag,~,~,grad ] = fminunc( fun, temp0, opt_unc );    
[ temp, residual, flag,~,~,grad ] = fmincon( fun, temp0, [], [], [], [], lb, ub, [], opt_con );    

% problem = createOptimProblem('fmincon','objective',fun,'x0',temp0,'lb',lb,'ub',ub,'options',opt_con);
% ms = MultiStart('Display','off','StartPointsToRun','all','XTolerance',1e-10,'UseParallel',true);
% [ temp, residual, flag, ms_output, ms_solutions ] = run(ms,problem,4000);
%gs = GlobalSearch('Display','off','StartPointsToRun','all','NumStageOnePoints',2000,'NumTrialPoints',4000,'XTolerance',1e-6);
%[ temp, residual(i), flag ] = run(gs,problem);
theta2 = temp(1:N_theta);
eta1_2 = temp(N_theta+1:2*N_theta);
eta2_2 = temp(2*N_theta+1:3*N_theta);

residual
flag

mu2    = X0;
sigma2 = zeros(cs.num_species);
for i = 2:N_steps+1
    [ mu2, sigma2 ] = solver_ss.update_mean_covariance(theta2,eta1_2,eta2_2,R2theta,mu2,sigma2,mu_l(:,i));
end
residual2 = norm( sigma2 - cov_l(:,:,end), 'fro' );

%%
fprintf('\n\nNumber of species: %d\n\n', xT);

for i = 1:numel(R2theta)
    fprintf('R2theta{%d}: ',i);
    fprintf(' %d ',R2theta{i});
    fprintf('\n');
end
fprintf('\n');

fprintf('Theta values:');
fprintf('\n\t1)'); fprintf(' %12.4e', theta1);
fprintf('\n\t2)'); fprintf(' %12.4e', theta2);
fprintf('\n\n');

fprintf('Eta1 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', eta1_2);
fprintf('\n\n');

fprintf('Eta2 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', eta2_2);
fprintf('\n\n');

fprintf('Mean values:');
fprintf('\n\t 0)'); fprintf(' %8.2f', mu_l(:,end)');  fprintf('\t|');  fprintf(' %8.2f', mu_l(:,end)');
fprintf('\n\t1a)'); fprintf(' %8.2f', mu1a(:,end)');
fprintf('\n\t 1)'); fprintf(' %8.2f', mu1(:,end)');
fprintf('\n\t 2)'); fprintf(' %8.2f', mu2(:,end)');
fprintf('\n\n');

fprintf('Mean value errors:');
fprintf('\n\t1a)'); fprintf(' %8.2f', abs(mu_l(:,end)'-mu1a(:,end)'));
fprintf('\n\t 1)'); fprintf(' %8.2f', abs(mu_l(:,end)'-mu1(:,end)'));
fprintf('\n\t 2)'); fprintf(' %8.2f', abs(mu_l(:,end)'-mu2(:,end)'));
fprintf('\n');

fprintf('\nCovariances:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', sigma1a(i,:));    fprintf('\t|');  
fprintf(' %8.2f', sigma1(i,:));     fprintf('\t|');  
fprintf(' %8.2f', sigma2(i,:));     fprintf('\t|');  
fprintf(' %8.2f', cov_l(i,:,end));
end
fprintf('\n');  

fprintf('\nCovariance absolute errors:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_l(i,:,end)));  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_l(i,:,end)));  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_l(i,:,end)));
end
fprintf('\n');  

fprintf('\nCovariance relative errors (%%):');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_l(i,:,end))./abs(cov_l(i,:,end))*100);  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_l(i,:,end))./abs(cov_l(i,:,end))*100);  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_l(i,:,end))./abs(cov_l(i,:,end))*100);
end
fprintf('\n'); 

fprintf('\nCovariance residuals:\n');
fprintf('\t1a) %10.2e\n', residual1a);
fprintf('\t 1) %10.2e\n', residual1);
fprintf('\t 2) %10.2e\n', residual2);

