clearvars
clear all
clc

% add all folders to the path variable
run('../pathlist');

%% system setup

% Chemical system
cs = cChemicalSystem('kinetic_law','combinatorial');
cs.add_reaction('A<->B');
cs.add_reaction('B<->C');
cs.add_reaction('C<->D');
cs.add_reaction('D<->E');

% cs.add_reaction('A<->B');
% cs.add_reaction('B->C');
% cs.add_reaction('C<->D');
% cs.add_reaction('D->A');

% cs.add_reaction('A + 2 B <-> 3 B');

% reaction rates
% cs.rates = 1e5*[ 1e0; 1e-2; 1e2; 1e-2; 1e0; 1e-1; 1e0; 1e0 ];
cs.rates = 1e7*[ 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0 ];
%cs.rates = 1e5*[ 1e0; 1e0; 1e0; 1e0; 1e0; 1e0 ];
% cs.rates = 1e-1*[ 1e1; 1e0; 1e1; 1e1; 1e0; 0.5e1 ];
% cs.rates = 1e5*[ 1e0; 1e1; 1e0; 1e0; 2e0; 1e1 ];
% cs.rates = 1e4*[ 1e0; 1e0 ];

cs.summary;

%% solver setup

% total number of species
xT = 1e2;

% initial state
X0    = zeros(cs.num_species,1);
X0(1) = xT;

Tfin      = 10;
N_steps   = Tfin;
tau       = 1;
solver_th = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',1,'eta1',1,'eta2',1);
solver_ss = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',-1);

%% analytical mean and covariance

p     = zeros(cs.num_species,1);
p(1)  = 1;
mu_a  = zeros(cs.num_species,N_steps+1);
cov_a = zeros(cs.num_species,cs.num_species,N_steps+1);
mu_a(:,1) = X0;
for i = 2:N_steps+1
    p = expm(tau*cs.A(X0)) * p;
    mu_a(:,i) = xT * p;
	for m = 1:cs.num_species
        for n = 1:cs.num_species
            cov_a(m,n,i) = -xT*p(m)*p(n);
        end
        cov_a(m,m,i) = xT*p(m)*(1-p(m));
    end
end


%% Lyapunov mean and covariance

theta_l = 1;

Ja = cs.propJacobian(X0);
d  = zeros(cs.num_reactions,1);
nu = cs.nu;
I  = eye(cs.num_species);

P1 = ( I - theta_l*tau*nu*Ja ) \ ( I + (1-theta_l)*tau*nu*Ja ) ;
p2 = ( I - theta_l*tau*nu*Ja ) \ ( nu * d );
P3 = 0.5*I - theta_l*tau * nu*Ja;
P4 = 0.5*I + (1-theta_l)*tau * nu*Ja;

mu_l(:,1)    = X0;
cov_l(:,:,1) = zeros(cs.num_species);
for i = 2:N_steps+1
    Q = P4 * cov_l(:,:,i-1) + cov_l(:,:,i-1) * P4' + tau * nu * diag(Ja*mu_l(:,i-1)) * nu';
    
    mu_l(:,i)    = P1 * mu_l(:,i-1) + tau * p2;
    cov_l(:,:,i) = lyap( P3, -Q );
end


%% theta = 1

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

mu1a    = X0;
sigma1a = zeros(cs.num_species);
theta1a = ones(1,cs.num_reactions);
for i = 2:N_steps
    [ mu1a, sigma1a ] = solver_th.update_mean_covariance(theta1a,ones(cs.num_reactions,1),ones(cs.num_reactions,1),R2theta,mu1a,sigma1a);
end
residual1a = norm( sigma1a - cov_a(:,:,end), 'fro' );


%% theta for reversible pairs

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

z   = cs.relaxRates(X0)';
ind = find(z<2);
theta1(:)   = sqrt(2./z) - 1./z;
theta1(ind) = 0.633975 - 0.0566243 * z(ind);

mu1    = X0;
sigma1 = zeros(cs.num_species);
eta12  = ones(cs.num_reactions,1);
for i = 1:N_steps
    [ mu1, sigma1 ] = solver_ss.update_mean_covariance(theta1,eta12,eta12,R2theta,mu1,sigma1);
end
residual1 = norm( sigma1 - cov_a(:,:,i), 'fro' );


%% optimal parameters

N_theta = cs.num_reactions;
R2theta = cell(1,N_theta);
for i = 1:N_theta
    R2theta{i} = [ i ];
end
for i = 1:N_theta
    theta0(i) = theta1(end,R2theta{i}(1));
end
temp0(1:N_theta)             = theta1;
temp0(N_theta+1:2*N_theta)   = 1;
temp0(2*N_theta+1:3*N_theta) = 1;

opt_unc = optimoptions(@fminunc,'Display','off','SpecifyObjectiveGradient',false,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10); %,'TypicalX',temp0);
opt_con = optimoptions(@fmincon,'Algorithm','sqp','SpecifyObjectiveGradient',false,'Display','off','DiffMaxChange',1e-2,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'TypicalX',temp0);
lb  = horzcat( 0.0*ones(1,N_theta), 0.7*ones(1,N_theta), 0.7*ones(1,N_theta) );
ub  = horzcat( 1.0*ones(1,N_theta), 1.1*ones(1,N_theta), 1.1*ones(1,N_theta) );


fun = @(X) solver_ss.Lyapunov_error_norm(X(1:N_theta),X(N_theta+1:2*N_theta),X(2*N_theta+1:3*N_theta),R2theta,mu_a(:,end),cov_a(:,:,end));
[ temp, residual, flag,~,~,grad ] = fminunc( fun, temp0, opt_unc );    
% [ temp, residual(i), flag,~,~,grad ] = fmincon( fun, temp0, [], [], [], [], lb, ub, [], opt_con );    
residual

% problem = createOptimProblem('fmincon','objective',fun,'x0',temp0,'lb',lb,'ub',ub,'options',opt_con);
% ms = MultiStart('Display','off','StartPointsToRun','all','XTolerance',1e-10,'UseParallel',true);
% [ temp, residual(i), flag, ms_output, ms_solutions ] = run(ms,problem,4000);
%gs = GlobalSearch('Display','off','StartPointsToRun','all','NumStageOnePoints',2000,'NumTrialPoints',4000,'XTolerance',1e-6);
%[ temp, residual(i), flag ] = run(gs,problem);
% theta2 = temp(1:N_theta);

theta2 = temp(1:N_theta);
eta1_2 = temp(N_theta+1:2*N_theta);
eta2_2 = temp(2*N_theta+1:3*N_theta);

mu2    = X0;
sigma2 = zeros(cs.num_species);
for i = 1:N_steps
    [ mu2, sigma2 ] = solver_ss.update_mean_covariance(temp(1:N_theta),temp(N_theta+1:2*N_theta),temp(2*N_theta+1:3*N_theta),R2theta,mu2,sigma2);
end
%residual2 = sqrt(residual(end));
residual2 = norm( sigma2 - cov_a(:,:,end), 'fro' );


%%
fprintf('\n\nNumber of species: %d\n\n', xT);

for i = 1:numel(R2theta)
    fprintf('R2theta{%d}: ',i);
    fprintf(' %d ',R2theta{i});
    fprintf('\n');
end
fprintf('\n');

fprintf('Theta values:');
fprintf('\n\t1)'); fprintf(' %12.4e', theta0);
fprintf('\n\t2)'); fprintf(' %12.4e', theta2);
fprintf('\n\n');

fprintf('Eta1 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', temp(N_theta+1:2*N_theta));
fprintf('\n\n');

fprintf('Eta2 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', temp(2*N_theta+1:3*N_theta));
fprintf('\n\n');

fprintf('Mean values:');
fprintf('\n\t 0)'); fprintf(' %8.2f', mu_a(:,end)');  fprintf('\t|');  fprintf(' %8.2f', mu_l(:,end)');
fprintf('\n\t1a)'); fprintf(' %8.2f', mu1a(:,end)');
fprintf('\n\t 1)'); fprintf(' %8.2f', mu1(:,end)');
fprintf('\n\t 2)'); fprintf(' %8.2f', mu2(:,end)');
fprintf('\n\n');

fprintf('Mean value errors:');
fprintf('\n\t1a)'); fprintf(' %8.2f', abs(mu_a(:,end)'-mu1a(:,end)'));
fprintf('\n\t 1)'); fprintf(' %8.2f', abs(mu_a(:,end)'-mu1(:,end)'));
fprintf('\n\t 2)'); fprintf(' %8.2f', abs(mu_a(:,end)'-mu2(:,end)'));
fprintf('\n');

fprintf('\nCovariances:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', sigma1a(i,:));    fprintf('\t|');  
fprintf(' %8.2f', sigma1(i,:));     fprintf('\t|');  
fprintf(' %8.2f', sigma2(i,:));     fprintf('\t|');  
fprintf(' %8.2f', cov_a(i,:,end));  fprintf('\t|');  fprintf(' %8.2f', cov_l(i,:,end));
end
fprintf('\n');  

fprintf('\nCovariance absolute errors:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_a(i,:,end)));  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_a(i,:,end)));  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_a(i,:,end)));
end
fprintf('\n');  

fprintf('\nCovariance relative errors (%%):');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_a(i,:,end))./abs(cov_a(i,:,end))*100);  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_a(i,:,end))./abs(cov_a(i,:,end))*100);  fprintf('\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_a(i,:,end))./abs(cov_a(i,:,end))*100);
end
fprintf('\n'); 

fprintf('\nCovariance residuals:\n');
fprintf('\t1a) %10.2e\n', residual1a);
fprintf('\t 1) %10.2e\n', residual1);
fprintf('\t 2) %10.2e\n', residual2);

