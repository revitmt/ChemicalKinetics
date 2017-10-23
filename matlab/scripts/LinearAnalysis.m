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

% reaction rates
% cs.rates = 1e5*[ 1e0; 1e-2; 1e2; 1e-2; 1e0; 1e-1; 1e0; 1e0 ];
cs.rates = 1e3*[ 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0 ];
% cs.rates = 1e4*[ 1e0; 1e0 ];
% cs.rates = 1e4*[ 1e0; 1e0 ];

cs.summary;

%% solver setup

% total number of species
xT = 1000;

% initial state
X0    = zeros(cs.num_species,1);
X0(1) = xT;

Tfin      = 5;
N_steps   = Tfin;
solver_th = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',1,'eta1',1,'eta2',1);
solver_ss = cSplitStepTauLeap('system',cs,'X0',X0,'T0',0,'Tfin',Tfin,'N_steps',N_steps,'theta',-1);

%% analytical mean and covariance

p     = zeros(cs.num_species,1);
p(1)  = 1;
mu_X  = zeros(cs.num_species,N_steps);
cov_X = zeros(cs.num_species,cs.num_species,N_steps);
for i = 1:N_steps
    p = expm(cs.A(X0)) * p;
    mu_X(:,i) = xT * p;
	for m = 1:cs.num_species
        for n = 1:cs.num_species
            cov_X(m,n,i) = -xT*p(m)*p(n);
        end
        cov_X(m,m,i) = xT*p(m)*(1-p(m));
    end
end


%% theta = 1

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

mu1a    = X0;
sigma1a = zeros(cs.num_species);
theta1a = zeros(1,cs.num_reactions);
for i = 1:N_steps
    mu_old  = mu1a;
    cov_old = sigma1a;

    theta1a(:) = 1;
    [ sigma1a, mu1a ] = solver_th.update_mean_covariance(theta1a,ones(cs.num_reactions,1),ones(cs.num_reactions,1),R2theta,mu_old,cov_old);
end
residual1a = norm( sigma1a - cov_X(:,:,i), 'fro' );


%% theta for reversible pairs

R2theta = cell(1,cs.num_reactions);
for i = 1:cs.num_reactions
    R2theta{i} = i;
end

mu1    = X0;
sigma1 = zeros(cs.num_species);
theta1 = zeros(N_steps,cs.num_reactions);
for i = 1:N_steps
    mu_old  = mu1;
    cov_old = sigma1;
    
    z   = cs.relaxRates(X0)';
    ind = find(z<2);
    theta1(i,:)   = sqrt(2./z) - 1./z;
    theta1(i,ind) = 0.633975 - 0.0566243 * z(ind);
    
    [ sigma1, mu1 ] = solver_ss.update_mean_covariance(theta1(i,:),ones(cs.num_reactions,1),ones(cs.num_reactions,1),R2theta,mu_old,cov_old);
end
residual1 = norm( sigma1 - cov_X(:,:,i), 'fro' );


%% optimal parameters

% N_rev   = size(cs.rev_reactions,1);
% N_irrev = size(cs.irrev_reactions,1);
% N_theta = N_rev + N_irrev;
% R2theta = cell(1,N_theta);
% for i = 1:N_rev
%     R2theta{i} = cs.rev_reactions(i,:);
% end
% for i = N_rev+1:N_theta
%     R2theta{i} = cs.irrev_reactions(i-N_rev);
% end

% N_theta = 2;
% R2theta = cell(1,N_theta);
% R2theta{1} = [ 1 3 5 ];
% R2theta{2} = [ 2 4 6 ];
%R2theta{3} = [ 4 ];

% N_theta = 1;
% R2theta = cell(1,N_theta);
% R2theta{1} = [ 1 2 3 4 5 6 ];

N_theta = cs.num_reactions;
R2theta = cell(1,N_theta);
for i = 1:N_theta
    R2theta{i} = [ i ];
end

mask = eye(cs.num_species);
mask = logical(mask);

for i = 1:N_theta
    theta0(i) = theta1(end,R2theta{i}(1));
end
temp0(1:N_theta)             = theta0;
temp0(N_theta+1:2*N_theta)   = 1;
temp0(2*N_theta+1:3*N_theta) = 1;

opt_unc = optimoptions(@fminunc,'Algorithm','trust-region','Display','off','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10); %,'TypicalX',temp0);
opt_con = optimoptions(@fmincon,'Algorithm','sqp','SpecifyObjectiveGradient',false,'Display','off','DiffMaxChange',1e-2,'OptimalityTolerance',1e-10,'StepTolerance',1e-10,'FunctionTolerance',1e-10,'TypicalX',temp0);

b   = 3;
lb  = horzcat( 0.0*ones(1,N_theta), 0.7*ones(1,N_theta), 0.7*ones(1,N_theta) );
ub  = horzcat( 1.0*ones(1,N_theta), 1.0*ones(1,N_theta), 1.0*ones(1,N_theta) );
Aeq = [];
beq = [];
% Aeq = diag( horzcat(zeros(1,N_theta),ones(1,2*N_theta)) );
% beq = vertcat(zeros(N_theta,1),ones(2*N_theta,1));

mu2      = X0;
sigma2   = zeros(cs.num_species);
theta2   = zeros(N_steps,N_theta);
residual = zeros(N_steps,1);
for i = 1:N_steps
    mu_old  = mu2;
    cov_old = sigma2;
    
    %[ temp, residual(i), flag ] = fminunc( fun, temp0, opt_unc );
    %fun = @(X) trace( ( solver.update_mean_covariance(X(1:N_theta),X(N_theta+1:2*N_theta),X(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old) - cov_X(:,:,i) ).^2 );
    %fun = @(X) norm( solver.update_mean_covariance(X(1:N_theta),X(N_theta+1:2*N_theta),X(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old) - cov_X(:,:,i), 'fro' )^2;
    
    %if norm( solver.update_mean_covariance(temp0(1:N_theta),temp0(N_theta+1:2*N_theta),temp0(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old) - cov_old, 'fro' ) > 1e-2
        fun = @(X) solver_ss.covariance_norm(X(1:N_theta),X(N_theta+1:2*N_theta),X(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old,mu_X(:,i),cov_X(:,:,i));%,mask);
        [ temp, residual(i), flag,~,~,grad ] = fmincon( fun, temp0, [], [], Aeq, beq, lb, ub, [], opt_con );    
%         if ( residual(i) >= 1 )
%             problem = createOptimProblem('fmincon','objective',fun,'x0',temp0,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub,'options',opt_con);
%             ms = MultiStart('Display','off','StartPointsToRun','all','XTolerance',1e-10,'UseParallel',true);
%             [ temp, residual(i), flag, ms_output, ms_solutions ] = run(ms,problem,4000);
%             %gs = GlobalSearch('Display','off','StartPointsToRun','all','NumStageOnePoints',2000,'NumTrialPoints',4000,'XTolerance',1e-6);
%             %[ temp, residual(i), flag ] = run(gs,problem);
%         end
    %else
    %    residual(i) = residual(i-1);
    %end
	fprintf('%3d %e ',i, residual(i));
    fprintf(' %e ',temp);
    fprintf('\n');
    temp0 = temp;
    theta2(i,:) = temp(1:N_theta);
    
    [ sigma2, mu2, A ] = solver_ss.update_mean_covariance(temp(1:N_theta),temp(N_theta+1:2*N_theta),temp(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old);
    %fprintf('%e ', eig(A));
    %fprintf('\n');
end
residual2 = sqrt(residual(end));


% [f,df] = solver.covariance_norm(temp(1:N_theta),temp(N_theta+1:2*N_theta),temp(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old,cov_X(:,:,i));
% f 
% df


for i = 1:N_theta
    theta0(i) = theta1(end,R2theta{i}(1));
end
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
fprintf('\n\t2)'); fprintf(' %12.4e', theta2(end,:));
fprintf('\n\n');

fprintf('Eta1 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', temp(N_theta+1:2*N_theta));
fprintf('\n\n');

fprintf('Eta2 values:');
fprintf('\n\t2)'); fprintf(' %12.4e', temp(2*N_theta+1:3*N_theta));
fprintf('\n\n');

fprintf('Mean values:');
fprintf('\n\t 0)'); fprintf(' %8.2f', mu_X(:,end)');
fprintf('\n\t1a)'); fprintf(' %8.2f', mu1a(:,end)');
fprintf('\n\t 1)'); fprintf(' %8.2f', mu1(:,end)');
fprintf('\n\t 2)'); fprintf(' %8.2f', mu2(:,end)');
fprintf('\n\n');

fprintf('Mean value errors:');
fprintf('\n\t0)'); fprintf(' %8.2f', abs(mu_X(:,end)'-mu1a(:,end)'));
fprintf('\n\t0)'); fprintf(' %8.2f', abs(mu_X(:,end)'-mu1(:,end)'));
fprintf('\n\t0)'); fprintf(' %8.2f', abs(mu_X(:,end)'-mu2(:,end)'));
fprintf('\n');

fprintf('\nCovariances:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', sigma1a(i,:));    fprintf('\t\t|');  
fprintf(' %8.2f', sigma1(i,:));     fprintf('\t\t|');  
fprintf(' %8.2f', sigma2(i,:));     fprintf('\t|');  
fprintf(' %8.2f', cov_X(i,:,end));  
end
fprintf('\n');  

fprintf('\nCovariance absolute errors:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_X(i,:,end)));  fprintf('\t\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_X(i,:,end)));  fprintf('\t\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_X(i,:,end)));
end
fprintf('\n');  

fprintf('\nCovariance relative errors (%%):');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma1a(i,:)-cov_X(i,:,end))./abs(cov_X(i,:,end))*100);  fprintf('\t\t|');  
fprintf(' %8.2f', abs(sigma1(i,:)-cov_X(i,:,end))./abs(cov_X(i,:,end))*100);  fprintf('\t\t|');  
fprintf(' %8.2f', abs(sigma2(i,:)-cov_X(i,:,end))./abs(cov_X(i,:,end))*100);
end
fprintf('\n'); 

fprintf('\nCovariance residuals:\n');
fprintf('\t1) %10.2e\n', residual1a);
fprintf('\t1) %10.2e\n', residual1);
fprintf('\t2) %10.2e\n', residual2);


%%

% fun = @(x,y) solver.covariance_norm([x y],temp(N_theta+1:2*N_theta),temp(2*N_theta+1:3*N_theta),R2theta,mu_old,cov_old,cov_X(:,:,end));
% 
% % fsurf(fun,[0 1],'MeshDensity',50);
% % zlim([1e2 1e6]);
% fcontour(fun,[0 1 0 1],'MeshDensity',200,'LevelList',linspace(1e4,5e4,50));
% 
% [f,df] = fun(temp(1),temp(2));
%  rectangle('Position',[ temp(1:2)-0.02 0.04 0.04])
%  line([temp(1) temp(1)+0.1*df(1)/norm(df(1:2))],[temp(2) temp(2)+0.1*df(2)/norm(df(1:2))])



%%

% for i = 1:numel(ms_solutions)
%     plot(i,ms_solutions(i).Fval,'o'); hold on
% end
% zlim([0 1e-4])
% hold off