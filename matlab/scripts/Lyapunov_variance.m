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
% cs.add_reaction('D<->A');

% reaction rates
cs.rates = 1e5*[ 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0; 1e0 ];
% cs.rates = 1e4*[ 1e0; 1e0 ];

cs.summary;

%% solver setup

% total number of species
xT = 1000;

% initial state
X0    = zeros(cs.num_species,1);
X0(1) = xT;

Tfin    = 1;
N_steps = Tfin;

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

%% theta

theta = 1;

Ja = cs.propJacobian(X0);
nu = cs.nu;
I  = eye(cs.num_species);

mu    = X0;
sigma = zeros(cs.num_species);

for i = 1:N_steps
    mu = ( I - theta*nu*Ja ) \ ( ( I + (1-theta)*nu*Ja ) * mu );

    A = 0.5*I - theta*nu*Ja;
    Q = ( 0.5*I + (1-theta)*nu*Ja ) * sigma + sigma * ( 0.5*I + (1-theta)*nu*Ja )' + nu * diag(Ja*mu) * nu';
    sigma = lyap( A, -Q );
end


%%
fprintf('\n\nNumber of species: %d\n\n', xT);

fprintf('Theta value: %12.4e\n\n', theta);

fprintf('Mean values:');
fprintf('\n\t 0)'); fprintf(' %8.2f', mu_X(:,end)');
fprintf('\n\t 1)'); fprintf(' %8.2f', mu');
fprintf('\n\n');

fprintf('Mean value errors:');
fprintf('\n\t0)'); fprintf(' %8.2f', abs(mu_X(:,end)'-mu'));
fprintf('\n');

fprintf('\nCovariances:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', sigma(i,:));     fprintf('\t\t|');  
fprintf(' %8.2f', cov_X(i,:,end));  
end
fprintf('\n');  

fprintf('\nCovariance absolute errors:');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma(i,:)-cov_X(i,:,end)));
end
fprintf('\n');  

fprintf('\nCovariance relative errors (%%):');
for i = 1:cs.num_species
fprintf('\n');  
fprintf(' %8.2f', abs(sigma(i,:)-cov_X(i,:,end))./abs(cov_X(i,:,end))*100); 
end
fprintf('\n'); 

fprintf('\nCovariance residual:\n');
fprintf('\t1) %10.2e\n', norm(sigma-cov_X(:,:,end),'fro'));


