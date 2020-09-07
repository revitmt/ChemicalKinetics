clearvars
clc

% add all folders to the path variable
run('../pathlist');

%cs = data_apnum2015_ex3('kinetic_law','combinatorial');
% cs = data_three_reversible_isomerizations('kinetic_law','combinatorial');

cs = cChemicalSystem('kinetic_law','combinatorial');

cs.add_reaction('A<->B');
cs.add_reaction('B<->C');
% cs.add_reaction('C->D');

% X = cs.X0;
cs.rates=[1;2;3;4];
X = [5;5;0];

cs.propensities(X);
cs.propJacobian(X);
cs.relaxRates(X);

cs.summary;


% %SSA = cSSA('system',cs,'X0',cs.X0,'T',linspace(0,1e-4,100));
% SSA = cSSA('system',cs,'X0',cs.X0,'Tfin',1e-5);
% 
% SSA.generate()
% 
% SSA
% 
% SSA.plot([1 2])
% %SSA.plot()

Tfin = 1e-3;
coef = 1.1;
solver1_1 = cSplitStepTauLeap('system',cs,'X0',cs.X0,'T0',0,             'Tfin',Tfin,     'N_steps',100,'theta',-1,'int_states',false);
solver2_1 = cSplitStepTauLeap('system',cs,'X0',cs.X0,'T0',0,             'Tfin',Tfin,     'N_steps',200,'theta',-1,'int_states',false);
solver1_2 = cThetaTauLeap('system',cs,'X0',cs.X0,'T0',solver1_1.Tfin,'Tfin',coef*Tfin,'N_steps',200,'theta',0,'int_states',false);
solver2_2 = cThetaTauLeap('system',cs,'X0',cs.X0,'T0',solver1_1.Tfin,'Tfin',coef*Tfin,'N_steps',200,'theta',0,'int_states',false);
%solver1.generate();
%solver1.plot([1 2])

solver_c = cCompositeIntegrator(solver1_1, solver1_2);
solver_f = cCompositeIntegrator(solver2_1, solver2_2);
% solver_c.add_integrator(solver1_2);
% solver_c.integrators{1} = solver1_c;
% solver_c.integrators{2} = solver2_c;

%solver

solver = cCoupledPaths(solver_c,solver_f);
% solver.generate()
% solver.plot([1])

simulation = cMonteCarlo('solver',solver,'functional',@(X)abs(X(1)-X(2)));
simulation.coupled_correlation(100);
% simulation.run(10);
% simulation.estimate_cost(100)
% simulation.stat


% simulation = cMLMC('system',cs,'X0',cs.X0,'Tfin',1e-3,'functional',@(X)abs(X(1)- X(2)),'Nlevels',13);
% simulation.choose_solvers();
% simulation.plot_params();

