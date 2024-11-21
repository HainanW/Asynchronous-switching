clear all 

%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-8;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on','Display', 'iter', ...
    'TolX', NLPtol, 'TolFun', NLPtol, 'TolCon', NLPtol,...
    'MaxFunEvals',10000,'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars

% number of control vars 
m.nu = 2;   
% number of FDs for two controls 
m.Nj = [3,3];
% number of state vars
m.nx = 3;    
% p: total number of time nodes 
m.p = sum(m.Nj) - m.nu; 
m.tF = 1.0;
x0.phi_u1 = [1 0 0 zeros(1, m.nx)]; % 6 = 3 + 3 
x0.phi_u2 = [1 0 0 zeros(1, m.nx)]; % 6 = 3 + 3
   x0.psi = [1 0 0 zeros(1, m.nx)]; % 8 = 3 + 5 
% Decision variable input sequence 
% u1 u2 t1vec t2vec tF 
lb = [-10*ones(1,sum(m.Nj)) zeros(1,sum(m.Nj)-m.nu) 1]; 
ub = [ 10*ones(1,sum(m.Nj))  ones(1,sum(m.Nj)-m.nu) 1];
%% Part 3.Run optimization with fmincon 
% Use this as initial guess 
dvar0 = [zeros(1,6) ... % u1 & u2 # 1-6
         0.2 0.8 ...     % t for u1 # 7-8 
         0.4 0.75 ...     % t for u2 # 9-10
         1.0];          % tF 
[Ctrl_matrix, pos_matrix, theta_vec, ...
          FHI_u1, FHI_u2] = Timehorizon_Sort(dvar0, m);
%-----------------------------------------------------------
load('functionHandle.mat')
[J0, grad0S] = costfun(dvar0,m,odefunvec,optODE,x0);
problem = createOptimProblem('fmincon','objective',...
    @(dvar)costfun(dvar,m,odefunvec,optODE,x0),...
    'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,...
    'nonlcon',  @(dvar)NLPcon(dvar,m));
tstart = tic ;
[dvarO,JO,exitflag,output,lambda,gradO]  = fmincon(problem);
dvarO
tend = toc(tstart);
ShowControl_VTNCVP(dvarO,m)
fprintf('Jmin = %9.6f\n', JO)
fprintf('CPU time = %4.3f s\n', tend)
ylim([-10.5 10.5])
print -dsvg Figure3.svg
return 
[J0, grad0S] = costfun(dvar0,m,odefunvec,optODE,x0);
[J1, grad1S] = costfun(dvar1,m,odefunvec,optODE,x0);
[J2, grad2S] = costfun(dvar2,m,odefunvec,optODE,x0);
% return 
% ms = MultiStart('StartPointsToRun','bounds');
% stpoints = RandomStartPointSet('NumStartPoints',10, ...
%     'ArtificialBound',10);
% [dvarO,JO,flag,outpt,allmins] = run(ms,problem,stpoints);

