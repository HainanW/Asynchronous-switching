clear all
%% Crane Problem in VTNCVP
%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-10;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on','Display', 'iter', ...
    'TolX', NLPtol, 'TolFun', NLPtol, 'TolCon', NLPtol,...
    'MaxFunEvals',20000, 'MaxIter',500, ...
    'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars
% number of control vars 
m.nu = 2;   
% number of FDs for two controls 
m.Nj = [6,3];
% number of state vars
m.nx = 7;    
% p: total number of time nodes 
m.p = sum(m.Nj) - m.nu; % p = 7 
% m1: total number of the inequality constraints
m.m1 = 2; 
m.x0 = [ 0 22 0 0.0 -1.0 0.0 0.0];
m.xF = [10 14 0 2.5 0.0 0.0];
x0.phi_u1 = [m.x0 0 0 zeros(1, m.nx+m.m1)]; 
x0.phi_u2 = [m.x0 0 0 zeros(1, m.nx+m.m1)];  
   x0.psi = [m.x0 0 0 zeros(1, m.nx+m.m1)];  
% Decision variable input sequence 
% u1 u2 t1vec t2vec tF 
lb = [-2.83374*ones(1,m.Nj(1)) ... u1 lb
      -0.80865*ones(1,m.Nj(2)) ... u2 lb
       zeros(1,m.p)        ... t^i_k lb
       9.0];                       %  tF lb
ub = [ 2.83374*ones(1,m.Nj(1)) ...
       0.71265*ones(1,m.Nj(2)) ...
       9*ones(1,m.p)  ...
       9.0];
sigma = 1e-6;
rho = 1;
epsilon = 2*sigma/((1+sqrt(5))*rho*m.m1*(m.p+1));

%% Part 3.Run optimization with fmincon 
% Use this as initial guess 
% dvar0 can not be a structure array 
dvar0 = [ 0.5*ones(1, m.Nj(1)) ...% u1 # 1-6 % Nj(1)
          0.5*ones(1, m.Nj(2)) ...% u2 # 7-9 % Nj(1) +1 : sum(Nj)
             9/6:9/6:9-9/6 ... % t_u1 # 10-14 % sum(Nj)+1 : sum(Nj)+N1-1
             9/3:9/3:9-9/3 ... % t_u2 # 15-16 % sum(Nj)+N1: 2*sum(Nj)-nu
                       9.0];   % tF # 17    % 
dvar0 = [ 1.0*ones(1, m.Nj(1)) ...
          0.5*ones(1, m.Nj(2)) ... 
          1.0 3.0 4.5 6.0 7.5  ... 
          3.0 6.0 ...
          9.0];

[Ctrl_matrix, pos_matrix, ... 
   theta_vec, FHI_u1, FHI_u2] = Timehorizon_Sort(dvar0, m);

load('functionHandle.mat')
% load('OptRes032723_Round1.mat', 'dvarOMat','minrow_idx','points');
% dvar0 = dvarOMat(minrow_idx,:);
[J0, grad0] = costfun(dvar0,m,odefunvec,optODE,x0,epsilon,rho);
% ShowControl_VTNCVP(points(minrow_idx,:), m)
% ShowControl_VTNCVP(dvarOMat(minrow_idx,:), m)
% 
% dvar0 = points(minrow_idx,:);
% dvar0 = dvar0;
dvarOvec = []; 
[JO,gradO]= costfun(dvar0,m,odefunvec,optODE, x0, epsilon, rho);
% consfunTS(dvar0, m, optODE, x0_orig)
JOvec = zeros(1,5);
counter = 1;
JOvec(counter) = JO; 
tendvec = [];
while counter <=10 
    %%
    tstart = tic;
    problem = createOptimProblem('fmincon','objective',...
        @(dvar)costfun(dvar,m,odefunvec,optODE, x0, epsilon, rho),...
        'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,'nonlcon',...
        @(dvar)NLPcon(dvar,m,optODE));
    [dvarO,JO,exitflag]= fmincon(problem);
    tend = toc(tstart);
    tendvec = [tendvec tend];
    fprintf('tf= %9.7f\n',JO)
    if abs(JO - JOvec(counter)) > sigma
        epsilon = epsilon*0.1; 
        rho = rho*10;
        JOvec(counter+1) = JO;
        dvar0 = dvarO;
        dvarOvec = [dvarOvec; dvarO];
    else 
        JOvec(counter+1) = JO;
        dvarOvec = [dvarOvec; dvarO];
        break 
    end 
    counter = counter +1; 
end
save('OptRes060723_Round3_singlestart.mat')

ShowControl_VTNCVP(dvarO, m)
return 
%% Part 4. Pre-preparation for Multi-Start 
ms_number = 50; % How many multi-start points 
% points: matrix store the ms intial points 
% Set up the first 10 columns
% For now, the points are from Round 6 

points = rand(ms_number, m.Nj(1)) * 5.66748 - 2.83374;
points(:, m.Nj(1) + 1 : sum(m.Nj)) = rand(ms_number, m.Nj(2)) * 1.52130 - 0.80865;
random_u1 = 0 + 9*rand(ms_number, m.Nj(1)-1);
random_u1 = sort(random_u1, 2);
random_u2 = 0 + 9*rand(ms_number, m.Nj(2)-1);
random_u2 = sort(random_u2, 2);

points(:, sum(m.Nj) + 1 : sum(m.Nj) + m.Nj(1) - 1) = random_u1; 
points(:, sum(m.Nj) + m.Nj(1) :  2*sum(m.Nj)-m.nu) = random_u2;
points(:, 2*sum(m.Nj)-m.nu + 1) = 9*ones(1, ms_number);

JOMat = zeros(ms_number,6); % 6 columns = 1(JOinitial)+while counter<5
dvarOMat = zeros(ms_number, 2*sum(m.Nj)-m.nu + 1);
ExitFlagMat = cell(ms_number,5);
tendMat = zeros(ms_number,1);% total time for each ms including all iterations
counterIdxVec = zeros(ms_number,1);


%% Part 5. 
for i = 1 :  ms_number
    dvar0 = points(i,:); 
    % Initialize the parameters 
    rho = 1;
    epsilon = 2*sigma/((1+sqrt(5))*rho*m.m1*(m.p+1));
    JO= costfun(dvar0,m,odefunvec,optODE,x0, epsilon, rho);
    counter = 1;
    JOMat(i, 1) = JO;
    tstart = tic;
    while counter <= 5
        problem = createOptimProblem('fmincon','objective',...
            @(dvar)costfun(dvar,m,odefunvec,optODE, x0, epsilon, rho),...
            'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,'nonlcon',...
            @(dvar)NLPcon(dvar, m, optODE));
        [dvarO,JO,ExitFlag]= fmincon(problem);
        % tend = toc(tstart);
        % tend = [tendvec tend];
        % fprintf('tf= %9.8f\n',JO)
        if abs(JO - JOMat(i, counter)) > sigma
            epsilon = epsilon*0.1;
            rho = rho*10;
            dvar0 = dvarO;
            JOMat(i, counter+1) = JO;
            ExitFlagMat(i,counter) = {ExitFlag};
            counter = counter+1;
        else
            JOMat(i, counter+1) = JO;
            ExitFlagMat(i,counter) = {ExitFlag}; 
            counter = counter+1;
            break
        end
    end
    tend = toc(tstart);
    tendMat(i) = tend;
    dvarOMat(i,:) = dvarO;
    counterIdxVec(i,:) = counter; 
    fprintf('Iteration %2.0d is done!\n', i) 
    fprintf('Jmin = %9.6f\n', JO)
    fprintf('CPU Time is %8.2f s\n', tend)
end

% Find the minimum of JOMat
JOmin_eachrow_vec = zeros(1,ms_number);
for i = 1 : ms_number 
    JOmin_eachrow_vec(i) = JOMat(i, counterIdxVec(i));
end 
[Jmin, minrow_idx] = min(JOmin_eachrow_vec);
% save('OptRes032723_Round1.mat')
ShowControl_VTNCVP(dvarOMat(minrow_idx), m)
% ms = MultiStart('StartPointsToRun','bounds');
% stpoints = RandomStartPointSet('NumStartPoints',10, ...
%     'ArtificialBound',10);
% [dvarO,JO,flag,outpt,allmins] = run(ms,problem,stpoints);

