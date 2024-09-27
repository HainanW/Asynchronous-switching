clear all 
%% Part 0.Crane Problem 
% [14, 7] Scenario 

%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-10;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on','Display', 'iter', ...
    'TolX', NLPtol, 'TolFun', NLPtol, 'TolCon', NLPtol,...
    'MaxFunEvals',20000, 'MaxIter',500,'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars

% number of control vars 
m.nu = 2;   
% number of FDs for two controls 
m.Nj = [14,7];
% number of state vars
m.nx = 7;    
% p: total number of time nodes 
m.p = sum(m.Nj) - m.nu; % p = 19 
% m1: total number of the inequality constraints
m.m1 = 2; 
m.x0 = [ 0 22 0 0.0 -1.0 0.0 0.0];
m.xF = [10 14 0 2.5  0.0 0.0];
x0.phi_u1 = [m.x0 0 0 zeros(1, m.nx+m.m1)]; 
x0.phi_u2 = [m.x0 0 0 zeros(1, m.nx+m.m1)];  
   x0.psi = [m.x0 0 0 zeros(1, m.nx+m.m1)];
% Decision variable input sequence 
% u1 u2 t1vec t2vec tF 
lb = [-2.83374*ones(1,sum(m.Nj(1))) ... u1 lb
      -0.80865*ones(1,sum(m.Nj(2))) ... u2 lb
      zeros(1,m.p)        ... delt lb
      9.0];                       %  tF lb
ub = [ 2.83374*ones(1,sum(m.Nj(1))) ...
       0.71265*ones(1,sum(m.Nj(2))) ...
       9*ones(1,m.p)  ...
       9.0];
sigma = 1e-6;
rho = 1;
epsilon = 2*sigma/((1+sqrt(5))*rho*m.m1*(m.p+1));

load('functionHandle.mat','odefunvec')
%% Part 3.Run Single optimization  
% Prepare for ms-start 
% Use this as initial guess 
% tart_time = 0;
% end_time = 9;
% 
% % 将时间段划分为14段
% num_segments = 14;
% segment_length = (end_time - start_time) / num_segments;
% 
% % 计算时间点
% time_points = start_time + segment_length : segment_length : end_time - segment_length;

% 0.010306
dvarO_opt = [2.8337 -0.0812 0.2813 0.5156 0.6707 0.8092 0.9051 1.0126 ...
   1.0640  1.1071   0.9635   0.5581  -0.1265  2.8337   ...
   0.0253 -0.0051   0.1421  -0.1369  -0.0179 -0.0001  0.7127 ...
   0.1368  0.7472   1.4971   2.2008   2.8930  3.6018  4.3498 ...
   5.0743  5.7594   6.5839   7.5830   8.2710  8.7714  ...
   1.0036  2.4421   3.6072   4.7734   6.0943  7.5968  9.0000];

dvar0 = [2.0*ones(1, m.Nj(1)) ... % u1 # 1-14 -14 pts 
         0.5*ones(1, m.Nj(2)) ... % u2 # 15-21 -7 pts 
         9/14 : 9/14 : 9-9/14 ...   % t for u1 # 22-34 -13 pts 
         9/7 : 9/7 : 9-9/7 ...   % t for u2 # 35-40 -6 pts 
         9.0];                 % tF   # 41
dvar0 = [2.8337 -0.0616:(1.0674--0.0616)/9:1.0674 0.7164 -0.0399 2.8337 ... % u1 # 1-14
         0.0006 0.0006 0.0697 0.0697 -0.1098  -0.0000 0.7127  ... % u2 # 15-21
         0.1357:(7.2281-.1357)/10:7.2281 8.1496  8.7770...   % t for u1 # 22-34 -13 pts 
         1.1971 2.3941 3.5245 4.6548  6.1026  7.5968 ...   % t for u2 # 35-40 -6 pts 
         9.0];                 % tF   # 41
ShowControl_VTNCVP(dvar0,m)
[Ctrl_matrix, pos_matrix, theta_vec, ...
    theta_u1, theta_u2] = Timehorizon_Sort(dvar0, m);
% Comment following line once the sensitivity function generated 
% sensitivity_gen(m)

load('functionHandle.mat')
load('OptRes053023_P_14_7_SS_Round1.mat')
ShowControl_VTNCVP(dvarO,m)
dvarOvec = []; 
[JO,gradO]= costfun(dvar0,m,odefunvec,optODE, x0, epsilon, rho);
% consfunTS(dvar0, m, optODE, x0_orig)
JOvec = zeros(1,5);
counter = 1;
JOvec(counter) = JO; 
tendvec = [];
while counter <= 5
    %%
    tstart = tic;
    problem = createOptimProblem('fmincon','objective',...
        @(dvar)costfun(dvar,m,odefunvec,optODE, x0, epsilon, rho),...
        'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,'nonlcon',...
        @(dvar)NLPcon(dvar,m,odefunvec, optODE));
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
% end % Comment on this.  
% save('OptRes053023_P_14_7_SS_Round1.mat')
ShowControl_VTNCVP(dvarO,m)
return 
%% Part 4. Pre-preparation for Multi-Start 
ms_number = 100; % How many multi-start points 
% points: matrix store the ms intial points 
% Set up the first 10 columns
% For now, the points are from Round 6 

points = rand(ms_number, m.Nj(1)) * 5.66748 - 2.83374;
points(:, m.Nj(1) + 1 : sum(m.Nj)) = rand(ms_number, m.Nj(2)) * 1.52130 - 0.80865;
random_u1 = 0 + 9*rand(ms_number, m.Nj(1)-1);
random_u1 = sort(random_u1, 2); % initial t^i_k needs to be in ascending order
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

%% Part 5. MS-start 
for i = 1 : ms_number 
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
            @(dvar)NLPcon(dvar,m,odefunvec, optODE));
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
save('OptRes052623_P_14_7_Round2.mat')
ShowControl_VTNCVP(dvarOMat(minrow_idx), m)