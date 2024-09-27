% What this file work: 
% 1. To check if the inequality constraints meet the requirements, 
% 2. To check if the termianl constraints meet the requirements, 
% 3. Be advised when running optimization in main_Case2_10_5.m,
% all warnings off 
% So re-run the optimal results again to see if any warnings or err msg. 

clear all 

m.nu = 2;   
m.Nj = [10,5]; % number of FEs for two controls 
m.nx = 7;  % number of state vars
m.p = sum(m.Nj) - m.nu; % p = 13: total number of time nodes 
m.m1 = 2; % m1: total number of the inequality constraints

m.x0 = [ 0 22 0 0.0 -1.0 0.0 0.0];
m.xF = [10 14 0 2.5  0.0 0.0];
x0.phi_u1 = [m.x0 0 0 zeros(1, m.nx+m.m1)]; 
x0.phi_u2 = [m.x0 0 0 zeros(1, m.nx+m.m1)];  
   x0.psi = [m.x0 0 0 zeros(1, m.nx+m.m1)];

sigma = 1e-6;
rho = 1;
epsilon = 2*sigma/((1+sqrt(5))*rho*m.m1*(m.p+1));
load("functionHandle.mat");

% % Round 4 
load('OptRes052923_P_10_5_Round5_MS_100.mat')
[Jmin_ascending, minrow_ascending] = sort(JOmin_eachrow_vec);
