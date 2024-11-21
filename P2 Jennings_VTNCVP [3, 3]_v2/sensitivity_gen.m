% function out = sentivity_gen(m)
%% 1 Equation & Parmeter Definition 
% 1.1 Parameter Definition 
clc 
Nj = [3,3]; % # of time stages 
nx = 3;
nu = 2; % nu = 3
P = sum(Nj(:)) - nu ; % P = 4
nx_t = nx + nx*( P + 1); % 3 + 3x5 = 18
xvarName = sprintfc('x%d(t)',1:nx_t); 
xvars = str2sym(xvarName);

syms(xvarName{:})
syms u1 u2 u3 

f_Orig(1,1) = u2;
f_Orig(2,1) = (-x1 + u1); 
f_Orig(3,1) = (-6*x1 - 12*x2 + 3*u1 + u2);
pfpx = Jacob_fun(f_Orig.*u3, nx); % pfpx = jacobian(f,sym('x',[1, nx]))
pfpu =  jacobian(f_Orig.*u3, [u1 u2 u3]);
pfpx_psi = pfpx; 
pfpu_psi = pfpu;

S_phi_u1 = sym(zeros(nx, Nj(1)));
S_phi_u2 = sym(zeros(nx, Nj(2)));
S_psi    = sym(zeros(nx, P + 1));
%% 2 Generate Sensitivity Matrix 
% Sensitivity Matrix of u1 
for i = 1 : nx*(Nj(1)) % 3x3 
    S_phi_u1(i) = str2sym(sprintfc('x%d(t)', i+nx));
end

% Sensitivity Matrix of u2
for i = 1 : nx*(Nj(2)) % 3x3 
    S_phi_u2(i) = str2sym(sprintfc('x%d(t)', i+nx));
end
 
% Sensitivity Matrix of theta 
for i = 1 : nx*(P + 1) % 3x5
    S_psi(i) = str2sym(sprintfc('x%d(t)', i+nx));
end 

pfpx_S_phi_u1 =  pfpx*S_phi_u1; % 3x3 
pfpx_S_phi_u2 =  pfpx*S_phi_u2; % 3x3 
pfpx_S_psi    =  pfpx_psi*S_psi; % 3x5 
 
odefunvec_phi_u1 = cell(1, Nj(1)); 
odefunvec_phi_u2 = cell(1, Nj(2));
odefunvec_psi    = cell(1, P+1);
%% 3 Generate function handle for phi_u1, phi_u2, and psi 

for i = 1 : Nj(1) % 1 : 3
    Fvec_phi_u1 = sym(zeros(nx*i+nx,1));
    pupw = zeros(nu+1, Nj(1)); % !3x3
    pupw(1, i) = 1;
    res_phi_u1 = pfpx_S_phi_u1 + pfpu*pupw; 
    Fvec_phi_u1(1:nx) = f_Orig.*u3;
    Fvec_phi_u1(nx*0*i+1+nx:nx*1*i+nx) = res_phi_u1(1:nx*i);
    odefunvec_phi_u1{i} = odeFunction(Fvec_phi_u1, ...
        xvars(1:nx + nx*Nj(1)),u1,u2,u3); % xvars = x1- x12 
end


for i = 1 : Nj(2) % 1 : 3
    Fvec_phi_u2 = sym(zeros(nx*i+nx,1));
    pupw = zeros(nu+1, Nj(2)); % !3x3
    pupw(2, i) = 1;
    res_phi_u2 = pfpx_S_phi_u2 + pfpu*pupw; 
    Fvec_phi_u2(1:nx) = f_Orig.*u3;
    Fvec_phi_u2(nx*0*i+1+nx:nx*1*i+nx) = res_phi_u2(1:nx*i);
    odefunvec_phi_u2{i} = odeFunction(Fvec_phi_u2, ...
        xvars(1:nx+nx*Nj(2)),u1,u2,u3);
end


for i = 1 : P + 1 % 1 : 5
    Fvec_psi = sym(zeros(nx*i+nx,1));
    pupw = zeros(nx, P+1);
    pupw(3, i) = 1;
    res_psi = pfpx_S_psi + pfpu_psi*pupw;
    Fvec_psi(1:nx) = f_Orig.*u3;
    Fvec_psi(nx*0*i+1+nx:nx*1*i+nx) = res_psi(1:nx*i);
    odefunvec_psi{i} = odeFunction(Fvec_psi, ...
        xvars(1:nx+nx*(P+1)), u1, u2, u3);
end

odefunvec.state = odeFunction(f_Orig.*u3, ...
        xvars(1:nx), u1, u2, u3); 
odefunvec.phi_u1 = odefunvec_phi_u1;
odefunvec.phi_u2 = odefunvec_phi_u2;
odefunvec.psi    = odefunvec_psi; 
%% 4 Save the function handle 
save('functionHandle.mat', "odefunvec")