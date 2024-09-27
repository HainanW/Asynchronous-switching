function sensitivity_gen(m)
%% Case 2 Crane Problem 
%% 1 Equation & Parmeter Definition 
% Now solving Crane Problem 
% 1.1 Parameter Definition 
 
Nj = m.Nj; % [14,7] % # of time stages 
nx = m.nx; % 7
nu = m.nu; % 2; % nu = 2
m1 = m.m1; % 2; 
P = m.p ; % P = 19
nx_t = nx + m1 + (nx+m1)*(P + 1);% 9+ 9*20=189
xvarName = sprintfc('x%d(t)',1:nx_t); 
xvars = str2sym(xvarName);

syms(xvarName{:})
% syms(thetavarName{:})
syms u1 u2 u3 epsilon 

%% 2 Constraint Related 
G1 = formula(x4^2-2.5^2);
G2 = formula(x5^2-1.0^2);
p1G1 = 0.5*(G1+(G1^2+4*epsilon^2)^(1/2));
p1G2 = 0.5*(G2+(G2^2+4*epsilon^2)^(1/2));

%% 3 Calculate Sensitivity eqn
f_Orig(1,1) = x4;
f_Orig(2,1) = x5; 
f_Orig(3,1) = x6;
f_Orig(4,1) = u1+17.2656*x3;
f_Orig(5,1) = u2;
f_Orig(6,1) = -(u1+27.0756*x3+2*x5*x6)/x2;
f_Orig(7,1) = x3^2 + x6^2; 
f_Orig(8,1) = p1G1 ; % x4 smooth penalty 
f_Orig(9,1) = p1G2 ; % x5 smooth penalty 

pfpx = Jacob_fun(f_Orig.*u3, nx+m1); % pfpx = jacobian(f,sym('x',[1, nx]))
pfpu =  jacobian(f_Orig.*u3, [u1 u2 u3]);

% Fvec_phi_u1 = sym(zeros(nx+m1+(nx+m1)*Nj(1), 1)); % 9+9*10 = 99
% Fvec_phi_u2 = sym(zeros(nx+m1+(nx+m1)*Nj(2), 1)); % 9+9*5 = 54
% Fvec_psi = sym(zeros(nx+m1+(nx+m1)*(P+1), 1)); % 9+9*14 = 135


S_phi_u1 = sym(zeros(nx+m1, Nj(1)));
S_phi_u2 = sym(zeros(nx+m1, Nj(2)));
S_psi    = sym(zeros(nx+m1, P+1));

%% 4 Generate Sensitivity Matrix 
% Sensitivity Matrix of u1 
for i = 1 : (nx+m1)*Nj(1) % 9x10 = 90
    S_phi_u1(i) = str2sym(sprintfc('x%d(t)', i+nx+m1));
end

% Sensitivity Matrix of u2
for i = 1 : (nx+m1)*Nj(2) % 9x5 = 45
    S_phi_u2(i) = str2sym(sprintfc('x%d(t)', i+nx+m1));
end
 
% Sensitivity Matrix of theta 
for i = 1 : (nx+m1)* (P+1) % 9x14 = 126
    S_psi(i) = str2sym(sprintfc('x%d(t)', i+nx+m1));
end 

pfpx_S_phi_u1 =  pfpx*S_phi_u1; % 9x10
pfpx_S_phi_u2 =  pfpx*S_phi_u2; % 9x5 
pfpx_S_psi    =  pfpx*S_psi;    % 9x14

odefunvec_phi_u1 = cell(1, Nj(1)); 
odefunvec_phi_u2 = cell(1, Nj(2));
odefunvec_psi    = cell(1, P+1);

%% 5 Generate function handle for phi_u1, phi_u2, and psi 

for i = 1 : Nj(1) % 1 : 14
    Fvec_phi_u1 = sym(zeros((nx+m1)*i + (nx+m1),1));
    pupw = zeros(nu+1, Nj(1)); % !3x10
    pupw(1, i) = 1;
    res_phi_u1 = pfpx_S_phi_u1 + pfpu*pupw; % 9x10
    
    Fvec_phi_u1(1:nx+m1) = f_Orig.*u3;
    Fvec_phi_u1(nx+m1+1:end) = res_phi_u1(1:(nx+m1)*i); 
    odefunvec_phi_u1{i} = odeFunction(Fvec_phi_u1, ...
        xvars(1:(nx+m1)+(nx+m1)*Nj(1)), u1, u2, u3, epsilon); 
end

for i = 1 : Nj(2) % 1 : 7
    Fvec_phi_u2 = sym(zeros((nx+m1)*i + (nx+m1),1));
    pupw = zeros(nu+1, Nj(2)); % !3x7
    pupw(2, i) = 1;
    res_phi_u2 = pfpx_S_phi_u2 + pfpu*pupw; % 9x7
    
    Fvec_phi_u2(1:nx+m1) = f_Orig.*u3;
    Fvec_phi_u2(nx+m1+1:end) = res_phi_u2(1:(nx+m1)*i); 
    odefunvec_phi_u2{i} = odeFunction(Fvec_phi_u2, ...
        xvars(1:(nx+m1)+(nx+m1)*Nj(2)), u1, u2, u3, epsilon);
end

for i = 1 : P+1 % 1 : 20
    Fvec_psi = sym(zeros((nx+m1)*i+(nx+m1),1));
    pupw = zeros(nu+1, P+1);
    pupw(3, i) = 1;
    res_psi = pfpx_S_psi + pfpu*pupw; % 9x20

    Fvec_psi(1:nx+m1) = f_Orig.*u3; 
    Fvec_psi(nx+m1+1:end) = res_psi(1:(nx+m1)*i);
    odefunvec_psi{i} = odeFunction(Fvec_psi, ...
        xvars(1: nx+m1+(nx+m1)*(P+1)), u1, u2, u3, epsilon);
end

odefun_state = odeFunction(f_Orig(1:nx).*u3,xvars(1:(nx+m1)),u1,u2,u3);

odefunvec.phi_u1 = odefunvec_phi_u1;
odefunvec.phi_u2 = odefunvec_phi_u2;
odefunvec.psi    = odefunvec_psi; 
odefunvec.state  = odefun_state; 

save('functionHandle.mat', "odefunvec")

end 