function [J, grad, pathviolation_1, pathviolation_2] = ...
costfun_special(dvar,m,odefunvec,optODE,x0,epsilon, rho)
%% 1. Sort the time sequence by input

[Ctrl_matrix,pos_matrix,theta_vec,...
    FHI_u1, FHI_u2] = Timehorizon_Sort(dvar, m);
tspan  = linspace(0, m.p+1, m.p+2);

%% 2. Calculate the objective and Gradient

z0_phi_u1 = x0.phi_u1;
odefunvec_phi_u1 = odefunvec.phi_u1;
z0_phi_u2 = x0.phi_u2;
odefunvec_phi_u2 = odefunvec.phi_u2;
z0_psi = x0.psi;
odefunvec_psi = odefunvec.psi;

for i = 1 : m.p + 1 % 14
    u1 = Ctrl_matrix(1,i);
    u2 = Ctrl_matrix(2,i);
    u3 = theta_vec(i);
    odefun_u1 = odefunvec_phi_u1{FHI_u1(i)};
    odefun_u2 = odefunvec_phi_u2{FHI_u2(i)};
    odefun_u3 = odefunvec_psi{i};
    odefun_u1 = @(t,Y)odefun_u1(t, Y, u1, u2, u3, epsilon);
    odefun_u2 = @(t,Y)odefun_u2(t, Y, u1, u2, u3, epsilon);
    odefun_u3 = @(t,Y)odefun_u3(t, Y, u1, u2, u3, epsilon);
    [~,z_u1] = ode45(odefun_u1,[tspan(i) tspan(i+1)],z0_phi_u1,optODE);
    z0_phi_u1 = z_u1(end,:)';
    [~,z_u2] = ode45(odefun_u2,[tspan(i) tspan(i+1)],z0_phi_u2,optODE);
    z0_phi_u2 = z_u2(end,:)';
    [~,z_u3] = ode45(odefun_u3,[tspan(i) tspan(i+1)],z0_psi,optODE);
       z0_psi = z_u3(end,:)';
    if i ~= m.p + 1
        % u1
        if FHI_u1(i+1) - FHI_u1(i) == 1
            z0_phi_u1 = [z0_phi_u1; zeros((m.nx+m.m1),1)];
        else % FHI_u1(i+1) - FHI_u1(i) == 0
            z0_phi_u1 = z0_phi_u1;
        end
        % u2
        if FHI_u2(i+1) - FHI_u2(i) == 1
            z0_phi_u2 = [z0_phi_u2; zeros((m.nx+m.m1),1)];
        else % FHI_u2(i+1) - FHI_u2(i) == 0
            z0_phi_u2 = z0_phi_u2;
        end
        % u3
        z0_psi = [z0_psi; zeros(m.nx+m.m1,1)];
    else % i == m.p + 1
        break
    end

end
%% 3. Gradient Conversion
% Be advised:
% pc_pt does not include partial c/partial t_F (pc_ptF)
pg1_pt = zeros(1,m.p); % x7_t
pg1_ptheta = z0_psi(16:9:end);
pg2_pt = zeros(1,m.p); % x8_t
pg2_ptheta = z0_psi(17:9:end);
pg3_pt = zeros(1,m.p); % x9_t
pg3_ptheta = z0_psi(18:9:end);

pos_sorted = sortrows(pos_matrix);
pos_vec = pos_sorted(:,3)';

for i = 1 : m.p
    pg1_pt(i) = pg1_ptheta(pos_vec(i)) - pg1_ptheta(pos_vec(i)+1);
    pg2_pt(i) = pg2_ptheta(pos_vec(i)) - pg2_ptheta(pos_vec(i)+1);
    pg3_pt(i) = pg3_ptheta(pos_vec(i)) - pg3_ptheta(pos_vec(i)+1);
end
%% 4. Output Obj and Gradient
J = z0_psi(7) + rho*sum(z0_psi(8:9));
pathviolation_1 = z0_psi(8);
pathviolation_2 = z0_psi(9);

if nargout > 1
    grad(           1 :   m.Nj(1)) = z0_phi_u1(16:9:end) + ...
        rho* (z0_phi_u1(17:9:end) + z0_phi_u1(18:9:end)) ; % u1 # 1-6 -6 values

    grad( m.Nj(1) + 1 : sum(m.Nj)) = z0_phi_u2(16:9:end) + ...
        rho* (z0_phi_u2(17:9:end) + z0_phi_u2(18:9:end))   ; % u2 # 7-9 -3 values

    grad(sum(m.Nj)+ 1 : sum(m.Nj) + m.p) = pg1_pt + ...
        rho*(pg2_pt + pg3_pt)  ;         % t for u1, u2 -7 pts

    grad(sum(m.Nj) + m.p + 1) = pg1_ptheta(end) + ...
        rho*(pg2_ptheta(end) + pg3_ptheta(end))    ; % tF #
else
    grad = [];
end
end