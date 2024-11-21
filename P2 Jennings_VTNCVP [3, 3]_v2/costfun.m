function [J, grad] = costfun(dvar, m, odefunvec, optODE, x0)


%% 1. Sort the time sequence by input

[Ctrl_matrix,pos_matrix,theta_vec,...
    FHI_u1, FHI_u2] = Timehorizon_Sort(dvar, m);
tspan  = linspace(0, sum(m.Nj) - m.nu + 1, sum(m.Nj) - m.nu + 2);
%% 2. Calculate the objective and Gradient

z0_phi_u1 = x0.phi_u1;
odefunvec_phi_u1 = odefunvec.phi_u1;
z0_phi_u2 = x0.phi_u2;
odefunvec_phi_u2 = odefunvec.phi_u2;
z0_psi = x0.psi;
odefunvec_psi = odefunvec.psi;

for i = 1 : m.p + 1

    u1 = Ctrl_matrix(1,i);
    u2 = Ctrl_matrix(2,i);
    u3 = theta_vec(i);
    odefun_u1 = odefunvec_phi_u1{FHI_u1(i)};
    odefun_u2 = odefunvec_phi_u2{FHI_u2(i)};
    odefun_u3 = odefunvec_psi{i};
    odefun_u1 = @(t,Y)odefun_u1(t, Y, u1, u2, u3);
    odefun_u2 = @(t,Y)odefun_u2(t, Y, u1, u2, u3);
    odefun_u3 = @(t,Y)odefun_u3(t, Y, u1, u2, u3);
    [~,z_u1] = ode45(odefun_u1,[tspan(i) tspan(i+1)],z0_phi_u1,optODE);
    z0_phi_u1 = z_u1(end,:)';
    [~,z_u2] = ode45(odefun_u2,[tspan(i) tspan(i+1)],z0_phi_u2,optODE);
    z0_phi_u2 = z_u2(end,:)';
    [~,z_u3] = ode45(odefun_u3,[tspan(i) tspan(i+1)],z0_psi,optODE);
    z0_psi = z_u3(end,:)';
    if i ~= m.p + 1
        % u1 
        if FHI_u1(i+1) - FHI_u1(i) == 1 
            z0_phi_u1 = [z0_phi_u1; zeros(3,1)];
        else % FHI_u1(i+1) - FHI_u1(i) == 0
            z0_phi_u1 = z0_phi_u1; 
        end 
        % u2 
        if FHI_u2(i+1) - FHI_u2(i) == 1
            z0_phi_u2 = [z0_phi_u2; zeros(3,1)];
        else % FHI_u2(i+1) - FHI_u2(i) == 0 
            z0_phi_u2 = z0_phi_u2; 
        end
        % u3 
        z0_psi = [z0_psi; zeros(3,1)];
    else % i == m.p + 1
        break 
    end

end
%% 3. Gradient Conversion
% Be advised:
% pc_pt does not include partial c/partial t_F (pc_ptF)
pg_pt = zeros(1,m.p); % 1x4
pg_ptheta = z0_psi(6:3:18); % 1x5
pos_sorted = sortrows(pos_matrix);
pos_vec = pos_sorted(:,3)';
for i = 1 : length(pg_pt) 
    pg_pt(i)= pg_ptheta(pos_vec(i)) - pg_ptheta(pos_vec(i)+1);
end
%% 4. Output Result
J = z0_psi(3);
if nargout > 1
    grad(1:3) = z0_phi_u1(6:3:end);
    grad(4:6) = z0_phi_u2(6:3:end);
    grad(7:10)= pg_pt';
    grad(11)  = pg_ptheta(end);
else
    grad = [];
end