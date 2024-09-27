function out = funint(dvar,m,optODE)

% integration function for piecewise constant control
[Ctrl_matrix, ~, theta_vec, ~, ~] = Timehorizon_Sort(dvar, m);
z0 = m.x0;
ts = linspace(0,m.p+1,m.p+2);
for ks = 1 : m.p + 1
    [~,z] = ode45(@(t,x)dyneqn(t,x,Ctrl_matrix(1,ks),...
        Ctrl_matrix(2,ks),theta_vec(ks)),...
        [ts(ks) ts(ks+1)],z0,optODE);
    z0 = z(end,:)';
end
out = z0; 
end