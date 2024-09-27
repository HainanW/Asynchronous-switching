function [c, ceq] = NLPcon(dvar, m, odefunvec, optODE)
% 
% u1   # 1-10 %          1 : Nj(1)
% u2   # 11-15 %    Nj(1)+1 : sum(Nj)
% t_u1 # 16-24 %    sum(Nj)+1 : sum(Nj)+1+Nj(1)-2=sum(Nj)+Nj(1)-1 
% t_u2 # 25-28 % sum(Nj)+Nj(1): 2*sum(Nj)-nu
% tF   # 29    % 2*sum(Nj)-nu+1


for i = 1 : m.Nj(1) - 2 % u1
    c(i) = dvar(i+sum(m.Nj)) - dvar(i+sum(m.Nj)+1);
end 

for i = 1 : m.Nj(2) - 2 % u2
    c(i+ m.Nj(1)-2) = dvar(i+sum(m.Nj)+m.Nj(1)-1) - dvar(i+sum(m.Nj)+m.Nj(1));
end 

%% 2 Terminal Equality Constraints 
[Ctrl_matrix, ~, theta_vec, ~, ~] = Timehorizon_Sort(dvar, m);
z0 = m.x0;
ts = linspace(0,m.p+1,m.p+2);
for ks = 1 : m.p + 1
    u1 = Ctrl_matrix(1,ks);
    u2 = Ctrl_matrix(2,ks);
    u3 = theta_vec(ks);
    odefun = odefunvec.state;
    odefun = @(t,Y) odefun(t,Y,u1,u2,u3);
    [~,z] = ode45(odefun,[ts(ks) ts(ks+1)],z0,optODE);
    z0 = z(end,:)';
end 
z_last = z0;

ceq(1) = z_last(1) - 10;
ceq(2) = z_last(2) - 14;
ceq(3) = z_last(3) - 0.0;
ceq(4) = z_last(4) - 2.5; 
ceq(5) = z_last(5);
ceq(6) = z_last(6);

end 