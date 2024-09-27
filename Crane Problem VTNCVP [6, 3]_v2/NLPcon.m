function [c, ceq] = NLPcon(dvar, m, optODE)
% 
% dvar0 =[ zeros(1, m.Nj(1)) ...% u1 # 1-6 % Nj(1)
%          zeros(1, m.Nj(2)) ...% u2 # 7-9 % Nj(1)+1 : sum(Nj)
%             9/6:9/6:9-9/6 ... % t_u1 # 10-14 % sum(Nj)+1 : sum(Nj)+N1-1
%             9/3:9/3:9-9/3 ... % t_u2 # 15-16 % sum(Nj)+N1: 2*sum(Nj)-nu
%                       9.0];   % tF # 17    % 

%% 1 Inequality constraints 
for i = 1 : m.Nj(1) - 2 % 1 : 4
    c(i) = dvar(i+9) - dvar(i+10);
end 

for i = 1 : m.Nj(2) - 2 % 1 : 1 
    c(i+ m.Nj(1)-2) = dvar(i+14) - dvar(i+15);
end 
%% 2 Terminal Equality Constraints 

f = funint(dvar,m,optODE);
ceq(1) = f(1) - 10;
ceq(2) = f(2) - 14;
ceq(3) = f(3) - 0.0;
ceq(4) = f(4) - 2.5; 
ceq(5) = f(5);
ceq(6) = f(6);

end 