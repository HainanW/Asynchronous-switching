function [c, ceq] = NLPcon(dvar, m, optODE)
% 
% u1   # 1-3 % Nj(1)
% u2   # 4-6 % Nj(1)+1 : sum(Nj)
% t_u1 # 7-8 % sum(Nj)+1 : sum(Nj)+N1-1
% t_u2 # 9-10 % sum(Nj)+N1: 2*sum(Nj)-nu
% tF   # 11   % 


for i = 1 : m.Nj(1) - 2 % u1
    c(i) = dvar(i+sum(m.Nj)) - dvar(i+sum(m.Nj)+1);
end 

for i = 1 : m.Nj(2) - 2 % u2
    c(i+ m.Nj(1)-2) = dvar(i+sum(m.Nj)+m.Nj(1)-1) - dvar(i+sum(m.Nj)+m.Nj(1));
end 

ceq = []; 

end 