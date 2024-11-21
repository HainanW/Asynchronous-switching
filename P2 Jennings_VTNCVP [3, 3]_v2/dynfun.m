function dydt = dynfun(t,x,u1,u2)

dydt(1) = u2;
dydt(2) = -x(1) + u1; 
dydt(3) = -6*x(1) - 12*x(2) + 3*u1 + u2;

end 
