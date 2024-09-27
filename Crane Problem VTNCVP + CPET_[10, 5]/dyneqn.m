function dx = dyneqn(t, x, u1, u2, u3)

dx = zeros(7,1);
dx(1) = u3*(x(4));
dx(2) = u3*(x(5));
dx(3) = u3*(x(6));
dx(4) = u3*(u1+17.2656*x(3));
dx(5) = u3*(u2);
dx(6) = u3*(-(u1+27.0756*x(3)+2*x(5)*x(6))/x(2));
dx(7) = u3*(x(3)^2+x(6)^2);

end