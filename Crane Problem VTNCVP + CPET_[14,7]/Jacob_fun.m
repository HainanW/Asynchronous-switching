function out = Jacob_fun(f,nx)
vars = str2sym(sprintfc('x%d(t)', 1:nx));
Jacob = sym(zeros(nx,nx));
for i = 1 : nx
    for j = 1 : nx
        Jacob(i,j) = diff(f(i),vars(j));
    end
end
out = Jacob;
end