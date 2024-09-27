function out = ShowControl_VTNCVP(dvarO, m)

% dvar0 = [ 0.5*ones(1, m.Nj(1)) ...% u1 # 1-6 % Nj(1)
%           0.5*ones(1, m.Nj(2)) ...% u2 # 7-9 % Nj(1) +1 : sum(Nj)
%              9/6:9/6:9-9/6 ... % t_u1 # 10-14 % sum(Nj)+1 : sum(Nj)+N1-1
%              9/3:9/3:9-9/3 ... % t_u2 # 15-16 % sum(Nj)+N1: 2*sum(Nj)-nu
%                        9.0];   % tF # 17    % 

% m: m.nu, m.Nj, m.nx, m.p, m.m1, m.x0, m.xF
tF = dvarO(end);
str = cell(1,2);
str{1} = 'k-';
str{2} = 'k--';

for i = 1 : m.nu %
    uvec = dvarO(sum(m.Nj(1:i-1))+1:sum(m.Nj(1:i)));
    tvec = dvarO(sum(m.Nj)+ (sum(m.Nj(1:i-1))+ 2-i : sum(m.Nj(1:i))-i ));
    ts_plot   = zeros(1, 2+2*(m.Nj(i)-1));
    uvec_plot = zeros(1, 2+2*(m.Nj(i)-1));
    ts_plot(1,1) = 0;
    for counter = 1 : m.Nj(i)
        if counter == m.Nj(i)
              ts_plot(2*counter) = tF;
            uvec_plot(2*counter-1:2*counter) = uvec(counter);
        else  
              ts_plot(2*counter  : 2*counter+1) = tvec(counter);
            uvec_plot(2*counter-1 : 2*counter)  = uvec(counter);
        end     
    end
    plot(ts_plot,uvec_plot,str{i},'LineWidth',1.5);
    hold on 
end

hold off 
xlabel('Time')
ylabel('Controls')
xlim([0 tF])
grid on 
legend({'\it{u}_{1}','\it{u}_{2}'},'Location','north')
end 