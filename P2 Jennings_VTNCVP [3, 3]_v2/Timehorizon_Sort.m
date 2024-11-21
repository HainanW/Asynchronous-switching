function [Ctrl_matrix, pos_matrix, theta_vec, ...
          FHI_u1, FHI_u2] = Timehorizon_Sort(dvar, m)

% u1   # 1-3 % Nj(1)
% u2   # 4-6 % Nj(1)+1 : sum(Nj)
% t_u1 # 7-8 % sum(Nj)+1 : sum(Nj)+N1-1
% t_u2 # 9-10 % sum(Nj)+N1: 2*sum(Nj)-nu
% tF   # 11   % 
% m.Nj = [3,3];
% m.nu = 2; % number of control vars
% m.nx = 2; % number of original state vars
% m.p = 4;
uvec_u1 = dvar(1:m.Nj(1));
uvec_u2 = dvar(m.Nj(1)+1 : sum(m.Nj));
tvec_u1 = dvar(sum(m.Nj)+1 : sum(m.Nj)+m.Nj(1)-1);
tvec_u2 = dvar(sum(m.Nj)+m.Nj(1) : 2*sum(m.Nj)-m.nu);
tF = dvar(end);

tvec_all = NaN(m.nu,max(m.Nj)-1);
tvec_all(1,1:m.Nj(1)-1) = tvec_u1;
tvec_all(2,1:m.Nj(2)-1) = tvec_u2; 

% 给定nu个行向量，每个行向量中的元素都是从小到大排列的
% 所有行向量的长度各不相同，每个行向量的长度加1后记录在Nj向量中
% 所有nu个数列元素数目为p=sum(Nj)-nu

% 假设nu个行向量存储在一个nu行，最多有max(Nj)列的矩阵tvec_all中
% 其中空缺的位置用NaN填充，例如，第i个行向量的元素存储在tvec_all(i, 1:Nj(i)-1)中，
% 其余位置为NaN
% 先将tvec_all转换为一个列向量vec_all，将NaN替换为Inf以便排序
% 转置成列向量,[] let reshape automatically calculate appropriate number of rows
colvec_all = reshape(tvec_all', [], 1); 
colvec_all(isnan(colvec_all)) = Inf; % 将NaN替换为Inf
[sorted_vec, idx] = sort(colvec_all); % 排序
% 因为ind2sub的排列顺序,调换了col_idx和row_idx的顺序
[col_idx, row_idx] = ind2sub(size(tvec_all'), idx); % 将线性索引转换为行列索引

% 创建pos_matrix矩阵
pos_matrix = [row_idx(1:m.p) col_idx(1:m.p) (1:m.p)'];
u1_ind = pos_matrix(pos_matrix(:, 1) == 1, 3);
u2_ind = pos_matrix(pos_matrix(:, 1) == 2, 3);

% Def: 
% FH_u1 At each stage which sensitivity _u1 to use 
% FH_u2 At each stage which sensitivity _u2 to use 
FHI_u1 = zeros(1, m.p+1); % FHI Function Handle Index _u1
FHI_u2 = zeros(1, m.p+1); 

Ctrl_matrix = zeros(m.nu,m.p + 1);

counter = 1; 
for i = 1 : m.Nj(1)
    if i == 1 % 1st Stage 
        for j = 1 : u1_ind(1)
            FHI_u1(counter) = i;
            counter = counter+1; 
            Ctrl_matrix(1,j) = uvec_u1(i);
        end 
    elseif i == m.Nj(1) % 6
        for j = u1_ind(end) + 1 : m.p + 1
            FHI_u1(counter) = i; 
            counter = counter+1;
            Ctrl_matrix(1,j) = uvec_u1(i);
        end 
    else % 2 3 4 5 
        for j = u1_ind(i-1) + 1 : u1_ind(i)
            FHI_u1(counter) = i;
            counter = counter+1;
            Ctrl_matrix(1,j) = uvec_u1(i);
        end 
    end 
end 

counter = 1; 
for i = 1 : m.Nj(2)
    if i == 1 % 1st Stage 
        for j = 1 : u2_ind(1)
            FHI_u2(counter) = i;
            counter = counter+1; 
            Ctrl_matrix(2,j) = uvec_u2(i);
        end 
    elseif i == m.Nj(2) % 3
        for j = u2_ind(end) + 1 : m.p + 1
            FHI_u2(counter) = i; 
            counter = counter+1;
            Ctrl_matrix(2,j) = uvec_u2(i);
        end 
    else  
        for j = u2_ind(i-1) + 1 : u2_ind(i)
            FHI_u2(counter) = i;
            counter = counter+1;
            Ctrl_matrix(2,j) = uvec_u2(i);
        end 
    end 
end 

tau_tspan = [0 sorted_vec(1:m.p)' tF];
theta_vec = zeros(1, m.p + 1);
for i = 1 : m.p + 1 % 5
    theta_vec(i) = tau_tspan(i+1) - tau_tspan(i) ;
end

end % Function Timehorizon_Sort End
%-------------------------------------------------------------------------