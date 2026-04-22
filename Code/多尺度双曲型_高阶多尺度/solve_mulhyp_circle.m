%% 二维多尺度双曲型方程的二阶渐近均匀化方法求解
% 方程: d*u'' - div(c*grad(u)) = f


function solve_mulhyp_circle(xyt_range,Nmmtp,equ_para,init,bc_type,bc_para)

%% ========== 用户可配置参数 ==========

% 空间区域设置
x_left = xyt_range(1);              % x方向左边界
x_right = xyt_range(2);             % x方向右边界
y_left = xyt_range(3);              % y方向左边界
y_right = xyt_range(4);             % y方向右边界

% 时间区域设置
t_start = xyt_range(5);             % 起始时间
t_end = xyt_range(6);             % 终止时间

% 网格设置
n_micro = Nmmtp(1);            % 微观网格点数（单胞）
n_macro = Nmmtp(2);            % 宏观网格点数
n_plot = Nmmtp(4);            % 细网格点数（用于绘图）
n_time = Nmmtp(3);            % 时间步数

% 尺度参数
epsilon = equ_para.epsilon*(x_right-x_left);           % 微观/宏观尺度比率

% 几何参数
R_circle = equ_para.len;         % 圆的实际半径

% 材料参数
c_in =equ_para.c_in;                % 圆内扩散系数
c_out = equ_para.c_ma;             % 圆外扩散系数
d_in = equ_para.d_in;                % 圆内密度
d_out = equ_para.d_ma;               % 圆外密度

%% 源项设置

f_source = equ_para.f;


%% 初始条件设置

u_initial = init.u0;

v_initial = init.v0;


%% 边界条件设置（统一使用函数句柄形式）
% 边界类型：'Dirichlet' 或 'Neumann'
% Dirichlet: h*u = r
% Neumann: n*c*grad(u) + q*u = g

% 左边界 (x = x_left)
if bc_type(1)==1
    BC_left.type = 'Dirichlet';
else
    BC_left.type='Neumann';
end

BC_left.h = bc_para.h_left;        
BC_left.r = bc_para.r_left;
BC_left.q=bc_para.q_left;
BC_left.g=bc_para.g_left;


% 右边界 (x = x_right)
if bc_type(2)==1
    BC_right.type = 'Dirichlet';
else
    BC_right.type = 'Neumann';
end

BC_right.h = bc_para.h_right;
BC_right.r = bc_para.r_right;
BC_right.q=bc_para.q_right;
BC_right.g=bc_para.g_right;

% 下边界 (y = y_left)
if bc_type(3)==1
    BC_bottom.type = 'Dirichlet';
else
    BC_bottom.type = 'Neumann';
end

BC_bottom.h = bc_para.h_bottom;
BC_bottom.r = bc_para.r_bottom;
BC_bottom.q=bc_para.q_bottom;
BC_bottom.g=bc_para.g_bottom;


% 上边界 (y = y_right)
if bc_type(4)==1
    BC_top.type = 'Dirichlet';
else
    BC_top.type = 'Neumann';
end

BC_top.h = bc_para.h_top;
BC_top.r = bc_para.r_top;
BC_top.q=bc_para.q_top;
BC_top.g=bc_para.g_top;


%% ========== 参数检查和计算 ==========
% 计算时间步长
T = t_end - t_start;
dt = T / n_time;
nt = n_time;

% 计算宏观区域尺寸
L_x = x_right - x_left;
L_y = y_right - y_left;

%% ========== 第一步：定义单胞材料参数 ==========

h_micro = 1/(n_micro-1);
[y1, y2] = meshgrid(linspace(0, 1, n_micro));

R_cell = R_circle/epsilon; 
dist_cell = sqrt((y1-0.5).^2 + (y2-0.5).^2);
c_cell = c_out * ones(n_micro, n_micro);
c_cell(dist_cell <= R_cell) = c_in;
d_cell = d_out * ones(n_micro, n_micro);
d_cell(dist_cell <= R_cell) = d_in;


%% ========== 第二步：求解一阶单胞问题 ==========

% 对于2D问题，需要求解两个单胞函数：N1 (对应x1方向) 和 N2 (对应x2方向)
N1_cell = solve_cell_problem_order1(c_cell, h_micro, n_micro, 1);
N2_cell = solve_cell_problem_order1(c_cell, h_micro, n_micro, 2);

%% ========== 第三步：计算均匀化系数 ==========

[dN1_dy1, dN1_dy2] = compute_gradient(N1_cell, h_micro);
[dN2_dy1, dN2_dy2] = compute_gradient(N2_cell, h_micro);

d_bar = mean(d_cell(:));

a_bar_11 = mean(mean(c_cell .* (1 + dN1_dy1)));
a_bar_22 = mean(mean(c_cell .* (1 + dN2_dy2)));
a_bar_12 = mean(mean(c_cell .* dN2_dy1));
a_bar_21 = mean(mean(c_cell .* dN1_dy2));


%% ========== 第四步：求解二阶单胞问题 ==========

N11_cell = solve_cell_problem_order2_N(c_cell, N1_cell, N1_cell, ...
    a_bar_11, h_micro, n_micro, 1, 1);

N12_cell = solve_cell_problem_order2_N(c_cell, N1_cell, N2_cell, ...
    a_bar_12, h_micro, n_micro, 1, 2);
N22_cell = solve_cell_problem_order2_N(c_cell, N2_cell, N2_cell, ...
    a_bar_22, h_micro, n_micro, 2, 2);

F_cell = solve_cell_problem_F(c_cell, d_cell, d_bar, h_micro, n_micro);

G_cell = zeros(n_micro, n_micro);

%% ========== 第五步：求解均匀化双曲型方程 ==========


% 宏观网格
h_macro_x = L_x/(n_macro-1);
h_macro_y = L_y/(n_macro-1);
[x1, x2] = meshgrid(linspace(x_left, x_right, n_macro), ...
                    linspace(y_left, y_right, n_macro));

% 处理初始条件
u0 = u_initial(x1, x2);
v0 = v_initial(x1, x2);


% 求解均匀化方程
[u_homo, v_homo, a_homo] = solve_homogenized_hyperbolic(...
    u0, v0, f_source, d_bar, a_bar_11, a_bar_22, a_bar_12, ...
    h_macro_x, h_macro_y, dt, nt, n_macro, x1, x2, t_start, ...
    BC_left, BC_right, BC_bottom, BC_top);


%% ========== 第六步：构造二阶双尺度解 ==========

h_plot = max(L_x, L_y)/(n_plot-1);
[X1, X2] = meshgrid(linspace(x_left, x_right, n_plot), ...
                    linspace(y_left, y_right, n_plot));

u_homo_fine = interp2_manual(x1, x2, u_homo, X1, X2, 'linear');

[du_dx1, du_dx2] = compute_gradient_2d(u_homo, h_macro_x, h_macro_y);
du_dx1_fine = interp2_manual(x1, x2, du_dx1, X1, X2, 'linear');
du_dx2_fine = interp2_manual(x1, x2, du_dx2, X1, X2, 'linear');

d2u_dx1dx1 = compute_second_derivative(u_homo, h_macro_x, 1, 1);
d2u_dx1dx2 = compute_mixed_derivative(u_homo, h_macro_x, h_macro_y);
d2u_dx2dx2 = compute_second_derivative(u_homo, h_macro_y, 2, 2);

d2u_dx1dx1_fine = interp2_manual(x1, x2, d2u_dx1dx1, X1, X2, 'linear');
d2u_dx1dx2_fine = interp2_manual(x1, x2, d2u_dx1dx2, X1, X2, 'linear');
d2u_dx2dx2_fine = interp2_manual(x1, x2, d2u_dx2dx2, X1, X2, 'linear');

v_homo_fine = interp2_manual(x1, x2, v_homo, X1, X2, 'linear');
a_homo_fine = interp2_manual(x1, x2, a_homo, X1, X2, 'linear');

Y1_local = mod(X1/epsilon, 1);
Y2_local = mod(X2/epsilon, 1);

N1_vals = interp2_manual(y1, y2, N1_cell, Y1_local, Y2_local, 'linear', 0);
N2_vals = interp2_manual(y1, y2, N2_cell, Y1_local, Y2_local, 'linear', 0);
N11_vals = interp2_manual(y1, y2, N11_cell, Y1_local, Y2_local, 'linear', 0);
N12_vals = interp2_manual(y1, y2, N12_cell, Y1_local, Y2_local, 'linear', 0);
N22_vals = interp2_manual(y1, y2, N22_cell, Y1_local, Y2_local, 'linear', 0);
F_vals = interp2_manual(y1, y2, F_cell, Y1_local, Y2_local, 'linear', 0);

term1 = epsilon * (N1_vals .* du_dx1_fine + N2_vals .* du_dx2_fine);

term2 = epsilon^2 * (N11_vals .* d2u_dx1dx1_fine + ...
                     N12_vals .* d2u_dx1dx2_fine + ...
                     N22_vals .* d2u_dx2dx2_fine + ...
                     F_vals .* a_homo_fine);

u_2scale = u_homo_fine + term1 + term2;


%% ========== 第七步：可视化结果 ==========

figure;

% 三维图
subplot(1, 2, 1);
surf(X1, X2, u_2scale);
xlabel('x'); ylabel('y'); zlabel('u');
title('二阶多尺度解（三维）');
colorbar;
view(3);
axis tight;
shading interp;

% 二维等高线图
subplot(1, 2, 2);
contourf(X1, X2, u_2scale, 20);
xlabel('x'); ylabel('y');
title('二阶多尺度解（二维）');
colorbar;
axis tight equal;
end

%% ========== 辅助函数 ==========

% 求解一阶单胞问题（五点差分格式）
function N = solve_cell_problem_order1(c_cell, h, n, direction)

    N_size = (n-2)^2;
    A = sparse(N_size, N_size);
    b = zeros(N_size, 1);
    
    % 内部节点编号
    node_map = zeros(n, n);
    idx = 1;
    for j = 2:n-1
        for i = 2:n-1
            node_map(i, j) = idx;
            idx = idx + 1;
        end
    end
    
    for j = 2:n-1
        for i = 2:n-1
            row = node_map(i, j);
            
            c_center = c_cell(i, j);
            c_right = 0.5 * (c_cell(i, j) + c_cell(i+1, j));
            c_left = 0.5 * (c_cell(i, j) + c_cell(i-1, j));
            c_up = 0.5 * (c_cell(i, j) + c_cell(i, j+1));
            c_down = 0.5 * (c_cell(i, j) + c_cell(i, j-1));
            
            coef_center = -(c_right + c_left + c_up + c_down) / h^2;
            A(row, row) = coef_center;
            
            if i < n-1
                col = node_map(i+1, j);
                A(row, col) = c_right / h^2;
            end
            
            if i > 2
                col = node_map(i-1, j);
                A(row, col) = c_left / h^2;
            end
            
            if j < n-1
                col = node_map(i, j+1);
                A(row, col) = c_up / h^2;
            end
            
            if j > 2
                col = node_map(i, j-1);
                A(row, col) = c_down / h^2;
            end
            
            if direction == 1
                b(row) = -(c_right - c_left) / (2*h);
            else
                b(row) = -(c_up - c_down) / (2*h);
            end
        end
    end
    
    N_vec = A \ b;
    
    N = zeros(n, n);
    for j = 2:n-1
        for i = 2:n-1
            N(i, j) = N_vec(node_map(i, j));
        end
    end
end

% 求解二阶单胞问题 N_αβ
function N_ab = solve_cell_problem_order2_N(c_cell, N_a, N_b, a_bar_ab, h, n, dir_a, dir_b)
    N_size = (n-2)^2;
    A = sparse(N_size, N_size);
    b = zeros(N_size, 1);
    
    node_map = zeros(n, n);
    idx = 1;
    for j = 2:n-1
        for i = 2:n-1
            node_map(i, j) = idx;
            idx = idx + 1;
        end
    end
    
    [dNa_dy1, dNa_dy2] = compute_gradient(N_a, h);
    [dNb_dy1, dNb_dy2] = compute_gradient(N_b, h);
    
    for j = 2:n-1
        for i = 2:n-1
            row = node_map(i, j);
            
            c_center = c_cell(i, j);
            c_right = 0.5 * (c_cell(i, j) + c_cell(i+1, j));
            c_left = 0.5 * (c_cell(i, j) + c_cell(i-1, j));
            c_up = 0.5 * (c_cell(i, j) + c_cell(i, j+1));
            c_down = 0.5 * (c_cell(i, j) + c_cell(i, j-1));
            
            coef_center = -(c_right + c_left + c_up + c_down) / h^2;
            A(row, row) = coef_center;
            
            if i < n-1
                A(row, node_map(i+1, j)) = c_right / h^2;
            end
            if i > 2
                A(row, node_map(i-1, j)) = c_left / h^2;
            end
            if j < n-1
                A(row, node_map(i, j+1)) = c_up / h^2;
            end
            if j > 2
                A(row, node_map(i, j-1)) = c_down / h^2;
            end

            delta_ab = (dir_a == dir_b);

            if dir_a == 1
                dNb_dya = dNb_dy1(i, j); 
            else
                dNb_dya = dNb_dy2(i, j);
            end

            if dir_b == 1
                div_term = (c_cell(i+1,j)*N_a(i+1,j) - c_cell(i-1,j)*N_a(i-1,j))/(2*h);
            else
                div_term = (c_cell(i,j+1)*N_a(i,j+1) - c_cell(i,j-1)*N_a(i,j-1))/(2*h);
            end

            b(row) = a_bar_ab - c_center*(delta_ab + dNb_dya) - div_term;
        end
    end

    N_vec = A \ b;

    N_ab = zeros(n, n);
    for j = 2:n-1
        for i = 2:n-1
            N_ab(i, j) = N_vec(node_map(i, j));
        end
    end
end

% 求解二阶单胞问题 F
function F = solve_cell_problem_F(c_cell, d_cell, d_bar, h, n)
    N_size = (n-2)^2;
    A = sparse(N_size, N_size);
    b = zeros(N_size, 1);
    
    node_map = zeros(n, n);
    idx = 1;
    for j = 2:n-1
        for i = 2:n-1
            node_map(i, j) = idx;
            idx = idx + 1;
        end
    end
    
    for j = 2:n-1
        for i = 2:n-1
            row = node_map(i, j);
            
            c_right = 0.5 * (c_cell(i, j) + c_cell(i+1, j));
            c_left = 0.5 * (c_cell(i, j) + c_cell(i-1, j));
            c_up = 0.5 * (c_cell(i, j) + c_cell(i, j+1));
            c_down = 0.5 * (c_cell(i, j) + c_cell(i, j-1));
            
            coef_center = -(c_right + c_left + c_up + c_down) / h^2;
            A(row, row) = coef_center;
            
            if i < n-1, A(row, node_map(i+1, j)) = c_right / h^2; end
            if i > 2, A(row, node_map(i-1, j)) = c_left / h^2; end
            if j < n-1, A(row, node_map(i, j+1)) = c_up / h^2; end
            if j > 2, A(row, node_map(i, j-1)) = c_down / h^2; end

            b(row) = d_cell(i, j) - d_bar;
        end
    end
    
    F_vec = A \ b;
    
    F = zeros(n, n);
    for j = 2:n-1
        for i = 2:n-1
            F(i, j) = F_vec(node_map(i, j));
        end
    end
end

function [du_dx1, du_dx2] = compute_gradient(u, h)
    [n1, n2] = size(u);
    du_dx1 = zeros(n1, n2);
    du_dx2 = zeros(n1, n2);

    for i = 2:n1-1
        for j = 2:n2-1
            du_dx1(i, j) = (u(i+1, j) - u(i-1, j)) / (2*h);
            du_dx2(i, j) = (u(i, j+1) - u(i, j-1)) / (2*h);
        end
    end

    du_dx1(1, :) = (u(2, :) - u(1, :)) / h;
    du_dx1(n1, :) = (u(n1, :) - u(n1-1, :)) / h;
    du_dx2(:, 1) = (u(:, 2) - u(:, 1)) / h;
    du_dx2(:, n2) = (u(:, n2) - u(:, n2-1)) / h;
end

function [du_dx1, du_dx2] = compute_gradient_2d(u, h_x, h_y)
    [du_dx2, du_dx1] = gradient(u, h_y, h_x);
end

function d2u = compute_second_derivative(u, h, dir1, dir2)
    [n1, n2] = size(u);
    d2u = zeros(n1, n2);
    
    if dir1 == 1 && dir2 == 1
        d2u(2:n1-1, :) = (u(3:n1, :) - 2*u(2:n1-1, :) + u(1:n1-2, :)) / h^2;
    elseif dir1 == 2 && dir2 == 2
        d2u(:, 2:n2-1) = (u(:, 3:n2) - 2*u(:, 2:n2-1) + u(:, 1:n2-2)) / h^2;
    end
end

function d2u = compute_mixed_derivative(u, h_x, h_y)
    [n1, n2] = size(u);
    d2u = zeros(n1, n2);

    d2u(2:n1-1, 2:n2-1) = (u(3:n1, 3:n2) - u(3:n1, 1:n2-2) - ...
                           u(1:n1-2, 3:n2) + u(1:n1-2, 1:n2-2)) / (4*h_x*h_y);
end

function [u, v, a] = solve_homogenized_hyperbolic(u0, v0, f_source, d_bar, ...
    a11, a22, a12, h_x, h_y, dt, nt, n, x1, x2, t_start, ...
    BC_left, BC_right, BC_bottom, BC_top)

    u = u0;        
    v = v0;        

    u = apply_boundary_conditions(u, n, x1, x2, t_start, ...
        BC_left, BC_right, BC_bottom, BC_top, h_x, h_y);

    f = f_source(x1, x2, t_start);
    if isscalar(f)
        f = f * ones(size(x1));
    end

    a = compute_acceleration(u, f, d_bar, a11, a22, a12, h_x, h_y, n, ...
        BC_left, BC_right, BC_bottom, BC_top, t_start, x1, x2);
    
    for step = 1:nt
        current_time = t_start + step * dt;
        
        f = f_source(x1, x2, current_time);
        if isscalar(f)
            f = f * ones(size(x1));
        end

        u_new = u + dt*v + 0.5*dt^2*a;
        
        u_new = apply_boundary_conditions(u_new, n, x1, x2, current_time, ...
            BC_left, BC_right, BC_bottom, BC_top, h_x, h_y);

        a_new = compute_acceleration(u_new, f, d_bar, a11, a22, a12, h_x, h_y, n, ...
            BC_left, BC_right, BC_bottom, BC_top, current_time, x1, x2);

        v_new = v + 0.5*dt*(a + a_new);

        v_new = apply_velocity_bc(v_new, n, BC_left, BC_right, BC_bottom, BC_top);

        u = u_new;
        v = v_new;
        a = a_new;
    end
end

function u = apply_boundary_conditions(u, n, x1, x2, current_time, ...
    BC_left, BC_right, BC_bottom, BC_top, h_x, h_y)

    persistent BC_cached BC_cache_valid

    if isempty(BC_cache_valid)
        BC_cache_valid = false;
    end
    
    if ~BC_cache_valid
        BC_cached = struct();
        
        % 检查左边界
        if strcmp(BC_left.type, 'Dirichlet')
            if ~check_time_dependency(BC_left.h, x1(1,:), x2(1,:)) && ...
               ~check_time_dependency(BC_left.r, x1(1,:), x2(1,:))
                BC_cached.left_h = BC_left.h(x1(1,:), x2(1,:), 0);
                BC_cached.left_r = BC_left.r(x1(1,:), x2(1,:), 0);
                BC_cached.left_time_indep = true;
            else
                BC_cached.left_time_indep = false;
            end
        else
            BC_cached.left_time_indep = false;
        end
        
        % 检查右边界
        if strcmp(BC_right.type, 'Dirichlet')
            if ~check_time_dependency(BC_right.h, x1(n,:), x2(n,:)) && ...
               ~check_time_dependency(BC_right.r, x1(n,:), x2(n,:))
                BC_cached.right_h = BC_right.h(x1(n,:), x2(n,:), 0);
                BC_cached.right_r = BC_right.r(x1(n,:), x2(n,:), 0);
                BC_cached.right_time_indep = true;
            else
                BC_cached.right_time_indep = false;
            end
        else
            BC_cached.right_time_indep = false;
        end
        
        % 检查下边界
        if strcmp(BC_bottom.type, 'Dirichlet')
            if ~check_time_dependency(BC_bottom.h, x1(:,1), x2(:,1)) && ...
               ~check_time_dependency(BC_bottom.r, x1(:,1), x2(:,1))
                BC_cached.bottom_h = BC_bottom.h(x1(:,1), x2(:,1), 0);
                BC_cached.bottom_r = BC_bottom.r(x1(:,1), x2(:,1), 0);
                BC_cached.bottom_time_indep = true;
            else
                BC_cached.bottom_time_indep = false;
            end
        else
            BC_cached.bottom_time_indep = false;
        end
        
        % 检查上边界
        if strcmp(BC_top.type, 'Dirichlet')
            if ~check_time_dependency(BC_top.h, x1(:,n), x2(:,n)) && ...
               ~check_time_dependency(BC_top.r, x1(:,n), x2(:,n))
                BC_cached.top_h = BC_top.h(x1(:,n), x2(:,n), 0);
                BC_cached.top_r = BC_top.r(x1(:,n), x2(:,n), 0);
                BC_cached.top_time_indep = true;
            else
                BC_cached.top_time_indep = false;
            end
        else
            BC_cached.top_time_indep = false;
        end
        
        BC_cache_valid = true;
    end
    
    % 左边界
    if strcmp(BC_left.type, 'Dirichlet')
        if BC_cached.left_time_indep
            u(1, :) = BC_cached.left_r ./ BC_cached.left_h;
        else
            h_val = BC_left.h(x1(1,:), x2(1,:), current_time);
            r_val = BC_left.r(x1(1,:), x2(1,:), current_time);
            u(1, :) = r_val ./ h_val;
        end
    elseif strcmp(BC_left.type, 'Neumann')
        q_val = BC_left.q(x1(1,:), x2(1,:), current_time);
        g_val = BC_left.g(x1(1,:), x2(1,:), current_time);
        c_eff = 1;  % 使用均匀化系数
        u(1, :) = (g_val + c_eff*u(2,:)/h_x) ./ (q_val + c_eff/h_x);
    end
    
    % 右边界
    if strcmp(BC_right.type, 'Dirichlet')
        if BC_cached.right_time_indep
            u(n, :) = BC_cached.right_r ./ BC_cached.right_h;
        else
            h_val = BC_right.h(x1(n,:), x2(n,:), current_time);
            r_val = BC_right.r(x1(n,:), x2(n,:), current_time);
            u(n, :) = r_val ./ h_val;
        end
    elseif strcmp(BC_right.type, 'Neumann')
        q_val = BC_right.q(x1(n,:), x2(n,:), current_time);
        g_val = BC_right.g(x1(n,:), x2(n,:), current_time);
        c_eff = 1;
        u(n, :) = (g_val + c_eff*u(n-1,:)/h_x) ./ (q_val + c_eff/h_x);
    end
    
    % 下边界
    if strcmp(BC_bottom.type, 'Dirichlet')
        if BC_cached.bottom_time_indep
            u(:, 1) = BC_cached.bottom_r ./ BC_cached.bottom_h;
        else
            h_val = BC_bottom.h(x1(:,1), x2(:,1), current_time);
            r_val = BC_bottom.r(x1(:,1), x2(:,1), current_time);
            u(:, 1) = r_val ./ h_val;
        end
    elseif strcmp(BC_bottom.type, 'Neumann')
        q_val = BC_bottom.q(x1(:,1), x2(:,1), current_time);
        g_val = BC_bottom.g(x1(:,1), x2(:,1), current_time);
        c_eff = 1;
        u(:, 1) = (g_val + c_eff*u(:,2)/h_y) ./ (q_val + c_eff/h_y);
    end
    
    % 上边界
    if strcmp(BC_top.type, 'Dirichlet')
        if BC_cached.top_time_indep
            u(:, n) = BC_cached.top_r ./ BC_cached.top_h;
        else
            h_val = BC_top.h(x1(:,n), x2(:,n), current_time);
            r_val = BC_top.r(x1(:,n), x2(:,n), current_time);
            u(:, n) = r_val ./ h_val;
        end
    elseif strcmp(BC_top.type, 'Neumann')
        q_val = BC_top.q(x1(:,n), x2(:,n), current_time);
        g_val = BC_top.g(x1(:,n), x2(:,n), current_time);
        c_eff = 1;
        u(:, n) = (g_val + c_eff*u(:,n-1)/h_y) ./ (q_val + c_eff/h_y);
    end
end

% 应用速度边界条件
function v = apply_velocity_bc(v, n, BC_left, BC_right, BC_bottom, BC_top)
    if strcmp(BC_left.type, 'Dirichlet')
        v(1, :) = 0;
    end
    if strcmp(BC_right.type, 'Dirichlet')
        v(n, :) = 0;
    end
    if strcmp(BC_bottom.type, 'Dirichlet')
        v(:, 1) = 0;
    end
    if strcmp(BC_top.type, 'Dirichlet')
        v(:, n) = 0;
    end
end

% 计算加速度 - 优化版本（向量化）
function a = compute_acceleration(u, f, d_bar, a11, a22, a12, h_x, h_y, n, ...
    BC_left, BC_right, BC_bottom, BC_top, current_time, x1, x2)
    
    a = zeros(n, n);

    i_inner = 2:n-1;
    j_inner = 2:n-1;

    d2u_dx1dx1 = zeros(n, n);
    d2u_dx1dx1(i_inner, :) = (u(i_inner+1, :) - 2*u(i_inner, :) + u(i_inner-1, :)) / h_x^2;

    d2u_dx2dx2 = zeros(n, n);
    d2u_dx2dx2(:, j_inner) = (u(:, j_inner+1) - 2*u(:, j_inner) + u(:, j_inner-1)) / h_y^2;

    d2u_dx1dx2 = zeros(n, n);
    d2u_dx1dx2(i_inner, j_inner) = (u(i_inner+1, j_inner+1) - u(i_inner+1, j_inner-1) - ...
                                     u(i_inner-1, j_inner+1) + u(i_inner-1, j_inner-1)) / (4*h_x*h_y);

    laplacian_u = a11*d2u_dx1dx1 + a22*d2u_dx2dx2 + 2*a12*d2u_dx1dx2;

    a(i_inner, j_inner) = (laplacian_u(i_inner, j_inner) + f(i_inner, j_inner)) / d_bar;

    if strcmp(BC_left.type, 'Dirichlet'), a(1, :) = 0; end
    if strcmp(BC_right.type, 'Dirichlet'), a(n, :) = 0; end
    if strcmp(BC_bottom.type, 'Dirichlet'), a(:, 1) = 0; end
    if strcmp(BC_top.type, 'Dirichlet'), a(:, n) = 0; end
end

function is_time_dependent = check_time_dependency(func, x_vals, y_vals)
    try
        val1 = func(x_vals, y_vals, 0);
        val2 = func(x_vals, y_vals, 1);
        is_time_dependent = ~isequal(val1, val2);
    catch
        is_time_dependent = true;
    end
end

function Vq = interp2_manual(X, Y, V, Xq, Yq, method, extrapval)
    
    if nargin < 7
        use_extrapval = false;
    else
        use_extrapval = true;
    end
    
    x_vec = X(1,:);
    y_vec = Y(:,1);
    
    xq_vec = Xq(1,:);
    yq_vec = Yq(:,1);
    
    V_intermediate = zeros(length(y_vec), length(xq_vec));
    for i = 1:length(y_vec)
        if use_extrapval
            V_intermediate(i, :) = interp1(x_vec, V(i, :), xq_vec, method, extrapval);
        else
            V_intermediate(i, :) = interp1(x_vec, V(i, :), xq_vec, method);
        end
    end
    
    Vq = zeros(length(yq_vec), length(xq_vec));
    for j = 1:length(xq_vec)
        if use_extrapval
            Vq(:, j) = interp1(y_vec, V_intermediate(:, j), yq_vec, method, extrapval);
        else
            Vq(:, j) = interp1(y_vec, V_intermediate(:, j), yq_vec, method);
        end
    end
end