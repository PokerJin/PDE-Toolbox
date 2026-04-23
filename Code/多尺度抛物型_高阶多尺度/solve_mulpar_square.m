% =========================================================================
% 使用二阶渐近均匀化方法求解二维多尺度抛物型偏微分方程
% 方程: d*u' - div(c*grad(u)) = f(x,y,t)
% 版本: 内含物已修改为正方形
% =========================================================================

function solve_mulpar_square(xyt_range,Nmmtp,equ_para,u0,bc_type,bc_para)

%% 1. 用户设定参数区域

% A. 微观结构与材料属性

% --- 修改 1: 参数更新 ---
square_half_side = equ_para.len;% 正方形夹杂物的实际半边长
% --- 修改结束 ---
c_inclusion = equ_para.c_in;        % 内含物热导率 c
c_matrix = equ_para.c_ma;         % 基体热导率 c
d_inclusion = equ_para.d_in;        % 内含物系数 d
d_matrix = equ_para.d_ma;           % 基体系数 d

% B. 时空域设定
x_min = xyt_range(1);            % x 轴左边界
x_max = xyt_range(2);            % x 轴右边界
y_min = xyt_range(3);            % y 轴左边界
y_max = xyt_range(4);            % y 轴右边界
T_start = xyt_range(5);          % 求解起始时间
T_final = xyt_range(6);          % 求解结束时间

epsilon = equ_para.epsilon*(x_max-x_min);          % 微观与宏观尺度比率
% C. 初始条件
IC_func = u0;

% D. 源项
f_func = equ_para.f;

% E. 边界条件
if bc_type(1)==1
    BC.left.type    = 'Dirichlet';
else
    BC.left.type    = 'Neumann';
end
BC.left.h       = bc_para.h_left;
BC.left.r       = bc_para.r_left;
BC.left.q=bc_para.q_left;
BC.left.g=bc_para.g_left;

if bc_type(2)==1
    BC.right.type   = 'Dirichlet';
else
    BC.right.type   = 'Neumann';
end

BC.right.h      = bc_para.h_right;
BC.right.r      = bc_para.r_right;
BC.right.q=bc_para.q_right;
BC.right.g=bc_para.g_right;

if bc_type(3)==1
    BC.bottom.type  = 'Dirichlet';
else
    BC.bottom.type  = 'Neumann';
end

BC.bottom.h     = bc_para.h_bottom;
BC.bottom.r     = bc_para.r_bottom;
BC.bottom.q=bc_para.q_bottom;
BC.bottom.g=bc_para.g_bottom;

if bc_type(4)==1
    BC.top.type='Dirichlet';
else
    BC.top.type  = 'Neumann';
end
BC.top.h        = bc_para.h_top;
BC.top.r        =bc_para.r_top;
BC.top.q=bc_para.q_top;
BC.top.g=bc_para.g_top;

% F. 数值参数
N_micro = Nmmtp(1);          % 微观单胞网格点数
N_macro = Nmmtp(2);           % 宏观均匀化问题网格点数
N_t = Nmmtp(3);               % 时间步数
N_fine = Nmmtp(4);           % 最终重构解的精细网格点数

%% 程序主执行部分
domain_bounds = [x_min, x_max, y_min, y_max];
time_interval = [T_start, T_final];

% --- 修改 2: 几何定义更新 ---
% 定义微观单胞内的材料属性函数
yc = [0.5, 0.5]; % 单胞中心
% 将宏观半边长映射到微观单胞坐标
s_cell = square_half_side / epsilon; 

% 定义一个匿名函数用于判断点是否在中心正方形内
is_inside = @(y1, y2) (abs(y1 - yc(1)) <= s_cell) & (abs(y2 - yc(2)) <= s_cell);

% 使用 is_inside 函数重新定义 c_func 和 d_func
c_func = @(y1, y2) is_inside(y1, y2) * c_inclusion + ~is_inside(y1, y2) * c_matrix;
d_func = @(y1, y2) is_inside(y1, y2) * d_inclusion + ~is_inside(y1, y2) * d_matrix;
% --- 修改结束 ---


% 步骤 1: 求解单胞问题和均匀化系数
[cell_sol, hom_coeffs] = solve_cell_problems(N_micro, c_func, d_func);

% 步骤 2: 求解宏观均匀化方程
u0 = solve_homogenized_pde(N_macro, N_t, domain_bounds, time_interval, hom_coeffs, f_func, IC_func, BC);

% 步骤 3: 重构二阶多尺度近似解
u_final = reconstruct_solution(u0, cell_sol, hom_coeffs, epsilon, N_fine, domain_bounds, time_interval);

% 步骤 4: 绘制结果
x_fine = linspace(x_min, x_max, N_fine);
y_fine = linspace(y_min, y_max, N_fine);
[XF, YF] = meshgrid(x_fine, y_fine);

figure;
subplot(1,2,1);
surf(XF, YF, u_final');
xlabel('x'); ylabel('y'); zlabel('u');
title('二阶多尺度解（三维）');
shading interp; colorbar; axis tight; view(3);

subplot(1,2,2);
contourf(XF, YF, u_final', 20);
xlabel('x'); ylabel('y');
title('二阶多尺度解（二维）');
axis tight equal; colorbar;
end

%% ========================================================================
%   辅助函数
% =========================================================================

function [cell_sol, hom_coeffs] = solve_cell_problems(N_micro, c_func, d_func)
    y = linspace(0, 1, N_micro);
    hy = y(2) - y(1);
    [Y1, Y2] = meshgrid(y, y);
    c = c_func(Y1, Y2);
    d = d_func(Y1, Y2);
    dc_dy1 = (circshift(c, [0, -1]) - circshift(c, [0, 1])) / (2 * hy);
    dc_dy2 = (circshift(c, [-1, 0]) - circshift(c, [1, 0])) / (2 * hy);
    
    rhs_N1 = -dc_dy1;
    rhs_N2 = -dc_dy2;
    cell_sol.N1 = solve_elliptic_5pt(N_micro, c, rhs_N1);
    cell_sol.N2 = solve_elliptic_5pt(N_micro, c, rhs_N2);
    hom_coeffs.d = mean(d(:));
    [dN1_dy1, dN1_dy2] = gradient(cell_sol.N1, hy);
    [dN2_dy1, dN2_dy2] = gradient(cell_sol.N2, hy);
    hom_coeffs.c = zeros(2, 2);
    hom_coeffs.c(1,1) = mean(mean(c .* (1 + dN1_dy1)));
    hom_coeffs.c(1,2) = mean(mean(c .* dN2_dy1));
    hom_coeffs.c(2,1) = mean(mean(c .* dN1_dy2));
    hom_coeffs.c(2,2) = mean(mean(c .* (1 + dN2_dy2)));
    rhs_S = d - hom_coeffs.d;
    cell_sol.S = solve_elliptic_5pt(N_micro, c, rhs_S);
    rhs_M11 = hom_coeffs.c(1,1) - c .* (1 + dN1_dy1);
    rhs_M12 = hom_coeffs.c(1,2) - c .* dN2_dy1;
    rhs_M21 = hom_coeffs.c(2,1) - c .* dN1_dy2;
    rhs_M22 = hom_coeffs.c(2,2) - c .* (1 + dN2_dy2);
    
    cell_sol.M11 = solve_elliptic_5pt(N_micro, c, rhs_M11);
    cell_sol.M12 = solve_elliptic_5pt(N_micro, c, rhs_M12);
    cell_sol.M21 = solve_elliptic_5pt(N_micro, c, rhs_M21);
    cell_sol.M22 = solve_elliptic_5pt(N_micro, c, rhs_M22);
    
    cell_sol.y_grid = y;
end

function U = solve_elliptic_5pt(N, c, F)
    h = 1 / (N - 1);
    n_pts = N * N;
    
    c_vec = reshape(c, n_pts, 1);
    c_E = (c_vec + circshift(c_vec, -N)) / 2;
    c_W = (c_vec + circshift(c_vec, N)) / 2;
    c_N = (c_vec + circshift(c_vec, -1)) / 2;
    c_S = (c_vec + circshift(c_vec, 1)) / 2;
    
    A = spdiags([-c_N/h^2, -c_S/h^2, (c_E+c_W+c_N+c_S)/h^2, -c_W/h^2, -c_E/h^2], ...
                [-1, 1, 0, N, -N], n_pts, n_pts);
    
    is_bnd = false(N, N);
    is_bnd(1,:)=true; is_bnd(end,:)=true; is_bnd(:,1)=true; is_bnd(:,end)=true;
    bnd_idx = find(reshape(is_bnd, n_pts, 1));
    
    A(bnd_idx, :) = 0;
    A(bnd_idx, bnd_idx) = speye(length(bnd_idx));
    
    f_vec = reshape(F, n_pts, 1);
    f_vec(bnd_idx) = 0;
    u_vec = A \ -f_vec;
    U = reshape(u_vec, N, N);
end

function u0 = solve_homogenized_pde(N_macro, N_t, domain_bounds, time_interval, hom_coeffs, f_func, IC_func, BC)
    x_min = domain_bounds(1); x_max = domain_bounds(2);
    y_min = domain_bounds(3); y_max = domain_bounds(4);
    Lx = x_max - x_min; Ly = y_max - y_min;
    
    T_start = time_interval(1); T_final = time_interval(2);
    
    hx = Lx / (N_macro - 1); hy = Ly / (N_macro - 1);
    t_steps = linspace(T_start, T_final, N_t + 1);
    dt = t_steps(2) - t_steps(1);
    
    x_vec = linspace(x_min, x_max, N_macro);
    y_vec = linspace(y_min, y_max, N_macro);
    [X, Y] = meshgrid(x_vec, y_vec);
    n_pts = N_macro * N_macro;
    
    c11 = hom_coeffs.c(1,1); c22 = hom_coeffs.c(2,2);
    d_hom = hom_coeffs.d;
    
    L_op = c11 * (1/hx^2) * spdiags([ones(n_pts,1) -2*ones(n_pts,1) ones(n_pts,1)], [-N_macro 0 N_macro], n_pts, n_pts) + ...
           c22 * (1/hy^2) * spdiags([ones(n_pts,1) -2*ones(n_pts,1) ones(n_pts,1)], [-1 0 1], n_pts, n_pts);
    
    [~, bnd_map] = get_boundary_indices(N_macro);
    
    G_vec = zeros(n_pts, 1);
    
    I = speye(n_pts);
    A_base = d_hom/dt * I - 0.5 * L_op;
    B_base = d_hom/dt * I + 0.5 * L_op;
    
    u0 = zeros(N_macro, N_macro, N_t + 1);
    IC_values = IC_func(X, Y);
    if isscalar(IC_values); IC_values = IC_values * ones(N_macro); end
    u0(:, :, 1) = IC_values;
    u_vec = reshape(u0(:, :, 1), n_pts, 1);
    
    for i = 1:N_t
        t_now = t_steps(i); t_next = t_steps(i+1);
        A = A_base; B = B_base;
        
        F_now = f_func(X, Y, t_now); F_next = f_func(X, Y, t_next);
        if isscalar(F_now); F_now = F_now * ones(N_macro); end
        if isscalar(F_next); F_next = F_next * ones(N_macro); end
        F_avg = (reshape(F_now, n_pts, 1) + reshape(F_next, n_pts, 1)) / 2;
        
        rhs = B * u_vec + F_avg + G_vec; 
        
        dirichlet_boundaries = {'left', 'right', 'bottom', 'top'};
        for k = 1:4
            bnd_name = dirichlet_boundaries{k};
            if strcmp(BC.(bnd_name).type, 'Dirichlet')
                indices = bnd_map.(bnd_name);
                x_coords = X(indices); y_coords = Y(indices);
                h_vals = BC.(bnd_name).h(x_coords, y_coords, t_next);
                r_vals = BC.(bnd_name).r(x_coords, y_coords, t_next);
                if isscalar(h_vals); h_vals = h_vals * ones(size(indices)); end
                if isscalar(r_vals); r_vals = r_vals * ones(size(indices)); end
                for j=1:length(indices)
                    idx = indices(j);
                    A(idx, :) = 0; A(idx, idx) = h_vals(j);
                    rhs(idx) = r_vals(j);
                end
            end
        end
        u_vec = A \ rhs;
        u0(:, :, i+1) = reshape(u_vec, N_macro, N_macro);
    end
end

function u_final = reconstruct_solution(u0, cell_sol, hom_coeffs, epsilon, N_fine, domain_bounds, time_interval)
    x_min = domain_bounds(1); x_max = domain_bounds(2);
    y_min = domain_bounds(3); y_max = domain_bounds(4);
    Lx = x_max - x_min; Ly = y_max - y_min;
    T_start = time_interval(1); T_final = time_interval(2);
    
    [N_macro, ~, N_t_plus_1] = size(u0);
    N_t = N_t_plus_1 - 1;
    hx = Lx / (N_macro - 1); hy = Ly / (N_macro - 1);
    dt = (T_final - T_start) / max(1, N_t);
    x_macro_vec = linspace(x_min, x_max, N_macro);
    y_macro_vec = linspace(y_min, y_max, N_macro);
    [X_macro, Y_macro] = meshgrid(x_macro_vec, y_macro_vec);
    
    u0_final = u0(:, :, end);
    u0_prev = (N_t > 0) * u0(:, :, end-1) + (N_t == 0) * u0_final;

    [u0_dx, u0_dy] = gradient(u0_final, hx, hy);
    [u0_dxdx, u0_dydx] = gradient(u0_dx, hx, hy);
    [u0_dxdy, u0_dydy] = gradient(u0_dy, hx, hy);
    
    if (N_t >= 2)
        u0_prev_2 = u0(:, :, end-2);
        u0_dt = (3*u0_final - 4*u0_prev + u0_prev_2) / (2 * dt);
    else
        u0_dt = (u0_final - u0_prev) / dt;
    end
    
    x_fine_vec = linspace(x_min, x_max, N_fine);
    y_fine_vec = linspace(y_min, y_max, N_fine);
    [XF, YF] = meshgrid(x_fine_vec, y_fine_vec); 
    U0_interp      = interp2_manual(X_macro, Y_macro, u0_final,  XF, YF, 'spline');
    U0_dx_interp   = interp2_manual(X_macro, Y_macro, u0_dx,     XF, YF, 'spline');
    U0_dy_interp   = interp2_manual(X_macro, Y_macro, u0_dy,     XF, YF, 'spline');
    U0_dxdx_interp = interp2_manual(X_macro, Y_macro, u0_dxdx,   XF, YF, 'spline');
    U0_dydy_interp = interp2_manual(X_macro, Y_macro, u0_dydy,   XF, YF, 'spline');
    U0_dxdy_interp = interp2_manual(X_macro, Y_macro, u0_dxdy,   XF, YF, 'spline');
    U0_dt_interp   = interp2_manual(X_macro, Y_macro, u0_dt,     XF, YF, 'spline');
    
    Y_coord1 = mod(XF / epsilon, 1);
    Y_coord2 = mod(YF / epsilon, 1);
    
    x_micro_vec = cell_sol.y_grid;
    y_micro_vec = cell_sol.y_grid;
    [X_micro_grid, Y_micro_grid] = meshgrid(x_micro_vec, y_micro_vec);
    N1_interp = interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.N1, Y_coord1, Y_coord2, 'spline');
    N2_interp = interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.N2, Y_coord1, Y_coord2, 'spline');
    S_interp  = interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.S,  Y_coord1, Y_coord2, 'spline');
    M11_interp= interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.M11,Y_coord1, Y_coord2, 'spline');
    M12_interp= interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.M12,Y_coord1, Y_coord2, 'spline');
    M21_interp= interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.M21,Y_coord1, Y_coord2, 'spline');
    M22_interp= interp2_manual(X_micro_grid, Y_micro_grid, cell_sol.M22,Y_coord1, Y_coord2, 'spline');
    
    u1 = N1_interp .* U0_dx_interp + N2_interp .* U0_dy_interp;
    u2 = S_interp .* U0_dt_interp + M11_interp .* U0_dxdx_interp ...
         + (M12_interp + M21_interp) .* U0_dxdy_interp + M22_interp .* U0_dydy_interp;
    u_final = U0_interp + epsilon * u1 + epsilon^2 * u2;
end

function [bnd_idx, bnd_map] = get_boundary_indices(N)
    is_bnd = false(N, N);
    is_bnd(1,:)=true; is_bnd(end,:)=true; is_bnd(:,1)=true; is_bnd(:,end)=true;
    bnd_idx = find(reshape(is_bnd, N*N, 1));
    map = zeros(N,N);
    map(:,1) = 1; map(:,N) = 2; map(1,:) = 3; map(N,:) = 4;
    map(1,1)=5; map(N,1)=6; map(1,N)=7; map(N,N)=8;
    bnd_map.left   = find(map==1 | map==5 | map==6);
    bnd_map.right  = find(map==2 | map==7 | map==8);
    bnd_map.bottom = find(map==3 | map==5 | map==7);
    bnd_map.top    = find(map==4 | map==6 | map==8);
end

function Vq = interp2_manual(X, Y, V, Xq, Yq, method)
    x_vec = X(1,:);
    y_vec = Y(:,1);

    xq_vec = Xq(1,:);
    yq_vec = Yq(:,1);

    V_intermediate = zeros(length(y_vec), length(xq_vec));
    for i = 1:length(y_vec)
        V_intermediate(i, :) = interp1(x_vec, V(i, :), xq_vec, method);
    end

    Vq = zeros(length(yq_vec), length(xq_vec));
    for j = 1:length(xq_vec)
        Vq(:, j) = interp1(y_vec, V_intermediate(:, j), yq_vec, method);
    end
end