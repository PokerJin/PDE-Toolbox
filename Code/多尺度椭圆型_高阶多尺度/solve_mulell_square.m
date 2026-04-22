%%-div(c*grad(u)) = f
% 此脚本使用二阶渐近均匀化方法(SOTS)和有限差分法求解二维多尺度椭圆型偏微分方程。
% 版本：内含物为正方形

function solve_mulell_square(xy_range,Nmmp,equ_para,bc_type,bc_para)

%% 1. 用户配置区
% --- 求解域设置 ---
x_left = xy_range(1);
x_right = xy_range(2);
y_left = xy_range(3);
y_right = xy_range(4);

% --- 源项 f(x, y) 设置 ---
f_handle =equ_para.f;

% --- 材料和几何参数设置 ---
epsilon = equ_para.epsilon*(x_right-x_left);           
square_side_length = equ_para.len/epsilon; % 正方形边长 (在[0,1]的单胞内)
c_inclusion = equ_para.c_in;         
c_matrix = equ_para.c_ma;          

% --- 网格点数设置 ---
N_micro = Nmmp(1);           % 细观单胞问题网格
N_macro = Nmmp(2);            % 宏观均匀化问题网格
N_fine = Nmmp(3);            % 最终重构解的精细网格

% --- 宏观问题边界条件配置 ---
% BC类型: 'Dirichlet' 或 'Neumann'
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

%% 2. 步骤一: 求解一阶单胞问题
hy_micro = 1 / N_micro;
hx_micro = 1 / N_micro;
y1_micro_vec = linspace(0, 1, N_micro + 1);
y2_micro_vec = linspace(0, 1, N_micro + 1);
[Y1_micro, Y2_micro] = meshgrid(y1_micro_vec, y2_micro_vec);

% ==================== 代码修改部分 (开始) ====================
% 定义正方形内含物
c_micro = c_matrix * ones(size(Y1_micro));
half_side = square_side_length / 2;
is_inside_square = (abs(Y1_micro - 0.5) <= half_side) & (abs(Y2_micro - 0.5) <= half_side);
c_micro(is_inside_square) = c_inclusion;
% ==================== 代码修改部分 (结束) ====================

[dc_dy1, dc_dy2] = gradient_fd(c_micro, hx_micro, hy_micro);

% 求解一阶辅助函数 N1 和 N2
N1 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, dc_dy1);
N2 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, dc_dy2);

%% 3. 步骤二: 计算均匀化系数 c_H
[dN1_dy1, dN1_dy2] = gradient_fd(N1, hx_micro, hy_micro);
[dN2_dy1, dN2_dy2] = gradient_fd(N2, hx_micro, hy_micro);

integrand_c11 = c_micro .* (1 + dN1_dy1);
integrand_c12 = c_micro .* dN2_dy1;
integrand_c21 = c_micro .* dN1_dy2;
integrand_c22 = c_micro .* (1 + dN2_dy2);

cH11 = mean(integrand_c11(:));
cH12 = mean(integrand_c12(:));
cH21 = mean(integrand_c21(:));
cH22 = mean(integrand_c22(:));
c_H = [cH11, cH12; cH21, cH22];

%% 4. 步骤三: 求解均匀化问题
Lx = x_right - x_left;
Ly = y_right - y_left;
hx_macro = Lx / N_macro;
hy_macro = Ly / N_macro;
x1_macro_vec = linspace(x_left, x_right, N_macro + 1);
x2_macro_vec = linspace(y_left, y_right, N_macro + 1);
[X1_macro, X2_macro] = meshgrid(x1_macro_vec, x2_macro_vec);

f_output = f_handle(X1_macro, X2_macro);
if isscalar(f_output)
    F_macro = f_output * ones(size(X1_macro));
else
    F_macro = f_output;
end

C11_macro = ones(N_macro + 1, N_macro + 1) * c_H(1,1);
C22_macro = ones(N_macro + 1, N_macro + 1) * c_H(2,2);

% --- 边界条件处理 ---
BC.x_vec = x1_macro_vec;
BC.y_vec = x2_macro_vec;
fields = {'left', 'right', 'bottom', 'top'};
for i = 1:length(fields)
    bnd_name = fields{i};
    if strcmpi(BC.(bnd_name).type, 'Dirichlet')
        if isnumeric(BC.(bnd_name).h); BC.(bnd_name).h = @(x,y) BC.(bnd_name).h; end
        if isnumeric(BC.(bnd_name).r); BC.(bnd_name).r = @(x,y) BC.(bnd_name).r; end
    else % Neumann
        if isnumeric(BC.(bnd_name).q); BC.(bnd_name).q = @(x,y) BC.(bnd_name).q; end
        if isnumeric(BC.(bnd_name).g); BC.(bnd_name).g = @(x,y) BC.(bnd_name).g; end
    end
end

u0 = solve_elliptic_fd_general_bc(hx_macro, hy_macro, C11_macro, C22_macro, F_macro, BC);

%% 5. 步骤四: 求解二阶单胞问题
[d_cN1_dy1, d_cN1_dy2] = gradient_fd(c_micro .* N1, hx_micro, hy_micro);
[d_cN2_dy1, d_cN2_dy2] = gradient_fd(c_micro .* N2, hx_micro, hy_micro);

F_N11 = -(cH11 - c_micro - c_micro.*dN1_dy1 - d_cN1_dy1);
F_N12 = -(cH12           - c_micro.*dN2_dy1 - d_cN1_dy2);
F_N21 = -(cH21           - c_micro.*dN1_dy2 - d_cN2_dy1);
F_N22 = -(cH22 - c_micro - c_micro.*dN2_dy2 - d_cN2_dy2);

N11 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, F_N11);
N12 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, F_N12);
N21 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, F_N21);
N22 = solve_elliptic_fd(hx_micro, hy_micro, c_micro, c_micro, F_N22);

%% 6. 步骤五: 重构二阶多尺度解
x1_fine_vec = linspace(x_left, x_right, N_fine + 1);
x2_fine_vec = linspace(y_left, y_right, N_fine + 1);
[X1_fine, X2_fine] = meshgrid(x1_fine_vec, x2_fine_vec);

% 计算 u0 在精细网格上的导数
[du0_dx1, du0_dx2] = gradient_fd(u0, hx_macro, hy_macro);
[d2u0_dx11, d2u0_dx12] = gradient_fd(du0_dx1, hx_macro, hy_macro);
[~, d2u0_dx22] = gradient_fd(du0_dx2, hx_macro, hy_macro);

% 将宏观解及其导数插值到精细网格
U0_fine = interp2_manual(X1_macro, X2_macro, u0, X1_fine, X2_fine, 'spline');
dU0_dx1_fine = interp2_manual(X1_macro, X2_macro, du0_dx1, X1_fine, X2_fine, 'spline');
dU0_dx2_fine = interp2_manual(X1_macro, X2_macro, du0_dx2, X1_fine, X2_fine, 'spline');
d2U0_dx11_fine = interp2_manual(X1_macro, X2_macro, d2u0_dx11, X1_fine, X2_fine, 'spline');
d2U0_dx12_fine = interp2_manual(X1_macro, X2_macro, d2u0_dx12, X1_fine, X2_fine, 'spline');
d2U0_dx22_fine = interp2_manual(X1_macro, X2_macro, d2u0_dx22, X1_fine, X2_fine, 'spline');

% 计算精细网格上的细观坐标
Y1_fine = mod(X1_fine / epsilon, 1);
Y2_fine = mod(X2_fine / epsilon, 1);

% 将单胞解插值到精细网格
N1_fine = interp2_manual(Y1_micro, Y2_micro, N1, Y1_fine, Y2_fine, 'spline');
N2_fine = interp2_manual(Y1_micro, Y2_micro, N2, Y1_fine, Y2_fine, 'spline');
N11_fine = interp2_manual(Y1_micro, Y2_micro, N11, Y1_fine, Y2_fine, 'spline');
N12_fine = interp2_manual(Y1_micro, Y2_micro, N12, Y1_fine, Y2_fine, 'spline');
N21_fine = interp2_manual(Y1_micro, Y2_micro, N21, Y1_fine, Y2_fine, 'spline');
N22_fine = interp2_manual(Y1_micro, Y2_micro, N22, Y1_fine, Y2_fine, 'spline');

% 根据二阶展开式重构解
u1_term = N1_fine .* dU0_dx1_fine + N2_fine .* dU0_dx2_fine;
u2_term = N11_fine .* d2U0_dx11_fine + (N12_fine + N21_fine) .* d2U0_dx12_fine + N22_fine .* d2U0_dx22_fine;
u2_eps = U0_fine + epsilon * u1_term + epsilon^2 * u2_term;

%% 7. 绘图
figure;
subplot(1,2,1);
surf(X1_fine, X2_fine, u2_eps);
colorbar;
xlabel('x');
ylabel('y');
zlabel('u');
title('二阶多尺度解（三维）');
view(3);
axis tight;
shading interp;

subplot(1,2,2);
contourf(X1_fine, X2_fine, u2_eps, 20);
colorbar;
axis tight equal;
xlabel('x');
ylabel('y');
title('二阶多尺度解（二维）');
end

% ------------------------- 函数定义 -------------------------
function u = solve_elliptic_fd_general_bc(hx, hy, C11, C22, F, BC)
% 使用五点中心差分法求解二维椭圆方程
% -div(C*grad(u)) = F
    [Ny_plus_1, Nx_plus_1] = size(F);
    Nx = Nx_plus_1 - 1;
    Ny = Ny_plus_1 - 1;
    
    N_total = (Nx + 1) * (Ny + 1); 
    L = spalloc(N_total, N_total, 5 * N_total);
    b = zeros(N_total, 1);
    x_vec = BC.x_vec;
    y_vec = BC.y_vec;
    for j = 0:Ny
        for i = 0:Nx
            k = j * (Nx + 1) + i + 1; 
            
            x = x_vec(i+1);
            y = y_vec(j+1);
            % 内部节点
            if i > 0 && i < Nx && j > 0 && j < Ny
                c_east  = 0.5 * (C11(j+1, i+1) + C11(j+1, i+2));
                c_west  = 0.5 * (C11(j+1, i+1) + C11(j+1, i));
                c_north = 0.5 * (C22(j+1, i+1) + C22(j+2, i+1));
                c_south = 0.5 * (C22(j+1, i+1) + C22(j, i+1));
                L(k, k)             = (c_east + c_west) / hx^2 + (c_north + c_south) / hy^2;
                L(k, k + 1)         = -c_east / hx^2;
                L(k, k - 1)         = -c_west / hx^2;
                L(k, k + (Nx + 1))  = -c_north / hy^2;
                L(k, k - (Nx + 1))  = -c_south / hy^2;
                b(k)                = F(j + 1, i + 1);
            else
                % --- 边界和角点节点 ---
                if i == 0 && j == 0 % 左下角
                    if strcmpi(BC.left.type, 'Dirichlet')
                        L(k,k) = BC.left.h(x,y); b(k) = BC.left.r(x,y);
                    elseif strcmpi(BC.bottom.type, 'Dirichlet')
                        L(k,k) = BC.bottom.h(x,y); b(k) = BC.bottom.r(x,y);
                    else 
                        L(k,k) = 1; b(k) = 0; 
                    end
                elseif i == Nx && j == 0 % 右下角
                    if strcmpi(BC.right.type, 'Dirichlet')
                        L(k,k) = BC.right.h(x,y); b(k) = BC.right.r(x,y);
                    elseif strcmpi(BC.bottom.type, 'Dirichlet')
                        L(k,k) = BC.bottom.h(x,y); b(k) = BC.bottom.r(x,y);
                    else
                        L(k,k) = 1; b(k) = 0;
                    end
                elseif i == 0 && j == Ny % 左上角
                    if strcmpi(BC.left.type, 'Dirichlet')
                        L(k,k) = BC.left.h(x,y); b(k) = BC.left.r(x,y);
                    elseif strcmpi(BC.top.type, 'Dirichlet')
                        L(k,k) = BC.top.h(x,y); b(k) = BC.top.r(x,y);
                    else
                        L(k,k) = 1; b(k) = 0;
                    end
                elseif i == Nx && j == Ny % 右上角
                    if strcmpi(BC.right.type, 'Dirichlet')
                        L(k,k) = BC.right.h(x,y); b(k) = BC.right.r(x,y);
                    elseif strcmpi(BC.top.type, 'Dirichlet')
                        L(k,k) = BC.top.h(x,y); b(k) = BC.top.r(x,y);
                    else
                        L(k,k) = 1; b(k) = 0;
                    end
                elseif i == 0 % 左边界 
                    if strcmpi(BC.left.type, 'Dirichlet')
                        L(k, k) = BC.left.h(x, y);
                        b(k) = BC.left.r(x, y);
                    else % Neumann
                        c_node = C11(j+1, i+1);
                        L(k, k) = BC.left.q(x, y) + c_node / hx; 
                        L(k, k+1) = -c_node / hx;                
                        b(k) = BC.left.g(x, y);
                    end
                elseif i == Nx % 右边界
                    if strcmpi(BC.right.type, 'Dirichlet')
                        L(k, k) = BC.right.h(x, y);
                        b(k) = BC.right.r(x, y);
                    else % Neumann
                        c_node = C11(j+1, i+1);
                        L(k, k) = BC.right.q(x, y) + c_node / hx;
                        L(k, k-1) = -c_node / hx;
                        b(k) = BC.right.g(x, y);
                    end
                elseif j == 0 % 下边界
                    if strcmpi(BC.bottom.type, 'Dirichlet')
                        L(k, k) = BC.bottom.h(x, y);
                        b(k) = BC.bottom.r(x, y);
                    else % Neumann
                        c_node = C22(j+1, i+1);
                        L(k, k) = BC.bottom.q(x, y) + c_node / hy; 
                        L(k, k + (Nx+1)) = -c_node / hy;           
                        b(k) = BC.bottom.g(x, y);
                    end
                elseif j == Ny % 上边界
                     if strcmpi(BC.top.type, 'Dirichlet')
                        L(k, k) = BC.top.h(x, y);
                        b(k) = BC.top.r(x, y);
                    else % Neumann
                        c_node = C22(j+1, i+1);
                        L(k, k) = BC.top.q(x, y) + c_node / hy;
                        L(k, k - (Nx+1)) = -c_node / hy;
                        b(k) = BC.top.g(x, y);
                    end
                end
            end
        end
    end
    % 求解线性系统
    u_vec = L \ b;
    u = reshape(u_vec, Nx + 1, Ny + 1)';
end

function u = solve_elliptic_fd(hx, hy, C11, C22, F)
% 使用五点中心差分法求解二维椭圆方程
% -div(C*grad(u)) = F
    [Ny_plus_1, Nx_plus_1] = size(F);
    Nx = Nx_plus_1 - 1;
    Ny = Ny_plus_1 - 1;

    N_internal = (Nx - 1) * (Ny - 1); 
    L = spalloc(N_internal, N_internal, 5 * N_internal);
    b = zeros(N_internal, 1);
    
    for j = 1:(Ny - 1)
        for i = 1:(Nx - 1)
            k = (j - 1) * (Nx - 1) + i;
            c_east  = 0.5 * (C11(j+1, i+1) + C11(j+1, i+2));
            c_west  = 0.5 * (C11(j+1, i+1) + C11(j+1, i));
            c_north = 0.5 * (C22(j+1, i+1) + C22(j+2, i+1));
            c_south = 0.5 * (C22(j+1, i+1) + C22(j, i+1));
            
            L(k, k) = (c_east + c_west) / hx^2 + (c_north + c_south) / hy^2;
            
            if i > 1
                L(k, k - 1) = -c_west / hx^2;
            end
            if i < Nx - 1
                L(k, k + 1) = -c_east / hx^2;
            end
            if j > 1
                L(k, k - (Nx - 1)) = -c_south / hy^2;
            end
            if j < Ny - 1
                L(k, k + (Nx - 1)) = -c_north / hy^2;
            end
            
            b(k) = F(j + 1, i + 1);
        end
    end

    u_vec = L \ b;
    u_internal = reshape(u_vec, Nx - 1, Ny - 1)';
    
    u = zeros(Ny + 1, Nx + 1);
    u(2:Ny, 2:Nx) = u_internal;
end

function [FX, FY] = gradient_fd(F, hx, hy)
% 使用中心差分计算梯度
    [m, n] = size(F);
    FX = zeros(m, n);
    FY = zeros(m, n);

    FX(:, 2:n-1) = (F(:, 3:n) - F(:, 1:n-2)) / (2 * hx);
    FX(:, 1) = (F(:, 2) - F(:, 1)) / hx;
    FX(:, n) = (F(:, n) - F(:, n-1)) / hx;

    FY(2:m-1, :) = (F(3:m, :) - F(1:m-2, :)) / (2 * hy);
    FY(1, :) = (F(2, :) - F(1, :)) / hy;
    FY(m, :) = (F(m, :) - F(m-1, :)) / hy;
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