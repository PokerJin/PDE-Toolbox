function solve_hyperbolic_leapfrog(xyt_range,Nxyt,equ_para,init,bc_type,bc_para)
% Solves the 2D hyperbolic PDE:
%   d*u_tt - div(c*grad(u)) + a*u = f
% using the Leapfrog finite difference method.

%% --- User-Defined Parameters ---

% --- Domain ---
x_left = xyt_range(1);
x_right = xyt_range(2);
y_left = xyt_range(3);
y_right = xyt_range(4);
t_start = xyt_range(5);
t_end = xyt_range(6);

% --- Discretization ---
Nx = Nxyt(1); 
Ny = Nxyt(2); 
Nt = Nxyt(3); 

% --- PDE Coefficients (Functions of x, y, t or constants) ---
d_coeff = equ_para.d;
c_coeff = equ_para.c;
a_coeff = equ_para.a;
f_source = equ_para.f;

% --- Initial Conditions ---
u0_func = init.u0;
v0_func = init.v0; 

% --- Boundary Conditions ---

% Left Boundary
if bc_type(1)==1
    left_type = 'dirichlet';
else
    left_type = 'neumann';
end


left_h_func = bc_para.h_left;
left_r_func = bc_para.r_left;
left_q_func = bc_para.q_left; 
left_g_func = bc_para.g_left; 

% Right Boundary
if bc_type(2)==1
    right_type = 'dirichlet';
else
    right_type = 'neumann';
end

right_h_func = bc_para.h_right;
right_r_func = bc_para.r_right;
right_q_func = bc_para.q_right; 
right_g_func = bc_para.g_right; 

% Bottom Boundary
if bc_type(3)==1
    bottom_type = 'dirichlet';
else
    bottom_type = 'neumann';
end

bottom_h_func = bc_para.h_bottom;
bottom_r_func = bc_para.r_bottom;
bottom_q_func = bc_para.q_bottom; 
bottom_g_func = bc_para.g_bottom;

% Top Boundary
if bc_type(4)==1
    top_type = 'dirichlet';
else
    top_type = 'neumann';
end

top_h_func = bc_para.h_top;
top_r_func = bc_para.r_top;
top_q_func = bc_para.q_top; 
top_g_func = bc_para.g_top; 


%% --- Setup Grid and Parameters ---

dx = (x_right - x_left) / (Nx - 1);
dy = (y_right - y_left) / (Ny - 1);
dt = (t_end - t_start) / (Nt - 1);

x = linspace(x_left, x_right, Nx);
y = linspace(y_left, y_right, Ny);
t = linspace(t_start, t_end, Nt);

[X, Y] = meshgrid(x, y); % Create spatial grid

c_center_val = eval_coeff(c_coeff, (x_left+x_right)/2, (y_left+y_right)/2, t_start);
d_center_val = eval_coeff(d_coeff, (x_left+x_right)/2, (y_left+y_right)/2, t_start);

if d_center_val <= 0
    error('Coefficient d must be positive for stability check.');
end
if c_center_val <= 0
    warning('Coefficient c is non-positive, CFL condition cannot be estimated reliably.');
    cfl_limit = Inf;
else
    speed_est = sqrt(c_center_val / d_center_val); % Approx wave speed
    cfl_limit = 1 / (speed_est * sqrt(1/dx^2 + 1/dy^2));
end

if dt > cfl_limit
    warning('dt = %.4f may exceed the CFL stability limit (%.4f). Solution might be unstable.', dt, cfl_limit);
end

% --- Initialize Solution Arrays ---
U_prev = zeros(Ny, Nx); 
U_curr = zeros(Ny, Nx); 
U_next = zeros(Ny, Nx); 

% --- Evaluate Initial Conditions ---
U_curr = u0_func(X, Y);
if isscalar(U_curr)
    U_curr = U_curr * ones(size(X));
end
V0 = v0_func(X, Y);
if isscalar(V0) 
    V0 = V0 * ones(size(X));
end


D0 = eval_coeff(d_coeff, X, Y, t_start);
C0 = eval_coeff(c_coeff, X, Y, t_start);
A0 = eval_coeff(a_coeff, X, Y, t_start);
F0 = eval_coeff(f_source, X, Y, t_start);

if any(D0(:) <= 0)
    error('Coefficient d must be strictly positive everywhere.');
end

%% --- First Time Step (Using Taylor Expansion) ---

Div_C_Grad_U0 = zeros(Ny, Nx);
for j = 2:Ny-1
    for i = 2:Nx-1
        c_iph = 0.5 * (C0(j, i+1) + C0(j, i));
        c_imh = 0.5 * (C0(j, i) + C0(j, i-1));
        c_jph = 0.5 * (C0(j+1, i) + C0(j, i));
        c_jmh = 0.5 * (C0(j, i) + C0(j-1, i));

        Fx_iph = c_iph * (U_curr(j, i+1) - U_curr(j, i)) / dx;
        Fx_imh = c_imh * (U_curr(j, i) - U_curr(j, i-1)) / dx;
        Fy_jph = c_jph * (U_curr(j+1, i) - U_curr(j, i)) / dy;
        Fy_jmh = c_jmh * (U_curr(j, i) - U_curr(j-1, i)) / dy;

        Div_C_Grad_U0(j, i) = (Fx_iph - Fx_imh) / dx + (Fy_jph - Fy_jmh) / dy;
    end
end

% Calculate u_tt at t=0
U_tt0 = (1 ./ D0) .* (Div_C_Grad_U0 - A0 .* U_curr + F0);

% Calculate U at the first time step (k=1)
U_next = U_curr + dt * V0 + 0.5 * dt^2 * U_tt0;

% --- Apply Boundary Conditions at First Time Step (t = dt) ---
t_next = t_start + dt;
y_vec = y(:); % Column vector for y coordinates
x_vec = x(:)'; % Row vector for x coordinates

eval_bc_param = @(func, x_arg, y_arg, t_arg) eval_coeff(func, x_arg, y_arg, t_arg);

% Left Boundary (i=1, x=x_left)
x_bnd = x_left;
x_grid_bnd = x_bnd * ones(size(y_vec)); 
if strcmp(left_type, 'dirichlet') 
    h_vals = eval_bc_param(left_h_func, x_grid_bnd, y_vec, t_next);
    r_vals = eval_bc_param(left_r_func, x_grid_bnd, y_vec, t_next);
    non_zero_h = abs(h_vals) > eps; 
    U_next(non_zero_h, 1) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
elseif strcmp(left_type, 'neumann')
    q_vals = eval_bc_param(left_q_func, x_grid_bnd, y_vec, t_next);
    g_vals = eval_bc_param(left_g_func, x_grid_bnd, y_vec, t_next);
    C_vals = eval_coeff(c_coeff, x_grid_bnd, y_vec, t_next); 
    denom = q_vals + C_vals / dx;
    non_zero_denom = abs(denom) > eps;
    U_next(non_zero_denom, 1) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(non_zero_denom, 2) / dx) ./ denom(non_zero_denom);
end

x_bnd = x_right;
x_grid_bnd = x_bnd * ones(size(y_vec));
if strcmp(right_type, 'dirichlet') 
    h_vals = eval_bc_param(right_h_func, x_grid_bnd, y_vec, t_next);
    r_vals = eval_bc_param(right_r_func, x_grid_bnd, y_vec, t_next);
    non_zero_h = abs(h_vals) > eps;
    U_next(non_zero_h, Nx) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
elseif strcmp(right_type, 'neumann')
    q_vals = eval_bc_param(right_q_func, x_grid_bnd, y_vec, t_next);
    g_vals = eval_bc_param(right_g_func, x_grid_bnd, y_vec, t_next);
    C_vals = eval_coeff(c_coeff, x_grid_bnd, y_vec, t_next);
    denom = q_vals + C_vals / dx;
    non_zero_denom = abs(denom) > eps;
    U_next(non_zero_denom, Nx) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(non_zero_denom, Nx-1) / dx) ./ denom(non_zero_denom);
end

% Bottom Boundary
y_bnd = y_left;
y_grid_bnd = y_bnd * ones(size(x_vec)); 
if strcmp(bottom_type, 'dirichlet') 
    h_vals = eval_bc_param(bottom_h_func, x_vec, y_grid_bnd, t_next);
    r_vals = eval_bc_param(bottom_r_func, x_vec, y_grid_bnd, t_next);
    non_zero_h = abs(h_vals) > eps;
    U_next(1, non_zero_h) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
elseif strcmp(bottom_type, 'neumann')
    q_vals = eval_bc_param(bottom_q_func, x_vec, y_grid_bnd, t_next);
    g_vals = eval_bc_param(bottom_g_func, x_vec, y_grid_bnd, t_next);
    C_vals = eval_coeff(c_coeff, x_vec, y_grid_bnd, t_next); 
    denom = q_vals + C_vals / dy;
    non_zero_denom = abs(denom) > eps;
    U_next(1, non_zero_denom) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(2, non_zero_denom) / dy) ./ denom(non_zero_denom);
end

% Top Boundary 
y_bnd = y_right;
y_grid_bnd = y_bnd * ones(size(x_vec));
if strcmp(top_type, 'dirichlet') 
    h_vals = eval_bc_param(top_h_func, x_vec, y_grid_bnd, t_next);
    r_vals = eval_bc_param(top_r_func, x_vec, y_grid_bnd, t_next);
    non_zero_h = abs(h_vals) > eps;
    U_next(Ny, non_zero_h) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
elseif strcmp(top_type, 'neumann')
    q_vals = eval_bc_param(top_q_func, x_vec, y_grid_bnd, t_next);
    g_vals = eval_bc_param(top_g_func, x_vec, y_grid_bnd, t_next);
    C_vals = eval_coeff(c_coeff, x_vec, y_grid_bnd, t_next); 
    denom = q_vals + C_vals / dy;
    non_zero_denom = abs(denom) > eps;
    U_next(Ny, non_zero_denom) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(Ny-1, non_zero_denom) / dy) ./ denom(non_zero_denom);
end

% --- Handle Corners---

if strcmp(left_type, 'dirichlet') && strcmp(bottom_type, 'dirichlet')
    x_c = x_left; y_c = y_left;
    h_l = eval_bc_param(left_h_func, x_c, y_c, t_next); r_l = eval_bc_param(left_r_func, x_c, y_c, t_next);
    h_b = eval_bc_param(bottom_h_func, x_c, y_c, t_next); r_b = eval_bc_param(bottom_r_func, x_c, y_c, t_next);
    if abs(h_b) > eps
        U_next(1, 1) = r_b / h_b;
    elseif abs(h_l) > eps
         U_next(1, 1) = r_l / h_l;
    end
end
if strcmp(right_type, 'dirichlet') && strcmp(bottom_type, 'dirichlet')
     x_c = x_right; y_c = y_left;
     h_r = eval_bc_param(right_h_func, x_c, y_c, t_next); r_r = eval_bc_param(right_r_func, x_c, y_c, t_next);
     h_b = eval_bc_param(bottom_h_func, x_c, y_c, t_next); r_b = eval_bc_param(bottom_r_func, x_c, y_c, t_next);
     if abs(h_b) > eps
        U_next(1, Nx) = r_b / h_b;
     elseif abs(h_r) > eps
         U_next(1, Nx) = r_r / h_r;
     end
end
if strcmp(left_type, 'dirichlet') && strcmp(top_type, 'dirichlet')
     x_c = x_left; y_c = y_right;
     h_l = eval_bc_param(left_h_func, x_c, y_c, t_next); r_l = eval_bc_param(left_r_func, x_c, y_c, t_next);
     h_t = eval_bc_param(top_h_func, x_c, y_c, t_next); r_t = eval_bc_param(top_r_func, x_c, y_c, t_next);
     if abs(h_t) > eps
        U_next(Ny, 1) = r_t / h_t;
     elseif abs(h_l) > eps
         U_next(Ny, 1) = r_l / h_l;
     end
end
if strcmp(right_type, 'dirichlet') && strcmp(top_type, 'dirichlet')
     x_c = x_right; y_c = y_right;
     h_r = eval_bc_param(right_h_func, x_c, y_c, t_next); r_r = eval_bc_param(right_r_func, x_c, y_c, t_next);
     h_t = eval_bc_param(top_h_func, x_c, y_c, t_next); r_t = eval_bc_param(top_r_func, x_c, y_c, t_next);
     if abs(h_t) > eps
        U_next(Ny, Nx) = r_t / h_t;
     elseif abs(h_r) > eps
         U_next(Ny, Nx) = r_r / h_r;
     end
end


% Update history
U_prev = U_curr;
U_curr = U_next;

%% --- Time Stepping Loop (Leapfrog) ---

for k = 2:Nt-1 
    t_curr = t(k); 
    t_next = t(k+1); 


    D = eval_coeff(d_coeff, X, Y, t_curr);
    C = eval_coeff(c_coeff, X, Y, t_curr);
    A = eval_coeff(a_coeff, X, Y, t_curr);
    F = eval_coeff(f_source, X, Y, t_curr);

    if any(D(:) <= 0)
        error('Coefficient d must be strictly positive everywhere at t = %.4f.', t_curr);
    end

    Div_C_Grad_U = zeros(Ny, Nx);
    for j = 2:Ny-1
        for i = 2:Nx-1
            c_iph = 0.5 * (C(j, i+1) + C(j, i));
            c_imh = 0.5 * (C(j, i) + C(j, i-1));
            c_jph = 0.5 * (C(j+1, i) + C(j, i));
            c_jmh = 0.5 * (C(j, i) + C(j-1, i));

            Fx_iph = c_iph * (U_curr(j, i+1) - U_curr(j, i)) / dx;
            Fx_imh = c_imh * (U_curr(j, i) - U_curr(j, i-1)) / dx;
            Fy_jph = c_jph * (U_curr(j+1, i) - U_curr(j, i)) / dy;
            Fy_jmh = c_jmh * (U_curr(j, i) - U_curr(j-1, i)) / dy;

            Div_C_Grad_U(j, i) = (Fx_iph - Fx_imh) / dx + (Fy_jph - Fy_jmh) / dy;
        end
    end

    U_next(2:Ny-1, 2:Nx-1) = 2 * U_curr(2:Ny-1, 2:Nx-1) - U_prev(2:Ny-1, 2:Nx-1) ...
        + (dt^2 ./ D(2:Ny-1, 2:Nx-1)) .* (Div_C_Grad_U(2:Ny-1, 2:Nx-1) ...
        - A(2:Ny-1, 2:Nx-1) .* U_curr(2:Ny-1, 2:Nx-1) + F(2:Ny-1, 2:Nx-1));


    % Left Boundary (i=1, x=x_left)
    x_bnd = x_left;
    x_grid_bnd = x_bnd * ones(size(y_vec));
    if strcmp(left_type, 'dirichlet') 
        h_vals = eval_bc_param(left_h_func, x_grid_bnd, y_vec, t_next);
        r_vals = eval_bc_param(left_r_func, x_grid_bnd, y_vec, t_next);
        non_zero_h = abs(h_vals) > eps;
        U_next(non_zero_h, 1) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
    elseif strcmp(left_type, 'neumann')
        q_vals = eval_bc_param(left_q_func, x_grid_bnd, y_vec, t_next);
        g_vals = eval_bc_param(left_g_func, x_grid_bnd, y_vec, t_next);
        C_vals = eval_coeff(c_coeff, x_grid_bnd, y_vec, t_next);
        denom = q_vals + C_vals / dx;
        non_zero_denom = abs(denom) > eps;
        U_next(non_zero_denom, 1) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(non_zero_denom, 2) / dx) ./ denom(non_zero_denom);
    end

    % Right Boundary 
    x_bnd = x_right;
    x_grid_bnd = x_bnd * ones(size(y_vec));
    if strcmp(right_type, 'dirichlet') 
        h_vals = eval_bc_param(right_h_func, x_grid_bnd, y_vec, t_next);
        r_vals = eval_bc_param(right_r_func, x_grid_bnd, y_vec, t_next);
        non_zero_h = abs(h_vals) > eps;
        U_next(non_zero_h, Nx) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
    elseif strcmp(right_type, 'neumann') 
        q_vals = eval_bc_param(right_q_func, x_grid_bnd, y_vec, t_next);
        g_vals = eval_bc_param(right_g_func, x_grid_bnd, y_vec, t_next);
        C_vals = eval_coeff(c_coeff, x_grid_bnd, y_vec, t_next); 
        denom = q_vals + C_vals / dx;
        non_zero_denom = abs(denom) > eps;
        U_next(non_zero_denom, Nx) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(non_zero_denom, Nx-1) / dx) ./ denom(non_zero_denom);
    end

    % Bottom Boundary 
    y_bnd = y_left;
    y_grid_bnd = y_bnd * ones(size(x_vec));
    if strcmp(bottom_type, 'dirichlet') 
        h_vals = eval_bc_param(bottom_h_func, x_vec, y_grid_bnd, t_next);
        r_vals = eval_bc_param(bottom_r_func, x_vec, y_grid_bnd, t_next);
        non_zero_h = abs(h_vals) > eps;
        U_next(1, non_zero_h) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
    elseif strcmp(bottom_type, 'neumann') 
        q_vals = eval_bc_param(bottom_q_func, x_vec, y_grid_bnd, t_next);
        g_vals = eval_bc_param(bottom_g_func, x_vec, y_grid_bnd, t_next);
        C_vals = eval_coeff(c_coeff, x_vec, y_grid_bnd, t_next); 
        denom = q_vals + C_vals / dy;
        non_zero_denom = abs(denom) > eps;
        U_next(1, non_zero_denom) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(2, non_zero_denom) / dy) ./ denom(non_zero_denom);
    end

    % Top Boundary 
    y_bnd = y_right;
    y_grid_bnd = y_bnd * ones(size(x_vec));
    if strcmp(top_type, 'dirichlet') 
        h_vals = eval_bc_param(top_h_func, x_vec, y_grid_bnd, t_next);
        r_vals = eval_bc_param(top_r_func, x_vec, y_grid_bnd, t_next);
        non_zero_h = abs(h_vals) > eps;
        U_next(Ny, non_zero_h) = r_vals(non_zero_h) ./ h_vals(non_zero_h);
    elseif strcmp(top_type, 'neumann')
        q_vals = eval_bc_param(top_q_func, x_vec, y_grid_bnd, t_next);
        g_vals = eval_bc_param(top_g_func, x_vec, y_grid_bnd, t_next);
        C_vals = eval_coeff(c_coeff, x_vec, y_grid_bnd, t_next); 
        denom = q_vals + C_vals / dy;
        non_zero_denom = abs(denom) > eps;
        U_next(Ny, non_zero_denom) = (g_vals(non_zero_denom) + C_vals(non_zero_denom) .* U_next(Ny-1, non_zero_denom) / dy) ./ denom(non_zero_denom);
    end

    % --- Handle Corners ---
    if strcmp(left_type, 'dirichlet') && strcmp(bottom_type, 'dirichlet')
        x_c = x_left; y_c = y_left;
        h_l = eval_bc_param(left_h_func, x_c, y_c, t_next); r_l = eval_bc_param(left_r_func, x_c, y_c, t_next);
        h_b = eval_bc_param(bottom_h_func, x_c, y_c, t_next); r_b = eval_bc_param(bottom_r_func, x_c, y_c, t_next);
        if abs(h_b) > eps
            U_next(1, 1) = r_b / h_b;
        elseif abs(h_l) > eps
             U_next(1, 1) = r_l / h_l;
        end
    end
    if strcmp(right_type, 'dirichlet') && strcmp(bottom_type, 'dirichlet')
         x_c = x_right; y_c = y_left;
         h_r = eval_bc_param(right_h_func, x_c, y_c, t_next); r_r = eval_bc_param(right_r_func, x_c, y_c, t_next);
         h_b = eval_bc_param(bottom_h_func, x_c, y_c, t_next); r_b = eval_bc_param(bottom_r_func, x_c, y_c, t_next);
         if abs(h_b) > eps
            U_next(1, Nx) = r_b / h_b;
         elseif abs(h_r) > eps
             U_next(1, Nx) = r_r / h_r;
         end
    end
    if strcmp(left_type, 'dirichlet') && strcmp(top_type, 'dirichlet')
         x_c = x_left; y_c = y_right;
         h_l = eval_bc_param(left_h_func, x_c, y_c, t_next); r_l = eval_bc_param(left_r_func, x_c, y_c, t_next);
         h_t = eval_bc_param(top_h_func, x_c, y_c, t_next); r_t = eval_bc_param(top_r_func, x_c, y_c, t_next);
         if abs(h_t) > eps
            U_next(Ny, 1) = r_t / h_t;
         elseif abs(h_l) > eps
             U_next(Ny, 1) = r_l / h_l;
         end
    end
    if strcmp(right_type, 'dirichlet') && strcmp(top_type, 'dirichlet')
         x_c = x_right; y_c = y_right;
         h_r = eval_bc_param(right_h_func, x_c, y_c, t_next); r_r = eval_bc_param(right_r_func, x_c, y_c, t_next);
         h_t = eval_bc_param(top_h_func, x_c, y_c, t_next); r_t = eval_bc_param(top_r_func, x_c, y_c, t_next);
         if abs(h_t) > eps
            U_next(Ny, Nx) = r_t / h_t;
         elseif abs(h_r) > eps
             U_next(Ny, Nx) = r_r / h_r;
         end
    end


    U_prev = U_curr;
    U_curr = U_next;


end 

U_final = U_curr;

%% --- Plotting Final Results ---

figure;

% Plot 3D Surface
subplot(1, 2, 1);
surf(X, Y, U_final);
title('数值解（三维图）');
xlabel('x');
ylabel('y');
zlabel('u');
colorbar;
axis tight;
view(3); 
shading interp;
grid on;

% Plot 2D Contour Map
subplot(1, 2, 2);
contourf(X, Y, U_final, 20); 
title('数值解（二维图）');
xlabel('x');
ylabel('y');
colorbar;
axis tight equal; 
grid on;
end

%% --- Helper function to ensure matrix output from coefficients ---
function coeff_mat = eval_coeff(coeff_func, x_grid, y_grid, t_val)
    val = coeff_func(x_grid, y_grid, t_val);
    if isscalar(val)
        coeff_mat = val * ones(size(x_grid));
    elseif isequal(size(val), size(x_grid))
        coeff_mat = val;
    else
        error('Coefficient function must return a scalar or a matrix matching the grid size.');
    end
end