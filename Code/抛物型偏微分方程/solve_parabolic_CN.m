function solve_parabolic_CN(xyt_range,Nxyt, equ_para,u0,bc_type,bc_para)
%% Parameters
% Domain parameters
x_min = xyt_range(1); x_max = xyt_range(2);
y_min = xyt_range(3); y_max = xyt_range(4);
t_min = xyt_range(5); t_max = xyt_range(6);

% Mesh parameters
Nx = Nxyt(1); 
Ny = Nxyt(2); 
Nt = Nxyt(3); 

dx = (x_max - x_min) / (Nx + 1);
dy = (y_max - y_min) / (Ny + 1);
dt = (t_max - t_min) / Nt;

% Create meshgrid
x = linspace(x_min, x_max, Nx + 2);
y = linspace(y_min, y_max, Ny + 2);
t = linspace(t_min, t_max, Nt + 1);
[X, Y] = meshgrid(x, y);

% d*u_t - div(c*grad(u)) + a*u = f
d = equ_para.d;                       % Coefficient of time derivative
c = equ_para.c;                       % Diffusion coefficient
a = equ_para.a;                       % Reaction coefficient
f = equ_para.f;                       % Source term

% Initial condition: u(x,y,0) = u0(x,y)
% u0 = @(x, y) sin(pi*x) .* sin(pi*y);

% Boundary conditions
% 1: Dirichlet (h*u = r)
% 2: Neumann (n*c*grad(u) + q*u = g)
left_boundary_type = bc_type(1);   
right_boundary_type = bc_type(2);   
bottom_boundary_type = bc_type(3);  
top_boundary_type = bc_type(4);    

% Dirichlet boundary parameters (h*u = r)
left_h = bc_para.h_left;
left_r = bc_para.r_left;
right_h = bc_para.h_right;
right_r = bc_para.r_right;
bottom_h = bc_para.h_bottom;
bottom_r = bc_para.r_bottom;
top_h = bc_para.h_top;
top_r = bc_para.r_top;

% Neumann boundary parameters (n*c*grad(u) + q*u = g)
left_q = bc_para.q_left;
left_g = bc_para.r_left;
right_q = bc_para.q_right;
right_g = bc_para.g_right;
bottom_q = bc_para.q_bottom;
bottom_g = bc_para.g_bottom;
top_q = bc_para.q_top;
top_g = bc_para.g_top;

%% Solve the PDE
[U_numerical, x_grid, y_grid] = solve_2d_parabolic_pde(d, c, a, f, u0, ...
    x_min, x_max, y_min, y_max, t_min, t_max, Nx, Ny, Nt, ...
    left_boundary_type, right_boundary_type, bottom_boundary_type, top_boundary_type, ...
    left_h, left_r, right_h, right_r, bottom_h, bottom_r, top_h, top_r, ...
    left_q, left_g, right_q, right_g, bottom_q, bottom_g, top_q, top_g);

%% Visualize the results
figure;

% Plot 3D surface plot
subplot(1, 2, 1);
surf(X, Y, U_numerical);
title('数值解（三维图）');
xlabel('x'); ylabel('y'); zlabel('u');
colorbar;
shading interp;
axis tight;
view(3);

% Plot contour map
subplot(1, 2, 2);
contourf(X, Y, U_numerical, 20);
title('数值解（二维图）');
xlabel('x'); ylabel('y');
colorbar;
axis equal tight; 

end
function [U, x, y] = solve_2d_parabolic_pde(d, c, a, f, u0, ...
    x_min, x_max, y_min, y_max, t_min, t_max, Nx, Ny, Nt, ...
    left_boundary_type, right_boundary_type, bottom_boundary_type, top_boundary_type, ...
    left_h, left_r, right_h, right_r, bottom_h, bottom_r, top_h, top_r, ...
    left_q, left_g, right_q, right_g, bottom_q, bottom_g, top_q, top_g)

% Initialize the grid
dx = (x_max - x_min) / (Nx + 1);
dy = (y_max - y_min) / (Ny + 1);
dt = (t_max - t_min) / Nt;

x = linspace(x_min, x_max, Nx + 2);
y = linspace(y_min, y_max, Ny + 2);
t = linspace(t_min, t_max, Nt + 1);

% Initialize solution arrays
U_prev = zeros(Ny + 2, Nx + 2);
U = zeros(Ny + 2, Nx + 2);

% Set initial condition
for i = 1:Ny+2
    for j = 1:Nx+2
        U_prev(i, j) = u0(x(j), y(i));
    end
end

theta = 0.5;

% Time stepping
for n = 1:Nt
    % Current time
    t_current = t(n);
    t_next = t(n+1);

    N = (Nx+2) * (Ny+2);

    A = sparse(N, N);
    b = zeros(N, 1);

    for i = 2:Ny+1
        for j = 2:Nx+1
            % Index in the linear system
            idx = (i-1)*(Nx+2) + j;

            xij = x(j);
            yij = y(i);

            if isa(d, 'function_handle')
                d_val = d(xij, yij, t_current);
            else
                d_val = d;
            end
            
            if isa(c, 'function_handle')
                c_val = c(xij, yij, t_current);
            else
                c_val = c;
            end
            
            if isa(a, 'function_handle')
                a_val = a(xij, yij, t_current);
            else
                a_val = a;
            end
            
            if isa(f, 'function_handle')
                f_curr = f(xij, yij, t_current);
                f_next = f(xij, yij, t_next);
            else
                f_curr = f;
                f_next = f;
            end

            diag_term = d_val/dt + theta*(2*c_val/dx^2 + 2*c_val/dy^2 + a_val);
            A(idx, idx) = diag_term;

            A(idx, (i-2)*(Nx+2) + j) = -theta * c_val/dy^2;
            
            A(idx, i*(Nx+2) + j) = -theta * c_val/dy^2;
            
            A(idx, (i-1)*(Nx+2) + j-1) = -theta * c_val/dx^2;
            
            A(idx, (i-1)*(Nx+2) + j+1) = -theta * c_val/dx^2;
            
            b(idx) = U_prev(i, j) * d_val/dt + (1-theta) * ( ...
                c_val/dx^2 * (U_prev(i, j-1) - 2*U_prev(i, j) + U_prev(i, j+1)) + ...
                c_val/dy^2 * (U_prev(i-1, j) - 2*U_prev(i, j) + U_prev(i+1, j)) - ...
                a_val * U_prev(i, j) ...
            ) + theta * f_next + (1-theta) * f_curr;
        end
    end
    
    % Apply boundary conditions
    % Left boundary 
    j = 1;
    for i = 1:Ny+2
        idx = (i-1)*(Nx+2) + j;
        yij = y(i);
        xij = x(j);
        
        if left_boundary_type == 1  % Dirichlet
            if isa(left_h, 'function_handle')
                h = left_h(xij, yij, t_next);
            else
                h = left_h;
            end
            
            if isa(left_r, 'function_handle')
                r = left_r(xij, yij, t_next);
            else
                r = left_r;
            end
            
            % h*u = r -> u = r/h
            A(idx, :) = 0;
            A(idx, idx) = h;
            b(idx) = r;
            
        elseif left_boundary_type == 2  % Neumann
            if isa(left_q, 'function_handle')
                q = left_q(xij, yij, t_next);
            else
                q = left_q;
            end
            
            if isa(left_g, 'function_handle')
                g = left_g(xij, yij, t_next);
            else
                g = left_g;
            end
            
            if isa(c, 'function_handle')
                c_val = c(xij, yij, t_next);
            else
                c_val = c;
            end
            
            A(idx, idx) = -1/dx + q/c_val;
            A(idx, idx+1) = 1/dx;
            b(idx) = g/c_val;
        end
    end
    
    % Right boundary
    j = Nx+2;
    for i = 1:Ny+2
        idx = (i-1)*(Nx+2) + j;
        yij = y(i);
        xij = x(j);
        
        if right_boundary_type == 1  % Dirichlet
            if isa(right_h, 'function_handle')
                h = right_h(xij, yij, t_next);
            else
                h = right_h;
            end
            
            if isa(right_r, 'function_handle')
                r = right_r(xij, yij, t_next);
            else
                r = right_r;
            end
            
            % h*u = r -> u = r/h
            A(idx, :) = 0;
            A(idx, idx) = h;
            b(idx) = r;
            
        elseif right_boundary_type == 2  % Neumann
            if isa(right_q, 'function_handle')
                q = right_q(xij, yij, t_next);
            else
                q = right_q;
            end
            
            if isa(right_g, 'function_handle')
                g = right_g(xij, yij, t_next);
            else
                g = right_g;
            end
            
            if isa(c, 'function_handle')
                c_val = c(xij, yij, t_next);
            else
                c_val = c;
            end

            A(idx, idx) = 1/dx + q/c_val;
            A(idx, idx-1) = -1/dx;
            b(idx) = g/c_val;
        end
    end
    
    % Bottom boundary (y = y_min)
    i = 1;
    for j = 1:Nx+2
        idx = (i-1)*(Nx+2) + j;
        xij = x(j);
        yij = y(i);
        
        if bottom_boundary_type == 1  % Dirichlet
            if isa(bottom_h, 'function_handle')
                h = bottom_h(xij, yij, t_next);
            else
                h = bottom_h;
            end
            
            if isa(bottom_r, 'function_handle')
                r = bottom_r(xij, yij, t_next);
            else
                r = bottom_r;
            end

            A(idx, :) = 0;
            A(idx, idx) = h;
            b(idx) = r;
            
        elseif bottom_boundary_type == 2  % Neumann
            if isa(bottom_q, 'function_handle')
                q = bottom_q(xij, yij, t_next);
            else
                q = bottom_q;
            end
            
            if isa(bottom_g, 'function_handle')
                g = bottom_g(xij, yij, t_next);
            else
                g = bottom_g;
            end
            
            if isa(c, 'function_handle')
                c_val = c(xij, yij, t_next);
            else
                c_val = c;
            end

            A(idx, idx) = -1/dy + q/c_val;
            A(idx, idx+(Nx+2)) = 1/dy;
            b(idx) = g/c_val;
        end
    end
    
    % Top boundary
    i = Ny+2;
    for j = 1:Nx+2
        idx = (i-1)*(Nx+2) + j;
        xij = x(j);
        yij = y(i);
        
        if top_boundary_type == 1  % Dirichlet
            if isa(top_h, 'function_handle')
                h = top_h(xij, yij, t_next);
            else
                h = top_h;
            end
            
            if isa(top_r, 'function_handle')
                r = top_r(xij, yij, t_next);
            else
                r = top_r;
            end
            
            A(idx, :) = 0;
            A(idx, idx) = h;
            b(idx) = r;
            
        elseif top_boundary_type == 2  % Neumann
            if isa(top_q, 'function_handle')
                q = top_q(xij, yij, t_next);
            else
                q = top_q;
            end
            
            if isa(top_g, 'function_handle')
                g = top_g(xij, yij, t_next);
            else
                g = top_g;
            end
            
            if isa(c, 'function_handle')
                c_val = c(xij, yij, t_next);
            else
                c_val = c;
            end
            
            A(idx, idx) = 1/dy + q/c_val;
            A(idx, idx-(Nx+2)) = -1/dy;
            b(idx) = g/c_val;
        end
    end
    
    % Solve the linear system using an iterative method
    
    u_flat=A\b;

    U = reshape(u_flat, Ny+2, Nx+2);
    
    U_prev = U;
end

end 