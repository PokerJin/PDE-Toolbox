function solve_parabolic_backword(xyt_range,Nxyt, equ_para,u0,bc_type,bc_para)
    % Solving 2D parabolic PDE: d*u' - div(c*grad(u)) + a*u = f
    % using backward difference for time and finite differences for space
    
    %% Define problem parameters
    % Domain boundaries
    x_left = xyt_range(1);
    x_right = xyt_range(2);
    y_left = xyt_range(3);
    y_right = xyt_range(4);
    t_start = xyt_range(5);
    t_end = xyt_range(6);
    
    % Grid resolution
    Nx = Nxyt(1); 
    Ny = Nxyt(2); 
    Nt = Nxyt(3); 
    
    % Calculate grid spacings
    dx = (x_right - x_left) / (Nx - 1);
    dy = (y_right - y_left) / (Ny - 1);
    dt = (t_end - t_start) / Nt;

    x = linspace(x_left, x_right, Nx);
    y = linspace(y_left, y_right, Ny);
    [X, Y] = meshgrid(x, y);
    
    %% Parameters for PDE: can be constants or function handles
    % d*u' - div(c*grad(u)) + a*u = f
    
    d = equ_para.d;     % Time coefficient 
    c = equ_para.c;     % Diffusion coefficient 
    a = equ_para.a;     % Reaction coefficient
    f = equ_para.f;     % Source term
    
    % Initial condition
    u_init = u0;
    
    %% Boundary condition settings
   

    map = {'dirichlet', 'neumann'};

    % 2. 直接使用 bc_type 作为索引来从 map 中选取元素
    bc_types = map(bc_type);

    % 显示结果
 
    % bc_types = {'dirichlet', 'dirichlet', 'dirichlet', 'dirichlet'};
    
    % Dirichlet: h*u = r
    left_h = bc_para.h_left;    
    left_r = bc_para.r_left;     
    right_h = bc_para.h_right;    
    right_r = bc_para.r_right;    
    bottom_h = bc_para.h_bottom;   
    bottom_r = bc_para.r_bottom;  
    top_h = bc_para.h_top;      
    top_r = bc_para.r_top;    
    
    % Neumann: n*c*grad(u) + q*u = g
    left_q = bc_para.q_left;   
    left_g = bc_para.g_left;     
    right_q = bc_para.q_right;   
    right_g = bc_para.g_right;    
    bottom_q = bc_para.q_bottom;   
    bottom_g = bc_para.g_bottom;   
    top_q = bc_para.q_top;      
    top_g = bc_para.g_top;     
    
    %% Initialize solution array
    u = zeros(Ny, Nx, Nt+1);
    
    % Set initial condition
    for i = 1:Ny
        for j = 1:Nx
            u(i, j, 1) = u_init(X(i, j), Y(i, j));
        end
    end
    
    % Apply boundary conditions for initial state
    u = apply_boundary_conditions(u, 1, x, y, dx, dy, 0, bc_types, ...
        left_h, left_r, right_h, right_r, bottom_h, bottom_r, top_h, top_r, ...
        left_q, left_g, right_q, right_g, bottom_q, bottom_g, top_q, top_g, c);
    
    %% Main time stepping loop
    for n = 1:Nt
        t = t_start + n * dt;

        u(:,:,n+1) = solve_time_step(u(:,:,n), Nx, Ny, dx, dy, dt, t, X, Y, ...
                                     d, c, a, f, bc_types, ...
                                     left_h, left_r, right_h, right_r, ...
                                     bottom_h, bottom_r, top_h, top_r, ...
                                     left_q, left_g, right_q, right_g, ...
                                     bottom_q, bottom_g, top_q, top_g);
    end
    
    %% Visualization
    figure;
    
    % 3D surface plot
    subplot(1,2,1);
    surf(X, Y, u(:,:,end));
    shading interp;
    title('数值解（三维图）');
    xlabel('x'); ylabel('y'); zlabel('u');
    view(3);
    colorbar;
    axis tight;
    
    % 2D contour map
    subplot(1,2,2);
    contourf(X, Y, u(:,:,end), 20);
    title('数值解（二维图）');
    xlabel('x'); ylabel('y');
    colorbar;
    axis tight equal;
end

function u_new = solve_time_step(u_old, Nx, Ny, dx, dy, dt, t, X, Y, d, c, a, f, bc_types, ...
                                 left_h, left_r, right_h, right_r, bottom_h, bottom_r, top_h, top_r, ...
                                 left_q, left_g, right_q, right_g, bottom_q, bottom_g, top_q, top_g)
    N = Nx * Ny;
    A = sparse(N, N);
    b = zeros(N, 1);
    
    % Process interior points
    for i = 2:Ny-1
        for j = 2:Nx-1
            row = (i-1)*Nx + j;
            x = X(i,j);
            y = Y(i,j);
            
            % Time derivative coefficient
            d_val = d(x, y, t);
            
            % Diffusion coefficient at various points
            c_ij = c(x, y, t);
            c_ip1j = c(x, y+dy/2, t); 
            c_im1j = c(x, y-dy/2, t);  
            c_ijp1 = c(x+dx/2, y, t);  
            c_ijm1 = c(x-dx/2, y, t); 
            
            % Reaction coefficient
            a_val = a(x, y, t);
            
            % Source term
            f_val = f(x, y, t);
 
            A(row, row) = d_val/dt + (c_ip1j + c_im1j)/dy^2 + (c_ijp1 + c_ijm1)/dx^2 + a_val;
            
            % Off-diagonal entries - coefficients of neighboring points
            A(row, row-Nx) = -c_im1j/dy^2; 
            A(row, row+Nx) = -c_ip1j/dy^2; 
            A(row, row-1) = -c_ijm1/dx^2;  
            A(row, row+1) = -c_ijp1/dx^2;  

            b(row) = f_val + d_val * u_old(i,j) / dt;
        end
    end
    

    j = 1;
    for i = 1:Ny
        row = (i-1)*Nx + j;
        y = Y(i,j);
        x = X(i,j);
        
        if strcmp(bc_types{1}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = left_h(x, y, t);
            r_val = left_r(x, y, t);
            
            A(row, :) = 0;
            A(row, row) = h_val;
            b(row) = r_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = left_q(x, y, t);
            g_val = left_g(x, y, t);
            c_val = c(x, y, t);
            
            A(row, row) = c_val/dx + q_val;
            A(row, row+1) = -c_val/dx;
            b(row) = g_val;
        end
    end
    
    % Right boundary (x=1)
    j = Nx;
    for i = 1:Ny
        row = (i-1)*Nx + j;
        y = Y(i,j);
        x = X(i,j);
        
        if strcmp(bc_types{2}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = right_h(x, y, t);
            r_val = right_r(x, y, t);
            
            A(row, :) = 0;
            A(row, row) = h_val;
            b(row) = r_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = right_q(x, y, t);
            g_val = right_g(x, y, t);
            c_val = c(x, y, t);
            
            A(row, row) = c_val/dx + q_val;
            A(row, row-1) = -c_val/dx;
            b(row) = g_val;
        end
    end
    
    % Bottom boundary (y=0)
    i = 1;
    for j = 2:Nx-1
        row = (i-1)*Nx + j;
        x = X(i,j);
        y = Y(i,j);
        
        if strcmp(bc_types{3}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = bottom_h(x, y, t);
            r_val = bottom_r(x, y, t);
            
            A(row, :) = 0;
            A(row, row) = h_val;
            b(row) = r_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = bottom_q(x, y, t);
            g_val = bottom_g(x, y, t);
            c_val = c(x, y, t);
            
            A(row, row) = c_val/dy + q_val;
            A(row, row+Nx) = -c_val/dy;
            b(row) = g_val;
        end
    end
    
    % Top boundary (y=1)
    i = Ny;
    for j = 2:Nx-1
        row = (i-1)*Nx + j;
        x = X(i,j);
        y = Y(i,j);
        
        if strcmp(bc_types{4}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = top_h(x, y, t);
            r_val = top_r(x, y, t);
            
            A(row, :) = 0;
            A(row, row) = h_val;
            b(row) = r_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = top_q(x, y, t);
            g_val = top_g(x, y, t);
            c_val = c(x, y, t);
            
            A(row, row) = c_val/dy + q_val;
            A(row, row-Nx) = -c_val/dy;
            b(row) = g_val;
        end
    end
    
    % Solve the linear system using an appropriate iterative solver
    % Check if matrix A is symmetric positive definite for PCG
    U=A\b;
    
    % Reshape solution vector to 2D grid
    u_new = reshape(U, Nx, Ny)';
end

function u = apply_boundary_conditions(u, n, x, y, dx, dy, t, bc_types, ...
                                      left_h, left_r, right_h, right_r, ...
                                      bottom_h, bottom_r, top_h, top_r, ...
                                      left_q, left_g, right_q, right_g, ...
                                      bottom_q, bottom_g, top_q, top_g, c)
    Nx = length(x);
    Ny = length(y);
    
    % Left boundary (x=0)
    j = 1;
    for i = 1:Ny
        if strcmp(bc_types{1}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = left_h(x(j), y(i), t);
            r_val = left_r(x(j), y(i), t);
            u(i, j, n) = r_val / h_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = left_q(x(j), y(i), t);
            g_val = left_g(x(j), y(i), t);
            c_val = c(x(j), y(i), t);
            
            % Forward difference for Neumann
            u(i, j, n) = (g_val + c_val*u(i, j+1, n)/dx) / (q_val + c_val/dx);
        end
    end
    
    % Right boundary (x=1)
    j = Nx;
    for i = 1:Ny
        if strcmp(bc_types{2}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = right_h(x(j), y(i), t);
            r_val = right_r(x(j), y(i), t);
            u(i, j, n) = r_val / h_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = right_q(x(j), y(i), t);
            g_val = right_g(x(j), y(i), t);
            c_val = c(x(j), y(i), t);
            
            % Backward difference for Neumann
            u(i, j, n) = (g_val + c_val*u(i, j-1, n)/dx) / (q_val + c_val/dx);
        end
    end
    
    % Bottom boundary (y=0)
    i = 1;
    for j = 2:Nx-1
        if strcmp(bc_types{3}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = bottom_h(x(j), y(i), t);
            r_val = bottom_r(x(j), y(i), t);
            u(i, j, n) = r_val / h_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = bottom_q(x(j), y(i), t);
            g_val = bottom_g(x(j), y(i), t);
            c_val = c(x(j), y(i), t);
            
            % Forward difference for Neumann
            u(i, j, n) = (g_val + c_val*u(i+1, j, n)/dy) / (q_val + c_val/dy);
        end
    end
    
    % Top boundary (y=1)
    i = Ny;
    for j = 2:Nx-1
        if strcmp(bc_types{4}, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = top_h(x(j), y(i), t);
            r_val = top_r(x(j), y(i), t);
            u(i, j, n) = r_val / h_val;
        else
            % Neumann: n*c*grad(u) + q*u = g
            q_val = top_q(x(j), y(i), t);
            g_val = top_g(x(j), y(i), t);
            c_val = c(x(j), y(i), t);
            
            % Backward difference for Neumann
            u(i, j, n) = (g_val + c_val*u(i-1, j, n)/dy) / (q_val + c_val/dy);
        end
    end
end 