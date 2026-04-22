function solve_parabolic_forword(xyt_range,Nxyt, equ_para,u0,bc_type,bc_para)
    % 2D Parabolic PDE Solver: d*u_t - div(c*grad(u)) + a*u = f
    % Using forward difference method in time
    
    % Parameters that can be modified
    x_min = xyt_range(1); x_max = xyt_range(2);  % Domain in x
    y_min = xyt_range(3); y_max = xyt_range(4);  % Domain in y
    t_min = xyt_range(5); t_max = xyt_range(6); % Time domain
    
    Nx = Nxyt(1); % Number of points in x
    Ny = Nxyt(2); % Number of points in y
    Nt = Nxyt(3); % Number of time steps
    
    dx = (x_max - x_min) / (Nx - 1);
    dy = (y_max - y_min) / (Ny - 1);
    dt = (t_max - t_min) / Nt;
    
    x = linspace(x_min, x_max, Nx);
    y = linspace(y_min, y_max, Ny);
    t = linspace(t_min, t_max, Nt+1);
    
    [X, Y] = meshgrid(x, y);

    d_coeff = equ_para.d;  % Time coefficient
    c_coeff = equ_para.c;  % Diffusion coefficient 
    a_coeff = equ_para.a;  % Reaction coefficient
    f_fun = equ_para.f; % Source term
    
    % Initial condition u(x,y,0)
    initial_condition = u0;
    
    % Boundary condition
    % Dirichlet: h*u = r
    % Neumann: n*c*grad(u) + q*u = g
    
    % Left boundary (x = x_min)
    if bc_type(1)==1
        left_bc_type = 'dirichlet';
    else
        left_bc_type='neumann';
    end
    left_h = bc_para.h_left;   
    left_r = bc_para.r_left;   
    left_q = bc_para.q_left;   
    left_g = bc_para.g_left;   
    
    % Right boundary (x = x_max)
    if bc_type(2)==1
        right_bc_type = 'dirichlet';
    else
        right_bc_type = 'neumann';
    end
    
    right_h = bc_para.h_right;  
    right_r = bc_para.r_right;  
    right_q = bc_para.q_right;  
    right_g = bc_para.g_right;  
    
    % Bottom boundary (y = y_min)
    if bc_type(3)==1
        bottom_bc_type = 'dirichlet';
    else
        bottom_bc_type = 'neumann';
    end
    
    bottom_h = bc_para.h_bottom; 
    bottom_r = bc_para.r_bottom; 
    bottom_q = bc_para.q_bottom; 
    bottom_g = bc_para.g_bottom; 
    
    % Top boundary (y = y_max)
    if bc_type(4)==1
        top_bc_type = 'dirichlet';
    else
        top_bc_type = 'neumann';
    end
    
    top_h = bc_para.h_top;    
    top_r = bc_para.r_top;    
    top_q = bc_para.q_top;    
    top_g = bc_para.g_top;    
    
    u_numerical = solve_pde(x, y, t, dx, dy, dt, d_coeff, c_coeff, a_coeff, f_fun, ...
        initial_condition, ...
        left_bc_type, left_h, left_r, left_q, left_g, ...
        right_bc_type, right_h, right_r, right_q, right_g, ...
        bottom_bc_type, bottom_h, bottom_r, bottom_q, bottom_g, ...
        top_bc_type, top_h, top_r, top_q, top_g);
    
    % Plot results
    figure;
    
    % 3D surface plot of numerical solution
    subplot(1, 2, 1);
    surf(X, Y, u_numerical);
    title('数值解（三维图）');
    view(3);
    shading interp;
    axis tight;
    xlabel('x'); ylabel('y'); zlabel('u');
    
    % 2D contour map of numerical solution
    subplot(1, 2, 2);
    contourf(X, Y, u_numerical, 20); % 20 contour levels
    colorbar;
    axis tight equal;
    title('数值解（二维图）');
    xlabel('x'); ylabel('y');
    
end

function u_numerical = solve_pde(x, y, t, dx, dy, dt, d_coeff, c_coeff, a_coeff, f_fun, ...
    initial_condition, ...
    left_bc_type, left_h, left_r, left_q, left_g, ...
    right_bc_type, right_h, right_r, right_q, right_g, ...
    bottom_bc_type, bottom_h, bottom_r, bottom_q, bottom_g, ...
    top_bc_type, top_h, top_r, top_q, top_g)
    
    Nx = length(x);
    Ny = length(y);
    Nt = length(t) - 1;
    
    % Initialize solution matrices
    u = zeros(Ny, Nx, Nt+1);
    
    % Initialize with the initial condition
    for i = 1:Nx
        for j = 1:Ny
            u(j, i, 1) = initial_condition(x(i), y(j));
        end
    end
    
    % Time stepping with forward difference
    for n = 1:Nt
        current_time = t(n);
        next_time = t(n+1);
        
        % Apply boundary conditions at current time
        u = apply_boundary_conditions(u, n, x, y, current_time, dx, dy, c_coeff, ...
            left_bc_type, left_h, left_r, left_q, left_g, ...
            right_bc_type, right_h, right_r, right_q, right_g, ...
            bottom_bc_type, bottom_h, bottom_r, bottom_q, bottom_g, ...
            top_bc_type, top_h, top_r, top_q, top_g);
        
        N = Nx * Ny;
        A = sparse(N, N);
        b = zeros(N, 1);
        
        % Construct the linear system
        for i = 2:(Nx-1)
            for j = 2:(Ny-1)
                idx = (j-1)*Nx + i;

                d_ij = d_coeff(x(i), y(j), current_time);
                c_ij = c_coeff(x(i), y(j), current_time);
                a_ij = a_coeff(x(i), y(j), current_time);
                f_ij = f_fun(x(i), y(j), current_time);

                A(idx, idx) = d_ij/dt + 2*c_ij*(1/dx^2 + 1/dy^2) + a_ij;

                A(idx, idx-1) = -c_ij/dx^2;
                A(idx, idx+1) = -c_ij/dx^2;

                A(idx, idx-Nx) = -c_ij/dy^2;
                A(idx, idx+Nx) = -c_ij/dy^2;

                b(idx) = d_ij*u(j,i,n)/dt + f_ij;
            end
        end

        for j = 1:Ny
            for i = [1, Nx]  
                idx = (j-1)*Nx + i;
                
                if i == 1  
                    if strcmp(left_bc_type, 'dirichlet')
                        h_val = left_h(x(i), y(j), current_time);
                        r_val = left_r(x(i), y(j), current_time);
                        
                        A(idx, :) = 0;
                        A(idx, idx) = h_val;
                        b(idx) = r_val;
                    elseif strcmp(left_bc_type, 'neumann')
                        c_val = c_coeff(x(i), y(j), current_time);
                        q_val = left_q(x(i), y(j), current_time);
                        g_val = left_g(x(i), y(j), current_time);
                        
                        A(idx, idx) = -3*c_val/(2*dx) + q_val;
                        A(idx, idx+1) = 2*c_val/dx;
                        A(idx, idx+2) = -c_val/(2*dx);
                        b(idx) = g_val;
                    end
                else  
                    if strcmp(right_bc_type, 'dirichlet')
                        h_val = right_h(x(i), y(j), current_time);
                        r_val = right_r(x(i), y(j), current_time);
                        
                        A(idx, :) = 0;
                        A(idx, idx) = h_val;
                        b(idx) = r_val;
                    elseif strcmp(right_bc_type, 'neumann')
                        c_val = c_coeff(x(i), y(j), current_time);
                        q_val = right_q(x(i), y(j), current_time);
                        g_val = right_g(x(i), y(j), current_time);
                        
                        A(idx, idx) = 3*c_val/(2*dx) + q_val;
                        A(idx, idx-1) = -2*c_val/dx;
                        A(idx, idx-2) = c_val/(2*dx);
                        b(idx) = g_val;
                    end
                end
            end
        end
        
        for i = 2:(Nx-1)
            for j = [1, Ny] 
                idx = (j-1)*Nx + i;
                
                if j == 1
                    if strcmp(bottom_bc_type, 'dirichlet')
                        % Dirichlet: h*u = r
                        h_val = bottom_h(x(i), y(j), current_time);
                        r_val = bottom_r(x(i), y(j), current_time);
                        
                        A(idx, :) = 0;
                        A(idx, idx) = h_val;
                        b(idx) = r_val;
                    elseif strcmp(bottom_bc_type, 'neumann')
                        c_val = c_coeff(x(i), y(j), current_time);
                        q_val = bottom_q(x(i), y(j), current_time);
                        g_val = bottom_g(x(i), y(j), current_time);
                        
                        A(idx, idx) = -3*c_val/(2*dy) + q_val;
                        A(idx, idx+Nx) = 2*c_val/dy;
                        A(idx, idx+2*Nx) = -c_val/(2*dy);
                        b(idx) = g_val;
                    end
                else  
                    if strcmp(top_bc_type, 'dirichlet')
                        % Dirichlet: h*u = r
                        h_val = top_h(x(i), y(j), current_time);
                        r_val = top_r(x(i), y(j), current_time);
                        
                        A(idx, :) = 0;
                        A(idx, idx) = h_val;
                        b(idx) = r_val;
                    elseif strcmp(top_bc_type, 'neumann')
                        c_val = c_coeff(x(i), y(j), current_time);
                        q_val = top_q(x(i), y(j), current_time);
                        g_val = top_g(x(i), y(j), current_time);
                        
                        A(idx, idx) = 3*c_val/(2*dy) + q_val;
                        A(idx, idx-Nx) = -2*c_val/dy;
                        A(idx, idx-2*Nx) = c_val/(2*dy);
                        b(idx) = g_val;
                    end
                end
            end
        end
        
        % Corner points
        corner_indices = [1, Nx, (Ny-1)*Nx+1, Ny*Nx];
        corner_names = {'bottom-left', 'bottom-right', 'top-left', 'top-right'};
        
        for k = 1:length(corner_indices)
            idx = corner_indices(k);

            if k == 1  % Bottom-left corner
                if strcmp(left_bc_type, 'dirichlet') || strcmp(bottom_bc_type, 'dirichlet')
                    A(idx, :) = 0;
                    A(idx, idx) = 1;
                    b(idx) = 0; 
                else
                    A(idx, idx) = 1;
                    A(idx, idx+1) = -0.5;
                    A(idx, idx+Nx) = -0.5;
                    b(idx) = 0;
                end
            elseif k == 2  
                if strcmp(right_bc_type, 'dirichlet') || strcmp(bottom_bc_type, 'dirichlet')
                    A(idx, :) = 0;
                    A(idx, idx) = 1;
                    b(idx) = 0;
                else
                    A(idx, idx) = 1;
                    A(idx, idx-1) = -0.5;
                    A(idx, idx+Nx) = -0.5;
                    b(idx) = 0;
                end
            elseif k == 3 
                if strcmp(left_bc_type, 'dirichlet') || strcmp(top_bc_type, 'dirichlet')
                    A(idx, :) = 0;
                    A(idx, idx) = 1;
                    b(idx) = 0;
                else
                    A(idx, idx) = 1;
                    A(idx, idx+1) = -0.5;
                    A(idx, idx-Nx) = -0.5;
                    b(idx) = 0;
                end
            else 
                if strcmp(right_bc_type, 'dirichlet') || strcmp(top_bc_type, 'dirichlet')
                    A(idx, :) = 0;
                    A(idx, idx) = 1;
                    b(idx) = 0;
                else
                    A(idx, idx) = 1;
                    A(idx, idx-1) = -0.5;
                    A(idx, idx-Nx) = -0.5;
                    b(idx) = 0;
                end
            end
        end
        
        % Solve the system using an iterative method
        % Check if A is symmetric positive definite for PCG
        u_next=A\b;
        u(:,:,n+1) = reshape(u_next, Nx, Ny)';
    end
    
    % Final numerical solution
    u_numerical = u(:,:,end);
end

function u = apply_boundary_conditions(u, time_idx, x, y, t, dx, dy, c_coeff, ...
    left_bc_type, left_h, left_r, left_q, left_g, ...
    right_bc_type, right_h, right_r, right_q, right_g, ...
    bottom_bc_type, bottom_h, bottom_r, bottom_q, bottom_g, ...
    top_bc_type, top_h, top_r, top_q, top_g)
    
    Nx = length(x);
    Ny = length(y);
    
    % Apply boundary conditions

    for j = 1:Ny
        if strcmp(left_bc_type, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = left_h(x(1), y(j), t);
            r_val = left_r(x(1), y(j), t);
            u(j, 1, time_idx) = r_val / h_val;
        elseif strcmp(left_bc_type, 'neumann')
            c_val = c_coeff(x(1), y(j), t);
            q_val = left_q(x(1), y(j), t);
            g_val = left_g(x(1), y(j), t);

            if abs(q_val) < 1e-10
                u(j, 1, time_idx) = (4*u(j, 2, time_idx) - u(j, 3, time_idx) - 2*dx*g_val/c_val) / 3;
            else
                u(j, 1, time_idx) = (4*u(j, 2, time_idx) - u(j, 3, time_idx) - 2*dx*g_val/c_val) / (3 + 2*dx*q_val/c_val);
            end
        end
    end
    
    % Right boundary (x = x_max)
    for j = 1:Ny
        if strcmp(right_bc_type, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = right_h(x(Nx), y(j), t);
            r_val = right_r(x(Nx), y(j), t);
            u(j, Nx, time_idx) = r_val / h_val;
        elseif strcmp(right_bc_type, 'neumann')
            % Neumann: n*c*grad(u) + q*u = g
            c_val = c_coeff(x(Nx), y(j), t);
            q_val = right_q(x(Nx), y(j), t);
            g_val = right_g(x(Nx), y(j), t);
            
            if abs(q_val) < 1e-10
                u(j, Nx, time_idx) = (4*u(j, Nx-1, time_idx) - u(j, Nx-2, time_idx) + 2*dx*g_val/c_val) / 3;
            else
                u(j, Nx, time_idx) = (4*u(j, Nx-1, time_idx) - u(j, Nx-2, time_idx) + 2*dx*g_val/c_val) / (3 + 2*dx*q_val/c_val);
            end
        end
    end
    
    % Bottom boundary (y = y_min)
    for i = 1:Nx
        if strcmp(bottom_bc_type, 'dirichlet')
            h_val = bottom_h(x(i), y(1), t);
            r_val = bottom_r(x(i), y(1), t);
            u(1, i, time_idx) = r_val / h_val;
        elseif strcmp(bottom_bc_type, 'neumann')
            c_val = c_coeff(x(i), y(1), t);
            q_val = bottom_q(x(i), y(1), t);
            g_val = bottom_g(x(i), y(1), t);
            
            if abs(q_val) < 1e-10
                u(1, i, time_idx) = (4*u(2, i, time_idx) - u(3, i, time_idx) - 2*dy*g_val/c_val) / 3;
            else
                u(1, i, time_idx) = (4*u(2, i, time_idx) - u(3, i, time_idx) - 2*dy*g_val/c_val) / (3 + 2*dy*q_val/c_val);
            end
        end
    end
    
    for i = 1:Nx
        if strcmp(top_bc_type, 'dirichlet')
            % Dirichlet: h*u = r
            h_val = top_h(x(i), y(Ny), t);
            r_val = top_r(x(i), y(Ny), t);
            u(Ny, i, time_idx) = r_val / h_val;
        elseif strcmp(top_bc_type, 'neumann')
            % Neumann: n*c*grad(u) + q*u = g
            c_val = c_coeff(x(i), y(Ny), t);
            q_val = top_q(x(i), y(Ny), t);
            g_val = top_g(x(i), y(Ny), t);
            
            if abs(q_val) < 1e-10
                u(Ny, i, time_idx) = (4*u(Ny-1, i, time_idx) - u(Ny-2, i, time_idx) + 2*dy*g_val/c_val) / 3;
            else
                u(Ny, i, time_idx) = (4*u(Ny-1, i, time_idx) - u(Ny-2, i, time_idx) + 2*dy*g_val/c_val) / (3 + 2*dy*q_val/c_val);
            end
        end
    end
end 