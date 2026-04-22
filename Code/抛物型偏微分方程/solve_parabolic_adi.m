function solve_parabolic_adi(xyt_range,Nxyt, equ_para,u0,bc_type,bc_para)
    % =================================================
    % Problem setup - modify these parameters as needed
    % =================================================
    
    % Domain parameters
    xmin = xyt_range(1); xmax = xyt_range(2);  
    ymin = xyt_range(1); ymax = xyt_range(2); 
    tmin = xyt_range(1); tmax = xyt_range(2); 
    
    Nx = Nxyt(1);  
    Ny = Nxyt(2);  
    Nt = Nxyt(3);
    
    dx = (xmax - xmin) / (Nx - 1);
    dy = (ymax - ymin) / (Ny - 1);
    dt = (tmax - tmin) / Nt;
    
    x = linspace(xmin, xmax, Nx);
    y = linspace(ymin, ymax, Ny);
    t = linspace(tmin, tmax, Nt+1);
    
    [X, Y] = meshgrid(x, y);
    
    % PDE coefficients (can be functions of x, y, t)
    d = equ_para.d;                % Coefficient of time derivative
    c = equ_para.c;                % Diffusion coefficient
    a = equ_para.a;                % Reaction coefficient
    f = equ_para.f;                % Source term
    
    % Boundary condition parameters - separate for each boundary
    % Dirichlet: h*u = r on boundary
    h_left = bc_para.h_left;     
    r_left = bc_para.r_left;       
    
    h_right = bc_para.h_right;        
    r_right = bc_para.r_right;          
    
    h_bottom = bc_para.h_bottom;     
    r_bottom = bc_para.r_bottom;    
    
    h_top = bc_para.h_top;         
    r_top = bc_para.r_top;        
    
    % Neumann: n*c*grad(u) + q*u = g on boundary
    q_left = bc_para.q_left;         
    g_left = bc_para.g_left;          
    
    q_right = bc_para.q_right;         
    g_right = bc_para.g_right;         
    
    q_bottom = bc_para.q_bottom;       
    g_bottom = bc_para.g_bottom;    
    
    q_top = bc_para.q_top;          
    g_top = bc_para.g_top;         
    
    % Create cell arrays for boundary condition functions
    h_funcs = {h_left, h_right, h_bottom, h_top};
    r_funcs = {r_left, r_right, r_bottom, r_top};
    q_funcs = {q_left, q_right, q_bottom, q_top};
    g_funcs = {g_left, g_right, g_bottom, g_top};
    
    % Initial condition
    u_init = u0; 
    
    % Boundary type flags: 0 for Dirichlet, 1 for Neumann
    bc_type = bc_type-1; % [left, right, bottom, top]
    
    % Iterative solver parameters
    tol = 1e-8;             % Tolerance for iterative solvers
    max_iter = 1000;        % Maximum iterations (will be limited by matrix size)
    
    % =================================================
    % Initialize solution and run solver
    % =================================================
    
    % Initialize solution matrix
    u = zeros(Ny, Nx, Nt+1);
    
    % Set initial condition
    for i = 1:Ny
        for j = 1:Nx
            u(i, j, 1) = u_init(x(j), y(i));
        end
    end
    
    % Solve using ADI method with iterative solvers
    u = adi_solver(u, x, y, t, dx, dy, dt, d, c, a, f, h_funcs, r_funcs, q_funcs, g_funcs, bc_type, tol, max_iter);


    figure;
    subplot(1,2,1);
    surf(X, Y, u(:,:,end));
    title('数值解（三维图）');
    xlabel('x'); ylabel('y'); zlabel('u');
    view(3);
    shading interp;
    colorbar;
    axis tight;

    subplot(1,2,2);
    contourf(X, Y, u(:,:,end), 20);
    title('数值解（二维图）');
    xlabel('x'); ylabel('y');
    colorbar;
    axis tight equal;

    % Add overall title
    %sgtitle(sprintf('The Solution of Parabolic PDE'));

    % Here are the commands to display the plot window full screen
    % Get a handle to the current graphics window
    % hFig = gcf;
    
    % Setting the graphics window to full screen
    % set(hFig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
end

function u = adi_solver(u, x, y, t, dx, dy, dt, d, c, a, f, h_funcs, r_funcs, q_funcs, g_funcs, bc_type, tol, max_iter)
    % ADI solver for 2D parabolic PDE
    % d*u_t - div(c*grad(u)) + a*u = f
    
    Nx = length(x);
    Ny = length(y);
    Nt = length(t) - 1;
    
    % Unpack boundary condition types
    bc_left = bc_type(1);    % x = xmin
    bc_right = bc_type(2);   % x = xmax
    bc_bottom = bc_type(3);  % y = ymin
    bc_top = bc_type(4);     % y = ymax
    
    % Unpack boundary condition function handles
    h_left = h_funcs{1};
    h_right = h_funcs{2};
    h_bottom = h_funcs{3};
    h_top = h_funcs{4};
    
    r_left = r_funcs{1};
    r_right = r_funcs{2};
    r_bottom = r_funcs{3};
    r_top = r_funcs{4};
    
    q_left = q_funcs{1};
    q_right = q_funcs{2};
    q_bottom = q_funcs{3};
    q_top = q_funcs{4};
    
    g_left = g_funcs{1};
    g_right = g_funcs{2};
    g_bottom = g_funcs{3};
    g_top = g_funcs{4};
    
    % Intermediate solution
    u_half = zeros(Ny, Nx);
    
    % Main time stepping loop
    for n = 1:Nt
        t_current = t(n);
        t_next = t(n+1);
        
        % First half step: implicit in x-direction, explicit in y-direction
        for i = 1:Ny
            % Set up tridiagonal system for this row
            A = zeros(Nx, Nx);
            b = zeros(Nx, 1);
            
            for j = 2:Nx-1
                % PDE coefficients at current position
                d_val = d(x(j), y(i), t_current + dt/2);
                c_val = c(x(j), y(i), t_current + dt/2);
                a_val = a(x(j), y(i), t_current + dt/2);
                f_val = f(x(j), y(i), t_current + dt/2);
                
                % Discrete operators in x-direction (implicit)
                A(j, j-1) = -c_val * dt/(2*dx^2);
                A(j, j) = d_val + c_val * dt/dx^2 + a_val * dt/2;
                A(j, j+1) = -c_val * dt/(2*dx^2);
                
                % RHS: includes explicit y-derivative terms
                y_diff_term = 0;
                if i > 1 && i < Ny
                    y_diff_term = c_val * (u(i+1, j, n) - 2*u(i, j, n) + u(i-1, j, n))/dy^2;
                end
                
                b(j) = d_val * u(i, j, n) + dt/2 * (y_diff_term - a_val*u(i, j, n) + f_val);
            end
            
            % Apply boundary conditions for x-direction
            % Left boundary (j = 1)
            if bc_left == 0  % Dirichlet
                h_val = h_left(x(1), y(i), t_current + dt/2);
                r_val = r_left(x(1), y(i), t_current + dt/2);
                A(1, 1) = h_val;
                b(1) = r_val;
            else  % Neumann
                c_val = c(x(1), y(i), t_current + dt/2);
                q_val = q_left(x(1), y(i), t_current + dt/2);
                g_val = g_left(x(1), y(i), t_current + dt/2);
                
                A(1, 1) = -3*c_val/(2*dx) + q_val;
                A(1, 2) = 4*c_val/(2*dx);
                A(1, 3) = -c_val/(2*dx);
                b(1) = g_val;
            end
            
            % Right boundary (j = Nx)
            if bc_right == 0  % Dirichlet
                h_val = h_right(x(Nx), y(i), t_current + dt/2);
                r_val = r_right(x(Nx), y(i), t_current + dt/2);
                A(Nx, Nx) = h_val;
                b(Nx) = r_val;
            else  % Neumann
                c_val = c(x(Nx), y(i), t_current + dt/2);
                q_val = q_right(x(Nx), y(i), t_current + dt/2);
                g_val = g_right(x(Nx), y(i), t_current + dt/2);
                
                A(Nx, Nx-2) = -c_val/(2*dx);
                A(Nx, Nx-1) = 4*c_val/(2*dx);
                A(Nx, Nx) = -3*c_val/(2*dx) + q_val;
                b(Nx) = g_val;
            end
            
            % Convert to sparse for efficiency
            A_sparse = sparse(A);
            
            u_row=A\b;
            u_half(i, :) = u_row';
        end
        
        % Second half step: explicit in x-direction, implicit in y-direction
        for j = 1:Nx
            % Set up tridiagonal system for this column
            A = zeros(Ny, Ny);
            b = zeros(Ny, 1);
            
            for i = 2:Ny-1
                % PDE coefficients at current position
                d_val = d(x(j), y(i), t_next);
                c_val = c(x(j), y(i), t_next);
                a_val = a(x(j), y(i), t_next);
                f_val = f(x(j), y(i), t_next);
                
                % Discrete operators in y-direction (implicit)
                A(i, i-1) = -c_val * dt/(2*dy^2);
                A(i, i) = d_val + c_val * dt/dy^2 + a_val * dt/2;
                A(i, i+1) = -c_val * dt/(2*dy^2);
                
                % RHS: includes explicit x-derivative terms from half step
                x_diff_term = 0;
                if j > 1 && j < Nx
                    x_diff_term = c_val * (u_half(i, j+1) - 2*u_half(i, j) + u_half(i, j-1))/dx^2;
                end
                
                b(i) = d_val * u_half(i, j) + dt/2 * (x_diff_term - a_val*u_half(i, j) + f_val);
            end
            
            % Apply boundary conditions for y-direction
            % Bottom boundary (i = 1)
            if bc_bottom == 0  % Dirichlet
                h_val = h_bottom(x(j), y(1), t_next);
                r_val = r_bottom(x(j), y(1), t_next);
                A(1, 1) = h_val;
                b(1) = r_val;
            else  % Neumann
                c_val = c(x(j), y(1), t_next);
                q_val = q_bottom(x(j), y(1), t_next);
                g_val = g_bottom(x(j), y(1), t_next);
                
                A(1, 1) = -3*c_val/(2*dy) + q_val;
                A(1, 2) = 4*c_val/(2*dy);
                A(1, 3) = -c_val/(2*dy);
                b(1) = g_val;
            end
            
            % Top boundary (i = Ny)
            if bc_top == 0  % Dirichlet
                h_val = h_top(x(j), y(Ny), t_next);
                r_val = r_top(x(j), y(Ny), t_next);
                A(Ny, Ny) = h_val;
                b(Ny) = r_val;
            else  % Neumann
                c_val = c(x(j), y(Ny), t_next);
                q_val = q_top(x(j), y(Ny), t_next);
                g_val = g_top(x(j), y(Ny), t_next);
                
                A(Ny, Ny-2) = -c_val/(2*dy);
                A(Ny, Ny-1) = 4*c_val/(2*dy);
                A(Ny, Ny) = -3*c_val/(2*dy) + q_val;
                b(Ny) = g_val;
            end
            
            % Convert to sparse for efficiency
            A_sparse = sparse(A);
            
            u_col=A_sparse\b;
            u(:, j, n+1) = u_col;
        end
    end
end