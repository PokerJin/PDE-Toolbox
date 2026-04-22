function solve_mulhyp_dir(xyt_range,Nxyt,equ_para,init,bc_type,bc_para)
    % Solves the 2D multiscale hyperbolic PDE:
    % d * u_tt - div(c * grad(u)) = f
    % using the Leapfrog finite difference method.

    % Domain [xL, xR] x [yL, yR] (must be square)
    xL = xyt_range(1);
    xR = xyt_range(2);
    yL = xyt_range(3);
    yR = xyt_range(4); % xR - xL must equal yR - yL

    % Time interval [t0, tF]
    t0 = xyt_range(5);
    tF = xyt_range(6);

    % Discretization steps
    dx = (xR-xL)/Nxyt(1);
    dy = (yR-yL)/Nxyt(2); 
    dt = (tF-t0)/Nxyt(3);

    % Multiscale parameters
    epsilon = equ_para.epsilon*(xR-xL);   % Microcell size (characteristic length of heterogeneity)
    r_circle = equ_para.len; % Radius of the circle in the center of each microcell (must be < epsilon/2)
    c_in = equ_para.c_in;      % Value of c inside the circles
    c_out = equ_para.c_ma;       % Value of c outside the circles
    d_in = equ_para.d_in;        % Value of d inside the circles
    d_out = equ_para.d_ma;     % Value of d outside the circles

    % Source term f(x, y, t)
    f_func = equ_para.f;

    % Initial conditions u(x, y, t0) and u_t(x, y, t0)
    u0_func = init.u0;
    ut0_func = init.v0;
    if bc_type(1)==1
        bc_para.g_left=@(x,y,t)0;
        bc_para.q_left=@(x,y,t)0;
    else
        bc_para.h_left=@(x,y,t)0;
        bc_para.r_left=@(x,y,t)0;
    end
    if bc_type(2)==1
        bc_para.g_right=@(x,y,t)0;
        bc_para.q_right=@(x,y,t)0;
    else
        bc_para.h_right=@(x,y,t)0;
        bc_para.r_right=@(x,y,t)0;
    end
    if bc_type(3)==1
        bc_para.g_bottom=@(x,y,t)0;
        bc_para.q_bottom=@(x,y,t)0;
    else
        bc_para.h_bottom=@(x,y,t)0;
        bc_para.r_bottom=@(x,y,t)0;
    end
    if bc_type(4)==1
        bc_para.g_top=@(x,y,t)0;
        bc_para.q_top=@(x,y,t)0;
    else
        bc_para.h_top=@(x,y,t)0;
        bc_para.r_top=@(x,y,t)0;
    end

    % Left Boundary (x = xL)
    if bc_type(1)==1
        BC.Left.type = 'Dirichlet';
    else
        BC.Left.type = 'Neumann';
    end
    
    BC.Left.h = bc_para.h_left; 
    BC.Left.r = bc_para.r_left; 
    BC.Left.q = bc_para.q_left;
    BC.Left.g = bc_para.g_left; 

    % Right Boundary (x = xR)
    if bc_type(2)==1
        BC.Right.type = 'Dirichlet';
    else
        BC.Right.type = 'Neumann';
    end
    
    BC.Right.h = bc_para.h_right; 
    BC.Right.r = bc_para.r_right; 
    BC.Right.q = bc_para.q_right;
    BC.Right.g = bc_para.g_right; 

    % Bottom Boundary (y = yL)
    if bc_type(3)==1
        BC.Bottom.type = 'Dirichlet';
    else
        BC.Bottom.type = 'Neumann';
    end
    BC.Bottom.h = bc_para.h_bottom; 
    BC.Bottom.r = bc_para.r_bottom; 
    BC.Bottom.q = bc_para.q_bottom;
    BC.Bottom.g = bc_para.g_bottom; 

    % Top Boundary (y = yR)
    if bc_type(4)==1
        BC.Top.type = 'Dirichlet';
    else
        BC.Top.type = 'Neumann';
    end
    
    BC.Top.h = bc_para.h_top; 
    BC.Top.r = bc_para.r_top; 
    BC.Top.q = bc_para.q_top;
    BC.Top.g = bc_para.g_top; 

    % --- End of User-Defined Parameters ---

    % --- Setup Grid ---
    if abs((xR - xL) - (yR - yL)) > 1e-10
        error('Domain must be square (xR - xL = yR - yL).');
    end
    x = (xL:dx:xR)';
    y = (yL:dy:yR)'; 
    t = t0:dt:tF;    

    Nx = length(x);
    Ny = length(y);
    Nt = length(t);

    [X, Y] = meshgrid(x, y);
    X = X'; 
    Y = Y';

    % --- Calculate Coefficients c and d ---
    C = zeros(Nx, Ny);
    D = zeros(Nx, Ny);

    for i = 1:Nx
        for j = 1:Ny
            xi = x(i);
            yj = y(j);

            cell_ix = floor((xi - xL + 1e-9*epsilon) / epsilon);
            cell_iy = floor((yj - yL + 1e-9*epsilon) / epsilon);

            xc = xL + epsilon * (cell_ix + 0.5);
            yc = yL + epsilon * (cell_iy + 0.5);

            dist_sq = (xi - xc)^2 + (yj - yc)^2;

            if dist_sq <= r_circle^2
                C(i, j) = c_in;
                D(i, j) = d_in;
            else
                C(i, j) = c_out;
                D(i, j) = d_out;
            end
        end
    end

    % --- Stability Check (Simplified CFL) ---
    max_wave_speed = sqrt(max(C(:) ./ D(:)));
    cfl_limit_x = max_wave_speed * dt / dx;
    cfl_limit_y = max_wave_speed * dt / dy;
    cfl_combined = sqrt(cfl_limit_x^2 + cfl_limit_y^2);

    if cfl_combined > 1.0
        warning('CFL condition (sqrt((c/d)_max) * dt * sqrt(1/dx^2 + 1/dy^2) <= 1) might be violated (%.4f > 1). Solution may be unstable.', cfl_combined);
    end
    if r_circle >= epsilon/2
         warning('Circle radius r_circle (%.4f) should be less than half the cell size epsilon/2 (%.4f) to avoid overlap.', r_circle, epsilon/2);
    end


    % --- Initialize Solution Arrays ---
    U_prev = zeros(Nx, Ny); 
    U_curr = zeros(Nx, Ny); 
    U_next = zeros(Nx, Ny); 

    % Set initial condition u(x, y, t0)
    U_prev = u0_func(X, Y);

    % Evaluate initial velocity and source term at t0
    Ut0 = ut0_func(X, Y);
    F0 = f_func(X, Y, t(1)); % t(1) = t0

    % Calculate spatial operator term L(U_prev) = div(C grad U_prev) at t=0
    Lu0 = compute_spatial_operator(U_prev, C, dx, dy, Nx, Ny, BC, t(1), x, y, xL, xR, yL, yR);

    % Calculate U_curr (solution at t=t1)
    U_curr = U_prev + dt * Ut0 + (dt^2 ./ (2 * D)) .* (Lu0 + F0);

    % Apply boundary conditions at t=t1 (n=1)
    U_curr = apply_boundary_conditions(U_curr, BC, t(2), x, y, Nx, Ny, xL, xR, yL, yR); % t(2) corresponds to n=1

    for n = 2:Nt-1 
        tn = t(n); 

        Fn = f_func(X, Y, tn);

        Lu_curr = compute_spatial_operator(U_curr, C, dx, dy, Nx, Ny, BC, tn, x, y, xL, xR, yL, yR);

        U_next = 2*U_curr - U_prev + (dt^2 ./ D) .* (Lu_curr + Fn);

        U_next = apply_boundary_conditions(U_next, BC, t(n+1), x, y, Nx, Ny, xL, xR, yL, yR);

        U_prev = U_curr;
        U_curr = U_next;
    end


    % --- Plotting Results ---
    figure;

    % 1. 3D Surface plot of the final solution U(tF)
    subplot(1, 2, 1);
    surf(y, x, U_curr); % Use y, x for correct orientation with meshgrid/matrix
    shading interp;
    colorbar;
    title('参考解（三维图）');
    xlabel('y');
    ylabel('x');
    zlabel('u');
    axis tight;
    view(3); % 3D view

    % 2. 2D Contour plot of the final solution U(tF)
    subplot(1, 2, 2);
    contourf(y, x, U_curr, 20); % Use y, x for correct orientation
    colorbar;
    title('参考解（二维图）');
    xlabel('y');
    ylabel('x');
    axis equal tight;

    % % 3. 2D plot of coefficient c(x, y)
    % subplot(1, 4, 3);
    % imagesc(y, x, C); % Use y, x for correct orientation
    % colorbar;
    % title('Coefficient c(x,y)');
    % xlabel('y');
    % ylabel('x');
    % axis equal tight;
    % set(gca,'YDir','normal'); % Ensure y-axis goes from bottom to top
    % 
    % % 4. 2D plot of coefficient d(x, y)
    % subplot(1, 4, 4);
    % imagesc(y, x, D); % Use y, x for correct orientation
    % colorbar;
    % title('Coefficient d(x,y)');
    % xlabel('y');
    % ylabel('x');
    % axis equal tight;
    % set(gca,'YDir','normal'); % Ensure y-axis goes from bottom to top

end

% =========================================================================
% Helper function to compute the spatial operator div(c * grad(u))
% =========================================================================
function Lu = compute_spatial_operator(U, C, dx, dy, Nx, Ny, BC, t_curr, x, y, xL, xR, yL, yR)
    Lu = zeros(Nx, Ny);

    % --- Compute C at half-grid points (arithmetic mean) ---
    Cx_half = zeros(Nx+1, Ny); % Staggered grid for x-fluxes
    Cy_half = zeros(Nx, Ny+1); % Staggered grid for y-fluxes

    % Interior points
    Cx_half(2:Nx, :) = (C(1:Nx-1, :) + C(2:Nx, :)) / 2;
    Cy_half(:, 2:Ny) = (C(:, 1:Ny-1) + C(:, 2:Ny)) / 2;

    % Boundary points (use value at boundary node)
    Cx_half(1, :) = C(1, :);
    Cx_half(Nx+1, :) = C(Nx, :);
    Cy_half(:, 1) = C(:, 1);
    Cy_half(:, Ny+1) = C(:, Ny);

    % --- Compute derivatives using ghost points implicitly for Neumann ---
    % Create padded U array to handle boundary conditions easily
    U_padded = zeros(Nx+2, Ny+2);
    U_padded(2:Nx+1, 2:Ny+1) = U;

    % --- Apply BCs to find ghost point values ---
    x_bnd_L = repmat(xL, size(y)); % x coordinates on Left boundary
    x_bnd_R = repmat(xR, size(y)); % x coordinates on Right boundary
    y_bnd_B = repmat(yL, size(x)); % y coordinates on Bottom boundary
    y_bnd_T = repmat(yR, size(x)); % y coordinates on Top boundary

    % Left Boundary (i=1)
    if strcmp(BC.Left.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Left.h, x_bnd_L, y, t_curr)'; 
        r_val = eval_bc_func(BC.Left.r, x_bnd_L, y, t_curr)'; 
        h_val(abs(h_val) < 1e-15) = 1; 
        U_padded(2, 2:Ny+1) = r_val ./ h_val; 
        U_padded(1, 2:Ny+1) = 2*U_padded(2, 2:Ny+1) - U_padded(3, 2:Ny+1);
    elseif strcmp(BC.Left.type, 'Neumann')
        q_val = eval_bc_func(BC.Left.q, x_bnd_L, y, t_curr)'; 
        g_val = eval_bc_func(BC.Left.g, x_bnd_L, y, t_curr)'; 
        U1j = U(1, :); 
        C1j = C(1, :); 
        U2j = U(2, :); 
        C1j_safe = C1j; C1j_safe(abs(C1j_safe)<1e-15) = 1e-15;
        U_padded(1, 2:Ny+1) = U2j - (2*dx ./ C1j_safe) .* (g_val - q_val .* U1j);
    end

    % Right Boundary (i=Nx)
    if strcmp(BC.Right.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Right.h, x_bnd_R, y, t_curr)'; 
        r_val = eval_bc_func(BC.Right.r, x_bnd_R, y, t_curr)'; 
        h_val(abs(h_val) < 1e-15) = 1;
        U_padded(Nx+1, 2:Ny+1) = r_val ./ h_val; % Row vector
        U_padded(Nx+2, 2:Ny+1) = 2*U_padded(Nx+1, 2:Ny+1) - U_padded(Nx, 2:Ny+1); 
    elseif strcmp(BC.Right.type, 'Neumann')
        q_val = eval_bc_func(BC.Right.q, x_bnd_R, y, t_curr)'; 
        g_val = eval_bc_func(BC.Right.g, x_bnd_R, y, t_curr)';
        UNxj = U(Nx, :);     
        CNxj = C(Nx, :);     
        UNx_1j = U(Nx-1, :); 
        CNxj_safe = CNxj; CNxj_safe(abs(CNxj_safe)<1e-15) = 1e-15;
        U_padded(Nx+2, 2:Ny+1) = UNx_1j + (2*dx ./ CNxj_safe) .* (g_val - q_val .* UNxj); 
    end

    % Bottom Boundary (j=1)
    if strcmp(BC.Bottom.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Bottom.h, x, y_bnd_B, t_curr); 
        r_val = eval_bc_func(BC.Bottom.r, x, y_bnd_B, t_curr); %
        h_val(abs(h_val) < 1e-15) = 1;
        U_padded(2:Nx+1, 2) = r_val ./ h_val;
        U_padded(2:Nx+1, 1) = 2*U_padded(2:Nx+1, 2) - U_padded(2:Nx+1, 3); 
    elseif strcmp(BC.Bottom.type, 'Neumann')
        q_val = eval_bc_func(BC.Bottom.q, x, y_bnd_B, t_curr);
        g_val = eval_bc_func(BC.Bottom.g, x, y_bnd_B, t_curr); 
        Ui1 = U(:, 1); % Column vector
        Ci1 = C(:, 1); % Column vector
        Ui2 = U(:, 2); % Column vector
        Ci1_safe = Ci1; Ci1_safe(abs(Ci1_safe)<1e-15) = 1e-15;
        U_padded(2:Nx+1, 1) = Ui2 - (2*dy ./ Ci1_safe) .* (g_val - q_val .* Ui1); 
    end

    % Top Boundary (j=Ny)
    if strcmp(BC.Top.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Top.h, x, y_bnd_T, t_curr);
        r_val = eval_bc_func(BC.Top.r, x, y_bnd_T, t_curr);
        h_val(abs(h_val) < 1e-15) = 1;
        U_padded(2:Nx+1, Ny+1) = r_val ./ h_val; % Column vector
        U_padded(2:Nx+1, Ny+2) = 2*U_padded(2:Nx+1, Ny+1) - U_padded(2:Nx+1, Ny); % Extrapolate ghost point
    elseif strcmp(BC.Top.type, 'Neumann')
        q_val = eval_bc_func(BC.Top.q, x, y_bnd_T, t_curr);
        g_val = eval_bc_func(BC.Top.g, x, y_bnd_T, t_curr);
        UiNy = U(:, Ny);    
        CiNy = C(:, Ny);     
        UiNy_1 = U(:, Ny-1); 
        CiNy_safe = CiNy; CiNy_safe(abs(CiNy_safe)<1e-15) = 1e-15;
        U_padded(2:Nx+1, Ny+2) = UiNy_1 + (2*dy ./ CiNy_safe) .* (g_val - q_val .* UiNy); 
    end

    term_x = (Cx_half(2:Nx+1, :) .* (U_padded(3:Nx+2, 2:Ny+1) - U_padded(2:Nx+1, 2:Ny+1)) / dx ...
            - Cx_half(1:Nx, :)   .* (U_padded(2:Nx+1, 2:Ny+1) - U_padded(1:Nx, 2:Ny+1)) / dx) / dx;

    term_y = (Cy_half(:, 2:Ny+1) .* (U_padded(2:Nx+1, 3:Ny+2) - U_padded(2:Nx+1, 2:Ny+1)) / dy ...
            - Cy_half(:, 1:Ny)   .* (U_padded(2:Nx+1, 2:Ny+1) - U_padded(2:Nx+1, 1:Ny)) / dy) / dy;

    Lu = term_x + term_y;

end

% =========================================================================
% Helper function to apply boundary conditions after Leapfrog step
% =========================================================================
function U = apply_boundary_conditions(U, BC, t_next, x, y, Nx, Ny, xL, xR, yL, yR)
    % Apply Dirichlet BCs directly. Neumann BCs were handled implicitly
    % via ghost points in the spatial operator calculation.

    x_bnd_L = repmat(xL, size(y)); 
    x_bnd_R = repmat(xR, size(y)); 
    y_bnd_B = repmat(yL, size(x)); 
    y_bnd_T = repmat(yR, size(x)); 

    % Left Boundary (i=1)
    if strcmp(BC.Left.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Left.h, x_bnd_L, y, t_next)'; 
        r_val = eval_bc_func(BC.Left.r, x_bnd_L, y, t_next)'; 
        h_val(abs(h_val) < 1e-15) = 1; 
        U(1, :) = r_val ./ h_val; 
    end

    % Right Boundary (i=Nx)
    if strcmp(BC.Right.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Right.h, x_bnd_R, y, t_next)'; 
        r_val = eval_bc_func(BC.Right.r, x_bnd_R, y, t_next)'; 
        h_val(abs(h_val) < 1e-15) = 1;
        U(Nx, :) = r_val ./ h_val; 
    end

    % Bottom Boundary (j=1)
    if strcmp(BC.Bottom.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Bottom.h, x, y_bnd_B, t_next); 
        r_val = eval_bc_func(BC.Bottom.r, x, y_bnd_B, t_next); 
        h_val(abs(h_val) < 1e-15) = 1;
        U(:, 1) = r_val ./ h_val; 
    end

    % Top Boundary (j=Ny)
    if strcmp(BC.Top.type, 'Dirichlet')
        h_val = eval_bc_func(BC.Top.h, x, y_bnd_T, t_next); 
        r_val = eval_bc_func(BC.Top.r, x, y_bnd_T, t_next); 
        h_val(abs(h_val) < 1e-15) = 1;
        U(:, Ny) = r_val ./ h_val; 
    end
end


% =========================================================================
% Helper function to evaluate boundary condition parameters (constants or functions)
% =========================================================================
function val = eval_bc_func(param, x_bnd, y_bnd, t)
    num_pts = numel(x_bnd);
    if isa(param, 'function_handle')
        try
            val = param(x_bnd, y_bnd, t);
        catch ME
            error('Error evaluating BC function handle: %s\n%s\nCheck function signature: @(x, y, t)', func2str(param), ME.message);
        end
    elseif isnumeric(param) && isscalar(param)
        val = repmat(param, num_pts, 1); % Repeat constant for all boundary points (column vector)
    else
        error('BC parameter must be a scalar constant or a function handle accepting (x, y, t).');
    end

    if isscalar(val)
        val = repmat(val, num_pts, 1);
    elseif ~iscolumn(val)
        if isrow(val) && numel(val) == num_pts
            val = val'; 
        elseif numel(val) == num_pts 
             warning('BC function handle %s returned an array of unexpected shape. Reshaping to column vector.', func2str(param));
             val = reshape(val, num_pts, 1);
        else 
             error('BC function handle %s returned an array of incorrect size (%d elements) for the boundary (%d points).', func2str(param), numel(val), num_pts);
        end
    elseif size(val, 1) ~= num_pts 
         error('BC function handle %s returned a column vector of incorrect size (%d rows) for the boundary (%d points).', func2str(param), size(val,1), num_pts);
    end

    if ~iscolumn(val) || size(val,1) ~= num_pts
         error('Failed to produce a column vector of size %d x 1 in eval_bc_func.', num_pts);
    end
end
