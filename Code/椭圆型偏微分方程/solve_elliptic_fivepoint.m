% Main script to define and solve a 2D elliptic PDE of the form
% -div(c * grad(u)) + a * u = f
% using a five-point finite difference method on a rectangular domain.
% =========================================================================
function solve_elliptic_fivepoint(xy_range,Nxy,equ_para,bc_type,bc_para)


    % -------------------------------------------------------------------------
    % 1. Define Domain and Grid
    % -------------------------------------------------------------------------
    x_left = xy_range(1);
    x_right = xy_range(2);
    y_left = xy_range(3);
    y_right = xy_range(4);
    Nx = Nxy(1); % Number of intervals in x-direction
    Ny = Nxy(2); % Number of intervals in y-direction
    
    % Calculate grid spacing
    dx = (x_right - x_left) / Nx;
    dy = (y_right - y_left) / Ny;
    
    % Create meshgrid (includes boundary points)
    % Node indices: i = 0...Nx, j = 0...Ny
    [X, Y] = meshgrid(linspace(x_left, x_right, Nx+1), ...
                      linspace(y_left, y_right, Ny+1));
    
    % -------------------------------------------------------------------------
    % 2. Define PDE Coefficients (can be constants or functions of x, y)
    % -------------------------------------------------------------------------
    % Diffusion coefficient c(x, y)
    c_func = equ_para.c; 
    % 
    % % Reaction coefficient a(x, y)
    a_func = equ_para.a; 
    % 
    % % Source term f(x, y)
    f_func = equ_para.f;
    
    % -------------------------------------------------------------------------
    % 3. Define Boundary Conditions
    % -------------------------------------------------------------------------
    % Specify type for each boundary: 'dirichlet' or 'neumann'
    if bc_type(1)==1
        left_type   = 'dirichlet';
    else
        left_type='neumann';
    end
    if bc_type(2)==1
        right_type  = 'dirichlet';
    else
        right_type  = 'neumann';
    end
    if bc_type(3)==1
        bottom_type='dirichlet';
    else
        bottom_type='neumann';
    end
    if bc_type(4)==1
        top_type    = 'dirichlet';
    else
        top_type    = 'neumann';
    end
    
    
    % --- Example: Homogeneous Dirichlet (u=0 on all boundaries) ---
    h_left_func   = bc_para.h_left; r_left_func   = bc_para.r_left;
    h_right_func  = bc_para.h_right; r_right_func  = bc_para.r_right;
    h_bottom_func = bc_para.h_bottom; r_bottom_func = bc_para.r_bottom;
    h_top_func    = bc_para.h_top; r_top_func    = bc_para.r_top;
    % Define dummy Neumann parameters (needed by solver function structure)
    q_left_func   = bc_para.q_left; g_left_func   =bc_para.g_left;
    q_right_func  = bc_para.q_right; g_right_func  = bc_para.g_left;
    q_bottom_func = bc_para.q_bottom; g_bottom_func = bc_para.g_bottom;
    q_top_func    = bc_para.q_top; g_top_func    = bc_para.g_top;
    % --- End Example ---
    
    % -------------------------------------------------------------------------
    % 4. Solve the PDE
    % -------------------------------------------------------------------------
    [u_numerical, X_solve, Y_solve] = solve_elliptic_pde_direct(...
        x_left, x_right, y_left, y_right, Nx, Ny, dx, dy, ... 
        c_func, a_func, f_func, ...                            
        left_type, right_type, bottom_type, top_type, ...      
        h_left_func, r_left_func, q_left_func, g_left_func, ... 
        h_right_func, r_right_func, q_right_func, g_right_func, ... 
        h_bottom_func, r_bottom_func, q_bottom_func, g_bottom_func, ... 
        h_top_func, r_top_func, q_top_func, g_top_func ...     
    );
    
    if ~isequal(X_solve, X) || ~isequal(Y_solve, Y)
        warning('Solver grid does not match input grid. Interpolating solution.');
        u_numerical = interp2(X_solve, Y_solve, u_numerical, X, Y, 'linear');
    end
    
    % -------------------------------------------------------------------------
    % 5. Plot Numerical Solution
    % -------------------------------------------------------------------------
    
    % figure('Name', 'Numerical Solution', 'Position', [100, 100, 1000, 450]);
    figure;
    % 3D Surface Plot
    subplot(1, 2, 1);
    surf(X, Y, u_numerical);
    shading interp;
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title('数值解（三维图）');
    axis tight;
    colorbar;
    view(3); % Standard 3D view
    
    % 2D Contour Plot
    subplot(1, 2, 2);
    contourf(X, Y, u_numerical, 20); % Filled contour plot with 20 levels
    xlabel('x');
    ylabel('y');
    title('数值解（二维图）');
    axis equal tight; % Ensure aspect ratio is 1 and axes fit data
    colorbar;
    grid on;
end

% =========================================================================
% solve_elliptic_pde_direct.m Function
% =========================================================================
function [u_numerical, X, Y] = solve_elliptic_pde_direct(...
    x_left, x_right, y_left, y_right, Nx, Ny, dx, dy, ... 
    c_func, a_func, f_func, ...                            
    left_type, right_type, bottom_type, top_type, ...      
    h_left_func, r_left_func, q_left_func, g_left_func, ... 
    h_right_func, r_right_func, q_right_func, g_right_func, ... 
    h_bottom_func, r_bottom_func, q_bottom_func, g_bottom_func, ... 
    h_top_func, r_top_func, q_top_func, g_top_func ...     
)

    N_total = (Nx + 1) * (Ny + 1);

    % Create meshgrid for calculations
    [X, Y] = meshgrid(linspace(x_left, x_right, Nx+1), ...
                      linspace(y_left, y_right, Ny+1));

    % Expected size for coefficient matrices
    expected_size = [Ny + 1, Nx + 1];

    % --- Evaluate PDE coefficients on the grid ---
    C_val = c_func(X, Y);
    if isscalar(C_val), C_val = C_val * ones(expected_size); end
    A_val = a_func(X, Y);
    if isscalar(A_val), A_val = A_val * ones(expected_size); end
    F_val = f_func(X, Y);
    if isscalar(F_val), F_val = F_val * ones(expected_size); end

    A_mat = spalloc(N_total, N_total, 5 * N_total);
    B_vec = zeros(N_total, 1);
    dx2 = dx^2; dy2 = dy^2;
    k_map = @(i, j) 1 + i + j * (Nx + 1);

    for j = 0:Ny
        for i = 0:Nx
            k = k_map(i, j); x = X(j+1, i+1); y = Y(j+1, i+1);
            is_left = (i == 0); is_right = (i == Nx);
            is_bottom = (j == 0); is_top = (j == Ny);

            coeff_center = A_val(j+1, i+1);
            coeff_left = 0; coeff_right = 0; coeff_bottom = 0; coeff_top = 0;
            rhs_val = F_val(j+1, i+1);

            if ~is_left,  c_imhalf_j = (C_val(j+1, i) + C_val(j+1, i+1)) / 2; coeff_center = coeff_center + c_imhalf_j / dx2; coeff_left   = -c_imhalf_j / dx2; end
            if ~is_right, c_iphalf_j = (C_val(j+1, i+2) + C_val(j+1, i+1)) / 2; coeff_center = coeff_center + c_iphalf_j / dx2; coeff_right  = -c_iphalf_j / dx2; end
            if ~is_bottom,c_i_jmhalf = (C_val(j, i+1) + C_val(j+1, i+1)) / 2; coeff_center = coeff_center + c_i_jmhalf / dy2; coeff_bottom = -c_i_jmhalf / dy2; end
            if ~is_top,   c_i_jphalf = (C_val(j+2, i+1) + C_val(j+1, i+1)) / 2; coeff_center = coeff_center + c_i_jphalf / dy2; coeff_top    = -c_i_jphalf / dy2; end

            boundary_applied = false;
            if is_left
                if strcmpi(left_type, 'dirichlet')
                    h_val = h_left_func(x, y); h_val = h_val(1); r_val = r_left_func(x, y); r_val = r_val(1);
                    if abs(h_val) < 1e-15, coeff_center=1; rhs_val=0; else, coeff_center=h_val; rhs_val=r_val; end
                    coeff_left=0; coeff_right=0; coeff_bottom=0; coeff_top=0; boundary_applied=true;
                elseif strcmpi(left_type, 'neumann')
                    q_val = q_left_func(x, y); q_val = q_val(1); g_val = g_left_func(x, y); g_val = g_val(1);
                    coeff_center = coeff_center + q_val / dx; rhs_val = rhs_val + g_val / dx;
                end
            end
            if is_right && ~boundary_applied
                if strcmpi(right_type, 'dirichlet')
                    h_val = h_right_func(x, y); h_val = h_val(1); r_val = r_right_func(x, y); r_val = r_val(1);
                    if abs(h_val) < 1e-15, coeff_center=1; rhs_val=0; else, coeff_center=h_val; rhs_val=r_val; end
                    coeff_left=0; coeff_right=0; coeff_bottom=0; coeff_top=0; boundary_applied=true;
                elseif strcmpi(right_type, 'neumann')
                    q_val = q_right_func(x, y); q_val = q_val(1); g_val = g_right_func(x, y); g_val = g_val(1);
                    coeff_center = coeff_center + q_val / dx; rhs_val = rhs_val - g_val / dx;
                end
            end
            if is_bottom && ~boundary_applied
                if strcmpi(bottom_type, 'dirichlet')
                    h_val = h_bottom_func(x, y); h_val = h_val(1); r_val = r_bottom_func(x, y); r_val = r_val(1);
                    if abs(h_val) < 1e-15, coeff_center=1; rhs_val=0; else, coeff_center=h_val; rhs_val=r_val; end
                    coeff_left=0; coeff_right=0; coeff_bottom=0; coeff_top=0; boundary_applied=true;
                elseif strcmpi(bottom_type, 'neumann')
                    q_val = q_bottom_func(x, y); q_val = q_val(1); g_val = g_bottom_func(x, y); g_val = g_val(1);
                    coeff_center = coeff_center + q_val / dy; rhs_val = rhs_val + g_val / dy;
                end
            end
            if is_top && ~boundary_applied
                if strcmpi(top_type, 'dirichlet')
                    h_val = h_top_func(x, y); h_val = h_val(1); r_val = r_top_func(x, y); r_val = r_val(1);
                    if abs(h_val) < 1e-15, coeff_center=1; rhs_val=0; else, coeff_center=h_val; rhs_val=r_val; end
                    coeff_left=0; coeff_right=0; coeff_bottom=0; coeff_top=0;
                elseif strcmpi(top_type, 'neumann')
                    q_val = q_top_func(x, y); q_val = q_val(1); g_val = g_top_func(x, y); g_val = g_val(1);
                    coeff_center = coeff_center + q_val / dy; rhs_val = rhs_val - g_val / dy;
                end
            end

            A_mat(k, k) = coeff_center; B_vec(k) = rhs_val;
            if coeff_left ~= 0,   A_mat(k, k_map(i-1, j)) = coeff_left;   end
            if coeff_right ~= 0,  A_mat(k, k_map(i+1, j)) = coeff_right;  end
            if coeff_bottom ~= 0, A_mat(k, k_map(i, j-1)) = coeff_bottom; end
            if coeff_top ~= 0,    A_mat(k, k_map(i, j+1)) = coeff_top;    end
        end
    end

    tol_sym = 1e-10 * nnz(A_mat);
    is_symmetric = nnz(A_mat - A_mat.') < tol_sym;
    tol_iter = 1e-8; maxit = min(N_total, 1000);
    flag = -1; U_vec = [];

    if is_symmetric
        is_diag_dominant = all(abs(diag(A_mat)) >= sum(abs(A_mat - spdiags(diag(A_mat),0,N_total,N_total)), 2));
        if is_diag_dominant
             solver_name = 'PCG';
             [U_vec, flag, relres, iter] = pcg(A_mat, B_vec, tol_iter, maxit);
        else
             solver_name = 'MINRES';
             [U_vec, flag, relres, iter] = minres(A_mat, B_vec, tol_iter, maxit);
        end
    else
        solver_name = 'BiCGSTAB';
        [U_vec, flag, relres, iter] = bicgstab(A_mat, B_vec, tol_iter, maxit);
    end

    if flag == 0
        solver_used = solver_name;
    else
        if flag ~= -1
        end
        U_vec = A_mat \ B_vec;
        solver_used = 'Direct (Fallback)';
    end

    u_numerical = reshape(U_vec, Ny + 1, Nx + 1);

end 
