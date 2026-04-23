% =========================================================================
% Solve 2D Multiscale Elliptic PDE using Finite Difference Method
% Equation: -div(c(x,y) * grad(u(x,y))) = f(x,y)
% Domain: [x_left, x_right] x [y_left, y_right] (Square)
% =========================================================================
function solve_mulell_dirfive(xy_range,Nxy,equ_para,bc_type,bc_para)

% --- User-Defined Parameters ---

% 1. Domain and Discretization
x_left = xy_range(1);
x_right = xy_range(2);
y_left = xy_range(3);
y_right = xy_range(4);
Nx = Nxy(1);
Ny = Nxy(2); 

% 2. Multiscale Coefficient Parameters
epsilon = equ_para.epsilon*(x_right-x_left);      % Microscale size (ratio to macroscale)
r_circle = equ_para.len;    % Radius of circular inclusions (must be < epsilon/2)
c_in = equ_para.c_in;           % Coefficient value inside circles
c_out = equ_para.c_ma;        % Coefficient value outside circles

% 3. Source Term Function f(x,y)
f_func = equ_para.f;

% 4. Boundary Conditions

% Left Boundary
if bc_type(1)==1
    bc.L.type = 'D';     % Type: 'D' or 'N'
else
    bc.L.type = 'N';     % Type: 'D' or 'N'
end

bc.L.h = bc_para.h_left;  
bc.L.r = bc_para.r_left;   
bc.L.q = bc_para.q_left;   
bc.L.g = bc_para.g_left;  

% Right Boundary 
if bc_type(2)==1
    bc.R.type = 'D';
else
    bc.R.type = 'N';
end

bc.R.h = bc_para.h_right;
bc.R.r = bc_para.r_right;
bc.R.q = bc_para.q_right;
bc.R.g = bc_para.g_right;

% Bottom Boundary
if bc_type(3)==1
    bc.B.type = 'D';
else
    bc.B.type = 'N';
end

bc.B.h = bc_para.h_bottom;
bc.B.r = bc_para.r_bottom;
bc.B.q = bc_para.q_bottom;
bc.B.g = bc_para.g_bottom;

% Top Boundary
if bc_type(4)==1
    bc.T.type='D';
else
    bc.T.type='N';
end
bc.T.h = bc_para.h_top;
bc.T.r = bc_para.r_top;
bc.T.q = bc_para.q_top;
bc.T.g = bc_para.g_top;


% --- Setup Grid ---
dx = (x_right - x_left) / (Nx - 1);
dy = (y_right - y_left) / (Ny - 1);
x = linspace(x_left, x_right, Nx);
y = linspace(y_left, y_right, Ny);
[X, Y] = meshgrid(x, y);
N_total = Nx * Ny; 

% --- Define Coefficient c(x,y) ---
C = define_coefficient(X, Y, epsilon, r_circle, c_in, c_out, x_left, y_left, x_right, y_right);

% --- Define Source Term f(x,y) on Grid ---
F_eval = f_func(X, Y);
if isscalar(F_eval)
    F_grid = ones(Ny, Nx) * F_eval;
else
    F_grid = F_eval;
    if ~isequal(size(F_grid), [Ny, Nx])
        error('Size of F_grid returned by f_func does not match grid dimensions [Ny, Nx].');
    end
end

% --- Assemble Linear System A*U = F_rhs ---
[A, F_rhs] = assemble_system(Nx, Ny, dx, dy, C, F_grid, bc, x, y, x_left, y_left, x_right, y_right);

% --- Solve the Linear System using the Direct Solver (\) ---
try
    U_vec = A \ F_rhs;
catch ME_direct
    error('Direct solver failed. Check matrix properties or system size.');
end

% --- Reshape Solution Vector to Grid ---
if isempty(U_vec)
    error('Solver failed to produce a solution vector.');
else
    U_sol = reshape(U_vec, Ny, Nx);
end

% --- Plot Results ---
plot_results(X, Y, U_sol, C, x, y);
end

% =========================================================================
% Helper Functions
% =========================================================================

function C = define_coefficient(X, Y, epsilon, r_circle, c_in, c_out, x_left, y_left, x_right, y_right)
    [Ny, Nx] = size(X);
    C = ones(Ny, Nx) * c_out;

    % Iterate through each grid point
    for j = 1:Ny
        for i = 1:Nx
            x_coord = X(j, i);
            y_coord = Y(j, i);

            % Find the microcell index (0-based for calculation)
            % Use small tolerance to handle points exactly on cell boundaries
            ix_calc = floor((x_coord - x_left + 1e-9*epsilon) / epsilon);
            iy_calc = floor((y_coord - y_left + 1e-9*epsilon) / epsilon);

            % Ensure indices stay within bounds (can happen at right/top edge)
            max_ix = floor((x_right - x_left - 1e-9*epsilon) / epsilon);
            max_iy = floor((y_right - y_left - 1e-9*epsilon) / epsilon);
            ix_calc = min(ix_calc, max_ix);
            iy_calc = min(iy_calc, max_iy);

            % Calculate the center of the current microcell
            xc = x_left + (ix_calc + 0.5) * epsilon;
            yc = y_left + (iy_calc + 0.5) * epsilon;

            % Check if the point is inside the circle within this microcell
            dist_sq = (x_coord - xc)^2 + (y_coord - yc)^2;
            if dist_sq <= r_circle^2 + 1e-9 
                C(j, i) = c_in; 
            end
        end
    end
end

% -------------------------------------------------------------------------

function [A, F_rhs] = assemble_system(Nx, Ny, dx, dy, C, F_grid, bc, x, y, x_left, y_left, x_right, y_right)
    N_total = Nx * Ny;
    A = spalloc(N_total, N_total, 5 * N_total);
    F_rhs = zeros(N_total, 1);

    dx2 = dx^2;
    dy2 = dy^2;

    harmonic_mean = @(c1, c2) ( (c1+c2) == 0 ) * 0 + ( (c1+c2) ~= 0 ) * (2 * c1 * c2 / (c1 + c2));

    for j = 1:Ny
        for i = 1:Nx 
            k = (j - 1) * Nx + i; 
            x_coord = x(i); 
            y_coord = y(j); 

            % --- Apply Boundary Conditions ---
            is_dirichlet = false;
            h_eff = 0; 
            r_eff = 0;

            % Check boundaries (corners prioritize Dirichlet if specified on either edge)
            if i == 1 && upper(bc.L.type) == 'D'
                h_eff = evaluate_bc_param(bc.L.h, x_coord, y_coord);
                r_eff = evaluate_bc_param(bc.L.r, x_coord, y_coord);
                is_dirichlet = true;
            elseif i == Nx && upper(bc.R.type) == 'D' 
                h_eff = evaluate_bc_param(bc.R.h, x_coord, y_coord);
                r_eff = evaluate_bc_param(bc.R.r, x_coord, y_coord);
                is_dirichlet = true;
            elseif j == 1 && upper(bc.B.type) == 'D'
                h_eff = evaluate_bc_param(bc.B.h, x_coord, y_coord);
                r_eff = evaluate_bc_param(bc.B.r, x_coord, y_coord);
                is_dirichlet = true;
            elseif j == Ny && upper(bc.T.type) == 'D'
                h_eff = evaluate_bc_param(bc.T.h, x_coord, y_coord);
                r_eff = evaluate_bc_param(bc.T.r, x_coord, y_coord);
                is_dirichlet = true;
            end

            % If it's a Dirichlet node:
            if is_dirichlet
                if abs(h_eff) < eps
                   A(k, k) = 1;
                   r_eff = 0;   
                   warning('Dirichlet BC h=0 encountered at (%.2f, %.2f). Setting A(k,k)=1, u=0.', x_coord, y_coord);
                else
                   A(k, k) = h_eff; % Set diagonal entry to h
                end
                F_rhs(k) = r_eff;
                continue; 
            end

            % --- Interior or Neumann/Robin Boundary Node ---
            term_center = 0; 
            rhs_contribution = F_grid(j, i); 

            % X-direction terms: -d/dx (c * du/dx)
            if i > 1 
                c_imh = harmonic_mean(C(j,i), C(j,i-1)); 
                term_left = c_imh / dx2;
                term_center = term_center + term_left; % Add to diagonal
                k_left = k - 1; % Global index of left neighbor
                A(k, k_left) = -term_left; % Off-diagonal entry
            else 
                if upper(bc.L.type) ~= 'N'
                    error('BC Type Mismatch: Left boundary node (%d,%d) is not Dirichlet but BC type is not N.',i,j);
                end
                q_val = evaluate_bc_param(bc.L.q, x_coord, y_coord); 
                g_val = evaluate_bc_param(bc.L.g, x_coord, y_coord);
                c_iph = harmonic_mean(C(j,1), C(j,2)); 
                term_center = term_center + 2*C(j,i)/dx2; 
                A(k,k) = A(k,k) + 2*q_val/dx; 
                rhs_contribution = rhs_contribution + 2*g_val/dx; 
                 term_center = 0; 
                 rhs_contribution = F_grid(j,i); 
                 c_iph = harmonic_mean(C(j,i), C(j,i+1));
                 term_right = c_iph / dx2;
                 term_center = term_center + term_right;
                 A(k, k+1) = -term_right;
                 c_iph = harmonic_mean(C(j,1), C(j,2));
                 term_right = c_iph / dx2;
                 term_center = term_center + term_right;
                 A(k, k+1) = -term_right;
                 term_center = term_center + (2*q_val*dx/ (2*dx*C(j,i))) * (C(j,i)/dx2); 
                 c_iph = harmonic_mean(C(j,1), C(j,2));
                 term_center = term_center + c_iph / dx2; 
                 A(k, k+1) = -c_iph / dx2; 
                 term_center = term_center + C(j,i)/dx2 + q_val/dx; 
                 rhs_contribution = rhs_contribution + g_val/dx; 

            end

            if i < Nx 
                c_iph = harmonic_mean(C(j,i), C(j,i+1));
                term_right = c_iph / dx2;
                term_center = term_center + term_right; 
                k_right = k + 1; 
                A(k, k_right) = -term_right;
            else 
                 if upper(bc.R.type) ~= 'N'
                     error('BC Type Mismatch: Right boundary node (%d,%d) is not Dirichlet but BC type is not N.',i,j);
                 end
                 q_val = evaluate_bc_param(bc.R.q, x_coord, y_coord);
                 g_val = evaluate_bc_param(bc.R.g, x_coord, y_coord);
                 c_imh = harmonic_mean(C(j,Nx), C(j,Nx-1)); 
                 term_center = term_center + c_imh / dx2; 
                 A(k, k-1) = -c_imh / dx2;
                 term_center = term_center + C(j,i)/dx2 + q_val/dx; 
                 rhs_contribution = rhs_contribution + g_val/dx; 
            end

             if j > 1 
                c_jmh = harmonic_mean(C(j,i), C(j-1,i));
                term_bottom = c_jmh / dy2;
                term_center = term_center + term_bottom; 
                k_bottom = k - Nx; 
                A(k, k_bottom) = -term_bottom; 
            else 
                 if upper(bc.B.type) ~= 'N'
                     error('BC Type Mismatch: Bottom boundary node (%d,%d) is not Dirichlet but BC type is not N.',i,j);
                 end
                 q_val = evaluate_bc_param(bc.B.q, x_coord, y_coord);
                 g_val = evaluate_bc_param(bc.B.g, x_coord, y_coord);
                 c_jph = harmonic_mean(C(1,i), C(2,i)); 
                 term_center = term_center + c_jph / dy2; 
                 A(k, k+Nx) = -c_jph / dy2; 
                 term_center = term_center + C(j,i)/dy2 + q_val/dy; 
                 rhs_contribution = rhs_contribution + g_val/dy; 
            end

            if j < Ny 
                c_jph = harmonic_mean(C(j,i), C(j+1,i)); 
                term_top = c_jph / dy2;
                term_center = term_center + term_top;
                k_top = k + Nx; 
                A(k, k_top) = -term_top; 
            else 
                 if upper(bc.T.type) ~= 'N'
                     error('BC Type Mismatch: Top boundary node (%d,%d) is not Dirichlet but BC type is not N.',i,j);
                 end
                 q_val = evaluate_bc_param(bc.T.q, x_coord, y_coord);
                 g_val = evaluate_bc_param(bc.T.g, x_coord, y_coord);
                 c_jmh = harmonic_mean(C(Ny,i), C(Ny-1,i)); 
                 term_center = term_center + c_jmh / dy2;
                 A(k, k-Nx) = -c_jmh / dy2; 
                 term_center = term_center + C(j,i)/dy2 + q_val/dy; 
                 rhs_contribution = rhs_contribution + g_val/dy; 
            end

            A(k, k) = term_center;
            F_rhs(k) = rhs_contribution;

        end 
    end 
end

% -------------------------------------------------------------------------

function val = evaluate_bc_param(param, x_coord, y_coord)
    if isa(param, 'function_handle')
        val_eval = param(x_coord, y_coord);
        if ~isscalar(val_eval) && ~isempty(val_eval)
             val = val_eval(1);
             elseif isempty(val_eval)
             error('BC function handle returned empty value for input (%.4f, %.4f).', x_coord, y_coord);
        else 
             val = val_eval;
        end
    elseif isscalar(param) 
        val = param;
    else 
         error('BC parameter is neither a scalar nor a function handle for coordinates (%.4f, %.4f)', x_coord, y_coord);
    end

    if ~isscalar(val)
        error('BC parameter evaluation did not result in a scalar value for coordinates (%.4f, %.4f)', x_coord, y_coord);
    end
end

% -------------------------------------------------------------------------

function plot_results(X, Y, U_sol, C, x, y)
    figure; 

    % 1. 3D Surface Plot of Solution U
    subplot(1, 2, 1); 
    surf(X, Y, U_sol); 
    colormap(gca, jet); 
    colorbar; 
    xlabel('x'); ylabel('y'); zlabel('u');
    title('参考解（三维图）');
    axis tight; 
    shading interp;
    view(3); 

    % 2. 2D Contour Plot of Solution U
    subplot(1, 2, 2); 
    contourf(X, Y, U_sol, 20);
    % colormap(gca, jet); 
    colorbar; 
    xlabel('x'); ylabel('y');
    title('参考解（二维图）');
    axis equal tight;

    % 3. 2D Plot of Coefficient C
    % subplot(1, 3, 3); 
    % imagesc(x, y, C); 
    % colormap(gca, copper); 
    % colorbar; 
    % axis xy; 
    % xlabel('x'); ylabel('y');
    % title('Coefficient c(x,y)');
    % axis equal tight; 

end
