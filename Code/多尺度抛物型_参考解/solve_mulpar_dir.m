function solve_mulpar_dir(xyt_range,Nxyt,equ_para,u0,bc_type,bc_para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve 2D Multiscale Parabolic PDE using Crank-Nicolson Finite Difference
% PDE: d*u_t - div(c*grad(u)) = f
% Domain: [xa, xb] x [ya, yb] (Square: xb-xa = yb-ya)
% Microstructure: Circles within microcells defined by epsilon
% Boundary Conditions: Dirichlet (h*u = r) or Neumann (n*c*grad(u) + q*u = g)
% MODIFIED: Uses only the direct solver (\) for the linear system.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. --- User-Defined Parameters ---

% Domain
xa = xyt_range(1); xb = xyt_range(2);
ya = xyt_range(3); yb = xyt_range(4); % Must have yb-ya = xb-xa

% Time
t_start = xyt_range(5); t_end = xyt_range(6);

% Discretization
dx = (xb-xa)/Nxyt(1); % Spatial step in x
dy = (yb-ya)/Nxyt(2);    % Spatial step in y (must be same as dx for this microstructure logic)
dt = (t_end-t_start)/Nxyt(3);  % Time step

% Microstructure Parameters
epsilon = equ_para.epsilon*(xb-xa);
circle_abs_radius = equ_para.len;

% Coefficients (Constant values)
c_in = equ_para.c_in;      % 'c' inside circles
c_out = equ_para.c_ma;    % 'c' outside circles
d_in = equ_para.d_in;      % 'd' inside circles
d_out = equ_para.d_ma;     % 'd' outside circles

% Source Term f(x, y, t)
f_func = equ_para.f;

% Initial Condition u(x, y, t_start)
u_initial_func = u0;

% --- Boundary Conditions ---
% Type: 'D' for Dirichlet (h*u = r)
%       'N' for Neumann (n*c*grad(u) + q*u = g)
if bc_type(1)==1
    bc_para.g_left=@(x,y)0;
    bc_para.q_left=@(x,y)0;
else
    bc_para.h_left=@(x,y)0;
    bc_para.r_left=@(x,y)0;
end
if bc_type(2)==1
    bc_para.g_right=@(x,y)0;
    bc_para.q_right=@(x,y)0;
else
    bc_para.h_right=@(x,y)0;
    bc_para.r_right=@(x,y)0;
end
if bc_type(3)==1
    bc_para.g_bottom=@(x,y)0;
    bc_para.q_bottom=@(x,y)0;
else
    bc_para.h_bottom=@(x,y)0;
    bc_para.r_bottom=@(x,y)0;
end
if bc_type(4)==1
    bc_para.g_top=@(x,y)0;
    bc_para.q_top=@(x,y)0;
else
    bc_para.h_top=@(x,y)0;
    bc_para.r_top=@(x,y)0;
end
% Left Boundary (x = xa)
if bc_type(1)==1
    bc_left.type = 'D';
else
    bc_left.type = 'N';
end

bc_left.h = bc_para.h_left;       % h(x,y,t) for Dirichlet
bc_left.r = bc_para.r_left;       % r(x,y,t) for Dirichlet
bc_left.q = bc_para.q_left;       % q(x,y,t) for Neumann
bc_left.g = bc_para.g_left;       % g(x,y,t) for Neumann

% Right Boundary (x = xb)
if bc_type(2)==1
    bc_right.type = 'D';
else
    bc_right.type = 'N';
end
bc_right.h = bc_para.h_right;      % h(x,y,t) for Dirichlet
bc_right.r = bc_para.r_right;      % r(x,y,t) for Dirichlet
bc_right.q = bc_para.q_right;      % q(x,y,t) for Neumann
bc_right.g = bc_para.g_right;      % g(x,y,t) for Neumann

% Bottom Boundary (y = ya)
if bc_type(3)==1
    bc_bottom.type = 'D';
else
    bc_bottom.type = 'N';
end

bc_bottom.h =bc_para.h_bottom;     % h(x,y,t) for Dirichlet
bc_bottom.r = bc_para.r_bottom;     % r(x,y,t) for Dirichlet
bc_bottom.q = bc_para.q_bottom;     % q(x,y,t) for Neumann
bc_bottom.g = bc_para.g_bottom;     % g(x,y,t) for Neumann

% Top Boundary (y = yb)
if bc_type(4)==1
    bc_top.type = 'D';
else
    bc_top.type = 'N';
end

bc_top.h = bc_para.h_top;        % h(x,y,t) for Dirichlet
bc_top.r = bc_para.r_top;        % r(x,y,t) for Dirichlet
bc_top.q = bc_para.q_top;        % q(x,y,t) for Neumann % Note: Original had q=1, g=0 for Neumann Top BC
bc_top.g = bc_para.g_top;        % g(x,y,t) for Neumann

%% 2. --- Grid Setup ---

x = xa:dx:xb;
y = ya:dy:yb;
t = t_start:dt:t_end;

Nx = length(x);
Ny = length(y);
Nt = length(t);

N = Nx * Ny; % Total number of grid points

[X, Y] = meshgrid(x, y);
X_vec = X(:);
Y_vec = Y(:);

%% 3. --- Determine Coefficients c and d based on Microstructure ---

C_mat = zeros(Ny, Nx);
D_mat = zeros(Ny, Nx);

% Pre-calculate squared absolute radius for efficiency
circle_abs_radius_sq = circle_abs_radius^2;

% Check if radius is valid
if circle_abs_radius > epsilon / 2.0
    warning('circle_abs_radius (%.4f) is larger than microcell half-size (%.4f). Circles may overlap or extend beyond cells.', circle_abs_radius, epsilon/2.0);
end


for j = 1:Ny % Loop over y
    for i = 1:Nx % Loop over x
        % Find the center of the microcell containing (x(i), y(j))
        cell_ix = floor((x(i) - xa) / epsilon);
        cell_iy = floor((y(j) - ya) / epsilon);

        % Adjust cell index for points exactly on the right/top boundary of a cell
        if abs(mod(x(i) - xa, epsilon)) < 1e-10 && x(i) > xa
            cell_ix = cell_ix - 1;
        end
         if abs(mod(y(j) - ya, epsilon)) < 1e-10 && y(j) > ya
            cell_iy = cell_iy - 1;
        end

        center_x = xa + (cell_ix + 0.5) * epsilon;
        center_y = ya + (cell_iy + 0.5) * epsilon;

        dist_sq = (x(i) - center_x)^2 + (y(j) - center_y)^2;

        is_inside = dist_sq <= circle_abs_radius_sq;

        if is_inside
            C_mat(j, i) = c_in;
            D_mat(j, i) = d_in;
        else
            C_mat(j, i) = c_out;
            D_mat(j, i) = d_out;
        end
    end
end

% --- Coefficient values at half-grid points (using arithmetic mean) ---
Cx_half = zeros(Ny, Nx-1); % c at (i+1/2, j)
for j = 1:Ny
    Cx_half(j,:) = (C_mat(j, 1:Nx-1) + C_mat(j, 2:Nx)) / 2;
end

Cy_half = zeros(Ny-1, Nx); % c at (i, j+1/2)
for i = 1:Nx
    Cy_half(:,i) = (C_mat(1:Ny-1, i) + C_mat(2:Ny, i)) / 2;
end

C_vec = reshape(C_mat, N, 1);
D_vec = reshape(D_mat, N, 1);
D_inv_dt = spdiags(D_vec / dt, 0, N, N); % Diagonal matrix D/dt

%% 4. --- Assemble Spatial Discretization Matrix L (for -div(c*grad(u))) ---

dx2 = dx^2;
dy2 = dy^2;

L = sparse(N, N);

% Loop over interior grid points (j=2:Ny-1, i=2:Nx-1)
for j = 2:Ny-1
    for i = 2:Nx-1
        k = i + (j-1)*Nx;

        c_e = Cx_half(j, i);  
        c_w = Cx_half(j, i-1); 
        c_n = Cy_half(j, i);   
        c_s = Cy_half(j-1, i); 

        % Finite difference approximation for -div(c*grad(u))
        L(k, k)   = (c_e + c_w)/dx2 + (c_n + c_s)/dy2;
        L(k, k+1) = -c_e / dx2; 
        L(k, k-1) = -c_w / dx2;
        L(k, k+Nx) = -c_n / dy2;
        L(k, k-Nx) = -c_s / dy2;
    end
end

%% 5. --- Initialize Solution Vector ---
U_init_raw = u_initial_func(X_vec, Y_vec);

if isscalar(U_init_raw)
    U = ones(N, 1) * U_init_raw;
elseif numel(U_init_raw) ~= N
     error('Initial condition function returned %d elements, expected N = %d.', numel(U_init_raw), N);
else
    U = reshape(U_init_raw, N, 1); 
end

if ~isequal(size(U), [N, 1])
    error('Initialization failed: U has size %s, expected [%d, 1].', mat2str(size(U)), N);
end


%% 6. --- Time Stepping Loop (Crank-Nicolson) ---

RHS_f = zeros(N, 1); % Preallocate source term vector

for n = 1:Nt-1 
    t_current = t(n);
    t_next = t(n+1);
    t_half = (t_current + t_next) / 2;

    if ~isequal(size(U), [N, 1])
        error('Error before RHS calculation at step %d: U has size %s, expected [%d, 1].', ...
              n, mat2str(size(U)), N);
    end

    % --- Define Crank-Nicolson matrices ---
    A = 2*D_inv_dt + L;
    B = 2*D_inv_dt - L;

    % --- Evaluate source term F at t_half ---
    F_half_raw = f_func(X_vec, Y_vec, t_half);

    if isscalar(F_half_raw)
        F_half = ones(N, 1) * F_half_raw;
    elseif numel(F_half_raw) ~= N
        error('Source function f_func returned %d elements at step %d, expected N = %d.', ...
              numel(F_half_raw), n, N);
    else
        F_half = reshape(F_half_raw, N, 1);
    end

    % --- Calculate Right Hand Side (RHS) ---
    RHS_f = 2 * F_half; % Scaled source term contribution
    RHS = B * U + RHS_f; % B*U(n) + 2*F(n+1/2)

    % Check RHS shape
     if ~isequal(size(RHS), [N, 1])
        error('Error after RHS calculation at step %d: RHS has size %s, expected [%d, 1]. U size was %s.', ...
              n, mat2str(size(RHS)), N, mat2str(size(U)));
    end

    % --- Create copies to modify for boundary conditions ---
    A_mod = A;
    RHS_mod = RHS;

    % --- Apply Boundary Conditions (Modify A_mod and RHS_mod) ---

    % --- Left Boundary  ---
    i = 1;
    current_x = x(i);
    for j = 1:Ny % Loop through y-indices
        k = i + (j-1)*Nx; % Global index
        current_y = y(j);
        if bc_left.type == 'D'
            h_val = bc_left.h(current_x, current_y, t_next); 
            r_val = bc_left.r(current_x, current_y, t_next); 
            if h_val ~= 0
                A_mod(k, :) = 0;    
                A_mod(k, k) = h_val; 
                RHS_mod(k) = r_val;
            else
                 A_mod(k, :) = 0;
                 A_mod(k, k) = 1; 
                 RHS_mod(k) = 0;
                 if abs(r_val) > 1e-10
                     warning('Dirichlet BC h=0 but r~=0 at x=xa, y=%.2f, t=%.2f. Enforcing u=0.', current_y, t_next);
                 end
            end
        elseif bc_left.type == 'N' 
            q_val = bc_left.q(current_x, current_y, t_next); 
            g_val = bc_left.g(current_x, current_y, t_next); 
            c_e = Cx_half(j, i); % c at (1.5, j)

            A_mod(k,:) = 0; 
            A_mod(k,k) = 2*D_vec(k)/dt; 

            A_mod(k, k)   = A_mod(k,k) + (c_e/dx2 + q_val/dx); 
            A_mod(k, k+1) = A_mod(k,k+1) - c_e/dx2; 

            % y-derivative part (if not a corner node)
            if j > 1 && j < Ny
                c_n = Cy_half(j, i);   % c at (1, j+1/2)
                c_s = Cy_half(j-1, i); % c at (1, j-1/2)
                A_mod(k, k)    = A_mod(k, k) + (c_n + c_s)/dy2; 
                A_mod(k, k+Nx) = A_mod(k, k+Nx) - c_n/dy2; 
                A_mod(k, k-Nx) = A_mod(k, k-Nx) - c_s/dy2; 
            elseif j == 1
                c_n = Cy_half(j, i); % c at (1, 1.5)
                 A_mod(k, k)    = A_mod(k, k) + c_n/dy2; 
                 A_mod(k, k+Nx) = A_mod(k, k+Nx) - c_n/dy2; 
            elseif j == Ny 
                 c_s = Cy_half(j-1, i); 
                 A_mod(k, k)    = A_mod(k, k) + c_s/dy2; 
                 A_mod(k, k-Nx) = A_mod(k, k-Nx) - c_s/dy2; 
            end

            RHS_mod(k) = RHS_mod(k) + 2*g_val/dx;
        end
    end

    % --- Right Boundary---
    i = Nx;
    current_x = x(i);
     for j = 1:Ny 
        k = i + (j-1)*Nx;
        current_y = y(j);
        if bc_right.type == 'D'
            h_val = bc_right.h(current_x, current_y, t_next);
            r_val = bc_right.r(current_x, current_y, t_next);
             if h_val ~= 0
                A_mod(k, :) = 0;
                A_mod(k, k) = h_val;
                RHS_mod(k) = r_val;
             else
                 A_mod(k, :) = 0;
                 A_mod(k, k) = 1;
                 RHS_mod(k) = 0;
                 if abs(r_val) > 1e-10
                     warning('Dirichlet BC h=0 but r~=0 at x=xb, y=%.2f, t=%.2f. Enforcing u=0.', current_y, t_next);
                 end
             end
        elseif bc_right.type == 'N' 
             q_val = bc_right.q(current_x, current_y, t_next);
             g_val = bc_right.g(current_x, current_y, t_next);
             c_w = Cx_half(j, i-1); 

             A_mod(k,:) = 0;
             A_mod(k,k) = 2*D_vec(k)/dt; 

             A_mod(k, k)   = A_mod(k,k) + (c_w/dx2 + q_val/dx); 
             A_mod(k, k-1) = A_mod(k,k-1) - c_w/dx2; 

             if j > 1 && j < Ny
                 c_n = Cy_half(j, i);  
                 c_s = Cy_half(j-1, i);
                 A_mod(k, k)    = A_mod(k, k) + (c_n + c_s)/dy2;
                 A_mod(k, k+Nx) = A_mod(k, k+Nx) - c_n/dy2;
                 A_mod(k, k-Nx) = A_mod(k, k-Nx) - c_s/dy2;
             elseif j == 1
                 c_n = Cy_half(j, i); 
                 A_mod(k, k)    = A_mod(k, k) + c_n/dy2;
                 A_mod(k, k+Nx) = A_mod(k, k+Nx) - c_n/dy2;
             elseif j == Ny
                 c_s = Cy_half(j-1, i); 
                 A_mod(k, k)    = A_mod(k, k) + c_s/dy2;
                 A_mod(k, k-Nx) = A_mod(k, k-Nx) - c_s/dy2;
             end

             RHS_mod(k) = RHS_mod(k) + 2*g_val/dx;
        end
    end

    % --- Bottom Boundary---
    j = 1;
    current_y = y(j);
    for i = 2:Nx-1 
        k = i + (j-1)*Nx; 
        current_x = x(i);
        if bc_bottom.type == 'D'
            h_val = bc_bottom.h(current_x, current_y, t_next);
            r_val = bc_bottom.r(current_x, current_y, t_next);
            if h_val ~= 0
                A_mod(k, :) = 0;
                A_mod(k, k) = h_val;
                RHS_mod(k) = r_val;
            else
                 A_mod(k, :) = 0;
                 A_mod(k, k) = 1;
                 RHS_mod(k) = 0;
                 if abs(r_val) > 1e-10
                     warning('Dirichlet BC h=0 but r~=0 at y=ya, x=%.2f, t=%.2f. Enforcing u=0.', current_x, t_next);
                 end
            end
        elseif bc_bottom.type == 'N' 
             q_val = bc_bottom.q(current_x, current_y, t_next);
             g_val = bc_bottom.g(current_x, current_y, t_next);
             c_n = Cy_half(j, i); % c at (i, 1.5)

             A_mod(k,:) = 0; 
             A_mod(k,k) = 2*D_vec(k)/dt; 

             A_mod(k, k)    = A_mod(k,k) + (c_n/dy2 + q_val/dy);
             A_mod(k, k+Nx) = A_mod(k,k+Nx) - c_n/dy2; 

             c_e = Cx_half(j, i);   % c at (i+1/2, 1)
             c_w = Cx_half(j, i-1); % c at (i-1/2, 1)
             A_mod(k, k)   = A_mod(k, k) + (c_e + c_w)/dx2;
             A_mod(k, k+1) = A_mod(k, k+1) - c_e/dx2;
             A_mod(k, k-1) = A_mod(k, k-1) - c_w/dx2;

             % Add g term modification to RHS_mod
             RHS_mod(k) = RHS_mod(k) + 2*g_val/dy;
        end
    end

    % --- Top Boundary---
    j = Ny;
    current_y = y(j);
     for i = 2:Nx-1 
        k = i + (j-1)*Nx; 
        current_x = x(i);
        if bc_top.type == 'D'
            h_val = bc_top.h(current_x, current_y, t_next);
            r_val = bc_top.r(current_x, current_y, t_next);
            if h_val ~= 0
                A_mod(k, :) = 0;
                A_mod(k, k) = h_val;
                RHS_mod(k) = r_val;
            else
                 A_mod(k, :) = 0;
                 A_mod(k, k) = 1;
                 RHS_mod(k) = 0;
                 if abs(r_val) > 1e-10
                     warning('Dirichlet BC h=0 but r~=0 at y=yb, x=%.2f, t=%.2f. Enforcing u=0.', current_x, t_next);
                 end
            end
        elseif bc_top.type == 'N' 
             q_val = bc_top.q(current_x, current_y, t_next);
             g_val = bc_top.g(current_x, current_y, t_next);
             c_s = Cy_half(j-1, i); % c at (i, Ny-0.5)

             A_mod(k,:) = 0; 
             A_mod(k,k) = 2*D_vec(k)/dt; 

             A_mod(k, k)    = A_mod(k,k) + (c_s/dy2 + q_val/dy); 
             A_mod(k, k-Nx) = A_mod(k,k-Nx) - c_s/dy2; 

              c_e = Cx_half(j, i);   
              c_w = Cx_half(j, i-1);
             A_mod(k, k)   = A_mod(k, k) + (c_e + c_w)/dx2;
             A_mod(k, k+1) = A_mod(k, k+1) - c_e/dx2;
             A_mod(k, k-1) = A_mod(k, k-1) - c_w/dx2;

             RHS_mod(k) = RHS_mod(k) + 2*g_val/dy;
        end
    end

    if ~isequal(size(A_mod), [N, N])
         error('Error before solver at step %d: A_mod has size %s, expected [%d, %d].', ...
              n, mat2str(size(A_mod)), N, N);
    end
    if ~isequal(size(RHS_mod), [N, 1])
         error('Error before solver at step %d: RHS_mod has size %s, expected [%d, 1].', ...
              n, mat2str(size(RHS_mod)), N);
    end

    try
        U_new_raw = A_mod \ RHS_mod;
    catch ME_direct
         error('Direct solver failed at step %d (t=%.3f): %s', n, t_next, ME_direct.message);
    end

    U_new = reshape(U_new_raw, N, 1);

     if ~isequal(size(U_new), [N, 1])
        error('Error after solver at step %d: U_new has size %s, expected [%d, 1]. Raw solver output size was %s.', ...
              n, mat2str(size(U_new)), N, mat2str(size(U_new_raw)));
    end

    U = U_new;

end

%% 7. --- Final Visualization ---

try
     if ~isequal(size(U), [N, 1])
        warning('Final solution vector U has size %s, expected [%d, 1]. Attempting reshape anyway.', ...
                 mat2str(size(U)), N);
    end

    U_final = reshape(U, Ny, []);
    if size(U_final, 2) ~= Nx
        warning('Final reshaped matrix column count (%d) does not match calculated Nx (%d). Check grid/solver consistency.', size(U_final, 2), Nx);
        if numel(U) == N
             U_final = reshape(U, Ny, Nx);
             fprintf('Reshaped using explicit [Ny, Nx] as fallback.\n');
        else
             error('Final solution vector U has %d elements, expected %d (Nx*Ny = %d*%d). Cannot reshape reliably.', numel(U), N, Nx, Ny);
        end
    end
catch ME
    error('Error reshaping final solution U (size %d x %d) into [%d x %d] matrix: %s\nCheck if N = Nx * Ny holds (%d vs %d*%d=%d).',...
          size(U,1), size(U,2), Ny, Nx, ME.message, numel(U), Nx, Ny, N);
end


% --- Create Plots ---
figure; % Wide figure

% Subplot 1: 3D Surface plot of the final solution
subplot(1, 2, 1);
surf(X, Y, U_final); 
title('参考解（三维）');
xlabel('x'); ylabel('y'); zlabel('u');
colorbar; 
axis tight; 
shading interp;
view(3); 

% Subplot 2: 2D Contour plot of the final solution
subplot(1, 2, 2);
contourf(X, Y, U_final, 20); % Filled contour plot with 20 levels
title('参考解（二维）');
xlabel('x'); ylabel('y');
colorbar;
axis equal tight; 

% Subplot 3: 2D plot of coefficient c
% subplot(1, 4, 3);
% imagesc(x, y, C_mat); 
% title('Coefficient c(x,y)');
% xlabel('x'); ylabel('y');
% colorbar;
% axis equal tight xy; 
% 
% % Subplot 4: 2D plot of coefficient d
% subplot(1, 4, 4);
% imagesc(x, y, D_mat);
% title('Coefficient d(x,y)');
% xlabel('x'); ylabel('y');
% colorbar;
% axis equal tight xy;


end 
