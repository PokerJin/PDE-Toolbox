% -------------------------------------------------------------------------
% MATLAB Code for Solving 2D Multiscale Elliptic PDE
% -div(c * grad(u)) = f
% using a Nine-Point Finite Difference Stencil
% CORRECTED BOUNDARY CONDITION HANDLING
% -------------------------------------------------------------------------

function solve_mulell_dirnine(xy_range,Nxy,equ_para,bc_type,bc_para)

% --- User-Defined Parameters ---

% Domain Boundaries (Square Domain Required: x_right-x_left = y_right-y_left)
x_left = xy_range(1);
x_right = xy_range(2);
y_left = xy_range(3);
y_right = xy_range(4);

% Grid Discretization
Nx = Nxy(1); % Number of points in x-direction 
Ny = Nxy(2); % Number of points in y-direction
dx = (x_right - x_left) / (Nx - 1);
dy = (y_right - y_left) / (Ny - 1);
if abs(dx - dy) > 1e-12 % Check if dx and dy are effectively equal
    warning('dx and dy are not equal. The 9-point stencil assumes dx=dy=h.');
end
h = dx; % Assuming dx = dy = h for the stencil

% Multiscale Coefficient c(x,y) Parameters
epsilon = equ_para.epsilon*(x_right-x_left);      % Microscale period (size of microcells)
r_circle = equ_para.len;    % Radius of circles within microcells (must be < epsilon/2)
c_in = equ_para.c_in;           % Value of c inside circles
c_out = equ_para.c_ma;        % Value of c outside circles

% Source Term f(x,y)
f_func = equ_para.f;

% Boundary Conditions (BC)
% --- Specify Type for Each Boundary ---
% 'D' for Dirichlet (h*u = r)
% 'N' for Neumann (n*c*grad(u) + q*u = g)
if bc_type(1)==1
    bc_type_left = 'D'; 
else 
    bc_type_left = 'N'; 
end
if bc_type(2)==1
    bc_type_right = 'D'; 
else
    bc_type_right = 'N'; 
end
% Right boundary (x = x_right)
if bc_type(3)==1
    bc_type_bottom = 'D';
else
    bc_type_bottom = 'N';
end
if bc_type(4)==1
    bc_type_top = 'D'; 
else
    bc_type_top = 'N';
end


% --- Define BC Parameters (can be constants or functions @(x,y)) ---

% Dirichlet Parameters (h*u = r)
h_left   = bc_para.h_left;      r_left   = bc_para.r_left;
h_right  = bc_para.h_right;      r_right  = bc_para.r_right;
h_bottom = bc_para.h_bottom;      r_bottom = bc_para.r_bottom;
h_top    = bc_para.h_top;      r_top    = bc_para.r_top;

% Neumann Parameters (n*c*grad(u) + q*u = g)

q_left   = bc_para.q_left;      g_left   =bc_para.g_left;
q_right  = bc_para.q_right;      g_right  = bc_para.g_left;
q_bottom = bc_para.q_bottom;      g_bottom = bc_para.g_bottom;
q_top    = bc_para.q_top;      g_top    = bc_para.g_top;


% --- Grid Generation ---
[X, Y] = meshgrid(linspace(x_left, x_right, Nx), linspace(y_left, y_right, Ny));
N_total = Nx * Ny;

% --- Calculate Coefficient Field c(x,y) ---
C = zeros(Ny, Nx);
for j = 1:Ny 
    for i = 1:Nx 
        x = X(j, i);
        y = Y(j, i);

        % Find center of the microcell containing (x,y)
        cell_idx_x = floor((x - x_left) / epsilon);
        cell_idx_y = floor((y - y_left) / epsilon);
        center_x = x_left + (cell_idx_x + 0.5) * epsilon;
        center_y = y_left + (cell_idx_y + 0.5) * epsilon;

        % Calculate distance from the center of the microcell
        dist = sqrt((x - center_x)^2 + (y - center_y)^2);

        % Check if inside or on the boundary of the circle
        if dist <= r_circle
            C(j, i) = c_in;
        else
            C(j, i) = c_out;
        end
    end
end

% --- Assemble Stiffness Matrix A and RHS Vector b ---

% Use sparse matrix for efficiency
A = spalloc(N_total, N_total, 9 * N_total); % Estimate non-zeros
b = zeros(N_total, 1);

% Helper function for 1D indexing
k = @(i, j) i + (j - 1) * Nx; 

% --- Step 1: Fill A and b assuming all points are interior ---
for j = 1:Ny
    for i = 1:Nx
        current_k = k(i, j);

        % Get coefficients at stencil points (handle boundaries by clamping indices)
        c_C = C(j,i);
        c_N = C(min(j+1, Ny), i);
        c_S = C(max(j-1, 1), i);
        c_E = C(j, min(i+1, Nx));
        c_W = C(j, max(i-1, 1));
        c_NE = C(min(j+1, Ny), min(i+1, Nx));
        c_NW = C(min(j+1, Ny), max(i-1, 1));
        c_SE = C(max(j-1, 1), min(i+1, Nx));
        c_SW = C(max(j-1, 1), max(i-1, 1));

        % Calculate effective coefficients between points
        ce = (c_E + c_C)/2; cw = (c_W + c_C)/2;
        cn = (c_N + c_C)/2; cs = (c_S + c_C)/2;
        cne = (c_NE + c_N + c_E + c_C)/4;
        cnw = (c_NW + c_N + c_W + c_C)/4;
        cse = (c_SE + c_S + c_E + c_C)/4;
        csw = (c_SW + c_S + c_W + c_C)/4;

        h2 = h*h;

        % 5-point stencil part coefficients
        A_C_5pt = (ce+cw+cn+cs)/h2;
        A_E_5pt = -ce/h2; A_W_5pt = -cw/h2;
        A_N_5pt = -cn/h2; A_S_5pt = -cs/h2;

        % Diagonal stencil part coefficients
        A_C_diag = (cne+cnw+cse+csw)/(2*h2);
        A_NE_diag = -cne/(2*h2); A_NW_diag = -cnw/(2*h2);
        A_SE_diag = -cse/(2*h2); A_SW_diag = -csw/(2*h2);

        % Combine for 9-point stencil (weighted average)
        A_C  = (2/3)*A_C_5pt + (1/3)*A_C_diag;
        A_E  = (2/3)*A_E_5pt;
        A_W  = (2/3)*A_W_5pt;
        A_N  = (2/3)*A_N_5pt;
        A_S  = (2/3)*A_S_5pt;
        A_NE = (1/3)*A_NE_diag;
        A_NW = (1/3)*A_NW_diag;
        A_SE = (1/3)*A_SE_diag;
        A_SW = (1/3)*A_SW_diag;

        % Assign stencil coefficients to matrix A
        % These will be overwritten later if the point is on a boundary
        A(current_k, current_k)   = A_C;
        if i+1 <= Nx, A(current_k, k(i+1, j))   = A_E; end
        if i-1 >= 1,  A(current_k, k(i-1, j))   = A_W; end
        if j+1 <= Ny, A(current_k, k(i, j+1))   = A_N; end
        if j-1 >= 1,  A(current_k, k(i, j-1))   = A_S; end
        if i+1 <= Nx && j+1 <= Ny, A(current_k, k(i+1, j+1)) = A_NE; end
        if i-1 >= 1  && j+1 <= Ny, A(current_k, k(i-1, j+1)) = A_NW; end
        if i+1 <= Nx && j-1 >= 1,  A(current_k, k(i+1, j-1)) = A_SE; end
        if i-1 >= 1  && j-1 >= 1,  A(current_k, k(i-1, j-1)) = A_SW; end

        % Assign source term to RHS vector b
        b(current_k) = f_func(X(j, i), Y(j, i));
    end % end i loop
end % end j loop

% --- Step 2: Apply Boundary Conditions (Overwrite rows for BCs) ---

% Bottom Boundary (j = 1)
j = 1;
for i = 1:Nx
    current_k = k(i, j);
    x = X(j, i); y = Y(j, i);
    if bc_type_bottom == 'D'
        A(current_k, :) = 0; % Clear the row
        A(current_k, current_k) = h_bottom(x, y);
        b(current_k) = r_bottom(x, y);
    elseif bc_type_bottom == 'N'
        % Using original Neumann logic (may need review for 9-point stencil)
        A(current_k, :) = 0; % Clear row first
        c_val = C(j,i);
        q_val = q_bottom(x,y);
        g_val = g_bottom(x,y);
        A(current_k, current_k) = c_val/dy + q_val; 
        if j+1 <= Ny, A(current_k, k(i, j+1)) = -c_val/dy; end 
        b(current_k) = g_val;
        if i>1, A(current_k, k(i-1,j))=0; end
        if i<Nx, A(current_k, k(i+1,j))=0; end
        if i>1 && j+1<=Ny, A(current_k, k(i-1,j+1))=0; end 
        if i<Nx && j+1<=Ny, A(current_k, k(i+1,j+1))=0; end
    end
end

% Top Boundary (j = Ny)
j = Ny;
for i = 1:Nx
    current_k = k(i, j);
    x = X(j, i); y = Y(j, i);
    if bc_type_top == 'D'
        A(current_k, :) = 0;
        A(current_k, current_k) = h_top(x, y);
        b(current_k) = r_top(x, y);
    elseif bc_type_top == 'N'
        % Using original Neumann logic
        A(current_k, :) = 0;
        c_val = C(j,i);
        q_val = q_top(x,y);
        g_val = g_top(x,y);
        A(current_k, current_k)   = c_val/dy + q_val; % Normal points outwards (+y), grad approx involves u(i,j-1)
        if j-1 >= 1, A(current_k, k(i, j-1)) = -c_val/dy; end % Coeff for u(i,j-1)
        b(current_k) = g_val;
        % Ensure other terms are zero
        if i>1, A(current_k, k(i-1,j))=0; end
        if i<Nx, A(current_k, k(i+1,j))=0; end
        if i>1 && j-1>=1, A(current_k, k(i-1,j-1))=0; end % SW relative to (i,j-1)
        if i<Nx && j-1>=1, A(current_k, k(i+1,j-1))=0; end % SE relative to (i,j-1)
    end
end

% Left Boundary 
i = 1;
% Avoid double-setting corners already handled by Top/Bottom Dirichlet
for j = (1 + (bc_type_bottom=='D')) : (Ny - (bc_type_top=='D'))
    current_k = k(i, j);
    x = X(j, i); y = Y(j, i);
    if bc_type_left == 'D'
        A(current_k, :) = 0;
        A(current_k, current_k) = h_left(x, y);
        b(current_k) = r_left(x, y);
    elseif bc_type_left == 'N'
        A(current_k, :) = 0;
        c_val = C(j,i);
        q_val = q_left(x,y);
        g_val = g_left(x,y);
        A(current_k, current_k)   = c_val/dx + q_val;
        if i+1 <= Nx, A(current_k, k(i+1, j)) = -c_val/dx; end 
        b(current_k) = g_val;
         if j>1, A(current_k, k(i,j-1))=0; end
         if j<Ny, A(current_k, k(i,j+1))=0; end
         if j>1 && i+1<=Nx, A(current_k, k(i+1,j-1))=0; end 
         if j<Ny && i+1<=Nx, A(current_k, k(i+1,j+1))=0; end
    end
end

% Right Boundary
i = Nx;
for j = (1 + (bc_type_bottom=='D')) : (Ny - (bc_type_top=='D'))
    current_k = k(i, j);
    x = X(j, i); y = Y(j, i);
     if bc_type_right == 'D'
        A(current_k, :) = 0;
        A(current_k, current_k) = h_right(x, y);
        b(current_k) = r_right(x, y);
     elseif bc_type_right == 'N'
         A(current_k, :) = 0;
         c_val = C(j,i);
         q_val = q_right(x,y);
         g_val = g_right(x,y);
         A(current_k, current_k)   = c_val/dx + q_val; 
         if i-1 >= 1, A(current_k, k(i-1, j)) = -c_val/dx; end 
         b(current_k) = g_val;
         % Ensure other terms are zero
         if j>1, A(current_k, k(i,j-1))=0; end
         if j<Ny, A(current_k, k(i,j+1))=0; end
         if j>1 && i-1>=1, A(current_k, k(i-1,j-1))=0; end 
         if j<Ny && i-1>=1, A(current_k, k(i-1,j+1))=0; end 
     end
end


% --- Solve the Linear System ---
U_vec = A \ b; 

% --- Reshape Solution and Plot ---
U = reshape(U_vec, Ny, Nx); 

% Create a single figure window for the required plots
figure; % [left, bottom, width, height]

% 1. 3D Surface Plot
subplot(1, 2, 1);
surf(X, Y, U);
% colormap(gca, jet);
% colorbar;
xlabel('x');
ylabel('y');
zlabel('u');
title('参考解（三维图）');
axis tight;
shading interp;
view(3);

% 2. 2D Contour Plot
subplot(1, 2, 2);
contourf(X, Y, U, 20); 
colorbar;
xlabel('x');
ylabel('y');
title('参考解（二维图）');
axis equal tight;

% % 3. Plot Coefficient c(x,y)
% subplot(1, 3, 3);
% imagesc(linspace(x_left, x_right, Nx), linspace(y_left, y_right, Ny), C);
% colorbar;
% xlabel('x');
% ylabel('y');
% title('Coefficient c(x,y)');
% axis equal tight xy; 
end