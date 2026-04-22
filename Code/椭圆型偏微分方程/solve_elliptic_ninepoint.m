% Simplified elliptic PDE solver without analytical solution comparison
% Solves 2D elliptic PDE: -div(c*grad(u))+a*u=f using nine-point difference scheme
function solve_elliptic_ninepoint(xy_range,Nxy,equ_para,bc_type,bc_para)
    % Define domain parameters
    x_left = xy_range(1);
    x_right =xy_range(2);
    y_bottom = xy_range(3);
    y_top = xy_range(4);
    nx = Nxy(1);
    ny = Nxy(2);
    dx = (x_right - x_left) / (nx - 1);
    dy = (y_top - y_bottom) / (ny - 1);
    
    tol=1e-8;
    maxit=1000;
    
    x = linspace(x_left, x_right, nx);
    y = linspace(y_bottom, y_top, ny);
    [X, Y] = meshgrid(x, y);
    
    c_func =equ_para.c;  % Diffusion coefficient
    a_func =equ_para.a; % Reaction term
    f_func =equ_para.f; % Source term
    
    % bc_type = [1, 1, 1, 1];
    
    % For Dirichlet: h*u = r
    h_funcs = {
        bc_para.h_left,
        bc_para.h_right,
        bc_para.h_bottom,
        bc_para.h_top
        };

    r_funcs = {
        bc_para.r_left,
        bc_para.r_right,
        bc_para.r_bottom,
        bc_para.r_top
        };

    % For Neumann: n*c*grad(u) + q*u = g
    q_funcs = {
        bc_para.q_left,
        bc_para.q_right,
        bc_para.q_bottom,
        bc_para.q_top
        };

    g_funcs = {
        bc_para.g_left,
        bc_para.g_right,
        bc_para.g_bottom,
        bc_para.g_top
        };
    
    % Solve the PDE
    u_numerical = solve_elliptic_9point(X, Y, dx, dy, c_func, a_func, f_func, ...
        h_funcs, r_funcs, q_funcs, g_funcs, bc_type,tol,maxit);
    
    % Visualize the numerical solution
    figure;
    subplot(1,2,1);
    surf(X, Y, u_numerical);
    title('数值解（三维图）');
    xlabel('x'); ylabel('y'); zlabel('u');
    shading interp;
    view(3);
    axis tight;
    colorbar;
    
    subplot(1,2,2);
    contourf(X, Y, u_numerical, 20);
    title('数值解（二维图）');
    xlabel('x'); ylabel('y');
    axis tight equal;
    colorbar;
    
end

%% Solver function
function u = solve_elliptic_9point(X, Y, dx, dy, c_func, a_func, f_func, ...
                                 h_funcs, r_funcs, q_funcs, g_funcs, bc_type,tol,maxit)

    [ny, nx] = size(X);
    N = nx * ny;  

    A = sparse(N, N);
    b = zeros(N, 1);

    for j = 1:ny
        for i = 1:nx
            idx = (j-1)*nx + i;

            xi = X(j,i);
            yi = Y(j,i);

            on_left = (i == 1);
            on_right = (i == nx);
            on_bottom = (j == 1);
            on_top = (j == ny);
            
            if on_left || on_right || on_bottom || on_top
                if on_left
                    normal = [-1, 0]; 
                    bc = bc_type(1);
                    boundary_idx = 1; 
                elseif on_right
                    normal = [1, 0];  
                    bc = bc_type(2);
                    boundary_idx = 2;
                elseif on_bottom
                    normal = [0, -1]; 
                    bc = bc_type(3);
                    boundary_idx = 3;  
                else % on_top
                    normal = [0, 1];   
                    bc = bc_type(4);
                    boundary_idx = 4;  
                end
                
                [A, b] = apply_boundary_condition(A, b, idx, i, j, xi, yi, nx, ny, ...
                                             dx, dy, normal, bc, c_func, ...
                                             h_funcs{boundary_idx}, r_funcs{boundary_idx}, ...
                                             q_funcs{boundary_idx}, g_funcs{boundary_idx});
            else
                [A, b] = apply_nine_point_stencil(A, b, idx, i, j, xi, yi, nx, ny, ...
                                             dx, dy, c_func, a_func, f_func);
            end
        end
    end
    
    
    u_vec = A \ b;
    u = reshape(u_vec, nx, ny)';
end

function [A, b] = apply_nine_point_stencil(A, b, idx, i, j, xi, yi, nx, ny, ...
                                       dx, dy, c_func, a_func, f_func)

    % Evaluate coefficients at current point
    c = c_func(xi, yi);
    a = a_func(xi, yi);
    f = f_func(xi, yi);
    
    % Grid spacings
    h = dx;
    k = dy;
    
    % Index for center point (i,j)
    center_idx = idx;
    
    % Cardinal direction indices
    east_idx = center_idx + 1;   
    west_idx = center_idx - 1;   
    north_idx = center_idx + nx; 
    south_idx = center_idx - nx; 
    
    % Diagonal direction indices
    ne_idx = north_idx + 1; 
    nw_idx = north_idx - 1; 
    se_idx = south_idx + 1; 
    sw_idx = south_idx - 1; 
    
    if abs(h - k) < 1e-10  
        A(center_idx, center_idx) = c * (20/6) / h^2;
 
        A(center_idx, east_idx) = c * (-4/6) / h^2;  % East
        A(center_idx, west_idx) = c * (-4/6) / h^2;  % West
        A(center_idx, north_idx) = c * (-4/6) / h^2; % North
        A(center_idx, south_idx) = c * (-4/6) / h^2; % South

        A(center_idx, ne_idx) = c * (-1/6) / h^2;  % Northeast
        A(center_idx, nw_idx) = c * (-1/6) / h^2;  % Northwest
        A(center_idx, se_idx) = c * (-1/6) / h^2;  % Southeast
        A(center_idx, sw_idx) = c * (-1/6) / h^2;  % Southwest
    else

        A(center_idx, center_idx) = 2*c*(1/h^2 + 1/k^2);

        A(center_idx, east_idx) = -c/h^2;  % East
        A(center_idx, west_idx) = -c/h^2;  % West
        A(center_idx, north_idx) = -c/k^2; % North
        A(center_idx, south_idx) = -c/k^2; % South
    end

    A(center_idx, center_idx) = A(center_idx, center_idx) + a;

    b(center_idx) = f;
end

function [A, b] = apply_boundary_condition(A, b, idx, i, j, xi, yi, nx, ny, ...
                                       dx, dy, normal, bc_type, c_func, ...
                                       h_func, r_func, q_func, g_func)

    A(idx, :) = 0;

    c = c_func(xi, yi);
    
    if bc_type == 1  
        h = h_func(xi, yi);
        r = r_func(xi, yi);

        A(idx, idx) = h;

        b(idx) = r;
    else 
        q = q_func(xi, yi);
        g = g_func(xi, yi);
        nx_val = normal(1);
        ny_val = normal(2);
        if abs(nx_val) > 0  
            if nx_val < 0  
                A(idx, idx) = -c/dx + q;
                A(idx, idx+1) = c/dx;
            else  
                A(idx, idx) = c/dx + q;
                A(idx, idx-1) = -c/dx;
            end
        else  
            if ny_val < 0  
                A(idx, idx) = -c/dy + q;
                A(idx, idx+nx) = c/dy;
            else  
                A(idx, idx) = c/dy + q;
                A(idx, idx-nx) = -c/dy;
            end
        end
        b(idx) = g;
    end
end