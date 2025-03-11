function sobolev_gradient_newton_method(nx, ny, Lx, Ly, A, sigma)
    % Create a regular grid of Morley elements
    [nodes, elements] = create_morley_grid(nx, ny, Lx, Ly);
    % Number of wavelengths in each direction
    nx_waves = 3;
    ny_waves = 2;

    % Wave numbers
    kx = 2 * pi * nx_waves / Lx;
    ky = 2 * pi * ny_waves / Ly;
  
    % Generate the initial condition for the stream function psi
    psi = generate_sinusoidal_initial_condition_with_nodes(nodes, A, kx, ky);  
    
    % Initial guess for the solution (flatten the psi matrix to a vector)
    u = psi;
    
    % Plot the initial condition
    figure;
    trisurf(elements, nodes(:,1), nodes(:,2), psi, 'EdgeColor', 'none');
    colorbar;
    title('Initial Condition for Stream Function \psi (Sinusoidal)');
    xlabel('x');
    ylabel('y');
    zlabel('psi');
    
    % Newton method parameters
    max_iter = 20;
    tol = 1e-6;
    max_grad_evals = 1000; % Maximum number of gradient evaluations
    alpha = 0.5; % Backtracking line search parameter
    beta = 0.8; % Backtracking line search parameter
    reg_param = 1e-4; % Regularization parameter
    
    % Perform the Newton method
    grad_evals = 0;
    coefficients = cell(size(elements, 1), 1); % Store coefficients for each element
    for elem = 1:size(elements, 1)
        element_nodes = elements(elem, :);
        element_coords = nodes(element_nodes, :);
        ue = u(element_nodes);
        coefficients{elem} = calculate_coefficients(element_coords, ue);
    end
    
    for iter = 1:max_iter
        % Compute the gradient and Hessian of the energy functional
        [grad, hess] = compute_gradient_hessian(u, nodes, elements, reg_param, coefficients);
        grad_evals = grad_evals + 1;
        
        % Check for convergence
        grad_norm = norm(grad);
        fprintf('Iteration %d: Gradient norm = %e\n', iter, grad_norm);
        if grad_norm < tol
            fprintf('Converged in %d iterations.\n', iter);
            break;
        end
        
        % Check if the maximum number of gradient evaluations is reached
        if grad_evals >= max_grad_evals
            fprintf('Maximum number of gradient evaluations reached. Stopping.\n');
            break;
        end
        
        % Compute the search direction
        delta_u = -hess \ grad;
        
        % Backtracking line search
        t = 1;
        while true
            u_new = u + t * delta_u;
            [grad_new, ~] = compute_gradient_hessian(u_new, nodes, elements, reg_param, coefficients);
            if norm(grad_new) <= (1 - alpha * t) * grad_norm
                break;
            end
            t = beta * t;
        end
        
        % Update the solution
        u = u_new;
        
        % Update coefficients with the new solution
        for elem = 1:size(elements, 1)
            element_nodes = elements(elem, :);
            element_coords = nodes(element_nodes, :);
            ue = u(element_nodes);
            coefficients{elem} = calculate_coefficients(element_coords, ue);
        end
    end
    
    % Reshape the solution to its original 2D form
    psi_solution = u;
    
    % Plot the gradient on the nodes
    figure;
    trisurf(elements, nodes(:,1), nodes(:,2), grad, 'EdgeColor', 'none');
    colorbar;
    title('Gradient on the Nodes');
    xlabel('x');
    ylabel('y');
    zlabel('Gradient');
    
    % Compute and plot the gradient on the edge elements
    edge_gradients = zeros(size(elements, 1), 1);
    for elem = 1:size(elements, 1)
        element_nodes = elements(elem, :);
        element_coords = nodes(element_nodes, :);
        ue = u(element_nodes);
        edge_gradients(elem) = norm(compute_edge_gradient(ue, element_coords, coefficients{elem}));
    end
    figure;
    trisurf(elements, nodes(:,1), nodes(:,2), edge_gradients, 'EdgeColor', 'none');
    colorbar;
    title('Gradient on the Edge Elements');
    xlabel('x');
    ylabel('y');
    zlabel('Gradient');
    
    % Plot the final solution
    figure;
    trisurf(elements, nodes(:,1), nodes(:,2), psi_solution, 'EdgeColor', 'none');
    colorbar;
    title('Final Solution of the PDE using Newton Method');
    xlabel('x');
    ylabel('y');
    zlabel('psi');
end

function [nodes, elements] = create_morley_grid(nx, ny, Lx, Ly)
    % Generate a regular grid of Morley elements for a domain of size Lx by Ly
    % with nx elements in the x-direction and ny elements in the y-direction.
    
    % Create grid points
    x = linspace(0, Lx, nx+1);
    y = linspace(0, Ly, ny+1);
    [X, Y] = meshgrid(x, y);
    nodes = [X(:), Y(:)];
    
    % Create elements with connectivity
    elements = [];
    for j = 1:ny
        for i = 1:nx
            % Indices of the four corners of the current rectangle
            ll = (j-1)*(nx+1) + i; % lower-left
            lr = ll + 1;           % lower-right
            ul = j*(nx+1) + i;     % upper-left
            ur = ul + 1;           % upper-right
            
            % Create two triangles for each rectangle
            elements = [elements; ll, lr, ur; ll, ur, ul]; % two triangles forming the rectangle
        end
    end
    
    % Remove duplicate nodes
    [nodes, ~, idx] = unique(nodes, 'rows');
    
    % Update the elements array with the new node indices
    elements = idx(elements);
end

function [grad, hess] = compute_gradient_hessian(u, nodes, elements, reg_param, coefficients)
    % Compute the gradient and Hessian of the energy functional for
    % stationary quasi-geostrophic flow in 2D

    num_nodes = size(nodes, 1);
    grad = zeros(num_nodes, 1);
    hess = sparse(num_nodes, num_nodes);
    
    % Loop over elements to assemble the gradient and Hessian
    for elem = 1:size(elements, 1)
        element_nodes = elements(elem, :);
        
        % Check if element_nodes exceeds the number of nodes
        if any(element_nodes > num_nodes)
            fprintf('Element nodes: %s\n', mat2str(element_nodes));
            error('Index exceeds the number of array elements. Index must not exceed %d.', num_nodes);
        end
        
        element_coords = nodes(element_nodes, :);
        ue = u(element_nodes);
        
        % Compute the element gradient and Hessian for nodes
        [elem_grad, elem_hess] = element_gradient_hessian_nodes(ue, element_coords);
        
        % Assemble into global gradient and Hessian
        grad(element_nodes) = grad(element_nodes) + elem_grad;
        hess(element_nodes, element_nodes) = hess(element_nodes, element_nodes) + elem_hess;
        
        % Compute the element gradient and Hessian for edges
        [elem_grad, elem_hess] = element_gradient_hessian_edges(ue, element_coords, coefficients{elem});
        
        % Assemble into global gradient and Hessian
        grad(element_nodes) = grad(element_nodes) + elem_grad;
        hess(element_nodes, element_nodes) = hess(element_nodes, element_nodes) + elem_hess;
    end
    
    % Add regularization term to the Hessian
    hess = hess + reg_param * speye(num_nodes);
end

function coefficients = calculate_coefficients(element_coords, ue)
    % Calculate the coefficients a, b, c, d, e, f for the Morley basis functions
    % element_coords: coordinates of the element nodes
    % ue: values of the degrees of freedom at the element nodes

    % Solve for the coefficients using the node data
    % Assume the basis function is of the form: phi = ax^2 + bxy + cy^2 + dx + ey + f
    % We need to solve the system of equations for the coefficients
    
    % Construct the matrix for the system of equations
    A = [element_coords(:,1).^2, element_coords(:,1).*element_coords(:,2), element_coords(:,2).^2, element_coords(:,1), element_coords(:,2), ones(size(element_coords, 1), 1)];
    
    % Solve for the coefficients
    coefficients = A \ ue;
end

function [elem_grad, elem_hess] = element_gradient_hessian_nodes(ue, element_coords)
    % Compute the gradient and Hessian for a single element based on nodes
    % ue: values of the degrees of freedom at the element nodes
    % element_coords: coordinates of the element nodes

    % Check dimensions for matrix operations
    if size(element_coords, 2) ~= 2
        error('Dimension mismatch: element_coords should have 2 columns.');
    end
    if length(ue) ~= size(element_coords, 1)
        error('Dimension mismatch: length of ue should match number of rows in element_coords.');
    end
    
    % Compute the basis functions and their gradients (linear for nodes)
    phi = ones(size(element_coords, 1), 1); % Placeholder linear basis function
    grad_phi = [1, 0, 0; 0, 1, 0]; % Gradient of linear basis function for nodes with correct dimensions (2x3)
    
    % Compute the Jacobian of the element transformation
    J = [element_coords(2, :) - element_coords(1, :); element_coords(3, :) - element_coords(1, :)];
    if size(J, 1) ~= size(J, 2)
        error('Jacobian matrix must be square.');
    end
    detJ = det(J);
    invJ = inv(J);
    
    % Compute the gradient of the stream function
    grad_psi = invJ * grad_phi * ue;
    
    % Compute the potential vorticity q
    q = sum(grad_psi .^ 2, 1);
    
    % Initialize the element gradient and Hessian
    elem_grad = zeros(3, 1);
    elem_hess = zeros(3, 3);
    
    % Compute the element gradient and Hessian
    for i = 1:3
        elem_grad(i) = detJ * (grad_psi' * grad_phi(:, i)) + phi(i) * q;
        for j = 1:3
            elem_hess(i, j) = detJ * (grad_phi(:, i)' * grad_phi(:, j));
        end
    end
end

function [elem_grad, elem_hess] = element_gradient_hessian_edges(ue, element_coords, coefficients)
    % Compute the gradient and Hessian for a single element based on edges
    % ue: values of the degrees of freedom at the element nodes
    % element_coords: coordinates of the element nodes
    % coefficients: coefficients a, b, c, d, e, f for the Morley basis functions

    % Check dimensions for matrix operations
    if size(element_coords, 2) ~= 2
        error('Dimension mismatch: element_coords should have 2 columns.');
    end
    if length(ue) ~= size(element_coords, 1)
        error('Dimension mismatch: length of ue should match number of rows in element_coords.');
    end
    
    % Extract the coefficients
    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);
    d = coefficients(4);
    e = coefficients(5);
    f = coefficients(6);
    
    % Compute the basis functions and their gradients
    phi = a * element_coords(:,1).^2 + b * element_coords(:,1) .* element_coords(:,2) + c * element_coords(:,2).^2 + d * element_coords(:,1) + e * element_coords(:,2) + f;
    grad_phi = [2 * a * element_coords(:,1) + b * element_coords(:,2) + d, b * element_coords(:,1) + 2 * c * element_coords(:,2) + e]'; % Transposed for correct dimensions
    
    % Compute the Jacobian of the element transformation
    J = [element_coords(2, :) - element_coords(1, :); element_coords(3, :) - element_coords(1, :)];
    if size(J, 1) ~= size(J, 2)
        error('Jacobian matrix must be square.');
    end
    detJ = det(J);
    invJ = inv(J);
    
    % Compute the gradient of the stream function
    grad_psi = invJ * (grad_phi * ue);
    
    % Compute the potential vorticity q
    q = sum(grad_psi .^ 2, 1);
    
    % Initialize the element gradient and Hessian
    elem_grad = zeros(3, 1);
    elem_hess = zeros(3, 3);
    
    % Compute the element gradient and Hessian
    for i = 1:3
        elem_grad(i) = detJ * (grad_psi' * grad_phi(:, i)) + phi(i) * q;
        for j = 1:3
            elem_hess(i, j) = detJ * (grad_phi(:, i)' * grad_phi(:, j));
        end
    end
end

function edge_gradient = compute_edge_gradient(ue, element_coords, coefficients)
    % Compute the gradient for edge elements
    % ue: values of the degrees of freedom at the element nodes
    % element_coords: coordinates of the element nodes
    % coefficients: coefficients a, b, c, d, e, f for the Morley basis functions

    % Check dimensions for matrix operations
    if size(element_coords, 2) ~= 2
        error('Dimension mismatch: element_coords should have 2 columns.');
    end
    if length(ue) ~= size(element_coords, 1)
        error('Dimension mismatch: length of ue should match number of rows in element_coords.');
    end
    
    % Extract the coefficients
    a = coefficients(1);
    b = coefficients(2);
    c = coefficients(3);
    d = coefficients(4);
    e = coefficients(5);
    f = coefficients(6);
    
    % Compute the basis functions and their gradients
    grad_phi = [2 * a * element_coords(:,1) + b * element_coords(:,2) + d, b * element_coords(:,1) + 2 * c * element_coords(:,2) + e]'; % Transposed for correct dimensions
    
    % Compute the Jacobian of the element transformation
    J = [element_coords(2, :) - element_coords(1, :); element_coords(3, :) - element_coords(1, :)];
    invJ = inv(J);
    
    % Compute the gradient of the stream function
    edge_gradient = invJ * (grad_phi * ue);
end

function psi = generate_sinusoidal_initial_condition_with_nodes(nodes, A, kx, ky)
    % Generate a sinusoidal initial condition for the stream function psi
    % nodes: Nx2 array of (x, y) coordinates
    % A: amplitude of the sinusoidal wave
    % kx, ky: wave numbers in x and y directions

    % Extract x and y from nodes
    x = nodes(:, 1);
    y = nodes(:, 2);
    
    % Compute the sinusoidal initial condition
    psi = A * sin(kx * x) .* sin(ky * y);
    
    % Plot the initial condition for visualization
    figure;
    scatter(x, y, 20, psi, 'filled');
    colorbar;
    title('Initial Condition for Stream Function \psi (Sinusoidal)');
    xlabel('x');
    ylabel('y');
end