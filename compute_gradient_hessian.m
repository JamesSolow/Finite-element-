function [grad, hess] = compute_gradient_hessian(u, nodes, elements)
    % Compute the gradient and Hessian of the energy functional for
    % stationary quasi-geostrophic flow in 2D

    num_nodes = size(nodes, 1);
    grad = zeros(num_nodes, 1);
    hess = sparse(num_nodes, num_nodes);
    
    % Loop over elements to assemble the gradient and Hessian
    for elem = 1:size(elements, 1)
        element_nodes = elements(elem, :);
        element_coords = nodes(element_nodes, :);
        ue = u(element_nodes);
        
        % Compute the element gradient and Hessian
        [elem_grad, elem_hess] = element_gradient_hessian(ue, element_coords);
        
        % Assemble into global gradient and Hessian
        grad(element_nodes) = grad(element_nodes) + elem_grad;
        hess(element_nodes, element_nodes) = hess(element_nodes, element_nodes) + elem_hess;
    end
end

function [elem_grad, elem_hess] = element_gradient_hessian(ue, element_coords)
    % Compute the gradient and Hessian for a single element
    % ue: values of the degrees of freedom at the element nodes
    % element_coords: coordinates of the element nodes

    % Define the basis functions and their gradients
    % (This will depend on the specific implementation of the Morley element)
    
    % Placeholder values for basis functions and their gradients
    % These need to be replaced with the actual basis functions
    phi = [1/3; 1/3; 1/3; 0; 0; 0]; % Example basis functions
    grad_phi = [1, 0; 0, 1; -1, -1; 0, 0; 0, 0; 0, 0]; % Example gradients of basis functions
    
    % Compute the Jacobian of the element transformation
    J = [element_coords(2, :) - element_coords(1, :); element_coords(3, :) - element_coords(1, :)];
    detJ = det(J);
    invJ = inv(J);
    
    % Compute the gradient of the stream function
    grad_psi = grad_phi * invJ * ue;
    
    % Compute the potential vorticity q
    q = sum(grad_psi .^ 2, 2);

    % Initialize the element gradient and Hessian
    elem_grad = zeros(6, 1);
    elem_hess = zeros(6, 6);
    
    % Compute the element gradient and Hessian
    for i = 1:6
        elem_grad(i) = detJ * (grad_psi' * grad_phi(:, i) + phi(i) * q);
        for j = 1:6
            elem_hess(i, j) = detJ * (grad_phi(:, i)' * grad_phi(:, j));
        end
    end
end