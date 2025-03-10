function create_morley_grid(nx, ny)
    % Create a regular grid of Morley elements
    % nx: number of elements in the x-direction
    % ny: number of elements in the y-direction

    % Define the vertices of the reference triangle
    vertices = [0, 0; 1, 0; 0, 1];
    
    % Define the midpoints of the edges of the reference triangle
    midpoints = (vertices([1 1 2], :) + vertices([2 3 3], :)) / 2;
    
    % Define the degrees of freedom (dofs) for the Morley element
    % 1-3: vertex values
    % 4-6: edge midpoints
    ref_nodes = [vertices; midpoints];
    
    % Initialize arrays for nodes and elements
    nodes = [];
    elements = [];
    
    % Generate the grid of elements
    for j = 0:ny
        for i = 0:nx
            % Calculate the coordinates of the lower-left vertex of the current element
            x0 = i;
            y0 = j;
            
            % Define the nodes of the current element
            element_nodes = ref_nodes + [x0, y0];
            
            % Append the nodes to the global nodes array
            nodes = [nodes; element_nodes];
            
            % Define the connectivity of the current element
            num_nodes = size(nodes, 1);
            element = num_nodes - 5:num_nodes;
            
            % Append the element to the global elements array
            elements = [elements; element];
        end
    end
    
    % Remove duplicate nodes
    [nodes, ~, idx] = unique(nodes, 'rows');
    
    % Update the elements array with the new node indices
    elements = idx(elements);
    
    % Plot the grid of Morley elements
    figure;
    hold on;
    for k = 1:size(elements, 1)
        fill(nodes(elements(k, :), 1), nodes(elements(k, :), 2), 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    end
    plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    xlim([-0.1 nx+1.1]);
    ylim([-0.1 ny+1.1]);
    axis equal;
    title('Regular Grid of Morley Elements');
    hold off;
end