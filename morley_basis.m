function [nodes, elements] = morley_basis()
    % Define the vertices of the reference triangle
    vertices = [0, 0; 1, 0; 0, 1];
    
    % Define the midpoints of the edges
    midpoints = (vertices([1 1 2], :) + vertices([2 3 3], :)) / 2;
    
    % Define the degrees of freedom (dofs) for the Morley element
    % 1-3: vertex values
    % 4-6: edge midpoints
    nodes = [vertices; midpoints];
    
    % Define the connectivity of the Morley element
    elements = [1 2 3 4 5 6];
    
    % Plot the reference triangle with nodes
    figure;
    hold on;
    fill(vertices(:,1), vertices(:,2), 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'k');
    plot(nodes(:,1), nodes(:,2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    text(nodes(:,1) + 0.02, nodes(:,2), arrayfun(@(n) sprintf('N%d', n), 1:6, 'UniformOutput', false));
    xlim([-0.1 1.1]);
    ylim([-0.1 1.1]);
    axis equal;
    title('Morley Element Nodes');
    hold off;
end