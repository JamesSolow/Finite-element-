function psi = generate_eddy_initial_condition(nx, ny, Lx, Ly, A, x0, y0, sigma)
    % Generate a Gaussian initial condition for the stream function psi
    % nx, ny: number of grid points in x and y directions
    % Lx, Ly: domain size in x and y directions
    % A: amplitude of the eddy
    % x0, y0: center of the eddy
    % sigma: standard deviation of the Gaussian

    % Create a grid
    x = linspace(0, Lx, nx);
    y = linspace(0, Ly, ny);
    [X, Y] = meshgrid(x, y);
    
    % Compute the Gaussian initial condition
    psi = A * exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2));
    
    % Plot the initial condition for visualization
    figure;
    contourf(X, Y, psi, 20);
    colorbar;
    title('Initial Condition for Stream Function \psi');
    xlabel('x');
    ylabel('y');
end