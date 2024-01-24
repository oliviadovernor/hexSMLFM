function plotLocalizationErrors(zValues, locs3D)
    plotXYError(zValues, locs3D);
    plotZError(zValues, locs3D);
end

function plotXYError(zValues, locs3D)
    % Extract positions
    xPositions = locs3D(:, 1);
    yPositions = locs3D(:, 2);

    % Calculate total errors in micrometers
    totalErrorX = abs(xPositions - 0) * 1000; % Convert to nanometers
    totalErrorY = abs(yPositions - 0) * 1000; % Convert to nanometers

    % Plotting
    figure;

    % Plotting Lines for Simulated Data
    plot(zValues, totalErrorX, 'o-', 'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', 'X Error', 'Color', [1, 0, 0]);
    hold on;
    plot(zValues, totalErrorY, 's-', 'LineWidth', 1, 'MarkerSize', 6, 'DisplayName', 'Y Error', 'Color', [0, 1, 0]);
    hold off;

    % Customize Plot
    xlabel('Axial Range (um)');
    ylabel('Lateral error (nm)');
    title('Accuracy of localisation: XY axes');
    legend('Location', 'Best');
    grid on;
end

function plotZError(zValues, locs3D)
    % Extract Z positions
    zPositions = locs3D(:, 3) - locs3D(60, 3);

    % Calculate total error in micrometers
    totalErrorZ = abs(zValues(:) - zPositions(:)) * 1000; % Convert to nanometers

    % Plotting
    figure;

    % Plotting Line for Simulated Data
    plot(zValues, totalErrorZ, '^-', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Z Error', 'Color', [0, 0, 1]);
    hold on;  % Allow multiple plots on the same figure

    % Find coefficients for shape-preserving interpolant
    fitResults = pchip(zValues, totalErrorZ);

    % Evaluate fit
    yplot = ppval(fitResults, zValues);

    % Plot the shape-preserving fit in blue
    plot(zValues, yplot, 'DisplayName', 'Shape-Preserving Fit', 'LineWidth', 2, 'Color', [1, 0, 0]);

    % Display the equation on the plot (slightly more right)
    xPosition = min(zValues) + 0.05 * (max(zValues) - min(zValues)); % Adjust the factor as needed
    yPosition = min(yplot);
    equation = sprintf('Shape-Preserving Fit');
    text(xPosition, yPosition, equation, 'FontSize', 10, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');

    % Customize Plot for Corrected Data
    xlabel('Axial Range (um)');
    ylabel('Corrected Axial error (nm)');
    title('Corrected Z-axis Localization Accuracy');
    legend('show');
    hold off;
end




