function plotLocalizationErrors(zValues, locs3D)
    plotXYError(zValues, locs3D);
end

function plotXYError(zValues, locs3D)
    % Extract positions
    xPositions = locs3D(:, 1);
    yPositions = locs3D(:, 2);

    % Calculate total errors in micrometers
    totalErrorX = abs(xPositions - 0) * 1000; % Convert to nanometers
    totalErrorY = abs(yPositions - 0) * 1000; % Convert to nanometers

    % Plotting
    figure(4);

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