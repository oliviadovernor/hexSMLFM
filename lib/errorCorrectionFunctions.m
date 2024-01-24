function correctedZPN = errorCorrectionFunctions(zValues, locs3D)

    % Extract Z positions
    zPositions = locs3D(:, 3) - locs3D(60, 3);

    % Calculate total error in micrometers
    totalErrorZ = abs(zValues(:) - zPositions(:)) * 1000; % Convert to nanometers
    ZErrorPN = (zValues(:) - zPositions(:)) * 1000; % Convert to nanometers
    disp(ZErrorPN);

    % Find coefficients for shape-preserving interpolant
    fitResults = pchip(zValues, totalErrorZ);
    fitResultsPN = pchip(zValues, ZErrorPN);
    
    % Evaluate piecewise polynomial
    fitArray(:, 1) = ppval(fitResults, zValues);

    %overall z position to be plotted
    fitArrayPN(:, 1) = ppval(fitResultsPN, zValues);

    % Correct the original Z positions
    correctedZPositions = (abs(zPositions) + abs(fitArray(:, 1)/1000));
    correctedZPN = (zPositions) + (fitArrayPN(:, 1)/1000);
    disp(correctedZPN);

    % Create a new corrected locs3D matrix
    locs3D_corrected = locs3D;
    
    % Replace Z positions in locs3D with corrected values
    locs3D_corrected(:, 3) = correctedZPN;

    % Extract corrected Z positions
    correctedZPN = locs3D_corrected(:, 3);

    % Calculate total error in micrometers for corrected data
    totalErrorZ_corrected = (zValues' - correctedZPN) * 1000;

    % Plotting (optional, you can remove this part if not needed)
    figure;
    % Plotting Line for Original Data
    plot(zValues, totalErrorZ, '^-', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Z Error (Original)', 'Color', [0, 0, 1]);
    hold on;
    % Plotting Line for Corrected Data
    plot(zValues, totalErrorZ_corrected, 'o-', 'LineWidth', 1, 'MarkerSize', 3, 'DisplayName', 'Z Error (Corrected)', 'Color', [1, 0, 0]);
    % Customize Plot
    xlabel('Axial Range (um)');
    ylabel('Axial error (nm)');
    title('Z-axis Localization Accuracy');
    % Customize Legend
    legend('show');
    hold off;
end
