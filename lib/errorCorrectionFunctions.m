function [z_calib, offset, CorrectedzPositions, calibration_errors] = errorCorrectionFunctions(zValues, sorted_locs3D, zFocal)

    % Extract Z positions
    zPositions = sorted_locs3D(:, 3) - sorted_locs3D(zFocal, 3);
    
    % plot of extracted vs intended z position 
    % determines calibration factor for original data 
    figure(3);
    plot(zPositions, zValues, '^-', 'LineWidth', 1, 'MarkerSize', 3);
    xlabel('Intended Range (um)');
    ylabel('Extracted Range (um)');
    title('Z-axis Localization Accuracy');
    
    % Perform linear regression
    coefficients = polyfit(zPositions, zValues, 1);
    slope = coefficients(1);
    intercept = coefficients(2);
    
    % Evaluate the linear fit
    fitLine = slope * zPositions + intercept;
    
    % Plot the linear fit
    hold on;
    plot(zPositions, fitLine, 'r-', 'LineWidth', 2);
    legend('Data', sprintf('Fit: y = %.4fx + %.4f', slope, intercept));
    hold off;

    z_calib = slope;
    offset = intercept; 

    % Apply correction
    CorrectedzPositions = zPositions * z_calib;
    figure(4);
    plot(CorrectedzPositions, zValues, '^-', 'LineWidth', 1, 'MarkerSize', 3);
    xlabel('Intended Range (um)');
    ylabel('Extracted Range (um)');
    title('Z-axis Localization Accuracy');

    % Calculate the errors
    calibration_errors = diff(zValues) - diff(CorrectedzPositions); 

    % Plot the errors
    figure(5);
    plot(zPositions(1:end-1), calibration_errors * 1000, 'o-', 'LineWidth', 1, 'MarkerSize', 3); % Convert to nanometers here
    xlabel('Intended Range (um)');
    ylabel('Calibration Errors (nm)');
    title('Calibration Errors');
        
    % Check if calibration error is acceptable
    max_error = max(abs(calibration_errors));
    if max_error < 0.1
        fprintf('Calibration successful. Maximum error: %.2f nm\n', max_error * 1000); % Convert to nanometers here
    else
        fprintf('Calibration error is too high (%.2f nm). Repeat calibration.\n', max_error * 1000); % Convert to nanometers here
    end
    
end
