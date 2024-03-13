function [z_calib, offset] = calibSMLFM(calibFilePath, calibFolderPath, calibStep, rangeToFit, type, save, locs_format, pixel_size_sample, type_mla, f_MLA, lens_pitch, mlaCentrePos, sizeOptic, NA, f_obj , f_TubeLens, f_FourierLens, pixel_size, nImmersion, nMedium)
    
    % SMLFM fitting of calibration
    % 1. Read 2D calibration localisation file
    locs_2d = File.readLocalisationFile(calibFilePath, locs_format, pixel_size_sample);
    outputFolder = calibFolderPath;
    disp(calibFolderPath)
    mkdir(outputFolder);

    % 2. Rotate x and y
    % Define theta based on the type of data; simulated or experimental
    if strcmp(type, 'experimental')
        theta = 2 * pi / 180; % change to match the orientation of the MLA for experimental data
    elseif strcmp(type, 'simulation')
        theta = 32 * pi / 180; % change to match the orientation of the MLA for simulation data
    else
        error('Invalid type specified.'); % Handle the case when 'type' is neither 'experimental' nor 'simulation'
    end
    x = locs_2d(:, 2);
    x = x - mean(x);
    y = locs_2d(:, 3);
    y = y - mean(y);
    locs_2d(:, 2) = (x .* cos(theta) - y .* sin(theta));
    locs_2d(:, 3) = (x .* sin(theta) + y .* cos(theta));

    % Scatter plot of rotated data at the image plane sample_pixel
    scatter(locs_2d(:, 2), locs_2d(:, 3), [], '.');
    figure(1);
    xlabel('x (\mu m)');
    ylabel('y (\mu m)');
    title('2D Localisations');

    % 3. Initialise MLA, Microscope, and LightFieldLocalisation objects
    % Create an instance of the MicroLensArray class
    microLensArray = Classes.MicroLensArray(type_mla, f_MLA, lens_pitch, mlaCentrePos, sizeOptic);    
    lfm = Classes.FourierLFM(NA, f_obj, f_TubeLens, f_FourierLens, pixel_size, nImmersion, nMedium, microLensArray);
    lfLocs = Classes.LightFieldLocalisations(locs_2d, microLensArray, lfm);
    
    % 4. Filter and rotate localisations
    lfLocs = lfLocs.resetFilteredLocs;
    lfLocs = lfLocs.filterRho([0 0.6]); % 0 0.6
 
    % Assign colours to microlenses - break code here to check if microlenses are correctly assigned
    % image plane sample_pixel using 114 nm - unsure units
    figure(2);
    lfLocs.plotXYs(3);
    xlabel('x (\mu m)');
    ylabel('y (\mu m)');
    title('Colour Assignment 2D Locs');
 
    % % Plot UV space with microlens array gives normalised pixel coordinates
    lfLocs.plotUVs;

    % 5. Find system aberrations
    fit_params = {};
    fit_params.frame_range = [lfLocs.minFrame, 1000];
    fit_params.max_disparity = 10; % find locs from -5 to 5 um
    fit_params.dist_search = 0.5;
    fit_params.angle_tol = 2 * pi / 180;
    fit_params.threshold = 3;
    fit_params.min_views = 4; % 5
    
    mla_rotation = -0 * pi / 180;
    lfLocs = lfLocs.resetFilteredLocs;
    
    % Plot XY localisations for further filtering
    % lfLocs = lfLocs.filterRho([0 0.6]);
    % lfLocs = lfLocs.filterSpotSize([0.1 1]);
    lfLocs.plotXYs(3)

    % Rotate UV localisations
    lfLocs = lfLocs.rotateUV(mla_rotation);
    lfLocs = lfLocs.setAlpha('integrateSphere');
    
    abberation_params = {};
    abberation_params.axial_window = 6; %sets z range for aberration
    abberation_params.min_views = 5;
    abberation_params.photon_thresh = 1;

    fprintf('Fitting first 1000 frames for aberration correction')
    numWorkers = 8;

    % Perform fitting for aberration correction
    [locs3D, fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs, lfm, fit_params, numWorkers);
    correction = Fitting.calculateViewError(fit_data, lfm, lfLocs.filteredLocs, abberation_params);

    lfLocs = lfLocs.correctUV(correction);
    fprintf('Global aberration calculated (um): views listed below (u, v, dx, dy)')
    correction(:, 1:4)

    % 6. Fit full data set on corrected localisations
    fit_params.frame_range = [lfLocs.minFrame, lfLocs.maxFrame];
    fit_params.max_disparity = 10;%8
    fit_params.dz = 0.5;%0.5 red increase for orange and other colours to 6 
    fit_params.angle_tol = 2 * pi / 180;
    fit_params.threshold = 3;
    fit_params.min_views = 4;
    fprintf('Fitting the whole data set')
    [locs3D, fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs, lfm, fit_params, numWorkers);
    fprintf('locs3D values:\n');
    fprintf(' %f\n', locs3D);   
    
    fprintf('Done!\n\n')
    
    fprintf('Calib: Total number of 2D localisations used for fitting: %d\n', sum(locs3D(:, end-2)))
    fprintf('Calib: Final number of 3D localisations: %d\n', size(locs3D, 1))

    % 7. Correcting intended vs extracted positions   
    sorted_locs3D = sortrows(locs3D, 3);  % Sort locs3D lowest to highest in order for calibraion correction
    zValues = rangeToFit(1):calibStep:rangeToFit(2);
    numSteps = (rangeToFit(2) - rangeToFit(1)) / calibStep + 1;
    zFocal = ceil((numSteps/2)); 

    % Call the error correction function to get corrected Z positions
    % Includes linear calibration plot, fitting and extraction of values 
    [z_calib, offset, CorrectedzPositions, calibration_errors] = errorCorrectionFunctions(zValues, sorted_locs3D, zFocal);
    % Update locs3D with corrected Z positions
    
    % plot XY
    plotLocalizationErrors(zValues, locs3D);
    

    % 8. Write the result
    % Save original and corrected Z positions

    switch save % if save = 'yes'
        case 'yes'
        fprintf('Writing localisations to output files...\n')
        
        % Original z positions from SMLFM algorithms 
        % VISP = [locs3D(:, 1) * 1000, locs3D(:, 2) * 1000, locs3D(:, 3) * 1000, locs3D(:, 7), locs3D(:, 8)];
        % dlmwrite([outputFolder 'VISP_(' datestr(now, 'dd-mm-yyyy_MM-HH') ').3d'], VISP, 'delimiter', '\t', 'precision', 7);
        writematrix(locs3D, [outputFolder 'original_locs3D.csv']);

        % Corrected z positions from linear fit
        % VISP = [locs3D(:, 1) * 1000, locs3D(:, 2) * 1000, CorrectedzPositions * 1000, locs3D(:, 7), locs3D(:, 8)];
        % dlmwrite([outputFolder 'VISP_(' datestr(now, 'dd-mm-yyyy_MM-HH') ').3d'], VISP, 'delimiter', '\t', 'precision', 7);
        locs3D(:, 3) = CorrectedzPositions;
        writematrix(locs3D, [outputFolder 'corrected_locs3D.csv']);
        % plots for calibration fit and accuracy error  
        fprintf('Saving plots...\n')
        saveas(figure(3), [outputFolder 'calibFit.fig']);
        saveas(figure(4), [outputFolder 'XYerror.fig']);
        saveas(figure(5), [outputFolder 'Zerror.fig']);
        fprintf('Calibration successfully analysed \n')

        case 'no'
    
    end
end

