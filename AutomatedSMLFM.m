% AutomatedSMLFM.m
% Authors: Olivia J. Dovernor (1)
% 1. Yusuf Hamied Department of Chemistry, Lensfield Road, University of
%    Cambridge, Cambridge, CB2 1EW, UK

% Purpose:
% The script allows the user to correct for a non-linear calibration. 
% Non-linear calibrations arise due to spherical aberration present in
% the system. 

% Inputs:
% Calibration Tiff.csv (10 frames at each step, averaged and analysed in FIJI) 
% Select csv to be analysed 

% Outputs:
% 3D localisations [x y z error (lateral) error (axial) views
% intensity frame] as a .csv, and a VISP-compatible .3d file
% Localisation error from by analysing intended and extracted Z positions
% Saves non-linear fit to be used to correct for Z position with other 
% data collected on same day of analysis. 

% Note: AutomatedSMLFM.m has been tested on MATLAB R2022b and requires the
% Statistics and Machine Learning Toolbox. This script is an extention of
% SMLFM.m

clear vars; close all; clc;
tic
addpath('lib')
save = 'yes'; % save output files: 'yes' or 'no'

% Folder with data
folder = 'F:\CAMRIG\3D_Characterisation\2024_03_12 3Dcalibration\AutomatedAnalysis\data\'; 
type = 'experimental'; %'simulation'
% Calibration data
calibFolder = 'F:\CAMRIG\3D_Characterisation\2024_03_12 3Dcalibration\AutomatedAnalysis\calib\';
calib = 'av_calib.csv';

calibStep = 1; %um
rangeToFit = [-4 4]; %range for fitting, nm for z and radians for angle

type_mla = 'hexagonal'; % 'hexagonal' or 'square'
locs_format = 'Peakfit'; % 'Peakfit', 'Thunderstorm', 'Picasso'
NA = 1.49; % numerical aperture of objective
nImmersion = 1.518; % immersion refractive index
nMedium = 1.33; % specimen/medium refractive index
f_obj = 200/100; % in mm
f_FourierLens = 175; % in mm
f_TubeLens = 200; % in mm
f_MLA = 100; % in mm
lens_pitch = 1400; % in microns
pixel_size = 6.5; % camera pixel size in microns
magnification = (f_TubeLens / f_obj * f_MLA / f_FourierLens); 
pixel_size_sample = (pixel_size / magnification); % pixel size in sample space (microns)
sizeOptic = 10000; % size of MLA optic (in microns)
mlaRotation = 0;
mlaCentrePos = ([0 0] * magnification); % in nm (adjust to suit data plotted on line 100)


%% Calibration %%
[z_calib, offset] = calibSMLFM(fullfile(calibFolder, calib), fullfile(calibFolder), calibStep, rangeToFit, type, save, locs_format, pixel_size_sample, type_mla, f_MLA, lens_pitch, mlaCentrePos, sizeOptic, NA, f_obj , f_TubeLens, f_FourierLens, pixel_size, nImmersion, nMedium);

%% Looping %%
files = dir(folder);  % Get info of files/folders in the current directory
filenames = {files.name}; 

subfolders = filenames([files.isdir]); % Directory names (including . and ..)
subfolders = subfolders(1,3:end);

for i = 1:length(subfolders) % Loop through subfolders
    subfolder = subfolders{i};
    disp(['Analyzing data in folder: ', subfolder]); 
    
    % Construct full path to the data folder
    dataFolderPath = fullfile(folder, subfolder);
    
    % List files in the data folder
    filesInSubfolder = dir(fullfile(dataFolderPath, '*.csv'));
    
    for k = 1:numel(filesInSubfolder) % Loop through each file in the subfolder
        fileName = filesInSubfolder(k).name;
        filePath = fullfile(dataFolderPath, fileName);
        disp(['Processing file: ', fileName]); 
        
        % 1. Read 2D calibration localisation file
        locs_2d = File.readLocalisationFile(filePath, locs_format, pixel_size_sample);
        outputFolder = fullfile(dataFolderPath, 'Results'); % Create an output folder within the subfolder
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
     
        % Assign colours to microlenses
        figure(2);
        lfLocs.plotXYs(3);
        xlabel('x (\mu m)');
        ylabel('y (\mu m)');
        title('Colour Assignment 2D Locs');
     
        % Plot UV space with microlens array gives normalised pixel coordinates
        lfLocs.plotUVs;
    
        % 5. Find system aberrations
        fit_params = {};
        fit_params.frame_range = [lfLocs.minFrame, 1000];
        fit_params.max_disparity = 10;
        fit_params.dist_search = 0.5;
        fit_params.angle_tol = 2 * pi / 180;
        fit_params.threshold = 3;
        fit_params.min_views = 4;
        
        mla_rotation = -0 * pi / 180;
        lfLocs = lfLocs.resetFilteredLocs;
        
        % Rotate UV localisations
        lfLocs = lfLocs.rotateUV(mla_rotation);
        lfLocs = lfLocs.setAlpha('integrateSphere');
        
        abberation_params = {};
        abberation_params.axial_window = 6;
        abberation_params.min_views = 5;
        abberation_params.photon_thresh = 1;
    
        fprintf('Fitting first 1000 frames for aberration correction\n')
        numWorkers = 8;
    
        % Perform fitting for aberration correction
        [locs3D, fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs, lfm, fit_params, numWorkers);
        correction = Fitting.calculateViewError(fit_data, lfm, lfLocs.filteredLocs, abberation_params);
    
        lfLocs = lfLocs.correctUV(correction);
        fprintf('Global aberration calculated (um): views listed below (u, v, dx, dy)\n')
        correction(:, 1:4)
    
        % 6. Fit full data set on corrected localisations
        fit_params.frame_range = [lfLocs.minFrame, lfLocs.maxFrame];
        fit_params.max_disparity = 10;
        fit_params.dz = 0.5;
        fit_params.angle_tol = 2 * pi / 180;
        fit_params.threshold = 3;
        fit_params.min_views = 4;
        fprintf('Fitting the whole data set\n')
        [locs3D, fit_data] = Fitting.lightfieldLocalisation(lfLocs.filteredLocs, lfm, fit_params, numWorkers);
        fprintf('locs3D values:\n');
        fprintf(' %f\n', locs3D);   
        
        fprintf('Done!\n\n')
        
        fprintf('Total number of 2D localisations used for fitting: %d\n', sum(locs3D(:, end-2)))
        fprintf('Final number of 3D localisations: %d\n', size(locs3D, 1))

        % Apply calibration correction to z coordinates 
        locs3D(:, 3) = (locs3D(:, 3) * z_calib) - offset;

        % 3D Scatter plot of localizations with corrected Z positions
        figure();
        x = locs3D(:, 1);
        y = locs3D(:, 2);
        z = locs3D(:, 3);
        n_views = locs3D(:, end-2);
        lateral_err = 1000 * locs3D(:, 4);
        axial_err = 1000 * locs3D(:, 5);
        
        keep = lateral_err < 200 & n_views > 3;
        
        figure(3); % Use a different figure for the scatter plot with corrected Z
        scatter3(x(keep), y(keep), z(keep), 80, z(keep), '.');
        xlabel('x \mum');
        ylabel('y \mum');
        zlabel('z \mum');
        c = colorbar;
        c.Label.String = 'Corrected Z (\mu m)';
        axis equal;
        % Set the color limits to cover the full range of corrected Z values
        clim([min(z(keep)), max(z(keep))]);

        % Histograms
        figure(4)
        subplot(1, 3, 1);
        histogram(lateral_err(keep), 1:5:200);
        xlabel('Lateral fit error (nm)');
        ylabel('Occurrence');
        subplot(1, 3, 2);
        histogram(axial_err(keep), 1:5:200);
        xlabel('Axial fit error (nm)');
        ylabel('Occurrence');
        subplot(1, 3, 3);
        histogram(locs3D(keep, end-1));
        xlabel('Number of photons');
        ylabel('Occurrence');
        set(gcf, 'Position', [100, 100, 1200, 300])
        
        % 2D Histogram of Photons vs. Axial Precision
        figure(5)
        histogram2(locs3D(keep, end-1), axial_err(keep), 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
        xlabel('Photons');
        ylabel('Axial precision');

        % 8. Save the results
        switch save % if save = 'yes'
            case 'yes'
                fprintf('Writing localisations to output files...\n')
                VISP = [locs3D(:, 1) * 1000, locs3D(:, 2) * 1000, locs3D(:, 3) * 1000, locs3D(:, 7), locs3D(:, 8)];
                dlmwrite(fullfile(outputFolder, ['VISP_' datestr(now, 'dd-mm-yyyy_MM-HH') '.3d']), VISP, 'delimiter', '\t', 'precision', 7);
                writematrix(locs3D, fullfile(outputFolder, 'locs3D.csv'));
                
                fprintf('Saving plots...\n')
                saveas(figure(3), fullfile(outputFolder, 'plot3D.fig'));
                saveas(figure(4), fullfile(outputFolder, 'histogram.fig'));
                saveas(figure(5), fullfile(outputFolder, 'precision.fig'));
            case 'no'
                % Do nothing if 'save' is set to 'no'
        end

        close all;
    end  
end
