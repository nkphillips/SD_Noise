%%% Initialize stimulus parameters

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define aperture size

p.aperture_height_px = w.screen_height_px;
p.aperture_width_px = w.screen_width_px;

p.aperture_radius_px = round(w.ppd * 3);

%% Define noise texture size

p.height_px = round(p.aperture_radius_px*2 + w.ppd);
p.width_px = round(p.aperture_radius_px*2 + w.ppd);

%% Define contrast and filter width levels

if toggles.training
    p.contrast = 0.7; % High enough to be clear, but not the absolute max
    p.orientation_bp_filter_width = 2; % Minimum noise for clearest orientation
elseif toggles.calibration
    p.num_levels = 8; % Default = 8
    p.contrast_min = 0.02;
    p.contrast_max = 0.9;
    p.contrast = round(logspace(log10(p.contrast_min), log10(p.contrast_max), p.num_levels),2);

    p.filter_width_min = 2;
    p.filter_width_max = 180;
    p.orientation_bp_filter_width = round(logspace(log10(p.filter_width_min), log10(p.filter_width_max), p.num_levels),2);

   
elseif toggles.level_type == 1
    p.contrast = [0.2 0.8]; 
    p.orientation_bp_filter_width = [2 80]; 


else
    calib_file = fullfile(dirs.data_dir, p.subj_ID, ['S' p.subj_ID '_calibrated_levels.mat']);
    if exist(calib_file, 'file')
        load(calib_file, 'calib');
        p.num_levels = length(calib.target_levels);
        p.contrast = calib.contrast_levels;
        p.orientation_bp_filter_width = calib.filter_width_levels;
        disp(['Loaded calibrated levels for S' p.subj_ID]);
    else
        error(['Calibration file not found for S' p.subj_ID '. Run fit_calibration.m first!']);
    end
end

p.sf_bp_filter_cutoffs = [1 4]; % unit: cycles/degree

%% Check that the number of levels between stimulus contrast and bp filter widths match

if length(p.orientation_bp_filter_width) == length(p.contrast)
    p.num_levels = length(p.contrast);
else
    error('Condition levels do not match in length!');
end

%% Define orientations

p.orientation_min = 0;
p.orientation_max = 179;

%% Define number of noise samples

p.num_noise_samples = 20;

%% Define probe

p.probe_offsets = round(linspace(0,15,7));

p.probe_length = round(2 * w.ppd);
if ~mod(p.probe_length, 2), p.probe_length = p.probe_length + 1; end

p.probe_thickness = round(w.ppd * 0.05);
if ~mod(p.probe_thickness, 2), p.probe_thickness = p.probe_thickness + 1; end


