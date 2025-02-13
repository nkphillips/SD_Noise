%%% Initialize stimulus parameters

%{   

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define size

stimuli.noise_height_px = w.screen_height_px / 2;
stimuli.noise_width_px = w.screen_width_px / 2;

%% Define contrasts

stimuli.contrast = [1 0.5 0.25]; % full contrast, medium, low 

%% Define orientation bandpass filter widths

stimuli.bp_filter_width = [0.1 5 20]; % "0" noise, medium, high ;unit: °

%% check num level

if length(stimuli.bp_filter_width) == length(stimuli.contrast)
    p.num_levels = length(stimuli.contrast);
else
    disp('Condition levels do not match in length!');
end

%% Define orientations

stimuli.orientation_min = 0;
stimuli.orientation_max = 359;

%% Define size
% 6° diameter

stimuli.aperture_height_px = w.ppd * 6;
stimuli.aperture_width_px = stimuli.aperture_height_px;

stimuli.aperture_radius_px = stimuli.aperture_height_px / 2;

stimuli.width_px = stimuli.aperture_height_px;
stimuli.height_px = stimuli.width_px;

%% Define number of noise samples

p.num_test_samples = 20; % 
p.num_mask_samples = 20;
