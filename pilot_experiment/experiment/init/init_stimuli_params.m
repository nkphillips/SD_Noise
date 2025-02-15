%%% Initialize stimulus parameters

%{   

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define aperture size
% 6° diameter

stimuli.aperture_height_px = w.screen_height_px / 2;
stimuli.aperture_width_px = stimuli.aperture_height_px;

stimuli.aperture_radius_px = w.ppd * 6;

%% Define noise texture size
% A bit larger than the aperture so that the noise is not cut off at the edges

stimuli.height_px = stimuli.aperture_height_px;
stimuli.width_px = stimuli.aperture_width_px;

%% Define contrasts

stimuli.contrast = [0.9 0.5 0.25]; % high contrast, medium, low 

%% Define orientation bandpass filter widths

stimuli.bp_filter_width = [0.1 5 20]; % "0" noise, medium, high ;unit: °

%% Check that the number of levels between stimulus contrast and bp filter widths match

if length(stimuli.bp_filter_width) == length(stimuli.contrast)
    p.num_levels = length(stimuli.contrast);
else
    disp('Condition levels do not match in length!');
end

%% Define orientations

stimuli.orientation_min = 0;
stimuli.orientation_max = 359;

%% Define number of noise samples

p.num_test_samples = 20; % 
p.num_mask_samples = 20;

%% Define probe line length

stimuli.probe_length = p.fixation_space_aperture_px;
stimuli.probe_thickness = w.ppd * 0.1;
stimuli.probe_color = w.black;

stimuli.probe_start = w.centerX - stimuli.probe_length;
stimuli.probe_end = w.centerX + stimuli.probe_length;

stimuli.probe_line_base = [w.centerX - stimuli.probe_length, w.centerX + stimuli.probe_length; ...
    w.centerY, w.centerY];
