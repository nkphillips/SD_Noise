%%% Initialize stimulus parameters

%{   

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define aperture size
% 6° diameter

p.aperture_height_px = round(w.screen_height_px / 2);
p.aperture_width_px = p.aperture_height_px;

p.aperture_radius_px = w.ppd * 6;

%% Define noise texture size
% A bit larger than the aperture so that the noise is not cut off at the edges

p.height_px = p.aperture_height_px;
p.width_px = p.aperture_width_px;

%% Define contrasts

if p.training
    p.contrast = 0.9; 
else
    p.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast
end

%% Define orientation and spatial frequency bandpass filter widths

if p.training
    p.orientation_bp_filter_width = 0.1;
else
    p.orientation_bp_filter_width = [0.1 5 20]; % low noise, medium, high ; unit: °
end

p.sf_bp_filter_cutoffs = [1 04]; % unit: cycles/degree

%% Check that the number of levels between stimulus contrast and bp filter widths match

if length(p.orientation_bp_filter_width) == length(p.contrast)
    p.num_levels = length(p.contrast);
else
    disp('Condition levels do not match in length!');
end

%% Define orientations

p.orientation_min = 0;
p.orientation_max = 179;

%% Define number of noise samples

p.num_test_samples = 20;
p.num_mask_samples = 20;

%% Define probe line 

p.probe_length = p.fixation_space_aperture_px;
p.probe_thickness = w.ppd * 0.1;
p.probe_color = w.black;

p.probe_start = w.centerX - p.probe_length;
p.probe_end = w.centerX + p.probe_length;

stimuli.probe_line_base = [w.centerX - p.probe_length, w.centerX + p.probe_length; ...
    w.centerY, w.centerY];
