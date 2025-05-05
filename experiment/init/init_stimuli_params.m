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

%% Define contrasts

p.contrast = [0.9 0.5 0.25]; % high contrast, medium, low; unit: % Michelson contrast

if p.training
    p.contrast = p.contrast(1); 
end

%% Define orientation and spatial frequency bandpass filter widths

% p.orientation_bp_filter_width = [1 5 10]; % low noise, medium, high ; unit: °
p.orientation_bp_filter_width = [2 40 80]; % low noise, medium, high ; unit: °

if p.training
    p.orientation_bp_filter_width = p.orientation_bp_filter_width(1);   
end

p.sf_bp_filter_cutoffs = [1 4]; % unit: cycles/degree

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

p.probe_length = round(1 * w.ppd);
if ~mod(p.probe_length, 2), p.probe_length = p.probe_length + 1; end

p.probe_thickness = round(w.ppd * 0.05);
if ~mod(p.probe_thickness, 2), p.probe_thickness = p.probe_thickness + 1; end

probe_line = ones(p.probe_length) * w.gray;

start_col = round(p.probe_length/2) - floor(p.probe_thickness/2);
end_col = round(p.probe_length/2) + floor(p.probe_thickness/2);
probe_line(:, start_col:end_col) = 0; 

stimuli.probe_line = probe_line;
