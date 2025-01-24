%%% Initialize stimuli parameters

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define range of orientation values

p.min_orientation = 0; % °
p.max_orientation = 180; % °

p.num_orientations = 100;

p.orientations = linspace(p.min_orientation, p.max_orientation,  p.num_orientations); % gen_lin_log_sf(p.min_sf, p.max_sf, p.num_sf);

p.num_noise_samples = 10;

%% Define size and contrast

stimuli.diameter_deg = 10;
stimuli.diameter_px = floor(stimuli.diameter_deg * w.ppd);

stimuli.contrast = 0.8;

%% Define aperture

% Aperture alpha
stimuli.aperture_alpha = 255; % 0 = completely transparent; 255 = completely opaque

% Aperture size
stimuli.aperture_diameter_px = stimuli.diameter_px;

stimuli.mask_radius_deg = 3; % This effectively defines how much of the stimuli behind the aperture is seen, so it should match the desired stimulus size
stimuli.mask_radius_px = stimuli.mask_radius_deg * w.ppd;

%% Define SF bandpass filter width and smoothening

stimuli.bp_filter_width = 0.1; % default = 0.1
stimuli.bp_filter_gauss_sd = 0.1 * w.ppd; % default = 0.1 * ppd
