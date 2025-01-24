%%% make_patches

%% Define fixation space and dot patches (how and where the fixation space and dot will be drawn)

fixation_dot_patch = CenterRectOnPoint([0 0 p.fixation_dot_px p.fixation_dot_px], w.centerX, w.centerY);
fixation_space_patch = CenterRectOnPoint([0 0 p.fixation_space_px p.fixation_space_px], w.centerX, w.centerY);

%% Aperture

aperture_patch = CenterRectOnPoint([0 0 stimuli.aperture_width_px stimuli.aperture_height_px], w.centerX, w.centerY);

%% Noise textures

% noise_L_patch = CenterRectOnPoint([0 0 stimuli.noise_diameter_px stimuli.noise_diameter_px], w.centerX - w.screen_width_px/4, w.centerY);
% noise_R_patch = CenterRectOnPoint([0 0 stimuli.noise_diameter_px stimuli.noise_diameter_px], w.centerX + w.screen_width_px/4, w.centerY);

w.noise_patch_shift = 0; % px

noise_L_patch = CenterRectOnPoint([0 0 stimuli.noise_width_px stimuli.noise_height_px], w.centerX - (w.screen_width_px/4), w.centerY);
noise_R_patch = CenterRectOnPoint([0 0 stimuli.noise_width_px stimuli.noise_height_px], w.centerX + (w.screen_width_px/4), w.centerY);
