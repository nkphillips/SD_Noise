%%% run_staircase

%% 
close all; clear all; clc

commandwindow;

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.demo_run = 0;
p.simulate_response = 0;

if p.simulate_response
    presumed_target = 0.707; % From Garcia-Perez 1998, two-down, one-up rule Table 2
end

p.which_setup = 2; % 0 = MacBook, 1 = 3329B, 2 = 3329C_ASUS, 3 = S32D850

% Sync Test
if sum(p.which_setup == [0 3]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end 

save_staircase_data = 1;

%% Set directories

p.subj_ID = '000'; % dummy subj = 999

dirs.script_dir = pwd;
dirs.init_dir = '../init'; addpath(dirs.init_dir);
dirs.modules_dir = '../script_modules'; addpath(dirs.modules_dir);
dirs.functions_dir = '../../../functions'; addpath(dirs.functions_dir);
dirs.texture_dir = '../textures'; addpath(dirs.texture_dir);
dirs.data_dir = '../../data'; addpath(dirs.data_dir);
dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/pclexp/Documents/Luis/MonitorCalibration'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window
init_fixation

%% Loading screen

loading_text = 'Loading calibration...';

loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Define stimulus

% Define size
stimuli.noise_width_px = w.screen_width_px;
stimuli.noise_height_px = w.screen_height_px;

stimuli.aperture_height_px = w.screen_height_px;
stimuli.aperture_width_px = w.screen_width_px;

stimuli.aperture_radius_px = w.ppd * 4;

stimuli.aperture_alpha = 255;

% Define SFs (must be identical to values used in main task)
p.reference_sfs = 2;
p.inducer_sfs = 2.^([-0.5 0 0.5]) * p.reference_sfs;
p.inducer_sfs = round(p.inducer_sfs,2);

stimuli.base_contrast = 0.7;

%% Initialize staircase

init_staircase

%% Create fixation

fixation_space_made = Screen('MakeTexture', w.window, fixation_space);

%% Create circular aperture
% alpha level for aperture:
% 0 = completely transparent (the texture of the aperture is invisible)
% 255 = completely opaque (the texture of the aperture dominates)

aperture = create_circular_aperture(stimuli.aperture_width_px, stimuli.aperture_height_px, stimuli.aperture_radius_px); % angle, texture size, radius of circle

aperture = imgaussfilt(aperture, 0.05 * w.ppd);

% First layer of aperture is just a gray rectangle the size of the screen
aperture_texture(:,:,1) = ones(size(aperture)) * w.gray;

% The second layer is the smoothened aperture made above
aperture_texture(:,:,2) = aperture * 255;

stimuli.aperture_made = Screen('MakeTexture', w.window, aperture_texture);

%% Create baseline noise textures

% Number of noise samples
p.num_samples = 10; % per SF

% Define SF bandpass filter width and smoothening
stimuli.bp_filter_width = 0.1; % default = 0.1

% Create textures
noise_textures = nan(stimuli.noise_height_px, stimuli.noise_width_px, length(p.inducer_sfs), p.num_samples);

for i = 1:size(noise_textures, 3) % inducer SFs
    for j = 1:size(noise_textures, 4) % Samples

        % Create noise texture
        noise_texture = create_noise_texture(stimuli.noise_height_px, stimuli.noise_width_px);

        % Apply SF bandpass filter
        noise_texture = make_sf_bp_filtered_img(noise_texture, p.inducer_sfs(i), stimuli.bp_filter_width, w.f_Nyquist, w.ppd);

        % Normalize and clip
        noise_texture = normalize_array(noise_texture, 'z-score');
        noise_texture(noise_texture < -2) = -2;
        noise_texture(noise_texture > 2) = 2;

        % Store
        noise_textures(:,:,i,j) = noise_texture / std2(noise_texture);

    end
end

%% Contrast normalization (jointly scale images linearly to visible range 0:255)

% Get min and max values
z_min = min(noise_textures(:));
z_max = max(noise_textures(:));

% Normalize and scale to visible range
textures = ((noise_textures - z_min) / (z_max - z_min) ) * 255;

% Store baseline textures
stimuli.full_contrast_textures = textures;
stimuli.base_contrast_textures = textures * stimuli.base_contrast;

% Make baseline textures for drawing
for i = 1:length(p.inducer_sfs)
    for j = 1:p.num_samples
        stimuli.textures_made(i,j) = Screen('MakeTexture', w.window, stimuli.base_contrast_textures(:,:,i,j));
    end
end

%% Define patches for drawing

make_patches

%% Define timing

% Define stimulus duration
t.stim_dur = 0.5;

% Define contrast delta duration
t.contrast_delta_dur = 0.1;

% Define inter-trial interval duration
t.iti_dur = 1;

% Define noise sample update rate
t.noise_sample_update_rate = 20; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

% Define frame rate and duration
[t.frame_dur, t.frame_rate] = get_framerate(t); % [s, Hz]

% Pre-allocate response duration vector
t.response_dur = nan(1, staircases.num_trials_per_sf * length(p.inducer_sfs));

%% Define frames

% Define total number of frames
frames.frames_count = round(t.stim_dur/t.frame_dur);

% Define number of frames for a noise sample update
frames.noise_sample_update_frames_count = round(t.frame_rate/t.noise_sample_update_rate);

% Define number of frames for a contrast delta update
frames.contrast_delta_frames_count = round(t.contrast_delta_dur/t.frame_dur);

% Define frame onsets
frames.frame_onsets = 0:t.frame_dur:t.stim_dur-t.frame_dur;

% Define frames of the noise sample updates
frames.noise_sample_update = zeros(1, frames.frames_count); 
frames.noise_sample_update(1:frames.noise_sample_update_frames_count:end) = 1; 

% Generate sequence of samples
num_sample_updates = sum(frames.noise_sample_update);
frames.update_noise_sample_seq = nan(staircases.max_trials_per_sc * staircases.num_staircases_per_sf, num_sample_updates, length(p.inducer_sfs));

for n_sf = 1:length(p.inducer_sfs)
    frames.update_noise_sample_seq(:,:,n_sf) = gen_unique_seq([staircases.num_trials_per_sf, num_sample_updates], 1:p.num_samples, p.num_samples/2);
end

% Generate contrast change frames
frames.contrast_delta = zeros(staircases.num_trials_per_sf, frames.frames_count, length(p.inducer_sfs));
for n_sf = 1:length(p.inducer_sfs)
    for n_trial = 1:staircases.num_trials_per_sf

        start_indx = randi([round(frames.frames_count/5), frames.frames_count - frames.contrast_delta_frames_count]);
        end_indx = start_indx + frames.contrast_delta_frames_count;

        frames.contrast_delta(n_trial, start_indx:end_indx, n_sf) = 1;

    end
end 

%% Loading complete

wait_for_trigger

%% Staircase 

detection_staircase

% Calculate the final final_contrast_deltas
 staircases.final_contrast_deltas = nan(1, length(p.inducer_sfs));
for n_sf = 1:length(p.inducer_sfs)

    % Pre-allocate vector to store the mean contrast deltas for a staircase
    reversal_contrast_deltas = nan(staircases.num_staircases_per_sf, 1);

    % Get the contrast levels at all the reversals from both staircases
    for n_sc = 1:staircases.num_staircases_per_sf
        
        % Get the number of reversals
        curr_num_reversals = nnz(staircases.reversal_indices(n_sc,:,n_sf) == 1);
        
        % Get the contrast levels at the reversals
        if curr_num_reversals > 0
            tmp_reversal_contrast_deltas = staircases.contrast_deltas(n_sc, staircases.reversal_indices(n_sc,:,n_sf) == 1, n_sf);
        end

        % Ignore the first reversals with respect to num_reversals_to_consider (counting back from the last reversal)
        tmp_reversal_contrast_deltas = tmp_reversal_contrast_deltas(end-staircases.num_reversals_to_consider+1:end);

        % Average all the reversals across a staircase
        reversal_contrast_deltas(n_sc) = mean(tmp_reversal_contrast_deltas, 'omitnan');

    end
    
    % Average the two levels across staircases
    staircases.final_contrast_deltas(n_sf) = mean(reversal_contrast_deltas);

end

%% Turn off Kb and restore display

KbQueueStop(p.device_number);

Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
if p.which_setup == 1, SetScreenDefault; end
Screen('CloseAll');
ShowCursor;

%% Save staircase data

staircases.inducer_sfs = p.inducer_sfs;

if save_staircase_data
    save_filename = ['staircase_data_S' p.subj_ID '_' t.the_date '_' t.the_time '.mat'];
    save([dirs.data_dir '/' p.subj_ID '/' save_filename], 'staircases');
    disp(['Staircase data saved to ' dirs.data_dir '/' p.subj_ID '/' save_filename]);
end

