% run_calibration

%% Clean and prepare workspace

clc;
clear all; %#ok<CLALL>
close all

commandwindow; % force cursor to command window
% Priority(1); % Set MATLAB/Psychtoolbox to "high" priority level

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.which_setup = 0; % 0 = MacBook, 1 = 3329C_ASUS, 2 = S32D850
p.disp_on = 0;
p.half_screen = 1;
p.simulate_response = 0;

% Sync Test
if sum(p.which_setup == [0 2]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

%% Set directories

p.subj_ID = '999';

dirs.project_dir = '../'; addpath(dirs.project_dir);
dirs.script_dir = pwd;
dirs.functions_dir = 'functions'; addpath(dirs.functions_dir);
dirs.data_dir = '../data';
dirs.texture_dir = 'textures';

if exist(dirs.data_dir,'dir') == 0
    mkdir(dirs.data_dir);
end

if exist(dirs.texture_dir,'dir') == 0
    mkdir(dirs.texture_dir);
end

addpath(dirs.data_dir);
addpath(dirs.texture_dir);

dirs.init_dir = 'init'; addpath(dirs.init_dir);
dirs.modules_dir = 'script_modules'; addpath(dirs.modules_dir);
dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/serenceslabexp/Desktop/MonitorCalibration/GammaTables'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

%init_device_input
%init_display
%open_window

%% Initialize device input

% Enable unified mode of KbName, so KbName accepts identical key names on
% all operating systems:
KbName('UnifyKeyNames');

%% Define device and relevant keys

p.device_number = 0;
[Kb_indices, product_names, ~] = GetKeyboardIndices;

switch p.which_setup

    case 0 % If using Macbook

        p.device_string = 'Apple Internal Keyboard / Trackpad'; % Macbook
        %     p.device_string =  'USB-HID Keyboard'; % external keyboard

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

    case 1 % If using 3329C

        p.device_string = 'LiteOn Lenovo Traditional USB Keyboard';

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

    case 2 % If using S32D850

        p.device_string =  'USB-HID Keyboard'; % external keyboard

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

end

%% Scan for device number

for i = 1:length(product_names)
    if strcmp(product_names{i}, p.device_string)
        p.device_number = Kb_indices(i);
        break;
    end
end

%% Error if no device is found

if p.device_number == 0
    error('No device by that name was detected');
end

%% Turn on keyboard input

KbQueueCreate(p.device_number);
KbQueueStart(p.device_number);

%% Initialize display parameters
% Ideal distance: 1 cm equals 1 visual degree at 57 cm

screens = Screen('Screens'); % Grab the available screens

%% Set screen parameters

if p.which_setup == 0 % Macbook

    p.display_setup = 'Macbook';

    w.use_screen = min(screens); % If there are two or more displays, 'min' should grab the most internal display.
    w.view_distance = 57; % in centimeters, cm
    w.screen_width = 30; % cm
    w.screen_width_px = 1512; % in pixels, px
    w.screen_height_px = 982; % px

elseif p.which_setup == 1 % 3329C_ASUS

    p.display_setup = '3329C_ASUS';

    w.use_screen = 0;
    w.gamma_correct = 1;
    w.view_distance = 42; % cm; default = 42
    w.screen_width = 58; % cm
    w.screen_width_px = 2560; % px
    w.screen_height_px = 1440; % px

elseif p.which_setup == 2 % S32D850

    p.display_setup = 'S32D850';

    w.use_screen = min(screens); % If there are two or more displays, 'max' should grab the most external display.
    w.view_distance = 57; % in centimeters, cm
    w.screen_width = 70.8; % cm
    w.screen_width_px = 2560; % in pixels, px
    w.screen_height_px = 1440; % px

end

%% Calculate visual angle of the display, pixels per degree of visual angle, size of a pixel, and Nyquist frequency

screen_length = w.screen_width;
screen_length_px = w.screen_width_px;

w.visual_angle = 2 * atan2d(screen_length/2,  w.view_distance); % Visual angle of the whole screen in degrees
w.ppd = screen_length_px/w.visual_angle; % Pixels per degree of visual angle
w.px_size = screen_length/screen_length_px; % size of pixel in cm
w.f_Nyquist = 1/(2*w.px_size);

%% Half screen

if p.half_screen
    w.view_distance = 57; % in centimeters, cm
    w.screen_width_px = w.screen_width_px/2; % in pixels, px
    w.screen_height_px = w.screen_height_px/2; % px
end

%% Define the background color

w.gray = 127; % halfway point between 0–254 (255 elements total)
w.bg_color = w.gray;

%% Define other colors

w.white = [255 255 255];
w.black = [0 0 0];
w.red = [255 0 0];
w.green = [0 255 0];
w.blue =  [0 0 255];

%% Open Window

w.window = PsychImaging('OpenWindow', w.use_screen, w.bg_color, [0 0 w.screen_width_px w.screen_height_px]);

if ~p.half_screen
    HideCursor;
    commandwindow;
end

%% Enable alpha blending

Screen('BlendFunction', w.window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

%% Load CLUT

w.DefaultCLUT = Screen('ReadNormalizedGammaTable', w.window);

if p.which_setup == 1 && w.gamma_correct

    load([dirs.monitor_cal_dir '/corrected_gamma_table_' p.display_setup '.mat'])
    w.CorrectedCLUT = corrected_gamma.table * 255;
    Screen('LoadCLUT', w.window, w.CorrectedCLUT);

end

%% Define center coordinates

w.centerX = w.screen_width_px/2;
w.centerY = w.screen_height_px/2;

%% Text settings

Screen('TextStyle', w.window, 1); % 0=normal, 1=bold, 2=italic
Screen('TextSize', w.window, 18);



%% Probe offset magnitudes

p.probe_offsets = round(linspace(0,15,7));

%% Initialize trackers

exit_session = 0;

%% Define stimuli parameters



p.num_levels = 10;
p.num_reps = 30;
p.num_trials_per_feature = p.num_levels * p.num_reps;

p.contrast_min = 0.05;
p.contrast_max = 0.9;
p.calibration_contrast_levels = round(logspace(log10(p.contrast_min), log10(p.contrast_max), p.num_levels),2);

p.filter_width_min = 2;
p.filter_width_max = 80;
p.calibration_filter_width_levels = round(logspace(log10(p.filter_width_min), log10(p.filter_width_max), p.num_levels),2);


p.num_noise_samples = 10;   % noise exemplars per level


%% Create stimuli textures

%init_calibration_textures

%%% init_textures

tic

%% Texture Toggles

textures_filename = ['SD_Noise_textures_' p.display_setup '.mat'];

if ~exist([dirs.texture_dir '/' textures_filename], 'file')
    generate_textures = 1;
    save_textures = 1;
else
    generate_textures = 0;
    save_textures = 0;
end

if p.training
    save_textures = 0;
end


% --- Calibration aliases so init_textures-style code runs ---
p.training = 0; % prevent undefined variable error

p.contrast = p.calibration_contrast_levels;
p.orientation_bp_filter_width = p.calibration_filter_width_levels;

p.num_test_samples = p.num_noise_samples;  % match naming used in init_textures



%% Aperture
% alpha level for aperture:
% 0 = completely transparent (the texture of the aperture is invisible)
% 255 = completely opaque (the texture of the aperture dominates)

aperture = create_circular_aperture(p.aperture_width_px, p.aperture_height_px, p.aperture_radius_px); % texture size, radius of circle
% figure, subplot(1,2,1), imshow(aperture)

aperture = imgaussfilt(aperture, 0.1 * w.ppd);
% subplot(1,2,2), imshow(aperture)
aperture_texture(:,:,1) = ones(size(aperture)) * w.gray;
aperture_texture(:,:,2) = aperture * 255; 

% figure, imshow(aperture_texture(:,:,2), [0 255])

%% Generate textures 

if generate_textures

    disp('Generating stimuli...')

    % Preallocate textures
    noise_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_test_samples);
    p.test_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_test_samples);
    p.mask_textures = nan(p.height_px, p.width_px, length(p.contrast), p.num_mask_samples);

    for i = 1:size(noise_textures, 3) % Contrasts
        for j = 1:size(noise_textures, 4) % Filter widths
            for k = 1:size(noise_textures, 5) % Samples
                
                base_noise = create_noise_texture(p.height_px, p.width_px);
                base_noise = bandpassFilterImg(base_noise, [0, 180], [0.5 6], w.ppd * 0.1, w.f_Nyquist);
                base_noise = centerTextureContrast(base_noise, p.contrast(i), w.gray);
                
                if j == 1
                    stimuli.mask_textures(:,:,i,k) = base_noise;
                end
                
                % Make orientation- and spatial frequency-bandpass filtered noise
                noise_texture = bandpassFilterImg(base_noise, [round(180 - p.orientation_bp_filter_width(j)/2), floor(180 + p.orientation_bp_filter_width(j)/2)], p.sf_bp_filter_cutoffs, w.ppd * 0.1, w.f_Nyquist);
                noise_texture = centerTextureContrast(noise_texture, p.contrast(i), w.gray);

                % Ignore certain combos
                if i > 1 && j > 1
                    continue
                end

                stimuli.test_textures(:,:,i,j,k) = noise_texture; % Convert to visible pixel values and scale by contrast

            end
        end
    end

else
    load([dirs.texture_dir '/' textures_filename]);
    stimuli = textures;
    clear textures;

end

disp(['Elapsed time: ' num2str(toc) ' s'])


%% CHAT RECOMMENDS ADDING THE FOLLOWING: Convert generated images to PTB textures (handles you can draw)

stimuli.test_tex = nan(size(stimuli.test_textures,3), size(stimuli.test_textures,4), size(stimuli.test_textures,5));

for i = 1:size(stimuli.test_textures,3) % contrast index
    for j = 1:size(stimuli.test_textures,4) % filter width index
        for k = 1:size(stimuli.test_textures,5) % sample index

            img = stimuli.test_textures(:,:,i,j,k);

            if any(isnan(img(:)))
                continue
            end

            % PTB expects uint8 or double in [0 1]. Your textures look like 0-255.
            img_uint8 = uint8(img);

            stimuli.test_tex(i,j,k) = Screen('MakeTexture', w.window, img_uint8);
        end
    end
end



%% Define stimulus sequences
% feature 1 = contrast
% feature 2 = filter width
%
% useful functions:
% datasample
% randi
% repmat
% nan

p.num_features = 2;
p.feature_name   = {'contrast','filter_width'};

presentation_order = nan(p.num_trials_per_feature, p.num_features); % matrix of indices. whatever approach you use, this var should store the order for both contrast and filter width sequences
responses = nan(p.num_reps, p.num_levels, p.num_features); % page 1: contrast, page 2: filter width

probe_offset_used = nan(p.num_trials_per_feature, p.num_features);


for feature = 1:p.num_features
    tmp = repmat(1:p.num_levels, p.num_reps, 1);
    presentation_order(:,feature) = datasample(tmp(:), p.num_trials_per_feature, 'Replace', false);
end

%% Make stimuli

stimuli_made = nan(p.num_levels, p.num_noise_samples, p.num_features);

for feature = 1:num_features

    % Screen('MakeTexture')

end

%% Calibration loop

trial_counter = zeros(1,p.num_levels);

for feature = 1:p.num_features
    for trial = 1:p.num_trials_per_feature

        % which level to index
        curr_lvl = presentation_order(trial,feature);

        % first indx in responses needs to point to which trial of the curr lvl it is
        trial_counter(curr_lvl) = trial_counter(curr_lvl) + 1;
        curr_lvl_count = trial_counter(curr_lvl);

        % use level index to pull the current stimulus we need to show
        % curr_stimulus = stimuli_made(curr_lvl, curr_noise_lvl, feature); % 3d array of textures / images

        % record response (responses = nan(p.num_reps, p.num_levels, 2);)
        curr_response = 1; % can be 1 or 2 (left or right; CCW or CW)
        responses(curr_lvl_count, curr_lvl, feature) = curr_response;


        % evaluate response

    end
end





% Might break feature 2 indexing, CHAT recommends...
for feature = 1:p.num_features

    trial_counter = zeros(1,p.num_levels);  % reset per feature

    for trial = 1:p.num_trials_per_feature
        curr_lvl = presentation_order(trial,feature);
        trial_counter(curr_lvl) = trial_counter(curr_lvl) + 1;
        curr_lvl_count = trial_counter(curr_lvl);


        curr_offset = datasample(p.probe_offsets, 1);
        probe_offset_used(trial, feature) = curr_offset;


        % store response
        curr_response = 1; % placeholder
        responses(curr_lvl_count, curr_lvl, feature) = curr_response;
    end
end

%% Estimate psychometric function



%% Extract constrast and filter width levels



%% Save data

