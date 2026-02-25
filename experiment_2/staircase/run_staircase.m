% run_staircase
% experiment 2

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
p.training = 0;
p.use_staircase = 0;
p.longer_stim_dur = 0;

% Sync Test
if sum(p.which_setup == [0 2]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end


%% Set directories

dirs.script_dir = pwd;
dirs.functions_dir = '../analyses/functions'; addpath(dirs.functions_dir);


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

%% Setup devices and display, open window

init_device_input
init_display
open_window


%% Loading screen

%% Define stimulus
% For precision (ie filter width), we could use a template noise sample.
%       staircase
% For contrast, we can just use template filtered noise texture ().

% max_contrast = 0.9; % Michelson contrast
% min_contrast = 0.05; % Michelson contrast
% max_filter_width = 85; % °
% min_filter_width = 2; % °

%%% Performance goals (orientation discrimination):
% lvl 1 = 90%
% lvl 2 = 75%
% lvl 3 = 65%
% (or we do linspace(lvl_1_perf, lvl_3_perf, num_lvls))


%%% Mechanics
% num_trials = 60;
% probe_offsets = round(linspace(0,15,7)); % ie  [0 3 5 8 10 13 15];



%% STEP 1: Simulate psychometric function (contrast)

contrast_vals = 0:0.01:1;
sim_psych_func_params = [];

free_params = sim_psych_func_params;
fixed_params = contrast_vals;
sim_psych_func = calcPsychometricFunction(free_params, fixed_params);

%% Init Textures: Toggles 

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

%% Calibration Loop

%Participants perform blocks of an orientation discrimination task.
%Loop through each feature, using feature as an index (e.g., 1 or 2)
%Within a block of trials, the contrast or filter width is randomly changing from trial to trial.
%For every trial, participants perform orientation discrimination task with a pseudorandom probe offset (see experiment_1/experiment/script_modules/experiment and trials_loop)
%Present test stimulus
%Present mask
%Present delay
%Present probe line
%Response period 
%Responses need to be stored with respect to feature and level 
%Inter-trial interval (ITI; i.e., small delay in between trials)


%% STEP 2: Response Simulation (contrast)

n_trials = 60;

probe_offsets = round(linspace(0,15,7));

offsets = datasample(probe_offsets, n_trials);
offsets = offsets .* datasample([-1 1], n_trials);
correct_ans = offsets > 0;

contrast = nan(n_trials,1);
resp     = nan(n_trials,1);   % 1 = correct, 0 = wrong
button_press = zeros(n_trials,1); % 1 = CW, 0 = CCW

min_contrast = 0.05;
max_contrast = 0.9;

contrast(1) = max_contrast;
step = 0.1;
min_step = 0.01;
max_step = 0.1;
step_tracker = nan;

% define observer bias
mu = 0;
sigma = 0.5;

num_correct = 0;


for t = 1:n_trials-1

    probe_offset = offsets(t);
    correct_resp = probe_offset > 0; % 1 = positive (CW), 0 = negative (CCW)

    p_CW = calc_pCW(probe_offset, mu, sigma, 0.1);
    button_press(t) = rand() < p_CW; % 1 = CW, 0 = CCW

    resp(t) = button_press(t) == correct_resp;

    if resp(t) == 1
        num_correct = num_correct + 1;
    elseif resp(t) == 0
        num_correct = 0;
        if ~isnan(step_tracker) && step_tracker == -1
            step = step/2;
            step = min(max(step,min_step),max_step);
        end
        step_tracker = 1;
        contrast(t+1) = contrast(t) + step;  % easier
    end

    % staircase rule (2-up-1-down)
    if num_correct == 2
        num_correct = 0;
        if ~isnan(step_tracker) && step_tracker == 1
            step = step/2;
            step = min(max(step,min_step),max_step);
        end
        step_tracker = -1;
        contrast(t+1) = contrast(t) - step;  % harder
    elseif num_correct == 1
        contrast(t+1) = contrast(t);
    end

    % staircase rule (3-up 1-down)

    % clamp
    contrast(t+1) = min(max(contrast(t+1),min_contrast),max_contrast);
end

figure;
yyaxis left
plot(1:n_trials,contrast,'LineWidth',1.5,'Marker','.')
ylabel('Contrast')
ylim([0 1])
yyaxis right
scatter(1:n_trials,resp,100,'LineWidth',1.5, 'Marker','.')
ylabel('Response (1=correct, 0=incorrect)')
xlabel('Trial')
title('Simulated contrast staircase')
box off;
set(gca, 'TickDir', 'out');

%% Estimate psychometric function

% Find all the unique contrast values, get the corresponding responses, then calculate performance for each value.
% This is what's fed into the psychometric function estimation.

unique_contrast_vals = unique(contrast);
performance = nan(length(unique_contrast_vals),1);
for i = 1:length(unique_contrast_vals)
    performance(i) = mean(resp(contrast == unique_contrast_vals(i)));
end

% Weibull psychometric function?
psych_func = @(free_params) calcPsychometricFunction(free_params, fixed_params);

init_params = [0.5, 1, 0.5, 0]; % this is what you're estimating (max performance, slope, C_50, baseline)
lower_bounds = [0, 0, 0, 0]; % (max performance, slope, C_50, baseline)
upper_bounds = [1, 10, 10, 1]; % (max performance, slope, C_50, baseline)

options = optimoptions('fmincon');
[params_est, sse, exit_flag] = fmincon(psych_func, init_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

est_psych_func = psych_func(contrast_vals, params_est);


perc_correct = [0.65 0.75 0.90];

% get the x-values that correspond to the performance goals
x_vals = interp1(contrast_vals, est_psych_func, perc_correct); % look into this...

% Plot psychometric function and data
% Should be s-shaped curve that fits the data.
figure('Color', 'w');
plot(contrast_vals, est_psych_func, 'LineWidth', 1);
hold on;
scatter(contrast, resp, 100, 'LineWidth', 1, 'Marker','.');
xlabel('Contrast');
ylabel('Response probability');
title('Psychometric function');
box off;
set(gca, 'TickDir', 'out');

%% Turn off Kb and restore display

KbQueueStop(p.device_number);
if p.which_setup == 1 && w.gamma_correct 
    Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
end
sca; ShowCursor;
