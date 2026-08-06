%%% simulate_calibration_subject.m
% Generates mock calibration runs for experiment_2

clear all; close all; clc;

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles
save_data = 1;
toggles.which_setup = 0;
toggles.show_debug_output = 1;
toggles.half_screen = 1;
toggles.simulate_response = 1;
toggles.training = 0;
toggles.calibration = 1;
toggles.level_type = 0;
toggles.simulation_mode = 1;
toggles.demo_run = 0;

%% Set directories
script_dir = fileparts(mfilename('fullpath'));
functions_dir = '../experiment/functions'; addpath(functions_dir);
data_dir = fullfile(script_dir, 'data');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
addpath(data_dir);

%% Set subjects
num_subjs = 20;

subj_IDs = cell(num_subjs,1);
for subj = 1:num_subjs
    subj_IDs{subj} = num2str(9900 + subj - 1);
end

%% Experiment parameters

p.feature_name = {'contrast', 'filter'};
p.num_features = numel(p.feature_name);
p.num_trials_per_lvl_per_block = 16; % matching recommended 10 trials per level per block for a 90 min session
p.num_blocks = 8;
p.num_blocks_per_feature = p.num_blocks / p.num_features;

%% Define stimulus parameters

p.num_levels = 7;
p.contrast_min = 0.10;
p.contrast_max = 0.9;
p.contrast = round(logspace(log10(p.contrast_min), log10(p.contrast_max), p.num_levels),2);

p.filter_width_min = 10;
p.filter_width_max = 80;
p.orientation_bp_filter_width = round(logspace(log10(p.filter_width_min), log10(p.filter_width_max), p.num_levels),2);

p.orientation_min = 0;
p.orientation_max = 179;
p.probe_offsets = round(linspace(0,15,7));

%% Simulate subjects
for subj = 1:num_subjs

    p.subj_ID = subj_IDs{subj};

    subj_dir = fullfile(data_dir, p.subj_ID);
    if ~exist(subj_dir, 'dir')
        mkdir(subj_dir);
    end

    % True parameters for this subject
    % Contrast (Weibull)
    true_alpha_c = 0.1 + (0.4 - 0.1) * rand();
    true_beta_c = 1.5 + (3.0 - 1.5) * rand();

    % Filter width (Weibull on precision = 1/fw)
    % precision = 1/fw ranges from ~0.0125 (fw=80) to ~0.1 (fw=10)
    true_alpha_fw = 0.02 + (0.07 - 0.02) * rand(); % threshold in precision units
    true_beta_fw  = 1.5 + (3.0 - 1.5) * rand();    % slope

    % Fixed parameters
    true_guess_rate = 0.5;
    true_lambda = 0.01;

    p.true_params.alpha_c = true_alpha_c;
    p.true_params.beta_c = true_beta_c;
    p.true_params.alpha_fw = true_alpha_fw;
    p.true_params.beta_fw = true_beta_fw;
    p.true_params.guess_rate = true_guess_rate;
    p.true_params.lambda = true_lambda;

    %% Define experiment

    p.feature_order = repmat(1:p.num_features, 1, p.num_blocks_per_feature);
    p.feature_order = Shuffle(p.feature_order);

    % Generate trial events
    for n_block = 1:p.num_blocks

        % If you are passing p.num_trials_per_lvl_per_block = 3 directly into BalanceFactors in init_experiment.m,
        % you are actually doing 3 * 10 levels = 30 trials *per block*.
        level_order = BalanceFactors(p.num_trials_per_lvl_per_block, 1, 1:p.num_levels);

        if n_block == 1
            p.num_trials_per_block = length(level_order);
            p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks);
            p.correct_response = nan(p.num_trials_per_block, p.num_blocks);
        end

        % Sample Test orientations
        test_orientations = sample_orientation(p.orientation_min, p.orientation_max, p.num_trials_per_block);

        % Probe orientations
        probe_offsets = datasample(p.probe_offsets, length(test_orientations));
        probe_orientations = calc_probe_orientation(test_orientations, probe_offsets');

        % Storing trial events
        p.trial_events(:,:,n_block) = [test_orientations, probe_orientations, level_order];

    end

    p.num_trials = p.num_trials_per_block * p.num_blocks;

    %% Simulate responses

    [behav_data.responses, behav_data.correct] = simulate_calibration_responses(p);

    %% Save the data

    save_filename = ['SD_Noise_Exp2_Calibration_S' p.subj_ID '_Run1.mat'];

    run_info.behav_data = behav_data;
    run_info.p = p;
    run_info.toggles = toggles;

    if save_data
        save(fullfile(subj_dir, save_filename), 'run_info', '-mat', '-v7.3');
        disp(['Saved ' save_filename]);
    end

end
