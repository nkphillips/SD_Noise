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

%% Set directories
script_dir = pwd;
functions_dir = '../experiment/functions'; addpath(functions_dir);
data_dir = '../data'; addpath(data_dir);

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

    if ~exist([data_dir '/' p.subj_ID],'dir')
        mkdir([data_dir '/' p.subj_ID])
    end

    % True parameters for this subject
    % Contrast (Weibull)
    true_alpha_c = 0.1 + (0.4 - 0.1) * rand();
    true_beta_c = 1.5 + (3.0 - 1.5) * rand();

    % Filter width (PTM; Lu & Dosher, 2008)
    % d' = Signal^g / sqrt((1+Nm^2)*fw^(2g) + Nm^2*Signal^(2g) + Na^2)
    true_signal_fw = 20 + (50 - 20) * rand();
    true_nmul_fw   = 0.10 + (0.45 - 0.10) * rand();
    true_nadd_fw   = 50 + (300 - 50) * rand();
    true_gamma_fw  = 2.0; % Fixed transducer exponent

    % Fixed parameters
    true_guess_rate = 0.5;
    true_lambda = 0.01;

    p.true_params.alpha_c = true_alpha_c;
    p.true_params.beta_c = true_beta_c;
    p.true_params.signal_fw = true_signal_fw;
    p.true_params.nmul_fw = true_nmul_fw;
    p.true_params.nadd_fw = true_nadd_fw;
    p.true_params.gamma_fw = true_gamma_fw;
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

    if save_data
        save([data_dir '/' p.subj_ID '/' save_filename], 'run_info', '-mat', '-v7.3');
        disp(['Saved ' save_filename]);
    end

end
