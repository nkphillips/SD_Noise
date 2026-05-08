% simulate_experiment.m
% Engine to simulate Experiment 2

clear all; close all; clc;

%% 1. Toggles & Setup

num_subjs = 20;
start_subj_id = 901;

% Hypothesis Configuration
% 'null': Flat SD across all levels
% 'contrast_only': Contrast scales SD, Precision is flat
% 'precision_only': Precision scales SD, Contrast is flat
% 'both_same': Both scale SD identically
% 'both_different': Both scale SD, but with different slopes
hypothesis = 'contrast_only';

% Simulation Features
inject_calib_correlation = true; % If true, makes baseline DoG amplitude dependent on psychometric fits

% Setup Directories
base_dir = fileparts(mfilename('fullpath'));
dirs.project_dir = '../../';
addpath(genpath('../experiment')); % Add experiment functions

dirs.data_dir = fullfile(base_dir, 'data');
if ~exist(dirs.data_dir, 'dir'), mkdir(dirs.data_dir); end

%% 2. Ground Truth Generation

gt_all = struct();

for i = 1:num_subjs
    subj_id = num2str(start_subj_id + i - 1);

    % --- Psychometric Parameters (Calibration) ---
    % Ranges validated by find_plausible_params.m to guarantee all inverted
    % calibration levels stay within physical bounds (contrast < 0.9,
    % filter width < 80 deg, all > 0).

    % Weibull parameters for Contrast
    gt.contrast_alpha = 0.4288 + 0.1017 * randn();
    gt.contrast_alpha = max(0.1236, min(0.7339, gt.contrast_alpha));
    gt.contrast_beta  = 2.7876 + 0.7375 * randn();
    gt.contrast_beta  = max(0.5753, min(5.0000, gt.contrast_beta));

    % Weibull parameters for Filter Width (on precision = 1/fw axis)
    % Same Weibull function as contrast, but the independent variable is
    % precision = 1/filter_width so that performance rises with x.
    % alpha_fw is the threshold in precision units (1/deg).
    % Filter widths range 10-80 deg -> precision range 0.0125-0.10.
    gt.filter_alpha = 0.0513 + 0.0119 * randn();
    gt.filter_alpha = max(0.0157, min(0.0870, gt.filter_alpha));
    gt.filter_beta  = 2.8027 + 0.7324 * randn();
    gt.filter_beta  = max(0.6054, min(5.0000, gt.filter_beta));

    gt.lambda = 0.01; % Lapse rate

    % --- Serial Dependence Parameters ---
    if inject_calib_correlation
        % Baseline SD amplitude correlated with psychometric parameters.
        % Higher alpha_fw (precision threshold) means the observer needs
        % narrower filters to perform well → more noise-sensitive → more SD.
        base_amp = 3 + 50.0 * gt.filter_alpha + 5.0 * gt.contrast_alpha + 0.5 * randn();
    else
        base_amp = 7.0 + 1.5 * randn();
    end
    base_amp = max(5, base_amp);
    base_w = 0.05; % ~33 deg FWHM

    gt.dog_amp = zeros(2, 3); % [Feature (1=Contrast, 2=Precision), Level (1,2,3)]
    gt.dog_w = zeros(2, 3);

    switch hypothesis
        case 'null'
            % Flat SD
            gt.dog_amp(1,:) = [base_amp, base_amp, base_amp];
            gt.dog_amp(2,:) = [base_amp, base_amp, base_amp];
            gt.dog_w(1,:) = [base_w, base_w, base_w];
            gt.dog_w(2,:) = [base_w, base_w, base_w];
        case 'contrast_only'
            % Contrast scales down as level increases (higher contrast = less SD)
            % Level 1 = lowest contrast (hardest), Level 3 = highest contrast (easiest)
            gt.dog_amp(1,:) = [base_amp+2, base_amp, base_amp-2];
            gt.dog_w(1,:)   = [base_w-0.015, base_w, base_w+0.015];

            gt.dog_amp(2,:) = [base_amp, base_amp, base_amp]; % Precision flat
            gt.dog_w(2,:)   = [base_w, base_w, base_w];
        case 'precision_only'
            % Level 1 = lowest noise (easiest), Level 3 = highest noise (hardest)
            % Therefore, higher level = more SD
            gt.dog_amp(1,:) = [base_amp, base_amp, base_amp]; % Contrast flat
            gt.dog_w(1,:)   = [base_w, base_w, base_w];

            gt.dog_amp(2,:) = [base_amp-2, base_amp, base_amp+2];
            gt.dog_w(2,:)   = [base_w+0.015, base_w, base_w-0.015];
        case 'both_same'
            gt.dog_amp(1,:) = [base_amp+2, base_amp, base_amp-2];
            gt.dog_w(1,:)   = [base_w-0.015, base_w, base_w+0.015];

            gt.dog_amp(2,:) = [base_amp-2, base_amp, base_amp+2];
            gt.dog_w(2,:)   = [base_w+0.015, base_w, base_w-0.015];
        case 'both_different'
            % Contrast has strong effect, Precision has weak effect
            gt.dog_amp(1,:) = [base_amp+4, base_amp, base_amp-3];
            gt.dog_w(1,:)   = [base_w-0.02, base_w, base_w+0.02];

            gt.dog_amp(2,:) = [base_amp-1, base_amp, base_amp+1];
            gt.dog_w(2,:)   = [base_w+0.005, base_w, base_w-0.005];
        case 'both_different_precision'
            % Precision has strong effect, Contrast has weak effect
            gt.dog_amp(1,:) = [base_amp+1, base_amp, base_amp-1];
            gt.dog_w(1,:)   = [base_w-0.005, base_w, base_w+0.005];

            gt.dog_amp(2,:) = [base_amp-4, base_amp, base_amp+3];
            gt.dog_w(2,:)   = [base_w+0.02, base_w, base_w-0.02];
    end

    gt_all(i).subj_id = subj_id;
    gt_all(i).gt = gt;
end

%% 3. Session 1: Calibration

disp('--- Starting Calibration Phase (Session 1) ---');

for i = 1:num_subjs
    subj_id = gt_all(i).subj_id;
    gt = gt_all(i).gt;

    disp(['Simulating Calibration for S' subj_id '...']);

    % Initialize experiment params
    p = struct();
    p.subj_ID = subj_id;
    p.calibration = 1;
    p.training = 0;
    p.simulation_mode = 1;
    p.display_setup = 'Macbook';

    % We need to mock 'w' to pass into init scripts
    w = struct();
    w.ppd = 40; % dummy ppd
    w.screen_width_px = 1920;
    w.screen_height_px = 1080;
    w.gray = 127; w.white = 255; w.black = 0;

    t = struct();

    % Call init scripts
    init_stimuli_params;

    % Calibration uses 8 blocks
    p.num_blocks = 8;
    p.feature_name = {'contrast', 'filter'};
    p.num_features = numel(p.feature_name);
    p.num_blocks_per_feature = p.num_blocks / p.num_features;
    p.feature_order = repmat(1:p.num_features, 1, p.num_blocks_per_feature);
    p.feature_order = Shuffle(p.feature_order);
    p.num_trials_per_lvl_per_block = 10;

    % Generate levels and events
    for n_block = 1:p.num_blocks
        level_order = BalanceFactors(p.num_trials_per_lvl_per_block, 1, 1:p.num_levels);
        if n_block == 1
            p.num_trials_per_block = length(level_order);
            p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks);
        end
        test_orientations = sample_orientation(p.orientation_min, p.orientation_max, p.num_trials_per_block);
        probe_offsets = datasample(p.probe_offsets, length(test_orientations));

        % Randomly flip probe offset sign
        signs = sign(randn(1, length(probe_offsets)));
        signs(signs == 0) = 1;
        probe_offsets = probe_offsets .* signs;

        probe_orientations = mod(test_orientations + probe_offsets', 180);
        p.trial_events(:,:,n_block) = [test_orientations, probe_orientations, level_order];
    end

    % Save ground truth into p so validate_calibration.m can access it
    p.true_params.alpha_c  = gt.contrast_alpha;
    p.true_params.beta_c   = gt.contrast_beta;
    p.true_params.alpha_fw = gt.filter_alpha;
    p.true_params.beta_fw  = gt.filter_beta;
    p.true_params.guess_rate = 0.5;
    p.true_params.lambda   = gt.lambda;

    % Simulate Responses
    [responses, correct] = simulate_responses(p, gt);

    % Save Run Info
    behav_data.responses = responses;
    behav_data.correct = correct;

    run_info.behav_data = behav_data;
    run_info.p = p;

    subj_dir = fullfile(dirs.data_dir, subj_id);
    if ~exist(subj_dir, 'dir'), mkdir(subj_dir); end

    save_filename = ['SD_Noise_Exp2_Calibration_S' subj_id '_Run1_' p.display_setup '.mat'];
    save(fullfile(subj_dir, save_filename), 'run_info');

end

disp('--- Running fit_calibration.m ---');

cd('../calibration');

calib_script = fileread('fit_calibration.m');
orig_subjs_line = regexp(calib_script, 'subj_IDs = \{[^}]+\};', 'match', 'once');
new_subjs_str = ['subj_IDs = {''' strjoin({gt_all.subj_id}, ''', ''') '''};'];

% Apply modifications to the calibration script
calib_script_mod = strrep(calib_script, 'clear all; close all; clc;', '% clear all; close all; clc;');
calib_script_mod = strrep(calib_script_mod, orig_subjs_line, new_subjs_str);
calib_script_mod = strrep(calib_script_mod, 'data_base_dir = ''../data'';', 'data_base_dir = fullfile(base_dir, ''data'');');

% Fix the save directory for figures to place them directly in the subject's simulation output folder
calib_script_mod = strrep(calib_script_mod, 'fig_dir = fullfile(script_dir, ''figures'');', 'fig_dir = fullfile(base_dir, ''experiment'', ''figures'', ''subject'', ''calibration'');');

% Ensure we close the figure at the end of the loop so we don't spawn 20+ windows
calib_script_mod = strrep(calib_script_mod, 'print(gcf, fig_filename, ''-dpdf'', ''-bestfit'');', 'print(gcf, fig_filename, ''-dpdf'', ''-bestfit''); close(gcf);');

fid = fopen('fit_calibration_sim.m', 'w');
fwrite(fid, calib_script_mod);
fclose(fid);

try
    fit_calibration_sim;
catch e
    disp('Error in calibration:');
    disp(e.message);
end

delete('fit_calibration_sim.m');
cd(base_dir); % Use absolute path to ensure we return to the simulation directory safely
disp('Calibration Phase Complete.');

%% 4. Sessions 2-5: Main Experiment

disp('--- Starting Main Experiment (Sessions 2-5) ---');

num_runs = 4; % equivalent to sessions 2-5

for i = 1:num_subjs
    subj_id = gt_all(i).subj_id;
    gt = gt_all(i).gt;

    disp(['Simulating Main Experiment for S' subj_id '...']);

    % Add path to simulate_responses
    addpath(base_dir);

    for run_num = 1:num_runs

        % Initialize experiment params
        p = struct();
        p.subj_ID = subj_id;
        p.calibration = 0;
        p.training = 0;
        p.simulation_mode = 1;
        p.display_setup = 'Macbook';

        % We need to mock 'w' to pass into init scripts
        w = struct();
        w.ppd = 40; % dummy ppd
        w.screen_width_px = 1920;
        w.screen_height_px = 1080;
        w.gray = 127; w.white = 255; w.black = 0;

        t = struct();

        % Initialize stimuli parameters (this will load S9XX_calibrated_levels.mat)
        % Note: init_stimuli_params requires dirs.data_dir
        dirs.data_dir = fullfile(base_dir, 'data');

        try
            init_stimuli_params;
        catch e
            disp(['Failed to load calibrated levels for S' subj_id '. Skipping...']);
            continue;
        end

        p.num_blocks = 6;
        p.feature_name = {'contrast', 'filter'};
        p.num_features = numel(p.feature_name);
        p.num_blocks_per_feature = p.num_blocks / p.num_features;
        p.feature_order = repmat(1:p.num_features, 1, p.num_blocks_per_feature);
        p.feature_order = Shuffle(p.feature_order);
        p.num_trials_per_lvl_per_block = 38;

        for n_block = 1:p.num_blocks
            level_order = BalanceFactors(p.num_trials_per_lvl_per_block, 1, 1:p.num_levels);
            if n_block == 1
                p.num_trials_per_block = length(level_order);
                p.trial_events = nan(p.num_trials_per_block, 3, p.num_blocks);
            end
            test_orientations = sample_orientation(p.orientation_min, p.orientation_max, p.num_trials_per_block);
            probe_offsets = datasample(p.probe_offsets, length(test_orientations));

            % Randomly flip probe offset sign
            signs = sign(randn(1, length(probe_offsets)));
            signs(signs == 0) = 1;
            probe_offsets = probe_offsets .* signs;

            probe_orientations = mod(test_orientations + probe_offsets', 180);
            p.trial_events(:,:,n_block) = [test_orientations, probe_orientations, level_order];
        end

        % Simulate Responses
        [responses, correct] = simulate_responses(p, gt);

        % Save Run Info
        behav_data.responses = responses;
        behav_data.correct = correct;

        run_info.behav_data = behav_data;
        run_info.p = p;

        subj_dir = fullfile(dirs.data_dir, subj_id);

        % Save with run number (2 to 5 conceptually, but we can just use 1 to 4 or 2 to 5)
        save_filename = ['SD_Noise_Exp2_S' subj_id '_Run' num2str(run_num) '_' p.display_setup '.mat'];
        save(fullfile(subj_dir, save_filename), 'run_info');
    end
end

disp('Simulation Complete!');

% Save ground truth so the dashboard and analysis scripts can check recovery accuracy
save('simulated_ground_truth.mat', 'gt_all', 'hypothesis');

