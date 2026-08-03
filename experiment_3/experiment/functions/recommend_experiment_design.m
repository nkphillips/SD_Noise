function [best_N, best_blocks] = recommend_experiment_design(t, p)
% recommend_experiment_design Calculates optimal block and trial
% distributions for the main experiment, assuming 2 fixed levels derived
% from the calibration phase.
%
% Usage:`
%   recommend_experiment_design()
%       Uses hardcoded default timing and features.
%
%   recommend_experiment_design(t, p)
%       Overrides all parameters and constraints.

%% Define Defaults if not provided

if nargin < 1 % nargin = number of arguments in input
    t.target_task_dur_mins = 60;
end

if nargin < 2
    t.max_task_dur_mins = 70;
end

if nargin < 3
    t.max_block_mins = 10;
end

if nargin < 4
    % Hardcoded timing defaults from init_timing.m
    t.test_dur = 0.5;
    t.mask_dur = 0.5;
    t.delay_dur = 2.0;
    t.probe_dur = 0.5;
    t.response_dur_est = 1;
    t.iti_min = 1.0;
    t.iti_max = 1.5;
    t.rest_dur = 10.0;
end

if nargin < 5
    % Hardcoded feature defaults
    p.num_features = 2; % Contrast and Filter Width
    p.num_levels = 2; % The main experiment always uses 2 calibrated levels
end


%%

% 1. Calculate Average Trial Time
avg_iti = mean([t.iti_min, t.iti_max]);
trial_dur = t.test_dur + t.mask_dur + t.delay_dur + t.probe_dur + t.response_dur_est + avg_iti;

% Constraints in seconds
target_task_secs = t.target_task_dur_mins * 60;
max_task_secs = t.max_task_dur_mins * 60;
max_block_secs = t.max_block_mins * 60;

% 2. Calculate optimal N (replicates per level per block)
% Each feature block tests its target feature at p.num_levels while keeping the other at baseline
trials_per_N = p.num_levels;
max_N_per_block = floor(max_block_secs / (trials_per_N * trial_dur));

if max_N_per_block < 1
    error('Max block duration is too short to fit even one replicate per level.');
end

best_N = max_N_per_block;
trials_per_block = best_N * p.num_levels;
block_dur = trials_per_block * trial_dur;

% 3. Calculate how many blocks fit in the target duration
% Block duration + rest duration (assume 1 rest after every block)
% Try to get as close to the target duration as possible without exceeding the max duration
best_blocks = 0;
best_diff_from_target = inf;

% Search for the best number of blocks
for b = p.num_features:p.num_features:100 % Check in multiples of num_features
    total_time_secs = (b * block_dur) + ((b - 1) * t.rest_dur);
    
    if total_time_secs <= max_task_secs
        diff = abs(total_time_secs - target_task_secs);
        if diff < best_diff_from_target
            best_diff_from_target = diff;
            best_blocks = b;
        end
    else
        break; % exceeded max duration
    end
end

if best_blocks < p.num_features
    error('Max task duration too short to run even one balanced set of %d blocks.', p.num_features);
end

% 5. Recalculate metrics for output
block_dur_mins = block_dur / 60;
total_trials = best_blocks * trials_per_block;
total_active_mins = (total_trials * trial_dur + t.rest_dur * (best_blocks - 1)) / 60;
total_trials_per_level = (best_blocks / p.num_features) * best_N;

fprintf('\n=== Searching for Optimal Main Experiment Design ===\n');
fprintf('Target Task Duration: %d mins (Max: %d mins)\n', t.target_task_dur_mins, t.max_task_dur_mins);
fprintf('Max Block Duration: %d mins\n', t.max_block_mins);
fprintf('Average Trial Time: %.2f seconds\n\n', trial_dur);

fprintf('=== Recommended Experiment Design ===\n');
fprintf('- Fixed Number of Levels: %d\n', p.num_levels);
fprintf('- Recommended Total Blocks: %d (multiple of %d features)\n', best_blocks, p.num_features);
fprintf('- Total Trials per Level (per feature): %d\n\n', total_trials_per_level);

fprintf('- Trials per Level per Block (N): %d\n', best_N);
fprintf('- Total Trials per Block: %d\n', trials_per_block);
fprintf('- Estimated Block Duration: %.1f minutes\n\n', block_dur_mins);

fprintf('- Total Task Trials: %d\n', total_trials);
fprintf('- Estimated Total Active Task Time: %.1f minutes\n', total_active_mins);
fprintf('=====================================\n\n');

end