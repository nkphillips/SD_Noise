function recommend_experiment_design(session_dur_mins, t, p)
% recommend_experiment_design Calculates optimal block and trial
% distributions based on session time constraints and trial durations.
%
% Usage:
%   recommend_experiment_design()
%       Uses hardcoded default timing and features.
%
%   recommend_experiment_design(session_dur_mins)
%       Overrides target session duration (default: 90 mins).
%
%   recommend_experiment_design(session_dur_mins, t, p)
%       Uses passed-in timing and parameter structs.

% Define Defaults if not provided
if nargin < 1
    session_dur_mins = 90;
end

if nargin < 2
    % Hardcoded timing defaults from init_timing.m
    t.test_dur = 0.5;
    t.mask_dur = 0.5;
    t.delay_dur = 2.0;
    t.probe_dur = 0.5;
    t.response_dur_est = 0.5;
    t.iti_min = 1.0;
    t.iti_max = 1.5;
    t.rest_dur = 10.0;
end

if nargin < 3
    % Hardcoded feature defaults
    p.num_features = 2; % Contrast and Filter Width
    p.num_levels = 7; % Reduced from 10 to allow more trials per level
    % We do not hardcode p.num_blocks so it can be dynamically optimized
end

% 1. Calculate Average Trial Time
avg_iti = mean([t.iti_min, t.iti_max]);
trial_dur = t.test_dur + t.mask_dur + t.delay_dur + t.probe_dur + t.response_dur_est + avg_iti;

% 2. Calculate Maximum Trials for Session
buffer_mins = 10; % Give 10 mins buffer for instructions, setup, etc
active_session_secs = (session_dur_mins - buffer_mins) * 60;

% 3. Determine Block Size Constraint (Max 10 mins per block)
max_block_secs = 10 * 60;
max_trials_per_block = floor(max_block_secs / trial_dur);

% 4. Determine Replicates Per Level
% We want BalanceFactors(N, 1, 1:p.num_levels) to yield a clean block.
% Find highest integer N (replicates) that fits in our 10 min block rule.
p_num_trials_per_lvl_per_block = floor(max_trials_per_block / p.num_levels);

if p_num_trials_per_lvl_per_block < 1
    error('Trial duration is too long to even fit 1 replicate per level in a 10-minute block.');
end

trials_per_block = p_num_trials_per_lvl_per_block * p.num_levels;
block_dur_mins = (trials_per_block * trial_dur) / 60;

% 5. Determine Session Capacity (Optimizing p.num_blocks)
% A 'run' (or session) can be thought of as a set of blocks.
% We want to find the maximum number of blocks that fit in the active session time,
% constrained to be a multiple of p.num_features.

% First, find absolute max number of blocks that fit
% Block duration + rest duration (assume 1 rest after every block)
block_with_rest_dur = (trials_per_block * trial_dur) + t.rest_dur;
max_possible_blocks = floor(active_session_secs / block_with_rest_dur);

% Floor to the nearest multiple of p.num_features to ensure balance
rec_num_blocks = floor(max_possible_blocks / p.num_features) * p.num_features;

if rec_num_blocks < p.num_features
    error('Session duration too short to run even one balanced set of %d blocks.', p.num_features);
end

% For output consistency, let's treat the entire set of blocks as 1 "run"
blocks_per_run = rec_num_blocks;
num_runs = 1;
total_blocks = rec_num_blocks;

total_trials = total_blocks * trials_per_block;
total_active_mins = (total_trials * trial_dur + t.rest_dur * (total_blocks - 1)) / 60;

% Print Recommendation
fprintf('\n=== Recommended Experiment Design ===\n');
fprintf('Target Session: %d minutes (Buffer: %d mins)\n', session_dur_mins, buffer_mins);
fprintf('Average Trial Time: %.2f seconds\n\n', trial_dur);

fprintf('- Recommended Trials per Level per Block (N in BalanceFactors): %d\n', p_num_trials_per_lvl_per_block);
fprintf('- Total Trials per Block: %d\n', trials_per_block);
fprintf('- Estimated Block Duration: %.1f minutes\n\n', block_dur_mins);

fprintf('- Recommended Total Blocks: %d (multiple of %d features)\n', total_blocks, p.num_features);
fprintf('- Total Session Trials: %d\n', total_trials);
fprintf('- Estimated Total Active Task Time: %.1f minutes\n', total_active_mins);
fprintf('=====================================\n\n');

end