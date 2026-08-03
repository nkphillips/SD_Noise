function recommend_calibration_design(target_trials_per_level, task_dur_mins, max_block_mins, t, p)
% recommend_calibration_design Calculates optimal block and trial
% distributions to reach a target number of trials per level while fitting
% within task duration and block duration constraints.
%
% Usage:
%   recommend_calibration_design()
%       Uses hardcoded default targets and timing.
%
%   recommend_calibration_design(target_trials_per_level)
%       Overrides target trials per level per feature (default: 40).
%
%   recommend_calibration_design(target_trials_per_level, task_dur_mins, max_block_mins, t, p)
%       Overrides all parameters and constraints.

% Define Defaults if not provided
if nargin < 1
    target_trials_per_level = 40; % 40-50 is standard for psychometric MLE fits
end

if nargin < 2
    task_dur_mins = 60;
end

if nargin < 3
    max_block_mins = 10;
end

if nargin < 4
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

if nargin < 5
    % Hardcoded feature defaults
    p.num_features = 2; % Contrast and Filter Width
end

% 1. Calculate Average Trial Time
avg_iti = mean([t.iti_min, t.iti_max]);
trial_dur = t.test_dur + t.mask_dur + t.delay_dur + t.probe_dur + t.response_dur_est + avg_iti;

% Constraints in seconds
active_task_secs = task_dur_mins * 60;
max_block_secs = max_block_mins * 60;

fprintf('\n=== Searching for Optimal Design ===\n');
fprintf('Target Trials Per Level: %d\n', target_trials_per_level);
fprintf('Max Task Duration: %d mins\n', task_dur_mins);
fprintf('Max Block Duration: %d mins\n', max_block_mins);
fprintf('Average Trial Time: %.2f seconds\n\n', trial_dur);

best_levels = 0;
best_N = 0;
best_blocks = 0;
best_actual_trials = 0;

% Search from high level counts down to 2
for num_levels = 15:-1:2
    
    max_N = floor(max_block_secs / (num_levels * trial_dur));
    
    if max_N < 1
        continue;
    end
    
    for N = max_N:-1:1
        trials_per_block = N * num_levels;
        block_dur = trials_per_block * trial_dur;
        
        max_possible_blocks = floor((active_task_secs + t.rest_dur) / (block_dur + t.rest_dur));
        
        % Ensure balanced number of blocks for all features
        B = floor(max_possible_blocks / p.num_features) * p.num_features;
        
        if B < p.num_features
            continue;
        end
        
        % Total trials for a given level of a specific feature
        actual_trials = N * (B / p.num_features);
        
        if actual_trials >= target_trials_per_level
            if num_levels > best_levels
                best_levels = num_levels;
                best_N = N;
                best_blocks = B;
                best_actual_trials = actual_trials;
            elseif num_levels == best_levels
                % If there is a tie in num_levels, prefer the configuration
                % that gets closest to the target without going vastly over.
                if actual_trials < best_actual_trials
                    best_levels = num_levels;
                    best_N = N;
                    best_blocks = B;
                    best_actual_trials = actual_trials;
                end
            end
        end
    end
end

if best_levels == 0
    error('Could not find any design that satisfies the constraints! Try lowering the target trials per level or increasing the task duration.');
end

% Recalculate metrics for best fit
trials_per_block = best_N * best_levels;
block_dur_mins = (trials_per_block * trial_dur) / 60;
total_trials = best_blocks * trials_per_block;
total_active_mins = (total_trials * trial_dur + t.rest_dur * (best_blocks - 1)) / 60;

fprintf('=== Recommended Calibration Design ===\n');
fprintf('- Recommended Number of Levels: %d\n', best_levels);
fprintf('- Recommended Total Blocks: %d (multiple of %d features)\n', best_blocks, p.num_features);
fprintf('- Total Trials per Level (per feature): %d (Target was %d)\n\n', best_actual_trials, target_trials_per_level);

fprintf('- Trials per Level per Block (N): %d\n', best_N);
fprintf('- Total Trials per Block: %d\n', trials_per_block);
fprintf('- Estimated Block Duration: %.1f minutes\n\n', block_dur_mins);

fprintf('- Total Task Trials: %d\n', total_trials);
fprintf('- Estimated Total Active Task Time: %.1f minutes\n', total_active_mins);
fprintf('=====================================\n\n');

end