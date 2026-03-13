function analyze_subject_behavior(subj_ID)
% analyze_subject_behavior(subj_ID)
% Calculates and prints overall performance, per-run performance,
% performance by feature level, reaction times, response bias, missed trials,
% and lvl-to-lvl transition performance for Experiment 2.
% Outputs to console and saves a .txt file in the subject's data directory.

if nargin < 1
    error('Please provide a subject ID (e.g., ''999'').');
end

% Define data directory
data_dir = fullfile('..', '..', 'data', subj_ID);

if ~exist(data_dir, 'dir')
    % try relative to root
    data_dir = fullfile(pwd, 'data', subj_ID);
    if ~exist(data_dir, 'dir')
        error(['Data directory for subject ' subj_ID ' not found.']);
    end
end

% Find all run files (excluding training and calibration runs)
file_pattern = sprintf('SD_Noise_Exp2_S%s_Run*.mat', subj_ID);
files = dir(fullfile(data_dir, file_pattern));

% Filter out calibration and training if any got caught (though naming usually prevents this)
keep_idx = cellfun(@isempty, regexp({files.name}, '(Training|Calibration)'));
files = files(keep_idx);

if isempty(files)
    error(['No main run data files found for subject ' subj_ID ' in ' data_dir]);
end

% Sort files to process in order
[~, sort_idx] = sort([files.datenum]);
files = files(sort_idx);

num_runs = length(files);

% Initialize variables for aggregating across runs
all_correct = [];
all_rt = [];
all_responses = [];

% Data structures for per-run and across-run metrics
run_correct = cell(num_runs, 1);
run_rt = cell(num_runs, 1);
run_missed = zeros(num_runs, 1);
run_fast_guesses = zeros(num_runs, 1);
run_num_trials = zeros(num_runs, 1);

% To hold feature-level data
feature_levels_correct = cell(2, 1); % 1: Contrast, 2: Filter Width
feature_levels_rt = cell(2, 1);

num_levels = NaN;
lvl_to_lvl_correct = cell(2, 1);
lvl_to_lvl_counts = cell(2, 1);

run_lvl_to_lvl_correct = cell(num_runs, 2);
run_lvl_to_lvl_counts = cell(num_runs, 2);

for r = 1:num_runs
    load(fullfile(data_dir, files(r).name), 'run_info');

    % Extract necessary variables
    behav = run_info.behav_data;
    p = run_info.p;

    if isnan(num_levels)
        num_levels = p.num_levels;
        for f = 1:2
            feature_levels_correct{f} = cell(num_levels, 1);
            feature_levels_rt{f} = cell(num_levels, 1);
            lvl_to_lvl_correct{f} = zeros(num_levels, num_levels);
            lvl_to_lvl_counts{f} = zeros(num_levels, num_levels);
        end
    end

    for f = 1:2
        run_lvl_to_lvl_correct{r, f} = zeros(num_levels, num_levels);
        run_lvl_to_lvl_counts{r, f} = zeros(num_levels, num_levels);
    end

    num_trials_per_block = size(behav.correct, 1);
    num_blocks = size(behav.correct, 2);

    r_correct = [];
    r_rt = [];
    r_responses = [];

    for b = 1:num_blocks
        feature = p.feature_order(b); % 1 or 2

        % Get block data
        block_correct = behav.correct(:, b);
        block_rt = behav.response_dur(:, b);
        block_response = behav.response(:, b);

        % Determine levels for each trial
        level_order_col = 3;
        levels = p.trial_events(:, level_order_col, b);

        r_correct = [r_correct; block_correct];
        r_rt = [r_rt; block_rt];
        r_responses = [r_responses; block_response];

        % Update feature-level stats
        for t = 1:num_trials_per_block
            lvl = levels(t);
            feature_levels_correct{feature}{lvl}(end+1) = block_correct(t);
            feature_levels_rt{feature}{lvl}(end+1) = block_rt(t);
        end

        % Update lvl-to-lvl transitions (within block)
        for t = 2:num_trials_per_block
            prev_lvl = levels(t-1);
            curr_lvl = levels(t);

            % Global
            lvl_to_lvl_correct{feature}(prev_lvl, curr_lvl) = lvl_to_lvl_correct{feature}(prev_lvl, curr_lvl) + block_correct(t);
            lvl_to_lvl_counts{feature}(prev_lvl, curr_lvl) = lvl_to_lvl_counts{feature}(prev_lvl, curr_lvl) + 1;

            % Run specific
            run_lvl_to_lvl_correct{r, feature}(prev_lvl, curr_lvl) = run_lvl_to_lvl_correct{r, feature}(prev_lvl, curr_lvl) + block_correct(t);
            run_lvl_to_lvl_counts{r, feature}(prev_lvl, curr_lvl) = run_lvl_to_lvl_counts{r, feature}(prev_lvl, curr_lvl) + 1;
        end
    end

    run_correct{r} = r_correct;
    run_rt{r} = r_rt;

    run_missed(r) = sum(r_responses == 0);
    run_fast_guesses(r) = sum(r_rt < 0.150 & r_responses ~= 0);
    run_num_trials(r) = length(r_responses);

    all_correct = [all_correct; r_correct];
    all_rt = [all_rt; r_rt];
    all_responses = [all_responses; r_responses];
end

%% Generate Report
report = {};
add_line = @(str) assignin('caller', 'report', [report; {str}]);

add_line(repmat('=', 1, 40));
add_line(sprintf('Subject %s Behavioral Summary', subj_ID));
add_line('Experiment: 2');
add_line(repmat('=', 1, 40));
add_line('');

% Overall Stats
overall_perf = mean(all_correct) * 100;
overall_rt_mean = mean(all_rt(all_responses ~= 0)) * 1000;
overall_rt_std = std(all_rt(all_responses ~= 0)) * 1000;

cw_resp = sum(all_responses == 2);
ccw_resp = sum(all_responses == 1);
valid_resp = cw_resp + ccw_resp;

total_missed = sum(run_missed);
total_fast = sum(run_fast_guesses);
total_trials = length(all_correct);

add_line('--- OVERALL STATISTICS ---');
add_line(sprintf('Overall Performance: %.1f%% (N=%d)', overall_perf, total_trials));
add_line(sprintf('Overall RT (Mean ± SD): %.0f ms ± %.0f ms', overall_rt_mean, overall_rt_std));
if valid_resp > 0
    add_line(sprintf('Response Bias: %.1f%% CW / %.1f%% CCW', (cw_resp/valid_resp)*100, (ccw_resp/valid_resp)*100));
else
    add_line('Response Bias: N/A');
end
add_line(sprintf('Missed Trials: %d (%.1f%%)', total_missed, (total_missed/total_trials)*100));
add_line(sprintf('Fast Guesses (<150ms): %d (%.1f%%)', total_fast, (total_fast/total_trials)*100));
add_line('');

% Per Run Performance
add_line('--- PER-RUN PERFORMANCE ---');
for r = 1:num_runs
    run_perf = mean(run_correct{r}) * 100;
    add_line(sprintf('Run %d: %.1f%% (N=%d)', r, run_perf, run_num_trials(r)));
end
add_line('');

% Performance By Feature Level
feature_names = {'Contrast', 'Filter Width'};
add_line('--- PERFORMANCE BY FEATURE LEVEL ---');
for f = 1:2
    add_line(sprintf('Feature: %s', feature_names{f}));
    for lvl = 1:num_levels
        correct_arr = feature_levels_correct{f}{lvl};
        rt_arr = feature_levels_rt{f}{lvl};
        if ~isempty(correct_arr)
            perf = mean(correct_arr) * 100;
            rt_m = mean(rt_arr) * 1000;
            add_line(sprintf('  Lvl %d: %.1f%% (Mean RT: %.0f ms, N=%d)', lvl, perf, rt_m, length(correct_arr)));
        end
    end
end
add_line('');

% Level to Level (Across Runs)
add_line('--- LEVEL -> LEVEL PERFORMANCE (Across Runs) ---');
for f = 1:2
    add_line(sprintf('Feature: %s', feature_names{f}));
    for prev_lvl = 1:num_levels
        for curr_lvl = 1:num_levels
            count = lvl_to_lvl_counts{f}(prev_lvl, curr_lvl);
            if count > 0
                perf = (lvl_to_lvl_correct{f}(prev_lvl, curr_lvl) / count) * 100;
                add_line(sprintf('  Lvl %d -> Lvl %d: %.1f%% (N=%d)', prev_lvl, curr_lvl, perf, count));
            end
        end
    end
end
add_line('');

% Level to Level (Within Runs)
add_line('--- LEVEL -> LEVEL PERFORMANCE (Within Runs) ---');
for r = 1:num_runs
    add_line(sprintf('Run %d:', r));
    for f = 1:2
        add_line(sprintf('  Feature: %s', feature_names{f}));
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                count = run_lvl_to_lvl_counts{r, f}(prev_lvl, curr_lvl);
                if count > 0
                    perf = (run_lvl_to_lvl_correct{r, f}(prev_lvl, curr_lvl) / count) * 100;
                    add_line(sprintf('    Lvl %d -> Lvl %d: %.1f%% (N=%d)', prev_lvl, curr_lvl, perf, count));
                end
            end
        end
    end
end

% Print to console
for i = 1:length(report)
    fprintf('%s\n', report{i});
end

% Save to text file
txt_filename = fullfile(data_dir, sprintf('S%s_behavioral_summary.txt', subj_ID));
fid = fopen(txt_filename, 'w');
if fid == -1
    warning('Could not open file %s for writing.', txt_filename);
else
    for i = 1:length(report)
        fprintf(fid, '%s\n', report{i});
    end
    fclose(fid);
    fprintf('\nReport saved to: %s\n', txt_filename);
end
end
