% summarizePerformance
% Creates a structure summarizing performance within and across conditions
% at subject-level (ind), group-level (grp), and super-subject level (all).
%
% Inputs:
%   delta_theta_windows: struct with performance data at all levels
%   p: struct with condition names and other parameters
%   num: struct with number of subjects, levels, conditions, etc.
%
% Output:
%   behav_data: struct with fields:
%       .ind: subject-level data
%           .overall: overall performance (collapsed across all conditions)
%           .by_condition: performance by condition (Contrast, Precision)
%       .grp: group-level data (same structure)
%       .all: super-subject level data (same structure)
%   Each level contains:
%       .overall: struct with fields
%           .percent_correct: scalar or vector (for ind, one per subject)
%           .percent_correct_sem: scalar or vector
%           .percent_ccw: scalar or vector
%           .percent_ccw_sem: scalar or vector
%       .by_condition: struct array (one per condition) with same fields

function behav_data = summarizePerformance(delta_theta_windows, p, num)

behav_data = struct('ind', [], 'grp', [], 'all', []);

%% Individual subject level
behav_data.ind = struct('overall', [], 'by_condition', []);

for subj = 1:num.subjs
    % Compute overall performance from raw data (collapse across all conditions and levels)
    all_responses = [];
    all_probe_offsets = [];

    for cond = 1:num.conds
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for i_window = 1:num.delta_theta_windows
                    resp = delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, i_window};
                    probe = delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                    if ~isempty(resp) && ~isempty(probe)
                        n = min(length(resp), length(probe));
                        all_responses = [all_responses; resp(1:n)];
                        all_probe_offsets = [all_probe_offsets; probe(1:n)];
                    end
                end
            end
        end
    end

    % Compute overall performance
    if ~isempty(all_responses) && ~isempty(all_probe_offsets)
        n = min(length(all_responses), length(all_probe_offsets));
        all_responses = all_responses(1:n);
        all_probe_offsets = all_probe_offsets(1:n);

        is_correct = (all_responses == (all_probe_offsets > 0));
        p_correct = mean(is_correct, 'omitnan');
        p_ccw = mean(1 - all_responses, 'omitnan');

        behav_data.ind.overall.percent_correct(subj) = 100 * p_correct;
        behav_data.ind.overall.percent_ccw(subj) = 100 * p_ccw;

        % SEM: binomial standard error
        n_trials = length(all_responses);
        if n_trials > 0
            behav_data.ind.overall.percent_correct_sem(subj) = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);
            behav_data.ind.overall.percent_ccw_sem(subj) = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
        else
            behav_data.ind.overall.percent_correct_sem(subj) = 0;
            behav_data.ind.overall.percent_ccw_sem(subj) = 0;
        end
    else
        behav_data.ind.overall.percent_correct(subj) = NaN;
        behav_data.ind.overall.percent_correct_sem(subj) = NaN;
        behav_data.ind.overall.percent_ccw(subj) = NaN;
        behav_data.ind.overall.percent_ccw_sem(subj) = NaN;
    end

    % Performance by condition
    for cond = 1:num.conds
        cond_responses = [];
        cond_probe_offsets = [];

        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for i_window = 1:num.delta_theta_windows
                    resp = delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, i_window};
                    probe = delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                    if ~isempty(resp) && ~isempty(probe)
                        n = min(length(resp), length(probe));
                        cond_responses = [cond_responses; resp(1:n)];
                        cond_probe_offsets = [cond_probe_offsets; probe(1:n)];
                    end
                end
            end
        end

        if ~isempty(cond_responses) && ~isempty(cond_probe_offsets)
            n = min(length(cond_responses), length(cond_probe_offsets));
            cond_responses = cond_responses(1:n);
            cond_probe_offsets = cond_probe_offsets(1:n);

            is_correct = (cond_responses == (cond_probe_offsets > 0));
            p_correct = mean(is_correct, 'omitnan');
            p_ccw = mean(1 - cond_responses, 'omitnan');

            behav_data.ind.by_condition(cond).percent_correct(subj) = 100 * p_correct;
            behav_data.ind.by_condition(cond).percent_ccw(subj) = 100 * p_ccw;

            n_trials = length(cond_responses);
            if n_trials > 0
                behav_data.ind.by_condition(cond).percent_correct_sem(subj) = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);
                behav_data.ind.by_condition(cond).percent_ccw_sem(subj) = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
            else
                behav_data.ind.by_condition(cond).percent_correct_sem(subj) = 0;
                behav_data.ind.by_condition(cond).percent_ccw_sem(subj) = 0;
            end
        else
            behav_data.ind.by_condition(cond).percent_correct(subj) = NaN;
            behav_data.ind.by_condition(cond).percent_correct_sem(subj) = NaN;
            behav_data.ind.by_condition(cond).percent_ccw(subj) = NaN;
            behav_data.ind.by_condition(cond).percent_ccw_sem(subj) = NaN;
        end
    end
end

%% Group level (averages across subjects)
behav_data.grp = struct('overall', [], 'by_condition', []);

% Overall: average across subjects
valid_subjs = ~isnan(behav_data.ind.overall.percent_correct);
if any(valid_subjs)
    behav_data.grp.overall.percent_correct = mean(behav_data.ind.overall.percent_correct(valid_subjs), 'omitnan');
    valid_subjs_pCW = ~isnan(behav_data.ind.overall.percent_ccw);
    if any(valid_subjs_pCW)
        behav_data.grp.overall.percent_ccw = mean(behav_data.ind.overall.percent_ccw(valid_subjs_pCW), 'omitnan');
    else
        behav_data.grp.overall.percent_ccw = NaN;
    end

    % SEM: standard error of the mean across subjects
    if sum(valid_subjs) > 1
        behav_data.grp.overall.percent_correct_sem = std(behav_data.ind.overall.percent_correct(valid_subjs), 'omitnan') / sqrt(sum(valid_subjs));
    else
        behav_data.grp.overall.percent_correct_sem = 0;
    end

    if sum(valid_subjs_pCW) > 1
        behav_data.grp.overall.percent_ccw_sem = std(behav_data.ind.overall.percent_ccw(valid_subjs_pCW), 'omitnan') / sqrt(sum(valid_subjs_pCW));
    else
        behav_data.grp.overall.percent_ccw_sem = 0;
    end
else
    behav_data.grp.overall.percent_correct = NaN;
    behav_data.grp.overall.percent_correct_sem = NaN;
    behav_data.grp.overall.percent_ccw = NaN;
    behav_data.grp.overall.percent_ccw_sem = NaN;
end

% By condition: average across subjects
for cond = 1:num.conds
    valid_subjs = ~isnan(behav_data.ind.by_condition(cond).percent_correct);
    if any(valid_subjs)
        behav_data.grp.by_condition(cond).percent_correct = mean(behav_data.ind.by_condition(cond).percent_correct(valid_subjs), 'omitnan');
        valid_subjs_pCW = ~isnan(behav_data.ind.by_condition(cond).percent_ccw);
        if any(valid_subjs_pCW)
            behav_data.grp.by_condition(cond).percent_ccw = mean(behav_data.ind.by_condition(cond).percent_ccw(valid_subjs_pCW), 'omitnan');
        else
            behav_data.grp.by_condition(cond).percent_ccw = NaN;
        end

        if sum(valid_subjs) > 1
            behav_data.grp.by_condition(cond).percent_correct_sem = std(behav_data.ind.by_condition(cond).percent_correct(valid_subjs), 'omitnan') / sqrt(sum(valid_subjs));
        else
            behav_data.grp.by_condition(cond).percent_correct_sem = 0;
        end

        if sum(valid_subjs_pCW) > 1
            behav_data.grp.by_condition(cond).percent_ccw_sem = std(behav_data.ind.by_condition(cond).percent_ccw(valid_subjs_pCW), 'omitnan') / sqrt(sum(valid_subjs_pCW));
        else
            behav_data.grp.by_condition(cond).percent_ccw_sem = 0;
        end
    else
        behav_data.grp.by_condition(cond).percent_correct = NaN;
        behav_data.grp.by_condition(cond).percent_correct_sem = NaN;
        behav_data.grp.by_condition(cond).percent_ccw = NaN;
        behav_data.grp.by_condition(cond).percent_ccw_sem = NaN;
    end
end

%% Super-subject level (from aggregated data)
behav_data.all = struct('overall', [], 'by_condition', []);

% Overall: collapse across all conditions and levels from raw data
all_responses = [];
all_probe_offsets = [];

for cond = 1:num.conds
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for i_window = 1:num.delta_theta_windows
                resp = delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window};
                probe = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                if ~isempty(resp) && ~isempty(probe)
                    n = min(length(resp), length(probe));
                    all_responses = [all_responses; resp(1:n)];
                    all_probe_offsets = [all_probe_offsets; probe(1:n)];
                end
            end
        end
    end
end

if ~isempty(all_responses) && ~isempty(all_probe_offsets)
    n = min(length(all_responses), length(all_probe_offsets));
    all_responses = all_responses(1:n);
    all_probe_offsets = all_probe_offsets(1:n);

    is_correct = (all_responses == (all_probe_offsets > 0));
    p_correct = mean(is_correct, 'omitnan');
    p_ccw = mean(1 - all_responses, 'omitnan');

    behav_data.all.overall.percent_correct = 100 * p_correct;
    behav_data.all.overall.percent_ccw = 100 * p_ccw;

    n_trials = length(all_responses);
    if n_trials > 0
        behav_data.all.overall.percent_correct_sem = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);
        behav_data.all.overall.percent_ccw_sem = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
    else
        behav_data.all.overall.percent_correct_sem = 0;
        behav_data.all.overall.percent_ccw_sem = 0;
    end
else
    behav_data.all.overall.percent_correct = NaN;
    behav_data.all.overall.percent_correct_sem = NaN;
    behav_data.all.overall.percent_ccw = NaN;
    behav_data.all.overall.percent_ccw_sem = NaN;
end

% By condition
for cond = 1:num.conds
    cond_responses = [];
    cond_probe_offsets = [];

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for i_window = 1:num.delta_theta_windows
                resp = delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, i_window};
                probe = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, i_window};

                if ~isempty(resp) && ~isempty(probe)
                    n = min(length(resp), length(probe));
                    cond_responses = [cond_responses; resp(1:n)];
                    cond_probe_offsets = [cond_probe_offsets; probe(1:n)];
                end
            end
        end
    end

    if ~isempty(cond_responses) && ~isempty(cond_probe_offsets)
        n = min(length(cond_responses), length(cond_probe_offsets));
        cond_responses = cond_responses(1:n);
        cond_probe_offsets = cond_probe_offsets(1:n);

        is_correct = (cond_responses == (cond_probe_offsets > 0));
        p_correct = mean(is_correct, 'omitnan');
        p_ccw = mean(1 - cond_responses, 'omitnan');

        behav_data.all.by_condition(cond).percent_correct = 100 * p_correct;
        behav_data.all.by_condition(cond).percent_ccw = 100 * p_ccw;

        n_trials = length(cond_responses);
        if n_trials > 0
            behav_data.all.by_condition(cond).percent_correct_sem = 100 * sqrt(p_correct * (1 - p_correct) / n_trials);
            behav_data.all.by_condition(cond).percent_ccw_sem = 100 * sqrt(p_ccw * (1 - p_ccw) / n_trials);
        else
            behav_data.all.by_condition(cond).percent_correct_sem = 0;
            behav_data.all.by_condition(cond).percent_ccw_sem = 0;
        end
    else
        behav_data.all.by_condition(cond).percent_correct = NaN;
        behav_data.all.by_condition(cond).percent_correct_sem = NaN;
        behav_data.all.by_condition(cond).percent_ccw = NaN;
        behav_data.all.by_condition(cond).percent_ccw_sem = NaN;
    end
end

end

