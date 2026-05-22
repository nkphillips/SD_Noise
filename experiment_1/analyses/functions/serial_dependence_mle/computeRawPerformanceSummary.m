function perf = computeRawPerformanceSummary(tbl_trials, subj_labels)
% computeRawPerformanceSummary  Raw 2AFC performance summaries.
%
% Correctness is defined from the probe sign: CW responses are correct when
% x_probe > 0 and CCW responses are correct when x_probe < 0. Zero-offset
% probes are excluded from accuracy because their sign is undefined.

    if nargin < 2 || isempty(subj_labels)
        subj_list = unique(tbl_trials.subject_id, 'stable');
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:numel(subj_list), 'UniformOutput', false);
    end

    required_vars = {'subject_id', 'cond_manipulation', 'cond_prev', 'cond_curr', 'x_probe', 'response'};
    missing_vars = setdiff(required_vars, tbl_trials.Properties.VariableNames);
    if ~isempty(missing_vars)
        error('computeRawPerformanceSummary:missingVariables', ...
            'tbl_trials is missing required variables: %s', strjoin(missing_vars, ', '));
    end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    manip_names = {'contrast', 'precision'};

    subject_overall = local_subjectOverall(tbl_trials, subj_list, subj_labels);
    subject_manipulation = local_subjectManipulation(tbl_trials, subj_list, subj_labels, manip_names);
    subject_level_pair = local_subjectLevelPair(tbl_trials, subj_list, subj_labels, manip_names);

    perf = struct();
    perf.subject_overall = subject_overall;
    perf.subject_manipulation = subject_manipulation;
    perf.subject_level_pair = subject_level_pair;
    perf.group_overall = local_groupSummary(subject_overall, {}, {});
    perf.group_manipulation = local_groupSummary(subject_manipulation, {'cond_manipulation'}, manip_names);
    perf.group_level_pair = local_groupSummary(subject_level_pair, ...
        {'cond_manipulation', 'cond_prev', 'cond_curr'}, manip_names);
    perf.super_subject_overall = local_superSubjectSummary(tbl_trials, {}, {}, subj_list);
    perf.super_subject_manipulation = local_superSubjectSummary(tbl_trials, {'cond_manipulation'}, manip_names, subj_list);
    perf.super_subject_level_pair = local_superSubjectSummary(tbl_trials, ...
        {'cond_manipulation', 'cond_prev', 'cond_curr'}, manip_names, subj_list);
    perf.manipulation_comparison = local_manipulationComparison(subject_manipulation);
end

function T = local_subjectOverall(tbl_trials, subj_list, subj_labels)
    rows = cell(numel(subj_list), 1);
    for s = 1:numel(subj_list)
        mask = tbl_trials.subject_id == subj_list(s);
        rows{s} = local_subjectRow(tbl_trials(mask, :), subj_list(s), subj_labels{s}, {});
    end
    T = vertcat(rows{:});
end

function T = local_subjectManipulation(tbl_trials, subj_list, subj_labels, manip_names)
    rows = cell(numel(subj_list) * numel(manip_names), 1);
    r = 0;
    for s = 1:numel(subj_list)
        for m = 1:numel(manip_names)
            r = r + 1;
            mask = tbl_trials.subject_id == subj_list(s) & tbl_trials.cond_manipulation == manip_names{m};
            rows{r} = local_subjectRow(tbl_trials(mask, :), subj_list(s), subj_labels{s}, ...
                {'cond_manipulation', string(manip_names{m})});
        end
    end
    T = vertcat(rows{:});
end

function T = local_subjectLevelPair(tbl_trials, subj_list, subj_labels, manip_names)
    rows = cell(numel(subj_list) * numel(manip_names) * 9, 1);
    r = 0;
    for s = 1:numel(subj_list)
        for m = 1:numel(manip_names)
            for prev_lvl = 1:3
                for curr_lvl = 1:3
                    r = r + 1;
                    mask = tbl_trials.subject_id == subj_list(s) & ...
                        tbl_trials.cond_manipulation == manip_names{m} & ...
                        tbl_trials.cond_prev == prev_lvl & tbl_trials.cond_curr == curr_lvl;
                    rows{r} = local_subjectRow(tbl_trials(mask, :), subj_list(s), subj_labels{s}, ...
                        {'cond_manipulation', string(manip_names{m}); ...
                         'cond_prev', prev_lvl; ...
                         'cond_curr', curr_lvl});
                end
            end
        end
    end
    T = vertcat(rows{:});
end

function row = local_subjectRow(tbl, subject_id, subject_label, extra_cols)
    stats = local_accuracyStats(tbl);
    row = table(subject_id, string(subject_label), stats.n_trials, stats.n_responses, ...
        stats.response_rate, stats.n_accuracy_trials, stats.n_correct, stats.accuracy, ...
        'VariableNames', {'subject_id', 'subject_label', 'n_trials', 'n_responses', ...
        'response_rate', 'n_accuracy_trials', 'n_correct', 'accuracy'});

    if ~isempty(extra_cols)
        for i = size(extra_cols, 1):-1:1
            row = addvars(row, extra_cols{i, 2}, 'Before', 'n_trials', ...
                'NewVariableNames', extra_cols{i, 1});
        end
    end
end

function T = local_groupSummary(subject_table, group_vars, manip_names)
    if isempty(group_vars)
        T = local_groupSummaryOne(subject_table, {}, {});
        return
    end

    rows = {};
    if isequal(group_vars, {'cond_manipulation'})
        for m = 1:numel(manip_names)
            mask = subject_table.cond_manipulation == manip_names{m};
            rows{end+1, 1} = local_groupSummaryOne(subject_table(mask, :), ...
                {'cond_manipulation'}, {string(manip_names{m})}); %#ok<AGROW>
        end
    else
        for m = 1:numel(manip_names)
            for prev_lvl = 1:3
                for curr_lvl = 1:3
                    mask = subject_table.cond_manipulation == manip_names{m} & ...
                        subject_table.cond_prev == prev_lvl & subject_table.cond_curr == curr_lvl;
                    rows{end+1, 1} = local_groupSummaryOne(subject_table(mask, :), group_vars, ...
                        {string(manip_names{m}), prev_lvl, curr_lvl}); %#ok<AGROW>
                end
            end
        end
    end
    T = vertcat(rows{:});
end

function row = local_groupSummaryOne(tbl, group_vars, group_vals)
    acc = tbl.accuracy;
    finite = isfinite(acc);
    n_subjects = sum(finite);
    mean_accuracy = mean(acc(finite), 'omitnan');
    sem_accuracy = local_sem(acc);
    [ci_lo, ci_hi] = local_meanCi(mean_accuracy, sem_accuracy, n_subjects);
    mean_response_rate = mean(tbl.response_rate(isfinite(tbl.response_rate)), 'omitnan');
    total_trials = sum(tbl.n_trials);
    total_accuracy_trials = sum(tbl.n_accuracy_trials);

    row = table(n_subjects, mean_accuracy, sem_accuracy, ci_lo, ci_hi, ...
        mean_response_rate, total_trials, total_accuracy_trials, ...
        'VariableNames', {'n_subjects', 'mean_accuracy', 'sem_accuracy', ...
        'ci95_lo', 'ci95_hi', 'mean_response_rate', 'total_trials', ...
        'total_accuracy_trials'});

    for i = numel(group_vars):-1:1
        row = addvars(row, group_vals{i}, 'Before', 'n_subjects', ...
            'NewVariableNames', group_vars{i});
    end
end

function T = local_superSubjectSummary(tbl_trials, group_vars, manip_names, subj_list)
    if isempty(group_vars)
        T = local_superSubjectSummaryOne(tbl_trials, {}, {}, subj_list);
        return
    end

    rows = {};
    if isequal(group_vars, {'cond_manipulation'})
        for m = 1:numel(manip_names)
            mask = tbl_trials.cond_manipulation == manip_names{m};
            rows{end+1, 1} = local_superSubjectSummaryOne(tbl_trials(mask, :), ...
                {'cond_manipulation'}, {string(manip_names{m})}, subj_list); %#ok<AGROW>
        end
    else
        for m = 1:numel(manip_names)
            for prev_lvl = 1:3
                for curr_lvl = 1:3
                    mask = tbl_trials.cond_manipulation == manip_names{m} & ...
                        tbl_trials.cond_prev == prev_lvl & tbl_trials.cond_curr == curr_lvl;
                    rows{end+1, 1} = local_superSubjectSummaryOne(tbl_trials(mask, :), group_vars, ...
                        {string(manip_names{m}), prev_lvl, curr_lvl}, subj_list); %#ok<AGROW>
                end
            end
        end
    end
    T = vertcat(rows{:});
end

function row = local_superSubjectSummaryOne(tbl, group_vars, group_vals, subj_list)
    stats = local_accuracyStats(tbl);
    sem_accuracy = sqrt(stats.accuracy * (1 - stats.accuracy) / max(stats.n_accuracy_trials, 1));
    if ~isfinite(stats.accuracy) || stats.n_accuracy_trials < 1
        sem_accuracy = NaN;
    end
    ci_lo = max(0, stats.accuracy - 1.96 * sem_accuracy);
    ci_hi = min(1, stats.accuracy + 1.96 * sem_accuracy);
    if ~isfinite(stats.accuracy)
        ci_lo = NaN;
        ci_hi = NaN;
    end
    n_subjects = numel(intersect(unique(tbl.subject_id), subj_list));

    row = table(n_subjects, stats.n_trials, stats.n_responses, stats.response_rate, ...
        stats.n_accuracy_trials, stats.n_correct, stats.accuracy, sem_accuracy, ci_lo, ci_hi, ...
        'VariableNames', {'n_subjects', 'n_trials', 'n_responses', 'response_rate', ...
        'n_accuracy_trials', 'n_correct', 'accuracy', 'sem_accuracy', 'ci95_lo', 'ci95_hi'});

    for i = numel(group_vars):-1:1
        row = addvars(row, group_vals{i}, 'Before', 'n_subjects', ...
            'NewVariableNames', group_vars{i});
    end
end

function T = local_manipulationComparison(subject_manipulation)
    contrast_rows = subject_manipulation(subject_manipulation.cond_manipulation == 'contrast', :);
    precision_rows = subject_manipulation(subject_manipulation.cond_manipulation == 'precision', :);
    [matched_ids, ia, ib] = intersect(contrast_rows.subject_id, precision_rows.subject_id, 'stable');

    contrast_acc = contrast_rows.accuracy(ia);
    precision_acc = precision_rows.accuracy(ib);
    valid = isfinite(contrast_acc) & isfinite(precision_acc);
    contrast_acc = contrast_acc(valid);
    precision_acc = precision_acc(valid);
    matched_ids = matched_ids(valid);
    diff_accuracy = precision_acc - contrast_acc;

    n_subjects = numel(diff_accuracy);
    mean_contrast = mean(contrast_acc, 'omitnan');
    mean_precision = mean(precision_acc, 'omitnan');
    mean_difference = mean(diff_accuracy, 'omitnan');
    sem_difference = local_sem(diff_accuracy);
    [ci_lo, ci_hi] = local_meanCi(mean_difference, sem_difference, n_subjects);

    if n_subjects > 1 && isfinite(sem_difference) && sem_difference > 0
        t_stat = mean_difference / sem_difference;
        df = n_subjects - 1;
        p_value = local_twoSidedTPValue(t_stat, df);
    else
        t_stat = NaN;
        df = max(n_subjects - 1, 0);
        p_value = NaN;
    end

    T = table(string('precision_minus_contrast'), n_subjects, numel(matched_ids), ...
        mean_contrast, mean_precision, mean_difference, sem_difference, ci_lo, ci_hi, ...
        t_stat, df, p_value, ...
        'VariableNames', {'comparison', 'n_subjects', 'n_matched_subjects', ...
        'contrast_mean_accuracy', 'precision_mean_accuracy', 'mean_difference', ...
        'sem_difference', 'ci95_lo', 'ci95_hi', 't_stat', 'df', 'p_value'});
end

function stats = local_accuracyStats(tbl)
    stats = struct();
    stats.n_trials = height(tbl);
    response_valid = isfinite(tbl.response) & (tbl.response == 0 | tbl.response == 1);
    stats.n_responses = sum(response_valid);
    if stats.n_trials > 0
        stats.response_rate = stats.n_responses / stats.n_trials;
    else
        stats.response_rate = NaN;
    end

    accuracy_valid = response_valid & isfinite(tbl.x_probe) & tbl.x_probe ~= 0;
    is_correct = (tbl.response == 1 & tbl.x_probe > 0) | ...
        (tbl.response == 0 & tbl.x_probe < 0);
    stats.n_accuracy_trials = sum(accuracy_valid);
    stats.n_correct = sum(is_correct & accuracy_valid);
    if stats.n_accuracy_trials > 0
        stats.accuracy = stats.n_correct / stats.n_accuracy_trials;
    else
        stats.accuracy = NaN;
    end
end

function sem = local_sem(x)
    x = x(isfinite(x));
    if numel(x) < 2
        sem = NaN;
    else
        sem = std(x, 0) / sqrt(numel(x));
    end
end

function [ci_lo, ci_hi] = local_meanCi(mu, sem, n)
    if n < 2 || ~isfinite(mu) || ~isfinite(sem)
        ci_lo = NaN;
        ci_hi = NaN;
        return
    end
    tcrit = local_tInv975(n - 1);
    ci_lo = mu - tcrit * sem;
    ci_hi = mu + tcrit * sem;
end

function p = local_twoSidedTPValue(t_stat, df)
    try
        p = 2 * (1 - tcdf(abs(t_stat), df));
    catch %#ok<CTCH>
        p = NaN;
    end
end

function tcrit = local_tInv975(df)
    try
        tcrit = tinv(0.975, df);
    catch %#ok<CTCH>
        tcrit = 1.96;
    end
end
