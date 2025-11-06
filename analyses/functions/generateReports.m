function generateReports(trial_counts, total_trials_per_cond, p, num, sd, reportMeta)
% generateReports
%
% Create text reports summarizing super-subject trial counts and DoG
% estimates with R^2 for every condition. Files are written to
% analyses/reports with a timestamp in both filename and header.

    % Resolve reports directory relative to this file location
    this_file_dir = fileparts(mfilename('fullpath'));
    reports_dir = fullfile(this_file_dir, '..', 'reports');
    if ~exist(reports_dir, 'dir')
        mkdir(reports_dir);
    end

    % Timestamp strings
    timestamp_human = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    timestamp_file = datestr(now, 'yyyymmdd_HHMMSS');

    % Suffix for filenames
    n_back_suffix = ['_' num2str(reportMeta.n_back) '_back'];

    %% Report 1: Trial counts (super subject)
    trial_report_path = fullfile(reports_dir, ['trial_counts_super_subject_' timestamp_file n_back_suffix '.txt']);
    fid = fopen(trial_report_path, 'w');
    fprintf(fid, 'Trial Counts (Super Subject)\n');
    fprintf(fid, 'Generated: %s\n', timestamp_human);
    fprintf(fid, 'N-Back: %d\n', reportMeta.n_back);
    fprintf(fid, '====================================\n\n');

    % Grand totals
    cond_totals = sum(total_trials_per_cond, 2, 'omitnan');
    grand_total = sum(cond_totals, 'omitnan');
    fprintf(fid, 'Grand total trials across all conditions: %d\n\n', grand_total);

    % Per-condition totals
    fprintf(fid, 'Trials per condition (across all subjects):\n');
    for cond = 1:num.conds
        curr_total = cond_totals(cond);
        fprintf(fid, '  %s: %d\n', p.cond_names{cond}, curr_total);
    end
    fprintf(fid, '\n');

    % Level-combination totals per condition
    fprintf(fid, 'Trials per level combination (across all subjects):\n');
    for cond = 1:num.conds
        fprintf(fid, '  %s:\n', p.cond_names{cond});
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                combo_total = sum(trial_counts(prev_lvl, curr_lvl, cond, :), 4, 'omitnan');
                if cond == 1
                    prev_label = p.contrast{prev_lvl};
                    curr_label = p.contrast{curr_lvl};
                else
                    prev_label = p.precision{prev_lvl};
                    curr_label = p.precision{curr_lvl};
                end
                fprintf(fid, '    %s -> %s: %d\n', prev_label, curr_label, combo_total);
            end
        end
    end
    fclose(fid);

    %% Report 2: DoG estimates and R^2 for every condition
    dog_report_path = fullfile(reports_dir, ['dog_estimates_r2_' timestamp_file n_back_suffix '.txt']);
    fid = fopen(dog_report_path, 'w');
    fprintf(fid, 'DoG Estimates and R^2 (Super Subject)\n');
    fprintf(fid, 'Generated: %s\n', timestamp_human);
    fprintf(fid, 'N-Back: %d\n', reportMeta.n_back);
    fprintf(fid, '====================================\n\n');

    % Note whether sigma is present (NLL objective)
    has_sigma = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'nll');

    for cond = 1:num.conds
        fprintf(fid, '%s\n', p.cond_names{cond});
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, :));
                r2_val = sd.all.r2(prev_lvl, curr_lvl, cond);

                if cond == 1
                    prev_label = p.contrast{prev_lvl};
                    curr_label = p.contrast{curr_lvl};
                else
                    prev_label = p.precision{prev_lvl};
                    curr_label = p.precision{curr_lvl};
                end

                if all(isnan(params))
                    fprintf(fid, '  %s -> %s: No data\n', prev_label, curr_label);
                    continue;
                end

                amp = params(1);
                w = params(2);
                base = params(3);
                if has_sigma
                    sigma = params(4);
                    fprintf(fid, '  %s -> %s: A = %.3f, w = %.3f, b = %.3f, sigma = %.3f, R^2 = %.3f\n', ...
                        prev_label, curr_label, amp, w, base, sigma, r2_val);
                else
                    fprintf(fid, '  %s -> %s: A = %.3f, w = %.3f, b = %.3f, R^2 = %.3f\n', ...
                        prev_label, curr_label, amp, w, base, r2_val);
                end
            end
        end
        fprintf(fid, '\n');
    end
    fclose(fid);

    % Console notice
    disp(['✓ Reports written to: ' reports_dir]);

    %% Report 3: Duration summaries
    if nargin >= 6 && ~isempty(reportMeta)
        durations_report_path = fullfile(reports_dir, ['analysis_durations_' timestamp_file n_back_suffix '.txt']);
        fid = fopen(durations_report_path, 'w');
        fprintf(fid, 'Analysis Duration Summary\n');
        fprintf(fid, 'Generated: %s\n', timestamp_human);
        fprintf(fid, 'N-Back: %d\n', reportMeta.n_back);
        fprintf(fid, '====================================\n\n');

        % Helper inline to print minutes and seconds
        fmt_duration = @(sec) sprintf('~%.1f min (%.1f s)', sec/60, sec);

        if isfield(reportMeta, 'overall_duration') && ~isempty(reportMeta.overall_duration)
            fprintf(fid, 'Overall analysis time: %s\n', fmt_duration(reportMeta.overall_duration));
        end
        if isfield(reportMeta, 'rb_duration') && ~isempty(reportMeta.rb_duration)
            if isfield(reportMeta, 'num_tasks') && ~isempty(reportMeta.num_tasks)
                fprintf(fid, 'Response bias: %d tasks, %s\n', reportMeta.num_tasks, fmt_duration(reportMeta.rb_duration));
            else
                fprintf(fid, 'Response bias: %s\n', fmt_duration(reportMeta.rb_duration));
            end
        end
        if isfield(reportMeta, 'sd_duration') && ~isempty(reportMeta.sd_duration)
            if isfield(reportMeta, 'num_sd_tasks') && ~isempty(reportMeta.num_sd_tasks)
                fprintf(fid, 'Serial dependence: %d conditions, %s\n', reportMeta.num_sd_tasks, fmt_duration(reportMeta.sd_duration));
            else
                fprintf(fid, 'Serial dependence: %s\n', fmt_duration(reportMeta.sd_duration));
            end
        end
        if isfield(reportMeta, 'bs_duration') && ~isempty(reportMeta.bs_duration) && isfinite(reportMeta.bs_duration)
            if isfield(reportMeta, 'bootstrap') && isfield(reportMeta.bootstrap, 'B')
                fprintf(fid, 'Bootstrapping: B = %d, %s\n', reportMeta.bootstrap.B, fmt_duration(reportMeta.bs_duration));
            else
                fprintf(fid, 'Bootstrapping: %s\n', fmt_duration(reportMeta.bs_duration));
            end
        end

        % Optional extra context
        if isfield(reportMeta, 'parallelization') && ~isempty(reportMeta.parallelization)
            fprintf(fid, '\nParallelization: %s\n', ternary(reportMeta.parallelization, 'Enabled', 'Disabled'));
        end
        if isfield(reportMeta, 'num_chunks') && ~isempty(reportMeta.num_chunks)
            fprintf(fid, 'Chunks used: %d\n', reportMeta.num_chunks);
        end

        fclose(fid);
    end
end

function out = ternary(cond, a, b)
    if cond
        out = a;
    else
        out = b;
    end
end


