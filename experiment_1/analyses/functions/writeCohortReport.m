function report_path = writeCohortReport(p, num, toggles, bootstrap, which_setup, n_back, total_trials_per_cond)
% writeCohortReport
%
% Writes a plain-text cohort audit log to analyses/reports/ documenting the
% exact subjects included in a run of SD_Noise_Analyses_And_Figures, along
% with per-subject run and trial counts and the key analysis settings. One
% file is produced per n-back iteration so subsequent figure / estimate
% auditing can recover the cohort without inspecting the saved .mat files.
%
% Inputs:
%   p                     - parameter struct (expects p.subj_IDs, p.cond_names,
%                           p.num_workers, p.num_chunks, p.guess_rate,
%                           optional p.sd_objective).
%   num                   - size struct (expects num.subjs, num.runs,
%                           num.conds, num.levels, num.delta_theta_windows).
%   toggles               - toggle struct (for parallelization, sd_objective,
%                           bootstrap flags).
%   bootstrap             - bootstrap struct (expects B, ci).
%   which_setup           - setup string used to filter data files.
%   n_back                - current n-back lag (integer).
%   total_trials_per_cond - [num.conds x num.subjs] matrix of per-subject,
%                           per-condition trial counts (as assembled in
%                           SD_Noise_Analyses_And_Figures).
%
% Output:
%   report_path - absolute path to the written file.

    % Resolve reports dir relative to this file (analyses/functions/ -> analyses/reports/)
    this_file_dir = fileparts(mfilename('fullpath'));
    reports_dir = fullfile(this_file_dir, '..', 'reports');
    if ~exist(reports_dir, 'dir')
        mkdir(reports_dir);
    end

    timestamp_human = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    timestamp_file  = datestr(now, 'yyyymmdd_HHMMSS');

    report_path = fullfile(reports_dir, ...
        ['subjects_' timestamp_file '_' num2str(n_back) '_back.txt']);

    fid = fopen(report_path, 'w');
    if fid == -1
        warning('writeCohortReport:cannotOpen', ...
            'Could not open report file for writing: %s', report_path);
        return
    end
    closer = onCleanup(@() fclose(fid));

    fprintf(fid, 'Subject Cohort Log (Super Subject)\n');
    fprintf(fid, 'Generated: %s\n', timestamp_human);
    fprintf(fid, 'N-Back:    %d\n', n_back);
    fprintf(fid, '====================================\n\n');

    % Run configuration
    fprintf(fid, 'Run configuration\n');
    fprintf(fid, '  which_setup       : %s\n', which_setup);
    if isfield(toggles, 'sd_objective')
        fprintf(fid, '  sd_objective      : %s\n', toggles.sd_objective);
    end
    if isfield(p, 'guess_rate')
        fprintf(fid, '  guess_rate        : %g\n', p.guess_rate);
    end
    if isfield(bootstrap, 'B_rb_perf') && ~isempty(bootstrap.B_rb_perf)
        fprintf(fid, '  bootstrap.B_rb_perf: %d\n', bootstrap.B_rb_perf);
    elseif isfield(bootstrap, 'B')
        fprintf(fid, '  bootstrap.B (legacy): %d\n', bootstrap.B);
    end
    if isfield(bootstrap, 'ci') && numel(bootstrap.ci) == 2
        fprintf(fid, '  bootstrap.ci      : [%g %g]\n', bootstrap.ci(1), bootstrap.ci(2));
    end
    if isfield(toggles, 'bootstrap_rb_perf')
        fprintf(fid, '  bootstrap_rb_perf : %d\n', toggles.bootstrap_rb_perf);
    end
    if isfield(toggles, 'bootstrap_sd_cluster')
        fprintf(fid, '  bootstrap_sd_cluster: %d\n', toggles.bootstrap_sd_cluster);
    end
    if isfield(bootstrap, 'B_subject_cluster_sd')
        fprintf(fid, '  B_subject_cluster_sd: %d\n', bootstrap.B_subject_cluster_sd);
    end
    if isfield(toggles, 'parallelization')
        fprintf(fid, '  parallelization   : %d\n', toggles.parallelization);
    end
    if isfield(p, 'num_workers')
        fprintf(fid, '  num_workers       : %d\n', p.num_workers);
    end
    if isfield(p, 'num_chunks')
        fprintf(fid, '  num_chunks        : %d\n', p.num_chunks);
    end
    fprintf(fid, '\n');

    % Cohort summary
    fprintf(fid, 'Cohort summary\n');
    fprintf(fid, '  num.subjs         : %d\n', num.subjs);
    fprintf(fid, '  num.conds         : %d\n', num.conds);
    fprintf(fid, '  num.levels        : %d\n', num.levels);
    if isfield(num, 'delta_theta_windows')
        fprintf(fid, '  num.delta_theta_windows : %d\n', num.delta_theta_windows);
    end
    fprintf(fid, '\n');

    % Ordered subject list
    fprintf(fid, 'Subjects (in analysis order)\n');
    for subj = 1:num.subjs
        if isfield(num, 'runs') && subj <= numel(num.runs)
            n_runs = num.runs(subj);
        else
            n_runs = NaN;
        end
        if ~isempty(total_trials_per_cond) && subj <= size(total_trials_per_cond, 2)
            subj_trials = sum(total_trials_per_cond(:, subj), 'omitnan');
        else
            subj_trials = NaN;
        end
        fprintf(fid, '  %02d. S%s  runs=%g  trials=%g\n', ...
            subj, p.subj_IDs{subj}, n_runs, subj_trials);
    end

    % Grand totals
    if ~isempty(total_trials_per_cond)
        grand_total = sum(total_trials_per_cond(:), 'omitnan');
        fprintf(fid, '\n  Grand total trials : %d\n', grand_total);
    end

end
