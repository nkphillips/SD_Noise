function [tbl_demeaned, baselines] = demeanTrialTable(tbl_trials, fit_opts, skip_at_bound_baseline)
% demeanTrialTable  Subtract a per-subject psychometric baseline mu_i from x_probe.
%
% For each subject i, fits a Gaussian psychometric with 25% lapse to all of that
% subject's trials (across all 18 cells), ignoring Delta-theta. The fitted mu_i
% is treated as the subject's overall response-bias baseline. Returns a copy of
% tbl_trials with x_probe replaced by x_probe - mu_i for each trial.
%
% When skip_at_bound_baseline is true, subjects whose baseline sigma_i saturates
% at either of its bounds are NOT demeaned (mu_i = 0 used instead). The rationale
% is that an at-bound sigma indicates the no-DoG, no-Delta-theta psychometric
% does not fit that subject's data; demeaning by an unidentified mu_i estimate
% would introduce more error than leaving the data untouched. The asymmetric
% inclusion of un-demeaned subjects in the bootstrap pool has small numerical
% effect (each un-demeaned subject contributes a small fraction of trials per
% cell, so per-cell beta absorbs only a fraction of their bias).
%
% Inputs:
%   tbl_trials               - trial table from buildTrialTableFromAllRuns.
%   fit_opts                 - struct passed to fitSubjectBaselineBias.
%   skip_at_bound_baseline   - logical (default true). If true, skip the demean
%                              for subjects whose sigma_baseline lands within
%                              1e-3 relative span of either bound.
%
% Outputs:
%   tbl_demeaned - copy of tbl_trials with x_probe replaced by x_probe - mu_i.
%   baselines    - per-subject diagnostic table:
%                    subject_id, n_trials, mu_baseline, sigma_baseline,
%                    exit_flag, fit_failed (logical), at_bound (logical),
%                    skipped (logical: true if mu_i = 0 was used instead).

    if nargin < 2 || isempty(fit_opts); fit_opts = struct(); end
    if nargin < 3 || isempty(skip_at_bound_baseline); skip_at_bound_baseline = true; end

    % Recover sigma bounds for the at-bound check
    sigma_lb = 1;   sigma_ub = 90;
    if isfield(fit_opts, 'lb') && ~isempty(fit_opts.lb), sigma_lb = fit_opts.lb(2); end
    if isfield(fit_opts, 'ub') && ~isempty(fit_opts.ub), sigma_ub = fit_opts.ub(2); end
    sigma_span = max(sigma_ub - sigma_lb, eps);
    bound_tol = 1e-3;

    subj_list = unique(tbl_trials.subject_id, 'stable');
    n_subj = numel(subj_list);

    rows = cell(n_subj, 1);
    mu_per_subj = nan(n_subj, 1);

    for i = 1:n_subj
        sid = subj_list(i);
        mask = tbl_trials.subject_id == sid;
        n_t = sum(mask);

        [mu_i, sigma_i, ef] = fitSubjectBaselineBias( ...
            tbl_trials.x_probe(mask), tbl_trials.response(mask), fit_opts);

        failed = false;
        at_bound = false;
        skipped = false;

        if ~isfinite(mu_i)
            mu_per_subj(i) = 0;
            failed = true;
            skipped = true;
            warning('demeanTrialTable:fitFailed', ...
                'Baseline fit failed for subject id %d (n_trials = %d); leaving mu = 0 (no demean).', ...
                sid, n_t);
        else
            % Check whether sigma_i landed at a bound
            if isfinite(sigma_i)
                at_bound = (abs(sigma_i - sigma_lb) / sigma_span < bound_tol) || ...
                           (abs(sigma_i - sigma_ub) / sigma_span < bound_tol);
            end

            if at_bound && skip_at_bound_baseline
                mu_per_subj(i) = 0;
                skipped = true;
                warning('demeanTrialTable:atBoundSkipped', ...
                    ['Subject id %d: baseline sigma at bound (%.2f); skipping demean ', ...
                     '(set skip_at_bound_baseline = false to override).'], sid, sigma_i);
            else
                mu_per_subj(i) = mu_i;
            end
        end

        rows{i} = {sid, n_t, mu_i, sigma_i, ef, failed, at_bound, skipped};
    end

    baselines = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'subject_id', 'n_trials', 'mu_baseline', 'sigma_baseline', ...
         'exit_flag', 'fit_failed', 'at_bound', 'skipped'});

    % Apply demean
    tbl_demeaned = tbl_trials;
    for i = 1:n_subj
        sid = subj_list(i);
        mask = tbl_demeaned.subject_id == sid;
        tbl_demeaned.x_probe(mask) = tbl_demeaned.x_probe(mask) - mu_per_subj(i);
    end
end
