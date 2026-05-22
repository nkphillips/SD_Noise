function [group_params, subject_params, prior_mu, subject_ids, exit_flags] = fitConditionHierarchical(tbl_cond, fit_opts)
% fitConditionHierarchical  Empirical-Bayes MAP fit for one condition table.
%
% First estimates an empirical prior from all condition trials using standard
% pooled MLE, then fits each subject with an L2 penalty around that prior.
%
% Outputs:
%   group_params   - 1x4 mean of subject MAP params [A, w, sigma, beta]
%   subject_params - nSubjects x 4 MAP params
%   prior_mu       - 1x4 pooled MLE params used as empirical prior
%   subject_ids    - stable unique subject IDs from tbl_cond
%   exit_flags     - nSubjects x 1 fmincon exit flags for subject MAP fits

    if nargin < 2 || isempty(fit_opts)
        fit_opts = struct();
    end

    group_params = nan(1, 4);
    subject_params = nan(0, 4);
    prior_mu = nan(1, 4);
    subject_ids = [];
    exit_flags = nan(0, 1);

    if isempty(tbl_cond) || height(tbl_cond) < 1
        return
    end

    required_vars = {'subject_id', 'delta_theta', 'x_probe', 'response'};
    missing_vars = setdiff(required_vars, tbl_cond.Properties.VariableNames);
    if ~isempty(missing_vars)
        error('fitConditionHierarchical:missingVariables', ...
            'tbl_cond is missing required variable(s): %s.', strjoin(missing_vars, ', '));
    end

    subject_ids = unique(tbl_cond.subject_id, 'stable');
    n_subj = numel(subject_ids);
    subject_params = nan(n_subj, 4);
    exit_flags = nan(n_subj, 1);

    pooled_opts = fit_opts;
    pooled_opts.use_map = false;
    [prior_mu, prior_exit] = fitConditionMLE( ...
        tbl_cond.delta_theta, tbl_cond.x_probe, tbl_cond.response, pooled_opts);

    if ~all(isfinite(prior_mu)) || prior_exit < 1
        return
    end

    map_opts = fit_opts;
    map_opts.use_map = true;
    map_opts.prior_means = prior_mu;

    for s = 1:n_subj
        mask = tbl_cond.subject_id == subject_ids(s);
        [subject_params(s, :), exit_flags(s)] = fitConditionMLE( ...
            tbl_cond.delta_theta(mask), tbl_cond.x_probe(mask), tbl_cond.response(mask), map_opts);
    end

    valid_subjects = all(isfinite(subject_params), 2);
    if any(valid_subjects)
        group_params = mean(subject_params(valid_subjects, :), 1);
    end
end
