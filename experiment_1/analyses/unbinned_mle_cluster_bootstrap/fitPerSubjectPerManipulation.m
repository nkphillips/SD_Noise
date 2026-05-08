function T = fitPerSubjectPerManipulation(tbl_trials, fit_opts, subj_labels)
% fitPerSubjectPerManipulation  Fit one DoG (A, w, sigma, beta) per subject per
% manipulation, pooling across the 9 prev x curr cells within that manipulation.
% Provides the S&S 2022 Fig 1H analog at the per-subject-per-manipulation grain.
%
% Inputs:
%   tbl_trials   - trial table from buildTrialTableFromAllRuns; subject_id is the
%                  numeric row index into subj_labels.
%   fit_opts     - struct passed to fitConditionMLE (e.g. .lb, .ub, .x0, .guess_rate).
%   subj_labels  - optional cell array of subject ID strings; if omitted, generated
%                  as 'S01', 'S02', ...
%
% Output:
%   T - long-form table:
%       subject_id (numeric), subject_label (string),
%       manipulation (string), n_trials,
%       A, w, fwhm, sigma, beta, exit_flag.

    if nargin < 2; fit_opts = struct(); end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    n_subj = numel(subj_list);

    if nargin < 3 || isempty(subj_labels)
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:n_subj, 'UniformOutput', false);
    end

    % Reproduce bounds used by fitConditionMLE so we can flag at-bound fits.
    fwhm_min_deg = 8;
    fwhm_max_deg = 120;
    w_lb_def = (2 * sqrt(log(2))) / fwhm_max_deg;
    w_ub_def = (2 * sqrt(log(2))) / fwhm_min_deg;
    lb = [-30; w_lb_def;  1; -10];
    ub = [ 30; w_ub_def; 90;  10];
    if isfield(fit_opts, 'lb') && ~isempty(fit_opts.lb), lb = fit_opts.lb(:); end
    if isfield(fit_opts, 'ub') && ~isempty(fit_opts.ub), ub = fit_opts.ub(:); end
    span = max(abs(ub - lb), eps);
    bound_tol = 1e-3;   % slightly looser than the bootstrap pipeline's 1e-4

    manips = {'contrast', 'precision'};

    rows = cell(n_subj * 2, 1);
    idx = 0;

    for s = 1:n_subj
        for m = 1:numel(manips)
            sid = subj_list(s);
            mask = (tbl_trials.subject_id == sid) & ...
                   (tbl_trials.cond_manipulation == manips{m});
            n_trials = sum(mask);

            A = NaN; w = NaN; sg = NaN; b = NaN; ef = NaN; fwhm = NaN;
            at_bound_any = false;
            at_bound_count = 0;

            if n_trials >= 50
                [pf, ef_local] = fitConditionMLE( ...
                    tbl_trials.delta_theta(mask), ...
                    tbl_trials.x_probe(mask), ...
                    tbl_trials.response(mask), fit_opts);

                if all(isfinite(pf))
                    A = pf(1); w = pf(2); sg = pf(3); b = pf(4);
                    fwhm = unbinnedWtoFwhm(w);

                    pf_col = pf(:);
                    bound_hits = (abs(pf_col - lb) ./ span < bound_tol) | ...
                                 (abs(pf_col - ub) ./ span < bound_tol);
                    at_bound_any = any(bound_hits);
                    at_bound_count = sum(bound_hits);
                end
                ef = ef_local;
            end

            idx = idx + 1;
            rows{idx} = {sid, string(subj_labels{s}), string(manips{m}), n_trials, ...
                         A, w, fwhm, sg, b, ef, at_bound_any, at_bound_count};
        end
    end

    T = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'subject_id', 'subject_label', 'manipulation', 'n_trials', ...
         'A', 'w', 'fwhm', 'sigma', 'beta', 'exit_flag', ...
         'at_bound_any', 'at_bound_count'});
end
