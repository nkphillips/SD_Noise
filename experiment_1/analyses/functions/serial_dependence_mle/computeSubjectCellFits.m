function T = computeSubjectCellFits(tbl_trials, fit_opts, subj_labels, varargin)
% computeSubjectCellFits  Fit one DoG per subject x manipulation x prev x curr cell.

    ip = inputParser;
    addParameter(ip, 'min_trials', 10, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(ip, 'bound_tol', 1e-3, @(x) isnumeric(x) && isscalar(x) && x > 0);
    parse(ip, varargin{:});

    if nargin < 2 || isempty(fit_opts)
        fit_opts = struct();
    end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    n_subj = numel(subj_list);
    if nargin < 3 || isempty(subj_labels)
        subj_labels = arrayfun(@(i) sprintf('S%02d', i), 1:n_subj, 'UniformOutput', false);
    end

    defaults = unbinnedMLEFitDefaults();
    lb = defaults.lb;
    ub = defaults.ub;
    if isfield(fit_opts, 'lb') && ~isempty(fit_opts.lb), lb = fit_opts.lb(:); end
    if isfield(fit_opts, 'ub') && ~isempty(fit_opts.ub), ub = fit_opts.ub(:); end
    span = max(abs(ub - lb), eps);

    rows = cell(n_subj * 18, 1);
    r = 0;
    for s = 1:n_subj
        sid = subj_list(s);
        for c = 1:18
            [m, prev_lvl, curr_lvl] = conditionSubscriptsFromIndex(c);
            if m == 1
                manip = 'contrast';
            else
                manip = 'precision';
            end

            mask = tbl_trials.subject_id == sid & ...
                tbl_trials.cond_manipulation == manip & ...
                tbl_trials.cond_prev == prev_lvl & tbl_trials.cond_curr == curr_lvl;
            n_trials = sum(mask);

            A = NaN; w = NaN; sg = NaN; beta = NaN; fwhm = NaN; ef = NaN;
            at_bound_any = false;
            at_bound_count = 0;
            if n_trials >= ip.Results.min_trials
                [pf, ef] = fitConditionMLE(tbl_trials.delta_theta(mask), ...
                    tbl_trials.x_probe(mask), tbl_trials.response(mask), fit_opts);
                if all(isfinite(pf))
                    A = pf(1);
                    w = pf(2);
                    sg = pf(3);
                    beta = pf(4);
                    fwhm = unbinnedWtoFwhm(w);
                    pf_col = pf(:);
                    bound_hits = (abs(pf_col - lb) ./ span < ip.Results.bound_tol) | ...
                        (abs(pf_col - ub) ./ span < ip.Results.bound_tol);
                    at_bound_any = any(bound_hits);
                    at_bound_count = sum(bound_hits);
                end
            end

            r = r + 1;
            rows{r} = {sid, string(subj_labels{s}), string(manip), prev_lvl, curr_lvl, ...
                n_trials, A, w, fwhm, sg, beta, ef, at_bound_any, at_bound_count};
        end
    end

    T = cell2table(vertcat(rows{:}), 'VariableNames', ...
        {'subject_id', 'subject_label', 'manipulation', 'cond_prev', 'cond_curr', ...
        'n_trials', 'A', 'w', 'fwhm', 'sigma', 'beta', 'exit_flag', ...
        'at_bound_any', 'at_bound_count'});
end
