function [params_row, curve_row, exit_row, at_bound_row] = bootstrapOneReplicate(tbl_trials, plan_row, grid, fit_opts, lb, ub, bound_tol, bootstrap_unit, fit_method)
% Single bootstrap iteration: 18 x 4 params, 18 x nGrid isolated DoG curves,
% 18 x 1 exit flags, 18 x 4 at-bound flags (for admission filtering downstream).

    if nargin < 5, lb = []; end
    if nargin < 6, ub = []; end
    if nargin < 7 || isempty(bound_tol), bound_tol = 1e-4; end
    if nargin < 8 || isempty(bootstrap_unit), bootstrap_unit = 'subject'; end
    if nargin < 9 || isempty(fit_method), fit_method = 'pooled'; end

    n_grid = numel(grid);
    num_conds = 18;
    params_row   = nan(num_conds, 4);
    curve_row    = nan(num_conds, n_grid);
    exit_row     = nan(num_conds, 1);
    at_bound_row = false(num_conds, 4);

    if ~isempty(lb) && ~isempty(ub)
        span = max(abs(ub(:)' - lb(:)'), eps);
    else
        span = [];
    end

    packed = generateBootstrapSample(tbl_trials, plan_row, bootstrap_unit);

    for c = 1:num_conds
        [m, prev, curr] = conditionSubscriptsFromIndex(c);
        mask = packed.manipulation == m & packed.cond_prev == prev & packed.cond_curr == curr;

        tbl_cond = packedConditionTable(packed, mask);
        [pfit, exitf] = fitConditionByMethod(tbl_cond, fit_opts, fit_method);
        params_row(c, :) = pfit;
        exit_row(c) = exitf;

        if all(isnan(pfit))
            curve_row(c, :) = nan(1, n_grid);
        else
            curve_row(c, :) = dogIsolated(grid, pfit(1), pfit(2))';
        end

        if ~isempty(span) && ~any(isnan(pfit))
            at_bound_row(c, :) = ...
                (abs(pfit - lb(:)') ./ span < bound_tol) | ...
                (abs(pfit - ub(:)') ./ span < bound_tol);
        end
    end

end

function tbl_cond = packedConditionTable(packed, mask)
    tbl_cond = table( ...
        packed.subject_id(mask), ...
        packed.delta_theta(mask), ...
        packed.x_probe(mask), ...
        packed.response(mask), ...
        'VariableNames', {'subject_id', 'delta_theta', 'x_probe', 'response'});
end
