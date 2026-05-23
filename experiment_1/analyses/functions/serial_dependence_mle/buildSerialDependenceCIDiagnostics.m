function diag_tbl = buildSerialDependenceCIDiagnostics(results)
% buildSerialDependenceCIDiagnostics  Per-cell diagnostics comparing percentile and BCa CIs.

    param_names = {'A', 'w', 'sigma', 'beta', 'FWHM'};
    num_conds = 18;
    n_rows = num_conds * numel(param_names);

    cond_manipulation = strings(n_rows, 1);
    cond_prev = zeros(n_rows, 1);
    cond_curr = zeros(n_rows, 1);
    parameter = strings(n_rows, 1);
    ci_method = repmat(string(results.ci_method), n_rows, 1);
    point = nan(n_rows, 1);
    boot_mean = nan(n_rows, 1);
    boot_median = nan(n_rows, 1);
    pc_lo = nan(n_rows, 1);
    pc_hi = nan(n_rows, 1);
    bca_lo = nan(n_rows, 1);
    bca_hi = nan(n_rows, 1);
    active_lo = nan(n_rows, 1);
    active_hi = nan(n_rows, 1);
    z0 = nan(n_rows, 1);
    acc = nan(n_rows, 1);
    admitted_fraction = nan(n_rows, 1);
    at_bound_fraction = nan(n_rows, 1);
    failed_fraction = nan(n_rows, 1);
    point_percentile_rank = nan(n_rows, 1);
    point_at_bound = false(n_rows, 1);
    extreme_z0 = false(n_rows, 1);
    extreme_acc = false(n_rows, 1);
    low_admitted_fraction = false(n_rows, 1);
    point_outside_percentile_ci = false(n_rows, 1);
    point_outside_bca_ci = false(n_rows, 1);
    bca_percentile_conflict = false(n_rows, 1);

    lb = results.lb(:)';
    ub = results.ub(:)';
    span = max(abs(ub - lb), eps);
    bound_tol = results.bound_tol;
    if isfield(results, 'failed_boot') && ~isempty(results.failed_boot)
        failed_boot = results.failed_boot;
    else
        failed_boot = all(isnan(results.params_boot), 3);
    end

    row = 0;
    for c = 1:num_conds
        [m, prev, curr] = conditionSubscriptsFromIndex(c);
        if m == 1
            manip = "contrast";
        else
            manip = "precision";
        end

        for ip = 1:numel(param_names)
            row = row + 1;
            pname = param_names{ip};
            cond_manipulation(row) = manip;
            cond_prev(row) = prev;
            cond_curr(row) = curr;
            parameter(row) = string(pname);
            failed_fraction(row) = mean(failed_boot(:, c));

            if ip <= 4
                k = ip;
                vals_all = squeeze(results.params_boot(:, c, k));
                if isfield(results, 'param_admitted') && ~isempty(results.param_admitted)
                    adm = squeeze(results.param_admitted(:, c, k));
                else
                    adm = results.admitted(:, c) & isfinite(vals_all);
                end
                vals = vals_all(adm);
                point(row) = results.overlay.params_point(c, k);
                pc_lo(row) = results.ci_percentile.param_lo(c, k);
                pc_hi(row) = results.ci_percentile.param_hi(c, k);
                bca_lo(row) = results.ci_bca.param_lo(c, k);
                bca_hi(row) = results.ci_bca.param_hi(c, k);
                active_lo(row) = results.ci_active.param_lo(c, k);
                active_hi(row) = results.ci_active.param_hi(c, k);
                if isfield(results.ci_bca, 'param_z0')
                    z0(row) = results.ci_bca.param_z0(c, k);
                    acc(row) = results.ci_bca.param_acc(c, k);
                end
                at_bound_fraction(row) = mean(squeeze(results.at_bound_boot(:, c, k)));
                point_at_bound(row) = isfinite(point(row)) && ...
                    ((abs(point(row) - lb(k)) / span(k) < bound_tol) || ...
                    (abs(point(row) - ub(k)) / span(k) < bound_tol));
            else
                vals_all_w = squeeze(results.params_boot(:, c, 2));
                if isfield(results, 'param_admitted') && ~isempty(results.param_admitted)
                    adm = squeeze(results.param_admitted(:, c, 2));
                else
                    adm = results.admitted(:, c) & isfinite(vals_all_w);
                end
                vals = serialDependenceWtoFwhm(vals_all_w(adm));
                point(row) = serialDependenceWtoFwhm(results.overlay.params_point(c, 2));
                pc_lo(row) = results.ci_percentile.fwhm_lo(c);
                pc_hi(row) = results.ci_percentile.fwhm_hi(c);
                bca_lo(row) = results.ci_bca.fwhm_lo(c);
                bca_hi(row) = results.ci_bca.fwhm_hi(c);
                active_lo(row) = results.ci_active.fwhm_lo(c);
                active_hi(row) = results.ci_active.fwhm_hi(c);
                if isfield(results.ci_bca, 'fwhm_z0')
                    z0(row) = results.ci_bca.fwhm_z0(c);
                    acc(row) = results.ci_bca.fwhm_acc(c);
                end
                at_bound_fraction(row) = mean(squeeze(results.at_bound_boot(:, c, 2)));
                w_point = results.overlay.params_point(c, 2);
                point_at_bound(row) = isfinite(w_point) && ...
                    ((abs(w_point - lb(2)) / span(2) < bound_tol) || ...
                    (abs(w_point - ub(2)) / span(2) < bound_tol));
            end

            vals = vals(isfinite(vals));
            admitted_fraction(row) = mean(adm);
            boot_mean(row) = mean(vals, 'omitnan');
            boot_median(row) = local_medianFinite(vals);
            if ~isempty(vals) && isfinite(point(row))
                point_percentile_rank(row) = mean(vals < point(row));
            end
            extreme_z0(row) = isfinite(z0(row)) && abs(z0(row)) > 2;
            extreme_acc(row) = isfinite(acc(row)) && abs(acc(row)) > 0.2;
            low_admitted_fraction(row) = isfinite(admitted_fraction(row)) && admitted_fraction(row) < 0.20;
            point_outside_percentile_ci(row) = local_outsideInterval(point(row), pc_lo(row), pc_hi(row));
            point_outside_bca_ci(row) = local_outsideInterval(point(row), bca_lo(row), bca_hi(row));
            bca_percentile_conflict(row) = local_ciConflict(pc_lo(row), pc_hi(row), bca_lo(row), bca_hi(row));
        end
    end

    diag_tbl = table(categorical(cond_manipulation), cond_prev, cond_curr, parameter, ci_method, ...
        point, boot_mean, boot_median, pc_lo, pc_hi, bca_lo, bca_hi, active_lo, active_hi, ...
        z0, acc, admitted_fraction, at_bound_fraction, failed_fraction, point_percentile_rank, ...
        point_at_bound, extreme_z0, extreme_acc, low_admitted_fraction, ...
        point_outside_percentile_ci, point_outside_bca_ci, bca_percentile_conflict, ...
        'VariableNames', {'cond_manipulation', 'cond_prev', 'cond_curr', 'parameter', 'ci_method', ...
        'point', 'boot_mean', 'boot_median', 'pc_lo', 'pc_hi', 'bca_lo', 'bca_hi', ...
        'active_lo', 'active_hi', 'z0', 'acc', 'admitted_fraction', 'at_bound_fraction', ...
        'failed_fraction', 'point_percentile_rank', 'point_at_bound', 'extreme_z0', ...
        'extreme_acc', 'low_admitted_fraction', 'point_outside_percentile_ci', ...
        'point_outside_bca_ci', 'bca_percentile_conflict'});
end

function m = local_medianFinite(values)
    values = values(isfinite(values));
    if isempty(values)
        m = NaN;
    else
        m = median(values);
    end
end

function tf = local_outsideInterval(x, lo, hi)
    if ~isfinite(x) || ~isfinite(lo) || ~isfinite(hi)
        tf = false;
        return
    end
    tf = x < min(lo, hi) || x > max(lo, hi);
end

function tf = local_ciConflict(pc_lo, pc_hi, bca_lo, bca_hi)
    if ~all(isfinite([pc_lo, pc_hi, bca_lo, bca_hi]))
        tf = false;
        return
    end
    pc_excludes = pc_lo > 0 || pc_hi < 0;
    bca_excludes = bca_lo > 0 || bca_hi < 0;
    if pc_excludes ~= bca_excludes
        tf = true;
    elseif pc_excludes && bca_excludes
        tf = sign(mean([pc_lo, pc_hi])) ~= sign(mean([bca_lo, bca_hi]));
    else
        tf = false;
    end
end
