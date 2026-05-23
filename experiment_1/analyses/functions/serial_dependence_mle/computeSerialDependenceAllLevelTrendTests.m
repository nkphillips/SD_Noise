function tt = computeSerialDependenceAllLevelTrendTests(results, n_back)
% computeSerialDependenceAllLevelTrendTests  Bootstrap trend tests over all 3 x 3 cells.
%
% Fits metric ~ previous_level + current_level + previous_level:current_level
% to each manipulation's 3 x 3 cell estimates, with levels coded [-1, 0, +1].
% A and FWHM use parameter-specific bootstrap/jackknife admission.

    if nargin < 2 || isempty(n_back)
        n_back = NaN;
    end

    ci_method = 'bca';
    if isfield(results, 'ci_method') && ~isempty(results.ci_method)
        ci_method = lower(char(results.ci_method));
    end
    bootstrap_unit = 'subject';
    if isfield(results, 'bootstrap_unit') && ~isempty(results.bootstrap_unit)
        bootstrap_unit = lower(char(results.bootstrap_unit));
    end
    use_bca = strcmp(ci_method, 'bca') && strcmp(bootstrap_unit, 'subject');

    required_fields = {'params_boot', 'admitted', 'overlay', 'ci_prctile'};
    for i = 1:numel(required_fields)
        if ~isfield(results, required_fields{i})
            tt = table();
            return
        end
    end
    if ~isfield(results.overlay, 'params_point')
        tt = table();
        return
    end
    if use_bca && (~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params'))
        tt = table();
        return
    end

    specs = local_buildSpecs();
    n_specs = numel(specs);
    alpha = 1 - (results.ci_prctile(2) - results.ci_prctile(1)) / 100;
    min_admit_fraction = 0.20;
    max_abs_z0 = 2;

    n_back_col = repmat(n_back, n_specs, 1);
    name = strings(n_specs, 1);
    parameter = strings(n_specs, 1);
    parameter_label = strings(n_specs, 1);
    manipulation = strings(n_specs, 1);
    term = strings(n_specs, 1);
    comparison = strings(n_specs, 1);
    estimate = nan(n_specs, 1);
    boot_median = nan(n_specs, 1);
    boot_mean = nan(n_specs, 1);
    bca_lo = nan(n_specs, 1);
    bca_hi = nan(n_specs, 1);
    pc_lo = nan(n_specs, 1);
    pc_hi = nan(n_specs, 1);
    p_bca = nan(n_specs, 1);
    n_admit = zeros(n_specs, 1);
    admit_fraction = nan(n_specs, 1);
    z0_col = nan(n_specs, 1);
    acc_col = nan(n_specs, 1);
    estimate_in_pc_ci = false(n_specs, 1);
    estimate_in_bca_ci = false(n_specs, 1);
    pc_excludes_zero = false(n_specs, 1);
    bca_excludes_zero = false(n_specs, 1);
    bca_pc_conflict = false(n_specs, 1);
    extreme_z0 = false(n_specs, 1);
    low_admit_fraction = true(n_specs, 1);
    valid_for_inference = false(n_specs, 1);
    ci_method_col = repmat(string(ci_method), n_specs, 1);

    coeffs = local_coeffDistributions(results);

    for i_spec = 1:n_specs
        spec = specs(i_spec);
        vals = coeffs.(spec.parameter);
        [est, boot_v, jk_v, n_good] = local_evalSpec(vals, spec, use_bca);

        name(i_spec) = string(spec.name);
        parameter(i_spec) = string(spec.parameter);
        parameter_label(i_spec) = string(spec.parameter_label);
        manipulation(i_spec) = string(spec.manipulation);
        term(i_spec) = string(spec.term);
        comparison(i_spec) = string(spec.comparison);
        estimate(i_spec) = est;
        n_admit(i_spec) = n_good;
        admit_fraction(i_spec) = n_good ./ results.B;
        low_admit_fraction(i_spec) = admit_fraction(i_spec) < min_admit_fraction;

        bv = boot_v(isfinite(boot_v));
        if numel(bv) < 10
            continue
        end

        boot_median(i_spec) = median(bv);
        boot_mean(i_spec) = mean(bv);
        pc_lo(i_spec) = prctile(bv, results.ci_prctile(1));
        pc_hi(i_spec) = prctile(bv, results.ci_prctile(2));
        estimate_in_pc_ci(i_spec) = local_inInterval(est, pc_lo(i_spec), pc_hi(i_spec));
        pc_excludes_zero(i_spec) = local_excludesZero(pc_lo(i_spec), pc_hi(i_spec));

        z0 = NaN;
        acc = NaN;
        if use_bca
            [a1, a2, z0, acc] = local_bcaQuantiles(boot_v, jk_v, est, alpha);
            [bca_lo(i_spec), bca_hi(i_spec)] = local_orderedPrctileCI(bv, 100 * a1, 100 * a2);
            p_bca(i_spec) = local_bcaTwoSidedPValue(boot_v, z0, acc);
        else
            bca_lo(i_spec) = pc_lo(i_spec);
            bca_hi(i_spec) = pc_hi(i_spec);
            p_bca(i_spec) = local_percentileTwoSidedPValue(boot_v);
        end
        z0_col(i_spec) = z0;
        acc_col(i_spec) = acc;
        estimate_in_bca_ci(i_spec) = local_inInterval(est, bca_lo(i_spec), bca_hi(i_spec));
        bca_excludes_zero(i_spec) = local_excludesZero(bca_lo(i_spec), bca_hi(i_spec));
        bca_pc_conflict(i_spec) = local_ciConflict(pc_lo(i_spec), pc_hi(i_spec), bca_lo(i_spec), bca_hi(i_spec));
        extreme_z0(i_spec) = isfinite(z0) && abs(z0) > max_abs_z0;
        valid_for_inference(i_spec) = estimate_in_pc_ci(i_spec) && ...
            ~bca_pc_conflict(i_spec) && ~extreme_z0(i_spec) && ~low_admit_fraction(i_spec);
    end

    p_holm = local_holmAdjust(p_bca);
    p_fdr_bh = local_bhAdjust(p_bca);
    p_bca_label = arrayfun(@local_formatPLabel, p_bca);
    p_holm_label = arrayfun(@local_formatPLabel, p_holm);
    p_fdr_bh_label = arrayfun(@local_formatPLabel, p_fdr_bh);
    supports_effect = valid_for_inference & bca_excludes_zero & pc_excludes_zero;

    tt = table(n_back_col, name, parameter, parameter_label, manipulation, term, ...
        comparison, estimate, boot_median, boot_mean, bca_lo, bca_hi, pc_lo, pc_hi, ...
        p_bca, p_bca_label, p_holm, p_holm_label, p_fdr_bh, p_fdr_bh_label, ...
        supports_effect, n_admit, admit_fraction, z0_col, acc_col, ...
        estimate_in_pc_ci, estimate_in_bca_ci, pc_excludes_zero, bca_excludes_zero, ...
        bca_pc_conflict, extreme_z0, low_admit_fraction, valid_for_inference, ci_method_col, ...
        'VariableNames', {'n_back', 'name', 'parameter', 'parameter_label', ...
        'manipulation', 'term', 'comparison', 'estimate', 'boot_median', 'boot_mean', ...
        'bca_lo', 'bca_hi', 'pc_lo', 'pc_hi', 'p_bca', 'p_bca_label', ...
        'p_holm', 'p_holm_label', 'p_fdr_bh', 'p_fdr_bh_label', ...
        'supports_effect', 'n_admit', 'admit_fraction', 'z0', 'acc', ...
        'estimate_in_pc_ci', 'estimate_in_bca_ci', 'pc_excludes_zero', ...
        'bca_excludes_zero', 'bca_pc_conflict', 'extreme_z0', ...
        'low_admit_fraction', 'valid_for_inference', 'ci_method'});
end

function specs = local_buildSpecs()
    specs = struct('name', {}, 'parameter', {}, 'parameter_label', {}, ...
        'manipulation', {}, 'term', {}, 'comparison', {});
    params = {'A', 'FWHM'};
    labels = {'Amplitude A (deg)', 'Range FWHM (deg)'};
    terms = {'previous_slope', 'current_slope', 'previous_current_interaction'};
    manips = {'contrast', 'precision'};

    for i_param = 1:numel(params)
        for i_manip = 1:numel(manips)
            for i_term = 1:numel(terms)
                specs(end+1) = local_spec(params{i_param}, labels{i_param}, ... %#ok<AGROW>
                    manips{i_manip}, terms{i_term}, 'within_manipulation');
            end
        end
        for i_term = 1:numel(terms)
            specs(end+1) = local_spec(params{i_param}, labels{i_param}, ... %#ok<AGROW>
                'precision_minus_contrast', terms{i_term}, 'between_manipulation');
        end
    end
end

function spec = local_spec(parameter, parameter_label, manipulation, term, comparison)
    spec = struct('name', sprintf('%s.%s.%s', parameter, manipulation, term), ...
        'parameter', parameter, ...
        'parameter_label', parameter_label, ...
        'manipulation', manipulation, ...
        'term', term, ...
        'comparison', comparison);
end

function coeffs = local_coeffDistributions(results)
    coeffs = struct();
    coeffs.A = local_oneParameterDistributions(results, 'A');
    coeffs.FWHM = local_oneParameterDistributions(results, 'FWHM');
end

function vals = local_oneParameterDistributions(results, parameter)
    B = results.B;
    n_subj = results.n_subj;
    vals = struct();
    vals.point = nan(2, 3);
    vals.boot = nan(B, 2, 3);
    vals.boot_good = false(B, 2);
    vals.jk = nan(n_subj, 2, 3);
    vals.jk_good = false(n_subj, 2);

    metric_boot = local_metricValues(results, parameter);
    metric_point = local_pointValues(results, parameter);
    metric_jk = local_jackknifeValues(results, parameter);
    adm = local_paramAdmission(results, parameter);
    jk_bound = local_jackknifeBound(results, parameter);
    X = local_designMatrix();

    for m = 1:2
        idx = local_manipIndices(m);
        vals.point(m, :) = local_fitTrend(metric_point(idx), X);

        boot_good = all(adm(:, idx), 2) & all(isfinite(metric_boot(:, idx)), 2);
        vals.boot_good(:, m) = boot_good;
        for b = find(boot_good)'
            vals.boot(b, m, :) = local_fitTrend(metric_boot(b, idx)', X);
        end

        jk_good = all(isfinite(metric_jk(:, idx)), 2) & ~any(jk_bound(:, idx), 2);
        vals.jk_good(:, m) = jk_good;
        for s = find(jk_good)'
            vals.jk(s, m, :) = local_fitTrend(metric_jk(s, idx)', X);
        end
    end
end

function metric = local_metricValues(results, parameter)
    switch parameter
        case 'A'
            metric = reshape(results.params_boot(:, :, 1), [results.B, 18]);
        case 'FWHM'
            w = reshape(results.params_boot(:, :, 2), [], 1);
            metric = reshape(serialDependenceWtoFwhm(w), [results.B, 18]);
        otherwise
            error('computeSerialDependenceAllLevelTrendTests:badParameter', 'Unknown parameter: %s', parameter);
    end
end

function metric = local_pointValues(results, parameter)
    switch parameter
        case 'A'
            metric = results.overlay.params_point(:, 1);
        case 'FWHM'
            metric = serialDependenceWtoFwhm(results.overlay.params_point(:, 2));
    end
    metric = metric(:);
end

function metric = local_jackknifeValues(results, parameter)
    if ~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params')
        metric = nan(results.n_subj, 18);
        return
    end
    switch parameter
        case 'A'
            metric = reshape(results.jackknife.params(:, :, 1), [results.n_subj, 18]);
        case 'FWHM'
            w = reshape(results.jackknife.params(:, :, 2), [], 1);
            metric = reshape(serialDependenceWtoFwhm(w), [results.n_subj, 18]);
    end
end

function adm = local_paramAdmission(results, parameter)
    switch parameter
        case 'A'
            k = 1;
        case 'FWHM'
            k = 2;
    end
    if isfield(results, 'param_admitted') && ~isempty(results.param_admitted)
        adm = reshape(results.param_admitted(:, :, k), [results.B, 18]);
    else
        adm = results.admitted;
    end
end

function jk_bound = local_jackknifeBound(results, parameter)
    if ~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'at_bound') || ...
            isempty(results.jackknife.at_bound) || ~isfield(results, 'discard_at_bound') || ...
            ~results.discard_at_bound
        jk_bound = false(results.n_subj, 18);
        return
    end
    switch parameter
        case 'A'
            k = 1;
        case 'FWHM'
            k = 2;
    end
    jk_bound = reshape(results.jackknife.at_bound(:, :, k), [results.n_subj, 18]);
end

function X = local_designMatrix()
    levels = [-1; 0; 1];
    prev = zeros(9, 1);
    curr = zeros(9, 1);
    r = 0;
    for p = 1:3
        for c = 1:3
            r = r + 1;
            prev(r) = levels(p);
            curr(r) = levels(c);
        end
    end
    X = [ones(9, 1), prev, curr, prev .* curr];
end

function idx = local_manipIndices(manipulation)
    idx = (manipulation - 1) * 9 + (1:9);
end

function beta = local_fitTrend(y, X)
    y = y(:);
    ok = isfinite(y);
    beta_full = nan(4, 1);
    if sum(ok) >= size(X, 2)
        beta_full = X(ok, :) \ y(ok);
    end
    beta = beta_full(2:4)';
end

function [estimate, boot_v, jk_v, n_good] = local_evalSpec(vals, spec, use_bca)
    term_idx = local_termIndex(spec.term);
    if strcmp(spec.manipulation, 'contrast')
        estimate = vals.point(1, term_idx);
        boot_v = squeeze(vals.boot(:, 1, term_idx));
        jk_v = squeeze(vals.jk(:, 1, term_idx));
        good = vals.boot_good(:, 1) & isfinite(boot_v);
    elseif strcmp(spec.manipulation, 'precision')
        estimate = vals.point(2, term_idx);
        boot_v = squeeze(vals.boot(:, 2, term_idx));
        jk_v = squeeze(vals.jk(:, 2, term_idx));
        good = vals.boot_good(:, 2) & isfinite(boot_v);
    else
        estimate = vals.point(2, term_idx) - vals.point(1, term_idx);
        boot_v = squeeze(vals.boot(:, 2, term_idx)) - squeeze(vals.boot(:, 1, term_idx));
        jk_v = squeeze(vals.jk(:, 2, term_idx)) - squeeze(vals.jk(:, 1, term_idx));
        good = vals.boot_good(:, 1) & vals.boot_good(:, 2) & isfinite(boot_v);
    end
    boot_v(~good) = NaN;
    n_good = sum(good);
    if ~use_bca
        jk_v(:) = NaN;
    end
end

function idx = local_termIndex(term)
    switch char(term)
        case 'previous_slope'
            idx = 1;
        case 'current_slope'
            idx = 2;
        case 'previous_current_interaction'
            idx = 3;
        otherwise
            error('computeSerialDependenceAllLevelTrendTests:badTerm', 'Unknown term: %s', char(term));
    end
end

function tf = local_inInterval(x, lo, hi)
    if ~isfinite(x) || ~isfinite(lo) || ~isfinite(hi)
        tf = false;
        return
    end
    tf = x >= min(lo, hi) && x <= max(lo, hi);
end

function tf = local_excludesZero(lo, hi)
    if ~isfinite(lo) || ~isfinite(hi)
        tf = false;
        return
    end
    tf = min(lo, hi) > 0 || max(lo, hi) < 0;
end

function tf = local_ciConflict(pc_lo, pc_hi, bca_lo, bca_hi)
    pc_sig = local_excludesZero(pc_lo, pc_hi);
    bca_sig = local_excludesZero(bca_lo, bca_hi);
    if pc_sig ~= bca_sig
        tf = true;
        return
    end
    if pc_sig && bca_sig
        tf = sign(mean([pc_lo, pc_hi])) ~= sign(mean([bca_lo, bca_hi]));
    else
        tf = false;
    end
end

function [a1, a2, z0, acc] = local_bcaQuantiles(boot_vals, jk_vals, point_est, alpha)
    boot_vals = boot_vals(isfinite(boot_vals));
    jk_vals = jk_vals(isfinite(jk_vals));
    a1 = alpha / 2; a2 = 1 - alpha / 2; z0 = NaN; acc = NaN;
    if numel(boot_vals) < 2 || numel(jk_vals) < 2 || ~isfinite(point_est)
        return
    end

    p_hat = mean(boot_vals < point_est);
    p_hat = min(max(p_hat, 1 / (numel(boot_vals) + 1)), 1 - 1 / (numel(boot_vals) + 1));
    z0 = norminv(p_hat);
    jk_mean = mean(jk_vals);
    num_acc = sum((jk_mean - jk_vals).^3);
    den_acc = 6 * (sum((jk_mean - jk_vals).^2)).^(3/2);
    if den_acc == 0
        acc = 0;
    else
        acc = num_acc / den_acc;
    end

    z_lo = norminv(alpha / 2);
    z_hi = norminv(1 - alpha / 2);
    a1 = normcdf(z0 + (z0 + z_lo) / (1 - acc * (z0 + z_lo)));
    a2 = normcdf(z0 + (z0 + z_hi) / (1 - acc * (z0 + z_hi)));
    if ~isfinite(a1) || a1 <= 0, a1 = alpha / 2; end
    if ~isfinite(a2) || a2 >= 1, a2 = 1 - alpha / 2; end
end

function [lo, hi] = local_orderedPrctileCI(values, prc_lo, prc_hi)
    q1 = prctile(values, prc_lo);
    q2 = prctile(values, prc_hi);
    lo = min(q1, q2);
    hi = max(q1, q2);
end

function p = local_bcaTwoSidedPValue(boot_v, z0, acc)
    bv = boot_v(isfinite(boot_v));
    n = numel(bv);
    if n < 10 || ~isfinite(z0) || ~isfinite(acc)
        p = NaN;
        return
    end

    q = mean(bv < 0);
    q = min(max(q, 1 / (n + 1)), 1 - 1 / (n + 1));
    p_perc = 2 * min(q, 1 - q);
    z_q = norminv(q);
    s = z_q - z0;
    denom = 1 + acc * s;
    if denom <= 0 || ~isfinite(denom)
        p = p_perc;
        return
    end

    f_bca_0 = normcdf(s / denom - z0);
    if ~isfinite(f_bca_0)
        p = p_perc;
    else
        p = min(max(2 * min(f_bca_0, 1 - f_bca_0), 0), 1);
    end
end

function p = local_percentileTwoSidedPValue(boot_v)
    bv = boot_v(isfinite(boot_v));
    n = numel(bv);
    if n < 10
        p = NaN;
        return
    end

    q = mean(bv < 0);
    q = min(max(q, 1 / (n + 1)), 1 - 1 / (n + 1));
    p = min(max(2 * min(q, 1 - q), 0), 1);
end

function p_adj = local_holmAdjust(p_raw)
    p_adj = nan(size(p_raw));
    valid = isfinite(p_raw);
    pv = p_raw(valid);
    if isempty(pv)
        return
    end
    [pv_sorted, sort_idx] = sort(pv, 'ascend');
    k = numel(pv_sorted);
    p_step = pv_sorted .* (k - (1:k)' + 1);
    p_step = min(cummax(p_step), 1);
    p_resorted = nan(k, 1);
    p_resorted(sort_idx) = p_step;
    p_adj(valid) = p_resorted;
end

function p_adj = local_bhAdjust(p_raw)
    p_adj = nan(size(p_raw));
    valid = isfinite(p_raw);
    pv = p_raw(valid);
    if isempty(pv)
        return
    end
    [pv_sorted, sort_idx] = sort(pv, 'ascend');
    k = numel(pv_sorted);
    q_step = pv_sorted .* (k ./ (1:k)');
    q_step = min(flipud(cummin(flipud(q_step))), 1);
    q_resorted = nan(k, 1);
    q_resorted(sort_idx) = q_step;
    p_adj(valid) = q_resorted;
end

function s = local_formatPLabel(p)
    if ~isfinite(p)
        s = "NA";
    elseif p < 0.001
        s = "<0.001***";
    elseif p < 0.01
        s = string(sprintf('%.3f**', p));
    elseif p < 0.05
        s = string(sprintf('%.3f*', p));
    else
        s = string(sprintf('%.3f', p));
    end
end
