function ht = computeTargetedDoGHypothesisTests(results, n_back)
% computeTargetedDoGHypothesisTests  Targeted BCa tests for A and DoG range.
%
% These tests ask whether the high-to-low information span within each
% manipulation changes DoG amplitude (A) or range (FWHM), and whether that
% span is larger for the hypothesized manipulation.

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
            ht = table();
            return
        end
    end
    if ~isfield(results.overlay, 'params_point')
        ht = table();
        return
    end
    if use_bca && (~isfield(results, 'jackknife') || ~isfield(results.jackknife, 'params'))
        ht = table();
        return
    end

    specs = local_buildSpecs();
    n_specs = numel(specs);

    n_back_col = repmat(n_back, n_specs, 1);
    name = strings(n_specs, 1);
    hypothesis = strings(n_specs, 1);
    parameter = strings(n_specs, 1);
    parameter_label = strings(n_specs, 1);
    axis = strings(n_specs, 1);
    test_family = strings(n_specs, 1);
    comparison = strings(n_specs, 1);
    expected_positive = false(n_specs, 1);
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
    estimate_boot_mismatch = false(n_specs, 1);
    valid_for_inference = false(n_specs, 1);
    ci_method_col = repmat(string(ci_method), n_specs, 1);

    alpha = 1 - (results.ci_prctile(2) - results.ci_prctile(1)) / 100;
    min_admit_fraction = 0.20;
    max_abs_z0 = 2;

    for i_spec = 1:n_specs
        spec = specs(i_spec);
        vals = local_values(results, spec.parameter);
        if use_bca
            jk_vals = local_jackknifeValues(results, spec.parameter);
        else
            jk_vals = nan(results.n_subj, 18);
        end
        point_vals = local_pointValues(results, spec.parameter);
        [est, boot_v, jk_v, n_good] = local_evalContrast(results, vals, jk_vals, point_vals, spec.weights, spec.parameter);

        name(i_spec) = string(spec.name);
        hypothesis(i_spec) = string(spec.hypothesis);
        parameter(i_spec) = string(spec.parameter);
        parameter_label(i_spec) = string(spec.parameter_label);
        axis(i_spec) = string(spec.axis);
        test_family(i_spec) = string(spec.test_family);
        comparison(i_spec) = string(spec.comparison);
        expected_positive(i_spec) = spec.expected_positive;
        estimate(i_spec) = est;
        n_admit(i_spec) = n_good;
        admit_fraction(i_spec) = n_good ./ size(vals, 1);
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
        if use_bca
            [a1, a2, z0, acc] = local_bcaQuantiles(boot_v, jk_v, est, alpha);
            [bca_lo(i_spec), bca_hi(i_spec)] = local_orderedPrctileCI(bv, 100 * a1, 100 * a2);
            z0_col(i_spec) = z0;
            acc_col(i_spec) = acc;
            p_bca(i_spec) = local_bcaTwoSidedPValue(boot_v, z0, acc);
        else
            bca_lo(i_spec) = pc_lo(i_spec);
            bca_hi(i_spec) = pc_hi(i_spec);
            p_bca(i_spec) = local_percentileTwoSidedPValue(boot_v);
        end
        estimate_in_bca_ci(i_spec) = local_inInterval(est, bca_lo(i_spec), bca_hi(i_spec));
        bca_excludes_zero(i_spec) = local_excludesZero(bca_lo(i_spec), bca_hi(i_spec));
        bca_pc_conflict(i_spec) = local_ciConflict(pc_lo(i_spec), pc_hi(i_spec), bca_lo(i_spec), bca_hi(i_spec));
        extreme_z0(i_spec) = isfinite(z0) && abs(z0) > max_abs_z0;
        estimate_boot_mismatch(i_spec) = ~estimate_in_pc_ci(i_spec);
        valid_for_inference(i_spec) = estimate_in_pc_ci(i_spec) && ...
            ~bca_pc_conflict(i_spec) && ~extreme_z0(i_spec) && ...
            ~low_admit_fraction(i_spec);
    end

    p_holm = local_holmAdjust(p_bca);
    p_fdr_bh = local_bhAdjust(p_bca);
    p_bca_label = arrayfun(@local_formatPLabel, p_bca);
    p_holm_label = arrayfun(@local_formatPLabel, p_holm);
    p_fdr_bh_label = arrayfun(@local_formatPLabel, p_fdr_bh);
    directed_support = expected_positive & estimate > 0 & ...
        bca_excludes_zero & pc_excludes_zero & bca_lo > 0 & pc_lo > 0;
    nondirectional_support = ~expected_positive & bca_excludes_zero & pc_excludes_zero;
    supports_hypothesis = valid_for_inference & (directed_support | nondirectional_support);

    ht = table(n_back_col, name, hypothesis, parameter, parameter_label, axis, ...
        test_family, comparison, expected_positive, estimate, boot_median, boot_mean, ...
        bca_lo, bca_hi, pc_lo, pc_hi, p_bca, p_bca_label, p_holm, p_holm_label, ...
        p_fdr_bh, p_fdr_bh_label, supports_hypothesis, n_admit, admit_fraction, ...
        z0_col, acc_col, estimate_in_pc_ci, estimate_in_bca_ci, pc_excludes_zero, ...
        bca_excludes_zero, bca_pc_conflict, extreme_z0, low_admit_fraction, ...
        estimate_boot_mismatch, valid_for_inference, ci_method_col, ...
        'VariableNames', {'n_back', 'name', 'hypothesis', 'parameter', ...
        'parameter_label', 'axis', 'test_family', 'comparison', ...
        'expected_positive', 'estimate', 'boot_median', 'boot_mean', ...
        'bca_lo', 'bca_hi', 'pc_lo', 'pc_hi', ...
        'p_bca', 'p_bca_label', 'p_holm', 'p_holm_label', ...
        'p_fdr_bh', 'p_fdr_bh_label', 'supports_hypothesis', ...
        'n_admit', 'admit_fraction', 'z0', 'acc', ...
        'estimate_in_pc_ci', 'estimate_in_bca_ci', 'pc_excludes_zero', ...
        'bca_excludes_zero', 'bca_pc_conflict', 'extreme_z0', ...
        'low_admit_fraction', 'estimate_boot_mismatch', 'valid_for_inference', ...
        'ci_method'});
end

function specs = local_buildSpecs()
    specs = struct('name', {}, 'hypothesis', {}, 'parameter', {}, ...
        'parameter_label', {}, 'axis', {}, 'test_family', {}, ...
        'comparison', {}, 'expected_positive', {}, 'weights', {});

    axes_to_use = {'average_prev_curr', 'previous_level', 'current_level'};
    for i_axis = 1:numel(axes_to_use)
        ax = axes_to_use{i_axis};
        wc = local_axisWeights(1, ax); % contrast
        wp = local_axisWeights(2, ax); % precision

        specs(end+1) = local_spec('contrast_A_high_to_low', ...
            'Contrast effect on bias strength', 'A', 'Amplitude A (deg)', ax, ...
            'within_manipulation', 'Does lower contrast change DoG amplitude?', true, wc); %#ok<AGROW>
        specs(end+1) = local_spec('precision_A_high_to_low', ...
            'Precision effect on bias strength', 'A', 'Amplitude A (deg)', ax, ...
            'within_manipulation', 'Does lower precision change DoG amplitude?', true, wp); %#ok<AGROW>
        specs(end+1) = local_spec('contrast_minus_precision_A', ...
            'Contrast effect > precision effect on bias strength', 'A', 'Amplitude A (deg)', ax, ...
            'between_manipulation', 'Is the contrast effect on amplitude larger than the precision effect?', true, wc - wp); %#ok<AGROW>

        specs(end+1) = local_spec('precision_FWHM_high_to_low', ...
            'Precision effect on bias range', 'FWHM', 'Range FWHM (deg)', ax, ...
            'within_manipulation', 'Does lower precision change DoG range?', true, wp); %#ok<AGROW>
        specs(end+1) = local_spec('contrast_FWHM_high_to_low', ...
            'Contrast effect on bias range', 'FWHM', 'Range FWHM (deg)', ax, ...
            'within_manipulation', 'Does lower contrast change DoG range?', true, wc); %#ok<AGROW>
        specs(end+1) = local_spec('precision_minus_contrast_FWHM', ...
            'Precision effect > contrast effect on bias range', 'FWHM', 'Range FWHM (deg)', ax, ...
            'between_manipulation', 'Is the precision effect on range larger than the contrast effect?', true, wp - wc); %#ok<AGROW>
    end

    for manipulation = 1:2
        if manipulation == 1
            prefix = 'contrast';
            manip_label = 'Contrast';
        else
            prefix = 'precision';
            manip_label = 'Precision';
        end
        wb = local_boundaryInteractionWeights(manipulation);
        specs(end+1) = local_spec(sprintf('%s_A_boundary_interaction', prefix), ...
            sprintf('%s boundary interaction on bias strength', manip_label), ...
            'A', 'Amplitude A (deg)', 'boundary_interaction', ...
            'within_manipulation', 'Does the L3-L1 current-level effect depend on previous-level endpoint?', false, wb); %#ok<AGROW>
        specs(end+1) = local_spec(sprintf('%s_FWHM_boundary_interaction', prefix), ...
            sprintf('%s boundary interaction on bias range', manip_label), ...
            'FWHM', 'Range FWHM (deg)', 'boundary_interaction', ...
            'within_manipulation', 'Does the L3-L1 current-level effect depend on previous-level endpoint?', false, wb); %#ok<AGROW>
    end
end

function spec = local_spec(name, hypothesis, parameter, parameter_label, axis, test_family, comparison, expected_positive, weights)
    spec = struct('name', sprintf('%s.%s', name, axis), ...
        'hypothesis', hypothesis, ...
        'parameter', parameter, ...
        'parameter_label', parameter_label, ...
        'axis', axis, ...
        'test_family', test_family, ...
        'comparison', comparison, ...
        'expected_positive', expected_positive, ...
        'weights', weights(:)');
end

function w = local_axisWeights(manipulation, axis_name)
    w_prev = local_prevWeights(manipulation, 3, 1);
    w_curr = local_currWeights(manipulation, 3, 1);
    switch axis_name
        case 'previous_level'
            w = w_prev;
        case 'current_level'
            w = w_curr;
        case 'average_prev_curr'
            w = 0.5 .* (w_prev + w_curr);
        otherwise
            error('computeTargetedDoGHypothesisTests:badAxis', 'Unknown axis: %s', axis_name);
    end
end

function w = local_prevWeights(manipulation, low_level, high_level)
    w = zeros(1, 18);
    for curr = 1:3
        w(local_cidx(manipulation, low_level, curr)) = 1/3;
        w(local_cidx(manipulation, high_level, curr)) = -1/3;
    end
end

function w = local_currWeights(manipulation, low_level, high_level)
    w = zeros(1, 18);
    for prev = 1:3
        w(local_cidx(manipulation, prev, low_level)) = 1/3;
        w(local_cidx(manipulation, prev, high_level)) = -1/3;
    end
end

function w = local_boundaryInteractionWeights(manipulation)
    w = zeros(1, 18);
    w(local_cidx(manipulation, 3, 3)) = 1;
    w(local_cidx(manipulation, 3, 1)) = -1;
    w(local_cidx(manipulation, 1, 3)) = -1;
    w(local_cidx(manipulation, 1, 1)) = 1;
end

function c = local_cidx(manipulation, prev, curr)
    c = (manipulation - 1) * 9 + (prev - 1) * 3 + curr;
end

function vals = local_values(results, parameter)
    switch parameter
        case 'A'
            vals = reshape(results.params_boot(:, :, 1), [results.B, 18]);
        case 'FWHM'
            w = reshape(results.params_boot(:, :, 2), [], 1);
            f = serialDependenceWtoFwhm(w);
            vals = reshape(f, [results.B, 18]);
        otherwise
            error('computeTargetedDoGHypothesisTests:badParameter', 'Unknown parameter: %s', parameter);
    end
end

function vals = local_jackknifeValues(results, parameter)
    n_subj = results.n_subj;
    switch parameter
        case 'A'
            vals = reshape(results.jackknife.params(:, :, 1), [n_subj, 18]);
        case 'FWHM'
            w = reshape(results.jackknife.params(:, :, 2), [], 1);
            f = serialDependenceWtoFwhm(w);
            vals = reshape(f, [n_subj, 18]);
    end
end

function vals = local_pointValues(results, parameter)
    switch parameter
        case 'A'
            vals = results.overlay.params_point(:, 1);
        case 'FWHM'
            vals = serialDependenceWtoFwhm(results.overlay.params_point(:, 2));
    end
    vals = vals(:);
end

function [estimate, boot_v, jk_v, n_good] = local_evalContrast(results, vals, jk_vals, point_vals, weights, parameter)
    w = weights(:);
    involved = abs(w) > sqrt(eps);
    adm = local_paramAdmission(results, parameter);
    good = all(adm(:, involved), 2) & all(isfinite(vals(:, involved)), 2);
    boot_v = nan(size(vals, 1), 1);
    boot_v(good) = vals(good, :) * w;
    n_good = sum(good);

    good_jk = all(isfinite(jk_vals(:, involved)), 2);
    if isfield(results, 'jackknife') && isfield(results.jackknife, 'at_bound') && ...
            ~isempty(results.jackknife.at_bound) && isfield(results, 'discard_at_bound') && results.discard_at_bound
        jk_bound = local_jackknifeBound(results, parameter);
        good_jk = good_jk & ~any(jk_bound(:, involved), 2);
    end
    jk_v = nan(size(jk_vals, 1), 1);
    jk_v(good_jk) = jk_vals(good_jk, :) * w;

    if all(isfinite(point_vals(involved)))
        estimate = point_vals' * w;
    else
        estimate = NaN;
    end
end

function adm = local_paramAdmission(results, parameter)
    switch parameter
        case 'A'
            k = 1;
        case 'FWHM'
            k = 2;
        otherwise
            k = 1;
    end
    if isfield(results, 'param_admitted') && ~isempty(results.param_admitted)
        adm = reshape(results.param_admitted(:, :, k), [results.B, 18]);
    else
        adm = results.admitted;
    end
end

function jk_bound = local_jackknifeBound(results, parameter)
    switch parameter
        case 'A'
            k = 1;
        case 'FWHM'
            k = 2;
        otherwise
            k = 1;
    end
    jk_bound = reshape(results.jackknife.at_bound(:, :, k), [results.n_subj, 18]);
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
