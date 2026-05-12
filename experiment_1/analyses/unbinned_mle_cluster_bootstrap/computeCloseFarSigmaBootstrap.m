function close_far = computeCloseFarSigmaBootstrap(tbl_trials, boot_plans, params_boot, admitted, jk_params, subj_list, opts)
% computeCloseFarSigmaBootstrap  Close-vs-far psychometric sigma diagnostic.
%
% Uses the fitted DoG+bias term as fixed mu_t within each condition, then
% estimates one sigma for close trials and one sigma for far trials.

    if nargin < 7 || isempty(opts)
        opts = struct();
    end
    opts = local_defaults(opts);

    if isempty(boot_plans)
        B = numel(opts.boot_plan_seeds);
    else
        B = size(boot_plans, 1);
    end
    num_conds = 18;
    n_subj = numel(subj_list);
    n_trials = height(tbl_trials);

    point_params = opts.point_params;
    sigma_point = nan(num_conds, 3);   % close, far, far-close
    n_close = nan(num_conds, 1);
    n_far = nan(num_conds, 1);
    for c = 1:num_conds
        [sigma_point(c, 1), sigma_point(c, 2), n_close(c), n_far(c)] = ...
            local_fitConditionCloseFar(tbl_trials, c, point_params(c, :), opts);
        sigma_point(c, 3) = sigma_point(c, 2) - sigma_point(c, 1);
    end

    sigma_boot = nan(B, num_conds, 3);
    if opts.use_parallel && B > 1
        parfor b = 1:B
            plan_row = local_planRow(boot_plans, b, opts.bootstrap_unit, n_subj, n_trials, opts.boot_plan_seeds);
            sb = local_bootOne(tbl_trials, plan_row, squeeze(params_boot(b, :, :)), opts);
            sigma_boot(b, :, :) = reshape(sb, 1, num_conds, 3);
        end
    else
        for b = 1:B
            plan_row = local_planRow(boot_plans, b, opts.bootstrap_unit, n_subj, n_trials, opts.boot_plan_seeds);
            sb = local_bootOne(tbl_trials, plan_row, squeeze(params_boot(b, :, :)), opts);
            sigma_boot(b, :, :) = reshape(sb, 1, num_conds, 3);
        end
    end

    sigma_jk = nan(n_subj, num_conds, 3);
    if opts.compute_jackknife && ~isempty(jk_params)
        if opts.use_parallel && n_subj > 1
            parfor k = 1:n_subj
                sj = local_jackOne(tbl_trials, subj_list, k, squeeze(jk_params(k, :, :)), opts);
                sigma_jk(k, :, :) = reshape(sj, 1, num_conds, 3);
            end
        else
            for k = 1:n_subj
                sj = local_jackOne(tbl_trials, subj_list, k, squeeze(jk_params(k, :, :)), opts);
                sigma_jk(k, :, :) = reshape(sj, 1, num_conds, 3);
            end
        end
    end

    alpha_lvl = 1 - (opts.ci_prctile(2) - opts.ci_prctile(1)) / 100;
    ci_pct = struct('lo', nan(num_conds, 3), 'hi', nan(num_conds, 3));
    ci_bca = ci_pct;
    ci_bca.z0 = nan(num_conds, 3);
    ci_bca.acc = nan(num_conds, 3);

    for c = 1:num_conds
        adm = admitted(:, c);
        for k = 1:3
            bv = squeeze(sigma_boot(:, c, k));
            bv = bv(adm);
            ci_pct.lo(c, k) = local_percentileFinite(bv, opts.ci_prctile(1));
            ci_pct.hi(c, k) = local_percentileFinite(bv, opts.ci_prctile(2));

            jv = squeeze(sigma_jk(:, c, k));
            [a1, a2, z0, acc] = local_bcaQuantiles(bv, jv, sigma_point(c, k), alpha_lvl);
            [ci_bca.lo(c, k), ci_bca.hi(c, k)] = local_orderedPercentileCI(bv, 100 * a1, 100 * a2);
            ci_bca.z0(c, k) = z0;
            ci_bca.acc(c, k) = acc;
        end
    end

    if strcmpi(opts.ci_method, 'bca') && opts.compute_jackknife
        active = ci_bca;
    else
        active = ci_pct;
    end

    [manip_label, prev_lvl, curr_lvl] = local_conditionLabels(num_conds);
    summary_table = table( ...
        categorical(manip_label), prev_lvl, curr_lvl, n_close, n_far, ...
        sigma_point(:, 1), sigma_point(:, 2), sigma_point(:, 3), ...
        active.lo(:, 1), active.hi(:, 1), ...
        active.lo(:, 2), active.hi(:, 2), ...
        active.lo(:, 3), active.hi(:, 3), ...
        'VariableNames', {'cond_manipulation', 'cond_prev', 'cond_curr', 'n_close', 'n_far', ...
        'sigma_close', 'sigma_far', 'delta_sigma_far_minus_close', ...
        'sigma_close_ci_lo', 'sigma_close_ci_hi', ...
        'sigma_far_ci_lo', 'sigma_far_ci_hi', ...
        'delta_sigma_ci_lo', 'delta_sigma_ci_hi'});

    close_far = struct();
    close_far.threshold_deg = opts.close_threshold_deg;
    close_far.min_trials = opts.min_trials;
    close_far.sigma_point = sigma_point;
    close_far.sigma_boot = sigma_boot;
    close_far.sigma_jackknife = sigma_jk;
    close_far.ci_percentile = ci_pct;
    close_far.ci_bca = ci_bca;
    close_far.ci_active = active;
    close_far.summary_table = summary_table;
end

function opts = local_defaults(opts)
    if ~isfield(opts, 'guess_rate') || isempty(opts.guess_rate), opts.guess_rate = 0.25; end
    if ~isfield(opts, 'close_threshold_deg') || isempty(opts.close_threshold_deg), opts.close_threshold_deg = 30; end
    if ~isfield(opts, 'min_trials') || isempty(opts.min_trials), opts.min_trials = 10; end
    if ~isfield(opts, 'sigma_bounds') || isempty(opts.sigma_bounds), opts.sigma_bounds = [1, 90]; end
    if ~isfield(opts, 'ci_prctile') || isempty(opts.ci_prctile), opts.ci_prctile = [2.5, 97.5]; end
    if ~isfield(opts, 'ci_method') || isempty(opts.ci_method), opts.ci_method = 'bca'; end
    if ~isfield(opts, 'compute_jackknife') || isempty(opts.compute_jackknife), opts.compute_jackknife = true; end
    if ~isfield(opts, 'use_parallel') || isempty(opts.use_parallel), opts.use_parallel = true; end
    if ~isfield(opts, 'bootstrap_unit') || isempty(opts.bootstrap_unit), opts.bootstrap_unit = 'subject'; end
    if ~isfield(opts, 'boot_plan_seeds') || isempty(opts.boot_plan_seeds), opts.boot_plan_seeds = []; end
    if ~isfield(opts, 'point_params') || isempty(opts.point_params)
        error('computeCloseFarSigmaBootstrap:missingPointParams', 'opts.point_params is required.');
    end
    opts.ci_method = lower(char(opts.ci_method));
    opts.bootstrap_unit = lower(char(opts.bootstrap_unit));
end

function sigma_row = local_bootOne(tbl_trials, plan_row, params_row, opts)
    packed = generateBootstrapSample(tbl_trials, plan_row, opts.bootstrap_unit);
    sigma_row = nan(18, 3);
    for c = 1:18
        [m, prev, curr] = conditionSubscriptsFromIndex(c);
        mask = packed.manipulation == m & packed.cond_prev == prev & packed.cond_curr == curr;
        [sg_close, sg_far] = local_fitCloseFarVectors( ...
            packed.delta_theta(mask), packed.x_probe(mask), packed.response(mask), params_row(c, :), opts);
        sigma_row(c, :) = [sg_close, sg_far, sg_far - sg_close];
    end
end

function plan_row = local_planRow(boot_plans, b, bootstrap_unit, n_subj, n_trials, boot_plan_seeds)
    if isempty(boot_plans)
        plan_row = makeBootstrapPlan(bootstrap_unit, n_subj, n_trials, boot_plan_seeds(b));
    else
        plan_row = boot_plans(b, :);
    end
end

function sigma_row = local_jackOne(tbl_trials, subj_list, k, params_row, opts)
    keep_mask = tbl_trials.subject_id ~= subj_list(k);
    tbl_jk = tbl_trials(keep_mask, :);
    sigma_row = nan(18, 3);
    for c = 1:18
        [sg_close, sg_far] = local_fitConditionCloseFar(tbl_jk, c, params_row(c, :), opts);
        sigma_row(c, :) = [sg_close, sg_far, sg_far - sg_close];
    end
end

function [sg_close, sg_far, n_close, n_far] = local_fitConditionCloseFar(tbl_trials, c, pfit, opts)
    [m, prev, curr] = conditionSubscriptsFromIndex(c);
    man = ones(height(tbl_trials), 1);
    man(tbl_trials.cond_manipulation == 'precision') = 2;
    mask = man == m & tbl_trials.cond_prev == prev & tbl_trials.cond_curr == curr;

    dt = tbl_trials.delta_theta(mask);
    xp = tbl_trials.x_probe(mask);
    rv = tbl_trials.response(mask);
    close_mask = abs(dt) < opts.close_threshold_deg;
    n_close = sum(close_mask & (rv == 0 | rv == 1));
    n_far = sum(~close_mask & (rv == 0 | rv == 1));
    [sg_close, sg_far] = local_fitCloseFarVectors(dt, xp, rv, pfit, opts);
end

function [sg_close, sg_far] = local_fitCloseFarVectors(delta_theta, x_probe, response, pfit, opts)
    if any(~isfinite(pfit(:))) || numel(pfit) < 4
        sg_close = NaN;
        sg_far = NaN;
        return
    end

    delta_theta = delta_theta(:);
    x_probe = x_probe(:);
    response = double(response(:));
    valid = response == 0 | response == 1;
    delta_theta = delta_theta(valid);
    x_probe = x_probe(valid);
    response = response(valid);

    close_mask = abs(delta_theta) < opts.close_threshold_deg;
    sg_close = local_fitSigmaScalar(delta_theta(close_mask), x_probe(close_mask), response(close_mask), pfit, opts);
    sg_far = local_fitSigmaScalar(delta_theta(~close_mask), x_probe(~close_mask), response(~close_mask), pfit, opts);
end

function sigma_hat = local_fitSigmaScalar(delta_theta, x_probe, response, pfit, opts)
    if numel(response) < opts.min_trials || numel(unique(response)) < 2
        sigma_hat = NaN;
        return
    end

    mu = dogIsolated(delta_theta, pfit(1), pfit(2)) + pfit(4);
    eps_log = eps;
    nll = @(sigma) -sum( ...
        response .* log(((1 - opts.guess_rate) .* normcdf(x_probe, mu, sigma) + 0.5 .* opts.guess_rate) + eps_log) + ...
        (1 - response) .* log(1 - ((1 - opts.guess_rate) .* normcdf(x_probe, mu, sigma) + 0.5 .* opts.guess_rate) + eps_log), ...
        'omitnan');

    try
        sigma_hat = fminbnd(nll, opts.sigma_bounds(1), opts.sigma_bounds(2));
        if ~isfinite(sigma_hat)
            sigma_hat = NaN;
        end
    catch
        sigma_hat = NaN;
    end
end

function [manip_label, prev_lvl, curr_lvl] = local_conditionLabels(num_conds)
    manip_label = cell(num_conds, 1);
    prev_lvl = zeros(num_conds, 1);
    curr_lvl = zeros(num_conds, 1);
    for c = 1:num_conds
        [m, prev_lvl(c), curr_lvl(c)] = conditionSubscriptsFromIndex(c);
        if m == 1
            manip_label{c} = 'contrast';
        else
            manip_label{c} = 'precision';
        end
    end
end

function [a1, a2, z0, a_acc] = local_bcaQuantiles(boot_vals, jk_vals, point_est, alpha)
    boot_vals = boot_vals(isfinite(boot_vals));
    jk_vals = jk_vals(isfinite(jk_vals));
    a1 = alpha / 2; a2 = 1 - alpha / 2; z0 = NaN; a_acc = NaN;
    if numel(boot_vals) < 2 || numel(jk_vals) < 2 || ~isfinite(point_est)
        return
    end

    p_hat = mean(boot_vals < point_est);
    p_hat = min(max(p_hat, 1 / (numel(boot_vals) + 1)), 1 - 1 / (numel(boot_vals) + 1));
    z0 = norminv(p_hat);

    jk_mean = mean(jk_vals);
    num_a = sum((jk_mean - jk_vals).^3);
    den_a = 6 * (sum((jk_mean - jk_vals).^2)).^(3/2);
    if den_a == 0
        a_acc = 0;
    else
        a_acc = num_a / den_a;
    end

    z_lo = norminv(alpha / 2);
    z_hi = norminv(1 - alpha / 2);
    a1 = normcdf(z0 + (z0 + z_lo) / (1 - a_acc * (z0 + z_lo)));
    a2 = normcdf(z0 + (z0 + z_hi) / (1 - a_acc * (z0 + z_hi)));

    if ~isfinite(a1) || a1 <= 0, a1 = alpha / 2; end
    if ~isfinite(a2) || a2 >= 1, a2 = 1 - alpha / 2; end
end

function q = local_percentileFinite(values, prc)
    values = values(isfinite(values));
    if isempty(values) || any(~isfinite(prc(:)))
        q = NaN;
    else
        q = prctile(values, prc);
    end
end

function [lo, hi] = local_orderedPercentileCI(values, prc_lo, prc_hi)
    q1 = local_percentileFinite(values, prc_lo);
    q2 = local_percentileFinite(values, prc_hi);
    lo = min(q1, q2);
    hi = max(q1, q2);
end
