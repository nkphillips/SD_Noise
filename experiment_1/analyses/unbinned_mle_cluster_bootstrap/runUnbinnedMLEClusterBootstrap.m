function results = runUnbinnedMLEClusterBootstrap(tbl_trials, varargin)
% runUnbinnedMLEClusterBootstrap  B-fold subject-cluster bootstrap, unbinned trial-level
% MLE per condition (S&S 2022 DoG + Gaussian psychometric + 25% guess rate).
% Returns BOTH percentile and BCa intervals for [A, w, sigma, beta], FWHM,
% and the smooth mu(delta_theta) = DoG + beta curve. Active CI is routed
% from ci_method (default 'bca').
%
% Name-value pairs:
%   'B' (2000), 'seed' (1), 'B_grid' (100), 'use_parallel' (true), 'num_workers' ([])
%   'fit_opts' (struct passed to fitConditionMLE; e.g. .lb, .ub, .x0, .guess_rate)
%   'ci_prctile' ([2.5 97.5])
%   'ci_method' ('bca' | 'percentile') -- routes active sd_ci_cluster.lo/hi
%   'compute_jackknife' (true) -- BCa requires leave-one-subject-out jackknife
%   'discard_at_bound' (true), 'discard_failed_fits' (true), 'bound_tol' (1e-4)
%   'guess_rate' (0.25) -- S&S 2022; passed through fit_opts
%   'contrast_labels', 'precision_labels' -- 1x3 cells for subplot titles
%   'mu_bounds' -- [lower upper] degrees for plot y-axis (default [-20 20])
%   'fig_subdir' -- e.g. '1_back'; figures and contrasts CSV go to figures/<subdir>/
%   'compute_contrasts' (true) -- auto-compute within-manipulation BCa contrasts;
%                                 results.contrast_table is added on success.
%   'contrast_params' ({'A'}) -- subset of {'A','w','sigma','beta'} to test contrasts on.
%   'contrast_specs' ([]) -- override specs; if empty, uses buildStandardContrasts.
%   'subject_influence' (true) -- per-subject leverage diagnostic from the jackknife.
%   'amplitude_sigma_correlation' (true) -- per-subject DoG fits + S&S Fig 1H scatter.
%   'demean_x_probe' (false) -- if true, subtract a per-subject Gaussian-psychometric
%                                baseline mu_i from x_probe before all SD fits.
%                                Computed once on the full data and reused across
%                                bootstrap and jackknife replicates. Saves a
%                                diagnostic figure and CSV under subject/baseline_bias/.
%   'skip_at_bound_baseline' (true) -- if true, subjects whose baseline sigma_i
%                                lands at either bound of the [1, 90] deg search
%                                range are NOT demeaned (mu_i = 0 used instead).
%                                Their trials still enter the bootstrap pool.
%                                Set to false to demean every subject regardless
%                                of fit quality.
%   'subject_data_quality' (true) -- per-subject diagnostic figure with trial
%                                counts per cell, empirical psychometric, and
%                                summary stats. Saved under subject/data_quality/.

    ip = inputParser;
    addParameter(ip, 'B', 2000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(ip, 'seed', 1, @(x) isnumeric(x) && isscalar(x));
    addParameter(ip, 'B_grid', 100, @(x) isnumeric(x) && isscalar(x) && x >= 5);
    addParameter(ip, 'use_parallel', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'num_workers', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'fit_opts', struct(), @isstruct);
    addParameter(ip, 'ci_prctile', [2.5, 97.5], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(ip, 'ci_method', 'bca', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'compute_jackknife', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'discard_at_bound', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'discard_failed_fits', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'bound_tol', 1e-4, @(x) isnumeric(x) && isscalar(x));
    addParameter(ip, 'guess_rate', 0.25, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x < 1);
    addParameter(ip, 'contrast_labels', {'90%', '50%', '25%'}, @(x) iscell(x) && numel(x) == 3);
    addParameter(ip, 'precision_labels', {'2°', '40°', '80°'}, @(x) iscell(x) && numel(x) == 3);
    addParameter(ip, 'mu_bounds', [-20, 20], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(ip, 'fig_subdir', '', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'compute_contrasts', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'contrast_params', {'A'}, @(x) iscell(x));   % which DoG params to test contrasts on
    addParameter(ip, 'contrast_specs', [], @(x) isempty(x) || isstruct(x));   % override; empty => use buildStandardContrasts
    addParameter(ip, 'subject_influence', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'amplitude_sigma_correlation', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'demean_x_probe', false, @(x) islogical(x) && isscalar(x));   % per-subject baseline demean before fitting
    addParameter(ip, 'skip_at_bound_baseline', true, @(x) islogical(x) && isscalar(x));  % skip demean for subjects whose sigma_baseline saturates
    addParameter(ip, 'subject_data_quality', true, @(x) islogical(x) && isscalar(x));   % per-subject trial counts + empirical psychometric figure
    addParameter(ip, 'subj_labels', {}, @(x) iscell(x));   % human-readable subject IDs for diagnostics
    parse(ip, varargin{:});

    B = ip.Results.B;
    seed = ip.Results.seed;
    n_grid = ip.Results.B_grid;
    use_parallel = ip.Results.use_parallel;
    num_workers = ip.Results.num_workers;
    fit_opts = ip.Results.fit_opts;
    ci_pct = ip.Results.ci_prctile;
    ci_method = lower(char(ip.Results.ci_method));
    compute_jackknife = ip.Results.compute_jackknife;
    discard_at_bound = ip.Results.discard_at_bound;
    discard_failed_fits = ip.Results.discard_failed_fits;
    bound_tol = ip.Results.bound_tol;
    guess_rate = ip.Results.guess_rate;

    fit_opts.guess_rate = guess_rate;

    if ~ismember(ci_method, {'bca', 'percentile'})
        error('runUnbinnedMLEClusterBootstrap:badCIMethod', ...
            'ci_method must be ''bca'' or ''percentile''.');
    end

    this_dir = fileparts(mfilename('fullpath'));
    fig_subdir = char(ip.Results.fig_subdir);
    if isempty(fig_subdir)
        fig_dir = fullfile(this_dir, 'figures');
    else
        fig_dir = fullfile(this_dir, 'figures', fig_subdir);
    end
    super_subj_dir = fullfile(fig_dir, 'super_subject');
    subj_dir = fullfile(fig_dir, 'subject');
    for d = {fig_dir, super_subj_dir, subj_dir}
        if ~exist(d{1}, 'dir'); mkdir(d{1}); end
    end

    addpath(this_dir);
    addpath(fullfile(this_dir, '..', 'functions'));

    % Capture the pre-demean trial table for the data-quality diagnostic
    tbl_trials_raw = tbl_trials;

    % -------- Optional per-subject demean of x_probe (D1-D8 confirmed plan) --------
    %   Compute mu_i once on the full data using fitSubjectBaselineBias (Gaussian
    %   psychometric, no DoG, ignoring Delta-theta), then subtract mu_i from
    %   x_probe per trial. The same mu_i is reused across all bootstrap and
    %   jackknife replicates. beta remains free in the SD model and is expected
    %   to fit near zero after demean.
    %
    %   When skip_at_bound_baseline = true (default), subjects whose baseline
    %   sigma_i saturates at its upper or lower bound are NOT demeaned (mu_i = 0
    %   used instead). Their data stays in the bootstrap pool with x_probe raw.
    baseline_bias_table = table();
    if ip.Results.demean_x_probe
        try
            baseline_fit_opts = struct('guess_rate', guess_rate);
            [tbl_trials, baseline_bias_table] = demeanTrialTable( ...
                tbl_trials_raw, baseline_fit_opts, ip.Results.skip_at_bound_baseline);

            bb_dir = fullfile(subj_dir, 'baseline_bias');
            plotSubjectBaselineBias(baseline_bias_table, bb_dir, ip.Results.subj_labels);
        catch demeanErr
            warning('runUnbinnedMLEClusterBootstrap:demeanFailed', ...
                'Demean step failed: %s. Proceeding with un-demeaned x_probe.', demeanErr.message);
            baseline_bias_table = table();
            tbl_trials = tbl_trials_raw;
        end
    end

    % -------- Per-subject data quality diagnostic (uses the raw, pre-demean table) --------
    if ip.Results.subject_data_quality
        try
            dq_dir = fullfile(subj_dir, 'data_quality');
            plotSubjectDataQuality(tbl_trials_raw, dq_dir, ip.Results.subj_labels, baseline_bias_table);
        catch dqErr
            warning('runUnbinnedMLEClusterBootstrap:dataQualityFailed', ...
                'Subject data quality figure failed: %s.', dqErr.message);
        end
    end

    subj_list = unique(tbl_trials.subject_id, 'stable');
    n_subj = numel(subj_list);
    if n_subj < 1
        error('runUnbinnedMLEClusterBootstrap:noSubjects', 'tbl_trials has no subjects.');
    end

    rng(seed, 'twister');
    boot_plans = zeros(B, n_subj);
    for b = 1:B
        boot_plans(b, :) = randi(n_subj, [1, n_subj]);
    end

    grid = linspace(-90, 90, n_grid)';
    num_conds = 18;

    % -------- Default bounds (must match fitConditionMLE defaults) --------
    fwhm_min_deg = 8;
    fwhm_max_deg = 120;
    w_lb_def = (2 * sqrt(log(2))) / fwhm_max_deg;
    w_ub_def = (2 * sqrt(log(2))) / fwhm_min_deg;
    lb_default = [-30; w_lb_def;   1; -10];   % beta_lb widened to -10 (was -5)
    ub_default = [ 30; w_ub_def;  90;  10];

    if isfield(fit_opts, 'lb') && ~isempty(fit_opts.lb)
        lb_used = fit_opts.lb(:);
    else
        lb_used = lb_default;
    end
    if isfield(fit_opts, 'ub') && ~isempty(fit_opts.ub)
        ub_used = fit_opts.ub(:);
    else
        ub_used = ub_default;
    end

    % -------- Bootstrap replicates (parfor over B) --------
    params_boot   = nan(B, num_conds, 4);
    curve_boot    = nan(B, num_conds, n_grid);   % isolated DoG only (no beta)
    exit_boot     = nan(B, num_conds);
    at_bound_boot = false(B, num_conds, 4);

    if use_parallel && B > 1
        pool = gcp('nocreate');
        if isempty(pool)
            if isempty(num_workers)
                parpool('local');
            else
                parpool('local', num_workers);
            end
        end

        t_bs = tic;
        parfor b = 1:B
            [pr, cr, er, ar] = bootstrapOneReplicate(tbl_trials, boot_plans(b, :), grid, fit_opts, lb_used, ub_used, bound_tol);
            params_boot(b, :, :)   = reshape(pr, 1, num_conds, 4);
            curve_boot(b, :, :)    = reshape(cr, 1, num_conds, n_grid);
            exit_boot(b, :)        = er(:)';
            at_bound_boot(b, :, :) = reshape(ar, 1, num_conds, 4);
        end
        boot_sec = toc(t_bs);
    else
        t_bs = tic;
        for b = 1:B
            [pr, cr, er, ar] = bootstrapOneReplicate(tbl_trials, boot_plans(b, :), grid, fit_opts, lb_used, ub_used, bound_tol);
            params_boot(b, :, :)   = reshape(pr, 1, num_conds, 4);
            curve_boot(b, :, :)    = reshape(cr, 1, num_conds, n_grid);
            exit_boot(b, :)        = er(:)';
            at_bound_boot(b, :, :) = reshape(ar, 1, num_conds, 4);
        end
        boot_sec = toc(t_bs);
    end

    % -------- Admission filtering --------
    admitted = true(B, num_conds);
    discard_reason = strings(B, num_conds);
    for b = 1:B
        for c = 1:num_conds
            pf = squeeze(params_boot(b, c, :));
            if all(isnan(pf))
                admitted(b, c) = false;
                discard_reason(b, c) = "failed_fit";
                continue
            end
            if discard_failed_fits && ~isnan(exit_boot(b, c)) && exit_boot(b, c) <= 0
                admitted(b, c) = false;
                discard_reason(b, c) = "failed_fit";
                continue
            end
            if discard_at_bound && any(squeeze(at_bound_boot(b, c, :)))
                admitted(b, c) = false;
                discard_reason(b, c) = "at_bound";
                continue
            end
        end
    end

    % -------- Pooled point estimate + figure-overlay metrics --------
    delta_centers_emp = -90:2:90;
    empirical_window_width = 32;
    overlay = computePooledFitsAndOverlay(tbl_trials, delta_centers_emp, empirical_window_width, guess_rate);

    % -------- Optional leave-one-subject-out jackknife (for BCa acceleration) --------
    jk = struct('params', nan(n_subj, num_conds, 4), ...
                'curves', nan(n_subj, num_conds, n_grid), ...
                'fwhm',   nan(n_subj, num_conds), ...
                'duration_sec', 0);
    if compute_jackknife
        t_jk = tic;
        jk_use_parallel = use_parallel && n_subj > 1;
        jk = runSubjectJackknife(tbl_trials, subj_list, grid, fit_opts, jk_use_parallel);
        jk.duration_sec = toc(t_jk);
    elseif strcmpi(ci_method, 'bca')
        warning('runUnbinnedMLEClusterBootstrap:noJackknife', ...
            'ci_method = ''bca'' but compute_jackknife = false. Falling back to percentile CI.');
        ci_method = 'percentile';
    end

    % -------- CI computation: percentile + BCa per cond --------
    alpha_lvl = 1 - (ci_pct(2) - ci_pct(1)) / 100;        % e.g. 0.05 for [2.5,97.5]

    % Per-parameter point estimates from pooled fit.
    pp = overlay.params_point;                            % [num_conds x 4]

    % Storage
    ci_perc = struct('param_lo', nan(num_conds, 4), 'param_hi', nan(num_conds, 4), ...
                     'fwhm_lo',  nan(num_conds, 1), 'fwhm_hi',  nan(num_conds, 1), ...
                     'curve_lo', nan(num_conds, n_grid), 'curve_hi', nan(num_conds, n_grid));
    ci_bca  = ci_perc;
    ci_bca.param_z0  = nan(num_conds, 4);
    ci_bca.param_acc = nan(num_conds, 4);

    % Per-parameter and FWHM
    for c = 1:num_conds
        adm = admitted(:, c);

        % Per-parameter
        for k = 1:4
            v = squeeze(params_boot(:, c, k));
            v = v(adm);
            ci_perc.param_lo(c, k) = percentileFinite(v, ci_pct(1));
            ci_perc.param_hi(c, k) = percentileFinite(v, ci_pct(2));

            jk_vals = squeeze(jk.params(:, c, k));
            theta_hat = pp(c, k);
            [a1, a2, z0, a_acc] = bcaQuantiles(v, jk_vals, theta_hat, alpha_lvl);
            ci_bca.param_lo(c, k)  = percentileFinite(v, 100 * a1);
            ci_bca.param_hi(c, k)  = percentileFinite(v, 100 * a2);
            ci_bca.param_z0(c, k)  = z0;
            ci_bca.param_acc(c, k) = a_acc;
        end

        % FWHM (derived from w via Gaussian-envelope mapping; monotonic in 1/w)
        w_vals = squeeze(params_boot(:, c, 2));
        w_vals = w_vals(adm);
        fwhm_vals = unbinnedWtoFwhm(w_vals);

        ci_perc.fwhm_lo(c) = percentileFinite(fwhm_vals, ci_pct(1));
        ci_perc.fwhm_hi(c) = percentileFinite(fwhm_vals, ci_pct(2));

        if isfinite(pp(c, 2)) && pp(c, 2) > 0
            fwhm_pt = unbinnedWtoFwhm(pp(c, 2));
        else
            fwhm_pt = NaN;
        end
        jk_w = squeeze(jk.params(:, c, 2));
        jk_fwhm = unbinnedWtoFwhm(jk_w);
        [a1, a2] = bcaQuantiles(fwhm_vals, jk_fwhm, fwhm_pt, alpha_lvl);
        ci_bca.fwhm_lo(c) = percentileFinite(fwhm_vals, 100 * a1);
        ci_bca.fwhm_hi(c) = percentileFinite(fwhm_vals, 100 * a2);

        % Smooth curve mu(delta_theta) = DoG(grid; A, w) + beta, pointwise CI
        beta_vals = squeeze(params_boot(:, c, 4));
        % full mu curve per replicate
        full_curve = squeeze(curve_boot(:, c, :)) + beta_vals;
        full_curve = full_curve(adm, :);

        if isfinite(pp(c, 1)) && isfinite(pp(c, 2)) && isfinite(pp(c, 4))
            curve_pt = (dogIsolated(grid, pp(c, 1), pp(c, 2)) + pp(c, 4))';
        else
            curve_pt = nan(1, n_grid);
        end

        % jackknife curve per held-out subject (DoG only) + beta
        if all(isfinite(pp(c, :)))
            jk_iso = squeeze(jk.curves(:, c, :));         % n_subj x n_grid
            jk_beta = squeeze(jk.params(:, c, 4));        % n_subj x 1
            jk_full_curve = jk_iso + jk_beta;
        else
            jk_full_curve = nan(n_subj, n_grid);
        end

        for ix = 1:n_grid
            cv = full_curve(:, ix);
            ci_perc.curve_lo(c, ix) = percentileFinite(cv, ci_pct(1));
            ci_perc.curve_hi(c, ix) = percentileFinite(cv, ci_pct(2));

            jkv = jk_full_curve(:, ix);
            [a1, a2] = bcaQuantiles(cv, jkv, curve_pt(ix), alpha_lvl);
            ci_bca.curve_lo(c, ix) = percentileFinite(cv, 100 * a1);
            ci_bca.curve_hi(c, ix) = percentileFinite(cv, 100 * a2);
        end
    end

    % -------- Route active CI --------
    if strcmpi(ci_method, 'bca')
        active = ci_bca;
    else
        active = ci_perc;
    end

    % -------- Per-cond medians (admitted only, finite-percentile-style) --------
    alpha_pt  = nan(num_conds, 1);
    w_pt_vec  = nan(num_conds, 1);
    sigma_pt  = nan(num_conds, 1);
    beta_pt   = nan(num_conds, 1);
    for c = 1:num_conds
        adm = admitted(:, c);
        alpha_pt(c) = medianFinite(params_boot(adm, c, 1));
        w_pt_vec(c) = medianFinite(params_boot(adm, c, 2));
        sigma_pt(c) = medianFinite(params_boot(adm, c, 3));
        beta_pt(c)  = medianFinite(params_boot(adm, c, 4));
    end

    % -------- Subscript labels --------
    manip = zeros(num_conds, 1);
    prev_lvl = zeros(num_conds, 1);
    curr_lvl = zeros(num_conds, 1);
    manip_label = cell(num_conds, 1);
    for c = 1:num_conds
        [manip(c), prev_lvl(c), curr_lvl(c)] = conditionSubscriptsFromIndex(c);
        if manip(c) == 1
            manip_label{c} = 'contrast';
        else
            manip_label{c} = 'precision';
        end
    end

    % -------- Summary table (active CI lo/hi for plot-script consumption) --------
    summary_table = table( ...
        categorical(manip_label), prev_lvl, curr_lvl, ...
        alpha_pt, w_pt_vec, sigma_pt, beta_pt, ...
        active.param_lo(:, 1), active.param_hi(:, 1), ...
        active.param_lo(:, 2), active.param_hi(:, 2), ...
        active.param_lo(:, 3), active.param_hi(:, 3), ...
        active.param_lo(:, 4), active.param_hi(:, 4), ...
        active.fwhm_lo,        active.fwhm_hi, ...
        'VariableNames', {'cond_manipulation', 'cond_prev', 'cond_curr', ...
        'alpha_median', 'w_median', 'sigma_median', 'beta_median', ...
        'alpha_prctile_lo', 'alpha_prctile_hi', 'w_prctile_lo', 'w_prctile_hi', ...
        'sigma_prctile_lo', 'sigma_prctile_hi', 'beta_prctile_lo', 'beta_prctile_hi', ...
        'fwhm_prctile_lo',  'fwhm_prctile_hi'});

    % -------- Discard summary --------
    discard_summary = struct();
    discard_summary.B_total            = B;
    discard_summary.admitted_per_cond  = sum(admitted, 1)';                  % num_conds x 1
    discard_summary.discarded_per_cond = sum(~admitted, 1)';
    discard_summary.global_admitted    = sum(admitted(:));
    discard_summary.global_discarded   = sum(~admitted(:));
    discard_summary.by_reason = struct( ...
        'failed_fit', sum(discard_reason(:) == "failed_fit"), ...
        'at_bound',   sum(discard_reason(:) == "at_bound"));

    % -------- Assemble results --------
    results = struct();
    results.params_boot      = params_boot;
    results.curve_boot       = curve_boot;
    results.exit_boot        = exit_boot;
    results.at_bound_boot    = at_bound_boot;
    results.admitted         = admitted;
    results.discard_reason   = discard_reason;
    results.discard_summary  = discard_summary;
    results.grid             = grid;
    results.B                = B;
    results.seed             = seed;
    results.n_subj           = n_subj;
    results.boot_plans       = boot_plans;
    results.ci_prctile       = ci_pct;
    results.ci_method        = ci_method;
    results.compute_jackknife = compute_jackknife;
    results.guess_rate       = guess_rate;
    results.demean_x_probe   = ip.Results.demean_x_probe;
    results.skip_at_bound_baseline = ip.Results.skip_at_bound_baseline;
    results.baseline_bias    = baseline_bias_table;
    results.lb               = lb_used(:)';
    results.ub               = ub_used(:)';
    results.bound_tol        = bound_tol;
    results.discard_at_bound = discard_at_bound;
    results.discard_failed_fits = discard_failed_fits;
    results.timing = struct('bootstrap_replicates_sec', boot_sec, ...
                            'jackknife_sec', jk.duration_sec);
    results.overlay  = overlay;
    results.empirical_delta_centers = delta_centers_emp;
    results.empirical_window_width  = empirical_window_width;
    results.jackknife = jk;
    results.ci_percentile = ci_perc;
    results.ci_bca = ci_bca;
    results.ci_active = active;
    results.summary_table = summary_table;

    % -------- Within-manipulation contrasts (BCa CIs on linear-combination contrasts) --------
    if ip.Results.compute_contrasts && compute_jackknife
        if isempty(ip.Results.contrast_specs)
            cspecs = buildStandardContrasts('params', ip.Results.contrast_params);
        else
            cspecs = ip.Results.contrast_specs;
        end
        try
            results.contrast_table = computeUnbinnedContrasts(results, cspecs);
            % Save CSV alongside the super-subject figures (it's a super-subject output)
            csv_path = fullfile(super_subj_dir, 'contrasts_bca.csv');
            try
                writetable(results.contrast_table, csv_path);
            catch %#ok<CTCH>
                % non-fatal -- the table is still in results
            end
        catch contrastErr
            warning('runUnbinnedMLEClusterBootstrap:contrastFailed', ...
                'Within-manipulation contrast computation failed: %s', contrastErr.message);
            results.contrast_table = table();
        end
    else
        if ip.Results.compute_contrasts && ~compute_jackknife
            warning('runUnbinnedMLEClusterBootstrap:noJackknifeForContrasts', ...
                'compute_contrasts requested but compute_jackknife = false; skipping (BCa requires jackknife).');
        end
        results.contrast_table = table();
    end

    % -------- Subject-influence diagnostic (uses leave-one-subject-out jackknife) --------
    if ip.Results.subject_influence
        if ~compute_jackknife
            warning('runUnbinnedMLEClusterBootstrap:noJackknifeForInfluence', ...
                'subject_influence requested but compute_jackknife = false; skipping.');
            results.subject_influence = struct();
        else
            try
                influence = computeSubjectInfluence(results);
                results.subject_influence = influence;
                inf_dir = fullfile(subj_dir, 'subject_influence');
                plotSubjectInfluence(influence, inf_dir, ip.Results.subj_labels);
            catch infErr
                warning('runUnbinnedMLEClusterBootstrap:subjectInfluenceFailed', ...
                    'Subject-influence diagnostic failed: %s', infErr.message);
                results.subject_influence = struct();
            end
        end
    end

    % -------- Per-subject per-manipulation DoG fits + S&S Fig 1H analog --------
    if ip.Results.amplitude_sigma_correlation
        try
            per_subj_fits = fitPerSubjectPerManipulation(tbl_trials, fit_opts, ip.Results.subj_labels);
            results.per_subject_per_manip_fits = per_subj_fits;
            corr_dir = fullfile(subj_dir, 'amplitude_sigma_correlation');
            r_summary = plotAmplitudeSigmaCorrelation(per_subj_fits, corr_dir);
            results.amplitude_sigma_correlation = r_summary;
        catch corrErr
            warning('runUnbinnedMLEClusterBootstrap:amplitudeSigmaFailed', ...
                'Amplitude vs sigma correlation diagnostic failed: %s', corrErr.message);
            results.amplitude_sigma_correlation = struct();
        end
    end

    % -------- Plots --------
    ps = plotSettings();
    plotDoGMLEBootstrapFigures(ps, super_subj_dir, curve_boot, grid, params_boot, ci_pct, ...
        'contrast_labels', ip.Results.contrast_labels, ...
        'precision_labels', ip.Results.precision_labels, ...
        'mu_bounds', ip.Results.mu_bounds, ...
        'overlay', overlay, ...
        'admitted', admitted, ...
        'curve_lo', active.curve_lo, ...
        'curve_hi', active.curve_hi, ...
        'ci_method', ci_method);

    plotUnbinnedSdScatterSummaries(super_subj_dir, summary_table, ...
        ip.Results.contrast_labels, ip.Results.precision_labels, ps, ci_pct);

end

% =========================================================================
% Helpers
% =========================================================================

function jk = runSubjectJackknife(tbl_trials, subj_list, grid, fit_opts, use_parallel)
    n_subj = numel(subj_list);
    n_grid = numel(grid);
    num_conds = 18;

    sd_p_cell  = cell(n_subj, 1);
    curve_cell = cell(n_subj, 1);

    if use_parallel
        parfor k = 1:n_subj
            [sp, cc] = jackknifeOneSubject(k, subj_list, tbl_trials, grid, fit_opts);
            sd_p_cell{k}  = sp;
            curve_cell{k} = cc;
        end
    else
        for k = 1:n_subj
            [sp, cc] = jackknifeOneSubject(k, subj_list, tbl_trials, grid, fit_opts);
            sd_p_cell{k}  = sp;
            curve_cell{k} = cc;
        end
    end

    jk = struct();
    jk.params = nan(n_subj, num_conds, 4);
    jk.curves = nan(n_subj, num_conds, n_grid);
    jk.fwhm   = nan(n_subj, num_conds);
    for k = 1:n_subj
        if ~isempty(sd_p_cell{k})
            jk.params(k, :, :) = sd_p_cell{k};
        end
        if ~isempty(curve_cell{k})
            jk.curves(k, :, :) = curve_cell{k};
        end
        for c = 1:num_conds
            wk = jk.params(k, c, 2);
            if isfinite(wk) && wk > 0
                jk.fwhm(k, c) = unbinnedWtoFwhm(wk);
            end
        end
    end
end

function [sd_params, curves] = jackknifeOneSubject(k, subj_list, tbl_trials, grid, fit_opts)
    num_conds = 18;
    n_grid = numel(grid);
    sd_params = nan(num_conds, 4);
    curves    = nan(num_conds, n_grid);

    keep_id = subj_list(k);
    keep_mask = tbl_trials.subject_id ~= keep_id;
    tbl_jk = tbl_trials(keep_mask, :);

    cm = tbl_jk.cond_manipulation;
    man = ones(height(tbl_jk), 1);
    man(cm == 'precision') = 2;

    for c = 1:num_conds
        [m, prev, curr] = conditionSubscriptsFromIndex(c);
        mask = man == m & tbl_jk.cond_prev == prev & tbl_jk.cond_curr == curr;
        pf = fitConditionMLE(tbl_jk.delta_theta(mask), tbl_jk.x_probe(mask), tbl_jk.response(mask), fit_opts);
        sd_params(c, :) = pf;
        if ~all(isnan(pf))
            curves(c, :) = dogIsolated(grid, pf(1), pf(2))';
        end
    end
end

function [a1, a2, z0, a_acc] = bcaQuantiles(boot_vals, jk_vals, point_est, alpha)
    boot_vals = boot_vals(isfinite(boot_vals));
    jk_vals   = jk_vals(isfinite(jk_vals));
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

function q = percentileFinite(values, prc)
    values = values(isfinite(values));
    if isempty(values) || any(~isfinite(prc(:)))
        q = NaN;
    else
        q = prctile(values, prc);
    end
end

function m = medianFinite(values)
    values = values(isfinite(values));
    if isempty(values)
        m = NaN;
    else
        m = median(values);
    end
end
