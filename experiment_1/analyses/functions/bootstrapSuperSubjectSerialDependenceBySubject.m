function [sd_ci_cluster, sd_boot_cluster] = bootstrapSuperSubjectSerialDependenceBySubject(delta_theta_windows, delta_theta_centers, num, p, bootstrap, toggles, sd_point, trial_pool)
% bootstrapSuperSubjectSerialDependenceBySubject
%
% Subject-cluster bootstrap CIs for pooled super-subject DoG fits. Each
% bootstrap replicate resamples subjects with replacement, rebuilds a pooled
% super-subject data set from the sampled subjects, and reruns the same
% response-bias-to-serial-dependence fitting path used for the point estimate.
%
% This function returns BOTH percentile and BCa intervals. The active CI
% fields (sd_ci_cluster.lo/hi/fwhm_lo/fwhm_hi/curve_lo/curve_hi) are set
% from sd_ci_cluster.(bootstrap.ci_method).* so downstream plotting code is
% unchanged.
%
% Inputs:
%   delta_theta_windows - per-subject windowed pool from makeDeltaThetaWindows.
%   delta_theta_centers - 1xW vector of window centers (deg).
%   num                 - struct with .levels, .conds, .subjs, .delta_theta_windows.
%   p                   - parameter struct (sd_bounds, sd_objective, sd_min_windows,
%                         sd_lobe_edges, fmincon_options, ...).
%   bootstrap           - settings struct: B_subject_cluster_sd, ci, ci_method,
%                         compute_jackknife, discard_at_bound, discard_failed_fits,
%                         bound_tol,
%                         subject_cluster_seed, subject_cluster_num_chunks,
%                         subject_cluster_curve_x.
%   toggles             - .disp_on, .parallelization.
%   sd_point            - point-estimate sd struct (sd.all.params_est used as
%                         theta_hat for BCa bias correction). When omitted or
%                         empty, the function falls back to percentile for
%                         sd_ci_cluster.lo/hi.
%   trial_pool          - REQUIRED for sd_objective='nll'. Per-subject un-windowed
%                         trial pool from buildPerSubjectCondTrialPool(all_runs, ...).
%                         SSE mode ignores it.

    if nargin < 6
        toggles = struct('disp_on', 0, 'parallelization', 0);
    end
    if nargin < 7
        sd_point = [];
    end
    if nargin < 8
        trial_pool = [];
    end

    if isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'nll') && isempty(trial_pool)
        error('bootstrapSuperSubjectSerialDependenceBySubject:missingTrialPool', ...
            ['NLL mode requires the un-windowed trial_pool argument ' ...
             '(buildPerSubjectCondTrialPool). Pass it as the 8th argument or run SSE.']);
    end

    % --- Defaults for new bootstrap fields ---
    if ~isfield(bootstrap, 'B_subject_cluster_sd') || isempty(bootstrap.B_subject_cluster_sd)
        bootstrap.B_subject_cluster_sd = 2000;
    end
    if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
        bootstrap.ci = [2.5, 97.5];
    end
    if ~isfield(bootstrap, 'ci_method') || isempty(bootstrap.ci_method)
        bootstrap.ci_method = 'bca';
    end
    if ~isfield(bootstrap, 'compute_jackknife') || isempty(bootstrap.compute_jackknife)
        bootstrap.compute_jackknife = strcmpi(bootstrap.ci_method, 'bca');
    end
    if ~isfield(bootstrap, 'discard_at_bound') || isempty(bootstrap.discard_at_bound)
        bootstrap.discard_at_bound = true;
    end
    if ~isfield(bootstrap, 'discard_failed_fits') || isempty(bootstrap.discard_failed_fits)
        bootstrap.discard_failed_fits = true;
    end
    if ~isfield(bootstrap, 'bound_tol') || isempty(bootstrap.bound_tol)
        bootstrap.bound_tol = 1e-4;
    end
    if ~isfield(bootstrap, 'subject_cluster_num_chunks') || isempty(bootstrap.subject_cluster_num_chunks)
        if isfield(p, 'num_chunks') && ~isempty(p.num_chunks)
            bootstrap.subject_cluster_num_chunks = min(p.num_chunks, bootstrap.B_subject_cluster_sd);
        else
            bootstrap.subject_cluster_num_chunks = min(1, bootstrap.B_subject_cluster_sd);
        end
    end
    if ~isfield(bootstrap, 'subject_cluster_curve_x') || isempty(bootstrap.subject_cluster_curve_x)
        bootstrap.subject_cluster_curve_x = linspace(-90, 90, 100);
    end

    B = bootstrap.B_subject_cluster_sd;
    prc = bootstrap.ci;
    % alpha = 1 - (CI width). e.g., prc = [2.5 97.5] -> alpha = 0.05
    alpha = 1 - (prc(2) - prc(1)) / 100;
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    curve_x = bootstrap.subject_cluster_curve_x(:)';
    n_curve = numel(curve_x);

    % --- Storage allocation ---
    sd_boot_cluster = struct();
    sd_boot_cluster.params         = nan(B, num.levels, num.levels, num.conds, num_params);
    sd_boot_cluster.fwhm           = nan(B, num.levels, num.levels, num.conds);
    sd_boot_cluster.curves         = nan(B, num.levels, num.levels, num.conds, n_curve);
    sd_boot_cluster.exit_flag      = nan(B, num.levels, num.levels, num.conds);
    sd_boot_cluster.at_bound       = false(B, num.levels, num.levels, num.conds, num_params);
    sd_boot_cluster.admitted       = false(B, num.levels, num.levels, num.conds);
    sd_boot_cluster.discard_reason = strings(B, num.levels, num.levels, num.conds);
    sd_boot_cluster.curve_x        = curve_x;
    sd_boot_cluster.B              = B;
    sd_boot_cluster.ci             = prc;
    sd_boot_cluster.ci_method         = bootstrap.ci_method;
    sd_boot_cluster.compute_jackknife = bootstrap.compute_jackknife;
    sd_boot_cluster.method            = 'subject_cluster_pooled_super_subject';

    sd_ci_cluster = struct();
    sd_ci_cluster.lo        = nan(num.levels, num.levels, num.conds, num_params);
    sd_ci_cluster.hi        = nan(num.levels, num.levels, num.conds, num_params);
    sd_ci_cluster.fwhm_lo   = nan(num.levels, num.levels, num.conds);
    sd_ci_cluster.fwhm_hi   = nan(num.levels, num.levels, num.conds);
    sd_ci_cluster.curve_lo  = nan(num.levels, num.levels, num.conds, n_curve);
    sd_ci_cluster.curve_hi  = nan(num.levels, num.levels, num.conds, n_curve);
    sd_ci_cluster.curve_x   = curve_x;
    sd_ci_cluster.B         = B;
    sd_ci_cluster.ci        = prc;
    sd_ci_cluster.ci_method = bootstrap.ci_method;
    sd_ci_cluster.method    = sd_boot_cluster.method;

    blank_ci = struct( ...
        'lo',       nan(num.levels, num.levels, num.conds, num_params), ...
        'hi',       nan(num.levels, num.levels, num.conds, num_params), ...
        'fwhm_lo',  nan(num.levels, num.levels, num.conds), ...
        'fwhm_hi',  nan(num.levels, num.levels, num.conds), ...
        'curve_lo', nan(num.levels, num.levels, num.conds, n_curve), ...
        'curve_hi', nan(num.levels, num.levels, num.conds, n_curve));
    sd_ci_cluster.percentile = blank_ci;
    sd_ci_cluster.bca = blank_ci;
    sd_ci_cluster.bca.z0           = nan(num.levels, num.levels, num.conds, num_params);
    sd_ci_cluster.bca.acceleration = nan(num.levels, num.levels, num.conds, num_params);

    if B < 1
        sd_boot_cluster.timing = struct('bootstrap_replicates_sec', NaN, 'jackknife_sec', NaN);
        sd_boot_cluster.compute_jackknife = bootstrap.compute_jackknife;
        sd_boot_cluster.discard_summary = buildDiscardSummary(sd_boot_cluster, B);
        return
    end

    % --- Subject sampling ---
    if isfield(bootstrap, 'subject_cluster_seed') && ~isempty(bootstrap.subject_cluster_seed)
        rng_state = rng;
        rng(bootstrap.subject_cluster_seed, 'twister');
        subject_samples = randi(num.subjs, B, num.subjs);
        rng(rng_state);
        sd_boot_cluster.seed = bootstrap.subject_cluster_seed;
    else
        subject_samples = randi(num.subjs, B, num.subjs);
        sd_boot_cluster.seed = [];
    end
    sd_boot_cluster.subject_samples = subject_samples;

    num_chunks = min(max(1, bootstrap.subject_cluster_num_chunks), B);
    chunk_size = ceil(B / num_chunks);
    actual_num_chunks = ceil(B / chunk_size);
    chunk_tasks = cell(actual_num_chunks, 1);
    for i_chunk = 1:actual_num_chunks
        chunk_start = (i_chunk - 1) * chunk_size + 1;
        chunk_end = min(i_chunk * chunk_size, B);
        chunk_tasks{i_chunk} = struct('replicate_indices', chunk_start:chunk_end);
    end

    inner_toggles = toggles;
    inner_toggles.parallelization = 0;
    inner_toggles.disp_on = 0;

    if isfield(toggles, 'disp_on') && toggles.disp_on
        disp(['Subject-cluster bootstrapping SD params: B = ' num2str(B) ...
            ', chunks = ' num2str(actual_num_chunks)]);
    end

    use_parallel = isfield(toggles, 'parallelization') && toggles.parallelization && actual_num_chunks > 1;
    chunk_results = cell(actual_num_chunks, 1);

    t_bs = tic;
    if use_parallel
        parfor i_chunk = 1:actual_num_chunks
            chunk_results{i_chunk} = processSubjectClusterBootstrapChunk( ...
                chunk_tasks{i_chunk}, subject_samples, delta_theta_windows, ...
                delta_theta_centers, curve_x, num, p, bootstrap, inner_toggles, trial_pool);
        end
    else
        for i_chunk = 1:actual_num_chunks
            chunk_results{i_chunk} = processSubjectClusterBootstrapChunk( ...
                chunk_tasks{i_chunk}, subject_samples, delta_theta_windows, ...
                delta_theta_centers, curve_x, num, p, bootstrap, inner_toggles, trial_pool);
        end
    end

    for i_chunk = 1:actual_num_chunks
        cr = chunk_results{i_chunk};
        reps = cr.replicate_indices;
        sd_boot_cluster.params(reps,:,:,:,:)         = cr.params;
        sd_boot_cluster.fwhm(reps,:,:,:)             = cr.fwhm;
        sd_boot_cluster.curves(reps,:,:,:,:)         = cr.curves;
        sd_boot_cluster.exit_flag(reps,:,:,:)        = cr.exit_flag;
        sd_boot_cluster.at_bound(reps,:,:,:,:)       = cr.at_bound;
        sd_boot_cluster.admitted(reps,:,:,:)         = cr.admitted;
        sd_boot_cluster.discard_reason(reps,:,:,:)   = cr.discard_reason;
    end
    sd_boot_cluster.timing.bootstrap_replicates_sec = toc(t_bs);

    % --- Jackknife (leave-one-subject-out) for BCa acceleration ---
    if bootstrap.compute_jackknife
        if isfield(toggles, 'disp_on') && toggles.disp_on
            disp(['  Computing leave-one-subject-out jackknife (N = ' num2str(num.subjs) ' fits)...']);
        end
        jk_use_parallel = isfield(toggles, 'parallelization') && toggles.parallelization;
        jk = runSubjectJackknife(delta_theta_windows, delta_theta_centers, curve_x, num, p, inner_toggles, jk_use_parallel, trial_pool);
    else
        jk = struct();
        jk.params       = nan(num.subjs, num.levels, num.levels, num.conds, num_params);
        jk.fwhm         = nan(num.subjs, num.levels, num.levels, num.conds);
        jk.curves       = nan(num.subjs, num.levels, num.levels, num.conds, n_curve);
        jk.duration_sec = 0;
        if strcmpi(bootstrap.ci_method, 'bca')
            warning('bootstrap:noJackknife', ...
                ['ci_method = ''bca'' but compute_jackknife = false. Falling back to percentile CI ' ...
                '(BCa requires leave-one-out jackknife).']);
            bootstrap.ci_method = 'percentile';
            sd_ci_cluster.ci_method = 'percentile';
            sd_boot_cluster.ci_method = 'percentile';
        end
    end
    sd_boot_cluster.jackknife = jk;
    sd_boot_cluster.timing.jackknife_sec = jk.duration_sec;

    % --- CI computation: percentile + BCa, then route active CI ---
    has_point_est = ~isempty(sd_point) && isfield(sd_point, 'all') && isfield(sd_point.all, 'params_est');

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                admitted_mask = squeeze(sd_boot_cluster.admitted(:, prev_lvl, curr_lvl, cond));

                % Per-parameter CIs
                for param = 1:num_params
                    v_all = squeeze(sd_boot_cluster.params(:, prev_lvl, curr_lvl, cond, param));
                    v = v_all(admitted_mask);

                    q = percentileFinite(v, prc);
                    sd_ci_cluster.percentile.lo(prev_lvl, curr_lvl, cond, param) = q(1);
                    sd_ci_cluster.percentile.hi(prev_lvl, curr_lvl, cond, param) = q(2);

                    if has_point_est
                        theta_hat = sd_point.all.params_est(prev_lvl, curr_lvl, cond, param);
                    else
                        theta_hat = NaN;
                    end
                    jk_vals = squeeze(jk.params(:, prev_lvl, curr_lvl, cond, param));
                    [a1, a2, z0, a_acc] = bcaQuantiles(v, jk_vals, theta_hat, alpha);
                    qb = percentileFinite(v, 100 * [a1, a2]);
                    sd_ci_cluster.bca.lo(prev_lvl, curr_lvl, cond, param) = qb(1);
                    sd_ci_cluster.bca.hi(prev_lvl, curr_lvl, cond, param) = qb(2);
                    sd_ci_cluster.bca.z0(prev_lvl, curr_lvl, cond, param) = z0;
                    sd_ci_cluster.bca.acceleration(prev_lvl, curr_lvl, cond, param) = a_acc;
                end

                % FWHM CI
                fwhm_all = squeeze(sd_boot_cluster.fwhm(:, prev_lvl, curr_lvl, cond));
                fwhm_v = fwhm_all(admitted_mask);
                qf = percentileFinite(fwhm_v, prc);
                sd_ci_cluster.percentile.fwhm_lo(prev_lvl, curr_lvl, cond) = qf(1);
                sd_ci_cluster.percentile.fwhm_hi(prev_lvl, curr_lvl, cond) = qf(2);

                if has_point_est
                    w_pt = max(sd_point.all.params_est(prev_lvl, curr_lvl, cond, 2), eps);
                    fwhm_pt = (2 * sqrt(log(2))) / w_pt;
                else
                    fwhm_pt = NaN;
                end
                jk_fwhm = squeeze(jk.fwhm(:, prev_lvl, curr_lvl, cond));
                [a1, a2] = bcaQuantiles(fwhm_v, jk_fwhm, fwhm_pt, alpha);
                qfb = percentileFinite(fwhm_v, 100 * [a1, a2]);
                sd_ci_cluster.bca.fwhm_lo(prev_lvl, curr_lvl, cond) = qfb(1);
                sd_ci_cluster.bca.fwhm_hi(prev_lvl, curr_lvl, cond) = qfb(2);

                % Smooth-curve pointwise CIs
                if has_point_est
                    dog_params_pt = squeeze(sd_point.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
                    if all(isfinite(dog_params_pt))
                        curve_pt_full = calcDoG(curve_x, dog_params_pt);
                    else
                        curve_pt_full = nan(1, n_curve);
                    end
                else
                    curve_pt_full = nan(1, n_curve);
                end

                for ix = 1:n_curve
                    c_all = squeeze(sd_boot_cluster.curves(:, prev_lvl, curr_lvl, cond, ix));
                    c_v = c_all(admitted_mask);

                    qc = percentileFinite(c_v, prc);
                    sd_ci_cluster.percentile.curve_lo(prev_lvl, curr_lvl, cond, ix) = qc(1);
                    sd_ci_cluster.percentile.curve_hi(prev_lvl, curr_lvl, cond, ix) = qc(2);

                    jk_c = squeeze(jk.curves(:, prev_lvl, curr_lvl, cond, ix));
                    [a1, a2] = bcaQuantiles(c_v, jk_c, curve_pt_full(ix), alpha);
                    qcb = percentileFinite(c_v, 100 * [a1, a2]);
                    sd_ci_cluster.bca.curve_lo(prev_lvl, curr_lvl, cond, ix) = qcb(1);
                    sd_ci_cluster.bca.curve_hi(prev_lvl, curr_lvl, cond, ix) = qcb(2);
                end
            end
        end
    end

    % --- Route active CI based on bootstrap.ci_method ---
    if strcmpi(bootstrap.ci_method, 'bca') && has_point_est
        active = sd_ci_cluster.bca;
    else
        % BCa requires a point estimate; fall back to percentile if not provided
        active = sd_ci_cluster.percentile;
        sd_ci_cluster.ci_method = 'percentile';
        sd_boot_cluster.ci_method = 'percentile';
    end
    sd_ci_cluster.lo       = active.lo;
    sd_ci_cluster.hi       = active.hi;
    sd_ci_cluster.fwhm_lo  = active.fwhm_lo;
    sd_ci_cluster.fwhm_hi  = active.fwhm_hi;
    sd_ci_cluster.curve_lo = active.curve_lo;
    sd_ci_cluster.curve_hi = active.curve_hi;

    % --- Discard summary ---
    sd_boot_cluster.discard_summary = buildDiscardSummary(sd_boot_cluster, B);
end

% =========================================================================
% Per-chunk bootstrap fitting
% =========================================================================

function chunk_result = processSubjectClusterBootstrapChunk(chunk_task, subject_samples, delta_theta_windows, delta_theta_centers, curve_x, num, p, bootstrap, toggles, trial_pool)
    reps = chunk_task.replicate_indices;
    if nargin < 10
        trial_pool = [];
    end
    n_reps = numel(reps);
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    n_curve = numel(curve_x);

    chunk_result = struct();
    chunk_result.replicate_indices = reps;
    chunk_result.params         = nan(n_reps, num.levels, num.levels, num.conds, num_params);
    chunk_result.fwhm           = nan(n_reps, num.levels, num.levels, num.conds);
    chunk_result.curves         = nan(n_reps, num.levels, num.levels, num.conds, n_curve);
    chunk_result.exit_flag      = nan(n_reps, num.levels, num.levels, num.conds);
    chunk_result.at_bound       = false(n_reps, num.levels, num.levels, num.conds, num_params);
    chunk_result.admitted       = false(n_reps, num.levels, num.levels, num.conds);
    chunk_result.discard_reason = strings(n_reps, num.levels, num.levels, num.conds);

    bounds_for_obj = p.sd_bounds(:, 1:num_params); % row 1 = upper, row 2 = lower
    ub = bounds_for_obj(1, :);
    lb = bounds_for_obj(2, :);
    span = max(abs(ub - lb), eps);

    for i_rep = 1:n_reps
        b = reps(i_rep);
        sampled_subjects = subject_samples(b, :);
        pooled_windows = poolSubjectDeltaThetaWindows(delta_theta_windows, sampled_subjects, num);
        rb_params = fitBootstrapResponseBias(pooled_windows, num, p);
        if ~isempty(trial_pool)
            pooled_trials = poolSubjectCondTrialPool(trial_pool, sampled_subjects, num);
        else
            pooled_trials = [];
        end
        [sd_params, sd_exit, sd_admit, sd_reason] = ...
            fitBootstrapSerialDependence(pooled_windows, rb_params, delta_theta_centers, num, p, toggles, pooled_trials);

        chunk_result.params(i_rep,:,:,:,:)       = sd_params;
        chunk_result.exit_flag(i_rep,:,:,:)      = sd_exit;
        chunk_result.discard_reason(i_rep,:,:,:) = sd_reason;

        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for cond = 1:num.conds
                    dog_params = squeeze(sd_params(prev_lvl, curr_lvl, cond, 1:3));
                    if all(isfinite(dog_params))
                        w = max(dog_params(2), eps);
                        chunk_result.fwhm(i_rep, prev_lvl, curr_lvl, cond) = (2 * sqrt(log(2))) / w;
                        chunk_result.curves(i_rep, prev_lvl, curr_lvl, cond, :) = calcDoG(curve_x, dog_params);
                    end

                    full_params = squeeze(sd_params(prev_lvl, curr_lvl, cond, 1:num_params))';
                    if any(isnan(full_params))
                        is_at_bound = false(1, num_params);
                    else
                        is_at_bound = ...
                            (abs(full_params - lb) ./ span < bootstrap.bound_tol) | ...
                            (abs(full_params - ub) ./ span < bootstrap.bound_tol);
                    end
                    chunk_result.at_bound(i_rep, prev_lvl, curr_lvl, cond, :) = is_at_bound;

                    is_admitted = sd_admit(prev_lvl, curr_lvl, cond);
                    if is_admitted && bootstrap.discard_failed_fits && ...
                            ~isnan(sd_exit(prev_lvl, curr_lvl, cond)) && ...
                            sd_exit(prev_lvl, curr_lvl, cond) <= 0
                        is_admitted = false;
                        chunk_result.discard_reason(i_rep, prev_lvl, curr_lvl, cond) = "failed_fit";
                    end
                    if is_admitted && bootstrap.discard_at_bound && any(is_at_bound)
                        is_admitted = false;
                        chunk_result.discard_reason(i_rep, prev_lvl, curr_lvl, cond) = "at_bound";
                    end
                    chunk_result.admitted(i_rep, prev_lvl, curr_lvl, cond) = is_admitted;
                end
            end
        end
    end
end

% =========================================================================
% Subject pooling for sampled (or jackknife) subject sets
% =========================================================================

function pooled_trials = poolSubjectCondTrialPool(trial_pool, sampled_subjects, num)
% Concatenate per-subject un-windowed trial vectors across the resampled subject set.
% trial_pool.ind(subj).{delta_thetas, probe_offsets, responses}{prev, curr, cond}.

    fields = {'delta_thetas', 'probe_offsets', 'responses'};
    for i_field = 1:numel(fields)
        pooled_trials.(fields{i_field}) = cell(num.levels, num.levels, num.conds);
    end

    for i_sample = 1:numel(sampled_subjects)
        subj = sampled_subjects(i_sample);
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for cond = 1:num.conds
                    for i_field = 1:numel(fields)
                        f = fields{i_field};
                        v = trial_pool.ind(subj).(f){prev_lvl, curr_lvl, cond};
                        if ~isempty(v)
                            pooled_trials.(f){prev_lvl, curr_lvl, cond} = ...
                                [pooled_trials.(f){prev_lvl, curr_lvl, cond}; v(:)];
                        end
                    end
                end
            end
        end
    end
end

function pooled_windows = poolSubjectDeltaThetaWindows(delta_theta_windows, sampled_subjects, num)
    fields = {'delta_thetas', 'probe_offsets', 'responses'};
    for i_field = 1:numel(fields)
        pooled_windows.(fields{i_field}) = cell(num.levels, num.levels, num.conds, num.delta_theta_windows);
    end

    for i_sample = 1:numel(sampled_subjects)
        subj = sampled_subjects(i_sample);
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for cond = 1:num.conds
                    for iw = 1:num.delta_theta_windows
                        for i_field = 1:numel(fields)
                            field = fields{i_field};
                            curr_vals = delta_theta_windows.ind(subj).(field){prev_lvl, curr_lvl, cond, iw};
                            if ~isempty(curr_vals)
                                pooled_windows.(field){prev_lvl, curr_lvl, cond, iw} = [ ...
                                    pooled_windows.(field){prev_lvl, curr_lvl, cond, iw}; curr_vals(:)];
                            end
                        end
                    end
                end
            end
        end
    end
end

% =========================================================================
% Probit per-window fit on the pooled super-subject
% =========================================================================

function rb_params = fitBootstrapResponseBias(pooled_windows, num, p)
    rb_params = nan(num.levels, num.levels, num.conds, num.delta_theta_windows, 2);

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                for iw = 1:num.delta_theta_windows
                    probe_offsets = pooled_windows.probe_offsets{prev_lvl, curr_lvl, cond, iw};
                    responses = pooled_windows.responses{prev_lvl, curr_lvl, cond, iw};
                    if isempty(probe_offsets) || isempty(responses)
                        continue
                    end

                    probe_offsets = probe_offsets(:);
                    responses = responses(:);
                    n = min(numel(probe_offsets), numel(responses));
                    if n == 0
                        continue
                    end

                    try
                        result = processResponseBiasTask({probe_offsets(1:n), responses(1:n)}, p);
                        rb_params(prev_lvl, curr_lvl, cond, iw, :) = result.params_est;
                    catch %#ok<CTCH>
                        % Leave failed bootstrap fits as NaN.
                    end
                end
            end
        end
    end
end

% =========================================================================
% DoG fit on pooled windows; returns admission flags + discard reasons
% =========================================================================

function [sd_params, sd_exit, sd_admit, sd_reason] = fitBootstrapSerialDependence(pooled_windows, rb_params, delta_theta_centers, num, p, toggles, pooled_trials)
    if nargin < 7
        pooled_trials = [];
    end
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    sd_params = nan(num.levels, num.levels, num.conds, num_params);
    sd_exit   = nan(num.levels, num.levels, num.conds);
    sd_admit  = false(num.levels, num.levels, num.conds);
    sd_reason = strings(num.levels, num.levels, num.conds);

    if ~isfield(p, 'sd_min_windows') || isempty(p.sd_min_windows)
        p.sd_min_windows = 3;
    end
    if ~isfield(p, 'sd_lobe_edges') || isempty(p.sd_lobe_edges)
        lobe_edges = [-90 -15 15 90];
    else
        lobe_edges = p.sd_lobe_edges;
    end

    dt_centers = delta_theta_centers(:)';

    for cond = 1:num.conds
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                if use_sse
                    mu_per_window = squeeze(rb_params(prev_lvl, curr_lvl, cond, :, 1));
                    has_data = false(num.delta_theta_windows, 1);
                    for iw = 1:num.delta_theta_windows
                        has_data(iw) = ~isempty(pooled_windows.delta_thetas{prev_lvl, curr_lvl, cond, iw});
                    end
                    valid_windows = has_data & isfinite(mu_per_window);

                    if sum(valid_windows) < p.sd_min_windows
                        sd_reason(prev_lvl, curr_lvl, cond) = "min_windows";
                        continue
                    end

                    valid_centers = dt_centers(valid_windows);
                    n_lobes = numel(lobe_edges) - 1;
                    has_lobe = false(1, n_lobes);
                    for k = 1:n_lobes
                        has_lobe(k) = any(valid_centers >= lobe_edges(k) & valid_centers <= lobe_edges(k+1));
                    end
                    if ~all(has_lobe)
                        sd_reason(prev_lvl, curr_lvl, cond) = "missing_lobe";
                        continue
                    end

                    task_data.probe_offsets = zeros(sum(valid_windows), 1);
                    task_data.responses     = mu_per_window(valid_windows);
                    task_data.delta_thetas  = dt_centers(valid_windows)';
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];
                else
                    % NLL mode: trial-level binary responses, each trial counted ONCE.
                    % Use pooled_trials (un-windowed); previously this branch vertcatted
                    % across overlapping Δθ windows and duplicated each trial ~32x.
                    if isempty(pooled_trials)
                        sd_reason(prev_lvl, curr_lvl, cond) = "no_trial_pool";
                        continue
                    end

                    probe_offsets = pooled_trials.probe_offsets{prev_lvl, curr_lvl, cond};
                    responses     = pooled_trials.responses{prev_lvl, curr_lvl, cond};
                    delta_thetas  = pooled_trials.delta_thetas{prev_lvl, curr_lvl, cond};
                    if isempty(probe_offsets) || isempty(responses) || isempty(delta_thetas)
                        sd_reason(prev_lvl, curr_lvl, cond) = "min_windows";
                        continue
                    end

                    task_data.probe_offsets = probe_offsets;
                    task_data.responses     = responses;
                    task_data.delta_thetas  = delta_thetas;
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];
                end

                try
                    result = processSerialDependenceTask(task_data, p);
                    sd_params(prev_lvl, curr_lvl, cond, 1:num_params) = result.params_est(1:num_params);
                    if isfield(result, 'exit_flag') && ~isempty(result.exit_flag)
                        sd_exit(prev_lvl, curr_lvl, cond) = result.exit_flag;
                    end
                    sd_admit(prev_lvl, curr_lvl, cond) = true;
                catch fitErr
                    sd_reason(prev_lvl, curr_lvl, cond) = "failed_fit";
                    if isfield(toggles, 'disp_on') && toggles.disp_on
                        warning('%s', ['Cluster bootstrap SD fit failed: ' fitErr.message]);
                    end
                end
            end
        end
    end
end

% =========================================================================
% Leave-one-subject-out jackknife (drives BCa acceleration)
% =========================================================================

function jk = runSubjectJackknife(delta_theta_windows, delta_theta_centers, curve_x, num, p, toggles, use_parallel, trial_pool)
    if nargin < 7 || isempty(use_parallel)
        use_parallel = false;
    end
    if nargin < 8
        trial_pool = [];
    end

    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    n_curve = numel(curve_x);
    N = num.subjs;

    sd_p_cell  = cell(N, 1);
    fwhm_cell  = cell(N, 1);
    curve_cell = cell(N, 1);

    pool = gcp('nocreate');
    run_parfor = use_parallel && ~isempty(pool) && N > 1;

    jk_elapsed = tic;
    if run_parfor
        parfor k = 1:N
            [sd_p, fk, ck] = jackknifeOneSubject(k, delta_theta_windows, delta_theta_centers, curve_x, num, p, toggles, trial_pool);
            sd_p_cell{k}  = sd_p;
            fwhm_cell{k}  = fk;
            curve_cell{k} = ck;
        end
    else
        for k = 1:N
            [sd_p, fk, ck] = jackknifeOneSubject(k, delta_theta_windows, delta_theta_centers, curve_x, num, p, toggles, trial_pool);
            sd_p_cell{k}  = sd_p;
            fwhm_cell{k}  = fk;
            curve_cell{k} = ck;
        end
    end

    jk = struct();
    jk.params = nan(N, num.levels, num.levels, num.conds, num_params);
    jk.fwhm   = nan(N, num.levels, num.levels, num.conds);
    jk.curves = nan(N, num.levels, num.levels, num.conds, n_curve);
    for k = 1:N
        if ~isempty(sd_p_cell{k})
            jk.params(k,:,:,:,:) = sd_p_cell{k};
        end
        if ~isempty(fwhm_cell{k})
            jk.fwhm(k,:,:,:) = fwhm_cell{k};
        end
        if ~isempty(curve_cell{k})
            jk.curves(k,:,:,:,:) = curve_cell{k};
        end
    end
    jk.duration_sec = toc(jk_elapsed);
end

function [sd_params, fwhm_k, curve_k] = jackknifeOneSubject(k, delta_theta_windows, delta_theta_centers, curve_x, num, p, toggles, trial_pool)
    if nargin < 8
        trial_pool = [];
    end
    n_curve = numel(curve_x);
    all_subjects = 1:num.subjs;
    sampled_subjects = all_subjects(all_subjects ~= k);
    pooled_windows = poolSubjectDeltaThetaWindows(delta_theta_windows, sampled_subjects, num);
    rb_params = fitBootstrapResponseBias(pooled_windows, num, p);
    if ~isempty(trial_pool)
        pooled_trials = poolSubjectCondTrialPool(trial_pool, sampled_subjects, num);
    else
        pooled_trials = [];
    end
    sd_params = fitBootstrapSerialDependence(pooled_windows, rb_params, delta_theta_centers, num, p, toggles, pooled_trials);

    fwhm_k  = nan(num.levels, num.levels, num.conds);
    curve_k = nan(num.levels, num.levels, num.conds, n_curve);
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                dog_params = squeeze(sd_params(prev_lvl, curr_lvl, cond, 1:3));
                if all(isfinite(dog_params))
                    w = max(dog_params(2), eps);
                    fwhm_k(prev_lvl, curr_lvl, cond)     = (2 * sqrt(log(2))) / w;
                    curve_k(prev_lvl, curr_lvl, cond, :) = calcDoG(curve_x, dog_params);
                end
            end
        end
    end
end

% =========================================================================
% BCa lower/upper percentiles on the unit interval
% =========================================================================

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

% =========================================================================
% Bootstrap admission summary
% =========================================================================

function summary = buildDiscardSummary(sd_boot_cluster, B)
    admitted_grid = sd_boot_cluster.admitted;
    discarded_grid = ~admitted_grid;

    summary = struct();
    summary.B_total             = B;
    summary.admitted_per_cell   = squeeze(sum(admitted_grid, 1));
    summary.discarded_per_cell  = squeeze(sum(discarded_grid, 1));
    summary.global_admitted     = sum(admitted_grid(:));
    summary.global_discarded    = sum(discarded_grid(:));

    reasons = sd_boot_cluster.discard_reason;
    summary.by_reason = struct( ...
        'min_windows',  sum(reasons(:) == "min_windows"), ...
        'missing_lobe', sum(reasons(:) == "missing_lobe"), ...
        'failed_fit',   sum(reasons(:) == "failed_fit"), ...
        'at_bound',     sum(reasons(:) == "at_bound") );
end

% =========================================================================
% Robust percentile that ignores NaNs / empty inputs
% =========================================================================

function q = percentileFinite(values, prc)
    values = values(isfinite(values));
    if isempty(values) || any(~isfinite(prc(:)))
        q = [NaN, NaN];
    else
        q = prctile(values, prc);
    end
end
