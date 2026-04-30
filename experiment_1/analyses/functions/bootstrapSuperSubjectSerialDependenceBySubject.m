function [sd_ci_cluster, sd_boot_cluster] = bootstrapSuperSubjectSerialDependenceBySubject(delta_theta_windows, delta_theta_centers, num, p, bootstrap, toggles)
% bootstrapSuperSubjectSerialDependenceBySubject
%
% Subject-cluster bootstrap CIs for pooled super-subject DoG fits. Each
% bootstrap replicate resamples subjects with replacement, rebuilds a pooled
% super-subject data set from the sampled subjects, and reruns the same
% response-bias-to-serial-dependence fitting path used for the point estimate.

    if nargin < 6
        toggles = struct('disp_on', 0, 'parallelization', 0);
    end
    if ~isfield(bootstrap, 'B_subject_cluster_sd') || isempty(bootstrap.B_subject_cluster_sd)
        bootstrap.B_subject_cluster_sd = 2000;
    end
    if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
        bootstrap.ci = [2.5, 97.5];
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
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    curve_x = bootstrap.subject_cluster_curve_x(:)';
    n_curve = numel(curve_x);

    sd_boot_cluster = struct();
    sd_boot_cluster.params = nan(B, num.levels, num.levels, num.conds, num_params);
    sd_boot_cluster.fwhm = nan(B, num.levels, num.levels, num.conds);
    sd_boot_cluster.curves = nan(B, num.levels, num.levels, num.conds, n_curve);
    sd_boot_cluster.curve_x = curve_x;
    sd_boot_cluster.B = B;
    sd_boot_cluster.ci = prc;
    sd_boot_cluster.method = 'subject_cluster_percentile_pooled_super_subject';

    sd_ci_cluster = struct();
    sd_ci_cluster.lo = nan(num.levels, num.levels, num.conds, num_params);
    sd_ci_cluster.hi = nan(num.levels, num.levels, num.conds, num_params);
    sd_ci_cluster.fwhm_lo = nan(num.levels, num.levels, num.conds);
    sd_ci_cluster.fwhm_hi = nan(num.levels, num.levels, num.conds);
    sd_ci_cluster.curve_lo = nan(num.levels, num.levels, num.conds, n_curve);
    sd_ci_cluster.curve_hi = nan(num.levels, num.levels, num.conds, n_curve);
    sd_ci_cluster.curve_x = curve_x;
    sd_ci_cluster.B = B;
    sd_ci_cluster.ci = prc;
    sd_ci_cluster.method = sd_boot_cluster.method;

    if B < 1
        return
    end

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

    if use_parallel
        parfor i_chunk = 1:actual_num_chunks
            chunk_results{i_chunk} = processSubjectClusterBootstrapChunk( ...
                chunk_tasks{i_chunk}, subject_samples, delta_theta_windows, ...
                delta_theta_centers, curve_x, num, p, inner_toggles);
        end
    else
        for i_chunk = 1:actual_num_chunks
            chunk_results{i_chunk} = processSubjectClusterBootstrapChunk( ...
                chunk_tasks{i_chunk}, subject_samples, delta_theta_windows, ...
                delta_theta_centers, curve_x, num, p, inner_toggles);
        end
    end

    for i_chunk = 1:actual_num_chunks
        reps = chunk_results{i_chunk}.replicate_indices;
        sd_boot_cluster.params(reps,:,:,:,:) = chunk_results{i_chunk}.params;
        sd_boot_cluster.fwhm(reps,:,:,:) = chunk_results{i_chunk}.fwhm;
        sd_boot_cluster.curves(reps,:,:,:,:) = chunk_results{i_chunk}.curves;
    end

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            for cond = 1:num.conds
                for param = 1:num_params
                    v = squeeze(sd_boot_cluster.params(:, prev_lvl, curr_lvl, cond, param));
                    q = percentileFinite(v, prc);
                    sd_ci_cluster.lo(prev_lvl, curr_lvl, cond, param) = q(1);
                    sd_ci_cluster.hi(prev_lvl, curr_lvl, cond, param) = q(2);
                end

                fwhm_vals = squeeze(sd_boot_cluster.fwhm(:, prev_lvl, curr_lvl, cond));
                q_fwhm = percentileFinite(fwhm_vals, prc);
                sd_ci_cluster.fwhm_lo(prev_lvl, curr_lvl, cond) = q_fwhm(1);
                sd_ci_cluster.fwhm_hi(prev_lvl, curr_lvl, cond) = q_fwhm(2);

                for ix = 1:n_curve
                    curve_vals = squeeze(sd_boot_cluster.curves(:, prev_lvl, curr_lvl, cond, ix));
                    q_curve = percentileFinite(curve_vals, prc);
                    sd_ci_cluster.curve_lo(prev_lvl, curr_lvl, cond, ix) = q_curve(1);
                    sd_ci_cluster.curve_hi(prev_lvl, curr_lvl, cond, ix) = q_curve(2);
                end
            end
        end
    end
end

function chunk_result = processSubjectClusterBootstrapChunk(chunk_task, subject_samples, delta_theta_windows, delta_theta_centers, curve_x, num, p, toggles)
    reps = chunk_task.replicate_indices;
    n_reps = numel(reps);
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    n_curve = numel(curve_x);

    chunk_result = struct();
    chunk_result.replicate_indices = reps;
    chunk_result.params = nan(n_reps, num.levels, num.levels, num.conds, num_params);
    chunk_result.fwhm = nan(n_reps, num.levels, num.levels, num.conds);
    chunk_result.curves = nan(n_reps, num.levels, num.levels, num.conds, n_curve);

    for i_rep = 1:n_reps
        b = reps(i_rep);
        sampled_subjects = subject_samples(b, :);
        pooled_windows = poolSubjectDeltaThetaWindows(delta_theta_windows, sampled_subjects, num);
        rb_params = fitBootstrapResponseBias(pooled_windows, num, p);
        sd_params = fitBootstrapSerialDependence(pooled_windows, rb_params, delta_theta_centers, num, p, toggles);

        chunk_result.params(i_rep,:,:,:,:) = sd_params(:,:,:,:);
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for cond = 1:num.conds
                    dog_params = squeeze(sd_params(prev_lvl, curr_lvl, cond, 1:3));
                    if all(isfinite(dog_params))
                        w = max(dog_params(2), eps);
                        chunk_result.fwhm(i_rep, prev_lvl, curr_lvl, cond) = (2 * sqrt(log(2))) / w;
                        chunk_result.curves(i_rep, prev_lvl, curr_lvl, cond, :) = calcDoG(curve_x, dog_params);
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

function sd_params = fitBootstrapSerialDependence(pooled_windows, rb_params, delta_theta_centers, num, p, toggles)
    use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
    num_params = 3 + (~use_sse);
    sd_params = nan(num.levels, num.levels, num.conds, num_params);

    if ~isfield(p, 'sd_min_windows') || isempty(p.sd_min_windows)
        p.sd_min_windows = 3;
    end

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
                        continue
                    end

                    task_data.probe_offsets = zeros(sum(valid_windows), 1);
                    task_data.responses = mu_per_window(valid_windows);
                    task_data.delta_thetas = delta_theta_centers(valid_windows)';
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];
                else
                    curr_probe_offsets = pooled_windows.probe_offsets(prev_lvl, curr_lvl, cond, :);
                    curr_responses = pooled_windows.responses(prev_lvl, curr_lvl, cond, :);
                    curr_delta_thetas = pooled_windows.delta_thetas(prev_lvl, curr_lvl, cond, :);

                    probe_offsets = vertcat(curr_probe_offsets{:});
                    responses = vertcat(curr_responses{:});
                    delta_thetas = vertcat(curr_delta_thetas{:});
                    if isempty(probe_offsets) || isempty(responses) || isempty(delta_thetas)
                        continue
                    end

                    task_data.probe_offsets = probe_offsets;
                    task_data.responses = responses;
                    task_data.delta_thetas = delta_thetas;
                    task_data.condition_info = [prev_lvl, curr_lvl, cond];
                end

                try
                    result = processSerialDependenceTask(task_data, p);
                    sd_params(prev_lvl, curr_lvl, cond, 1:num_params) = result.params_est(1:num_params);
                catch fitErr
                    if isfield(toggles, 'disp_on') && toggles.disp_on
                        warning('%s', ['Cluster bootstrap SD fit failed: ' fitErr.message]);
                    end
                end
            end
        end
    end
end

function q = percentileFinite(values, prc)
    values = values(isfinite(values));
    if isempty(values)
        q = [NaN, NaN];
    else
        q = prctile(values, prc);
    end
end
