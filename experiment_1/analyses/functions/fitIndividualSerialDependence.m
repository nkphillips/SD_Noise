function [sd, ind_sd_meta] = fitIndividualSerialDependence(sd, rb, delta_theta_windows, delta_theta_centers, num, p, toggles, trial_pool)
% fitIndividualSerialDependence
%
% Populate sd.ind by fitting the DoG separately for each subject and
% previous/current level pair. In SSE mode this uses each subject's fitted
% response-bias mu values across delta-theta windows. In NLL mode this uses
% trial_pool.ind(subj) (un-windowed trials, each counted once).

    if nargin < 7
        toggles = struct('disp_on', 0, 'parallelization', 0);
    end
    if nargin < 8
        trial_pool = [];
    end
    if ~isfield(p, 'sd_min_windows') || isempty(p.sd_min_windows)
        p.sd_min_windows = 3;
    end

    if strcmp(p.sd_objective, 'nll') && isempty(trial_pool)
        error('fitIndividualSerialDependence:missingTrialPool', ...
            ['NLL mode requires the un-windowed trial_pool (8th arg). ' ...
             'Build it with buildPerSubjectCondTrialPool(all_runs, p, num, n_back).']);
    end

    start_time = tic;
    task_list = {};
    task_indices = [];
    skipped_too_few_windows = 0;

    for subj = 1:num.subjs
        for cond = 1:num.conds
            for prev_lvl = 1:num.levels
                for curr_lvl = 1:num.levels
                    if strcmp(p.sd_objective, 'sse')
                        mu_per_window = squeeze(rb.ind(subj).params_est(prev_lvl, curr_lvl, cond, :, 1));
                        has_data = false(num.delta_theta_windows, 1);
                        for iw = 1:num.delta_theta_windows
                            has_data(iw) = ~isempty(delta_theta_windows.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond, iw});
                        end

                        valid_windows = has_data & isfinite(mu_per_window);
                        if sum(valid_windows) < p.sd_min_windows
                            skipped_too_few_windows = skipped_too_few_windows + 1;
                            continue
                        end

                        x_dt = delta_theta_centers(valid_windows)';
                        y_mu = mu_per_window(valid_windows)';

                        task_data.probe_offsets = zeros(size(x_dt));
                        task_data.responses = y_mu;
                        task_data.delta_thetas = x_dt;
                        task_data.condition_info = [prev_lvl, curr_lvl, cond];
                    else
                        % NLL mode: trial-level binary responses, each trial counted ONCE.
                        % (Previously vertcatted across overlapping windows -> ~32x duplication.)
                        probe_offsets = trial_pool.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond};
                        responses     = trial_pool.ind(subj).responses{prev_lvl, curr_lvl, cond};
                        delta_thetas  = trial_pool.ind(subj).delta_thetas{prev_lvl, curr_lvl, cond};

                        if isempty(probe_offsets) || isempty(responses) || isempty(delta_thetas)
                            skipped_too_few_windows = skipped_too_few_windows + 1;
                            continue
                        end

                        task_data.probe_offsets = probe_offsets;
                        task_data.responses = responses;
                        task_data.delta_thetas = delta_thetas;
                        task_data.condition_info = [prev_lvl, curr_lvl, cond];
                    end

                    task_list{end+1} = task_data; %#ok<AGROW>
                    task_indices = [task_indices; subj, prev_lvl, curr_lvl, cond]; %#ok<AGROW>
                end
            end
        end
    end

    num_tasks = numel(task_list);
    if isfield(toggles, 'disp_on') && toggles.disp_on
        disp(['  - Created ' num2str(num_tasks) ' individual serial-dependence tasks']);
        disp(['  - Skipped ' num2str(skipped_too_few_windows) ' cells with too few valid windows/trials']);
    end

    if num_tasks == 0
        ind_sd_meta = struct('duration', toc(start_time), 'num_tasks', 0, ...
            'skipped_too_few_windows', skipped_too_few_windows);
        return
    end

    use_parallel = isfield(toggles, 'parallelization') && toggles.parallelization && num_tasks > 1;
    num_chunks = min(p.num_chunks, num_tasks);
    [results, ~] = processTasksInChunks(task_list, num_chunks, use_parallel, ...
        @processSerialDependenceTask, toggles, p);

    for it = 1:num_tasks
        subj = task_indices(it, 1);
        prev_lvl = task_indices(it, 2);
        curr_lvl = task_indices(it, 3);
        cond = task_indices(it, 4);

        r = results{it};
        if strcmp(p.sd_objective, 'sse')
            sd.ind(subj).start_params(prev_lvl, curr_lvl, cond, 1:3) = r.start_params;
            sd.ind(subj).params_est(prev_lvl, curr_lvl, cond, 1:3) = r.params_est;
        else
            sd.ind(subj).start_params(prev_lvl, curr_lvl, cond, :) = r.start_params;
            sd.ind(subj).params_est(prev_lvl, curr_lvl, cond, :) = r.params_est;
        end
        sd.ind(subj).start_nll(prev_lvl, curr_lvl, cond) = r.start_metric;
        sd.ind(subj).nll(prev_lvl, curr_lvl, cond) = r.final_metric;
        sd.ind(subj).exit_flag(prev_lvl, curr_lvl, cond) = r.exit_flag;

        y_observed = task_list{it}.responses(:);
        dt = task_list{it}.delta_thetas(:);
        if strcmp(p.sd_objective, 'sse')
            y_fitted = calcDoG(dt, r.params_est(1:3));
        else
            po = task_list{it}.probe_offsets(:);
            mu = calcDoG(dt, r.params_est(1:3));
            y_fitted = (1 - p.guess_rate) * normcdf(po, mu, r.params_est(4)) + 0.5 * p.guess_rate;
        end
        sd.ind(subj).r2(prev_lvl, curr_lvl, cond) = calcR2(y_observed, y_fitted);
    end

    ind_sd_meta = struct();
    ind_sd_meta.duration = toc(start_time);
    ind_sd_meta.num_tasks = num_tasks;
    ind_sd_meta.skipped_too_few_windows = skipped_too_few_windows;
end
