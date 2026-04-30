function [rb, ind_rb_meta] = fitIndividualResponseBias(rb, delta_theta_windows, num, p, toggles)
% fitIndividualResponseBias
%
% Populate rb.ind by fitting the response-bias CDF separately for each
% subject, previous/current level pair, condition, and delta-theta window.
% This reuses the same task wrapper as the super-subject path.

    if nargin < 5
        toggles = struct('disp_on', 0, 'parallelization', 0);
    end

    start_time = tic;
    task_list = {};
    task_indices = [];

    for subj = 1:num.subjs
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for cond = 1:num.conds
                    for iw = 1:num.delta_theta_windows
                        probe_offsets = delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, iw};
                        responses = delta_theta_windows.ind(subj).responses{prev_lvl, curr_lvl, cond, iw};

                        if ~isempty(probe_offsets) && ~isempty(responses)
                            task_list{end+1} = {probe_offsets, responses}; %#ok<AGROW>
                            task_indices = [task_indices; subj, prev_lvl, curr_lvl, cond, iw]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end

    num_tasks = numel(task_list);
    if isfield(toggles, 'disp_on') && toggles.disp_on
        disp(['  - Created ' num2str(num_tasks) ' individual response-bias tasks']);
    end

    if num_tasks == 0
        ind_rb_meta = struct('duration', toc(start_time), 'num_tasks', 0);
        return
    end

    use_parallel = isfield(toggles, 'parallelization') && toggles.parallelization && num_tasks > 1;
    num_chunks = min(p.num_chunks, num_tasks);
    [results, ~] = processTasksInChunks(task_list, num_chunks, use_parallel, ...
        @processResponseBiasTask, toggles, p);

    for it = 1:num_tasks
        subj = task_indices(it, 1);
        prev_lvl = task_indices(it, 2);
        curr_lvl = task_indices(it, 3);
        cond = task_indices(it, 4);
        iw = task_indices(it, 5);

        r = results{it};
        probe_offsets = delta_theta_windows.ind(subj).probe_offsets{prev_lvl, curr_lvl, cond, iw};

        rb.ind(subj).start_params(prev_lvl, curr_lvl, cond, iw, :) = r.start_params;
        rb.ind(subj).start_nll(prev_lvl, curr_lvl, cond, iw) = r.start_nll;
        rb.ind(subj).params_est(prev_lvl, curr_lvl, cond, iw, :) = r.params_est;
        rb.ind(subj).nll(prev_lvl, curr_lvl, cond, iw) = r.nll;
        rb.ind(subj).null_nll(prev_lvl, curr_lvl, cond, iw) = -numel(probe_offsets) * log(0.5);
        rb.ind(subj).better_than_null(prev_lvl, curr_lvl, cond, iw) = ...
            r.nll < rb.ind(subj).null_nll(prev_lvl, curr_lvl, cond, iw);
        rb.ind(subj).exit_flag(prev_lvl, curr_lvl, cond, iw) = r.exit_flag;
    end

    ind_rb_meta = struct();
    ind_rb_meta.duration = toc(start_time);
    ind_rb_meta.num_tasks = num_tasks;
end
