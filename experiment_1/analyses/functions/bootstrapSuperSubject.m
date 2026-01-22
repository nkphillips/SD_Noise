function [rb_ci, perf_ci] = bootstrapSuperSubject(delta_theta_windows, num, p, bootstrap, toggles)
% bootstrapSuperSubject
%
% Compute bootstrap confidence intervals for super-subject windowed metrics.
% - Response bias (mu per \Delta\theta window) via refitting the RB model
% - Performance metrics (Percent Correct, Percent CCW) per window
%
% Inputs:
%   delta_theta_windows: struct with .all fields containing trial-level cells
%   num: struct with levels/conds/windows
%   p: parameter struct (used for RB fitting)
%   bootstrap: struct with fields:
%       .B  - number of bootstrap iterations (e.g., 200)
%       .ci - 1x2 vector of percentiles (e.g., [2.5 97.5])
%   toggles: struct for printing/parallelization flags
%
% Outputs:
%   rb_ci.mu_lo/mu_hi: [prev, curr, cond, window] lower/upper CI for mu
%   perf_ci.(pc_lo|pc_hi|pccw_lo|pccw_hi): same shape, CIs for performance

    if nargin < 6
        toggles = struct('disp_on', 0);
    end

    if ~isfield(bootstrap, 'B') || isempty(bootstrap.B)
        bootstrap.B = 200;
    end
    if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
        bootstrap.ci = [2.5, 97.5];
    end

    B = bootstrap.B;
    prc = bootstrap.ci;

    % Preallocate CI outputs with NaNs
    sz = [num.levels, num.levels, num.conds, num.delta_theta_windows];
    rb_ci.mu_lo = nan(sz);
    rb_ci.mu_hi = nan(sz);
    perf_ci.pc_lo = nan(sz);
    perf_ci.pc_hi = nan(sz);
    perf_ci.pccw_lo = nan(sz);
    perf_ci.pccw_hi = nan(sz);

    % Build a task list of windows with data for efficiency
    tasks = struct('prev', {}, 'curr', {}, 'cond', {}, 'win', {}, 'po', {}, 'y', {});
    for cond = 1:num.conds
        for prev_lvl = 1:num.levels
            for curr_lvl = 1:num.levels
                for iw = 1:num.delta_theta_windows
                    po = delta_theta_windows.all.probe_offsets{prev_lvl, curr_lvl, cond, iw};
                    y  = delta_theta_windows.all.responses{prev_lvl, curr_lvl, cond, iw};
                    if ~isempty(po) && ~isempty(y)
                        % Ensure column vectors and equal length
                        po = po(:); y = y(:);
                        n = min(numel(po), numel(y));
                        po = po(1:n); y = y(1:n);
                        tasks(end+1) = struct('prev', prev_lvl, 'curr', curr_lvl, 'cond', cond, 'win', iw, 'po', po, 'y', y); %#ok<AGROW>
                    end
                end
            end
        end
    end

    if toggles.disp_on
        disp(['Bootstrapping super-subject (windows with data): ' num2str(numel(tasks)) ' windows, B = ' num2str(B)]);
    end

    % Prepare task list for parallel/sequential processing
    num_tasks = numel(tasks);
    results = cell(num_tasks, 1);

    if isfield(toggles,'parallelization') && toggles.parallelization && num_tasks > 1
        % Use the unified chunking infrastructure for consistency
        task_list = cell(num_tasks,1);
        for it = 1:num_tasks
            task_list{it} = tasks(it);
        end
        num_chunks = min(p.num_chunks, num_tasks);
        [results, ~] = processTasksInChunks(task_list, num_chunks, true, @processBootstrapWindowTask, toggles, p, bootstrap, toggles);
    else
        for it = 1:num_tasks
            results{it} = processBootstrapWindowTask(tasks(it), p, bootstrap, toggles);
        end
    end

    % Write results back into CI arrays
    for it = 1:num_tasks
        prev_lvl = tasks(it).prev;
        curr_lvl = tasks(it).curr;
        cond     = tasks(it).cond;
        iw       = tasks(it).win;
        r = results{it};
        if ~isempty(r)
            rb_ci.mu_lo(prev_lvl, curr_lvl, cond, iw) = r.mu_lo;
            rb_ci.mu_hi(prev_lvl, curr_lvl, cond, iw) = r.mu_hi;
            perf_ci.pc_lo(prev_lvl, curr_lvl, cond, iw) = r.pc_lo;
            perf_ci.pc_hi(prev_lvl, curr_lvl, cond, iw) = r.pc_hi;
            perf_ci.pccw_lo(prev_lvl, curr_lvl, cond, iw) = r.pccw_lo;
            perf_ci.pccw_hi(prev_lvl, curr_lvl, cond, iw) = r.pccw_hi;
        end
    end

    % Optional summary
    if isfield(toggles,'disp_on') && toggles.disp_on
        n_mu = sum(isfinite(rb_ci.mu_lo(:)) & isfinite(rb_ci.mu_hi(:)));
        n_pc = sum(isfinite(perf_ci.pc_lo(:)) & isfinite(perf_ci.pc_hi(:)));
        n_pccw = sum(isfinite(perf_ci.pccw_lo(:)) & isfinite(perf_ci.pccw_hi(:)));
        disp(['  - Bootstrap CIs computed for windows (mu/pc/pccw): ' num2str(n_mu) ' / ' num2str(n_pc) ' / ' num2str(n_pccw)]);
    end

end


