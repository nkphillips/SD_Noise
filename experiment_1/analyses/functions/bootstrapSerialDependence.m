function sd_ci = bootstrapSerialDependence(delta_theta_windows, delta_theta_centers, num, p, rb, bootstrap, toggles)
% bootstrapSerialDependence
%
% Compute bootstrap confidence intervals for super-subject serial dependence
% parameters across all [prev, curr, cond] combinations.
%
% Supports both SSE (window-level mu) and NLL (trial-level) objectives:
%   - SSE: resamples windows with replacement among valid windows
%   - NLL: resamples trials with replacement
%
% Inputs:
%   delta_theta_windows: struct with .all fields containing trial-level cells
%   delta_theta_centers: vector of window centers (for SSE mode)
%   num: struct with levels/conds/windows
%   p: parameter struct (expects p.sd_objective)
%   rb: response bias struct (expects rb.all.params_est with mu per window)
%   bootstrap: struct with fields:
%       .B             - number of bootstrap iterations (e.g., 200)
%       .ci            - 1x2 vector percentiles (e.g., [2.5 97.5])
%       .min_trials    - optional, minimum trials for NLL mode
%       .min_windows   - optional, minimum windows for SSE mode
%   toggles: struct for printing/parallelization flags
%
% Output:
%   sd_ci.lo / sd_ci.hi: [prev, curr, cond, param] lower/upper CIs
%       param order: [A, w, b] for SSE; [A, w, b, sigma] for NLL

if nargin < 7
    toggles = struct('disp_on', 0);
end
if ~isfield(bootstrap, 'B') || isempty(bootstrap.B)
    bootstrap.B = 10;
end
if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
    bootstrap.ci = [2.5, 97.5];
end

use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
num_params = 3 + (~use_sse);

% Preallocate CI arrays
sz = [num.levels, num.levels, num.conds, num_params];
sd_ci.lo = nan(sz);
sd_ci.hi = nan(sz);

% Build task list
tasks = struct('prev', {}, 'curr', {}, 'cond', {}, 'mode', {}, ...
    'x_dt', {}, 'y_mu', {}, 'probe_offsets', {}, 'responses', {}, 'delta_thetas', {});

for cond = 1:num.conds
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            if use_sse
                % Use observed mu per window from RB fits
                mu_per_window = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1));
                has_data = false(num.delta_theta_windows, 1);
                for iw = 1:num.delta_theta_windows
                    has_data(iw) = ~isempty(delta_theta_windows.all.delta_thetas{prev_lvl, curr_lvl, cond, iw});
                end
                valid_windows = has_data & ~isnan(mu_per_window);
                if any(valid_windows)
                    x_dt = delta_theta_centers(valid_windows)';
                    y_mu = mu_per_window(valid_windows)';
                    tasks(end+1) = struct('prev', prev_lvl, 'curr', curr_lvl, 'cond', cond, ...
                        'mode', 'sse', 'x_dt', x_dt, 'y_mu', y_mu, ...
                        'probe_offsets', [], 'responses', [], 'delta_thetas', []); %#ok<AGROW>
                end
            else
                % NLL: aggregate trials
                curr_probe_offsets = delta_theta_windows.all.probe_offsets(prev_lvl, curr_lvl, cond, :);
                curr_responses = delta_theta_windows.all.responses(prev_lvl, curr_lvl, cond, :);
                curr_delta_thetas = delta_theta_windows.all.delta_thetas(prev_lvl, curr_lvl, cond, :);
                po = vertcat(curr_probe_offsets{:});
                y  = vertcat(curr_responses{:});
                dt = vertcat(curr_delta_thetas{:});
                if ~isempty(po) && ~isempty(y) && ~isempty(dt)
                    tasks(end+1) = struct('prev', prev_lvl, 'curr', curr_lvl, 'cond', cond, ...
                        'mode', 'nll', 'x_dt', [], 'y_mu', [], ...
                        'probe_offsets', po, 'responses', y, 'delta_thetas', dt); %#ok<AGROW>
                end
            end
        end
    end
end

if toggles.disp_on
    disp(['Bootstrapping SD params across conditions: ' num2str(numel(tasks)) ' tasks, B = ' num2str(bootstrap.B)]);
end

% Prepare tasks for chunked processing
num_tasks = numel(tasks);
results = cell(num_tasks, 1);

if isfield(toggles,'parallelization') && toggles.parallelization && num_tasks > 1
    task_list = cell(num_tasks,1);
    for it = 1:num_tasks
        task_list{it} = tasks(it);
    end
    num_chunks = min(p.num_chunks, num_tasks);
    [results, ~] = processTasksInChunks(task_list, num_chunks, true, @processBootstrapSerialDependenceTask, toggles, p, bootstrap, toggles);
else
    for it = 1:num_tasks
        results{it} = processBootstrapSerialDependenceTask(tasks(it), p, bootstrap, toggles);
    end
end

% Write results back into CI arrays
for it = 1:num_tasks
    prev_lvl = tasks(it).prev;
    curr_lvl = tasks(it).curr;
    cond     = tasks(it).cond;
    r = results{it};
    if ~isempty(r)
        sd_ci.lo(prev_lvl, curr_lvl, cond, 1:num_params) = r.lo;
        sd_ci.hi(prev_lvl, curr_lvl, cond, 1:num_params) = r.hi;
    end
end

end


