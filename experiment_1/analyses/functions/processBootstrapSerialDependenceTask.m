function result = processBootstrapSerialDependenceTask(task_data, p, bootstrap, toggles)
% processBootstrapSerialDependenceTask
% Bootstrap one condition's serial dependence parameters and return CIs.
%
% Inputs (task_data depends on objective):
%   - For SSE (window-level fit):
%       task_data.mode = 'sse'
%       task_data.x_dt: [n_windows x 1] delta-theta centers
%       task_data.y_mu: [n_windows x 1] observed mu per window (from RB fits)
%   - For NLL (trial-level fit):
%       task_data.mode = 'nll'
%       task_data.probe_offsets: [n_trials x 1]
%       task_data.responses:     [n_trials x 1]
%       task_data.delta_thetas:  [n_trials x 1]
%
% Outputs:
%   result: struct with fields .lo and .hi (1xP) for parameters
%           P = 3 if SSE objective; P = 4 if NLL objective
%
% Notes:
% - Uses p.sd_objective to determine parameter count and bounds.
% - bootstrap.B iterations; bootstrap.ci = [lo hi] percentiles.

if nargin < 4
    toggles = struct('disp_on', 0);
end
if ~isfield(bootstrap, 'B') || isempty(bootstrap.B)
    bootstrap.B = 10;
end
if ~isfield(bootstrap, 'ci') || isempty(bootstrap.ci)
    bootstrap.ci = [2.5, 97.5];
end
if ~isfield(bootstrap, 'min_trials') || isempty(bootstrap.min_trials)
    bootstrap.min_trials = 0;
end
if ~isfield(bootstrap, 'min_windows') || isempty(bootstrap.min_windows)
    bootstrap.min_windows = 3;
end

B = bootstrap.B;
prc = bootstrap.ci;

use_sse = isfield(p, 'sd_objective') && strcmp(p.sd_objective, 'sse');
num_params = 3 + (~use_sse);

params_b = nan(B, num_params);

switch task_data.mode
    case 'sse'
        x_dt = task_data.x_dt(:);
        y_mu = task_data.y_mu(:);
        n = min(numel(x_dt), numel(y_mu));
        x_dt = x_dt(1:n);
        y_mu = y_mu(1:n);
        if n < bootstrap.min_windows
            result.lo = nan(1, num_params);
            result.hi = nan(1, num_params);
            return
        end
        for b = 1:B
            try
                idx = randi(n, n, 1);
                x_b = x_dt(idx);
                y_b = y_mu(idx);
                raw_data = [zeros(size(x_b)), y_b, x_b]; % [probe_offsets, responses(=mu), delta_thetas]
                init_params = p.sd_init_params;
                [~, ~, params_est_b] = estimateSerialDependence(raw_data, p, init_params);
                params_b(b,1:3) = params_est_b(1:3);
            catch %#ok<CTCH>
                % leave as NaN
            end
        end

    case 'nll'
        po = task_data.probe_offsets(:);
        y  = task_data.responses(:);
        dt = task_data.delta_thetas(:);
        n = min([numel(po), numel(y), numel(dt)]);
        po = po(1:n); y = y(1:n); dt = dt(1:n);
        if n <= bootstrap.min_trials
            result.lo = nan(1, num_params);
            result.hi = nan(1, num_params);
            return
        end
        for b = 1:B
            try
                idx = randi(n, n, 1);
                po_b = po(idx);
                y_b  = y(idx);
                dt_b = dt(idx);
                raw_data = [po_b, y_b, dt_b];
                init_params = p.sd_init_params;
                [~, ~, params_est_b] = estimateSerialDependence(raw_data, p, init_params);
                if use_sse
                    params_b(b,1:3) = params_est_b(1:3);
                else
                    params_b(b,1:4) = params_est_b(1:4);
                end
            catch %#ok<CTCH>
                % leave as NaN
            end
        end
    otherwise
        error('Unknown task_data.mode: %s', task_data.mode);
end

% Compute percentiles per parameter
lo = nan(1, num_params);
hi = nan(1, num_params);
for k = 1:num_params
    v = params_b(:,k);
    v = v(~isnan(v));
    if isempty(v)
        lo(k) = NaN; hi(k) = NaN;
    else
        q = prctile(v, prc);
        lo(k) = q(1); hi(k) = q(2);
    end
end
result.lo = lo;
result.hi = hi;

end


