function result = processSerialDependenceTask(task_data, p, previous_results)
    % processSerialDependenceTask - Process a single serial dependence estimation task
    %
    % This function is a wrapper that takes task data and processes it using
    % the optimized serial dependence estimation function with smart initialization.
    %
    % Inputs:
    %   task_data - struct containing {probe_offsets, responses, delta_thetas, condition_info}
    %   p - parameter struct for the estimation
    %   previous_results - previous fitting results for smart initialization (optional)
    %
    % Outputs:
    %   result - struct containing estimation results
    
    if nargin < 3
        previous_results = [];
    end
    
    % Extract data from task
    probe_offsets = task_data.probe_offsets;
    responses = task_data.responses;
    delta_thetas = task_data.delta_thetas;
    condition_info = task_data.condition_info;
    
    % Ensure column vectors and consistent lengths
    probe_offsets = probe_offsets(:);
    responses = responses(:);
    delta_thetas = delta_thetas(:);
    n = min([numel(probe_offsets), numel(responses), numel(delta_thetas)]);
    probe_offsets = probe_offsets(1:n);
    responses = responses(1:n);
    delta_thetas = delta_thetas(1:n);

    % Create raw_data matrix
    raw_data = [probe_offsets, responses, delta_thetas];
    
    % Get smart initial parameters if previous results available
    if ~isempty(previous_results)
        init_params = getSmartInitialParams(raw_data, p, previous_results);
    else
        init_params = p.sd_init_params;
    end
    
    % Run estimation with smart initialization
    [start_params, start_metric, params_est, final_metric, exit_flag] = estimateSerialDependence(raw_data, p, init_params);
    
    % Package results
    result.start_params = start_params;
    result.start_metric = start_metric;
    result.params_est = params_est;
    result.final_metric = final_metric;
    result.exit_flag = exit_flag;
    result.condition_info = condition_info;
    
end 