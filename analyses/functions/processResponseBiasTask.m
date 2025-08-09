function result = processResponseBiasTask(task_data, p)
    % processResponseBiasTask - Process a single response bias estimation task
    %
    % This function is a wrapper that takes task data and processes it using
    % the optimized response bias estimation function.
    %
    % Inputs:
    %   task_data - cell array containing {probe_offsets, responses}
    %   p - parameter struct for the estimation
    %
    % Outputs:
    %   result - struct containing estimation results
    
    % Extract data from task (task_data is already a cell array with {probe_offsets, responses})
    probe_offsets = task_data{1};
    responses = task_data{2};
    
    % Create fixed_params
    fixed_params = [probe_offsets, responses];
    
    % Run estimation (silent - progress is tracked at higher level)
    [start_params, start_nll, params_est, nll, exit_flag] = estimateResponseBias_optimized(fixed_params, p);
    
    % Package results
    result.start_params = start_params;
    result.start_nll = start_nll;
    result.params_est = params_est;
    result.nll = nll;
    result.exit_flag = exit_flag;
end 