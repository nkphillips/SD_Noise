function [best_params, best_sse] = gridSearchSerialDependence(init_params, fixed_params, which_search, param_bounds)

    %% Define search steps

    switch which_search
        case 'coarse'
            search_steps = 20;
            amplitude_vals = linspace(param_bounds(2,1), param_bounds(1,1), search_steps); % param_bounds(2,1) is lower_amp, param_bounds(1,1) is upper_amp
            % param_bounds(2,2) is lower_width_param, param_bounds(1,2) is upper_width_param
            if param_bounds(2,2) <= 0 || param_bounds(1,2) <= 0 % logspace requires positive bounds
                % Fallback to linspace if bounds are not suitable for logspace (e.g., include zero or negative)
                width_param_vals = linspace(param_bounds(2,2), param_bounds(1,2), search_steps);
            else
                width_param_vals = logspace(log10(param_bounds(2,2)), log10(param_bounds(1,2)), search_steps);
            end
            
        case 'fine'
            search_steps = 100; % Keep fine search dense

            % Amplitude bounds for fine search
            if init_params(1)/2 < param_bounds(2,1)
                amplitude_lb = param_bounds(2,1);
            else
                amplitude_lb = init_params(1)/2;
            end

            if init_params(1)*1.5 > param_bounds(1,1)
                amplitude_ub = param_bounds(1,1);
            else
                amplitude_ub = init_params(1)*1.5;
            end
            amplitude_vals = linspace(amplitude_lb, amplitude_ub, search_steps);

            % Width_param bounds for fine search
            % Ensure init_params(2) (initial width_param) is positive for logspace calculations
            current_init_width = init_params(2);
            if current_init_width <= 0
                current_init_width = param_bounds(2,2); % Fallback to lower bound if init is not positive
            end
            
            width_param_lb_candidate = current_init_width / 2;
            width_param_ub_candidate = current_init_width * 1.5;

            if width_param_lb_candidate < param_bounds(2,2)
                width_param_lb = param_bounds(2,2);
            else
                width_param_lb = width_param_lb_candidate;
            end

            if width_param_ub_candidate > param_bounds(1,2)
                width_param_ub = param_bounds(1,2);
            else
                width_param_ub = width_param_ub_candidate;
            end
            
            % Ensure fine search bounds for width_param are positive for logspace
            if width_param_lb <= 0 || width_param_ub <= 0
                 % Fallback to linspace if bounds are not suitable for logspace
                width_param_vals = linspace(width_param_lb, width_param_ub, search_steps);
            else
                % Ensure lb is not greater than ub for logspace
                if width_param_lb > width_param_ub
                    width_param_lb = width_param_ub / 2; % Adjust if lb > ub to prevent error, or use global bound
                     if width_param_lb < param_bounds(2,2)
                         width_param_lb = param_bounds(2,2);
                     end
                     if width_param_lb > width_param_ub % if still problematic, default to narrower range
                         width_param_lb = width_param_ub;
                     end
                end
                if width_param_lb == width_param_ub % if they become equal, logspace will error
                    width_param_vals = repmat(width_param_lb, 1, search_steps); % or linspace
                else
                    width_param_vals = logspace(log10(width_param_lb), log10(width_param_ub), search_steps);
                end
            end
    end

    %% Extract parameters

    delta_theta = fixed_params(:,1);
    measured_bias = fixed_params(:,2);

    %% Populate sse grid

    sse_grid = nan(length(width_param_vals), length(amplitude_vals));

    for n = 1:length(amplitude_vals)
        for i = 1:length(width_param_vals)
            % Pass amplitude_vals(n) and width_param_vals(i) to gaussianPrime
            estimated_bias = gaussianPrime([amplitude_vals(n), width_param_vals(i)], delta_theta);
            sse_grid(i,n) = calcSSE(measured_bias, estimated_bias);
        end
    end

    %% Find best parameters
    
    % sse_grid has dimensions (length(width_param_vals), length(amplitude_vals))
    % Each column corresponds to an amplitude, each row to a width_param.
    
    [min_sse_each_column, row_indices_of_min_sse] = min(sse_grid, [], 1); % Finds min SSE in each column (best width_param for each amplitude)
                                                                          % row_indices_of_min_sse contains the index of the best width_param for each amplitude.
    
    [best_sse, col_index_of_best_sse] = min(min_sse_each_column); % Finds the overall best SSE and which amplitude column it occurred in.
    
    best_amplitude_val = amplitude_vals(col_index_of_best_sse);
    best_width_param_idx = row_indices_of_min_sse(col_index_of_best_sse); % Get the row index for width_param corresponding to the best amplitude
    best_width_param_val = width_param_vals(best_width_param_idx);
    
    best_params = [best_amplitude_val, best_width_param_val];

end