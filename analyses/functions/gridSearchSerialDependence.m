function [best_params, best_sse] = gridSearchSerialDependence(init_params, fixed_params, which_search, param_bounds)

    %% Define search steps

    switch which_search
        case 'coarse'
            search_steps = 10;
            amplitude_vals = linspace(param_bounds(1,1), param_bounds(2,1), search_steps);
            w_vals = linspace(param_bounds(1,2), param_bounds(2,2), search_steps);
        case 'fine'
            
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

            if init_params(2)/2 < param_bounds(2,2)
                w_lb = param_bounds(2,2);
            else
                w_lb = init_params(2)/2;
            end

            if init_params(2)*1.5 > param_bounds(1,2)
                w_ub = param_bounds(1,2);
            else
                w_ub = init_params(2)*1.5;
            end

            search_steps = 100;

            amplitude_vals = linspace(amplitude_lb, amplitude_ub, search_steps);
            w_vals = linspace(w_lb, w_ub, search_steps);
    end

    %% Extract parameters

    delta_theta = fixed_params(:,1);
    measured_bias = fixed_params(:,2);

    %% Populate sse grid

    sse_grid = nan(length(w_vals), length(amplitude_vals));

    for n = 1:length(amplitude_vals)
        for i = 1:length(w_vals)
            estimated_bias = gaussianPrime([amplitude_vals(n), w_vals(i)], delta_theta);
            sse_grid(i,n) = calcSSE(measured_bias, estimated_bias);
        end
    end

    %% Find best parameters
    
    [best_sse_per_amplitude, best_w_per_amplitude_index] = min(sse_grid);
    [best_sse, best_amplitude_index] = min(best_sse_per_amplitude);
    
    best_params = [mean(amplitude_vals(best_amplitude_index)), mean(w_vals(best_w_per_amplitude_index(best_amplitude_index)))];

end