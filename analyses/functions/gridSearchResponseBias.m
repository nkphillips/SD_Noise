function [best_params, best_sse] = gridSearchResponseBias(init_params, fixed_params, which_search, param_bounds)

    %% Define search steps

    switch which_search
        case 'coarse'
            search_steps = 10;
            mu_vals = linspace(param_bounds(1,1), param_bounds(2,1), search_steps);
            sigma_vals = linspace(param_bounds(1,2), param_bounds(2,2), search_steps);
        case 'fine'
            
            if init_params(1)/2 < param_bounds(2,1)
                mu_lb = param_bounds(2,1);
            else
                mu_lb = init_params(1)/2;
            end

            if init_params(1)*1.5 > param_bounds(1,1)
                mu_ub = param_bounds(1,1);
            else
                mu_ub = init_params(1)*1.5;
            end

            if init_params(2)/2 < param_bounds(2,2)
                sigma_lb = param_bounds(2,2);
            else
                sigma_lb = init_params(2)/2;
            end

            if init_params(2)*1.5 > param_bounds(1,2)
                sigma_ub = param_bounds(1,2);
            else
                sigma_ub = init_params(2)*1.5;
            end

            search_steps = 100;

            mu_vals = linspace(mu_lb, mu_ub, search_steps);
            sigma_vals = linspace(sigma_lb, sigma_ub, search_steps);
    end

    %% Extract parameters

    probe_offsets = fixed_params{1}(:,1);
    measured_p_CW = fixed_params{1}(:,2);
    guess_rate = fixed_params{2};

    %% Populate sse grid

    sse_grid = nan(length(sigma_vals), length(mu_vals));

    for n = 1:length(mu_vals)
        for i = 1:length(sigma_vals)
            estimated_p_CW = calc_pCW(probe_offsets, mu_vals(n), sigma_vals(i), guess_rate);
            sse_grid(i,n) = calcSSE(measured_p_CW, estimated_p_CW);
        end
    end

    %% Find best parameters

    [best_sse_per_mu, best_sigma_per_mu_index] = min(sse_grid);
    [best_sse, best_mu_index] = min(best_sse_per_mu);
    
    best_params = [ mean(mu_vals(best_mu_index)), mean(sigma_vals(best_sigma_per_mu_index(best_mu_index)))];

end

