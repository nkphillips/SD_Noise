function [best_params, best_nll] = gridSearchResponseBias(init_params, fixed_params, which_search, p)

    %% Define search steps

    switch which_search
        case 'coarse'
            search_steps = 10;
            mu_vals = linspace(p.rb_bounds(1,1), p.rb_bounds(2,1), search_steps);
            sigma_vals = linspace(p.rb_bounds(1,2), p.rb_bounds(2,2), search_steps);
        
        case 'fine'
            if init_params(1)/2 < p.rb_bounds(2,1)
                mu_lb = p.rb_bounds(2,1);
            else
                mu_lb = init_params(1)/2;
            end

            if init_params(1)*1.5 > p.rb_bounds(1,1)
                mu_ub = p.rb_bounds(1,1);
            else
                mu_ub = init_params(1)*1.5;
            end

            if init_params(2)/2 < p.rb_bounds(2,2)
                sigma_lb = p.rb_bounds(2,2);
            else
                sigma_lb = init_params(2)/2;
            end

            if init_params(2)*1.5 > p.rb_bounds(1,2)
                sigma_ub = p.rb_bounds(1,2);
            else
                sigma_ub = init_params(2)*1.5;
            end

            search_steps = 50;
            mu_vals = linspace(mu_lb, mu_ub, search_steps);
            sigma_vals = linspace(sigma_lb, sigma_ub, search_steps);
            
    end

    %% Vectorized evaluation

    % Create parameter grid
    [mu_grid, sigma_grid] = meshgrid(mu_vals, sigma_vals);
    params_matrix = [mu_grid(:), sigma_grid(:)];
    
    % Extract data
    probe_offsets = fixed_params(:,1);
    responses = fixed_params(:,2);
    
    % Vectorized evaluation
    nll_vec = zeros(size(params_matrix, 1), 1);
    
    for i = 1:size(params_matrix, 1)
        mu = params_matrix(i, 1);
        sigma = params_matrix(i, 2);        
        p_CW = calc_pCW(probe_offsets, mu, sigma, p.guess_rate);
        nll_vec(i) = calcNLL(responses, p_CW);
    end
    
    % Reshape back to grid
    nll_grid = reshape(nll_vec, size(mu_grid));

    %% Find best parameters

    [best_nll, best_linear_index] = min(nll_grid(:));
    [best_sigma_index, best_mu_index] = ind2sub(size(nll_grid), best_linear_index);
    
    best_params = [mu_vals(best_mu_index), sigma_vals(best_sigma_index)];

end 