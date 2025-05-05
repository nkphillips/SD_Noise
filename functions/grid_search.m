function [best_params, best_sse] = grid_search(offsets, pCW, guess_rate)

    mu_vals = -10:0.5:10;
    sigma_vals = 0.5:0.5:10;
    sse_grid = nan(length(sigma_vals), length(mu_vals));

    for n = 1:length(mu_vals)
        for i = 1:length(sigma_vals)
            estimate_pCW = calc_pCW(offsets, mu_vals(n), sigma_vals(i), guess_rate);
            sse_grid(i,n) = sum((pCW - estimate_pCW).^2);
        end
    end

    [best_sse_per_mu, best_sigma_per_mu_index] = min(sse_grid);
    [best_sse, best_mu_index] = min(best_sse_per_mu);
    
    best_params = [ mean(mu_vals(best_mu_index)), mean(sigma_vals(best_sigma_per_mu_index(best_mu_index)))];

end

