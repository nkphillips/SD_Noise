function responses = simulate_responses(p)

    % Calculate orientation difference
    orient_diff = p.trial_events(:, 2) - p.trial_events(:, 1);

    % Calculate difference between previous and current trial orientation
    delta_theta = p.trial_events(1:end-1,1) - p.trial_events(2:end,1);

    % Get response bias from serial dependence
    bias = nan(length(delta_theta),1);
    for n = 1:length(delta_theta)
        curr_level = p.trial_events(n,end);
        bias(n) = gaussian_prime(p.DoG_params(curr_level,:), delta_theta(n));
    end
    bias = [0; bias];

    % Get internal orientation difference (actual diff + bias + noise)
    internal_orient_diff = orient_diff + bias + p.DoG_params(3) .* randn(length(orient_diff), 1);

    % Calculate probability of correct discrimination using cumulative normal
    mu = p.psychometric_params(1); % ie threshold
    sigma = p.psychometric_params(2); % ie slope
    guess_rate = p.psychometric_params(3);

    p_CW = calc_pCW(internal_orient_diff, mu, sigma, guess_rate);
    
    % Simulate responses
    responses = rand(length(p_CW), 1) < p_CW;


end



% Determine correct vs incorrect responses based on sign of offset
    is_CW = orient_diff > 0;
    is_CCW = orient_diff < 0;

    correct = (responses == 1 & is_CW) | (responses == 0 & is_CCW);

% Store all relevant trial-level data in a debug matrix
    delta_theta_padded = [delta_theta; NaN];
    p.debug_matrix = [ ...
        orient_diff, ...          % 1
        delta_theta_padded, ...   % 2
        bias, ...                 % 3
        internal_orient_diff, ...% 4
        p_CW, ...                 % 5
        responses, ...            % 6
        correct ...               % 7
    ];

end

% Probability CW

function p_CW = calc_pCW(internal_diff, psychometric_params)
    % Calculate the probability of CW response from internal orientation difference
    mu = psychometric_params(1);
    sigma = psychometric_params(2);
    lapse = psychometric_params(3);

    p_CW = (1 - lapse) * normcdf(internal_diff, mu, sigma) + 0.5 * lapse;
end

% Evaluate correctness

function correct = correct_response(responses, probe_offsets)
    % Determine correctness based on probe offset direction
    correct = (responses == 1 & probe_offsets > 0) | ...
              (responses == 0 & probe_offsets < 0);
end





%% Gridsearch

function [best_params, best_sse] = grid_search_mu_sigma(offsets, pCW)

    mu_vals = -10:0.5:10;
    sigma_vals = 0.5:0.5:10;
    best_sse = inf;

    for n = mu_vals
        for i = sigma_vals
            predicted_values = normcdf(offsets, n, i);
            sse = sum((pCW - predicted_values).^2);
            if sse < best_sse
                best_sse = sse;
                best_params = [n, i];
            end
        end
    end
end

