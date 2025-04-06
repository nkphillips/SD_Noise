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
    p_correct = (1 - guess_rate) * normcdf(internal_orient_diff, mu, sigma) + 0.5 * guess_rate;
    
    % Simulate response based on discrimination probability
    responses = rand() < p_correct;

end






