function [responses, correct] = simulate_responses(p)

    mu = p.psychometric_params(1); % ie threshold
    sigma = p.psychometric_params(2); % ie slope
    guess_rate = p.psychometric_params(3);

    responses = nan(p.num_trials_per_block, p.num_blocks);
    correct = nan(p.num_trials_per_block, p.num_blocks);

    for n_block = 1:p.num_blocks

        probe_offset = p.trial_events(:, 2, n_block) - p.trial_events(:, 1, n_block);

        % Get response bias from serial dependence 
        % delta_theta = p.trial_events(1:end-1,1,n_block) - p.trial_events(2:end,1,n_block);
        % bias = nan(length(delta_theta),1);
        % for n = 1:length(delta_theta)
        %     curr_level = p.trial_events(n+1,end,n_block);
        %     bias(n) = gaussian_prime(p.DoG_params(curr_level,:), delta_theta(n));
        % end
        % bias = [0; bias];

        % Get internal orientation difference (actual diff + bias + noise)
    %     bias_from_sd = bias + p.DoG_params(3) .* randn(length(orient_diff), 1);
        bias_from_sd = 0;
        internal_orient_diff = probe_offset + bias_from_sd;

        % Calculate probability of pCW
        p_CW = calc_pCW(internal_orient_diff, mu, sigma, guess_rate);
        
        % Simulate responses (1 = CW)
        responses(:,n_block) = rand(length(p_CW), 1) < p_CW;
        
        % Determine correct vs incorrect responses based on sign of offset
        is_CW = probe_offset > 0;
        is_CCW = probe_offset < 0;

        correct(:,n_block) = (responses(:,n_block) == 1 & is_CW) | (responses(:,n_block) == 0 & is_CCW);
    end

end

