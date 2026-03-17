function [responses, correct] = simulate_calibration_responses(p)
% simulate_calibration_responses Simulates responses for the calibration phase

    alpha_c = p.true_params.alpha_c;
    beta_c = p.true_params.beta_c;
    signal_fw = p.true_params.signal_fw;
    nmul_fw = p.true_params.nmul_fw;
    nadd_fw = p.true_params.nadd_fw;
    ptm_g = p.true_params.gamma_fw;
    guess_rate = p.true_params.guess_rate;
    lambda = p.true_params.lambda;

    responses = nan(p.num_trials_per_block, p.num_blocks);
    correct = nan(p.num_trials_per_block, p.num_blocks);

    for n_block = 1:p.num_blocks
        curr_feature = p.feature_order(n_block);
        
        test_orientations = p.trial_events(:, 1, n_block);
        probe_orientations = p.trial_events(:, 2, n_block);
        levels = p.trial_events(:, 3, n_block);
        
        % Calculate wrap-around offset
        probe_offset = probe_orientations - test_orientations;
        probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
        probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;
        
        % We need to incorporate the stimulus level into the p_CW calculation 
        % to properly simulate the psychometric function.
        % First calculate the probability of being correct for each trial
        P_correct = nan(p.num_trials_per_block, 1);
        
        if curr_feature == 1 % Contrast feature
            actual_levels = p.contrast(levels)';
            % For contrast blocks, filter width is easiest (index 1)
            % Weibull
            P_correct = guess_rate + (1 - guess_rate - lambda) * (1 - exp(-(actual_levels./alpha_c).^beta_c));
        elseif curr_feature == 2 % Filter Width feature
            actual_levels = p.orientation_bp_filter_width(levels)';
            % PTM (Perceptual Template Model)
            d_prime = signal_fw.^ptm_g ./ sqrt((1 + nmul_fw.^2) .* actual_levels.^(2*ptm_g) ...
                + nmul_fw.^2 .* signal_fw.^(2*ptm_g) + nadd_fw.^2);
            P_correct = (1 - lambda) * normcdf(d_prime / sqrt(2)) + lambda * 0.5;
        end
        
        % Now determine the probability of responding CW (1) based on what the correct answer is
        p_CW = nan(p.num_trials_per_block, 1);
        
        is_CW = probe_offset > 0;
        is_CCW = probe_offset < 0;
        is_zero = probe_offset == 0;
        
        % If the correct answer is CW, the probability of saying CW is P_correct
        p_CW(is_CW) = P_correct(is_CW);
        
        % If the correct answer is CCW, the probability of saying CW is (1 - P_correct)
        % (i.e. the probability of being incorrect)
        p_CW(is_CCW) = 1 - P_correct(is_CCW);
        
        % For 0-offset, it's ambiguous, so 50% chance
        p_CW(is_zero) = 0.5;
        
        % Simulate responses (1 = CW, 0 = CCW) using the calculated probabilities
        responses(:, n_block) = rand(p.num_trials_per_block, 1) < p_CW;
        
        % Determine correctness
        correct(:, n_block) = (responses(:, n_block) == 1 & is_CW) | (responses(:, n_block) == 0 & is_CCW);
        
        % For ambiguous trials, assign correctness by chance or consider them incorrect (following experiment 1 convention)
        correct(is_zero, n_block) = rand(sum(is_zero), 1) < 0.5;
    end
end