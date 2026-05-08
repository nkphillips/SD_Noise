function [responses, correct] = simulate_responses(p, gt)
% SIMULATE_RESPONSES Generates simulated 2AFC responses for Experiment 2
%   p:  experiment parameters (contains trial_events, num_blocks, etc.)
%   gt: ground truth parameters for the simulated subject
%       gt.contrast_alpha, gt.contrast_beta (Weibull)
%       gt.filter_alpha, gt.filter_beta (Weibull on precision = 1/fw)
%       gt.lambda (lapse rate)
%       gt.dog_amp, gt.dog_w (Serial dependence parameters, e.g. matrices [num_features x num_levels])

    responses = nan(p.num_trials_per_block, p.num_blocks);
    correct = nan(p.num_trials_per_block, p.num_blocks);

    for n_block = 1:p.num_blocks
        
        test_ori = p.trial_events(:, 1, n_block);
        probe_ori = p.trial_events(:, 2, n_block);
        levels = p.trial_events(:, 3, n_block);
        
        probe_offset = probe_ori - test_ori;
        % Wrap offset
        probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
        probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;
        
        is_CW = probe_offset > 0;
        is_CCW = probe_offset < 0;
        is_zero = probe_offset == 0;

        curr_feature = p.feature_order(n_block);
        
        %% 1. Calculate base sensitivity (sigma) or probability correct for each trial
        
        p_correct_base = zeros(length(levels), 1);
        
        if curr_feature == 1 % Contrast
            c_vals = p.contrast(levels)';
            p_correct_base = 0.5 + (0.5 - gt.lambda) * (1 - exp(-(c_vals ./ gt.contrast_alpha).^gt.contrast_beta));
        else % Filter Width (Weibull on precision = 1/fw)
            fw_vals = p.orientation_bp_filter_width(levels)';
            precision_vals = 1 ./ fw_vals;
            p_correct_base = 0.5 + (0.5 - gt.lambda) * (1 - exp(-(precision_vals ./ gt.filter_alpha).^gt.filter_beta));
        end
        
        %% 2. Calibration: direct binomial sampling from p_correct_base
        % For calibration blocks there is no serial dependence, so we sample
        % correct/incorrect directly from the psychometric function. This
        % avoids the Jensen's inequality compression that arises from mapping
        % p_correct through variable probe offsets.
        
        if p.calibration
            correct_trial = rand(length(levels), 1) < p_correct_base;
            responses(is_CW, n_block)  = correct_trial(is_CW);
            responses(is_CCW, n_block) = ~correct_trial(is_CCW);
            responses(is_zero, n_block) = rand(sum(is_zero), 1) > 0.5;
            correct(:, n_block) = correct_trial;
            correct(is_zero, n_block) = rand(sum(is_zero), 1) > 0.5;
            continue;
        end
        
        %% 3. Main experiment: convert P_correct to internal noise (sigma)
        % We model the probability of a CW response as a cumulative normal
        % over the probe offset. 
        % P_correct = normcdf(mean_abs_offset / sigma). We can invert this to find a working sigma.
        mean_abs_offset = mean(abs(probe_offset(~is_zero)));
        if isnan(mean_abs_offset) || mean_abs_offset == 0
            mean_abs_offset = 5; % fallback
        end
        
        % Invert: P_correct = normcdf(x / sigma)  =>  sigma = x / norminv(P_correct)
        adj_p_correct = max(0.51, min(0.99, p_correct_base)); 
        sigma_trial = mean_abs_offset ./ norminv(adj_p_correct);
        
        %% 4. Calculate Serial Dependence (Bias)
        
        bias = zeros(length(levels), 1);
        
        delta_theta = [0; test_ori(1:end-1) - test_ori(2:end)];
        delta_theta(delta_theta > 90) = delta_theta(delta_theta > 90) - 180;
        delta_theta(delta_theta < -90) = delta_theta(delta_theta < -90) + 180;
        
        prev_levels = [1; levels(1:end-1)];
        
        for t_idx = 2:length(levels)
            p_lvl = prev_levels(t_idx);
            
            amp = gt.dog_amp(curr_feature, p_lvl);
            w_val = gt.dog_w(curr_feature, p_lvl);
            
            % DoG: y = A * c * w * x * exp(-(w*x)^2), c = sqrt(e) so peak = A
            c = sqrt(exp(1)); 
            dt = delta_theta(t_idx);
            bias(t_idx) = amp * c * w_val * dt * exp(-(w_val * dt)^2);
        end
        
        %% 5. Calculate P(CW) and simulate responses
        
        % Subtract bias so that a positive bias (attractive SD) shifts the
        % perceived test orientation toward the previous stimulus.
        internal_offset = probe_offset - bias;
        
        p_CW = (1 - gt.lambda) * normcdf(internal_offset ./ sigma_trial) + gt.lambda * 0.5;
        
        rand_vals = rand(length(p_CW), 1);
        responses(:, n_block) = rand_vals < p_CW;
        
        correct(:, n_block) = (responses(:, n_block) == 1 & is_CW) | (responses(:, n_block) == 0 & is_CCW);
        correct(is_zero, n_block) = rand(sum(is_zero), 1) > 0.5;
        
    end

end