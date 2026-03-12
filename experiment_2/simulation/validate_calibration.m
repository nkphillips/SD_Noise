% validate_calibration.m
% Validates the calibration fits using the simulated data

clear all; close all; clc;

%% Set directories
script_dir = pwd;
data_dir = '../data';

num_subjs = 5;
subj_IDs = cell(num_subjs,1);
for subj = 1:num_subjs
    subj_IDs{subj} = num2str(9900 + subj - 1);
end

target_levels = [0.65, 0.75, 0.85];
gamma = 0.5;
lambda = 0.01;
options = optimset('Display', 'off');

%% Validate each subject
for subj = 1:num_subjs

    subj_ID = subj_IDs{subj};

    % Load all calibration runs for this subject
    file_pattern = [data_dir '/' subj_ID '/SD_Noise_Exp2_Calibration_S' subj_ID '_Run*.mat'];
    files = dir(file_pattern);

    if isempty(files)
        warning('No calibration data found for subject %s. Skipping.', subj_ID);
        continue;
    end

    all_contrast_correct = [];
    all_contrast_levels = [];
    all_filter_correct = [];
    all_filter_levels = [];

    true_params = [];

    for f = 1:length(files)
        load([data_dir '/' subj_ID '/' files(f).name], 'run_info');
        p = run_info.p;
        behav_data = run_info.behav_data;

        if isempty(true_params) && isfield(p, 'true_params')
            true_params = p.true_params;
        end

        for n_block = 1:p.num_blocks
            curr_feature = p.feature_order(n_block);

            % Calculate probe offset to identify 0-offset trials
            test_ori = p.trial_events(:, 1, n_block);
            probe_ori = p.trial_events(:, 2, n_block);
            probe_offset = probe_ori - test_ori;

            % Wrap offset
            probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
            probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;

            % We exclude 0-offset trials because there is no true "correct" answer.
            % An observer will guess 50% of the time regardless of contrast/noise,
            % which artificially pulls down the asymptote from 1.0 to ~0.92,
            % causing the Weibull/EqNoise fits (which assume lambda=0.01) to underestimate slope.
            valid_trials = probe_offset ~= 0;

            levels = p.trial_events(valid_trials, 3, n_block);
            correct = behav_data.correct(valid_trials, n_block);

            if curr_feature == 1 % Contrast feature
                actual_levels = p.contrast(levels)';
                all_contrast_levels = [all_contrast_levels; actual_levels];
                all_contrast_correct = [all_contrast_correct; correct];
            elseif curr_feature == 2 % Filter Width feature
                actual_levels = p.orientation_bp_filter_width(levels)';
                all_filter_levels = [all_filter_levels; actual_levels];
                all_filter_correct = [all_filter_correct; correct];
            end
        end
    end

    % Compute proportions
    [unique_c, ~, idx_c] = unique(all_contrast_levels);
    n_c = accumarray(idx_c, 1);
    k_c = accumarray(idx_c, all_contrast_correct, [], @sum);
    prop_c = k_c ./ n_c;

    [unique_fw, ~, idx_fw] = unique(all_filter_levels);
    n_fw = accumarray(idx_fw, 1);
    k_fw = accumarray(idx_fw, all_filter_correct, [], @sum);
    prop_fw = k_fw ./ n_fw;

    %% Fit Models
    % --- A. Contrast (Weibull) ---
    weibull_prob = @(x, alpha, beta) gamma + (1 - gamma - lambda) * (1 - exp(-(x./alpha).^beta));

    nll_weibull = @(params) -sum( ...
        k_c .* log(max(eps, weibull_prob(unique_c, params(1), params(2)))) + ...
        (n_c - k_c) .* log(max(eps, 1 - weibull_prob(unique_c, params(1), params(2)))) ...
        );

    % Coarse grid search for Contrast
    alpha_grid = linspace(0.01, 1.0, 20);
    beta_grid = linspace(0.5, 5.0, 20);
    min_nll_c = inf;
    best_guess_c = [mean(unique_c), 2];
    for a = alpha_grid
        for b = beta_grid
            curr_nll = nll_weibull([a, b]);
            if curr_nll < min_nll_c
                min_nll_c = curr_nll;
                best_guess_c = [a, b];
            end
        end
    end

    lb_c = [0.001, 0.5];
    ub_c = [1.0, 5.0];
    best_params_c = fmincon(nll_weibull, best_guess_c, [], [], [], [], lb_c, ub_c, [], options);
    alpha_c = best_params_c(1);
    beta_c = best_params_c(2);

    % --- B. Filter Width (Equivalent Noise) ---
    en_prob = @(x, sig, sig_int) (1 - lambda) * normcdf( (sig ./ sqrt(x.^2 + sig_int.^2)) / sqrt(2) ) + lambda * 0.5;

    nll_en = @(params) -sum( ...
        k_fw .* log(max(eps, en_prob(unique_fw, params(1), params(2)))) + ...
        (n_fw - k_fw) .* log(max(eps, 1 - en_prob(unique_fw, params(1), params(2)))) ...
        );

    % Coarse grid search for Filter Width
    sig_grid = linspace(1, 50, 20);
    sig_int_grid = linspace(1, 30, 20);
    min_nll_fw = inf;
    best_guess_fw = [15, 10];
    for s = sig_grid
        for si = sig_int_grid
            curr_nll = nll_en([s, si]);
            if curr_nll < min_nll_fw
                min_nll_fw = curr_nll;
                best_guess_fw = [s, si];
            end
        end
    end

    lb_fw = [1.0, 1.0];
    ub_fw = [50.0, 30.0];
    best_params_fw = fmincon(nll_en, best_guess_fw, [], [], [], [], lb_fw, ub_fw, [], options);
    signal_fw = best_params_fw(1);
    sigma_int_fw = best_params_fw(2);

    %% Calculate Errors
    fprintf('Subject %s Parameters Validation:\n', subj_ID);
    fprintf('  Contrast (Alpha): True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.alpha_c, alpha_c, abs(true_params.alpha_c - alpha_c));
    fprintf('  Contrast (Beta) : True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.beta_c, beta_c, abs(true_params.beta_c - beta_c));
    fprintf('  EqNoise (Signal): True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.signal_fw, signal_fw, abs(true_params.signal_fw - signal_fw));
    fprintf('  EqNoise (Sig_int):True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.sigma_int_fw, sigma_int_fw, abs(true_params.sigma_int_fw - sigma_int_fw));
    fprintf('--------------------------------------------------\n');

    %% Extract 3 Levels
    K = (target_levels - gamma) ./ (1 - gamma - lambda);
    calib_contrast = alpha_c .* (-log(1 - K)).^(1/beta_c);
    calib_contrast = max(0.01, min(1.0, calib_contrast));

    P_target_adj = (target_levels - lambda * 0.5) ./ (1 - lambda);
    d_target = sqrt(2) .* norminv(P_target_adj);
    calib_filter = zeros(1, length(target_levels));

    for i = 1:length(target_levels)
        val = (signal_fw / d_target(i))^2 - sigma_int_fw^2;
        if val < 0
            calib_filter(i) = 0;
        else
            calib_filter(i) = sqrt(val);
        end
    end
    
    % The arrays are naturally created ordered by target_levels: 
    % contrast: [65% val, 75% val, 85% val] -> ascending 
    % filter: [65% val, 75% val, 85% val] -> descending (wider filter = harder)
    % We want to reverse the filter array so that it is ALSO strictly ascending physical values
    % (i.e., narrowest/easiest filter first, widest/hardest filter last)
    calib_filter = sort(calib_filter, 'ascend');
    
    % We must also reverse the target levels when plotting so that they map correctly to the reversed filter array
    filter_target_levels = flip(target_levels);

    %% Plot
    figure('Color', 'w', 'Position', [100, 100, 1000, 400], 'Name', ['S' subj_ID ' Simulation Fits']);

    % Contrast
    subplot(1,2,1);
    plot(unique_c, prop_c, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
    x_fit_c = linspace(0, max(unique_c)*1.2, 100);
    plot(x_fit_c, weibull_prob(x_fit_c, alpha_c, beta_c), 'r-', 'LineWidth', 2);
    plot(x_fit_c, weibull_prob(x_fit_c, true_params.alpha_c, true_params.beta_c), 'g--', 'LineWidth', 2);
    
    for i = 1:3
        plot([0 calib_contrast(i)], [target_levels(i) target_levels(i)], 'b--', 'HandleVisibility', 'off');
        plot([calib_contrast(i) calib_contrast(i)], [0 target_levels(i)], 'b--', 'HandleVisibility', 'off');
        plot(calib_contrast(i), target_levels(i), 'b*', 'MarkerSize', 8);
    end
    
    xlim([0 max(max(unique_c), max(calib_contrast))*1.1]);
    ylim([0.4 1.0]);
    xlabel('Contrast'); ylabel('Proportion Correct');
    title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f)', alpha_c, beta_c));
    legend('Simulated Data', 'Fit', 'True Generative', 'Location', 'SouthEast');
    set(gca, 'TickDir', 'out', 'TickLength', [0.015, 0.015]);
    box off; axis square;

    % Filter Width
    subplot(1,2,2);
    plot(unique_fw, prop_fw, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
    x_fit_fw = linspace(0, max(unique_fw)*1.2, 100);
    plot(x_fit_fw, en_prob(x_fit_fw, signal_fw, sigma_int_fw), 'r-', 'LineWidth', 2);
    plot(x_fit_fw, en_prob(x_fit_fw, true_params.signal_fw, true_params.sigma_int_fw), 'g--', 'LineWidth', 2);
    
    for i = 1:3
        plot([0 calib_filter(i)], [filter_target_levels(i) filter_target_levels(i)], 'b--', 'HandleVisibility', 'off');
        plot([calib_filter(i) calib_filter(i)], [0 filter_target_levels(i)], 'b--', 'HandleVisibility', 'off');
        plot(calib_filter(i), filter_target_levels(i), 'b*', 'MarkerSize', 8);
    end
    
    xlim([0 max(max(unique_fw), max(calib_filter))*1.1]);
    ylim([0.4 1.0]);
    xlabel('Filter Width (\sigma_{ext})'); ylabel('Proportion Correct');
    title(sprintf('Eq Noise Fit (Sig=%.2f, \\sigma_{int}=%.2f)', signal_fw, sigma_int_fw));
    legend('Simulated Data', 'Fit', 'True Generative', 'Location', 'SouthEast');
    set(gca, 'TickDir', 'out', 'TickLength', [0.015, 0.015]);
    box off; axis square;
    
    %% Save Figure
    fig_dir = [script_dir '/figures'];
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    
    % Save as PDF
    fig_filename = [fig_dir '/S' subj_ID '_Calibration_Simulation.pdf'];
    exportgraphics(gcf, fig_filename, 'ContentType', 'vector');
        
end
