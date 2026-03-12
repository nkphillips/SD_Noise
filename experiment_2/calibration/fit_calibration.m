% fit_calibration.m
% Extracts 3 subject-specific levels from calibration data.

clear all; close all; clc;

%% Set directories
addpath('../analyses/plotting');

% Load plot settings
ps = plotSettings();

%% 1. Set subject IDs and directories
% You can add multiple subject IDs to this cell array
subj_IDs = {'999'};

data_base_dir = '../../data/';

for s = 1:length(subj_IDs)
    subj_ID = subj_IDs{s};
    data_dir = [data_base_dir subj_ID];

    disp(['Processing calibration for Subject ' subj_ID '...']);

    % Load all calibration runs for this subject
    file_pattern = [data_dir '/SD_Noise_Exp2_Calibration_S' subj_ID '_Run*.mat'];
    files = dir(file_pattern);

    if isempty(files)
        warning('No calibration data found for subject %s. Skipping...', subj_ID);
        continue;
    end

    %% 2. Aggregate Data
    all_contrast_correct = [];
    all_contrast_levels = [];

    all_filter_correct = [];
    all_filter_levels = [];

    for f = 1:length(files)
        load([data_dir '/' files(f).name], 'run_info');
        p = run_info.p;
        behav_data = run_info.behav_data;

        for n_block = 1:p.num_blocks
            curr_feature = p.feature_order(n_block);

            % Calculate probe offset to exclude 0-offset trials from proportion correct
            test_ori = p.trial_events(:, 1, n_block);
            probe_ori = p.trial_events(:, 2, n_block);
            probe_offset = probe_ori - test_ori;

            % Wrap offset
            probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
            probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;

            % Exclude 0-offset trials because guessing brings the asymptote down to ~0.92,
            % which severely biases the Weibull/EqNoise fits (which expect a 1.0 asymptote).
            valid_trials = probe_offset ~= 0;

            % p.trial_events structure: [test_orientations, probe_orientations, level_order]
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


    % Compute proportions and standard errors
    [unique_c, ~, idx_c] = unique(all_contrast_levels);
    n_c = accumarray(idx_c, 1);
    k_c = accumarray(idx_c, all_contrast_correct, [], @sum);
    prop_c = k_c ./ n_c;
    se_c = sqrt(prop_c .* (1 - prop_c) ./ n_c);

    [unique_fw, ~, idx_fw] = unique(all_filter_levels);
    n_fw = accumarray(idx_fw, 1);
    k_fw = accumarray(idx_fw, all_filter_correct, [], @sum);
    prop_fw = k_fw ./ n_fw;
    se_fw = sqrt(prop_fw .* (1 - prop_fw) ./ n_fw);

    %% 3. Fit Models (Maximum Likelihood Estimation)

    % Fixed parameters
    gamma = 0.5;   % Guess rate (2AFC)
    lambda = 0.01; % Lapse rate
    options = optimset('Display', 'off');

    % --- A. Contrast (Weibull) ---
    % P(x) = gamma + (1 - gamma - lambda) * (1 - exp(-(x/alpha)^beta))
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

    % fmincon for Contrast
    % bounds: alpha (threshold) in [0.001, 1.0], beta (slope) in [0.5, 5.0]
    lb_c = [0.001, 0.5];
    ub_c = [1.0, 5.0];
    best_params_c = fmincon(nll_weibull, best_guess_c, [], [], [], [], lb_c, ub_c, [], options);
    alpha_c = best_params_c(1);
    beta_c = best_params_c(2);

    % --- B. Filter Width (Equivalent Noise) ---
    % d' = Signal / sqrt(x^2 + sigma_int^2)
    % P(x) = (1 - lambda) * normcdf(d' / sqrt(2)) + lambda * 0.5
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

    % fmincon for Filter Width
    % bounds: signal in [1.0, 50.0], sigma_int in [1.0, 30.0]
    lb_fw = [1.0, 1.0];
    ub_fw = [50.0, 30.0];
    best_params_fw = fmincon(nll_en, best_guess_fw, [], [], [], [], lb_fw, ub_fw, [], options);
    signal_fw = best_params_fw(1);
    sigma_int_fw = best_params_fw(2);

    %% 4. Inverse Steps to Find 3 Levels

    target_levels = [0.65, 0.75, 0.85];

    % --- A. Contrast Inverse ---
    % K = (P_target - gamma) / (1 - gamma - lambda)
    % x = alpha * [-ln(1 - K)]^(1/beta)

    K = (target_levels - gamma) ./ (1 - gamma - lambda);
    calib_contrast = alpha_c .* (-log(1 - K)).^(1/beta_c);

    % Restrict to reasonable values
    calib_contrast = max(0.01, min(1.0, calib_contrast));

    % --- B. Filter Width Inverse ---
    % P_target = (1 - lambda) * normcdf(...) + lambda * 0.5
    % P_target_adj = (P_target - lambda * 0.5) / (1 - lambda)
    % d'_target = sqrt(2) * norminv(P_target_adj)
    % x = sqrt( (Signal / d'_target)^2 - sigma_int^2 )

    P_target_adj = (target_levels - lambda * 0.5) ./ (1 - lambda);
    d_target = sqrt(2) .* norminv(P_target_adj);

    calib_filter = zeros(1, length(target_levels));
    for i = 1:length(target_levels)
        val = (signal_fw / d_target(i))^2 - sigma_int_fw^2;
        if val < 0
            warning('Target performance %.2f unreachable (internal noise too high). Clamping to 0.', target_levels(i));
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

    %% 5. Save Levels

    calib.target_levels = target_levels;
    calib.contrast_levels = calib_contrast;
    calib.filter_width_levels = calib_filter;
    calib.fit_params.contrast_alpha = alpha_c;
    calib.fit_params.contrast_beta = beta_c;
    calib.fit_params.filter_signal = signal_fw;
    calib.fit_params.filter_sigma_int = sigma_int_fw;

    save([data_dir '/S' subj_ID '_calibrated_levels.mat'], 'calib');
    disp(['Successfully saved calibrated levels to: ' data_dir '/S' subj_ID '_calibrated_levels.mat']);

    %% 6. Plot the Fits

    figure('Color', 'w', 'Position', [100, 100, 1000, 400], 'Name', ['S' subj_ID ' Calibration Fits']);

    % Plot Contrast
    subplot(1,2,1);
    errorbar(unique_c, prop_c, se_c, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    x_fit_c = linspace(0, max(unique_c), 100);
    plot(x_fit_c, weibull_prob(x_fit_c, alpha_c, beta_c), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    for i = 1:3
        plot([0 calib_contrast(i)], [target_levels(i) target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([calib_contrast(i) calib_contrast(i)], [0 target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(calib_contrast(i), target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end
    xlim([0 max(max(unique_c), max(calib_contrast))*1.1]);
    ylim([0.4 1.0]);
    xlabel('Contrast', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f)', alpha_c, beta_c), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    legend('Data', 'Fit', 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    % Plot Filter Width
    subplot(1,2,2);
    errorbar(unique_fw, prop_fw, se_fw, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    x_fit_fw = linspace(0, max(unique_fw), 100);
    plot(x_fit_fw, en_prob(x_fit_fw, signal_fw, sigma_int_fw), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    for i = 1:3
        plot([0 calib_filter(i)], [filter_target_levels(i) filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([calib_filter(i) calib_filter(i)], [0 filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(calib_filter(i), filter_target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end
    xlim([0 max(max(unique_fw), max(calib_filter))*1.1]);
    ylim([0.4 1.0]);
    xlabel('Filter Width (\sigma_{ext})', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Eq Noise Fit (Sig=%.2f, \\sigma_{int}=%.2f)', signal_fw, sigma_int_fw), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;
    
    drawnow; % ensure figures are rendered before continuing the loop
    
    %% Save Figure
    fig_dir = [script_dir '/figures'];
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    
    % Save as PDF
    fig_filename = [fig_dir '/S' subj_ID '_Calibration_Fit.pdf'];
    
    % Use exportgraphics to prevent clipping and ensure consistent sizing
    if exist('exportgraphics', 'file')
        exportgraphics(gcf, fig_filename, 'ContentType', 'vector');
    else
        set(gcf, 'PaperPositionMode', 'auto');
        set(gcf, 'PaperOrientation', 'landscape');
        print(gcf, fig_filename, '-dpdf', '-bestfit');
    end

end % End of subject loop
