% validate_calibration.m
% Validates the calibration fits using the simulated data

clear all; close all; clc;

%% Set directories
script_dir = fileparts(mfilename('fullpath'));
data_dir = 'data';
addpath('../analyses/plotting');

% Load plot settings
ps = plotSettings();

num_subjs = 20;
subj_IDs = cell(num_subjs,1);
for subj = 1:num_subjs
    subj_IDs{subj} = num2str(9900 + subj - 1);
end

target_levels = [0.65, 0.75, 0.85];
guess_rate = 0.5;
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

    % Compute proportions and standard errors
    [unique_c, ~, idx_c] = unique(all_contrast_levels);
    n_c = accumarray(idx_c, 1);
    k_c = accumarray(idx_c, all_contrast_correct, [], @sum);
    prop_c = k_c ./ n_c;
    se_c = sqrt(prop_c .* (1 - prop_c) ./ n_c); % Binomial standard error

    [unique_fw, ~, idx_fw] = unique(all_filter_levels);
    n_fw = accumarray(idx_fw, 1);
    k_fw = accumarray(idx_fw, all_filter_correct, [], @sum);
    prop_fw = k_fw ./ n_fw;
    se_fw = sqrt(prop_fw .* (1 - prop_fw) ./ n_fw); % Binomial standard error

    %% Fit Models
    % Both features use Weibull. For filter width, we transform x -> 1/x.
    weibull_prob = @(x, alpha, beta) guess_rate + (1 - guess_rate - lambda) * (1 - exp(-(x./alpha).^beta));

    % --- A. Contrast (Weibull) ---
    nll_contrast = @(params) -sum( ...
        k_c .* log(max(eps, weibull_prob(unique_c, params(1), params(2)))) + ...
        (n_c - k_c) .* log(max(eps, 1 - weibull_prob(unique_c, params(1), params(2)))) ...
        );

    lb_c = [0.001, 0.5];
    ub_c = [1.0, 5.0];
    best_guess_c = gridSearch(nll_contrast, lb_c, ub_c);
    best_params_c = fmincon(nll_contrast, best_guess_c, [], [], [], [], lb_c, ub_c, [], options);
    alpha_c = best_params_c(1);
    beta_c = best_params_c(2);

    % --- B. Filter Width (Weibull on precision = 1/fw) ---
    unique_precision = 1 ./ unique_fw;

    nll_filter = @(params) -sum( ...
        k_fw .* log(max(eps, weibull_prob(unique_precision, params(1), params(2)))) + ...
        (n_fw - k_fw) .* log(max(eps, 1 - weibull_prob(unique_precision, params(1), params(2)))) ...
        );

    lb_fw = [0.005, 0.5];
    ub_fw = [0.20, 5.0];
    best_guess_fw = gridSearch(nll_filter, lb_fw, ub_fw);
    best_params_fw = fmincon(nll_filter, best_guess_fw, [], [], [], [], lb_fw, ub_fw, [], options);
    alpha_fw = best_params_fw(1);
    beta_fw  = best_params_fw(2);

    %% Calculate Errors
    fprintf('Subject %s Parameters Validation:\n', subj_ID);
    fprintf('  Contrast (Alpha):   True=%.4f, Fit=%.4f, Err=%.4f\n', true_params.alpha_c, alpha_c, abs(true_params.alpha_c - alpha_c));
    fprintf('  Contrast (Beta) :   True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.beta_c, beta_c, abs(true_params.beta_c - beta_c));
    fprintf('  Filter (Alpha_fw):  True=%.4f, Fit=%.4f, Err=%.4f\n', true_params.alpha_fw, alpha_fw, abs(true_params.alpha_fw - alpha_fw));
    fprintf('  Filter (Beta_fw) :  True=%.3f, Fit=%.3f, Err=%.3f\n', true_params.beta_fw, beta_fw, abs(true_params.beta_fw - beta_fw));
    fprintf('--------------------------------------------------\n');

    %% Extract 3 Levels

    K = (target_levels - guess_rate) ./ (1 - guess_rate - lambda);
    calib_contrast = alpha_c .* (-log(1 - K)).^(1/beta_c);
    calib_contrast = max(0.01, min(1.0, calib_contrast));

    % Filter width: invert Weibull in precision domain, then convert back
    precision_targets = alpha_fw .* (-log(1 - K)).^(1/beta_fw);
    calib_filter = 1 ./ precision_targets;

    % Clamp filter width to physical limits
    calib_filter = max(10.0, min(90.0, calib_filter));
    precision_targets = 1 ./ calib_filter;

    calib_filter = sort(calib_filter, 'ascend');
    filter_target_levels = flip(target_levels);

    %% Plot
    figure('Color', 'w', 'Position', [100, 100, 1000, 400], 'Name', ['S' subj_ID ' Simulation Fits']);

    % Contrast (log-scaled)
    subplot(1,2,1);
    errorbar(unique_c, prop_c, se_c, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    x_fit_c = logspace(log10(min(unique_c)*0.7), log10(max(unique_c)*1.3), 100);
    plot(x_fit_c, weibull_prob(x_fit_c, alpha_c, beta_c), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    plot(x_fit_c, weibull_prob(x_fit_c, true_params.alpha_c, true_params.beta_c), '--', 'Color', ps.colors.green, 'LineWidth', 2);

    xl_c = xlim;
    for i = 1:3
        plot([xl_c(1) calib_contrast(i)], [target_levels(i) target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([calib_contrast(i) calib_contrast(i)], [0.0 target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(calib_contrast(i), target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end

    set(gca, 'XScale', 'log');
    xlim([min(unique_c)*0.7 max(unique_c)*1.3]);
    ylim([0.0 1.0]);
    xticks([0.1 0.2 0.4 0.6 0.8 1.0]);
    xticklabels({'0.1', '0.2', '0.4', '0.6', '0.8', '1.0'});
    xlabel('Contrast', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f)', alpha_c, beta_c), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    legend('Simulated Data', 'Fit', 'True Generative', 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    % Filter Width (Weibull on precision = 1/fw, log-scaled)
    subplot(1,2,2);
    errorbar(unique_precision, prop_fw, se_fw, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    x_fit_prec = logspace(log10(min(unique_precision)*0.7), log10(max(unique_precision)*1.3), 100);
    plot(x_fit_prec, weibull_prob(x_fit_prec, alpha_fw, beta_fw), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    plot(x_fit_prec, weibull_prob(x_fit_prec, true_params.alpha_fw, true_params.beta_fw), '--', 'Color', ps.colors.green, 'LineWidth', 2);

    xl = xlim;
    for i = 1:3
        prec_i = 1 / calib_filter(i);
        plot([xl(1) prec_i], [filter_target_levels(i) filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([prec_i prec_i], [0.0 filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(prec_i, filter_target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end

    set(gca, 'XScale', 'log');
    xlim([min(unique_precision)*0.7 max(unique_precision)*1.3]);
    ylim([0.0 1.0]);

    prec_ticks = [1/80, 1/40, 1/20, 1/10];
    xticks(prec_ticks);
    xticklabels({'1/80', '1/40', '1/20', '1/10'});

    xlabel('Precision (1/filter width)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.4f, \\beta=%.2f)', alpha_fw, beta_fw), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    %% Save Figure
    fig_dir = fullfile(script_dir, 'experiment', 'figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    % Save as PDF
    fig_filename = fullfile(fig_dir, sprintf('S%s_Calibration_Simulation.pdf', subj_ID));

    % Use print to maintain original figure size and margins
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperOrientation', 'landscape');
    print(gcf, fig_filename, '-dpdf', '-bestfit');

end
