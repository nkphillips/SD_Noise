% fit_calibration.m
% Extracts 3 subject-specific levels from calibration data.

clear all; close all; clc;

%% Set directories
script_dir = pwd;
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

            % Exclude 0-offset trials because guessing brings the asymptote down
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
    n_c = accumarray(idx_c, 1); % count number of trials per contrast lvl
    k_c = accumarray(idx_c, all_contrast_correct, [], @sum); % count number of correct trials per contrast lvl
    prop_c = k_c ./ n_c; % calculate proportion of correct trials
    se_c = sqrt(prop_c .* (1 - prop_c) ./ n_c); %

    [unique_fw, ~, idx_fw] = unique(all_filter_levels);
    n_fw = accumarray(idx_fw, 1);
    k_fw = accumarray(idx_fw, all_filter_correct, [], @sum);
    prop_fw = k_fw ./ n_fw;
    se_fw = sqrt(prop_fw .* (1 - prop_fw) ./ n_fw);

    %% 3. Fit Models (Maximum Likelihood Estimation)
    % 1. Contrast: Weibull function
    % 2. Filter Width: Perceptual Template Model (PTM; Lu & Dosher, 2008)

    % Fixed parameters
    gamma = 0.5;      % Guess rate (2AFC)
    lambda = 0.01;    % Lapse rate
    ptm_gamma = 2.0;  % Fixed transducer exponent (Lu & Dosher, 1999)
    g = ptm_gamma;
    options = optimset('Display', 'off');

    % --- A. Contrast (Weibull) ---
    % P(x) = gamma + (1 - gamma - lambda) * (1 - exp(-(x/alpha)^beta))
    weibull_prob = @(x, alpha, beta) gamma + (1 - gamma - lambda) * (1 - exp(-(x./alpha).^beta));

    nll_weibull = @(params) -sum( ...
        k_c .* log(max(eps, weibull_prob(unique_c, params(1), params(2)))) + ...
        (n_c - k_c) .* log(max(eps, 1 - weibull_prob(unique_c, params(1), params(2)))) ...
        );

    lb_c = [0.001, 0.5];
    ub_c = [1.0, 5.0];
    best_guess_c = gridSearch(nll_weibull, lb_c, ub_c);
    best_params_c = fmincon(nll_weibull, best_guess_c, [], [], [], [], lb_c, ub_c, [], options);
    alpha_c = best_params_c(1);
    beta_c = best_params_c(2);

    % --- B. Filter Width (PTM) ---
    % d' = Signal^g / sqrt((1+Nm^2)*fw^(2g) + Nm^2*Signal^(2g) + Na^2)
    % P(x) = (1-lambda)*normcdf(d'/sqrt(2)) + lambda*0.5
    % Params: [Signal, N_mul, N_add]
    ptm_dprime = @(x, S, Nm, Na) S.^g ./ sqrt((1+Nm.^2).*x.^(2*g) + Nm.^2.*S.^(2*g) + Na.^2);
    ptm_prob = @(x, S, Nm, Na) (1 - lambda) * normcdf(ptm_dprime(x, S, Nm, Na) / sqrt(2)) + lambda * 0.5;

    nll_ptm = @(params) -sum( ...
        k_fw .* log(max(eps, ptm_prob(unique_fw, params(1), params(2), params(3)))) + ...
        (n_fw - k_fw) .* log(max(eps, 1 - ptm_prob(unique_fw, params(1), params(2), params(3)))) ...
        );

    % bounds: [Signal, N_mul, N_add]
    lb_fw = [5.0, 0.001, 1.0];
    ub_fw = [80.0, 0.65, 600.0];
    best_guess_fw = gridSearch(nll_ptm, lb_fw, ub_fw);
    best_params_fw = fmincon(nll_ptm, best_guess_fw, [], [], [], [], lb_fw, ub_fw, [], options);
    signal_fw = best_params_fw(1);
    nmul_fw   = best_params_fw(2);
    nadd_fw   = best_params_fw(3);

    %% 4. Inverse Steps to Find 3 Levels

    target_levels = [0.65, 0.75, 0.85];

    % --- A. Contrast Inverse ---
    K = (target_levels - gamma) ./ (1 - gamma - lambda);
    calib_contrast = alpha_c .* (-log(1 - K)).^(1/beta_c);
    calib_contrast = max(0.01, min(1.0, calib_contrast));

    % --- B. Filter Width Inverse (PTM) ---
    % fw^(2g) = (S^(2g) * (1/d'^2 - Nm^2) - Na^2) / (1 + Nm^2)
    P_target_adj = (target_levels - lambda * 0.5) ./ (1 - lambda);
    d_target = sqrt(2) .* norminv(P_target_adj);

    S2g = signal_fw^(2*g);

    calib_filter = zeros(1, length(target_levels));
    for i = 1:length(target_levels)
        if 1/d_target(i)^2 <= nmul_fw^2
            warning('Target %.2f unreachable (d'' exceeds 1/N_mul). Clamping to 0.', target_levels(i));
            calib_filter(i) = 0;
            continue;
        end
        numer = S2g * (1/d_target(i)^2 - nmul_fw^2) - nadd_fw^2;
        denom = 1 + nmul_fw^2;
        if numer <= 0
            warning('Target %.2f unreachable (noise too high). Clamping to 0.', target_levels(i));
            calib_filter(i) = 0;
        else
            calib_filter(i) = (numer / denom)^(1/(2*g));
        end
    end

    calib_filter = sort(calib_filter, 'ascend');
    filter_target_levels = flip(target_levels);

    %% 5. Save Levels

    calib.target_levels = target_levels;
    calib.contrast_levels = calib_contrast;
    calib.filter_width_levels = calib_filter;
    calib.fit_params.contrast_alpha = alpha_c;
    calib.fit_params.contrast_beta = beta_c;
    calib.fit_params.filter_signal = signal_fw;
    calib.fit_params.filter_nmul = nmul_fw;
    calib.fit_params.filter_nadd = nadd_fw;
    calib.fit_params.filter_gamma = ptm_gamma;

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
    xlim([0 max(unique_c)*1.1]);
    ylim([0.4 1.0]);
    xlabel('Contrast', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f)', alpha_c, beta_c), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    legend('Data', 'Fit', 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    % Plot Filter Width
    subplot(1,2,2);
    errorbar(unique_fw, prop_fw, se_fw, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    x_fit_fw = linspace(0.1, max(unique_fw), 100);
    plot(x_fit_fw, ptm_prob(x_fit_fw, signal_fw, nmul_fw, nadd_fw), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    for i = 1:3
        plot([0 calib_filter(i)], [filter_target_levels(i) filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([calib_filter(i) calib_filter(i)], [0 filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(calib_filter(i), filter_target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end
    xlim([0 max(unique_fw)*1.1]);
    ylim([0.4 1.0]);
    xlabel('Filter Width (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('PTM Fit (S=%.1f, N_m=%.2f, N_a=%.0f)', signal_fw, nmul_fw, nadd_fw), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
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

end
