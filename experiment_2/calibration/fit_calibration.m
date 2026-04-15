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
subj_IDs = {'777', '000'};

data_base_dir = '../data';

for s = 1:length(subj_IDs)
    subj_ID = subj_IDs{s};
    data_dir = fullfile(data_base_dir, subj_ID);

    disp(['Processing calibration for Subject ' subj_ID '...']);

    % Load all calibration runs for this subject
    file_pattern = fullfile(data_dir, ['SD_Noise_Exp2_Calibration_S' subj_ID '_Run*.mat']);
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
        load(fullfile(data_dir, files(f).name), 'run_info');
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
    %
    % Both features use the same Weibull psychometric function:
    %   P(x) = gamma + (1 - gamma - lambda) * (1 - exp(-(x/alpha)^beta))
    %
    % For CONTRAST this is straightforward: higher contrast -> better
    % performance, so we fit the Weibull directly on contrast values.
    %
    % For FILTER WIDTH there is a complication: wider filters inject more
    % orientation noise, so performance DROPS as filter width INCREASES.
    % The Weibull expects performance to RISE with x.
    %
    % THE FIX: We transform the x-axis from filter_width to "precision",
    % defined as:
    %       precision = 1 / filter_width
    %
    % Precision is high when the filter is narrow (easy, high accuracy)
    % and low when the filter is wide (hard, low accuracy). This makes
    % precision monotonically increasing with accuracy, exactly what the
    % Weibull needs.
    %
    % After fitting and inverting the Weibull in the precision domain, we
    % convert back to physical filter widths via:
    %       filter_width = 1 / precision

    % Fixed parameters (same for both features)
    gamma = 0.5;   % Guess rate for a 2-AFC task (chance = 50%)
    lambda = 0.01; % Lapse rate (caps upper asymptote at 99%)
    options = optimset('Display', 'off');

    % The Weibull function shared by both features.
    % x = contrast (Feature 1) or precision = 1/fw (Feature 2).
    % alpha = threshold, beta = slope.
    weibull_prob = @(x, alpha, beta) gamma + (1 - gamma - lambda) * (1 - exp(-(x./alpha).^beta));

    % ---------------------------------------------------------------
    % --- A. Contrast (Weibull on raw contrast values) ---
    % ---------------------------------------------------------------
    nll_contrast = @(params) -sum( ...
        k_c .* log(max(eps, weibull_prob(unique_c, params(1), params(2)))) + ...
        (n_c - k_c) .* log(max(eps, 1 - weibull_prob(unique_c, params(1), params(2)))) ...
        );

    % Bounds: alpha (threshold) in [0.0100, 1.0000], beta (slope) in [0.5, 5.0]
    lb_c = [0.0100, 0.5000];
    ub_c = [1.0000, 5.0000];
    best_guess_c = gridSearch(nll_contrast, lb_c, ub_c);
    best_params_c = fmincon(nll_contrast, best_guess_c, [], [], [], [], lb_c, ub_c, [], options);
    alpha_c = best_params_c(1);
    beta_c  = best_params_c(2);

    % ---------------------------------------------------------------
    % --- B. Filter Width (Weibull on transformed precision = 1/fw) ---
    % ---------------------------------------------------------------

    % STEP 1: Transform filter widths to precision.
    % Example: fw = [10, 20, 40, 60, 80] deg
    %       -> precision = [0.100, 0.050, 0.025, 0.017, 0.013]
    % Now larger values = narrower filter = easier = higher accuracy.
    unique_precision = 1 ./ unique_fw;

    % The NLL uses the SAME Weibull function, just evaluated on precision
    % instead of filter width.
    nll_filter = @(params) -sum( ...
        k_fw .* log(max(eps, weibull_prob(unique_precision, params(1), params(2)))) + ...
        (n_fw - k_fw) .* log(max(eps, 1 - weibull_prob(unique_precision, params(1), params(2)))) ...
        );

    % Bounds for the precision-domain Weibull.
    % alpha_fw: threshold in precision units. Filter widths range from
    %   2 to 180 deg, so precision ranges from 1/180 = 0.0056 to 1/2 = 0.5.
    %   A reasonable threshold sits somewhere in that range.
    % beta_fw: slope, same general range as contrast.
    lb_fw = [0.0050, 0.5000];
    ub_fw = [0.6000, 5.0000];
    best_guess_fw = gridSearch(nll_filter, lb_fw, ub_fw);
    best_params_fw = fmincon(nll_filter, best_guess_fw, [], [], [], [], lb_fw, ub_fw, [], options);
    alpha_fw = best_params_fw(1);
    beta_fw  = best_params_fw(2);

    %% 4. Inverse Steps to Find 3 Levels
    %
    % We analytically invert the Weibull to find the stimulus value that
    % produces each target accuracy. The inverse Weibull is:
    %
    %   K = (P_target - gamma) / (1 - gamma - lambda)
    %   x_target = alpha * [-ln(1 - K)]^(1/beta)
    %
    % For contrast, x_target IS the physical contrast.
    % For filter width, x_target is in precision units, so we convert:
    %   fw_target = 1 / x_target

    target_levels = [0.65, 0.75, 0.85]; % [0.65, 0.75, 0.85]

    % Intermediate quantity K (same formula for both features)
    K = (target_levels - gamma) ./ (1 - gamma - lambda);

    % --- A. Contrast Inverse ---
    calib_contrast = alpha_c .* (-log(1 - K)).^(1/beta_c);
    calib_contrast = max(0.01, min(1.0, calib_contrast));

    % --- B. Filter Width Inverse ---
    % STEP 1: Invert the Weibull in the precision domain.
    precision_targets = alpha_fw .* (-log(1 - K)).^(1/beta_fw);

    % STEP 2: Convert precision back to physical filter width.
    %   precision = 1/fw  =>  fw = 1/precision
    calib_filter = 1 ./ precision_targets;

    % --- SANITY CHECKS & CLAMPING FOR REAL DATA ---
    % Warn if the subject's parameters exceed physical monitor/stimulus bounds
    if any(calib_contrast >= 1.0)
        warning('S%s: Contrast target >= 1.0 required. Clamping to 1.0. Subject may have performed poorly.', subj_ID);
    end
    if any(calib_filter >= 180)
        warning('S%s: Filter width target >= 180 deg required. Clamping to 180. Subject may have performed poorly.', subj_ID);
    end
    if any(calib_filter <= 2)
        warning('S%s: Filter width target <= 2 deg required. Clamping to 2. Subject performed extremely well.', subj_ID);
    end

    % Clamp filter width to physically meaningful limits [2 deg, 180 deg]
    % (We clamp to 2 minimum because extremely narrow filters
    % can cause artifacting/aliasing in the spatial frequency domain).
    calib_filter = max(2.0, min(180.0, calib_filter));

    % The result naturally comes out DESCENDING (highest precision = narrowest
    % filter comes first for the hardest target). Sort ascending so that
    % index 1 = narrowest (easiest), index 3 = widest (hardest).
    calib_filter = sort(calib_filter, 'ascend');

    % Re-calculate clamped precision targets so the plot markers stay on the bounds
    % and correspond to the newly sorted calib_filter indices
    precision_targets = 1 ./ calib_filter;

    % The target_levels [0.65, 0.75, 0.85] map to descending filter widths,
    % so after sorting we flip the labels for plotting.
    filter_target_levels = flip(target_levels);

    % Calculate R^2 for Contrast and Filter Width fits
    % Note: R^2 can be negative if the constrained Weibull fit (forced to 0.5
    % guess rate) has a higher sum of squared errors than a flat line at the mean.
    % We clamp it to 0 in these cases.
    ss_res_c = sum((prop_c - weibull_prob(unique_c, alpha_c, beta_c)).^2);
    ss_tot_c = sum((prop_c - mean(prop_c)).^2);
    r2_c = max(0, 1 - (ss_res_c / ss_tot_c));

    ss_res_fw = sum((prop_fw - weibull_prob(unique_precision, alpha_fw, beta_fw)).^2);
    ss_tot_fw = sum((prop_fw - mean(prop_fw)).^2);
    r2_fw = max(0, 1 - (ss_res_fw / ss_tot_fw));

    %% 5. Save Levels

    calib.target_levels = target_levels;
    calib.contrast_levels = calib_contrast;
    calib.filter_width_levels = calib_filter;
    calib.fit_params.contrast_alpha = alpha_c;
    calib.fit_params.contrast_beta = beta_c;
    calib.fit_params.filter_alpha = alpha_fw;  % threshold in precision (1/deg) units
    calib.fit_params.filter_beta = beta_fw;

    save(fullfile(data_dir, ['S' subj_ID '_calibrated_levels.mat']), 'calib');
    disp(['Successfully saved calibrated levels to: ' fullfile(data_dir, ['S' subj_ID '_calibrated_levels.mat'])]);

    %% 6. Plot the Fits

    figure('Color', 'w', 'Position', [100, 100, 1000, 400], 'Name', ['S' subj_ID ' Calibration Fits']);

    % Plot Contrast (log-scaled)
    subplot(1,2,1);
    errorbar(unique_c, prop_c, se_c, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    % Determine axis limits dynamically to include both data and extrapolated targets
    min_x_c = min([unique_c; calib_contrast(:)]);
    max_x_c = max([unique_c; calib_contrast(:)]);
    
    x_fit_c = logspace(log10(min_x_c*0.7), log10(max_x_c*1.3), 100);
    plot(x_fit_c, weibull_prob(x_fit_c, alpha_c, beta_c), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    xl_c = xlim;
    for i = 1:3
        plot([xl_c(1) calib_contrast(i)], [target_levels(i) target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([calib_contrast(i) calib_contrast(i)], [0.0 target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(calib_contrast(i), target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end
    set(gca, 'XScale', 'log');
    xlim([min_x_c*0.7 max_x_c*1.3]);
    ylim([0.0 1.0]);
    xticks([0.1 0.2 0.4 0.6 0.8 1.0]);
    xticklabels({'0.1', '0.2', '0.4', '0.6', '0.8', '1.0'});
    xlabel('Contrast', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f, R^2=%.3f)', alpha_c, beta_c, r2_c), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    legend('Data', 'Fit', 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    % Plot Filter Width (on the transformed precision = 1/fw axis, log-scaled)
    subplot(1,2,2);
    errorbar(unique_precision, prop_fw, se_fw, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.black, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
    % Determine axis limits dynamically to include both data and extrapolated targets
    min_x_prec = min([unique_precision; precision_targets(:)]);
    max_x_prec = max([unique_precision; precision_targets(:)]);

    x_fit_prec = logspace(log10(min_x_prec*0.7), log10(max_x_prec*1.3), 100);
    plot(x_fit_prec, weibull_prob(x_fit_prec, alpha_fw, beta_fw), '-', 'Color', ps.colors.red, 'LineWidth', 2);
    xl = xlim;
    for i = 1:3
        prec_i = precision_targets(i);
        plot([xl(1) prec_i], [filter_target_levels(i) filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot([prec_i prec_i], [0.0 filter_target_levels(i)], '--', 'Color', ps.colors.blue, 'HandleVisibility', 'off');
        plot(prec_i, filter_target_levels(i), '*', 'Color', ps.colors.blue, 'MarkerSize', 8);
    end
    set(gca, 'XScale', 'log');
    xlim([min_x_prec*0.7 max_x_prec*1.3]);
    ylim([0.0 1.0]);

    % Dynamically generate xticks for precision based on the plotted range
    possible_fws = [180, 140, 100, 80, 60, 40, 20, 10, 5, 2];
    valid_fws = possible_fws(1./possible_fws >= min_x_prec*0.7 & 1./possible_fws <= max_x_prec*1.3);
    
    % If for some reason none fell perfectly in range, grab the min and max bounds
    if isempty(valid_fws)
        valid_fws = [round(1/min_x_prec), round(1/max_x_prec)];
    end
    
    prec_ticks = sort(1 ./ valid_fws, 'ascend');
    
    % Sort valid_fws descending to match the ascending precision ticks
    valid_fws = sort(valid_fws, 'descend');
    
    tick_labels = cell(1, length(valid_fws));
    for i = 1:length(valid_fws)
        tick_labels{i} = ['1/' num2str(valid_fws(i))];
    end
    
    xticks(prec_ticks);
    xticklabels(tick_labels);

    xlabel('Precision (1/filter width)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type); ylabel('Proportion Correct', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('Weibull Fit (\\alpha=%.4f, \\beta=%.2f, R^2=%.2f)', alpha_fw, beta_fw, r2_fw), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;

    drawnow; % ensure figures are rendered before continuing the loop

    %% Save Figure
    fig_dir = fullfile(script_dir, 'figures');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    % Save as PDF
    fig_filename = fullfile(fig_dir, ['S' subj_ID '_Calibration_Fit.pdf']);

    % Use print to maintain original figure size and margins
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperOrientation', 'landscape');
    print(gcf, fig_filename, '-dpdf', '-bestfit');

end
