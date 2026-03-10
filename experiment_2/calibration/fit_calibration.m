% fit_calibration.m
% Extracts 3 subject-specific levels from calibration data.

clear all; close all; clc;

%% 1. Set subject ID and directories
subj_ID = '999';
data_dir = ['../../data/' subj_ID];

% Load all calibration runs for this subject
file_pattern = [data_dir '/SD_Noise_Exp2_Calibration_S' subj_ID '_Run*.mat'];
files = dir(file_pattern);

if isempty(files)
    error('No calibration data found for subject %s.', subj_ID);
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
        
        % p.trial_events structure: [test_orientations, probe_orientations, level_order]
        levels = p.trial_events(:, 3, n_block);
        correct = behav_data.correct(:, n_block);
        
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

init_guess_c = [mean(unique_c), 2]; % initial alpha and beta
best_params_c = fminsearch(nll_weibull, init_guess_c, options);
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

init_guess_fw = [15, 10]; % initial Signal and internal noise
best_params_fw = fminsearch(nll_en, init_guess_fw, options);
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
plot(unique_c, prop_c, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
x_fit_c = linspace(0, max(unique_c), 100);
plot(x_fit_c, weibull_prob(x_fit_c, alpha_c, beta_c), 'r-', 'LineWidth', 2);
for i = 1:3
    plot([0 calib_contrast(i)], [target_levels(i) target_levels(i)], 'b--', 'HandleVisibility', 'off');
    plot([calib_contrast(i) calib_contrast(i)], [0 target_levels(i)], 'b--', 'HandleVisibility', 'off');
    plot(calib_contrast(i), target_levels(i), 'b*', 'MarkerSize', 8);
end
xlim([0 max(max(unique_c), max(calib_contrast))*1.1]);
ylim([0.4 1.0]);
xlabel('Contrast'); ylabel('Proportion Correct');
title(sprintf('Weibull Fit (\\alpha=%.2f, \\beta=%.2f)', alpha_c, beta_c));
box off; axis square;

% Plot Filter Width
subplot(1,2,2);
plot(unique_fw, prop_fw, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); hold on;
x_fit_fw = linspace(0, max(unique_fw), 100);
plot(x_fit_fw, en_prob(x_fit_fw, signal_fw, sigma_int_fw), 'r-', 'LineWidth', 2);
for i = 1:3
    plot([0 calib_filter(i)], [target_levels(i) target_levels(i)], 'b--', 'HandleVisibility', 'off');
    plot([calib_filter(i) calib_filter(i)], [0 target_levels(i)], 'b--', 'HandleVisibility', 'off');
    plot(calib_filter(i), target_levels(i), 'b*', 'MarkerSize', 8);
end
xlim([0 max(max(unique_fw), max(calib_filter))*1.1]);
ylim([0.4 1.0]);
xlabel('Filter Width (\sigma_{ext})'); ylabel('Proportion Correct');
title(sprintf('Eq Noise Fit (Sig=%.2f, \\sigma_{int}=%.2f)', signal_fw, sigma_int_fw));
box off; axis square;
