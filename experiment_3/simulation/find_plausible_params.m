% find_plausible_params.m
% Grid search over Weibull parameter spaces for both features (contrast
% and filter width) to find ranges that guarantee all calibrated levels
% fall within the physical stimulus bounds (contrast < 0.9, filter
% width in [0, 80] deg).
%
% Both features use the same Weibull psychometric function. For filter
% width, the independent variable is precision = 1/filter_width.
%
% The output informs the ground truth parameter ranges used in
% simulate_experiment.m and the fitting bounds in fit_calibration.m.

clear all; close all; clc;

addpath('../analyses/plotting');
ps = plotSettings();

%% Physical stimulus bounds (from init_stimuli_params / simulate_calibration_subject)

contrast_max = 0.9;
contrast_min = 0.10;
filter_width_max = 80;
filter_width_min = 10;

%% Fixed psychometric parameters

guess_rate = 0.5; % 2AFC guess rate
lambda = 0.01;    % Lapse rate

%% Target performance levels for calibration inversion

target_levels = [0.65, 0.75, 0.85];

% Weibull inversion constants (same formula for both features)
K = (target_levels - guess_rate) ./ (1 - guess_rate - lambda);
neg_log_1mK = -log(1 - K);

%% ============================================================
%  1. WEIBULL (Contrast) Grid Search
%  ============================================================

n_grid = 300;
alpha_range = linspace(0.01, 0.80, n_grid);
beta_range  = linspace(0.5, 5.0, n_grid);
[A_grid, B_grid] = meshgrid(alpha_range, beta_range);

weibull_valid = true(size(A_grid));
c_levels = zeros([size(A_grid), 3]);

for t = 1:3
    c_levels(:,:,t) = A_grid .* neg_log_1mK(t).^(1./B_grid);
    weibull_valid = weibull_valid & c_levels(:,:,t) >= contrast_min & c_levels(:,:,t) <= contrast_max;
end

weibull_prob = @(x, a, b) guess_rate + (1 - guess_rate - lambda) .* (1 - exp(-(x./a).^b));
P_at_cmin = weibull_prob(contrast_min, A_grid, B_grid);
P_at_cmax = weibull_prob(contrast_max, A_grid, B_grid);

weibull_quality = P_at_cmin > 0.55 & P_at_cmin < 0.90 & P_at_cmax > 0.85;
weibull_feasible = weibull_valid & weibull_quality;

[r, c_idx] = find(weibull_feasible);
alpha_feasible_min = alpha_range(min(c_idx));
alpha_feasible_max = alpha_range(max(c_idx));
beta_feasible_min  = beta_range(min(r));
beta_feasible_max  = beta_range(max(r));

fprintf('\n=== WEIBULL (Contrast) Feasible Region ===\n');
fprintf('alpha: [%.3f, %.3f]\n', alpha_feasible_min, alpha_feasible_max);
fprintf('beta:  [%.2f, %.2f]\n', beta_feasible_min, beta_feasible_max);

%% ============================================================
%  2. WEIBULL on PRECISION (Filter Width) Grid Search
%  ============================================================
% Same Weibull function as contrast, but the independent variable is
% precision = 1/filter_width. This transforms the decreasing
% performance-vs-fw relationship into a standard increasing one.
%
% Inversion: precision_target = alpha_fw * (-log(1-K))^(1/beta_fw)
%            fw_target = 1 / precision_target
%
% Feasibility: all 3 inverted fw must be in (0, filter_width_max].

% Precision range corresponding to fw = [10, 80] deg
prec_min = 1 / filter_width_max; % 0.0125
prec_max = 1 / filter_width_min; % 0.1

alpha_fw_range = linspace(0.005, 0.15, n_grid);
beta_fw_range  = linspace(0.5, 5.0, n_grid);
[A_fw_grid, B_fw_grid] = meshgrid(alpha_fw_range, beta_fw_range);

fw_valid = true(size(A_fw_grid));
fw_levels = zeros([size(A_fw_grid), 3]);

for t = 1:3
    % Invert the Weibull in the precision domain
    prec_level = A_fw_grid .* neg_log_1mK(t).^(1./B_fw_grid);
    % Convert precision back to filter width
    fw_level = 1 ./ prec_level;
    fw_levels(:,:,t) = fw_level;
    fw_valid = fw_valid & fw_level >= filter_width_min & fw_level <= filter_width_max;
end

% Quality: psychometric function spans a meaningful range across the
% tested precision values [1/80, 1/10]
weibull_prob_fw = @(x, a, b) guess_rate + (1 - guess_rate - lambda) .* (1 - exp(-(x./a).^b));
P_at_prec_max = weibull_prob_fw(prec_max, A_fw_grid, B_fw_grid); % narrowest filter (easiest)
P_at_prec_min = weibull_prob_fw(prec_min, A_fw_grid, B_fw_grid); % widest filter (hardest)

fw_quality = P_at_prec_max > 0.75 & P_at_prec_min > 0.52 & P_at_prec_min < 0.80;
fw_feasible = fw_valid & fw_quality;

[r_fw, c_fw] = find(fw_feasible);
alpha_fw_feasible_min = alpha_fw_range(min(c_fw));
alpha_fw_feasible_max = alpha_fw_range(max(c_fw));
beta_fw_feasible_min  = beta_fw_range(min(r_fw));
beta_fw_feasible_max  = beta_fw_range(max(r_fw));

fprintf('\n=== WEIBULL on PRECISION (Filter Width) Feasible Region ===\n');
fprintf('alpha_fw (precision threshold): [%.4f, %.4f]\n', alpha_fw_feasible_min, alpha_fw_feasible_max);
fprintf('beta_fw  (slope):               [%.2f, %.2f]\n', beta_fw_feasible_min, beta_fw_feasible_max);

%% ============================================================
%  3. RECOMMENDED RANGES for simulate_experiment.m
%  ============================================================

% Contrast Weibull
rec_alpha_mean = (alpha_feasible_min + alpha_feasible_max) / 2;
rec_alpha_std  = (alpha_feasible_max - alpha_feasible_min) / 6;
rec_beta_mean  = (beta_feasible_min + beta_feasible_max) / 2;
rec_beta_std   = (beta_feasible_max - beta_feasible_min) / 6;

% Filter Width Weibull (precision domain)
rec_alpha_fw_mean = (alpha_fw_feasible_min + alpha_fw_feasible_max) / 2;
rec_alpha_fw_std  = (alpha_fw_feasible_max - alpha_fw_feasible_min) / 6;
rec_beta_fw_mean  = (beta_fw_feasible_min + beta_fw_feasible_max) / 2;
rec_beta_fw_std   = (beta_fw_feasible_max - beta_fw_feasible_min) / 6;

fprintf('\n=== RECOMMENDED simulate_experiment.m PARAMETERS ===\n');
fprintf('\n--- Contrast Weibull ---\n');
fprintf('gt.contrast_alpha = %.2f + %.2f * randn();\n', rec_alpha_mean, rec_alpha_std);
fprintf('gt.contrast_alpha = max(%.3f, min(%.3f, gt.contrast_alpha));\n', alpha_feasible_min, alpha_feasible_max);
fprintf('gt.contrast_beta  = %.2f + %.2f * randn();\n', rec_beta_mean, rec_beta_std);
fprintf('gt.contrast_beta  = max(%.2f, min(%.2f, gt.contrast_beta));\n', beta_feasible_min, beta_feasible_max);

fprintf('\n--- Filter Width Weibull (precision = 1/fw) ---\n');
fprintf('gt.filter_alpha = %.4f + %.4f * randn();\n', rec_alpha_fw_mean, rec_alpha_fw_std);
fprintf('gt.filter_alpha = max(%.4f, min(%.4f, gt.filter_alpha));\n', alpha_fw_feasible_min, alpha_fw_feasible_max);
fprintf('gt.filter_beta  = %.2f + %.2f * randn();\n', rec_beta_fw_mean, rec_beta_fw_std);
fprintf('gt.filter_beta  = max(%.2f, min(%.2f, gt.filter_beta));\n', beta_fw_feasible_min, beta_fw_feasible_max);

fprintf('\n--- fit_calibration.m bounds ---\n');
fprintf('lb_fw = [%.4f, %.2f];  ub_fw = [%.4f, %.2f];\n', ...
    alpha_fw_feasible_min * 0.5, max(0.5, beta_fw_feasible_min * 0.5), ...
    alpha_fw_feasible_max * 1.5, beta_fw_feasible_max * 1.5);

%% ============================================================
%  4. VISUALIZATIONS
%  ============================================================

fig = figure('Color', 'w', 'Position', [50, 50, 1400, 600], 'Name', 'Plausible Parameter Search');

% --- Contrast Weibull feasibility map ---
subplot(1,2,1);
imagesc(alpha_range, beta_range, double(weibull_feasible));
set(gca, 'YDir', 'normal');
colormap(gca, [1 1 1; ps.colors.blue]);
hold on;

alpha_boundary = contrast_max ./ neg_log_1mK(3).^(1./beta_range);
plot(alpha_boundary, beta_range, '-', 'Color', ps.colors.red, 'LineWidth', 2);

rectangle('Position', [alpha_feasible_min, beta_feasible_min, ...
    alpha_feasible_max - alpha_feasible_min, beta_feasible_max - beta_feasible_min], ...
    'EdgeColor', ps.colors.green, 'LineWidth', 2, 'LineStyle', '--');

xlabel('\alpha (Contrast Threshold)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('\beta (Slope)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title('Contrast Weibull: Feasible Region', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
legend({'c_{85} = 0.9 boundary', 'Recommended'}, 'Location', 'northeast', 'FontSize', ps.axes_tick_font_size);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

% --- Filter Width Weibull (precision) feasibility map ---
subplot(1,2,2);
imagesc(alpha_fw_range, beta_fw_range, double(fw_feasible));
set(gca, 'YDir', 'normal');
colormap(gca, [1 1 1; ps.colors.blue]);
hold on;

% Boundary where the hardest target (85%) reaches fw_max (precision_min)
% precision_target = alpha * neg_log_1mK(3)^(1/beta)
% Boundary: alpha = prec_min / neg_log_1mK(3)^(1/beta)
alpha_fw_boundary = prec_min ./ neg_log_1mK(3).^(1./beta_fw_range);
plot(alpha_fw_boundary, beta_fw_range, '-', 'Color', ps.colors.red, 'LineWidth', 2);

rectangle('Position', [alpha_fw_feasible_min, beta_fw_feasible_min, ...
    alpha_fw_feasible_max - alpha_fw_feasible_min, beta_fw_feasible_max - beta_fw_feasible_min], ...
    'EdgeColor', ps.colors.green, 'LineWidth', 2, 'LineStyle', '--');

xlabel('\alpha_{fw} (Precision Threshold)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('\beta_{fw} (Slope)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title('Filter Width Weibull (1/fw): Feasible', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
legend({'fw_{85} = 80 boundary', 'Recommended'}, 'Location', 'northeast', 'FontSize', ps.axes_tick_font_size);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

%% Save figure
fig_dir = fullfile(fileparts(mfilename('fullpath')), 'experiment', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print(fig, fullfile(fig_dir, 'plausible_parameter_search.pdf'), '-dpdf', '-bestfit');
disp(['Figure saved to ' fullfile(fig_dir, 'plausible_parameter_search.pdf')]);
