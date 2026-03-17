% find_plausible_params.m
% Grid search over Weibull and PTM (Perceptual Template Model) parameter
% spaces to find ranges that guarantee all calibrated levels fall within
% the physical stimulus bounds (contrast < 0.9, filter width < 80 deg,
% all > 0).
%
% PTM model (Lu & Dosher, 2008, Eq. 13; gamma fixed at 2.0):
%   d' = Signal^g / sqrt((1+Nm^2)*fw^(2g) + Nm^2*Signal^(2g) + Na^2)
%   P  = (1-lambda)*normcdf(d'/sqrt(2)) + lambda*0.5
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
ptm_gamma = 2.0;  % Fixed transducer exponent (Lu & Dosher, 1999; Ramirez et al., 2021 n~1.78)
g = ptm_gamma;

%% Target performance levels for calibration inversion

target_levels = [0.65, 0.75, 0.85];

% Weibull inversion constants
K = (target_levels - guess_rate) ./ (1 - guess_rate - lambda);
neg_log_1mK = -log(1 - K);

% PTM inversion constants
P_target_adj = (target_levels - lambda * 0.5) ./ (1 - lambda);
d_target = sqrt(2) .* norminv(P_target_adj); % [0.551, 0.965, 1.485]

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
    weibull_valid = weibull_valid & c_levels(:,:,t) > 0 & c_levels(:,:,t) < contrast_max;
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
%  2. PTM (Filter Width) Grid Search
%  ============================================================
% 3D search over (Signal, N_mul, N_add) with gamma fixed.
% Inversion: fw^(2g) = (S^(2g)*(1/d'^2 - Nm^2) - Na^2) / (1+Nm^2)

n_s = 120;  n_m = 40;  n_a = 120;
signal_range = linspace(5, 80, n_s);
nmul_range   = linspace(0.01, 0.65, n_m);
nadd_range   = linspace(1, 600, n_a);

ptm_feasible_3d = false(n_m, n_s, n_a); % (nmul, signal, nadd)

ptm_dprime = @(fw, S, Nm, Na) S.^g ./ sqrt((1+Nm.^2).*fw.^(2*g) + Nm.^2.*S.^(2*g) + Na.^2);
ptm_prob   = @(fw, S, Nm, Na) (1-lambda).*normcdf(ptm_dprime(fw,S,Nm,Na)./sqrt(2)) + lambda*0.5;

for im = 1:n_m
    Nm = nmul_range(im);
    for is = 1:n_s
        S = signal_range(is);
        S2g = S^(2*g);
        for ia = 1:n_a
            Na = nadd_range(ia);

            valid = true;
            for t = 1:3
                denom = 1 + Nm^2;
                numer = S2g * (1/d_target(t)^2 - Nm^2) - Na^2;
                if numer <= 0 || (1/d_target(t)^2 <= Nm^2)
                    valid = false; break;
                end
                fw_val = (numer / denom)^(1/(2*g));
                if fw_val <= 0 || fw_val >= filter_width_max
                    valid = false; break;
                end
            end

            if valid
                P_fwmin = ptm_prob(filter_width_min, S, Nm, Na);
                P_fwmax = ptm_prob(filter_width_max, S, Nm, Na);
                if P_fwmin > 0.75 && P_fwmax > 0.55 && P_fwmax < 0.80
                    ptm_feasible_3d(im, is, ia) = true;
                end
            end
        end
    end
end

ptm_any_feasible = squeeze(any(ptm_feasible_3d, 1)); % (signal, nadd) -- any N_mul works

% Extract marginal feasible ranges
[fs, fa, fm] = ind2sub(size(ptm_feasible_3d), find(ptm_feasible_3d));
signal_feasible_min = signal_range(min(fa));
signal_feasible_max = signal_range(max(fa));
nmul_feasible_min   = nmul_range(min(fs));
nmul_feasible_max   = nmul_range(max(fs));
nadd_feasible_min   = nadd_range(min(fm));
nadd_feasible_max   = nadd_range(max(fm));

fprintf('\n=== PTM (Filter Width) Feasible Region (gamma = %.1f) ===\n', ptm_gamma);
fprintf('Signal: [%.1f, %.1f]\n', signal_feasible_min, signal_feasible_max);
fprintf('N_mul:  [%.3f, %.3f]\n', nmul_feasible_min, nmul_feasible_max);
fprintf('N_add:  [%.1f, %.1f]\n', nadd_feasible_min, nadd_feasible_max);

% Project feasible region onto Signal x N_add plane at a representative N_mul
mid_nmul_idx = round(n_m / 2);
ptm_slice = squeeze(ptm_feasible_3d(mid_nmul_idx, :, :)); % (signal, nadd)

%% ============================================================
%  3. RECOMMENDED RANGES for simulate_experiment.m
%  ============================================================

% Weibull
rec_alpha_mean = (alpha_feasible_min + alpha_feasible_max) / 2;
rec_alpha_std  = (alpha_feasible_max - alpha_feasible_min) / 6;
rec_beta_mean  = (beta_feasible_min + beta_feasible_max) / 2;
rec_beta_std   = (beta_feasible_max - beta_feasible_min) / 6;

% PTM
rec_signal_mean = (signal_feasible_min + signal_feasible_max) / 2;
rec_signal_std  = (signal_feasible_max - signal_feasible_min) / 6;
rec_nmul_mean   = (nmul_feasible_min + nmul_feasible_max) / 2;
rec_nmul_std    = (nmul_feasible_max - nmul_feasible_min) / 6;
rec_nadd_mean   = (nadd_feasible_min + nadd_feasible_max) / 2;
rec_nadd_std    = (nadd_feasible_max - nadd_feasible_min) / 6;

fprintf('\n=== RECOMMENDED simulate_experiment.m PARAMETERS ===\n');
fprintf('\n--- Weibull ---\n');
fprintf('gt.contrast_alpha = %.2f + %.2f * randn();\n', rec_alpha_mean, rec_alpha_std);
fprintf('gt.contrast_alpha = max(%.3f, min(%.3f, gt.contrast_alpha));\n', alpha_feasible_min, alpha_feasible_max);
fprintf('gt.contrast_beta  = %.2f + %.2f * randn();\n', rec_beta_mean, rec_beta_std);
fprintf('gt.contrast_beta  = max(%.2f, min(%.2f, gt.contrast_beta));\n', beta_feasible_min, beta_feasible_max);

fprintf('\n--- PTM (gamma = %.1f fixed) ---\n', ptm_gamma);
fprintf('gt.filter_signal = %.1f + %.1f * randn();\n', rec_signal_mean, rec_signal_std);
fprintf('gt.filter_signal = max(%.1f, min(%.1f, gt.filter_signal));\n', signal_feasible_min, signal_feasible_max);
fprintf('gt.filter_nmul   = %.3f + %.3f * randn();\n', rec_nmul_mean, rec_nmul_std);
fprintf('gt.filter_nmul   = max(%.3f, min(%.3f, gt.filter_nmul));\n', nmul_feasible_min, nmul_feasible_max);
fprintf('gt.filter_nadd   = %.1f + %.1f * randn();\n', rec_nadd_mean, rec_nadd_std);
fprintf('gt.filter_nadd   = max(%.1f, min(%.1f, gt.filter_nadd));\n', nadd_feasible_min, nadd_feasible_max);

fprintf('\n--- fit_calibration.m bounds ---\n');
fprintf('lb_fw = [%.1f, %.3f, %.1f];  ub_fw = [%.1f, %.3f, %.1f];\n', ...
    signal_feasible_min * 0.5, max(0.001, nmul_feasible_min * 0.5), nadd_feasible_min * 0.5, ...
    signal_feasible_max * 1.5, nmul_feasible_max * 1.5, nadd_feasible_max * 1.5);

%% ============================================================
%  4. VISUALIZATIONS
%  ============================================================

fig = figure('Color', 'w', 'Position', [50, 50, 1400, 600], 'Name', 'Plausible Parameter Search');

% --- Weibull feasibility map ---
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
title('Weibull: Feasible Region', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
legend({'c_{85} = 0.9 boundary', 'Recommended'}, 'Location', 'northeast', 'FontSize', ps.axes_tick_font_size);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

% --- PTM feasibility map (Signal x N_add slice at median N_mul) ---
subplot(1,2,2);
imagesc(signal_range, nadd_range, double(ptm_slice'));
set(gca, 'YDir', 'normal');
colormap(gca, [1 1 1; ps.colors.blue]);
hold on;

rectangle('Position', [signal_feasible_min, nadd_feasible_min, ...
    signal_feasible_max - signal_feasible_min, nadd_feasible_max - nadd_feasible_min], ...
    'EdgeColor', ps.colors.green, 'LineWidth', 2, 'LineStyle', '--');

xlabel('Signal', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('N_{add}', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title(sprintf('PTM: Feasible (\\gamma=%.0f, N_{mul}=%.2f slice)', ptm_gamma, nmul_range(mid_nmul_idx)), ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
legend({'Recommended'}, 'Location', 'northeast', 'FontSize', ps.axes_tick_font_size);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

%% Save figure
fig_dir = 'experiment/figures';
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
exportgraphics(fig, fullfile(fig_dir, 'plausible_parameter_search.pdf'), 'ContentType', 'vector');
disp(['Figure saved to ' fullfile(fig_dir, 'plausible_parameter_search.pdf')]);
