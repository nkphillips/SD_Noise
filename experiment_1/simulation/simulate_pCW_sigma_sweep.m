%%% simulate_pCW_sigma_sweep
%
% Generates figures illustrating the response bias model (normal CDF +
% guessing lapse) and corresponding normal PDF under a "no encoding bias"
% scenario (mu = 0) while sigma (internal representational noise) increases.
%
% Companion to simulate_pCW.m, which sweeps mu (encoding bias) at fixed sigma.
% Here we instead fix mu = 0 and vary sigma to visualize how increased noise
% flattens the CDF / broadens the PDF, without introducing directional bias.
%
% Figures are saved as vector PDFs to simulation/figures/.

%% Prepare workspace

clear all; close all; clc;

script_dir = pwd;
functions_dir = '../analyses/functions';
plotting_dir  = '../experiment/functions';
addpath(functions_dir);
addpath(plotting_dir);

ps = plotSettings();

%% Output directory

fig_dir = fullfile(script_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Model inputs

guess_rate = 0;       % noise-free CDF (match simulate_pCW.m for consistency)
mu         = 0;       % fixed: no encoding bias

max_offset    = 20;
min_offset    = -max_offset;
probe_offset  = min_offset:0.1:max_offset;

sigma_levels = [2, 5, 10];          % increasing sensory noise (deg)
sigma_labels = arrayfun(@(s) sprintf('\\sigma = %g', s), sigma_levels, 'UniformOutput', false);

%% Colors (low -> high sigma: darker -> lighter gray)

n_sig = numel(sigma_levels);
color_map = zeros(n_sig, 3);
for i = 1:n_sig
    shade = 0.05 + 0.55 * (i-1) / max(n_sig-1, 1);
    color_map(i,:) = [shade shade shade];
end

%% Per-sigma individual figures (matches layout of sfn reference PDFs)

for i = 1:n_sig

    sigma = sigma_levels(i);

    p_cw     = calc_pCW(probe_offset, mu, sigma, guess_rate);
    norm_pdf = normpdf(probe_offset, mu, sigma);
    norm_pdf = norm_pdf / max(norm_pdf);

    fig_name = sprintf('No Encoding Bias (sigma=%g)', sigma);
    fg = figure('Color', 'w', 'Name', fig_name);

    plot(probe_offset, p_cw, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
    hold on;
    plot(probe_offset, norm_pdf, 'LineStyle', '--', 'Color', ps.colors.black, 'LineWidth', ps.line_width);
    xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);

    xlim([min_offset max_offset]);
    ylim([0 1]);
    xticks(min_offset:5:max_offset);
    yticks(0:0.1:1);

    xlabel('\delta\theta (°)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    ylabel('P(RespCW)',        'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title(sprintf('\\mu = 0, \\sigma = %g', sigma), ...
        'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');

    set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
             'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
    box off;
    axis square;

    exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

end

%% Combined overlay figure

fg_all = figure('Color', 'w', 'Name', 'Response Bias Model - sigma sweep (mu=0)');

% CDFs
subplot(1,2,1); hold on;
h_cdf = gobjects(n_sig, 1);
for i = 1:n_sig
    sigma = sigma_levels(i);
    p_cw = calc_pCW(probe_offset, mu, sigma, guess_rate);
    h_cdf(i) = plot(probe_offset, p_cw, 'Color', color_map(i,:), 'LineWidth', ps.line_width);
end
xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
xlim([min_offset max_offset]); ylim([0 1]);
xticks(min_offset:5:max_offset); yticks(0:0.1:1);
xlabel('\delta\theta (°)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('P(RespCW)',        'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title('CDF: \mu = 0, varying \sigma', 'FontSize', ps.axes_label_font_size, ...
    'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_cdf, sigma_labels, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);
set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
         'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
box off; axis square;

% PDFs (normalized to unit peak for visual comparison)
subplot(1,2,2); hold on;
h_pdf = gobjects(n_sig, 1);
for i = 1:n_sig
    sigma = sigma_levels(i);
    norm_pdf = normpdf(probe_offset, mu, sigma);
    norm_pdf = norm_pdf / max(norm_pdf);
    h_pdf(i) = plot(probe_offset, norm_pdf, '--', 'Color', color_map(i,:), 'LineWidth', ps.line_width);
end
xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
xlim([min_offset max_offset]); ylim([0 1]);
xticks(min_offset:5:max_offset); yticks(0:0.1:1);
xlabel('\delta\theta (°)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Normalized density', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title('PDF: \mu = 0, varying \sigma', 'FontSize', ps.axes_label_font_size, ...
    'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_pdf, sigma_labels, 'Location', 'northeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);
set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
         'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
box off; axis square;

exportgraphics(fg_all, fullfile(fig_dir, 'Response Bias Model - sigma sweep (mu=0).pdf'), ...
    'ContentType', 'vector');

%% Done

disp(' ');
disp(['✓ Saved ' num2str(n_sig + 1) ' figures to: ' fig_dir]);
for i = 1:n_sig
    disp(['  - No Encoding Bias (sigma=' num2str(sigma_levels(i)) ').pdf']);
end
disp('  - Response Bias Model - sigma sweep (mu=0).pdf');
