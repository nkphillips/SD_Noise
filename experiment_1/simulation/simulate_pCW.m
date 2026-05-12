%%% simulate_pCW
%
% Generates figures illustrating how an encoding bias shifts the response
% bias model (normal CDF + guessing lapse) and corresponding normal PDF.
%
% Figures are saved as vector PDFs to simulation/figures/.

%% Prepare workspace

clear all; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
functions_dir = fullfile(script_dir, '..', 'analyses', 'functions');
plotting_dir  = fullfile(script_dir, '..', 'experiment', 'functions');
addpath(functions_dir);
addpath(plotting_dir);

ps = plotSettings();

%% Output directory

fig_dir = fullfile(script_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Model inputs

guess_rate = 0;
sigma      = 5;       % fixed representational noise (deg)

max_offset   = 20;
min_offset   = -max_offset;
probe_offset = min_offset:0.1:max_offset;

x_lim   = [min_offset max_offset];
y_lim   = [0 1];
x_ticks = min_offset:5:max_offset;
y_ticks = 0:0.1:1;

bias_levels = [-5, 0, 5];
bias_names  = {'Counter Clockwise Encoding Bias', 'No Encoding Bias', 'Clockwise Encoding Bias'};
bias_labels = arrayfun(@(m) sprintf('\\mu = %g', m), bias_levels, 'UniformOutput', false);

color_map = [
    ps.colors.blue
    ps.colors.black
    ps.colors.red
];

%% Per-bias individual figures

for i = 1:numel(bias_levels)

    mu = bias_levels(i);

    p_cw     = calc_pCW(probe_offset, mu, sigma, guess_rate);
    norm_pdf = normpdf(probe_offset, mu, sigma);
    norm_pdf = norm_pdf / max(norm_pdf);

    fig_name = bias_names{i};
    fg = figure('Color', 'w', 'Name', fig_name);

    h_cdf = plot(probe_offset, p_cw, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
    hold on;
    h_pdf = plot(probe_offset, norm_pdf, 'LineStyle', '--', 'Color', ps.colors.black, 'LineWidth', ps.line_width);
    xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);

    xlim(x_lim);
    ylim(y_lim);
    xticks(x_ticks);
    yticks(y_ticks);

    xlabel('\delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    ylabel('P(RespCW)',          'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    title({fig_name, sprintf('\\mu = %g, \\sigma = %g', mu, sigma)}, ...
        'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
    legend([h_cdf h_pdf], {'P(RespCW)', 'Normalized PDF'}, 'Location', 'northwest', 'Box', 'off', ...
        'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);

    set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
             'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
    box off;
    axis square;

    exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

end

%% Combined overlay figure

fg_all = figure('Color', 'w', 'Name', sprintf('Response Bias Model - encoding bias sweep (sigma=%g)', sigma));

% CDFs
subplot(1,2,1); hold on;
h_cdf = gobjects(numel(bias_levels), 1);
for i = 1:numel(bias_levels)
    mu = bias_levels(i);
    p_cw = calc_pCW(probe_offset, mu, sigma, guess_rate);
    h_cdf(i) = plot(probe_offset, p_cw, 'Color', color_map(i,:), 'LineWidth', ps.line_width);
end
xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
xlim(x_lim); ylim(y_lim);
xticks(x_ticks); yticks(y_ticks);
xlabel('\delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('P(RespCW)',          'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({'CDF: varying \mu', sprintf('\\sigma = %g', sigma)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_cdf, bias_labels, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);
set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
         'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
box off; axis square;

% PDFs (normalized to unit peak for visual comparison)
subplot(1,2,2); hold on;
h_pdf = gobjects(numel(bias_levels), 1);
for i = 1:numel(bias_levels)
    mu = bias_levels(i);
    norm_pdf = normpdf(probe_offset, mu, sigma);
    norm_pdf = norm_pdf / max(norm_pdf);
    h_pdf(i) = plot(probe_offset, norm_pdf, '--', 'Color', color_map(i,:), 'LineWidth', ps.line_width);
end
xline(0, 'Color', ps.colors.black, 'LineWidth', ps.line_width);
xlim(x_lim); ylim(y_lim);
xticks(x_ticks); yticks(y_ticks);
xlabel('\delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Normalized density', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({'PDF: varying \mu', sprintf('\\sigma = %g', sigma)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_pdf, bias_labels, 'Location', 'northeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);
set(gca, 'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
         'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length]);
box off; axis square;

exportgraphics(fg_all, fullfile(fig_dir, sprintf('Response Bias Model - encoding bias sweep (sigma=%g).pdf', sigma)), ...
    'ContentType', 'vector');

%% Done

disp(' ');
disp(['Saved ' num2str(numel(bias_levels) + 1) ' figures to: ' fig_dir]);
for i = 1:numel(bias_names)
    disp(['  - ' bias_names{i} '.pdf']);
end
disp(['  - Response Bias Model - encoding bias sweep (sigma=' num2str(sigma) ').pdf']);
