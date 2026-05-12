%%% simulate_DoG
%
% Generates figures illustrating the derivative-of-Gaussian (DoG) serial
% dependence model, plus amplitude and width sweeps.
%
% Figures are saved as vector PDFs to simulation/figures/.

%% Prepare workspace

clear all; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
plotting_dir = fullfile(script_dir, '..', 'experiment', 'functions');
addpath(plotting_dir);

ps = plotSettings();

%% Output directory

fig_dir = fullfile(script_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Define universal parameters

delta_thetas = -90:0.01:90;
c = sqrt(2)/exp(-0.5);
b = 0;

x_lim   = [-90 90];
y_lim   = [-8 8];
x_ticks = -90:45:90;
y_ticks = -8:4:8;

color_map = [
    ps.colors.blue
    ps.colors.black
    ps.colors.red
    ];

%% Baseline DoG

A = 5;
w = 1/20;

y = (delta_thetas * A * w * c) .* exp(-(w * delta_thetas).^2) + b;

fig_name = 'Baseline DoG';
fg = figure('Color', 'w', 'Name', fig_name);

h_dog = plot(delta_thetas, y, 'LineWidth', ps.line_width, 'Color', ps.colors.black);
hold on;
xline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);
yline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);

xlim(x_lim);
ylim(y_lim);
xticks(x_ticks);
yticks(y_ticks);

xlabel('\Delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Bias (deg)',         'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({fig_name, sprintf('A = %g deg, w = 1/%g, b = %g', A, 1/w, b)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_dog, {'DoG'}, 'Location', 'northeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);

set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
box off;
axis square;

exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

%% Annotated DoG width schematic

A = 5;
w = 1/20;

y = (delta_thetas * A * w * c) .* exp(-(w * delta_thetas).^2) + b;
x_peak = 1 / (w * sqrt(2));
y_peak = A + b;
fwhm = 2 * sqrt(log(2)) / w;
half_fwhm = fwhm / 2;

fig_name = 'DoG Width Annotation';
fg = figure('Color', 'w', 'Name', fig_name);

h_dog = plot(delta_thetas, y, 'LineWidth', ps.line_width, 'Color', ps.colors.black);
hold on;
xline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black, 'HandleVisibility', 'off');
yline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black, 'HandleVisibility', 'off');

% w is an inverse-width parameter; annotate intuitive distances in degrees.
scatter(x_peak, y_peak, 45, ps.colors.red, 'filled', 'MarkerEdgeColor', ps.colors.white);
xline(x_peak, ':', 'Color', ps.colors.red, 'LineWidth', ps.line_width, 'HandleVisibility', 'off');
xline(-half_fwhm, ':', 'Color', ps.colors.blue, 'LineWidth', ps.line_width, 'HandleVisibility', 'off');
xline(half_fwhm, ':', 'Color', ps.colors.blue, 'LineWidth', ps.line_width, 'HandleVisibility', 'off');

y_peak_annot = 0.8;
plot([0 x_peak], [y_peak_annot y_peak_annot], '-', ...
    'Color', ps.colors.red, 'LineWidth', 2 * ps.line_width);
text(x_peak / 2, y_peak_annot + 0.45, sprintf('1/(w\\surd2) = %.1f°', x_peak), ...
    'HorizontalAlignment', 'center', 'Color', ps.colors.red, ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, 'Interpreter', 'tex');

y_fwhm_annot = -5.8;
plot([-half_fwhm half_fwhm], [y_fwhm_annot y_fwhm_annot], '-', ...
    'Color', ps.colors.blue, 'LineWidth', 2 * ps.line_width);
plot([-half_fwhm -half_fwhm], y_fwhm_annot + [-0.25 0.25], '-', ...
    'Color', ps.colors.blue, 'LineWidth', 2 * ps.line_width);
plot([half_fwhm half_fwhm], y_fwhm_annot + [-0.25 0.25], '-', ...
    'Color', ps.colors.blue, 'LineWidth', 2 * ps.line_width);
text(0, y_fwhm_annot - 0.65, sprintf('Gaussian-envelope FWHM = %.1f°', fwhm), ...
    'HorizontalAlignment', 'center', 'Color', ps.colors.blue, ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);

text(x_peak + 3, y_peak, 'lobe peak = A', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'Color', ps.colors.red, 'FontName', ps.font_type, ...
    'FontSize', ps.axes_tick_font_size);

xlim(x_lim);
ylim([-5 5]);
xticks(x_ticks);
yticks([-5 0 5]);

xlabel('\Delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Bias (deg)',         'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({fig_name, sprintf('A = %g deg, w = 1/%g, b = %g', A, 1/w, b)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_dog, {'DoG'}, 'Location', 'northeast', 'Box', 'off', ...
    'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);

set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
box off;
axis square;

exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

%% DoG amplitude sweep

A_levels = [2 4 6];
w = 1/20;
y = zeros(numel(A_levels), numel(delta_thetas));

for i = 1:numel(A_levels)
    A = A_levels(i);
    y(i,:) = (delta_thetas * A * w * c) .* exp(-(w * delta_thetas).^2) + b;
end

fig_name = 'DoG Amplitude Sweep';
fg = figure('Color', 'w', 'Name', fig_name);

h_amp = gobjects(numel(A_levels), 1);
hold on;
for i = 1:numel(A_levels)
    h_amp(i) = plot(delta_thetas, y(i,:), 'LineWidth', ps.line_width, 'Color', color_map(i,:));
end
xline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);
yline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);

xlim(x_lim);
ylim(y_lim);
xticks(x_ticks);
yticks(y_ticks);

xlabel('\Delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Bias (deg)',         'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({fig_name, sprintf('w = 1/%g, b = %g', 1/w, b)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_amp, arrayfun(@(A) sprintf('A = %g deg', A), A_levels, 'UniformOutput', false), ...
    'Location', 'northeast', 'Box', 'off', 'FontName', ps.font_type, ...
    'FontSize', ps.axes_tick_font_size);

set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
box off;
axis square;

exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

%% DoG width sweep

A = 5;
w_levels = [1/10 1/20 1/30];
y = zeros(numel(w_levels), numel(delta_thetas));

for i = 1:numel(w_levels)
    w = w_levels(i);
    y(i,:) = (delta_thetas * A * w * c) .* exp(-(w * delta_thetas).^2) + b;
end

fig_name = 'DoG Width Sweep';
fg = figure('Color', 'w', 'Name', fig_name);

h_width = gobjects(numel(w_levels), 1);
hold on;
for i = 1:numel(w_levels)
    h_width(i) = plot(delta_thetas, y(i,:), 'LineWidth', ps.line_width, 'Color', color_map(i,:));
end
xline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);
yline(0, 'LineWidth', ps.line_width, 'Color', ps.colors.black);

xlim(x_lim);
ylim(y_lim);
xticks(x_ticks);
yticks(y_ticks);

xlabel('\Delta\theta (deg)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('Bias (deg)',         'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
title({fig_name, sprintf('A = %g deg, b = %g', A, b)}, ...
    'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type, 'FontWeight', 'normal');
legend(h_width, arrayfun(@(w) sprintf('w = 1/%g', 1/w), w_levels, 'UniformOutput', false), ...
    'Location', 'northeast', 'Box', 'off', 'FontName', ps.font_type, ...
    'FontSize', ps.axes_tick_font_size);

set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length], ...
    'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
box off;
axis square;

exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

%% Done

disp(' ');
disp(['Saved 4 figures to: ' fig_dir]);
disp('  - Baseline DoG.pdf');
disp('  - DoG Width Annotation.pdf');
disp('  - DoG Amplitude Sweep.pdf');
disp('  - DoG Width Sweep.pdf');
