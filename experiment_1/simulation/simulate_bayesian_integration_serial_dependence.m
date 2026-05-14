%%% simulate_bayesian_integration_serial_dependence
%
% Generates a two-panel schematic of Bayesian integration for serial
% dependence. The prior is centered on a previous stimulus, the likelihood is
% centered on the current stimulus, and the posterior mode gives the
% attractive bias toward the prior.
%
% Figure is saved as a vector PDF to simulation/figures/.

%% Prepare workspace

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
plotting_dir = fullfile(script_dir, '..', 'experiment', 'functions');
addpath(plotting_dir);

ps = plotSettings();

%% Output directory

fig_dir = fullfile(script_dir, 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Define orientation domain

x = linspace(-90, 90, 1000);

%% Define Bayesian integration parameters

prior_mu = 40;
prior_sigma = 15;
likelihood_mu = 0;

conditions = struct( ...
    'title', {'Baseline', 'High Sensory Uncertainty'}, ...
    'likelihood_sigma', {10, 30});

prior = normpdf(x, prior_mu, prior_sigma);

%% Plot figure

fig_name = 'Bayesian Integration Serial Dependence';
fg = figure('Color', 'w', 'Name', fig_name, 'Position', [100 100 1000 400]);

prior_color = 0.45 * [1 1 1];
likelihood_color = ps.colors.blue;
posterior_color = ps.colors.black;
line_width = 2;

posterior_modes = nan(1, numel(conditions));
posterior_peak_ys = nan(1, numel(conditions));
all_y = prior;

for i_condition = 1:numel(conditions)
    likelihood_sigma = conditions(i_condition).likelihood_sigma;
    likelihood = normpdf(x, likelihood_mu, likelihood_sigma);
    posterior = prior .* likelihood;
    posterior = posterior ./ trapz(x, posterior);

    posterior_mode = exactGaussianProductMode( ...
        prior_mu, prior_sigma, likelihood_mu, likelihood_sigma);
    posterior_peak_y = interp1(x, posterior, posterior_mode, 'pchip');

    posterior_modes(i_condition) = posterior_mode;
    posterior_peak_ys(i_condition) = posterior_peak_y;
    all_y = [all_y likelihood posterior]; %#ok<AGROW>

    ax = subplot(1, 2, i_condition);
    hold(ax, 'on');

    h_prior = plot(ax, x, prior, '--', ...
        'Color', prior_color, 'LineWidth', line_width);
    h_likelihood = plot(ax, x, likelihood, '-', ...
        'Color', likelihood_color, 'LineWidth', line_width);
    h_posterior = plot(ax, x, posterior, '-', ...
        'Color', posterior_color, 'LineWidth', line_width);

    plot(ax, [posterior_mode posterior_mode], [0 posterior_peak_y], ':', ...
        'Color', posterior_color, 'LineWidth', line_width, ...
        'HandleVisibility', 'off');

    attractive_bias = posterior_mode - likelihood_mu;
    text(ax, posterior_mode + 4, 0.82 * posterior_peak_y, ...
        sprintf('Attractive Bias = %.1f deg', attractive_bias), ...
        'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size, ...
        'Color', posterior_color, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');

    xlim(ax, [-90 90]);
    xticks(ax, -90:45:90);
    xlabel(ax, 'Orientation ($^\circ$)', ...
        'Interpreter', 'latex', 'FontSize', ps.axes_label_font_size);
    ylabel(ax, 'Probability Density', ...
        'FontName', ps.font_type, 'FontSize', ps.axes_label_font_size);
    title(ax, conditions(i_condition).title, ...
        'FontName', ps.font_type, 'FontSize', ps.axes_label_font_size, ...
        'FontWeight', 'normal');

    set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length ps.tick_length], ...
        'FontName', ps.font_type, 'FontSize', 14);
    box(ax, 'off');

    if i_condition == 1
        legend(ax, [h_prior h_likelihood h_posterior], ...
            {'Prior', 'Likelihood', 'Posterior'}, ...
            'Location', 'northeast', 'Box', 'off', ...
            'FontName', ps.font_type, 'FontSize', ps.axes_tick_font_size);
    end
end

y_lim = [0 1.12 * max(all_y)];
for i_condition = 1:numel(conditions)
    ylim(subplot(1, 2, i_condition), y_lim);
end

exportgraphics(fg, fullfile(fig_dir, [fig_name '.pdf']), 'ContentType', 'vector');

fprintf('Baseline posterior mode: %.2f deg\n', posterior_modes(1));
fprintf('High uncertainty posterior mode: %.2f deg\n', posterior_modes(2));
fprintf('Saved figure: %s\n', fullfile(fig_dir, [fig_name '.pdf']));

%% Local functions

function posterior_mode = exactGaussianProductMode(mu_1, sigma_1, mu_2, sigma_2)
% Mode of the product of two Gaussian densities.
precision_1 = 1 / sigma_1^2;
precision_2 = 1 / sigma_2^2;
posterior_mode = (mu_1 * precision_1 + mu_2 * precision_2) / ...
    (precision_1 + precision_2);
end
