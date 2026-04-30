%%% rerender_super_bias_arrows
% Post-hoc rerender of the super-subject response-bias grids so subplot
% titles use the same \rightarrow arrow style as the poster Fig6 panels.
%
% This uses saved 04.24.2026 estimates and does not rerun model fitting.

clear; clc; close all;

addpath('functions');
addpath('plotting');

analysis_date = '04.24.2026';
n_back_list = [1 2 3];
estimates = loadEstimates(analysis_date, n_back_list);

p = struct();
p.cond_names = {'Contrast','Precision'};
p.contrast = {'90%','50%','25%'};
p.precision = {'2°','40°','80°'};
p.rb_bounds = [20, 90; -20, 0.1];

plt_opts = posterPlotOpts();
plt_opts.save_sup_figures = 1;
plt_opts.fg_type = 'pdf';
plt_opts.marker_size_scatter = 50;

for n_back = n_back_list
    key = sprintf('n%d', n_back);
    est = estimates.(key);
    if isempty(est)
        warning('rerenderSuperBiasArrows:noData', ...
            'Skipping n_back=%d (no estimates).', n_back);
        continue
    end

    delta_theta_centers = est.delta_theta_centers;
    if isempty(delta_theta_centers)
        delta_theta_centers = -90:1:90;
    end

    rb = est.rb;
    sd = est.sd;
    rb_ci = est.rb_ci;

    plt_opts.sup_figure_path = fullfile('figures', 'super', ...
        analysis_date, sprintf('%d_back', n_back));
    if ~exist(plt_opts.sup_figure_path, 'dir')
        mkdir(plt_opts.sup_figure_path);
    end

    disp(sprintf('Rerendering super-subject bias arrows for %d-back...', n_back)); %#ok<DSPS>

    renderBiasGrid(rb, sd, rb_ci, delta_theta_centers, p, plt_opts, false);
    renderBiasGrid(rb, sd, rb_ci, delta_theta_centers, p, plt_opts, true);
end

disp('Done.');

function renderBiasGrid(rb, sd, rb_ci, delta_theta_centers, p, plt_opts, with_dog)
    num_levels = 3;
    num_conds = 2;

    for cond = 1:num_conds
        if with_dog
            fg_name = ['Super Subj Bias with DoG ' p.cond_names{cond}];
        else
            fg_name = ['Super Subj Bias ' p.cond_names{cond}];
        end

        fg = figure('Visible','off','Color','w');

        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);
                hold on;

                if cond == 1
                    fg_title = sprintf('%s \\rightarrow %s', ...
                        p.contrast{prev_lvl}, p.contrast{curr_lvl});
                else
                    fg_title = sprintf('%s \\rightarrow %s', ...
                        p.precision{prev_lvl}, p.precision{curr_lvl});
                end

                sd_params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
                baseline = [];
                if ~isempty(sd_params) && ~any(isnan(sd_params))
                    baseline = sd_params(3);
                end

                mu = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1));
                if isfield(rb_ci, 'mu_lo')
                    mu_lo = squeeze(rb_ci.mu_lo(prev_lvl, curr_lvl, cond, :));
                    mu_hi = squeeze(rb_ci.mu_hi(prev_lvl, curr_lvl, cond, :));
                    plotResponseBias(delta_theta_centers, mu, plt_opts, cond, mu_lo, mu_hi, baseline);
                else
                    plotResponseBias(delta_theta_centers, mu, plt_opts, cond, [], [], baseline);
                end

                if with_dog && ~isempty(sd_params)
                    dog_params = sd_params(1:3);
                    delta_smooth = linspace(-90, 90, 100);
                    dog_fit = calcDoG(delta_smooth, dog_params);
                    if isfield(plt_opts, 'rb_subtract_baseline') && plt_opts.rb_subtract_baseline
                        dog_fit = dog_fit - dog_params(3);
                    end
                    plot(delta_smooth, dog_fit, '--', ...
                        'LineWidth', plt_opts.line_width, ...
                        'Color', plt_opts.colors.black);

                    amp_est = dog_params(1);
                    w_est = max(dog_params(2), eps);
                    fwhm_est = (2 * sqrt(log(2))) / w_est;
                    text(0.95, 0.10, sprintf('A = %.2f°\nFWHM = %.1f°\nb = %.2f°', ...
                        amp_est, fwhm_est, dog_params(3)), ...
                        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                        'VerticalAlignment', 'bottom', 'FontSize', 8, ...
                        'Color', plt_opts.colors.black);

                    curr_r2 = sd.all.r2(prev_lvl, curr_lvl, cond);
                    if ~isnan(curr_r2)
                        text(0.05, 0.90, sprintf('R^2 = %.2f', curr_r2), ...
                            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'top', 'FontSize', 8, ...
                            'Color', plt_opts.colors.black);
                    end
                end

                axis square;
                title(fg_title);
                if prev_lvl == num_levels
                    xlabel('\Delta\theta (°)');
                    xticks(-90:45:90);
                    xtickangle(0);
                else
                    xticks(-90:45:90);
                    xtickangle(0);
                    set(gca, 'XTickLabel', []);
                end
                if curr_lvl == 1
                    ylabel('Bias (°)');
                end
                set(gca, 'TickDir', 'out', ...
                    'TickLength', [plt_opts.tick_length, plt_opts.tick_length]);
                xlim([-90 90]);
                line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');
                line([0, 0], [p.rb_bounds(2,1), p.rb_bounds(1,1)], ...
                    'LineWidth', 1, 'Color', 'k');
                ylim([p.rb_bounds(2,1) p.rb_bounds(1,1)]);
                yticks(p.rb_bounds(2,1):5:p.rb_bounds(1,1));
                box off;
            end
        end

        saveas(fg, fullfile(plt_opts.sup_figure_path, [fg_name '.' plt_opts.fg_type]));
        close(fg);
    end
end
