function plotBiasDogGrid(analysis_date, n_back_list, opts)
% plotBiasDogGrid
%
% Poster Figure 6 - Full previous-by-current bias + DoG panels. For each
% n-back and condition, renders a 3x3 grid where rows are previous level,
% columns are current level, and each cell shows:
%   - mu(\Delta\theta) with rb_ci.mu_* bands
%   - DoG overlay from sd.all.params_est
%   - R^2, A, FWHM, b annotations
%
% Usage:
%   plotBiasDogGrid('04.24.2026', [1 2 3]);

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'objective'), opts.objective = 'sse'; end
    if ~isfield(opts, 'save'),      opts.save      = true;  end

    plt_opts = posterPlotOpts();
    p        = defaultLabels();

    this_dir = fileparts(mfilename('fullpath'));
    analyses_dir = fullfile(this_dir, '..');
    if ~isfield(opts, 'fig_dir') || isempty(opts.fig_dir)
        opts.fig_dir = fullfile(analyses_dir, 'figures', 'poster', analysis_date);
    end
    if opts.save && ~exist(opts.fig_dir, 'dir'); mkdir(opts.fig_dir); end

    if isfield(opts, 'estimates') && ~isempty(opts.estimates)
        estimates = opts.estimates;
    else
        estimates = loadEstimates(analysis_date, n_back_list, opts);
    end

    num_levels = 3;
    num_conds = 2;
    rb_ylim = [-20 20];

    for i = 1:numel(n_back_list)
        n_back = n_back_list(i);
        key = sprintf('n%d', n_back);
        est = estimates.(key);
        if isempty(est)
            warning('plotBiasDogGrid:noData', ...
                'Skipping n_back=%d (no estimates).', n_back);
            continue
        end

        delta = est.delta_theta_centers;
        if isempty(delta), delta = -90:1:90; end
        delta_smooth = linspace(-90, 90, 200);

        rb = est.rb;
        sd = est.sd;
        rb_ci = est.rb_ci;
        has_mu_ci = ~isempty(rb_ci) && isfield(rb_ci, 'mu_lo') && isfield(rb_ci, 'mu_hi');

        for cond = 1:num_conds
            fg = figure('Visible','off','Color',plt_opts.figure_color, ...
                'Name', sprintf('Poster Fig6 Bias DoG Grid %s %d-back', ...
                p.cond_names{cond}, n_back));

            for prev_lvl = 1:num_levels
                for curr_lvl = 1:num_levels
                    subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);
                    hold on;

                    mu = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 1));
                    if has_mu_ci
                        mu_lo = squeeze(rb_ci.mu_lo(prev_lvl, curr_lvl, cond, :));
                        mu_hi = squeeze(rb_ci.mu_hi(prev_lvl, curr_lvl, cond, :));
                    else
                        mu_lo = [];
                        mu_hi = [];
                    end

                    sd_params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
                    baseline = [];
                    if ~any(isnan(sd_params))
                        baseline = sd_params(3);
                    end

                    plotResponseBias(delta, mu, plt_opts, cond, mu_lo, mu_hi, baseline);

                    if ~any(isnan(sd_params))
                        dog_fit = calcDoG(delta_smooth, sd_params(1:3));
                        if isfield(plt_opts, 'rb_subtract_baseline') && plt_opts.rb_subtract_baseline
                            dog_fit = dog_fit - sd_params(3);
                        end
                        plot(delta_smooth, dog_fit, '--', ...
                            'LineWidth', plt_opts.line_width, ...
                            'Color', plt_opts.colors.black);

                        amp_est  = sd_params(1);
                        w_est    = max(sd_params(2), eps);
                        fwhm_est = (2 * sqrt(log(2))) / w_est;
                        text(0.95, 0.10, sprintf('A = %.2f°\nFWHM = %.1f°\nb = %.2f°', ...
                            amp_est, fwhm_est, sd_params(3)), ...
                            'Units','normalized','HorizontalAlignment','right', ...
                            'VerticalAlignment','bottom','FontSize',8, ...
                            'Color', plt_opts.colors.black);

                        curr_r2 = sd.all.r2(prev_lvl, curr_lvl, cond);
                        if ~isnan(curr_r2)
                            text(0.05, 0.90, sprintf('R^2 = %.2f', curr_r2), ...
                                'Units','normalized','HorizontalAlignment','left', ...
                                'VerticalAlignment','top','FontSize',8, ...
                                'Color', plt_opts.colors.black);
                        end
                    end

                    if cond == 1
                        fg_title = sprintf('%s \\rightarrow %s', ...
                            p.contrast{prev_lvl}, p.contrast{curr_lvl});
                    else
                        fg_title = sprintf('%s \\rightarrow %s', ...
                            p.precision{prev_lvl}, p.precision{curr_lvl});
                    end
                    title(fg_title);

                    axis square;
                    xlim([-90 90]); ylim(rb_ylim);
                    xticks(-90:45:90); xtickangle(0);
                    if prev_lvl == num_levels
                        xlabel('\Delta\theta (°)');
                    else
                        set(gca, 'XTickLabel', []);
                    end
                    if curr_lvl == 1
                        ylabel('Bias (°)');
                    end
                    set(gca, 'TickDir','out', ...
                        'TickLength', [plt_opts.tick_length, plt_opts.tick_length]);
                    line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');
                    line([0, 0], rb_ylim, 'LineWidth', 1, 'Color', 'k');
                    yticks(rb_ylim(1):5:rb_ylim(2));
                    box off;
                end
            end

            sgtitle(sprintf('%s full 3x3 - %d-back', p.cond_names{cond}, n_back));

            if opts.save
                fname = sprintf('Poster Fig6 Bias DoG Grid %s %d_back.%s', ...
                    p.cond_names{cond}, n_back, plt_opts.fg_type);
                saveas(fg, fullfile(opts.fig_dir, fname));
            end
            close(fg);
        end
    end
end

function p = defaultLabels()
    p.cond_names = {'Contrast','Precision'};
    p.contrast   = {'90%','50%','25%'};
    p.precision  = {'2°','40°','80°'};
end
