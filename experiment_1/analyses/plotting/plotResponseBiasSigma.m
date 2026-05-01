function plotResponseBiasSigma(analysis_date, n_back_list, opts)
% plotResponseBiasSigma
%
% Poster Figure 1 - Response uncertainty sigma(\Delta\theta) grids.
% For each n-back in n_back_list and each condition (Contrast, Precision),
% renders a 3x3 grid (prev_lvl x curr_lvl) of the super-subject sigma as a
% function of delta_theta, overlaying the bootstrap sigma CI band from
% rb_ci.sigma_lo / rb_ci.sigma_hi (computed in Phase 0).
%
% Usage:
%   plotResponseBiasSigma('11.24.2025', [1 2 3]);
%   plotResponseBiasSigma(date, n_backs, struct('save', true));
%
% Optional opts:
%   opts.objective - 'sse' (default) or 'nll'
%   opts.save      - true/false (default true)
%   opts.fig_dir   - override output directory

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

    for i = 1:numel(n_back_list)
        n_back = n_back_list(i);
        key = sprintf('n%d', n_back);
        est = estimates.(key);
        if isempty(est)
            warning('plotResponseBiasSigma:noData', ...
                'Skipping n_back=%d (no estimates).', n_back);
            continue
        end

        delta = est.delta_theta_centers;
        if isempty(delta)
            delta = -90:1:90;
        end

        rb = est.rb;
        rb_ci = est.rb_ci;
        has_ci = ~isempty(rb_ci) && isfield(rb_ci, 'sigma_lo') && isfield(rb_ci, 'sigma_hi');

        [num_prev, num_curr, num_conds, ~, ~] = size(rb.all.params_est);

        for cond = 1:num_conds
            fg = figure('Visible', 'off', 'Color', plt_opts.figure_color, ...
                        'Name', sprintf('Super Subj Sigma %s (%d-back)', p.cond_names{cond}, n_back));

            for prev_lvl = 1:num_prev
                for curr_lvl = 1:num_curr
                    subplot(num_prev, num_curr, curr_lvl + (prev_lvl-1)*num_curr);

                    sigma_vals = squeeze(rb.all.params_est(prev_lvl, curr_lvl, cond, :, 2));
                    lo = [];
                    hi = [];
                    if cond == 1
                        lineColor = plt_opts.colors.blue;
                    else
                        lineColor = plt_opts.colors.green;
                    end

                    plot_band = false;
                    if has_ci
                        lo = squeeze(rb_ci.sigma_lo(prev_lvl, curr_lvl, cond, :));
                        hi = squeeze(rb_ci.sigma_hi(prev_lvl, curr_lvl, cond, :));
                        idx = isfinite(sigma_vals) & isfinite(lo) & isfinite(hi);
                        if any(idx)
                            plot_band = true;
                            err_upper = hi(idx) - sigma_vals(idx);
                            err_lower = sigma_vals(idx) - lo(idx);
                            err = [err_upper(:)'; err_lower(:)'];
                            hold on;
                            xlim([-90 90]);
                            sigma_ylim = getSigmaYLim(sigma_vals, lo, hi);
                            ylim(sigma_ylim);
                            line([0 0], sigma_ylim, ...
                                'LineWidth', plt_opts.line_width, ...
                                'Color', 0.8 * plt_opts.colors.white, ...
                                'HandleVisibility', 'off');
                            shadedErrorBar(delta(idx), sigma_vals(idx), err, ...
                                'lineProps', {'-', 'Color', lineColor, 'LineWidth', plt_opts.line_width}, ...
                                'transparent', true, 'patchSaturation', 0.4);
                        end
                    end
                    if ~plot_band
                        hold on;
                        xlim([-90 90]);
                        sigma_ylim = getSigmaYLim(sigma_vals, lo, hi);
                        ylim(sigma_ylim);
                        line([0 0], sigma_ylim, ...
                            'LineWidth', plt_opts.line_width, ...
                            'Color', 0.8 * plt_opts.colors.white, ...
                            'HandleVisibility', 'off');
                        plot(delta, sigma_vals, 'LineWidth', plt_opts.line_width, 'Color', lineColor);
                    end

                    axis square;
                    if cond == 1
                        fg_title = [p.contrast{prev_lvl} ' \rightarrow ' p.contrast{curr_lvl}];
                    else
                        fg_title = [p.precision{prev_lvl} ' \rightarrow ' p.precision{curr_lvl}];
                    end
                    title(fg_title);

                    xlim([-90 90]); ylim(sigma_ylim);
                    xticks(-90:45:90); xtickangle(0);
                    if prev_lvl == num_prev
                        xlabel('\Delta\theta (°)');
                    else
                        set(gca, 'XTickLabel', []);
                    end
                    if curr_lvl == 1
                        ylabel('\sigma (°)');
                    end
                    set(gca, 'TickDir', 'out', 'TickLength', ...
                        [plt_opts.tick_length, plt_opts.tick_length]);
                    box off;
                    hold on;
                end
            end

            sgtitle(sprintf('%s, n-back = %d', p.cond_names{cond}, n_back));

            if opts.save
                fname = sprintf('Super Subj Sigma %s %d_back.%s', ...
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

function y_lim = getSigmaYLim(sigma_vals, lo, hi)
    y_vals = sigma_vals(:);
    if ~isempty(lo)
        y_vals = [y_vals; lo(:); hi(:)];
    end
    y_vals = y_vals(isfinite(y_vals));

    if isempty(y_vals)
        y_lim = [0 1];
        return
    end

    y_min = min(y_vals);
    y_max = max(y_vals);
    y_range = y_max - y_min;
    if y_range <= 0
        y_range = max(abs(y_max), 1);
    end

    pad = 0.05 * y_range;
    y_lim = [max(0, y_min - pad), y_max + pad];
    if y_lim(1) == y_lim(2)
        y_lim = [max(0, y_lim(1) - 0.5), y_lim(2) + 0.5];
    end
end
