function plotSigmaSummary(analysis_date, n_back_list, opts)
% plotSigmaSummary
%
% Poster Figure 2 - Response uncertainty sigma summary across the full
% previous-by-current 3x3 design. One figure is written per condition
% (Contrast, Precision). Rows are previous level, columns are current level,
% and each cell plots mean sigma across delta_theta windows as a function of
% n-back. This preserves all 3x3 comparisons instead of averaging over
% previous level or showing only the diagonal.
%
% Usage:
%   plotSigmaSummary('11.24.2025', [1 2 3]);

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'objective'),     opts.objective = 'sse'; end
    if ~isfield(opts, 'save'),          opts.save      = true;  end
    if ~isfield(opts, 'max_abs_delta'), opts.max_abs_delta = 60; end

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

    num_conds = 2;
    num_levels = 3;
    num_nback = numel(n_back_list);

    % Preallocate: [prev_lvl, curr_lvl, cond, nback]
    mean_sigma = nan(num_levels, num_levels, num_conds, num_nback);

    for i = 1:num_nback
        n_back = n_back_list(i);
        key = sprintf('n%d', n_back);
        est = estimates.(key);
        if isempty(est); continue; end

        delta = est.delta_theta_centers;
        if isempty(delta), delta = -90:1:90; end
        keep_win = abs(delta) <= opts.max_abs_delta;

        sigma_arr = squeeze(est.rb.all.params_est(:,:,:,:,2));
        % sigma_arr: [prev_lvl, curr_lvl, cond, window]

        for cond = 1:num_conds
            for prev_lvl = 1:num_levels
                for curr_lvl = 1:num_levels
                    slab = squeeze(sigma_arr(prev_lvl, curr_lvl, cond, keep_win));
                    mean_sigma(prev_lvl, curr_lvl, cond, i) = mean(slab(:), 'omitnan');
                end
            end
        end
    end

    for cond = 1:num_conds
        if cond == 1
            line_color = plt_opts.colors.blue;
            level_labels = p.contrast;
        else
            line_color = plt_opts.colors.green;
            level_labels = p.precision;
        end

        fg = figure('Visible','off','Color',plt_opts.figure_color, ...
                    'Name', sprintf('Poster Fig 2 - Sigma summary %s', p.cond_names{cond}));

        all_y = squeeze(mean_sigma(:,:,cond,:));
        ymax_ = max(all_y(:), [], 'omitnan');
        if isnan(ymax_) || ymax_ <= 0, ymax_ = 10; end
        y_lim = [0 ceil(ymax_ * 1.15)];

        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);
                hold on;

                y = squeeze(mean_sigma(prev_lvl, curr_lvl, cond, :));
                plot(n_back_list(:), y, '-', 'Color', line_color, ...
                    'LineWidth', plt_opts.line_width, 'HandleVisibility', 'off');
                scatter(n_back_list(:), y, plt_opts.marker_size_scatter, ...
                    'MarkerFaceColor', line_color, 'MarkerEdgeColor', [1 1 1], ...
                    'MarkerFaceAlpha', 0.9, 'LineWidth', plt_opts.line_width);

                title(sprintf('%s \\rightarrow %s', level_labels{prev_lvl}, level_labels{curr_lvl}));
                xlabel('n-back');
                if curr_lvl == 1
                    ylabel('Mean \sigma (°)');
                end
                xlim([min(n_back_list)-0.5, max(n_back_list)+0.5]);
                ylim(y_lim);
                xticks(n_back_list);
                axis square; box off;
                set(gca, 'TickDir','out', 'LineWidth', plt_opts.line_width);
            end
        end

        sgtitle(sprintf('\\sigma across n-back - %s', p.cond_names{cond}));

        if opts.save
            fname = sprintf('Poster Fig2 Sigma Summary %s.%s', ...
                p.cond_names{cond}, plt_opts.fg_type);
            saveas(fg, fullfile(opts.fig_dir, fname));
        end
        close(fg);
    end
end

function p = defaultLabels()
    p.cond_names = {'Contrast','Precision'};
    p.contrast   = {'90%','50%','25%'};
    p.precision  = {'2°','40°','80°'};
end
