function plotDoGPastLevelGrid(sd, contrast_lbl, precision_lbl, opts)
% plotDoGPastLevelGrid
%
% 2 x 3 DoG-only grid: rows = contrast vs precision; columns = past level.
% Each panel draws the super-subject DoG for all three current levels (uses
% sd.all.params_est(prev, curr, cond, 1:3) = [A, w, baseline] with calcDoG).
%
% No file I/O: pass sd (and labels) from the caller. For serial-dependence
% MLE outputs, use
% sd.all.params_est = packSummaryTableToSdParamsEst(results.summary_table).
%
% Inputs:
%   sd              - struct with field sd.all.params_est, size [3 x 3 x 2 x 3..4]
%   contrast_lbl    - 1x3 cell, x / legend labels for manipulation 1
%   precision_lbl   - 1x3 cell, labels for manipulation 2
%
% opts fields (all optional except fig layout):
%   .n_back         - scalar for figure title (default 1)
%   .cond_names     - 1x2 cell (default {'Contrast','Precision'})
%   .title_tag      - short suffix for sgtitle, e.g. 'serial-dependence MLE (BCa)'
%   .save           - logical (default false)
%   .fig_dir        - if .save true, directory for PDF
%   .fname          - if nonempty and .save, filename; default DoG Past Level Grid <n_back>.pdf
%   .plt_opts       - struct from posterPlotOpts (default internal posterPlotOpts)

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'n_back') || isempty(opts.n_back)
        opts.n_back = 1;
    end
    if ~isfield(opts, 'cond_names') || isempty(opts.cond_names)
        opts.cond_names = {'Contrast', 'Precision'};
    end
    if ~isfield(opts, 'title_tag'), opts.title_tag = ''; end
    if ~isfield(opts, 'save'), opts.save = false; end
    if ~isfield(opts, 'fig_dir'), opts.fig_dir = ''; end
    if ~isfield(opts, 'fname'), opts.fname = ''; end
    if ~isfield(opts, 'plt_opts') || isempty(opts.plt_opts)
        plt_opts = posterPlotOpts();
    else
        plt_opts = opts.plt_opts;
    end

    if isempty(sd) || ~isfield(sd, 'all') || ~isfield(sd.all, 'params_est')
        warning('plotDoGPastLevelGrid:noSD', 'Missing sd.all.params_est; nothing to plot.');
        return
    end

    num_levels = 3;
    num_conds = 2;

    delta_smooth = linspace(-90, 90, 300);
    lev_scales = [1;0.70;0.25];

    y_all = [];
    fg = figure('Visible', 'off', 'Color', plt_opts.figure_color, ...
        'Units', 'inches', 'Position', [1 1 12 6.5], ...
        'PaperUnits', 'inches', 'PaperSize', [12 6.5], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 12 6.5], ...
        'Name', sprintf('DoG by past level %d-back', opts.n_back));

    n_panels = num_conds * num_levels;
    ax_handles = gobjects(n_panels, 1);

    for cond = 1:num_conds
        if cond == 1
            base_color = plt_opts.colors.blue;
            curr_lbl = contrast_lbl;
            past_lbl = contrast_lbl;
        else
            base_color = plt_opts.colors.green;
            curr_lbl = precision_lbl;
            past_lbl = precision_lbl;
        end
        curve_colors = repmat(base_color, num_levels, 1) .* lev_scales;

        for prev_lvl = 1:num_levels
            idx_ax = (cond - 1) * num_levels + prev_lvl;
            ax = subplot(num_conds, num_levels, idx_ax);
            ax_handles(idx_ax) = ax;
            hold(ax, 'on');

            for curr_lvl = 1:num_levels
                sd_params = squeeze(sd.all.params_est(prev_lvl, curr_lvl, cond, 1:3));
                if any(isnan(sd_params(:)))
                    continue
                end
                dog_fit = calcDoG(delta_smooth, sd_params(1:3));
                if isfield(plt_opts, 'rb_subtract_baseline') && plt_opts.rb_subtract_baseline
                    dog_fit = dog_fit - sd_params(3);
                end
                y_all = [y_all; dog_fit(:)]; %#ok<AGROW>
                plot(ax, delta_smooth, dog_fit, '-', ...
                    'Color', curve_colors(curr_lvl, :), ...
                    'LineWidth', plt_opts.line_width, ...
                    'DisplayName', sprintf('Current %s', curr_lbl{curr_lvl}));
            end

            axis(ax, 'square');
            xlim(ax, [-90 90]);
            set(ax, 'TickDir', 'out', ...
                'TickLength', [plt_opts.tick_length, plt_opts.tick_length], ...
                'LineWidth', plt_opts.line_width, 'Box', 'off', 'FontSize', 13);
            xticks(ax, -90:45:90);

            title(ax, sprintf('Past %s', past_lbl{prev_lvl}), 'FontSize', 11);

            if prev_lvl == 1
                ylabel(ax, 'Bias (°)', 'FontSize', 12);
            end
            if cond == num_conds
                xlabel(ax, '\Delta\theta (°)', 'FontSize', 12);
            else
                set(ax, 'XTickLabel', []);
            end

            if prev_lvl == num_levels
                legend(ax, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 9);
            end
        end
    end

    if ~isempty(y_all) && all(isfinite(y_all(:)))
        y_min = min(y_all(:));
        y_max = max(y_all(:));
        if y_min == y_max
            pad = max(abs(y_min), 1) * 0.1;
        else
            pad = 0.08 * (y_max - y_min);
        end
        y_lo = y_min - pad;
        y_hi = y_max + pad;
    else
        y_lo = [];
        y_hi = [];
    end

    for idx_ax = 1:n_panels
        ax = ax_handles(idx_ax);
        if ~isgraphics(ax)
            continue
        end
        if ~isempty(y_lo)
            ylim(ax, [y_lo, y_hi]);
        end
        line(ax, [-90, 90], [0, 0], 'LineWidth', 1, ...
            'Color', plt_opts.colors.black, 'HandleVisibility', 'off');
        yy = ylim(ax);
        line(ax, [0, 0], yy, 'LineWidth', 1, ...
            'Color', plt_opts.colors.black, 'HandleVisibility', 'off');
    end

    if ~isempty(opts.title_tag)
        sg = sprintf('%s | %s; %d-back', opts.title_tag, ...
            strjoin(opts.cond_names, ' / '), opts.n_back);
    else
        sg = sprintf('DoG by past level: %s; %d-back', ...
            strjoin(opts.cond_names, ' / '), opts.n_back);
    end
    sgtitle(fg, sg, 'FontSize', 13);

    if opts.save
        if isempty(opts.fig_dir)
            warning('plotDoGPastLevelGrid:noFigDir', ...
                'opts.save is true but opts.fig_dir is empty; figure not saved.');
        else
            if ~exist(opts.fig_dir, 'dir')
                mkdir(opts.fig_dir);
            end
            if ~isempty(opts.fname)
                out_pdf = fullfile(opts.fig_dir, opts.fname);
            else
                out_pdf = fullfile(opts.fig_dir, sprintf('DoG Past Level Grid %d_back.%s', ...
                    opts.n_back, plt_opts.fg_type));
            end
            saveas(fg, out_pdf);
        end
    end
    close(fg);
end
