function plotAmplitudeWidthScatter(analysis_date, n_back_list, opts)
% plotAmplitudeWidthScatter
%
% Poster Figure 5 - Amplitude vs FWHM scatter across cells and n-backs.
% Two subplots (Contrast, Precision). Each point = one
% (prev_lvl, curr_lvl, n_back) cell; x = FWHM, y = amplitude.
% Marker size encodes n-back; color encodes current level; marker shape
% encodes previous level.
%
% Usage:
%   plotAmplitudeWidthScatter('11.24.2025', [1 2 3]);

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

    estimates = loadEstimates(analysis_date, n_back_list, opts);

    num_conds = 2;
    num_levels = 3;
    num_nback = numel(n_back_list);
    c = 2 * sqrt(log(2));

    prev_markers = {'o','s','^'}; % 3 marker shapes for prev_lvl
    nback_sizes  = linspace(40, 140, max(num_nback, 1));

    fg = figure('Visible','off','Color',plt_opts.figure_color, ...
                'Name','Poster Fig5 Amplitude vs FWHM', ...
                'Position', [100 100 1100 700]);
    tl = tiledlayout(fg, 2, num_conds, 'TileSpacing', 'compact', 'Padding', 'compact');

    for cond = 1:num_conds
        ax = nexttile(tl, cond);
        hold(ax, 'on');

        if cond == 1
            base_color = plt_opts.colors.blue;
            level_labels = p.contrast;
        else
            base_color = plt_opts.colors.green;
            level_labels = p.precision;
        end

        level_shades = linspace(0.25, 0.75, num_levels);

        for i = 1:num_nback
            key = sprintf('n%d', n_back_list(i));
            est = estimates.(key);
            if isempty(est); continue; end
            pe = est.sd.all.params_est;

            for prev_lvl = 1:num_levels
                for curr_lvl = 1:num_levels
                    A = pe(prev_lvl, curr_lvl, cond, 1);
                    w = pe(prev_lvl, curr_lvl, cond, 2);
                    if ~isfinite(A) || ~isfinite(w) || w <= 0
                        continue
                    end
                    fwhm = c / w;
                    shade = level_shades(curr_lvl);
                    mk_color = shade * base_color + (1 - shade) * [1 1 1];
                    scatter(ax, fwhm, A, nback_sizes(i), ...
                        'Marker', prev_markers{prev_lvl}, ...
                        'MarkerFaceColor', mk_color, ...
                        'MarkerEdgeColor', [0.3 0.3 0.3], ...
                        'MarkerFaceAlpha', 0.85, ...
                        'LineWidth', plt_opts.line_width);
                end
            end
        end

        title(ax, p.cond_names{cond});
        xlabel(ax, 'FWHM (°)');
        ylabel(ax, 'Amplitude (°)');
        axis(ax, 'square'); box(ax, 'off');
        set(ax, 'TickDir','out', 'LineWidth', plt_opts.line_width);

        legend_ax = nexttile(tl, cond + num_conds);
        drawLegendPanel(legend_ax, base_color, level_shades, level_labels, ...
            prev_markers, n_back_list, nback_sizes, p.cond_names{cond}, plt_opts);
    end

    title(tl, 'DoG amplitude vs FWHM');

    if opts.save
        fname = sprintf('Poster Fig5 Amplitude vs FWHM Scatter.%s', plt_opts.fg_type);
        exportgraphics(fg, fullfile(opts.fig_dir, fname), 'ContentType', 'vector');
    end
    close(fg);
end

function p = defaultLabels()
    p.cond_names = {'Contrast','Precision'};
    p.contrast   = {'90%','50%','25%'};
    p.precision  = {'2°','40°','80°'};
end

function drawLegendPanel(ax, base_color, level_shades, level_labels, prev_markers, n_back_list, nback_sizes, cond_name, plt_opts)
% drawLegendPanel
% Custom, in-bounds legend panel for the three encodings used in Fig 5.

    axes(ax); %#ok<LAXES>
    cla(ax);
    hold(ax, 'on');
    axis(ax, [0 1 0 1]);
    axis(ax, 'off');

    marker_edge = [0.3 0.3 0.3];
    text_x = 0.18;
    marker_x = 0.08;
    y = 0.90;
    dy = 0.12;

    text(ax, 0.02, y, sprintf('Current level (%s color)', cond_name), ...
        'FontWeight', 'bold', 'FontSize', 9, 'Units', 'normalized');
    y = y - dy;
    for lvl = 1:numel(level_labels)
        shade = level_shades(lvl);
        mk_color = shade * base_color + (1 - shade) * [1 1 1];
        scatter(ax, marker_x, y, 70, 'o', ...
            'MarkerFaceColor', mk_color, 'MarkerEdgeColor', marker_edge, ...
            'MarkerFaceAlpha', 0.85, 'LineWidth', plt_opts.line_width);
        text(ax, text_x, y, level_labels{lvl}, ...
            'VerticalAlignment', 'middle', 'FontSize', 9);
        y = y - dy;
    end

    y = y - 0.03;
    text(ax, 0.02, y, 'Previous level (marker shape)', ...
        'FontWeight', 'bold', 'FontSize', 9, 'Units', 'normalized');
    y = y - dy;
    neutral_color = [0.75 0.75 0.75];
    for lvl = 1:numel(level_labels)
        scatter(ax, marker_x, y, 70, prev_markers{lvl}, ...
            'MarkerFaceColor', neutral_color, 'MarkerEdgeColor', marker_edge, ...
            'MarkerFaceAlpha', 0.85, 'LineWidth', plt_opts.line_width);
        text(ax, text_x, y, level_labels{lvl}, ...
            'VerticalAlignment', 'middle', 'FontSize', 9);
        y = y - dy;
    end

    y = 0.90;
    marker_x = 0.58;
    text_x = 0.70;
    text(ax, 0.52, y, 'n-back (marker size)', ...
        'FontWeight', 'bold', 'FontSize', 9, 'Units', 'normalized');
    y = y - dy;
    for i = 1:numel(n_back_list)
        scatter(ax, marker_x, y, nback_sizes(i), 'o', ...
            'MarkerFaceColor', neutral_color, 'MarkerEdgeColor', marker_edge, ...
            'MarkerFaceAlpha', 0.85, 'LineWidth', plt_opts.line_width);
        text(ax, text_x, y, sprintf('%d-back', n_back_list(i)), ...
            'VerticalAlignment', 'middle', 'FontSize', 9);
        y = y - dy;
    end
end
