function out = plotUnbinnedSimpleSlopeTrendTestsByNBack(fig_dir, simple_tests, ps, ci_method, contrast_lbl, precision_lbl)
% plotUnbinnedSimpleSlopeTrendTestsByNBack  Forest plots for subplot-level slopes.

    if nargin < 4 || isempty(ci_method)
        ci_method = local_tableCIMethod(simple_tests);
    end
    if nargin < 5 || isempty(contrast_lbl)
        contrast_lbl = {'L1', 'L2', 'L3'};
    end
    if nargin < 6 || isempty(precision_lbl)
        precision_lbl = {'L1', 'L2', 'L3'};
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    out = struct();
    out.amplitude_by_past_pdf = local_render(fig_dir, simple_tests, 'A', 'previous_level', ...
        'Current-level slope within fixed previous level', 'A slope (deg / current-level step)', ...
        'simple_slope_trends_amplitude_by_past_level.pdf', ps, ci_method, contrast_lbl, precision_lbl);
    out.amplitude_by_current_pdf = local_render(fig_dir, simple_tests, 'A', 'current_level', ...
        'Previous-level slope within fixed current level', 'A slope (deg / previous-level step)', ...
        'simple_slope_trends_amplitude_by_current_level.pdf', ps, ci_method, contrast_lbl, precision_lbl);
    out.fwhm_by_past_pdf = local_render(fig_dir, simple_tests, 'FWHM', 'previous_level', ...
        'Current-level slope within fixed previous level', 'FWHM slope (deg / current-level step)', ...
        'simple_slope_trends_fwhm_by_past_level.pdf', ps, ci_method, contrast_lbl, precision_lbl);
    out.fwhm_by_current_pdf = local_render(fig_dir, simple_tests, 'FWHM', 'current_level', ...
        'Previous-level slope within fixed current level', 'FWHM slope (deg / previous-level step)', ...
        'simple_slope_trends_fwhm_by_current_level.pdf', ps, ci_method, contrast_lbl, precision_lbl);
end

function out_pdf = local_render(fig_dir, T, parameter, fixed_axis, title_str, x_label, fname, ps, ci_method, contrast_lbl, precision_lbl)
    keep = strcmp(cellstr(T.parameter), parameter) & strcmp(cellstr(T.fixed_axis), fixed_axis);
    T = T(keep, :);
    if isempty(T) || height(T) == 0
        out_pdf = '';
        return
    end

    n_back_list = unique(T.n_back(:))';
    manips = {'contrast', 'precision'};
    manip_labels = {'Contrast', 'Precision'};
    fixed_levels = 1:3;
    colors = local_colors(ps);
    plot_opts = local_plotOpts(ps);

    x_vals = [T.estimate; T.bca_lo; T.bca_hi];
    x_vals = x_vals(isfinite(x_vals));
    if isempty(x_vals)
        x_lims = [-1, 1];
    else
        pad = 0.12 * max(max(x_vals) - min(x_vals), eps);
        x_lims = [min(x_vals) - pad, max(x_vals) + pad];
        if x_lims(1) > 0, x_lims(1) = 0; end
        if x_lims(2) < 0, x_lims(2) = 0; end
    end

    fig_h = max(5.5, 2.0 * numel(n_back_list));
    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 fig_h], ...
        'PaperUnits', 'inches', 'PaperSize', [11 fig_h], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 fig_h]);
    tl = tiledlayout(numel(n_back_list), numel(manips), 'Padding', 'compact', 'TileSpacing', 'compact');

    for i_nb = 1:numel(n_back_list)
        for i_m = 1:numel(manips)
            ax = nexttile(tl);
            level_labels = local_levelLabels(manips{i_m}, contrast_lbl, precision_lbl);
            local_plotPanel(ax, T, n_back_list(i_nb), manips{i_m}, fixed_levels, ...
                level_labels, colors(i_m, :), x_lims, x_label, plot_opts);
            if i_nb == 1
                title(ax, manip_labels{i_m}, 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
            end
            if i_m == 1
                ylabel(ax, sprintf('%d-back', n_back_list(i_nb)), ...
                    'FontSize', plot_opts.axes_label_font_size, 'FontWeight', 'bold');
            end
            if i_nb < numel(n_back_list)
                set(ax, 'XTickLabel', []);
            else
                xlabel(ax, x_label, 'FontSize', plot_opts.axes_label_font_size, 'Interpreter', 'none');
            end
        end
    end

    title(tl, sprintf('%s: %s (95%% %s CI)', local_parameterLabel(parameter), title_str, local_ciLabel(ci_method)), ...
        'FontSize', plot_opts.axes_label_font_size + 1, 'Interpreter', 'none');

    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function local_plotPanel(ax, T, n_back, manipulation, fixed_levels, level_labels, color, x_lims, x_label, plot_opts)
    hold(ax, 'on');
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    y = 1:numel(fixed_levels);
    labels = local_fixedLevelTickLabels(fixed_levels, level_labels);
    for i_level = 1:numel(fixed_levels)
        keep = T.n_back == n_back & strcmp(cellstr(T.manipulation), manipulation) & ...
            T.fixed_level == fixed_levels(i_level);
        if ~any(keep)
            continue
        end
        row = T(find(keep, 1), :);
        est = row.estimate(1);
        lo = row.bca_lo(1);
        hi = row.bca_hi(1);
        ok = local_isValid(row);
        shade = 0.45 + 0.22 * i_level;
        c = min(color .* shade, 1);
        if ok
            marker_face = c;
            line_style = '-';
        else
            marker_face = [1 1 1];
            line_style = '--';
        end
        if isfinite(lo) && isfinite(hi)
            plot(ax, [lo, hi], [y(i_level), y(i_level)], line_style, ...
                'Color', c * 0.85, 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
        end
        if isfinite(est)
            scatter(ax, est, y(i_level), 42, 'MarkerFaceColor', marker_face, ...
                'MarkerEdgeColor', c, 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
            sig_label = local_sigLabel(row);
            if ~isempty(sig_label)
                dx = 0.02 * diff(x_lims);
                text(ax, est + dx, y(i_level) + 0.16, sig_label, ...
                    'FontSize', 8, 'Color', c * 0.75, ...
                    'HorizontalAlignment', 'left', 'Interpreter', 'none');
            end
        end
    end
    set(ax, 'YTick', y, 'YTickLabel', labels);
    xlim(ax, x_lims);
    ylim(ax, [0.4, numel(fixed_levels) + 0.6]);
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
        'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
    xlabel(ax, x_label, 'FontSize', plot_opts.axes_label_font_size, 'Interpreter', 'none');
end

function labels = local_levelLabels(manipulation, contrast_lbl, precision_lbl)
    if strcmpi(manipulation, 'contrast')
        labels = contrast_lbl;
    else
        labels = precision_lbl;
    end
    if isstring(labels)
        labels = cellstr(labels);
    end
end

function labels = local_fixedLevelTickLabels(fixed_levels, level_labels)
    labels = cell(numel(fixed_levels), 1);
    for i = 1:numel(fixed_levels)
        lvl = fixed_levels(i);
        if numel(level_labels) >= lvl
            labels{i} = char(level_labels{lvl});
        else
            labels{i} = sprintf('L%d', lvl);
        end
    end
end

function tf = local_isValid(row)
    tf = true;
    if ismember('valid_for_inference', row.Properties.VariableNames)
        tf = tf && logical(row.valid_for_inference(1));
    end
    if ismember('supports_effect', row.Properties.VariableNames)
        tf = tf && logical(row.supports_effect(1));
    end
end

function s = local_sigLabel(row)
    s = '';
    if ismember('supports_effect', row.Properties.VariableNames) && ~logical(row.supports_effect(1))
        return
    end
    if ismember('p_fdr_bh_label', row.Properties.VariableNames)
        val = char(row.p_fdr_bh_label(1));
        if contains(val, '*')
            s = ['q=' val];
            return
        end
    end
    if ismember('p_bca_label', row.Properties.VariableNames)
        val = char(row.p_bca_label(1));
        if contains(val, '*')
            s = ['p=' val];
        end
    end
end

function c = local_colors(ps)
    blue = [0.20 0.45 0.75];
    green = [0.30 0.65 0.30];
    try
        blue = ps.colors.blue;
        green = ps.colors.green;
    catch
    end
    c = [blue; green];
end

function plot_opts = local_plotOpts(ps)
    plot_opts.figure_color = 'w';
    plot_opts.line_width = 1;
    plot_opts.tick_length = 0.020;
    plot_opts.axes_label_font_size = 14;
    plot_opts.axes_tick_font_size = 13;
    try
        plot_opts.figure_color = ps.figure_color;
        plot_opts.line_width = ps.line_width;
        plot_opts.tick_length = ps.tick_length;
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    catch
    end
end

function label = local_parameterLabel(parameter)
    switch parameter
        case 'A'
            label = 'Amplitude';
        case 'FWHM'
            label = 'FWHM';
        otherwise
            label = parameter;
    end
end

function ci_method = local_tableCIMethod(T)
    if ~isempty(T) && height(T) > 0 && ismember('ci_method', T.Properties.VariableNames)
        ci_method = char(T.ci_method(1));
    else
        ci_method = 'bootstrap';
    end
end

function s = local_ciLabel(ci_method)
    ci_method = lower(char(ci_method));
    switch ci_method
        case 'bca'
            s = 'BCa';
        case 'percentile'
            s = 'percentile';
        otherwise
            s = char(ci_method);
    end
end
