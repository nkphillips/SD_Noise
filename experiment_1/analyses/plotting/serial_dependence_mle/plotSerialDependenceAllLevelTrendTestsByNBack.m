function out = plotSerialDependenceAllLevelTrendTestsByNBack(fig_dir, trend_tests, ps, ci_method)
% plotSerialDependenceAllLevelTrendTestsByNBack  Forest plots for all-level trend tests.

    if nargin < 4 || isempty(ci_method)
        ci_method = local_tableCIMethod(trend_tests);
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    out = struct();
    out.amplitude_pdf = local_renderParameter(fig_dir, trend_tests, 'A', ...
        'Amplitude trend (deg/level)', 'all_level_trend_tests_amplitude_by_nback.pdf', ps, ci_method);
    out.fwhm_pdf = local_renderParameter(fig_dir, trend_tests, 'FWHM', ...
        'FWHM trend (deg / level)', 'all_level_trend_tests_fwhm_by_nback.pdf', ps, ci_method);
end

function out_pdf = local_renderParameter(fig_dir, T, parameter, x_label, fname, ps, ci_method)
    T = T(strcmp(cellstr(T.parameter), parameter), :);
    if isempty(T) || height(T) == 0
        out_pdf = '';
        return
    end

    n_back_list = unique(T.n_back(:))';
    terms = {'previous_slope', 'current_slope', 'previous_current_interaction'};
    term_labels = {'Previous-level slope', 'Current-level slope', 'Previous x current'};
    row_names = {'precision', 'contrast', 'precision_minus_contrast'};
    row_labels = {'Precision', 'Contrast', 'Precision - contrast'};
    colors = local_colors(ps);
    plot_opts = local_plotOpts(ps);

    x_vals = [T.estimate; T.bca_lo; T.bca_hi];
    x_vals = x_vals(isfinite(x_vals));
    if isempty(x_vals)
        x_lims = [-1, 1];
    else
        pad = 0.12 * max(diff([min(x_vals), max(x_vals)]), eps);
        x_lims = [min(x_vals) - pad, max(x_vals) + pad];
        if x_lims(1) == x_lims(2)
            x_lims = x_lims + [-1, 1];
        end
        if x_lims(1) > 0, x_lims(1) = 0; end
        if x_lims(2) < 0, x_lims(2) = 0; end
    end

    fig_h = max(5.5, 2.0 * numel(n_back_list));
    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 12 fig_h], ...
        'PaperUnits', 'inches', 'PaperSize', [12 fig_h], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 12 fig_h]);
    tl = tiledlayout(numel(n_back_list), numel(terms), 'Padding', 'compact', 'TileSpacing', 'compact');

    for i_nb = 1:numel(n_back_list)
        for i_term = 1:numel(terms)
            ax = nexttile(tl);
            local_plotOnePanel(ax, T, n_back_list(i_nb), terms{i_term}, row_names, row_labels, ...
                colors, x_lims, x_label, plot_opts);
            if i_nb == 1
                title(ax, term_labels{i_term}, 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
            end
            if i_term == 1
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

    title(tl, sprintf('All-level %s trend tests across n-back (95%% %s CI)', ...
        local_parameterLabel(parameter), local_ciLabel(ci_method)), ...
        'FontSize', plot_opts.axes_label_font_size + 1, 'Interpreter', 'none');

    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function local_plotOnePanel(ax, T, n_back, term, row_names, row_labels, colors, x_lims, ~, plot_opts)
    hold(ax, 'on');
    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    y = numel(row_names):-1:1;
    for i_row = 1:numel(row_names)
        keep = T.n_back == n_back & strcmp(cellstr(T.term), term) & ...
            strcmp(cellstr(T.manipulation), row_names{i_row});
        if ~any(keep)
            continue
        end
        row = T(find(keep, 1), :);
        est = row.estimate(1);
        lo = row.bca_lo(1);
        hi = row.bca_hi(1);
        if local_isSignificant(row)
            marker_face = colors(i_row, :);
        else
            marker_face = [1 1 1];
        end
        if isfinite(lo) && isfinite(hi)
            plot(ax, [lo, hi], [y(i_row), y(i_row)], '-', ...
                'Color', colors(i_row, :) * 0.85, 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
        end
        if isfinite(est)
            scatter(ax, est, y(i_row), 42, 'MarkerFaceColor', marker_face, ...
                'MarkerEdgeColor', colors(i_row, :), 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
            sig_label = local_sigLabel(row);
            if ~isempty(sig_label)
                dx = 0.02 * diff(x_lims);
                text(ax, est + dx, y(i_row) + 0.16, sig_label, ...
                    'FontSize', 8, 'Color', colors(i_row, :) * 0.75, ...
                    'HorizontalAlignment', 'left', 'Interpreter', 'none');
            end
        end
    end
    set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(row_labels));
    xlim(ax, x_lims);
    ylim(ax, [0.4, numel(row_names) + 0.6]);
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
        'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
end

function s = local_sigLabel(row)
    s = '';
    if ismember('p_fdr_bh_label', row.Properties.VariableNames)
        val = char(row.p_fdr_bh_label(1));
        if ~isempty(strtrim(val))
            s = ['q=' val];
            return
        end
    end
    if ismember('p_bca_label', row.Properties.VariableNames)
        val = char(row.p_bca_label(1));
        if ~isempty(strtrim(val))
            s = ['p=' val];
        end
    end
end

function tf = local_isSignificant(row)
    tf = false;
    if ismember('p_fdr_bh_label', row.Properties.VariableNames)
        tf = tf || contains(char(row.p_fdr_bh_label(1)), '*');
    end
    if ismember('p_bca_label', row.Properties.VariableNames)
        tf = tf || contains(char(row.p_bca_label(1)), '*');
    end
end

function c = local_colors(ps)
    blue = [0.20 0.45 0.75];
    green = [0.30 0.65 0.30];
    red = [0.75 0.25 0.20];
    try
        blue = ps.colors.blue;
        green = ps.colors.green;
        red = ps.colors.red;
    catch
    end
    c = [green; blue; red];
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
            label = 'amplitude';
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
