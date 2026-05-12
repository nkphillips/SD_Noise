function out = plotTargetedDoGEndpointEffectsByNBack(fig_dir, results, n_back_list, contrast_lbl, precision_lbl, ps, varargin)
% plotTargetedDoGEndpointEffectsByNBack  Across-n-back dumbbell + forest summary.
%
% Rows are n-back values. The left column shows absolute endpoint values for
% L1/L3; the right column shows the corresponding targeted delta tests.

    ip = inputParser;
    addParameter(ip, 'show_middle_level', false, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'axis_name', 'average_prev_curr', @(x) ischar(x) || isstring(x));
    parse(ip, varargin{:});

    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    axis_name = char(ip.Results.axis_name);
    plot_opts = local_buildPlotOpts(ps);
    out = struct();
    out.amplitude_pdf = local_renderMetric(fig_dir, results, n_back_list, 'A', ...
        'Amplitude A (deg)', 'targeted_dog_endpoint_effects_by_nback_amplitude.pdf', ...
        contrast_lbl, precision_lbl, plot_opts, axis_name, ip.Results.show_middle_level);
    out.fwhm_pdf = local_renderMetric(fig_dir, results, n_back_list, 'FWHM', ...
        'FWHM (deg)', 'targeted_dog_endpoint_effects_by_nback_fwhm.pdf', ...
        contrast_lbl, precision_lbl, plot_opts, axis_name, ip.Results.show_middle_level);
end

function out_pdf = local_renderMetric(fig_dir, results, n_back_list, parameter, x_label, fname, ...
    contrast_lbl, precision_lbl, plot_opts, axis_name, show_middle_level)

    n_back_list = n_back_list(:)';
    packs = cell(numel(n_back_list), 1);
    tests = cell(numel(n_back_list), 1);
    abs_vals = [];
    delta_vals = [];
    ci_method = local_resultsCIMethod(results);

    for i_nb = 1:numel(n_back_list)
        res_i = local_getResult(results, n_back_list(i_nb), i_nb);
        if isempty(res_i) || ~isfield(res_i, 'summary_table') || isempty(res_i.summary_table)
            continue
        end
        packs{i_nb} = packUnbinnedScatterParams(res_i.summary_table);
        tests{i_nb} = computeTargetedDoGHypothesisTests(res_i, n_back_list(i_nb));
        abs_vals = [abs_vals; local_collectAbsoluteValues(packs{i_nb}, parameter, axis_name, show_middle_level)]; %#ok<AGROW>
        delta_vals = [delta_vals; local_collectDeltaValues(tests{i_nb}, parameter)]; %#ok<AGROW>
    end

    fig_h = max(5.5, 2.3 * numel(n_back_list));
    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 fig_h], ...
        'PaperUnits', 'inches', 'PaperSize', [11 fig_h], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 fig_h]);
    tl = tiledlayout(numel(n_back_list), 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    abs_lims = local_limits(abs_vals, false);
    delta_lims = local_limits(delta_vals, true);

    for i_nb = 1:numel(n_back_list)
        ax_abs = nexttile(tl);
        if ~isempty(packs{i_nb})
            vals = local_metricValues(packs{i_nb}, parameter);
            local_plotDumbbell(ax_abs, vals, axis_name, abs_lims, x_label, ...
                contrast_lbl, precision_lbl, plot_opts, show_middle_level, n_back_list(i_nb));
        end
        title(ax_abs, 'Absolute endpoints', 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');

        ax_delta = nexttile(tl);
        if ~isempty(tests{i_nb})
            local_plotForest(ax_delta, tests{i_nb}, parameter, axis_name, delta_lims, plot_opts);
        end
        title(ax_delta, 'Tested deltas', 'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
    end

    title(tl, sprintf('Targeted %s endpoint effects across n-back (%s; 95%% %s CI)', ...
        parameter, local_axisLabel(axis_name), local_ciLabel(ci_method)), ...
        'FontSize', plot_opts.axes_label_font_size + 1, 'Interpreter', 'none');

    out_pdf = fullfile(fig_dir, fname);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
end

function local_plotDumbbell(ax, vals, axis_name, x_lims, x_label, contrast_lbl, precision_lbl, plot_opts, show_middle_level, n_back)
    hold(ax, 'on');
    y = [2, 1];
    labels = {'Contrast', 'Precision'};
    colors = [plot_opts.colors.blue; plot_opts.colors.green];
    level_labels = {contrast_lbl, precision_lbl};

    for m = 1:2
        v1 = local_axisLevelValue(vals(:, :, m), axis_name, 1);
        v2 = local_axisLevelValue(vals(:, :, m), axis_name, 2);
        v3 = local_axisLevelValue(vals(:, :, m), axis_name, 3);
        plot(ax, [v1, v3], [y(m), y(m)], '-', 'Color', colors(m, :) * 0.85, ...
            'LineWidth', plot_opts.line_width + 1, 'HandleVisibility', 'off');
        scatter(ax, v1, y(m), 55, 'MarkerFaceColor', colors(m, :), ...
            'MarkerEdgeColor', [1 1 1], 'LineWidth', plot_opts.line_width, ...
            'DisplayName', sprintf('%s high (%s)', labels{m}, level_labels{m}{1}));
        scatter(ax, v3, y(m), 55, 'MarkerFaceColor', colors(m, :) * 0.45, ...
            'MarkerEdgeColor', [1 1 1], 'LineWidth', plot_opts.line_width, ...
            'DisplayName', sprintf('%s low (%s)', labels{m}, level_labels{m}{3}));
        if show_middle_level
            scatter(ax, v2, y(m), 40, 'MarkerFaceColor', colors(m, :) * 0.65 + 0.35, ...
                'MarkerEdgeColor', colors(m, :) * 0.9, 'MarkerFaceAlpha', 0.35, ...
                'LineWidth', plot_opts.line_width, ...
                'DisplayName', sprintf('%s middle (%s)', labels{m}, level_labels{m}{2}));
        end
    end

    text(ax, x_lims(1), 2.45, sprintf('%d-back', n_back), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', plot_opts.axes_tick_font_size, 'FontWeight', 'bold');
    set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(labels));
    xlabel(ax, x_label, 'FontSize', plot_opts.axes_label_font_size);
    xlim(ax, x_lims);
    ylim(ax, [0.4, 2.6]);
    local_styleAxis(ax, plot_opts);
end

function local_plotForest(ax, tests, parameter, axis_name, x_lims, plot_opts)
    hold(ax, 'on');
    y = [5, 4, 3, 2, 1];
    row_names = local_effectNames(parameter);
    labels = local_effectLabels(parameter);
    colors = local_forestColors(plot_opts, numel(row_names));

    xline(ax, 0, 'k-', 'HandleVisibility', 'off');
    for i = 1:numel(row_names)
        row = local_getTestRow(tests, row_names{i}, parameter, axis_name);
        if isempty(row)
            continue
        end
        est = row.estimate(1);
        lo = row.bca_lo(1);
        hi = row.bca_hi(1);
        ok = local_isValid(row);
        if ok
            marker_face = colors(i, :);
            line_style = '-';
        else
            marker_face = [1 1 1];
            line_style = '--';
        end
        if isfinite(lo) && isfinite(hi)
            plot(ax, [lo, hi], [y(i), y(i)], line_style, ...
                'Color', colors(i, :) * 0.85, 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
        end
        if isfinite(est)
            scatter(ax, est, y(i), 40, 'MarkerFaceColor', marker_face, ...
                'MarkerEdgeColor', colors(i, :), 'LineWidth', plot_opts.line_width, ...
                'HandleVisibility', 'off');
            sig_label = local_sigLabel(row);
            if ~isempty(sig_label)
                dx = 0.02 * diff(x_lims);
                text(ax, est + dx, y(i) + 0.18, sig_label, 'FontSize', 8, ...
                    'Color', colors(i, :) * 0.75, 'HorizontalAlignment', 'left', ...
                    'Interpreter', 'none');
            end
        end
    end

    set(ax, 'YTick', fliplr(y), 'YTickLabel', fliplr(labels));
    xlabel(ax, sprintf('Delta %s (L3 - L1)', parameter), ...
        'FontSize', plot_opts.axes_label_font_size, 'Interpreter', 'none');
    xlim(ax, x_lims);
    ylim(ax, [0.4, 5.6]);
    local_styleAxis(ax, plot_opts);
end

function vals = local_metricValues(pack, parameter)
    switch parameter
        case 'A'
            vals = pack.amp;
        case 'FWHM'
            vals = pack.fwhm;
        otherwise
            vals = nan(3, 3, 2);
    end
end

function vals_out = local_collectAbsoluteValues(pack, parameter, axis_name, show_middle_level)
    vals = local_metricValues(pack, parameter);
    levels = [1, 3];
    if show_middle_level
        levels = [1, 2, 3];
    end
    vals_out = [];
    for m = 1:2
        for i_level = 1:numel(levels)
            vals_out(end+1, 1) = local_axisLevelValue(vals(:, :, m), axis_name, levels(i_level)); %#ok<AGROW>
        end
    end
end

function vals_out = local_collectDeltaValues(tests, parameter)
    vals_out = [];
    if isempty(tests) || height(tests) == 0
        return
    end
    keep = strcmp(cellstr(tests.parameter), parameter);
    rows = tests(keep, :);
    for i = 1:height(rows)
        vals_out = [vals_out; rows.estimate(i); rows.bca_lo(i); rows.bca_hi(i)]; %#ok<AGROW>
    end
end

function v = local_axisLevelValue(mat, axis_name, level)
    switch axis_name
        case 'previous_level'
            v = mean(mat(level, :), 'omitnan');
        case 'current_level'
            v = mean(mat(:, level), 'omitnan');
        case 'average_prev_curr'
            v = 0.5 .* (mean(mat(level, :), 'omitnan') + mean(mat(:, level), 'omitnan'));
        otherwise
            v = NaN;
    end
end

function row = local_getTestRow(tests, effect_name, parameter, axis_name)
    if contains(effect_name, 'boundary_interaction')
        target = sprintf('%s.%s', effect_name, 'boundary_interaction');
    else
        target = sprintf('%s.%s', effect_name, axis_name);
    end
    keep = strcmp(cellstr(tests.name), target) & strcmp(cellstr(tests.parameter), parameter);
    if any(keep)
        row = tests(find(keep, 1), :);
    else
        row = table();
    end
end

function names = local_effectNames(parameter)
    switch parameter
        case 'A'
            names = {'contrast_A_high_to_low', 'precision_A_high_to_low', ...
                'contrast_minus_precision_A', 'contrast_A_boundary_interaction', ...
                'precision_A_boundary_interaction'};
        case 'FWHM'
            names = {'contrast_FWHM_high_to_low', 'precision_FWHM_high_to_low', ...
                'precision_minus_contrast_FWHM', 'contrast_FWHM_boundary_interaction', ...
                'precision_FWHM_boundary_interaction'};
        otherwise
            names = {};
    end
end

function labels = local_effectLabels(parameter)
    switch parameter
        case 'A'
            labels = {'Contrast', 'Precision', 'Contrast - precision', ...
                'Contrast boundary', 'Precision boundary'};
        case 'FWHM'
            labels = {'Contrast', 'Precision', 'Precision - contrast', ...
                'Contrast boundary', 'Precision boundary'};
        otherwise
            labels = {};
    end
end

function colors = local_forestColors(plot_opts, n)
    base = [plot_opts.colors.blue; plot_opts.colors.green; plot_opts.colors.purple; ...
        plot_opts.colors.orange; plot_opts.colors.gray];
    if n <= size(base, 1)
        colors = base(1:n, :);
    else
        colors = repmat(base, ceil(n / size(base, 1)), 1);
        colors = colors(1:n, :);
    end
end

function res_i = local_getResult(results, n_back, i_nb)
    if isstruct(results) && isfield(results, sprintf('n%d', n_back))
        res_i = results.(sprintf('n%d', n_back));
    elseif isstruct(results) && isfield(results, 'summary_table') && i_nb == 1
        res_i = results;
    else
        res_i = [];
    end
end

function tf = local_isValid(row)
    if ismember('valid_for_inference', row.Properties.VariableNames)
        tf = logical(row.valid_for_inference(1));
    else
        tf = true;
    end
end

function s = local_sigLabel(row)
    s = '';
    if ~local_isValid(row) || ~ismember('p_holm', row.Properties.VariableNames) || ...
            ~isfinite(row.p_holm(1)) || row.p_holm(1) >= 0.05
        return
    end
    if ismember('p_holm_label', row.Properties.VariableNames)
        s = sprintf('Holm p %s', char(row.p_holm_label(1)));
    else
        s = sprintf('Holm p %.3f', row.p_holm(1));
    end
end

function lims = local_limits(vals, include_zero)
    vals = vals(isfinite(vals));
    if isempty(vals)
        lims = [-1 1];
        return
    end
    if include_zero
        vals = [vals; 0];
    end
    vmin = min(vals);
    vmax = max(vals);
    if vmin == vmax
        pad = max(abs(vmin), 1) * 0.1;
    else
        pad = 0.12 * (vmax - vmin);
    end
    lims = [vmin - pad, vmax + pad];
end

function s = local_axisLabel(axis_name)
    switch axis_name
        case 'previous_level'
            s = 'previous level';
        case 'current_level'
            s = 'current level';
        case 'average_prev_curr'
            s = 'average prev/current';
        otherwise
            s = axis_name;
    end
end

function local_styleAxis(ax, plot_opts)
    axis(ax, 'square');
    box(ax, 'off');
    set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
        'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);
end

function ci_method = local_resultsCIMethod(results)
    ci_method = 'bootstrap';
    if isstruct(results) && isfield(results, 'ci_method') && ~isempty(results.ci_method)
        ci_method = results.ci_method;
        return
    end
    if isstruct(results)
        fields = fieldnames(results);
        for i = 1:numel(fields)
            val = results.(fields{i});
            if isstruct(val) && isfield(val, 'ci_method') && ~isempty(val.ci_method)
                ci_method = val.ci_method;
                return
            end
        end
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

function plot_opts = local_buildPlotOpts(ps)
    plot_opts.line_width = ps.line_width;
    plot_opts.colors = ps.colors;
    plot_opts.figure_color = ps.figure_color;
    plot_opts.tick_length = ps.tick_length;
    plot_opts.axes_label_font_size = 14;
    plot_opts.axes_tick_font_size = 13;
    if isfield(ps, 'axes_label_font_size')
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
    end
    if isfield(ps, 'axes_tick_font_size')
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    end
    if ~isfield(plot_opts.colors, 'purple')
        plot_opts.colors.purple = [102 51 204] ./ 255;
    end
    if ~isfield(plot_opts.colors, 'gray')
        plot_opts.colors.gray = [128 128 128] ./ 255;
    end
end
