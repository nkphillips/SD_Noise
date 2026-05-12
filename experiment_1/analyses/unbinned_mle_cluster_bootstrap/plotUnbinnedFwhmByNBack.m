function plotUnbinnedFwhmByNBack(fig_dir, results, n_back_list, contrast_lbl, precision_lbl, ps, ci_pct, ci_method)
% plotUnbinnedFwhmByNBack  DoG-lobe FWHM decay across n-back.
%
% Renders two 2 x 3 grids:
%   1) columns are past level; within each panel, lines are current levels.
%   2) columns are current level; within each panel, lines are past levels.
% Points are horizontally dodged at each n-back to keep bootstrap CIs readable.

    if nargin < 7 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end
    if nargin < 8 || isempty(ci_method)
        ci_method = local_resultsCIMethod(results);
    end
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    plot_opts = local_buildPlotOpts(ps);
    n_back_list = n_back_list(:)';
    n_nb = numel(n_back_list);
    num_levels = 3;
    num_conds = 2;

    fwhm = nan(num_levels, num_levels, num_conds, n_nb);
    fwhm_lo = fwhm;
    fwhm_hi = fwhm;
    for i_nb = 1:n_nb
        res_i = local_getResult(results, n_back_list(i_nb), i_nb);
        if isempty(res_i) || ~isfield(res_i, 'summary_table') || isempty(res_i.summary_table)
            continue
        end
        pack = packUnbinnedScatterParams(res_i.summary_table);
        fwhm(:, :, :, i_nb) = pack.fwhm;
        fwhm_lo(:, :, :, i_nb) = pack.fwhm_lo;
        fwhm_hi(:, :, :, i_nb) = pack.fwhm_hi;
    end

    all_y = [fwhm(:); fwhm_lo(:); fwhm_hi(:)];
    all_y = all_y(isfinite(all_y));
    if isempty(all_y)
        y_lims = [0 1];
    else
        y_min = min(all_y);
        y_max = max(all_y);
        if y_min == y_max
            y_pad = max(abs(y_min), 1) * 0.1;
        else
            y_pad = 0.10 * (y_max - y_min);
        end
        y_lims = [max(0, y_min - y_pad), y_max + y_pad];
    end

    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s CI', ci_nominal_pct, local_ciLabel(ci_method));

    local_renderByFixedAxis(fig_dir, fwhm, fwhm_lo, fwhm_hi, n_back_list, y_lims, ...
        contrast_lbl, precision_lbl, plot_opts, ci_label, ...
        'past', 'current', 'unbinned_mle_sd_fwhm_by_nback_by_past_level.pdf');

    local_renderByFixedAxis(fig_dir, fwhm, fwhm_lo, fwhm_hi, n_back_list, y_lims, ...
        contrast_lbl, precision_lbl, plot_opts, ci_label, ...
        'current', 'past', 'unbinned_mle_sd_fwhm_by_nback_by_current_level.pdf');

end

function local_renderByFixedAxis(fig_dir, fwhm, fwhm_lo, fwhm_hi, n_back_list, y_lims, ...
    contrast_lbl, precision_lbl, plot_opts, ci_label, fixed_axis, line_axis, out_name)

    num_levels = 3;
    num_conds = 2;
    fig = figure('Color', plot_opts.figure_color, 'Visible', 'off', ...
        'Units', 'inches', 'Position', [1 1 11 7], ...
        'PaperUnits', 'inches', 'PaperSize', [11 7], ...
        'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 7]);
    tl = tiledlayout(num_conds, num_levels, 'Padding', 'compact', 'TileSpacing', 'compact');

    dodge = linspace(-0.16, 0.16, num_levels);

    for cond = 1:num_conds
        if cond == 1
            base_color = plot_opts.colors.blue;
            labels = contrast_lbl;
            row_lbl = 'Contrast';
        else
            base_color = plot_opts.colors.green;
            labels = precision_lbl;
            row_lbl = 'Precision';
        end
        line_colors = bsxfun(@times, repmat(base_color, num_levels, 1), [1; 0.70; 0.25]);

        for fixed_level = 1:num_levels
            ax = nexttile(tl);
            hold(ax, 'on');

            for line_level = 1:num_levels
                if strcmpi(fixed_axis, 'past')
                    prev = fixed_level;
                    curr = line_level;
                else
                    prev = line_level;
                    curr = fixed_level;
                end

                y = squeeze(fwhm(prev, curr, cond, :))';
                lo = squeeze(fwhm_lo(prev, curr, cond, :))';
                hi = squeeze(fwhm_hi(prev, curr, cond, :))';
                ok = isfinite(y);
                if ~any(ok)
                    continue
                end
                x = n_back_list + dodge(line_level);
                yneg = max(0, y - lo);
                ypos = max(0, hi - y);
                errorbar(ax, x(ok), y(ok), yneg(ok), ypos(ok), '-o', ...
                    'Color', line_colors(line_level, :) * 0.85, ...
                    'MarkerFaceColor', line_colors(line_level, :), ...
                    'MarkerEdgeColor', [1 1 1], ...
                    'MarkerSize', 5, ...
                    'LineWidth', plot_opts.line_width, ...
                    'CapSize', 0, ...
                    'DisplayName', sprintf('%s %s', local_titleCase(line_axis), labels{line_level}));
            end

            xlim(ax, [min(n_back_list) - 0.45, max(n_back_list) + 0.45]);
            ylim(ax, y_lims);
            xticks(ax, n_back_list);
            axis(ax, 'square');
            box(ax, 'off');
            set(ax, 'TickDir', 'out', 'TickLength', [plot_opts.tick_length, plot_opts.tick_length], ...
                'LineWidth', plot_opts.line_width, 'FontSize', plot_opts.axes_tick_font_size);

            title(ax, sprintf('%s: %s %s', row_lbl, fixed_axis, labels{fixed_level}), ...
                'FontSize', plot_opts.axes_tick_font_size, 'Interpreter', 'none');
            if cond == num_conds
                xlabel(ax, 'n-back', 'FontSize', plot_opts.axes_label_font_size);
            else
                set(ax, 'XTickLabel', []);
            end
            if fixed_level == 1
                ylabel(ax, 'Width FWHM (deg)', 'FontSize', plot_opts.axes_label_font_size);
            end
            if fixed_level == num_levels
                legend(ax, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 9);
            end
        end
    end

    title(tl, sprintf('Unbinned MLE: DoG-lobe FWHM across n-back, subplots by %s level (%s)', ...
        fixed_axis, ci_label), 'FontSize', plot_opts.axes_label_font_size + 1);

    out_pdf = fullfile(fig_dir, out_name);
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig);
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
    try
        plot_opts.axes_label_font_size = ps.axes_label_font_size;
        plot_opts.axes_tick_font_size = ps.axes_tick_font_size;
    catch
    end
end

function s = local_titleCase(s)
    if isempty(s)
        return
    end
    s = char(s);
    s(1) = upper(s(1));
end
