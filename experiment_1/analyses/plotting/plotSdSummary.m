function plotSdSummary(estimates, n_back_list, plt_opts, p, opts, param_index, ylabel_str, fname_stem, transform, ci_transform)
% plotSdSummary
%
% Shared renderer for Figure 3 (amplitude) and Figure 4 (FWHM/width).
% Renders one 3x3 grid per condition. Rows are previous level, columns are
% current level, and each cell plots the transformed DoG parameter across
% n-back with CI bars. This preserves all previous-by-current comparisons.
%
% Inputs:
%   estimates    - struct returned by loadEstimates.
%   n_back_list  - vector of n-backs to plot.
%   plt_opts     - struct from posterPlotOpts.
%   p            - label struct with p.cond_names, p.contrast, p.precision.
%   opts         - struct with optional fields:
%                    .fig_dir (mandatory for save)
%                    .save (default true)
%   param_index  - 1 = amplitude, 2 = width.
%   ylabel_str   - y-axis label.
%   fname_stem   - file stem (without extension).
%   transform    - function handle applied to raw parameter value.
%   ci_transform - function handle ci_transform(lo_raw, hi_raw) returning
%                  [lo_transformed, hi_transformed]. For amplitude this is
%                  identity; for FWHM it inverts (c/w) because FWHM is a
%                  monotone-decreasing function of w.

    num_conds = 2;
    num_levels = 3;
    num_nback = numel(n_back_list);

    vals = nan(num_levels, num_levels, num_conds, num_nback);
    lo   = nan(num_levels, num_levels, num_conds, num_nback);
    hi   = nan(num_levels, num_levels, num_conds, num_nback);

    for i = 1:num_nback
        key = sprintf('n%d', n_back_list(i));
        est = estimates.(key);
        if isempty(est); continue; end
        sd = est.sd;
        if isfield(est, 'sd_ci_cluster') && ~isempty(est.sd_ci_cluster) && ...
                isfield(est.sd_ci_cluster, 'lo') && ~isempty(est.sd_ci_cluster.lo)
            sd_ci = est.sd_ci_cluster;
            ci_label = 'subject-cluster 95% CI';
        else
            sd_ci = est.sd_ci;
            ci_label = '95% CI';
        end
        has_ci = ~isempty(sd_ci) && isfield(sd_ci, 'lo') && isfield(sd_ci, 'hi');

        pe = sd.all.params_est; % [prev, curr, cond, param]

        for cond = 1:num_conds
            for prev_lvl = 1:num_levels
                for curr_lvl = 1:num_levels
                    v = pe(prev_lvl, curr_lvl, cond, param_index);
                    vals(prev_lvl, curr_lvl, cond, i) = transform(v);

                    if has_ci
                        lo_raw = sd_ci.lo(prev_lvl, curr_lvl, cond, param_index);
                        hi_raw = sd_ci.hi(prev_lvl, curr_lvl, cond, param_index);
                        [lo_v, hi_v] = ci_transform(lo_raw, hi_raw);
                        lo(prev_lvl, curr_lvl, cond, i) = lo_v;
                        hi(prev_lvl, curr_lvl, cond, i) = hi_v;
                    end
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
            'Name', sprintf('%s %s', fname_stem, p.cond_names{cond}));

        cond_vals = squeeze(vals(:,:,cond,:));
        y_min = min(cond_vals(:), [], 'omitnan');
        y_max = max(cond_vals(:), [], 'omitnan');
        if isnan(y_min) || isnan(y_max)
            y_min = 0; y_max = 1;
        end
        y_pad = max((y_max - y_min) * 0.15, 0.1);
        y_lim = [y_min - y_pad, y_max + y_pad];
        if strcmp(ylabel_str, 'Amplitude (°)')
            y_lim(1) = min(0, y_lim(1));
        elseif strcmp(ylabel_str, 'FWHM (°)')
            y_lim(1) = max(0, y_lim(1));
        end

        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);
                hold on;

                y = squeeze(vals(prev_lvl, curr_lvl, cond, :));
                lo_curr = squeeze(lo(prev_lvl, curr_lvl, cond, :));
                hi_curr = squeeze(hi(prev_lvl, curr_lvl, cond, :));

                plot(n_back_list(:), y, '-', 'Color', line_color, ...
                     'LineWidth', plt_opts.line_width, 'HandleVisibility', 'off');
                scatter(n_back_list(:), y, plt_opts.marker_size_scatter, ...
                        'MarkerFaceColor', line_color, 'MarkerEdgeColor', [1 1 1], ...
                        'MarkerFaceAlpha', 0.9, 'LineWidth', plt_opts.line_width);

                if any(isfinite(lo_curr)) && any(isfinite(hi_curr))
                    yneg = max(0, y - lo_curr);
                    ypos = max(0, hi_curr - y);
                    errorbar(n_back_list(:), y, yneg, ypos, ...
                             'Color', line_color * 0.8, 'CapSize', 0, ...
                             'LineStyle', 'none', 'LineWidth', plt_opts.line_width, ...
                             'HandleVisibility', 'off');
                end

                title(sprintf('%s \\rightarrow %s', level_labels{prev_lvl}, level_labels{curr_lvl}));
                xlabel('n-back');
                if curr_lvl == 1
                    ylabel(ylabel_str);
                end
                xlim([min(n_back_list)-0.5, max(n_back_list)+0.5]);
                ylim(y_lim);
                xticks(n_back_list);
                axis square; box off;
                set(gca, 'TickDir','out', 'LineWidth', plt_opts.line_width);
            end
        end

        sgtitle(sprintf('%s - %s', fname_stem, p.cond_names{cond}));
        if exist('ci_label', 'var') && ~isempty(ci_label)
            annotation(fg, 'textbox', [0.01 0.01 0.4 0.04], ...
                'String', ci_label, 'EdgeColor', 'none', ...
                'FontSize', 8, 'Color', plt_opts.colors.gray);
        end

        if isfield(opts, 'save') && opts.save
            fname = sprintf('%s %s.%s', fname_stem, p.cond_names{cond}, plt_opts.fg_type);
            saveas(fg, fullfile(opts.fig_dir, fname));
        end
        close(fg);
    end
end
