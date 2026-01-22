% plotPerformance
% Plots percent correct and percent CCW as a function of delta_theta arranged as a
% num.levels x num.levels grid (rows = previous level, cols = current level), for
% each condition (Contrast, Precision). Uses pre-computed performance arrays.

function plotPerformance(delta_theta_centers, performance_array, pCW_array, p, plt_opts, fg_prefix, perf_ci, save_path, save_flag)


metrics = { 'Percent Correct', 'Percent CCW' };
num_levels = size(performance_array, 1);
num_conds = size(performance_array, 3);

have_ci = (nargin >= 7) && ~isempty(perf_ci);

% For each metric, compute a global y-limit across ALL conditions and subplots,
% then plot each condition using that shared y-limit. This allows direct
% comparison/overlay between conditions.
for m_idx = 1:numel(metrics)
    metric_name = metrics{m_idx};

    % Compute global y-range across conditions for this metric
    y_min_global = inf;
    y_max_global = -inf;
    for cond_tmp = 1:num_conds
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                y_vals_tmp = nan(size(delta_theta_centers));
                for i_window = 1:numel(delta_theta_centers)
                    % Use pre-computed performance values
                    switch m_idx
                        case 1 % Percent Correct
                            y_vals_tmp(i_window) = performance_array(prev_lvl, curr_lvl, cond_tmp, i_window);
                        case 2 % Percent CCW
                            y_vals_tmp(i_window) = pCW_array(prev_lvl, curr_lvl, cond_tmp, i_window);
                    end
                end
                finite_y = y_vals_tmp(~isnan(y_vals_tmp));
                if ~isempty(finite_y)
                    y_min_global = min(y_min_global, min(finite_y));
                    y_max_global = max(y_max_global, max(finite_y));
                end
            end
        end
    end

    % Finalize global y-limits with padding and clamping
    if isfinite(y_min_global) && isfinite(y_max_global)
        y_min = max(0, y_min_global - 10);
        y_max = min(100, y_max_global + 10);
        if y_max - y_min < 20
            pad = (20 - (y_max - y_min)) / 2;
            y_min = max(0, y_min - pad);
            y_max = min(100, y_max + pad);
        end
        y_lim_metric = [y_min, y_max];
    else
        y_lim_metric = [0, 100];
    end

    % Now plot per condition using the shared metric y-limits
    for cond = 1:num_conds
        fg_name = [fg_prefix ' ' metric_name ' ' p.cond_names{cond}];

        clf(gcf);
        set(gcf, 'Color', 'w');

        % Columns = current level; Rows = previous level
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels

                subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);

                y_vals = nan(size(delta_theta_centers));
                for i_window = 1:numel(delta_theta_centers)
                    % Use pre-computed performance values
                    switch m_idx
                        case 1 % Percent Correct
                            y_vals(i_window) = performance_array(prev_lvl, curr_lvl, cond, i_window);
                        case 2 % Percent CCW
                            y_vals(i_window) = pCW_array(prev_lvl, curr_lvl, cond, i_window);
                    end
                end

                % Choose color by condition
                if cond == 1
                    plot_color = plt_opts.colors.blue;
                else
                    plot_color = plt_opts.colors.green;
                end

                if have_ci
                    % Select corresponding CI arrays
                    pc_lo = []; pc_hi = []; pccw_lo = []; pccw_hi = [];
                    if m_idx == 1 && isfield(perf_ci, 'pc_lo') && isfield(perf_ci, 'pc_hi')
                        pc_lo = squeeze(perf_ci.pc_lo(prev_lvl, curr_lvl, cond, :));
                        pc_hi = squeeze(perf_ci.pc_hi(prev_lvl, curr_lvl, cond, :));
                    elseif m_idx == 2 && isfield(perf_ci, 'pccw_lo') && isfield(perf_ci, 'pccw_hi')
                        pccw_lo = squeeze(perf_ci.pccw_lo(prev_lvl, curr_lvl, cond, :));
                        pccw_hi = squeeze(perf_ci.pccw_hi(prev_lvl, curr_lvl, cond, :));
                    end

                    if m_idx == 1 && ~isempty(pc_lo)
                        err_upper = pc_hi(:)' - y_vals(:)';
                        err_lower = y_vals(:)' - pc_lo(:)';
                        err = [err_upper; err_lower];
                        shadedErrorBar(delta_theta_centers(:)', y_vals(:)', err, 'lineProps', {'-','Color', plot_color, 'LineWidth', plt_opts.line_width}, 'transparent', true, 'patchSaturation', 0.5);
                    elseif m_idx == 2 && ~isempty(pccw_lo)
                        err_upper = pccw_hi(:)' - y_vals(:)';
                        err_lower = y_vals(:)' - pccw_lo(:)';
                        err = [err_upper; err_lower];
                        shadedErrorBar(delta_theta_centers(:)', y_vals(:)', err, 'lineProps', {'-','Color', plot_color, 'LineWidth', plt_opts.line_width}, 'transparent', true, 'patchSaturation', 0.5);
                    else
                        plot(delta_theta_centers, y_vals, 'LineWidth', plt_opts.line_width, 'Color', plot_color);
                    end
                else
                    plot(delta_theta_centers, y_vals, 'LineWidth', plt_opts.line_width, 'Color', plot_color);
                end

                % Title with actual values
                if cond == 1
                    fg_title = [p.contrast{prev_lvl} ' -> ' p.contrast{curr_lvl}];
                else
                    fg_title = [p.precision{prev_lvl} ' -> ' p.precision{curr_lvl}];
                end
                title(fg_title);

                axis square;
                xlim([-90 90]);
                ylim(y_lim_metric);

                % Reference lines
                xl = xlim;
                line([0 0], y_lim_metric, 'Color', 'k', 'LineWidth', 1);
                line(xl, [50 50], 'Color', 'k', 'LineWidth', 1, 'LineStyle', '--');

                % Only show x tick labels on bottom row
                if prev_lvl == num_levels
                    xlabel('\Delta\theta (°)');
                    xticks(-90:45:90);
                    xtickangle(0);
                else
                    xticks(-90:45:90);
                    xtickangle(0);
                    set(gca, 'XTickLabel', []);
                end

                % Only show y label on left column
                if curr_lvl == 1
                    if m_idx == 1
                        ylabel('% Correct');
                    else
                        ylabel('% CCW');
                    end
                end

                set(gca, 'TickDir', 'out', 'TickLength', [plt_opts.tick_length, plt_opts.tick_length]);
                box off;
                hold on;
            end
        end

        if save_flag
            saveas(gcf, fullfile(save_path, [fg_name '.' plt_opts.fg_type]));
        end
    end
end

end


