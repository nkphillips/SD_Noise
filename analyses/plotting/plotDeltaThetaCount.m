function plotDeltaThetaCount(delta_thetas_cell, p, plt_settings)
% plotDeltaThetaCount - Plot delta-theta count histograms for any analysis level
%
% This function accepts either:
%   - A 4D cell array of delta-thetas per window: {prev_lvl, curr_lvl, cond, window}
%   - A 3D cell array already aggregated over windows: {prev_lvl, curr_lvl, cond}
%
% And produces, for each condition, a num.levels x num.levels grid of histograms
% showing the distribution of delta-theta values for each (prev->curr) level pair.
%
% Inputs:
%   delta_thetas_cell - 3D or 4D cell of delta-thetas, as described above
%   p                 - parameter struct; expects fields: num.levels, num.conds,
%                       contrast, precision, cond_names
%   plt_settings      - plotting settings; expects fields: colors, line_width,
%                       tick_length

    % Determine dimensions and create a function to fetch concatenated delta-thetas per pair
    dims = ndims(delta_thetas_cell);
    if dims == 4
        % {prev, curr, cond, window}: concatenate across windows
        get_pair_deltas = @(i,j,c) vertcat(delta_thetas_cell{i,j,c,:});
    elseif dims == 3
        % {prev, curr, cond}: already aggregated
        get_pair_deltas = @(i,j,c) delta_thetas_cell{i,j,c};
    else
        error('delta_thetas_cell must be a 3D or 4D cell array');
    end

    % Infer dimensions from the provided cell
    num_levels = size(delta_thetas_cell, 1);
    num_conds = size(delta_thetas_cell, 3);

    % Define consistent bin edges across [-90, 90]
    edges = -90:5:90;

    % Find global max count across all subplots and conditions to unify y-limits
    max_y = 0;
    for cond = 1:num_conds
        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels
                deltas = get_pair_deltas(prev_lvl, curr_lvl, cond);
                if ~isempty(deltas)
                    counts = histcounts(deltas, edges);
                    if ~isempty(counts)
                        max_y = max(max_y, max(counts));
                    end
                end
            end
        end
    end
    if max_y == 0
        max_y = 1; % avoid zero-height axes if no data
    end

    % Plot per condition
    for cond = 1:num_conds
        clf(gcf);
        set(gcf, 'Color', 'w');

        for prev_lvl = 1:num_levels
            for curr_lvl = 1:num_levels

                subplot(num_levels, num_levels, curr_lvl + (prev_lvl-1)*num_levels);

                deltas = get_pair_deltas(prev_lvl, curr_lvl, cond);

                if ~isempty(deltas)
                    histogram(deltas, 'BinEdges', edges, 'Normalization', 'count', ...
                        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');

                    % Formatting
                    axis square;
                    xlim([-90 90]);
                    ylim([0 max_y]);
                    xticks([-90 -45 0 45 90]);
                    line([0 0], [0 max_y], 'LineWidth', plt_settings.line_width, 'Color', [0 0 0]);
                    set(gca, 'TickDir', 'out', 'TickLength', [plt_settings.tick_length, plt_settings.tick_length], 'Box', 'off');

                    % Bottom row x-labels only
                    if prev_lvl == num_levels
                        xlabel('\Delta\theta (°)');
                    else
                        set(gca, 'XTickLabel', []);
                    end

                    % Left column y-label only
                    if curr_lvl == 1
                        ylabel('Count');
                    end

                    hold on;
                else
                    % No data for this cell
                    axis off;
                end
            end
        end

    end
end


