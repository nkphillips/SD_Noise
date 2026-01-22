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
%   plt_opts      - plotting settings; expects fields: colors, line_width,
%                       tick_length
function plotDeltaThetaCount(all_delta_thetas, num, p, plt_opts)
%plotDeltaThetaCount Plots histograms of delta thetas for super subject.

    %% Generate histogram of delta thetas

    %%% Super subject level %%%

    % Define consistent bin edges and find max count across subplots (ignoring windows)
    edges = -90:5:90; % 5-degree bins across [-90, 90]
    max_y = 0;
    for cond = 1:num.conds
        for prev_lvl = 1:num.levels 
            for curr_lvl = 1:num.levels
                counts = histcounts(all_delta_thetas{prev_lvl, curr_lvl, cond}, edges);
                if ~isempty(counts)
                    max_y = max(max_y, max(counts));
                end
            end
        end
    end

    % Plot
    if plt_opts.plot_sup_figures
        for cond = 1:num.conds

            fg_name = ['Super Subj Delta Thetas ' p.cond_names{cond}];

            for prev_lvl = 1:num.levels
                for curr_lvl = 1:num.levels

                    % Create title with actual values
                    if cond == 1
                        % Contrast condition
                        fg_title = [p.contrast{prev_lvl} ' -> ' p.contrast{curr_lvl}];
                    else
                        % Precision condition
                        fg_title = [p.precision{prev_lvl} ' -> ' p.precision{curr_lvl}];
                    end

                    subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);

                    % Plot histogram
                    histogram(all_delta_thetas{prev_lvl, curr_lvl, cond}, 'BinEdges', edges, 'Normalization', 'count');

                    % Format figure
                    axis square;
                    title(fg_title);
                    
                    % Only show x-axis labels on bottom row
                    if curr_lvl == num.levels
                        xlabel('\Delta\theta (°)');
                        xticks(-90:45:90);
                        xtickangle(0);
                    else
                        xticks(-90:45:90);
                        xtickangle(0);
                        set(gca, 'XTickLabel', []);
                    end
                    
                    % Only show y-axis label on bottom left subplot
                    if prev_lvl == 1 && curr_lvl == num.levels
                        ylabel('Count');
                    end
                    
                    set(gca, 'TickDir', 'out', 'TickLength', [plt_opts.tick_length, plt_opts.tick_length]);
                    xlim([-90 90]);
                    ylim([0, max_y]);
                    line([0, 0], [0, max_y], 'LineWidth', plt_opts.line_width, 'Color', plt_opts.colors.black);
                    box off;
                    hold on;

                end
            end

            % Save figure
            if plt_opts.save_sup_figures
                saveas(gcf, fullfile(plt_opts.sup_figure_path, [fg_name '.' plt_opts.fg_type]));
            end

            clf(gcf);

        end
    end

end


