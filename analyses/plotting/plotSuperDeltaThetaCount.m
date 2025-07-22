function plotSuperDeltaThetaCount(num, analyzed_data, cond_names, figure_path, save_figures, fg)
% plotSuperDeltaThetaCount - Plots delta theta count histograms for the super subject
%
% Inputs:
%   num - Structure containing experiment parameters
%   analyzed_data - Data structure containing all analyzed data
%   cond_names - Cell array of condition names
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

% Use the provided figure handle instead of creating a new one
set(0, 'CurrentFigure', fg);
clf(fg); % Clear the figure before plotting

figure_name = 'Super Delta Theta Counts';

% Initialize bins for histograms
binEdges = -180:10:180;

% Create subplots for each level pair and condition
for cond = 1:num.conds
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            % Calculate subplot index
            subplot_idx = (prev_lvl - 1) * num.levels * num.conds + (curr_lvl - 1) * num.conds + cond;
            subplot(num.levels*num.levels, num.conds, subplot_idx);
            
            % Plot histogram if we have data
            if ~isempty(analyzed_data.super.all_delta_thetas{prev_lvl, curr_lvl, cond})
                histogram(analyzed_data.super.all_delta_thetas{prev_lvl, curr_lvl, cond}, binEdges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
                
                % Add title for the first row only
                if prev_lvl == 1 && curr_lvl == 1
                    title(cond_names{cond});
                end
                
                % Add labels for the rows and columns
                if cond == 1
                    ylabel(['L' num2str(prev_lvl) ' \rightarrow L' num2str(curr_lvl)]);
                end
                
                % Add x-axis label for the bottom row
                if prev_lvl == num.levels && curr_lvl == num.levels
                    xlabel('\Delta \theta (°)');
                end
                
                % Format axes
                xlim([-180 180]);
                xticks(-180:90:180);
                box off;
                
                text(0.1, 0.9, ['n = ' num2str(length(analyzed_data.super.all_delta_thetas{prev_lvl, curr_lvl, cond}))], 'Units', 'normalized');
            else
                % No data available
                text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
                axis off;
            end
        end
    end
end

% Save the figure if requested
if save_figures
    saveas(fg, [figure_path '/' figure_name '.png']); clf;   
    disp(['Saved ' figure_name '.png']);
end

end 