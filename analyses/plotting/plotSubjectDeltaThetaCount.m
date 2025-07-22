function plotSubjectDeltaThetaCount(num, analyzed_data, subj, subj_id, cond_names, figure_path, save_figures, fg)
% plotSubjectDeltaThetaCount - Plots delta theta count histograms for an individual subject
%
% Inputs:
%   num - Structure containing experiment parameters
%   analyzed_data - Data structure containing all analyzed data
%   subj - Subject index
%   subj_id - Subject ID string
%   cond_names - Cell array of condition names
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

% Use the provided figure handle instead of creating a new one
set(0, 'CurrentFigure', fg);
clf(fg); % Clear the figure before plotting

figure_name = ['S' subj_id ' Delta Theta Counts'];

% Initialize bins for histograms
binEdges = -180:10:180;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;

% Create subplots for each level pair and condition
for cond = 1:num.conds
    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            % Calculate subplot index
            subplot_idx = (prev_lvl - 1) * num.levels * num.conds + (curr_lvl - 1) * num.conds + cond;
            subplot(num.levels*num.levels, num.conds, subplot_idx);
            
            % Initialize counts array
            counts_array = [];
            
            % Loop through all runs for this subject
            for n_run = 1:length(analyzed_data.ind(subj).counts)
                if ~isempty(analyzed_data.ind(subj).counts(n_run).delta_theta) && ...
                   prev_lvl <= size(analyzed_data.ind(subj).counts(n_run).delta_theta, 1) && ...
                   curr_lvl <= size(analyzed_data.ind(subj).counts(n_run).delta_theta, 2) && ...
                   cond <= size(analyzed_data.ind(subj).counts(n_run).delta_theta, 3)
                    
                    % Get delta theta values for this level pair and condition
                    delta_theta_vals = analyzed_data.ind(subj).counts(n_run).delta_theta{prev_lvl, curr_lvl, cond};
                    if ~isempty(delta_theta_vals)
                        counts_array = [counts_array, delta_theta_vals];
                    end
                end
            end
            
            % Plot histogram if we have data
            if ~isempty(counts_array)
                histogram(counts_array, binEdges, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
                
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
                set(gca, 'TickDir', 'out');
                
                % Display the total count
                text(0.1, 0.9, ['n = ' num2str(length(counts_array))], 'Units', 'normalized');
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