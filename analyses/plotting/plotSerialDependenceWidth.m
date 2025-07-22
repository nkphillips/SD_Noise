function plotSerialDependenceWidth(data, subject_label, num, cond_names, figure_path, save_figures, fg)
% plotSerialDependenceWidth - Plots serial dependence width for each level pair
%
% Inputs:
%   data - Data structure containing serial dependence parameters
%   subject_label - Label for the subject (e.g., 'super' or 'S1')
%   num - Structure containing experiment parameters
%   cond_names - Cell array of condition names
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

% Use the provided figure handle instead of creating a new one
set(0, 'CurrentFigure', fg);
clf(fg); % Clear the figure before plotting

figure_name = [subject_label ' Serial Dependence Width per Level Pair'];

for cond = 1:num.conds
    subplot(1, num.conds, cond);
    
    % Check if data has the required field
    if isfield(data, 'serial_dependence_params_est') && ~isempty(data.serial_dependence_params_est)
        % Extract width values (second parameter of serial dependence)
        widths_matrix = squeeze(data.serial_dependence_params_est(:, :, cond, 2))';
        
        % Plot as a bar graph
        bar(widths_matrix);
        
        % Format plot
        title(cond_names{cond});
        xlabel('Current Level');
        ylabel('Width');
        xticks(1:num.levels);
        
        % Set y-axis limits - typically width is around 0-150 degrees
        ylim([0 150]);
        
        box off; 
        axis square;
        set(gca, 'TickDir', 'out');
    else
        % Display a message if no data available
        text(0.5, 0.5, 'No serial dependence data available', 'HorizontalAlignment', 'center');
        axis off;
    end
end

% Save the figure if requested
if save_figures
    saveas(fg, [figure_path '/' figure_name '.png']); clf;
    disp(['Saved ' figure_name '.png']);
end

end 