function plotSerialDependenceAmplitude(data, subject_label, num, cond_names, p, figure_path, save_figures, fg)
% plotSerialDependenceAmplitude - Plots serial dependence amplitude for each level pair
%
% Inputs:
%   data - Data structure containing serial dependence parameters
%   subject_label - Label for the subject (e.g., 'super' or 'S1')
%   num - Structure containing experiment parameters
%   cond_names - Cell array of condition names
%   p - Structure containing parameters
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

% Use the provided figure handle instead of creating a new one
set(0, 'CurrentFigure', fg);
clf(fg); % Clear the figure before plotting

figure_name = [subject_label ' Serial Dependence Amplitude per Level Pair'];

for cond = 1:num.conds
    subplot(1, num.conds, cond);
    
    % Check if data has the required field
    if isfield(data, 'serial_dependence_params_est') && ~isempty(data.serial_dependence_params_est)
        % Extract amplitude values (first parameter of serial dependence)
        amplitudes_matrix = squeeze(data.serial_dependence_params_est(:, :, cond, 1))';
        
        % Plot as a bar graph
        bar(amplitudes_matrix);
        
        % Format plot
        title(cond_names{cond});
        xlabel('Current Level');
        ylabel('Amplitude');
        xticks(1:num.levels);
        
        % Set y-axis limits
        if isfield(p, 'serial_dependence_bounds') && ~isempty(p.serial_dependence_bounds)
            ylim([0 p.serial_dependence_bounds(1,1)]);
        else
            ylim([0 10]); % Default if bounds not available
        end
        
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