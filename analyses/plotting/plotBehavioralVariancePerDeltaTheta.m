function plotBehavioralVariancePerDeltaTheta(delta_theta_bin_centers, data, subject_label, num, cond_names, p, figure_path, save_figures, fg)
% plotBehavioralVariancePerDeltaTheta - Plots behavioral variance as a function of delta theta for each level pair
%
% Inputs:
%   delta_theta_bin_centers - Centers of the delta theta bins
%   data - Data structure containing behavioral variance metrics
%   subject_label - Label for the subject (e.g., 'super' or 'S1')
%   num - Structure containing experiment parameters
%   cond_names - Cell array of condition names
%   p - Structure containing parameters
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

for cond = 1:num.conds
    % Use the provided figure handle instead of creating a new one
    set(0, 'CurrentFigure', fg);
    clf(fg); % Clear the figure before plotting
    
    figure_name = [subject_label ' Behavioral Variance ' cond_names{cond} ' per Level Pair'];

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);
            
            % Check if using super subject data or individual subject data
            if strcmp(subject_label, 'super')
                sigma_values = squeeze(data.sigma_per_delta_theta_bin(prev_lvl, curr_lvl, cond, :));
            else
                % For individual subjects
                if isfield(data, 'sigma_per_delta_theta') && ~isempty(data.sigma_per_delta_theta)
                    sigma_values = squeeze(data.sigma_per_delta_theta(prev_lvl, curr_lvl, cond, :));
                else
                    sigma_values = nan(size(delta_theta_bin_centers));
                end
            end
            
            % Plot sigma values
            plot(delta_theta_bin_centers, sigma_values, 'LineWidth', 1, 'Color', 'k');
            
            title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
            if curr_lvl == num.levels, xlabel('\Delta \theta (°)'); end
            if prev_lvl == 1, ylabel('\sigma'); end
            ylim([0 p.response_bias_bounds(1,2)]);
            xlim([-90 90]);
            xticks(min(xlim):45:max(xlim));
            line([min(xlim), max(xlim)], [0, 0], 'LineWidth', 1, 'Color', 'k');
            box off;
            set(gca, 'TickDir', 'out');
        end
    end

    % Save the figure if requested
    if save_figures
        saveas(fg, [figure_path '/' figure_name '.png']); clf;
        disp(['Saved ' figure_name '.png']);
    end
    
end
end 