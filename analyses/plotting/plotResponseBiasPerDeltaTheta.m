function plotResponseBiasPerDeltaTheta(delta_theta_bin_centers, data, subject_label, num, cond_names, p, green, figure_path, save_figures, fg)
% plotResponseBiasPerDeltaTheta - Plots response bias as a function of delta theta for each level pair
%
% Inputs:
%   delta_theta_bin_centers - Centers of the delta theta bins
%   data - Data structure containing response bias metrics
%   subject_label - Label for the subject (e.g., 'super' or 'S1')
%   num - Structure containing experiment parameters
%   cond_names - Cell array of condition names
%   p - Structure containing parameters
%   green - Color for mu plot
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

for cond = 1:num.conds
    % Use the provided figure handle instead of creating a new one
    set(0, 'CurrentFigure', fg);
    clf(fg); % Clear the figure before plotting
    
    figure_name = [subject_label ' Response Bias ' cond_names{cond} ' per Level Pair'];

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);
            
            % Check if using super subject data or individual subject data
            if strcmp(subject_label, 'super')
                mu_values = squeeze(data.mu_per_delta_theta_bin(prev_lvl, curr_lvl, cond, :));
                
                % Plot mu values
                plot(delta_theta_bin_centers, mu_values, 'LineWidth', 1, 'Color', green);
                hold on;
                
                % Plot serial dependence estimated bias if available
                if isfield(data, 'serial_dependence_estimated_bias') && ~isempty(data.serial_dependence_estimated_bias)
                    serial_dep_bias = squeeze(data.serial_dependence_estimated_bias(prev_lvl, curr_lvl, cond, :));
                    plot(-180:180, serial_dep_bias, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-');
                    
                    % Add R² to title if available
                    if isfield(data, 'serial_dependence_r2') && ~isempty(data.serial_dependence_r2)
                        r2 = data.serial_dependence_r2(prev_lvl, curr_lvl, cond);
                        title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl) ' (R^2 = ' num2str(round(r2, 2)) ')'];
                    else
                        title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                    end
                else
                    title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                end
            else
                % For individual subjects
                if isfield(data, 'mu_per_delta_theta') && ~isempty(data.mu_per_delta_theta)
                    mu_values = squeeze(data.mu_per_delta_theta(prev_lvl, curr_lvl, cond, :));
                    
                    % Plot mu values
                    plot(delta_theta_bin_centers, mu_values, 'LineWidth', 1, 'Color', green);
                    hold on;
                    
                    % Plot serial dependence estimated bias if available
                    if isfield(data, 'serial_dependence_estimated_bias') && ~isempty(data.serial_dependence_estimated_bias) && ...
                       prev_lvl <= size(data.serial_dependence_estimated_bias, 1) && ...
                       curr_lvl <= size(data.serial_dependence_estimated_bias, 2) && ...
                       cond <= size(data.serial_dependence_estimated_bias, 3)
                        serial_dep_bias = squeeze(data.serial_dependence_estimated_bias(prev_lvl, curr_lvl, cond, :));
                        if ~isempty(serial_dep_bias) && ~all(isnan(serial_dep_bias))
                            plot(-180:180, serial_dep_bias, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-');
                        end
                        
                        % Add R² to title if available
                        if isfield(data, 'serial_dependence_r2') && ~isempty(data.serial_dependence_r2) && ...
                           prev_lvl <= size(data.serial_dependence_r2, 1) && ...
                           curr_lvl <= size(data.serial_dependence_r2, 2) && ...
                           cond <= size(data.serial_dependence_r2, 3)
                            r2 = data.serial_dependence_r2(prev_lvl, curr_lvl, cond);
                            if ~isnan(r2)
                                title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl) ' (R^2 = ' num2str(round(r2, 2)) ')'];
                            else
                                title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                            end
                        else
                            title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                        end
                    else
                        title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                    end
                else
                    mu_values = nan(size(delta_theta_bin_centers));
                    title_str = [num2str(prev_lvl) ' -> ' num2str(curr_lvl)];
                end
            end
            
            title(title_str);
            if curr_lvl == num.levels, xlabel('\Delta \theta (°)'); end
            if prev_lvl == 1, ylabel('\mu'); end
            ylim([p.response_bias_bounds(2,1) p.response_bias_bounds(1,1)]);
            xlim([-180 180]);
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