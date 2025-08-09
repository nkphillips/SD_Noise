%% plotSerialDependence   
% Bar plots of serial dependence estimates (amplitude and width)
% Each subplot in the figure is the noise type, bar group is the previous level, while each bar is the current level. 

% Inputs:
% p: parameters
% plt_settings: plotting settings
% sd: serial dependence estimates

% Outputs:
% fg: figure handle

function plotSerialDependence(sd_data, param_index, param_name, p, plt_settings)
    % plotSerialDependence - Create heatmap of serial dependence parameters
    %
    % This function creates a heatmap visualization of serial dependence parameters
    % across all experimental conditions (3×3×2 structure).
    %
    % Inputs:
    %   sd_data - 4D array [prev_lvl, curr_lvl, cond, param] from sd.all.params_est
    %   param_index - which parameter to plot (1=amplitude, 2=width, 3=baseline, 4=sigma)
    %   param_name - string name of parameter for title
    %   p - parameter struct
    %   plt_settings - plotting settings struct
    
    % Extract the parameter data
    param_data = squeeze(sd_data(:,:,:,param_index)); % [prev_lvl, curr_lvl, cond]

    % If plotting width (param_index==2), convert from model parameter w (1/deg)
    % to Gaussian FWHM in degrees: FWHM = (2 * sqrt(ln 2)) / w ≈ 1.6651 / w
    if param_index == 2
        with_epsilon = max(param_data, eps);
        param_data = (2 * sqrt(log(2))) ./ with_epsilon;
        % Clarify label for width
        param_name = [param_name ' (FWHM, °)'];
    end
    
    % Create figure with subplots for each condition
    num_conds = size(param_data, 3);
    
    % Calculate global min and max for consistent colorbar range
    global_min = min(param_data(:));
    global_max = max(param_data(:));
    
    for cond = 1:num_conds
        subplot(1, num_conds, cond);
        
        % Get data for current condition
        cond_data = param_data(:,:,cond);
        
        % Create heatmap with consistent range
        imagesc(cond_data, [global_min, global_max]);

        % Condition-specific colormap: low values = dark hue, high values = light of same hue
        if cond == 1
            base_color = plt_settings.colors.blue; % contrast
        else
            base_color = plt_settings.colors.green; % precision
        end
        start_color = 0.35 * base_color; % dark tone
        white_color = [1 1 1];
        steps = 256;
        cm = [linspace(start_color(1), white_color(1), steps)', ...
              linspace(start_color(2), white_color(2), steps)', ...
              linspace(start_color(3), white_color(3), steps)'];
        colormap(gca, cm);
        
        % Add colorbar
        colorbar;
        
        % Set axis labels and ticks
        xlabel('Current Level');
        ylabel('Previous Level');
        
        % Use actual contrast/precision values for tick labels
        if cond == 1
            % Contrast condition
            set(gca, 'XTick', 1:3, 'XTickLabel', p.contrast);
            set(gca, 'YTick', 1:3, 'YTickLabel', p.contrast);
        else
            % Precision condition
            set(gca, 'XTick', 1:3, 'XTickLabel', p.precision);
            set(gca, 'YTick', 1:3, 'YTickLabel', p.precision);
        end
        
        % Add value annotations
        for i = 1:3
            for j = 1:3
                text(j, i, sprintf('%.2f', cond_data(i,j)), ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Color', plt_settings.colors.gray - 0.1, 'FontWeight', 'bold', ...
                    'FontSize', 10);
            end
        end
        
        % Set title
        if cond == 1
            title('Contrast');
        else
            title('Precision');
        end
        
        % Set axis properties
        axis square;
        set(gca, 'TickDir', 'out', 'LineWidth', plt_settings.line_width, 'Box', 'off');
    end
    
end