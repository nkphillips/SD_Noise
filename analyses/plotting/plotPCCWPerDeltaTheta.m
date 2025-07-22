function plotPCCWPerDeltaTheta(delta_theta_bin_centers, data, subject_label, num, cond_names, figure_path, save_figures, fg)
% plotPCCWPerDeltaTheta - Plots probability of counterclockwise response as a function of delta theta for each level pair
%
% Inputs:
%   delta_theta_bin_centers - Centers of the delta theta bins
%   data - Data structure containing pCW metrics
%   subject_label - Label for the subject (e.g., 'super' or 'S1')
%   num - Structure containing experiment parameters
%   cond_names - Cell array of condition names
%   figure_path - Path to save figures
%   save_figures - Boolean flag to determine if figures should be saved
%   fg - Figure handle to use for plotting

for cond = 1:num.conds
    % Use the provided figure handle instead of creating a new one
    set(0, 'CurrentFigure', fg);
    clf(fg); % Clear the figure before plotting
    
    figure_name = [subject_label ' pCCW ' cond_names{cond} ' per Level Pair'];

    for prev_lvl = 1:num.levels
        for curr_lvl = 1:num.levels
            subplot(num.levels, num.levels, prev_lvl + (curr_lvl-1)*num.levels);
            
            % Get pCW mean and SEM
            if isfield(data, 'pCW_per_delta_theta_bin_mean') && ~isempty(data.pCW_per_delta_theta_bin_mean)
                pCW_mean = squeeze(data.pCW_per_delta_theta_bin_mean(prev_lvl, curr_lvl, cond, :));
                pCW_sem = squeeze(data.pCW_per_delta_theta_bin_sem(prev_lvl, curr_lvl, cond, :));
            else
                % If data not available, use NaNs
                pCW_mean = nan(size(delta_theta_bin_centers));
                pCW_sem = nan(size(delta_theta_bin_centers));
            end
            
            % Plot 1 - pCW_mean for pCCW
            shadedErrorBar(delta_theta_bin_centers, 1 - pCW_mean, pCW_sem, 'lineprops', '-b');
            
            title([num2str(prev_lvl) ' -> ' num2str(curr_lvl)]);
            if curr_lvl == num.levels, xlabel('\Delta \theta (°)'); end
            if prev_lvl == 1, ylabel('p(Resp|CCW)'); end
            ylim([0, 1]);
            xlim([-180 180]);
            xticks(min(xlim):45:max(xlim));
            line([min(xlim), max(xlim)], [0.5, 0.5], 'LineWidth', 1, 'Color', 'k');
            box off;
            set(gca, 'TickDir', 'out');
        end
    end

    % Save the figure
    if save_figures
        saveas(fg, [figure_path '/' figure_name '.png']); clf
        disp(['Saved ' figure_name '.png']);
    end
    
end
end 