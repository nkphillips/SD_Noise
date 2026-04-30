function plotSerialDependenceSubjectSummary(sd, param_index, param_name, p, plt_opts, fg, error_mode)
% plotSerialDependenceSubjectSummary
%
% Scatter plots of subject-level serial-dependence parameter summaries.
% Central points are across-subject means from sd.ind, with either SEM or
% t-based 95% CIs. For width, each subject's w is converted to FWHM before
% computing the mean/error bars.
%
% Inputs:
%   sd          - serial-dependence struct with populated sd.ind
%   param_index - 1=amplitude, 2=width, 3=baseline
%   param_name  - display name
%   p           - label/config struct
%   plt_opts    - plotting options
%   fg          - figure handle
%   error_mode  - 'sem' or 't95'

    if nargin < 7 || isempty(error_mode)
        error_mode = 'sem';
    end

    num_subjects = numel(sd.ind);
    num_conds = numel(p.cond_names);
    num_levels = numel(p.contrast);

    subj_values = nan(num_subjects, num_levels, num_levels, num_conds);
    for subj = 1:num_subjects
        raw = squeeze(sd.ind(subj).params_est(:, :, :, param_index));
        if param_index == 2
            raw = (2 * sqrt(log(2))) ./ max(raw, eps);
        end
        subj_values(subj, :, :, :) = raw;
    end

    mean_data = squeeze(mean(subj_values, 1, 'omitnan'));
    sem_data = squeeze(std(subj_values, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(subj_values), 1)));
    n_data = squeeze(sum(isfinite(subj_values), 1));

    if param_index == 2
        param_name = [param_name ' (FWHM, °)'];
    end

    switch lower(error_mode)
        case 'sem'
            err_data = sem_data;
            error_label = 'mean ± SEM';
            file_label = 'Subject Mean SEM';
        case 't95'
            tcrit = tinv(0.975, max(n_data - 1, 1));
            tcrit(n_data <= 1) = NaN;
            err_data = sem_data .* tcrit;
            error_label = 'mean ± t_{12} 95% CI';
            file_label = 'Subject Mean t95';
        otherwise
            error('Unknown error_mode: %s', error_mode);
    end

    clf(fg);

    for cond = 1:num_conds
        subplot(1, num_conds, cond);
        hold on;

        cond_data = squeeze(mean_data(:, :, cond));
        cond_err = squeeze(err_data(:, :, cond));

        if cond == 1
            base_color = plt_opts.colors.blue;
            legend_vals = p.contrast;
            x_labels = fliplr(p.contrast);
        else
            base_color = plt_opts.colors.green;
            legend_vals = p.precision;
            x_labels = fliplr(p.precision);
        end
        marker_colors = repmat(base_color, size(cond_data, 1), 1) .* [1 0.70 0.25]';

        x = 1:num_levels;
        y = fliplr(cond_data)';
        err = fliplr(cond_err)';

        for i = 1:size(y, 2)
            plot(x, y(:, i), '-', 'Color', marker_colors(i, :), ...
                'LineWidth', plt_opts.line_width, 'HandleVisibility', 'off');
            scatter(x, y(:, i), 50, ...
                'MarkerFaceColor', marker_colors(i, :), ...
                'MarkerEdgeColor', [1 1 1], ...
                'MarkerFaceAlpha', 0.75, ...
                'LineWidth', plt_opts.line_width);
            errorbar(x, y(:, i), err(:, i), err(:, i), ...
                'Color', marker_colors(i, :) * 0.8, ...
                'CapSize', 0, ...
                'LineStyle', 'none', ...
                'LineWidth', plt_opts.line_width, ...
                'HandleVisibility', 'off');
        end

        title(p.cond_names{cond});
        set(gca, 'XTick', 1:num_levels, 'XTickLabel', x_labels);
        xlabel('Current level');
        ylabel(param_name);
        text(0.02, 0.02, error_label, 'Units', 'normalized', ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
            'FontSize', 7, 'Color', plt_opts.colors.gray, ...
            'Interpreter', 'tex');

        if strcmp(param_name, 'Amplitude (°)')
            ylim([0 10]);
        end
        axis square;
        set(gca, 'TickDir', 'out', 'LineWidth', plt_opts.line_width, 'Box', 'off');
        legend(legend_vals);
    end

    if plt_opts.save_sup_figures
        fg_name = ['Super Subj SD ' file_label ' ' param_name ' Scatter'];
        saveas(fg, fullfile(plt_opts.sup_figure_path, [fg_name '.' plt_opts.fg_type]));
    end
    clf(fg);
end
