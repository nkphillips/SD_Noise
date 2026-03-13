function plot_frame_timing_stats(exe_timing, p, dirs, t, ps)
% plot_frame_timing_stats generates summary statistic figures of frame timing
% vs intended timing and saves them as a PDF.

% Set figure path and check existence
figure_path = fullfile(dirs.data_dir, p.subj_ID);
if ~exist(figure_path, 'dir')
    mkdir(figure_path)
    disp([figure_path ' created.'])
end

% Open figure
fg = figure('Color', [1 1 1], 'Position', [100 100 1200 600]);
set(0, 'CurrentFigure', fg)

if p.training
    exp_phase = 'Training_';
elseif p.calibration
    exp_phase = 'Calibration_';
else
    exp_phase = '';
end
figure_name = ['SD_Noise_Exp2_' exp_phase 'S' p.subj_ID '_Run' num2str(p.run_num) '_' p.display_setup '_FrameTiming'];

sgtitle(figure_name, 'Interpreter', 'none');

% Set paper size to match figure dimensions to prevent clipping
set(fg, 'PaperPositionMode', 'auto');
set(fg, 'PaperOrientation', 'landscape');

phases = {'test', 'mask', 'delay', 'probe'};

for i = 1:length(phases)
    phase = phases{i};

    % --- Row 1: Missed Deadlines ---
    missed_field = [phase '_Missed'];
    if isfield(exe_timing, missed_field) && ~isempty(exe_timing.(missed_field))
        missed_data = cat(1, exe_timing.(missed_field){:});
        frames_missed = sum(missed_data(:) > 0);
        pct_missed = round(mean(missed_data(:) > 0) * 100);

        subplot(2, 4, i)
        plot(missed_data', '.', 'MarkerSize', 6); hold on;
        line([0 size(missed_data,2)+1], [0 0], 'Color', 'k', 'LineStyle', '-')

        title_str = [capitalize(phase) ' missed: ' num2str(frames_missed) ' (' num2str(pct_missed) '%)'];
        title(title_str)
        box off; axis square;
        xlabel('Frame #')
        ylabel('Deadline offset (s)')
    end

    % --- Row 2: Actual Frame Durations ---
    stim_onset_field = [phase '_StimulusOnsetTime'];
    if isfield(exe_timing, stim_onset_field) && ~isempty(exe_timing.(stim_onset_field))
        stim_onset_data = cat(1, exe_timing.(stim_onset_field){:});

        % Calculate diff along frame dimension (dim 2)
        frame_durs = diff(stim_onset_data, 1, 2);

        subplot(2, 4, i + 4)
        if ~isempty(frame_durs) && size(frame_durs, 2) > 0
            histogram(frame_durs(:), 50); hold on;
            yl = ylim;
            line([t.frame_dur t.frame_dur], yl, 'Color', ps.colors.red, 'LineWidth', 1, 'LineStyle', '--');

            title([capitalize(phase) ' Frame Durs (s)'])
            box off; axis square;
            xlabel('Duration (s)')
            ylabel('Count')
        else
            title([capitalize(phase) ' Frame Durs (N/A)'])
            axis off;
        end
    end
end

% Save and close
set(fg, 'PaperUnits', 'inches');
set(fg, 'PaperSize', [12 6]);
set(fg, 'PaperPosition', [0 0 12 6]);
saveas(gcf, fullfile(figure_path, [figure_name '.pdf']));
close(fg)

end

function str = capitalize(str)
if ~isempty(str)
    str(1) = upper(str(1));
end
end