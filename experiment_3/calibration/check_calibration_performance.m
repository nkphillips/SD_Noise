%%% check_calibration_performance
% Independently checks raw calibration performance against the physical
% contrast and orientation-filter-width values used on each trial.
% This script does not fit a psychometric function or modify saved data.

clear;
close all;
clc;

%% Settings

subj_IDs = {'001'};
exclude_zero_offset_trials = true;
connect_points = true;

script_dir = fileparts(mfilename('fullpath'));
data_base_dir = fullfile(script_dir, '..', 'data');

%% Load and validate calibration trials

for i_subj = 1:numel(subj_IDs)

    subj_ID = subj_IDs{i_subj};
    data_dir = fullfile(data_base_dir, subj_ID);

    file_patterns = { ...
        fullfile(data_dir, ['SD_Noise_Exp3_Calibration_S' subj_ID '_Run*.mat']), ...
        fullfile(data_dir, ['SD_Noise_Exp2_Calibration_S' subj_ID '_Run*.mat'])};

    files = [];
    for i_pattern = 1:numel(file_patterns)
        files = [files; dir(file_patterns{i_pattern})]; %#ok<AGROW>
    end

    if isempty(files)
        warning('No calibration runs found for Subject %s.', subj_ID);
        continue;
    end

    % Remove duplicate paths if both filename patterns find the same file.
    full_paths = fullfile({files.folder}, {files.name});
    [~, unique_idx] = unique(full_paths, 'stable');
    files = files(unique_idx);

    feature = strings(0, 1);
    stimulus_value = zeros(0, 1);
    correct = zeros(0, 1);
    run_number = zeros(0, 1);
    block_number = zeros(0, 1);
    level_index = zeros(0, 1);
    abs_probe_offset = zeros(0, 1);

    n_scoring_mismatches = 0;
    n_scored_trials = 0;

    for i_file = 1:numel(files)

        loaded = load(fullfile(files(i_file).folder, files(i_file).name), 'run_info');
        run_info = loaded.run_info;
        p = run_info.p;
        behav_data = run_info.behav_data;

        if isfield(run_info, 'toggles') && run_info.toggles.simulate_response
            warning('%s used simulated responses.', files(i_file).name);
        end

        for n_block = 1:p.num_blocks

            curr_feature = p.feature_order(n_block);
            test_orientation = p.trial_events(:, 1, n_block);
            probe_orientation = p.trial_events(:, 2, n_block);
            curr_level_index = p.trial_events(:, 3, n_block);
            curr_response = behav_data.response(:, n_block);
            saved_correct = behav_data.correct(:, n_block);

            % Axial signed difference: probe minus test, wrapped to [-90, 90].
            probe_offset = probe_orientation - test_orientation;
            probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
            probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;

            % Reconstruct the response scoring used during the experiment:
            % response 1 = CCW (negative offset), response 2 = CW.
            expected_response = 2 * ones(size(probe_offset));
            expected_response(probe_offset < 0) = 1;
            recomputed_correct = expected_response == curr_response;

            scored = ~isnan(curr_response) & ~isnan(saved_correct);
            n_scored_trials = n_scored_trials + nnz(scored);
            n_scoring_mismatches = n_scoring_mismatches + ...
                nnz(saved_correct(scored) ~= recomputed_correct(scored));

            valid = scored;
            if exclude_zero_offset_trials
                valid = valid & probe_offset ~= 0;
            end

            valid_levels = curr_level_index(valid);
            n_valid = nnz(valid);

            if curr_feature == 1
                curr_feature_name = repmat("Contrast", n_valid, 1);
                curr_stimulus_value = p.contrast(valid_levels)';
            elseif curr_feature == 2
                curr_feature_name = repmat("Filter width", n_valid, 1);
                curr_stimulus_value = p.orientation_bp_filter_width(valid_levels)';
            else
                error('Unexpected feature index %d in block %d.', curr_feature, n_block);
            end

            feature = [feature; curr_feature_name]; %#ok<AGROW>
            stimulus_value = [stimulus_value; curr_stimulus_value]; %#ok<AGROW>
            correct = [correct; saved_correct(valid)]; %#ok<AGROW>
            run_number = [run_number; repmat(i_file, n_valid, 1)]; %#ok<AGROW>
            block_number = [block_number; repmat(n_block, n_valid, 1)]; %#ok<AGROW>
            level_index = [level_index; valid_levels]; %#ok<AGROW>
            abs_probe_offset = [abs_probe_offset; abs(probe_offset(valid))]; %#ok<AGROW>
        end
    end

    trials = table(feature, stimulus_value, level_index, abs_probe_offset, correct, ...
        run_number, block_number, 'VariableNames', ...
        {'Feature', 'StimulusValue', 'LevelIndex', 'AbsProbeOffset', ...
        'Correct', 'Run', 'Block'});

    if isempty(trials)
        warning('No valid calibration trials found for Subject %s.', subj_ID);
        continue;
    end

    %% Summarize raw performance at each physical stimulus value

    summary_table = groupsummary(trials, {'Feature', 'StimulusValue'}, ...
        {'mean', 'std'}, 'Correct');
    summary_table.SE = summary_table.std_Correct ./ sqrt(summary_table.GroupCount);
    summary_table.Properties.VariableNames{'mean_Correct'} = 'ProportionCorrect';

    fprintf('\nSubject %s calibration performance\n', subj_ID);
    fprintf('Files loaded: %d\n', numel(files));
    fprintf('Scored trials checked: %d\n', n_scored_trials);
    fprintf('Saved-versus-recomputed scoring mismatches: %d\n\n', ...
        n_scoring_mismatches);
    disp(summary_table(:, {'Feature', 'StimulusValue', 'GroupCount', ...
        'ProportionCorrect', 'SE'}));

    %% Summarize performance separately for each absolute probe offset

    offset_summary = groupsummary(trials, ...
        {'Feature', 'StimulusValue', 'AbsProbeOffset'}, ...
        {'mean', 'std'}, 'Correct');
    offset_summary.SE = offset_summary.std_Correct ./ ...
        sqrt(offset_summary.GroupCount);
    offset_summary.Properties.VariableNames{'mean_Correct'} = ...
        'ProportionCorrect';

    fprintf('\nPerformance separated by absolute probe offset\n\n');
    disp(offset_summary(:, {'Feature', 'StimulusValue', 'AbsProbeOffset', ...
        'GroupCount', 'ProportionCorrect', 'SE'}));

    %% Plot raw performance only (no Weibull fit)

    figure('Color', 'w', 'Position', [100 100 1000 420], ...
        'Name', ['S' subj_ID ' Raw Calibration Performance']);

    feature_names = ["Contrast", "Filter width"];
    for i_feature = 1:numel(feature_names)

        subplot(1, 2, i_feature);
        rows = summary_table.Feature == feature_names(i_feature);
        feature_summary = sortrows(summary_table(rows, :), 'StimulusValue');

        x = feature_summary.StimulusValue;
        y = feature_summary.ProportionCorrect;
        y_error = feature_summary.SE;

        if connect_points
            plot(x, y, '-', 'Color', [0.45 0.45 0.45], ...
                'LineWidth', 1, 'HandleVisibility', 'off');
            hold on;
        end
        errorbar(x, y, y_error, 'o', 'Color', [0 0 0], ...
            'MarkerFaceColor', [0 0 0], 'MarkerSize', 6, ...
            'LineWidth', 1, 'CapSize', 0);

        set(gca, 'XScale', 'log', 'XMinorTick', 'off', ...
            'TickDir', 'out', 'TickLength', [0.020 0.020], ...
            'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1);
        ylim([0 1]);
        ylabel('Proportion Correct', 'FontName', 'Helvetica', 'FontSize', 14);
        box off;
        axis square;

        if i_feature == 1
            xlabel('Contrast', 'FontName', 'Helvetica', 'FontSize', 14);
            title('Raw Performance vs. Contrast');
        else
            xlabel('Orientation Filter Width (degrees)', ...
                'FontName', 'Helvetica', 'FontSize', 14);
            title('Raw Performance vs. Filter Width');
            set(gca, 'XDir', 'reverse');
        end
    end

    sgtitle(['Subject ' subj_ID ': Calibration Data Check']);

    %% Plot performance by absolute probe offset

    figure('Color', 'w', 'Position', [100 100 1100 450], ...
        'Name', ['S' subj_ID ' Calibration Performance by Probe Offset']);

    offset_values = unique(trials.AbsProbeOffset);
    offset_colors = lines(numel(offset_values));

    for i_feature = 1:numel(feature_names)

        subplot(1, 2, i_feature);
        hold on;

        for i_offset = 1:numel(offset_values)
            rows = offset_summary.Feature == feature_names(i_feature) & ...
                offset_summary.AbsProbeOffset == offset_values(i_offset);
            curr_summary = sortrows(offset_summary(rows, :), 'StimulusValue');

            if isempty(curr_summary)
                continue;
            end

            errorbar(curr_summary.StimulusValue, ...
                curr_summary.ProportionCorrect, curr_summary.SE, '-o', ...
                'Color', offset_colors(i_offset, :), ...
                'MarkerFaceColor', offset_colors(i_offset, :), ...
                'MarkerSize', 5, 'LineWidth', 1, 'CapSize', 0, ...
                'DisplayName', [num2str(offset_values(i_offset)) char(176)]);
        end

        yline(0.5, ':', 'Chance', 'HandleVisibility', 'off');
        set(gca, 'XScale', 'log', 'XMinorTick', 'off', ...
            'TickDir', 'out', 'TickLength', [0.020 0.020], ...
            'FontName', 'Helvetica', 'FontSize', 13, 'LineWidth', 1);
        ylim([0 1]);
        ylabel('Proportion Correct', 'FontName', 'Helvetica', 'FontSize', 14);
        box off;
        axis square;

        if i_feature == 1
            xlabel('Contrast', 'FontName', 'Helvetica', 'FontSize', 14);
            title('Contrast Performance by Probe Offset');
        else
            xlabel('Orientation Filter Width (degrees)', ...
                'FontName', 'Helvetica', 'FontSize', 14);
            title('Filter-Width Performance by Probe Offset');
            set(gca, 'XDir', 'reverse');
        end

        legend('Location', 'bestoutside');
    end

    sgtitle(['Subject ' subj_ID ': Absolute Probe-Offset Check']);

end
