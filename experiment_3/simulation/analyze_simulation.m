% analyze_simulation.m
% Analyzes the simulated data, fits DoG models to recover parameters,
% compares against ground-truth, and outputs figures and reports.

clear all; close all; clc;

% Directories
sim_dir = fileparts(mfilename('fullpath'));
exp_dir = fullfile(sim_dir, 'experiment');
figs_grp_dir = fullfile(exp_dir, 'figures', 'group');
data_dir = fullfile(sim_dir, 'data');

if ~exist(figs_grp_dir, 'dir'), mkdir(figs_grp_dir); end

addpath(fullfile(sim_dir, '../analyses/plotting'));
ps = plotSettings();

% Load ground truth
load('simulated_ground_truth.mat', 'gt_all', 'hypothesis');
num_subjs = length(gt_all);

%% Data Extraction & Subject-Level Fitting
% For each subject, calculate bias per delta_theta window for each feature/level combo
delta_theta_edges = -90:20:90;
delta_theta_centers = delta_theta_edges(1:end-1) + diff(delta_theta_edges)/2;

recovered_all = struct();

% DoG function: y = A * c * w * x * exp(-(w*x)^2)
% Normalizing constant c = exp(0.5)
dog_func = @(params, x) params(1) * exp(0.5) * params(2) * x .* exp(-(params(2) * x).^2);

for i = 1:num_subjs
    subj_id = gt_all(i).subj_id;
    disp(['Analyzing S' subj_id '...']);

    figs_subj_dir = fullfile(exp_dir, 'figures', 'subject', 'sd');
    if ~exist(figs_subj_dir, 'dir'), mkdir(figs_subj_dir); end

    % Aggregate data for the subject (Sessions 2-5 are saved as SD_Noise_Exp2_S9XX_RunX.mat)
    % Specifically avoiding 'Calibration' and 'Training'
    files = dir(fullfile(data_dir, subj_id, 'SD_Noise_Exp2_S*_Run*.mat'));
    files = files(cellfun(@isempty, regexp({files.name}, 'Calibration|Training')));

    if isempty(files)
        disp(['No mock data found for S' subj_id]);
        continue;
    end

    % Data structure for fitting: cell array for [Feature, Level]
    subj_data = cell(2, 3);

    for f = 1:length(files)
        load(fullfile(files(f).folder, files(f).name), 'run_info');

        % If run_info somehow doesn't exist or is corrupted, skip it
        if ~exist('run_info', 'var') || isempty(run_info)
            continue;
        end
        p = run_info.p;
        responses = run_info.behav_data.responses;

        for n_block = 1:p.num_blocks
            feature = p.feature_order(n_block);
            test_ori = p.trial_events(:, 1, n_block);
            probe_ori = p.trial_events(:, 2, n_block);
            levels = p.trial_events(:, 3, n_block);

            probe_offset = probe_ori - test_ori;
            probe_offset(probe_offset > 90) = probe_offset(probe_offset > 90) - 180;
            probe_offset(probe_offset < -90) = probe_offset(probe_offset < -90) + 180;

            % delta theta
            delta_theta = [0; test_ori(1:end-1) - test_ori(2:end)];
            delta_theta(delta_theta > 90) = delta_theta(delta_theta > 90) - 180;
            delta_theta(delta_theta < -90) = delta_theta(delta_theta < -90) + 180;

            prev_levels = [1; levels(1:end-1)];

            for t_idx = 2:length(levels)
                p_lvl = prev_levels(t_idx);
                dt = delta_theta(t_idx);
                resp = responses(t_idx, n_block); % 1 = CW, 0 = CCW
                po = probe_offset(t_idx);

                % Store trial info
                if isempty(subj_data{feature, p_lvl})
                    subj_data{feature, p_lvl} = [dt, po, resp];
                else
                    subj_data{feature, p_lvl}(end+1, :) = [dt, po, resp];
                end
            end
        end
    end

    % Fit response bias per window, then DoG
    rec_amp = zeros(2, 3);
    rec_w = zeros(2, 3);

    fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 600], 'Name', ['S' subj_id ' Recovered SD'], 'Color', 'w');
    tiledlayout(fig, 2, 3, 'Padding', 'compact');

    for feat = 1:2
        for lvl = 1:3
            data = subj_data{feat, lvl};
            if isempty(data), continue; end

            dt = data(:,1);
            po = data(:,2);
            resp = data(:,3);

            % Simple windowed bias estimation (average offset where P(CW) = 0.5)
            % Or just mean error. Since it's binary, fit logistic per window.
            bias_windowed = zeros(size(delta_theta_centers));
            for w_idx = 1:length(delta_theta_centers)
                idx = dt >= delta_theta_edges(w_idx) & dt < delta_theta_edges(w_idx+1);
                if sum(idx) > 15
                    try
                        mdl = fitglm(po(idx), resp(idx), 'Distribution', 'binomial', ...
                            'Link', 'probit', 'LikelihoodPenalty', 'jeffreys-prior');
                        b0 = mdl.Coefficients.Estimate(1);
                        b1 = mdl.Coefficients.Estimate(2);
                        if abs(b1) < 1e-4
                            bias_windowed(w_idx) = NaN;
                        else
                            bias_windowed(w_idx) = -b0 / b1;
                        end
                    catch
                        bias_windowed(w_idx) = NaN;
                    end
                else
                    bias_windowed(w_idx) = NaN;
                end
            end

            % Fit DoG
            valid = ~isnan(bias_windowed) & ~isinf(bias_windowed) & abs(bias_windowed) < 15;
            if sum(valid) > 3
                x_val = delta_theta_centers(valid);
                y_val = bias_windowed(valid);

                options = optimoptions('lsqcurvefit', 'Display', 'off');
                lb = [0, 0.01]; ub = [15, 0.2];
                init_guess = [5, 0.05];
                try
                    best_params = lsqcurvefit(dog_func, init_guess, x_val, y_val, lb, ub, options);
                    rec_amp(feat, lvl) = best_params(1);
                    rec_w(feat, lvl) = best_params(2);
                catch
                    rec_amp(feat, lvl) = NaN;
                    rec_w(feat, lvl) = NaN;
                end
            else
                rec_amp(feat, lvl) = NaN;
                rec_w(feat, lvl) = NaN;
            end

            % Plot subject figure
            nexttile;
            plot(delta_theta_centers, bias_windowed, 'o', 'MarkerFaceColor', ps.colors.black, 'Color', ps.colors.black, 'MarkerSize', 6); hold on;
            x_smooth = linspace(-90, 90, 100);
            if ~isnan(rec_amp(feat, lvl))
                plot(x_smooth, dog_func([rec_amp(feat, lvl), rec_w(feat, lvl)], x_smooth), '-', 'Color', ps.colors.blue, 'LineWidth', ps.line_width);
            end
            % Ground truth
            gt_a = gt_all(i).gt.dog_amp(feat, lvl);
            gt_w = gt_all(i).gt.dog_w(feat, lvl);
            plot(x_smooth, dog_func([gt_a, gt_w], x_smooth), '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
            title(sprintf('F%d L%d', feat, lvl), 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
            xlim([-90 90]); ylim([-10 10]);
            set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
            box off; axis square;
            if feat == 1 && lvl == 1, legend('Data', 'Fit', 'GT', 'Location', 'best', 'FontSize', ps.axes_tick_font_size); end
        end
    end

    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperOrientation', 'landscape');
    print(fig, fullfile(figs_subj_dir, sprintf('S%s_SD_recovery.pdf', subj_id)), '-dpdf', '-bestfit');
    close(fig);

    recovered_all(i).subj_id = subj_id;
    recovered_all(i).rec_amp = rec_amp;
    recovered_all(i).rec_w = rec_w;
end

%% Group Level Analysis & RMSE/R2 Calculation
% Aggregate params
valid_subjs = find(~cellfun(@isempty, {recovered_all.subj_id}));
n_valid = length(valid_subjs);

gt_amps = reshape([gt_all.gt], num_subjs, 1);
gt_amps_matrix = zeros(n_valid, 2, 3);
rec_amps_matrix = zeros(n_valid, 2, 3);

gt_w_matrix = zeros(n_valid, 2, 3);
rec_w_matrix = zeros(n_valid, 2, 3);

for idx = 1:n_valid
    i = valid_subjs(idx);
    gt_amps_matrix(idx, :, :) = gt_all(i).gt.dog_amp;
    rec_amps_matrix(idx, :, :) = recovered_all(i).rec_amp;

    gt_w_matrix(idx, :, :) = gt_all(i).gt.dog_w;
    rec_w_matrix(idx, :, :) = recovered_all(i).rec_w;
end

% Flatten to calculate global R^2 and RMSE for Amplitude
gt_flat = gt_amps_matrix(:);
rec_flat = rec_amps_matrix(:);
valid = ~isnan(rec_flat);

RMSE = sqrt(mean((gt_flat(valid) - rec_flat(valid)).^2));
corr_matrix = corrcoef(gt_flat(valid), rec_flat(valid));
if numel(corr_matrix) > 1
    R2 = corr_matrix(1,2)^2;
else
    R2 = NaN;
end

mean_gt = squeeze(mean(gt_amps_matrix, 1));
mean_rec = squeeze(mean(rec_amps_matrix, 1, 'omitnan'));
sem_rec = squeeze(std(rec_amps_matrix, 0, 1, 'omitnan')) / sqrt(n_valid);

% Flatten to calculate global R^2 and RMSE for Width
gt_w_flat = gt_w_matrix(:);
rec_w_flat = rec_w_matrix(:);
valid_w = ~isnan(rec_w_flat);

RMSE_w = sqrt(mean((gt_w_flat(valid_w) - rec_w_flat(valid_w)).^2));
corr_w_matrix = corrcoef(gt_w_flat(valid_w), rec_w_flat(valid_w));
if numel(corr_w_matrix) > 1
    R2_w = corr_w_matrix(1,2)^2;
else
    R2_w = NaN;
end

mean_gt_w = squeeze(mean(gt_w_matrix, 1));
mean_rec_w = squeeze(mean(rec_w_matrix, 1, 'omitnan'));
sem_rec_w = squeeze(std(rec_w_matrix, 0, 1, 'omitnan')) / sqrt(n_valid);

%% Plot Group Recovery
y_upper = ceil(max([mean_gt(:); mean_rec(:)] + max(sem_rec(:))) + 2);
y_upper = max(y_upper, 10);

fig_grp = figure('Visible', 'off', 'Position', [200, 200, 800, 400], 'Name', 'Group Average Amplitudes', 'Color', 'w');
tiledlayout(fig_grp, 1, 2, 'Padding', 'compact');

nexttile;
errorbar(1:3, mean_rec(1,:), sem_rec(1,:), 'o-', 'Color', ps.colors.blue, 'MarkerFaceColor', ps.colors.blue, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
plot(1:3, mean_gt(1,:), '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
title('Contrast (F1)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
xlabel('Level', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('DoG Amplitude', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
legend('Recovered', 'Ground Truth', 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
ylim([0 y_upper]);
xticks([1 2 3]);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

nexttile;
errorbar(1:3, mean_rec(2,:), sem_rec(2,:), 'o-', 'Color', ps.colors.green, 'MarkerFaceColor', ps.colors.green, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width); hold on;
plot(1:3, mean_gt(2,:), '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
title('Precision (F2)', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
xlabel('Level', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylabel('DoG Amplitude', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
ylim([0 y_upper]);
xticks([1 2 3]);
set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
box off; axis square;

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print(fig_grp, fullfile(figs_grp_dir, 'Group_Amplitude_Recovery.pdf'), '-dpdf', '-bestfit');
close(fig_grp);

%% Calibration vs. SD Correlation Analysis
calib_params = nan(n_valid, 4); % [contrast_alpha, contrast_beta, filter_alpha, filter_beta]
mean_rec_amp_subj = nan(n_valid, 1);

for idx = 1:n_valid
    i = valid_subjs(idx);
    sid = recovered_all(i).subj_id;
    calib_file = fullfile(data_dir, sid, ['S' sid '_calibrated_levels.mat']);
    if exist(calib_file, 'file')
        tmp = load(calib_file, 'calib');
        calib_params(idx, 1) = tmp.calib.fit_params.contrast_alpha;
        calib_params(idx, 2) = tmp.calib.fit_params.contrast_beta;
        calib_params(idx, 3) = tmp.calib.fit_params.filter_alpha;
        calib_params(idx, 4) = tmp.calib.fit_params.filter_beta;
    end
    mean_rec_amp_subj(idx) = mean(rec_amps_matrix(idx, :), 'all', 'omitnan');
end

% Plot correlation figure (4 panels: 2 contrast + 2 filter Weibull params)
fig_corr = figure('Visible', 'off', 'Position', [300, 300, 1200, 500], 'Name', 'Calibration vs SD', 'Color', 'w');
tiledlayout(fig_corr, 1, 4, 'Padding', 'compact');
labels = {'Contrast \alpha', 'Contrast \beta', 'Filter \alpha_{fw}', 'Filter \beta_{fw}'};
corr_r = zeros(1,4);
corr_p = zeros(1,4);

for p_idx = 1:4
    nexttile;
    x_data = calib_params(:, p_idx);
    y_data = mean_rec_amp_subj;
    valid_idx = ~isnan(x_data) & ~isnan(y_data);
    x_data = x_data(valid_idx);
    y_data = y_data(valid_idx);

    plot(x_data, y_data, 'o', 'Color', ps.colors.black, 'MarkerFaceColor', ps.colors.gray, 'MarkerSize', 6); hold on;

    if length(x_data) > 2
        [R, P] = corrcoef(x_data, y_data);
        corr_r(p_idx) = R(1,2);
        corr_p(p_idx) = P(1,2);

        mdl = polyfit(x_data, y_data, 1);
        x_fit = [min(x_data) max(x_data)];
        y_fit = polyval(mdl, x_fit);
        plot(x_fit, y_fit, '-', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
        title(sprintf('%s\nr=%.2f, p=%.3f', labels{p_idx}, corr_r(p_idx), corr_p(p_idx)), 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
    else
        title(labels{p_idx}, 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
    end
    xlabel(labels{p_idx}, 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type);
    if p_idx == 1
        ylabel('Mean SD Amplitude', 'FontSize', ps.axes_label_font_size, 'FontName', ps.font_type);
    end
    set(gca, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', ps.axes_tick_font_size, 'FontName', ps.font_type, 'LineWidth', ps.line_width);
    box off; axis square;
end
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperOrientation', 'landscape');
print(fig_corr, fullfile(figs_grp_dir, 'Group_Calibration_SD_Correlation.pdf'), '-dpdf', '-bestfit');
close(fig_corr);

%% Write Outputs

% 1. TXT stats
fid = fopen(fullfile(exp_dir, 'simulation_summary_stats.txt'), 'w');
fprintf(fid, 'Simulation Summary Statistics\n');
fprintf(fid, '=============================\n');
fprintf(fid, 'Hypothesis Injected: %s\n', hypothesis);
fprintf(fid, 'Number of Subjects: %d\n', num_subjs);
fprintf(fid, '\n--- Amplitude Recovery Metrics ---\n');
fprintf(fid, 'Overall Amplitude R^2: %.3f\n', R2);
fprintf(fid, 'Overall Amplitude RMSE: %.3f degrees\n', RMSE);
fprintf(fid, '\nMean Ground Truth Amplitudes (Contrast): %.2f, %.2f, %.2f\n', mean_gt(1,1), mean_gt(1,2), mean_gt(1,3));
fprintf(fid, 'Mean Recovered Amplitudes (Contrast): %.2f, %.2f, %.2f\n', mean_rec(1,1), mean_rec(1,2), mean_rec(1,3));
fprintf(fid, '\nMean Ground Truth Amplitudes (Precision): %.2f, %.2f, %.2f\n', mean_gt(2,1), mean_gt(2,2), mean_gt(2,3));
fprintf(fid, 'Mean Recovered Amplitudes (Precision): %.2f, %.2f, %.2f\n', mean_rec(2,1), mean_rec(2,2), mean_rec(2,3));

fprintf(fid, '\n--- Width Recovery Metrics ---\n');
fprintf(fid, 'Overall Width R^2: %.3f\n', R2_w);
fprintf(fid, 'Overall Width RMSE: %.3f\n', RMSE_w);
fprintf(fid, '\nMean Ground Truth Widths (Contrast): %.3f, %.3f, %.3f\n', mean_gt_w(1,1), mean_gt_w(1,2), mean_gt_w(1,3));
fprintf(fid, 'Mean Recovered Widths (Contrast): %.3f, %.3f, %.3f\n', mean_rec_w(1,1), mean_rec_w(1,2), mean_rec_w(1,3));
fprintf(fid, '\nMean Ground Truth Widths (Precision): %.3f, %.3f, %.3f\n', mean_gt_w(2,1), mean_gt_w(2,2), mean_gt_w(2,3));
fprintf(fid, 'Mean Recovered Widths (Precision): %.3f, %.3f, %.3f\n', mean_rec_w(2,1), mean_rec_w(2,2), mean_rec_w(2,3));

fprintf(fid, '\n--- Calibration vs SD Correlations ---\n');
fprintf(fid, 'Contrast Alpha vs Mean SD:    r = %.3f, p = %.3f\n', corr_r(1), corr_p(1));
fprintf(fid, 'Contrast Beta vs Mean SD:     r = %.3f, p = %.3f\n', corr_r(2), corr_p(2));
fprintf(fid, 'Filter Alpha_fw vs Mean SD:   r = %.3f, p = %.3f\n', corr_r(3), corr_p(3));
fprintf(fid, 'Filter Beta_fw vs Mean SD:    r = %.3f, p = %.3f\n', corr_r(4), corr_p(4));
fclose(fid);

% 2. Markdown report
fid = fopen(fullfile(exp_dir, 'simulation_report.md'), 'w');
fprintf(fid, '# Simulation Report: Experiment 2\n\n');
fprintf(fid, '## Configuration\n');
fprintf(fid, '- **Hypothesis Tested**: `%s`\n', hypothesis);
fprintf(fid, '- **Number of Subjects**: %d\n', num_subjs);
fprintf(fid, '- **Generated Output**: Stored in `data/9XX` (matching real experiment format)\n\n');

fprintf(fid, '## Parameter Recovery Accuracy\n');
fprintf(fid, 'We fit a Derivative of Gaussian (DoG) model to the simulated subject responses (binned into 15° \\Delta\\theta windows) using a probit link function.\n\n');
fprintf(fid, '- **Overall Amplitude $R^2$**: %.3f\n', R2);
fprintf(fid, '- **Overall Amplitude RMSE**: %.3f degrees\n', RMSE);
fprintf(fid, '- **Overall Width $R^2$**: %.3f\n', R2_w);
fprintf(fid, '- **Overall Width RMSE**: %.3f\n\n', RMSE_w);

fprintf(fid, '### Contrast Modulation (Level 1 $\\rightarrow$ Level 3)\n');
fprintf(fid, '- Amplitude Ground Truth: `[%.2f, %.2f, %.2f]`\n', mean_gt(1,1), mean_gt(1,2), mean_gt(1,3));
fprintf(fid, '- Amplitude Recovered: `[%.2f, %.2f, %.2f]`\n', mean_rec(1,1), mean_rec(1,2), mean_rec(1,3));
fprintf(fid, '- Width Ground Truth: `[%.3f, %.3f, %.3f]`\n', mean_gt_w(1,1), mean_gt_w(1,2), mean_gt_w(1,3));
fprintf(fid, '- Width Recovered: `[%.3f, %.3f, %.3f]`\n\n', mean_rec_w(1,1), mean_rec_w(1,2), mean_rec_w(1,3));

fprintf(fid, '### Precision Modulation (Level 1 $\\rightarrow$ Level 3)\n');
fprintf(fid, '- Amplitude Ground Truth: `[%.2f, %.2f, %.2f]`\n', mean_gt(2,1), mean_gt(2,2), mean_gt(2,3));
fprintf(fid, '- Amplitude Recovered: `[%.2f, %.2f, %.2f]`\n', mean_rec(2,1), mean_rec(2,2), mean_rec(2,3));
fprintf(fid, '- Width Ground Truth: `[%.3f, %.3f, %.3f]`\n', mean_gt_w(2,1), mean_gt_w(2,2), mean_gt_w(2,3));
fprintf(fid, '- Width Recovered: `[%.3f, %.3f, %.3f]`\n\n', mean_rec_w(2,1), mean_rec_w(2,2), mean_rec_w(2,3));

fprintf(fid, '## Calibration vs Serial Dependence\n');
fprintf(fid, 'Correlations between individual differences in psychometric parameters and baseline Serial Dependence (Mean Amplitude):\n\n');
fprintf(fid, '- **Contrast Threshold ($\\alpha$)**: $r = %.3f$, $p = %.3f$\n', corr_r(1), corr_p(1));
fprintf(fid, '- **Contrast Slope ($\\beta$)**: $r = %.3f$, $p = %.3f$\n', corr_r(2), corr_p(2));
fprintf(fid, '- **Filter Threshold ($\\alpha_{fw}$)**: $r = %.3f$, $p = %.3f$\n', corr_r(3), corr_p(3));
fprintf(fid, '- **Filter Slope ($\\beta_{fw}$)**: $r = %.3f$, $p = %.3f$\n\n', corr_r(4), corr_p(4));

fprintf(fid, '## Verification of Statistical Tests\n');
fprintf(fid, 'If the `%s` hypothesis was injected, you should see:\n', hypothesis);
fprintf(fid, '1. A robust RM-ANOVA interaction if `both_different` or either `_only` was chosen.\n');
fprintf(fid, '2. Strong fixed-effects parameters in the LME matching the generated slope.\n');
fclose(fid);

save('simulated_results.mat', 'recovered_all', 'gt_amps_matrix', 'rec_amps_matrix', 'gt_w_matrix', 'rec_w_matrix', 'mean_gt', 'mean_rec', 'sem_rec', 'mean_gt_w', 'mean_rec_w', 'sem_rec_w', 'R2', 'RMSE', 'R2_w', 'RMSE_w', 'hypothesis');
disp('Analysis complete! Outputs generated in experiment_2/simulation/experiment/');
