classdef SimulationDashboard < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        TabGroup             matlab.ui.container.TabGroup
        SubjectTab           matlab.ui.container.Tab
        GroupTab             matlab.ui.container.Tab

        % Data Source
        LoadDataButton       matlab.ui.control.Button
        DataSourceLabel      matlab.ui.control.Label

        % Subject Tab Components
        SubjectDropDown      matlab.ui.control.DropDown
        SubjectDropDownLabel matlab.ui.control.Label
        AxContrast           matlab.ui.control.UIAxes
        AxPrecision          matlab.ui.control.UIAxes

        % Group Tab Components
        GroupToggle          matlab.ui.control.Switch
        AxGrpContrast        matlab.ui.control.UIAxes
        AxGrpPrecision       matlab.ui.control.UIAxes
        StatsTextArea        matlab.ui.control.TextArea

        % Internal Data
        isSimulated          logical
        simData              struct
        realData             struct
        subjIDs              cell
    end

    methods (Access = private)

        % Button pushed function: LoadDataButton
        function loadData(app, ~)
            [file, path] = uigetfile('*.mat', 'Select Results Data File');
            if isequal(file, 0)
                return;
            end

            filepath = fullfile(path, file);
            data = load(filepath);

            % Auto-detect simulated data
            if isfield(data, 'hypothesis') || isfield(data, 'gt_amps_matrix') || isfield(data, 'recovered_all')
                app.isSimulated = true;
                app.simData = data;
                app.DataSourceLabel.Text = ['Data: Simulated (' data.hypothesis ')'];

                % Populate subject dropdown
                app.subjIDs = {data.recovered_all.subj_id};
                app.SubjectDropDown.Items = app.subjIDs;

                app.GroupToggle.Value = 'Simulated';
                app.GroupToggle.Enable = 'off'; % Only sim data available
            elseif isfield(data, 'bae2024_results')
                app.isSimulated = true;
                app.simData = data;
                app.DataSourceLabel.Text = 'Data: Simulated (Bae2024 categories)';
                app.subjIDs = {data.bae2024_results.subject.subj_id};
                if isempty(app.subjIDs)
                    app.SubjectDropDown.Items = {'None'};
                else
                    app.SubjectDropDown.Items = app.subjIDs;
                    app.SubjectDropDown.Value = app.subjIDs{1};
                end
                app.GroupToggle.Value = 'Simulated';
                app.GroupToggle.Enable = 'off';
            else
                app.isSimulated = false;
                app.realData = data;
                app.DataSourceLabel.Text = 'Data: Real Experiment';
                % E.g., load real subject IDs here if available
                app.GroupToggle.Value = 'Real';
            end

            % Update plots
            updateSubjectPlots(app);
            updateGroupPlots(app);
        end

        % Value changed function: SubjectDropDown
        function subjectChanged(app, ~)
            updateSubjectPlots(app);
        end

        function updateSubjectPlots(app)
            if isempty(app.subjIDs)
                return;
            end

            subj_idx = find(strcmp(app.subjIDs, app.SubjectDropDown.Value));
            if isempty(subj_idx), return; end

            if app.isSimulated
                ps = plotSettings();
                if isfield(app.simData, 'bae2024_results')
                    cla(app.AxContrast);
                    cla(app.AxPrecision);
                    if isempty(app.subjIDs) || strcmp(app.SubjectDropDown.Value, 'None')
                        return;
                    end
                    subj_idx = find(strcmp(app.subjIDs, app.SubjectDropDown.Value), 1);
                    if isempty(subj_idx), return; end
                    subj = app.simData.bae2024_results.subject(subj_idx);
                    amp_vals = nan(1,4);
                    labels = {'Cardinal','45-like','Congruent','Incongruent'};
                    for k = 1:4
                        if ~isempty(subj.category_fit{k}) && isfield(subj.category_fit{k}, 'params')
                            amp_vals(k) = subj.category_fit{k}.params(1);
                        end
                    end
                    bar(app.AxContrast, 1:4, amp_vals, 'FaceColor', ps.colors.blue); hold(app.AxContrast, 'on');
                    app.AxContrast.Title.String = 'Bae2024 Subject Amplitudes';
                    app.AxContrast.XTick = 1:4;
                    app.AxContrast.XTickLabel = labels;
                    app.AxContrast.YLabel.String = 'DoG Amplitude (deg)';

                    side_vals = [subj.intercept_cw, subj.intercept_ccw];
                    bar(app.AxPrecision, 1:2, side_vals, 'FaceColor', ps.colors.green); hold(app.AxPrecision, 'on');
                    app.AxPrecision.Title.String = 'Bae2024 Subject Intercepts';
                    app.AxPrecision.XTick = [1 2];
                    app.AxPrecision.XTickLabel = {'CW side','CCW side'};
                    app.AxPrecision.YLabel.String = 'DoG Intercept (deg)';
                    return;
                end
                % Plot Recovered vs GT
                % Contrast
                cla(app.AxContrast);
                rec_c = app.simData.recovered_all(subj_idx).rec_amp(1,:);
                gt_c = app.simData.gt_amps_matrix(subj_idx, 1, :);
                gt_c = squeeze(gt_c)';

                plot(app.AxContrast, 1:3, rec_c, 'o-', 'Color', ps.colors.blue, 'MarkerFaceColor', ps.colors.blue, 'MarkerSize', 6, 'LineWidth', ps.line_width); hold(app.AxContrast, 'on');
                plot(app.AxContrast, 1:3, gt_c, '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
                app.AxContrast.Title.String = 'Contrast Modulation (DoG Amplitude)';
                app.AxContrast.Title.FontSize = ps.axes_label_font_size;
                app.AxContrast.XLabel.String = 'Feature Level (1=Hard, 3=Easy)';
                app.AxContrast.XLabel.FontSize = ps.axes_label_font_size;
                app.AxContrast.YLabel.String = 'Amplitude (deg)';
                app.AxContrast.YLabel.FontSize = ps.axes_label_font_size;
                app.AxContrast.XTick = [1 2 3];
                legend(app.AxContrast, {'Recovered', 'Ground Truth'}, 'Location', 'best', 'FontSize', ps.axes_tick_font_size);

                % Precision
                cla(app.AxPrecision);
                rec_p = app.simData.recovered_all(subj_idx).rec_amp(2,:);
                gt_p = app.simData.gt_amps_matrix(subj_idx, 2, :);
                gt_p = squeeze(gt_p)';

                plot(app.AxPrecision, 1:3, rec_p, 'o-', 'Color', ps.colors.green, 'MarkerFaceColor', ps.colors.green, 'MarkerSize', 6, 'LineWidth', ps.line_width); hold(app.AxPrecision, 'on');
                plot(app.AxPrecision, 1:3, gt_p, '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
                app.AxPrecision.Title.String = 'Precision Modulation (DoG Amplitude)';
                app.AxPrecision.Title.FontSize = ps.axes_label_font_size;
                app.AxPrecision.XLabel.String = 'Feature Level (1=Hard, 3=Easy)';
                app.AxPrecision.XLabel.FontSize = ps.axes_label_font_size;
                app.AxPrecision.XTick = [1 2 3];
                legend(app.AxPrecision, {'Recovered', 'Ground Truth'}, 'Location', 'best', 'FontSize', ps.axes_tick_font_size);
            else
                % Real data plotting logic here
            end
        end

        function updateGroupPlots(app)
            if app.isSimulated
                ps = plotSettings();
                if isfield(app.simData, 'bae2024_results')
                    cla(app.AxGrpContrast);
                    cla(app.AxGrpPrecision);
                    grp = app.simData.bae2024_results.group;
                    summ = app.simData.bae2024_results.summary;
                    labels = {'Cardinal','45-like','Congruent','Incongruent'};
                    mean_amp = summ.mean_amp;
                    sem_amp = summ.sem_amp;
                    errorbar(app.AxGrpContrast, 1:4, mean_amp, sem_amp, 'o-', 'Color', ps.colors.blue, ...
                        'MarkerFaceColor', ps.colors.blue, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width);
                    app.AxGrpContrast.Title.String = 'Bae2024 Group Amplitudes';
                    app.AxGrpContrast.XTick = 1:4;
                    app.AxGrpContrast.XTickLabel = labels;
                    app.AxGrpContrast.YLabel.String = 'DoG Amplitude (deg)';

                    mean_side = mean(grp.intercept_side, 1, 'omitnan');
                    sem_side = std(grp.intercept_side, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(grp.intercept_side), 1));
                    errorbar(app.AxGrpPrecision, 1:2, mean_side, sem_side, 'o-', 'Color', ps.colors.green, ...
                        'MarkerFaceColor', ps.colors.green, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width);
                    app.AxGrpPrecision.Title.String = 'Bae2024 Intercept Shift';
                    app.AxGrpPrecision.XTick = [1 2];
                    app.AxGrpPrecision.XTickLabel = {'CW side','CCW side'};
                    app.AxGrpPrecision.YLabel.String = 'DoG Intercept (deg)';

                    s = summ.stats;
                    stats_str = sprintf(['Congruent A>0: p=%.4f\n' ...
                        'Incongruent A>0: p=%.4f\n' ...
                        'Congruent-Incongruent: p=%.4f\n' ...
                        'CW-CCW intercept: p=%.4f'], ...
                        s.amp_congruent_gt0.p, s.amp_incongruent_gt0.p, ...
                        s.amp_congruent_vs_incongruent.p, s.intercept_cw_vs_ccw.p);
                    app.StatsTextArea.Value = stats_str;
                    return;
                end
                cla(app.AxGrpContrast);
                errorbar(app.AxGrpContrast, 1:3, app.simData.mean_rec(1,:), app.simData.sem_rec(1,:), 'o-', 'Color', ps.colors.blue, 'MarkerFaceColor', ps.colors.blue, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width);
                hold(app.AxGrpContrast, 'on');
                plot(app.AxGrpContrast, 1:3, app.simData.mean_gt(1,:), '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
                app.AxGrpContrast.Title.String = 'Group Average: Contrast';
                app.AxGrpContrast.Title.FontSize = ps.axes_label_font_size;
                app.AxGrpContrast.XTick = [1 2 3];

                cla(app.AxGrpPrecision);
                errorbar(app.AxGrpPrecision, 1:3, app.simData.mean_rec(2,:), app.simData.sem_rec(2,:), 'o-', 'Color', ps.colors.green, 'MarkerFaceColor', ps.colors.green, 'MarkerSize', 6, 'CapSize', 0, 'LineWidth', ps.line_width);
                hold(app.AxGrpPrecision, 'on');
                plot(app.AxGrpPrecision, 1:3, app.simData.mean_gt(2,:), '--', 'Color', ps.colors.red, 'LineWidth', ps.line_width);
                app.AxGrpPrecision.Title.String = 'Group Average: Precision';
                app.AxGrpPrecision.Title.FontSize = ps.axes_label_font_size;
                app.AxGrpPrecision.XTick = [1 2 3];

                % Stats text
                stats_str = sprintf('Amplitude R^2: %.3f | RMSE: %.3f deg\nWidth R^2: %.3f | RMSE: %.4f', ...
                    app.simData.R2, app.simData.RMSE, app.simData.R2_w, app.simData.RMSE_w);
                app.StatsTextArea.Value = stats_str;
            end
        end

        % Create UI Figure and components
        function createComponents(app)
            ps = plotSettings();
            app.UIFigure = uifigure('Position', [100 100 800 600], 'Name', 'Experiment 2 Dashboard', 'Color', 'w');

            % Load Data
            app.LoadDataButton = uibutton(app.UIFigure, 'push', 'Position', [20, 550, 100, 30], ...
                'Text', 'Load Data', 'ButtonPushedFcn', @(~,~) loadData(app, []));
            app.DataSourceLabel = uilabel(app.UIFigure, 'Position', [140, 550, 300, 30], 'Text', 'Data: None');

            % Tabs
            app.TabGroup = uitabgroup(app.UIFigure, 'Position', [20, 20, 760, 500]);
            app.SubjectTab = uitab(app.TabGroup, 'Title', 'Subject View');
            app.GroupTab = uitab(app.TabGroup, 'Title', 'Group View');

            % --- Subject Tab ---
            app.SubjectDropDownLabel = uilabel(app.SubjectTab, 'Position', [20, 430, 60, 22], 'Text', 'Subject:');
            app.SubjectDropDown = uidropdown(app.SubjectTab, 'Position', [80, 430, 100, 22], ...
                'Items', {'None'}, 'ValueChangedFcn', @(~,~) subjectChanged(app, []));

            app.AxContrast = uiaxes(app.SubjectTab, 'Position', [20, 50, 350, 350]);
            app.AxPrecision = uiaxes(app.SubjectTab, 'Position', [390, 50, 350, 350]);

            % --- Group Tab ---
            app.GroupToggle = uiswitch(app.GroupTab, 'Position', [50, 430, 50, 22], 'Items', {'Real', 'Simulated'});

            app.AxGrpContrast = uiaxes(app.GroupTab, 'Position', [20, 50, 350, 350]);
            app.AxGrpPrecision = uiaxes(app.GroupTab, 'Position', [390, 50, 350, 350]);

            app.StatsTextArea = uitextarea(app.GroupTab, 'Position', [200, 410, 400, 60], 'Editable', 'off');

            % Apply styles to all axes
            axesList = [app.AxContrast, app.AxPrecision, app.AxGrpContrast, app.AxGrpPrecision];
            for i = 1:length(axesList)
                ax = axesList(i);
                ax.FontName = ps.font_type;
                ax.FontSize = ps.axes_tick_font_size;
                ax.TickDir = 'out';
                ax.TickLength = [ps.tick_length, ps.tick_length];
                ax.LineWidth = ps.line_width;
                ax.Box = 'off';
                ax.XGrid = 'off'; ax.YGrid = 'off';
            end
        end
    end

    methods (Access = public)
        % Construct app
        function app = SimulationDashboard
            % Ensure plotting scripts are on path
            addpath(fullfile(fileparts(mfilename('fullpath')), '../analyses/plotting'));
            createComponents(app);
            app.isSimulated = false;
        end

        % Delete app
        function delete(app)
            delete(app.UIFigure);
        end
    end
end