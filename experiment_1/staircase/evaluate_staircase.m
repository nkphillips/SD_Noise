% evaluate_staircase

clear all; close all; clc;

%% Prepare workspace

subj_IDs = {'001'};

script_dir = pwd;
staircase_data_dir = '../data'; addpath(staircase_data_dir);

%% Load staircase data

% Prepare a structure to store all the staircase data from each subject
all_staircases = cell(length(subj_IDs), 1);

for subj = 1:length(subj_IDs)

    save_filename = [staircase_data_dir '/' subj_IDs{subj} '/staircase_data_S' subj_IDs{subj} '_*.mat'];
    files_found = dir([staircase_data_dir '/' subj_IDs{subj}]);

    % Keep only the file names
    files_found = {files_found.name};

    % Keep only the file names that match staircase filename
    files_found = files_found(contains(files_found, 'staircase'));

    % Keep only the most recent file
    files_found = files_found{end};

    % Check if the file exists
    if ~isempty(files_found)
        
        % Store loaded staircase data struct into all_staircases
        loaded_staircase = load([staircase_data_dir '/' subj_IDs{subj} '/' files_found]);
        all_staircases{subj} = loaded_staircase.staircases;

        disp(['Loaded staircase data for subject ' subj_IDs{subj}]);

    else
        disp(['No staircase data found for subject ' subj_IDs{subj}]);
    end

end

% Convert all_staircases to a single structure
all_staircases = [all_staircases{:}];

num_conds = size(all_staircases(1).probe_offsets, 4);
num_levels = size(all_staircases(1).probe_offsets, 3);
cond_names = {'contrast', 'filter width'};
colors = {'b', 'r'}; % different colors for each staircase

%% Plot final probe offsets

figure_name = ['Final probe offsets S' subj_IDs{subj}];
figure('Name', figure_name, 'Color', 'w');

for subj = 1:length(subj_IDs)

    x = 1:num_conds;
    y = all_staircases(subj).final_probe_offsets;

    bar(x, y, 'k');

    % Format
    xticks(1:num_conds);
    xticklabels(cond_names);
    xlabel('Condition');
    ylabel('Final probe offset (°)');
    box off;
    set(gca, 'TickDir', 'out');

end

%% Plot staircases

for subj = 1:length(subj_IDs)

    % Figures for each condition and subplots for each level
    for cond = 1:num_conds

        figure_name = ['Staircase Trials S' subj_IDs{subj} ' ' cond_names{cond} ];
        figure('Name', figure_name, 'Color', 'w');

        for lvl = 1:num_levels

            subplot(1, num_levels, lvl);
            hold on;

            % Plot probe offsets for each staircase
            for n_sc = 1:size(all_staircases(subj).probe_offsets, 1)
            
                % Plot probe offset trajectory
                plot(all_staircases(subj).probe_offsets(n_sc,:,lvl,cond), colors{n_sc}, 'LineWidth', 1);
                
                % Mark reversals
                reversal_idx = find(all_staircases(subj).reversal_indices(n_sc,:,lvl,cond) == 1);
                reversal_values = all_staircases(subj).probe_offsets(n_sc,reversal_idx,lvl,cond);
                plot(reversal_idx, reversal_values, 'k.', 'MarkerSize', 15);
            end
            hold off;

            % Add convergence metrics
            for n_sc = 1:size(all_staircases(subj).probe_offsets, 1)
                reversal_idx = find(all_staircases(subj).reversal_indices(n_sc,:,lvl,cond) == 1);
                if length(reversal_idx) >= all_staircases(subj).num_reversals_to_consider
                    last_reversals = all_staircases(subj).probe_offsets(n_sc,reversal_idx(end-all_staircases(subj).num_reversals_to_consider+1:end),lvl,cond);
                    cv = std(last_reversals) / mean(last_reversals);
                    mean_threshold = mean(last_reversals);
                end
            end

            % Format figure
            title(['Level ' num2str(lvl)]);
            xlabel('Trial');
            ylabel('Probe offset (°)');
            box off; 
            set(gca, 'TickDir', 'out');
            ylim([0 all_staircases(subj).max_probe_offset]);
            grid on;
            hold off;
            axis square;
        end

    end
end

%% Plot convergence analysis

for subj = 1:length(subj_IDs)
    for cond = 1:num_conds

        figure_name = ['Convergence S' subj_IDs{subj} ' ' cond_names{cond}];
        figure('Name', figure_name, 'Color', 'w');

        for lvl = 1:num_levels
            subplot(1, num_levels, lvl);
            hold on;
        
            for n_sc = 1:size(all_staircases(subj).probe_offsets, 1)
                % Get all reversal values
                reversal_idx = find(all_staircases(subj).reversal_indices(n_sc,:,lvl,cond) == 1);
                reversal_values = all_staircases(subj).probe_offsets(n_sc,reversal_idx,lvl,cond);
                
            % Plot running mean of reversals
            running_mean = movmean(reversal_values, 4);
            plot(1:length(running_mean), running_mean, colors{n_sc}, 'LineWidth', 2);
            
            % Plot individual reversals
            scatter(1:length(reversal_values), reversal_values, 50, colors{n_sc}, 'filled', 'MarkerFaceAlpha', 0.3);
            end
        
            % Format figure
            title(['Level ' num2str(lvl)]);
            text(0.5, 0.9, sprintf('Final threshold: %.1f°', all_staircases(subj).final_probe_offsets(cond,lvl)), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center');
            xlabel('Reversal Number');
            ylabel('Probe offset (°)');
            box off;
            set(gca, 'TickDir', 'out');
            ylim([0 all_staircases(subj).max_probe_offset]);
            grid on;
            hold off;
            axis square;
        
        end
    end
end

%% Performance analysis

for subj = 1:length(subj_IDs)
    for cond = 1:num_conds

        figure_name = ['Performance at threshold S' subj_IDs{subj} ' ' cond_names{cond}];
        figure('Name', figure_name, 'Color', 'w');

        for lvl = 1:num_levels
        
            subplot(1, num_levels, lvl);
            hold on;
            
            performances = nan(1, size(all_staircases(subj).probe_offsets, 1));
            num_trials = nan(1, size(all_staircases(subj).probe_offsets, 1));
            mean_deltas = nan(1, size(all_staircases(subj).probe_offsets, 1));
            
            for n_sc = 1:size(all_staircases(subj).probe_offsets, 1)
                
                % Get number of trials performed
                num_trials = sum(~isnan(all_staircases(subj).responses(n_sc,:,lvl,cond)));

                % Find reversal indices
                reversal_idx = find(all_staircases(subj).reversal_indices(n_sc,:,lvl,cond) == 1);
                            
                if length(reversal_idx) >= all_staircases(subj).num_reversals_to_consider
                    
                    % Get the last N reversals
                    last_reversals_idx = reversal_idx(end-all_staircases(subj).num_reversals_to_consider+1:end);
                    
                    % Get the trials between first last half reversals and last trial
                    trial_range = last_reversals_idx(1):num_trials;
                    
                    % Calculate mean probe offsets for these reversals
                    mean_delta = mean(all_staircases(subj).probe_offsets(n_sc, last_reversals_idx, lvl, cond));
                    mean_deltas(n_sc) = mean_delta;
                    
                    % Find trials near this probe offsets (±1%)
                    delta_window = 1 + mean_delta;
                    deltas = all_staircases(subj).probe_offsets(n_sc, trial_range, lvl, cond);
                    responses = all_staircases(subj).responses(n_sc, trial_range, lvl, cond);
                    relevant_trials = abs(deltas - mean_delta) <= delta_window;
                    
                    % Calculate and store performance
                    performances(n_sc) = mean(responses(relevant_trials), 'omitnan') * 100;
                    num_trials(n_sc) = sum(relevant_trials);
                    
                else
                    disp('  Not enough reversals');
                end

            end
        
            % Only create plot if we have valid data
            if any(~isnan(performances))
                % Plot bars
                h = scatter(1:2, performances, 50, 'filled', 'MarkerFaceAlpha', 0.8);
                % Set colors for each marker
                h.CData = [0 0 1; 1 0 0]; % Blue for SC1, Red for SC2
                
                % Add target line
                yline(70.7, 'k--', 'Alpha', 0.3);
                
                % Format figure
                title(['Level ' num2str(lvl)]);
                ylabel('Discrimination Rate (%)');
                ylim([0 100]);
                xlim([0.5 2.5]);
                xticks(1:2);
                xticklabels({'SC1', 'SC2'});
                set(gca, 'TickDir', 'out');
                axis square;
            else
                text(0.5, 0.5, 'No data (insufficient reversals)', ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
            end
            
            hold off;

        end

    end
end