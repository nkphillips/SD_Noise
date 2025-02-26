% evaluate_staircase

clear all; close all; clc;

%% Prepare workspace

subj_IDs = {'000'};

script_dir = pwd;
staircase_data_dir = '../data'; addpath(staircase_data_dir);

%% Load staircase data

% Prepare a structure to store all the staircase data from each subject
all_staircases = {};

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

%% Plot staircases

for subj = 1:length(subj_IDs)

    % Subplots of all contrast deltas per SF
    figure_name = ['Staircase Trials S' subj_IDs{subj}];
    figure('Name', figure_name, 'Color', 'w');

    for n_sf = 1:size(all_staircases(subj).contrast_deltas, 3)

        subplot(size(all_staircases(subj).contrast_deltas, 3), 1, n_sf);
        plot(100*all_staircases(subj).contrast_deltas(:,:,n_sf)');

        % Format figure
        % title(['SF ' num2str(all_staircases(subj).inducer_sfs(n_sf))]);
        if n_sf == 3, xlabel('Trial'); end
        if n_sf == 1, ylabel('Contrast delta (%)'); end
        box off; set(gca, 'TickDir', 'out');

    end

    % Plot final thresholds
    figure_name = ['Final Thresholds S' subj_IDs{subj}];
    figure('Name', figure_name, 'Color', 'w');
    bar(100*all_staircases(subj).final_contrast_deltas);
    
    % Format figure
    xlabel('Inducer SF');
    ylabel('Final contrast delta (%)');

    % xticklabels(num2str(round(all_staircases(subj).inducer_sfs, 2)));

    box off;set(gca, 'TickDir', 'out');

end