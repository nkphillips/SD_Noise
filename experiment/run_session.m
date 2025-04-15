%%% run_session
% Skeleton script for executing multiple runs of experiment (num_runs_to_complete).

%% Clean and prepare workspace

clc;
clear all; %#ok<CLALL>
close all

commandwindow; % force cursor to command window
% Priority(1); % Set MATLAB/Psychtoolbox to "high" priority level

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.disp_on = 1;
p.half_screen = 1;
p.simulate_response = 1;
p.training = 0;

p.which_setup = 2; % 0 = MacBook, 1 = 3329C_ASUS, 2 = S32D850

% Sync Test
if sum(p.which_setup == [0 3]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

%% Set directories

p.subj_ID = '999';

dirs.project_dir = '../'; addpath(dirs.project_dir);
dirs.script_dir = pwd;
dirs.functions_dir = '../functions'; addpath(dirs.functions_dir);


dirs.data_dir = '../data'; 
dirs.texture_dir = 'textures';

if exist(dirs.data_dir,'dir') == 0
    mkdir(dirs.data_dir);
end

if exist(dirs.texture_dir,'dir') == 0
    mkdir(dirs.texture_dir);
end

addpath(dirs.data_dir);
addpath(dirs.texture_dir);

dirs.init_dir = 'init'; addpath(dirs.init_dir);
dirs.modules_dir = 'script_modules'; addpath(dirs.modules_dir);
dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/serenceslabexp/Desktop/MonitorCalibration'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window

%% Load probe offset magnitudes

if ~p.training

    save_filename = [dirs.data_dir '/' 'staircase_data_S' p.subj_ID '_*.mat'];
    
    files_found = dir([dirs.data_dir '/' p.subj_ID]);

    % Keep only the file names
    files_found = {files_found.name};

    % Keep only the file names that match staircase filename
    files_found = files_found(contains(files_found, 'staircase'));

    % Keep only the most recent file
    files_found = files_found{end};

    % Check if the file exists
    if ~isempty(files_found)

        % Load staircase data
        load([dirs.data_dir '/' p.subj_ID '/' files_found]);

        p.staircases = staircases;

        % Get final probe offsets
        p.probe_offsets = staircases.final_probe_offsets; % num_levels x num_conds

        disp(['Loaded staircase data for subject ' p.subj_ID]);

    else
        error('No staircase data found for subject %s', p.subj_ID);
    end

elseif p.training

    p.probe_offsets = 5 * ones(3, 2);

end

%% Initialize trackers

num_runs_to_complete = 1;

n_run = 1;

exit_session = 0; % default: 0

%% Session loop

if p.training
    dirs.save_filename_template = ['SD_Noise_Pilot_Training_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
else
    dirs.save_filename_template = ['SD_Noise_Pilot_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
end

while ~exit_session

    % Enter training/experiment
    if p.training
        training_info = run_experiment(p, w, dirs);
    else
        run_info = run_experiment(p, w, dirs);
    end

    % Rest between runs / All done screen
    if n_run ~= num_runs_to_complete

        % Rest text
        rest_text = ['Run ' num2str(n_run) ' of ' num2str(num_runs_to_complete) ' done! ' ...
            'When ready for the next run, press ' KbName(p.trigger_key) ' to begin.'];
        
        rest_text_boundary = Screen('TextBounds', w.window, rest_text);
        rest_text_patch = CenterRectOnPoint(rest_text_boundary, w.centerX, w.centerY);

        % Draw and flip
        Screen('DrawText', w.window, rest_text, rest_text_patch(1), rest_text_patch(2), w.white);
        Screen('Flip', w.window);

        % Let user progress to next run when ready
        KbQueueFlush(p.device_number);
        ready_for_next_run = 0;

        while ~ready_for_next_run

            % Check for response
            if ~p.simulate_response
                [key_pressed, first_press] = KbQueueCheck(p.device_number);
                which_press = find(first_press);
            else
                key_pressed = 1;
                which_press = p.trigger_key;
            end
            
            % Evaluate response
            if key_pressed && ~isempty(which_press)
                if which_press(1) == p.keypress_numbers(end) 
                    % escape session
                    exit_session = 1;
                elseif which_press(1) == p.trigger_key
                    % update flags for progression
                    n_run = n_run + 1;
                    ready_for_next_run = 1;
                end
            end

        end


    else

        exit_session = 1;

        all_done_text = 'Session complete!';
        
        all_done_text_boundary = Screen('TextBounds', w.window, all_done_text);
        all_done_text_patch = CenterRectOnPoint(all_done_text_boundary, w.centerX, w.centerY);

        Screen('DrawText', w.window, all_done_text, all_done_text_patch(1), all_done_text_patch(2), w.white);
        Screen('Flip', w.window);
        
        WaitSecs(2);

    end

end

%% Turn off Kb and restore display

KbQueueStop(p.device_number);

Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
if p.which_setup == 1 && w.gamma_correct, SetScreenDefault; end
Screen('CloseAll');
ShowCursor;