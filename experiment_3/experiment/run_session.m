% run_session

%% Clean and prepare workspace

clc;
clear all; %#ok<CLALL>
close all;

commandwindow; % force cursor to command window
Priority(1); % Set MATLAB/Psychtoolbox to "high" priority level

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

toggles.which_setup = 2; % 0 = MacBook, 1 = 3329B_ASUS, 2 = 3329C_ASUS, 3 = 3329D_ASUS, 4 = S32D850
toggles.show_debug_output = 1; % Display trial information and diagnostic messages in the Command Window 
toggles.half_screen = 1; 
toggles.simulate_response = 0;
toggles.training = 0;
toggles.calibration = 0;
toggles.level_type = 1; % 1 = fixed , 0 = calibrated 
toggles.simulation_mode = 0;
toggles.demo_run = 0;


% Sync Test
if ~any(toggles.which_setup == 1:3)
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

% Verify toggles
while toggles.training && toggles.calibration
    x = input('Training and calibration are both toggled to on. For training, enter 1. For calibration, enter 2: ');
    if x == 1
        toggles.calibration = 0;
    elseif x == 2
        toggles.training = 0;
    end
end

%% Set directories

p.subj_ID = '888';

dirs.project_dir = '../'; addpath(dirs.project_dir); % Set the project directory to the parent folder and add to MATLAB's search path
dirs.script_dir = pwd; % Store current path of experiment
dirs.functions_dir = 'functions'; addpath(dirs.functions_dir); % Create functions folder; add functions folder to MATLAB's search path
dirs.data_dir = '../data'; % Create data folder located one level above
dirs.texture_dir = [dirs.script_dir '/textures']; % Create textures folder path in current folder

if exist(dirs.data_dir,'dir') == 0 % 0= does not exist
    mkdir(dirs.data_dir);
end

if exist(dirs.texture_dir,'dir') == 0 
    mkdir(dirs.texture_dir);
end

addpath(dirs.data_dir); % Add data folder to MATLAB's search path
addpath(dirs.texture_dir); % Add textures folder to MATLAB's search path

dirs.init_dir = 'init'; addpath(dirs.init_dir); % Add experiment-initialization scripts to MATLAB's search path
dirs.modules_dir = 'script_modules'; addpath(dirs.modules_dir);  % Add main experiment modules to MATLAB's search path
dirs.logs_dir = fullfile(dirs.data_dir, p.subj_ID, 'logs'); % Build current participant's log-folder path

% Add monitor gamma-calibration folder when using a lab ASUS display
if any(toggles.which_setup == 1:3)
    dirs.monitor_cal_dir = '/home/serenceslabexp/Documents/MonitorCalibration/GammaTables'; addpath(dirs.monitor_cal_dir); % Store gamma-table folder path and add to MATLAB's search path
end

%% Set device and display; open window

init_device_input % Configure the keyboard or participant response device
init_display % Define selected monitor and display properties
open_window % Open Psychtoolbox experiment window

%% Initialize trackers

exit_session = 0; % default: 0

%% Session loop

if toggles.training
    dirs.save_filename_template = ['SD_Noise_Exp2_Training_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
elseif toggles.calibration
    dirs.save_filename_template = ['SD_Noise_Exp2_Calibration_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
else
    dirs.save_filename_template = ['SD_Noise_Exp2_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
end

while ~exit_session

    % Enter experiment
    run_info = run_experiment(p, w, dirs, toggles);

    % All done screen
    if ~exit_session

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
if toggles.which_setup == 1 && w.gamma_correct
    Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
end

Priority(0);
Screen('CloseAll');
ShowCursor;
