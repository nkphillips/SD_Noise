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

p.which_setup = 3; % 0 = MacBook, 1 = 3329B_ASUS, 2 = 3329C_ASUS, 3 = 3329D_ASUS, 4 = S32D850
p.disp_on = 1;
p.half_screen = 1;
p.simulate_response = 0;
p.training = 0;
p.calibration = 1;

% Sync Test
if ~any(p.which_setup == 1:3)
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

% Verify toggles
while p.training && p.calibration
    x = input('Training and calibration are both toggled to on. For training, enter 1. For calibration, enter 2: ');
    if x == 1
        p.calibration = 0;
    elseif x == 2
        p.training = 0;
    end
end

%% Set directories

p.subj_ID = '999';

dirs.project_dir = '../'; addpath(dirs.project_dir);
dirs.script_dir = pwd;
dirs.functions_dir = 'functions'; addpath(dirs.functions_dir);
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

if any(p.which_setup == 1:3)
    dirs.monitor_cal_dir = '/home/serenceslabexp/Documents/MonitorCalibration/GammaTables'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window

%% Initialize trackers

exit_session = 0; % default: 0

%% Session loop

if p.training
    dirs.save_filename_template = ['SD_Noise_Exp2_Training_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
elseif p.calibration
    dirs.save_filename_template = ['SD_Noise_Exp2_Calibration_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
else
    dirs.save_filename_template = ['SD_Noise_Exp2_S' p.subj_ID '_Run*_' p.display_setup '.mat'];
end

while ~exit_session

    % Enter experiment
    run_info = run_experiment(p, w, dirs);

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
if p.which_setup == 1 && w.gamma_correct
    Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
end

Priority(0);
Screen('CloseAll');
ShowCursor;
