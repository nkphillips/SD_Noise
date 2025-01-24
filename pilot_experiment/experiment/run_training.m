%%% run_training
% Phase 1: Adjustment of SF
% Phase 2: Practice main task

%% Clean and prepare workspace

clc;
clear all; %#ok<CLALL> 
close all

commandwindow; % force cursor to command window
Priority(1); % Set MATLAB/Psychtoolbox to "high" priority level

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.demo_run = 0;
p.simulate_response = 0;

p.which_setup = 1; % 0 = MacBook, 1 = 3329D, 2 = Scanner, 3 = S32D850

% Sync Test
if sum(p.which_setup == [0 3]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

%% Set directories

p.subj_ID = '003';

dirs.script_dir = pwd;
dirs.init_dir = 'init'; addpath(dirs.init_dir);
dirs.modules_dir = 'script_modules'; addpath(dirs.modules_dir);
dirs.functions_dir = '../../functions'; addpath(dirs.functions_dir);
dirs.texture_dir = 'textures'; addpath(dirs.texture_dir);
dirs.data_dir = '../data'; addpath(dirs.data_dir);

dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];
dirs.training_dir = [dirs.data_dir '/' p.subj_ID '/training'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/pclexp/Documents/Luis/MonitorCalibration'; addpath(dirs.monitor_cal_dir);
end

%% Initialize device input

% Enable unified mode of KbName, so KbName accepts identical key names on all operating systems:
KbName('UnifyKeyNames');

p.device_number = 0;
[Kb_indices, product_names, ~] = GetKeyboardIndices;

switch p.which_setup

    case 0 % If using Macbook

        p.device_string = 'Apple Internal Keyboard / Trackpad'; % Macbook
        %     p.device_string =  'USB-HID Keyboard'; % external keyboard

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@') KbName('Space') KbName('Escape')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

    case 1 % If using 3329D

        p.device_string = 'Dell Dell USB Keyboard';

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@') KbName('Space') KbName('Escape')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

    case 3 % If using S32D850

        p.device_string =  'USB-HID Keyboard'; % external keyboard

        % Define keypress numbers
        p.keypress_numbers = [KbName('1!') KbName('2@') KbName('Space') KbName('Escape')];

        % Assign trigger key
        p.trigger_key = KbName('Space');

end

% Scan for device number
for i = 1:length(product_names)
    if strcmp(product_names{i}, p.device_string)
        p.device_number = Kb_indices(i);
        break;
    end
end

% Error if no device is found
if p.device_number == 0
    error('No device by that name was detected');
end

%% Set display settings and open window

init_display
init_fixation
open_window

% Define fixation space and dot patches (how and where the fixation space and dot will be drawn)
created_fixation_space = Screen('MakeTexture', w.window, fixation_space);
fixation_dot_patch = CenterRectOnPoint([0 0 p.fixation_dot_px p.fixation_dot_px], w.centerX, w.centerY);
fixation_space_patch = CenterRectOnPoint([0 0 p.fixation_space_px p.fixation_space_px], w.centerX, w.centerY);

%% Loading screen for phase 1

loading_text = 'Loading adjustment ...'; % standby text for loading experiment and stimuli
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Phase 1 - Explore SF

% Get start timestamp
t.explore_sf_start = GetSecs; % seconds

% Step into adjustment script
adjust_sf

% Get end timestamp
t.explore_sf_end = GetSecs; % seconds
t.explore_sf_dur = t.explore_sf_end - t.explore_sf_start; % seconds

%% Pause before phase 2

while KbCheck; end
WaitSecs(0.5);

% Rest text
rest_text = ['When ready for practice, press ' KbName(p.trigger_key) ' to begin.'];

rest_text_boundary = Screen('TextBounds', w.window, rest_text);
rest_text_patch = CenterRectOnPoint(rest_text_boundary, w.centerX, w.centerY);

% Draw and flip
Screen('DrawText', w.window, rest_text, rest_text_patch(1), rest_text_patch(2), w.white);
Screen('Flip', w.window);

% Let user progress to next run when ready
ready_for_practice = 0;

while ~ready_for_practice

    % Check for response
    if ~p.simulate_response

        [key_pressed, ~, first_press] = KbCheck(p.device_number);
        which_press = find(first_press);

    else

        key_pressed = 1;
        which_press = p.keypress_numbers(3);

    end

    % Evaluate response
    if key_pressed & ~isempty(which_press)
        
        disp(['Key Pressed: ', num2str(key_pressed), ' Which Press: ', num2str(which_press)]);

        if which_press(1) == p.keypress_numbers(end)
            
            % escape session
            exit_session = 1;
       
        elseif which_press(1) == p.trigger_key
            
            % update flags for progression
            ready_for_practice = 1;
       
        end
    end

end


%% Phase 2 - Practice Task

while KbCheck; end
WaitSecs(0.5);
KbEventFlush(p.device_number);
KbQueueCreate(p.device_number);
KbQueueStart(p.device_number);

p.run_num = 1;
p.fixation_dot_color = w.black;

init_practice_experiment
practice_experiment

%% Training complete!

loading_text = 'Training Complete!';
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);
WaitSecs(1);

%% Save training info

loading_text = 'Saving training info...';
loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

training_info.behav_data = behav_data;
training_info.p = p;
training_info.t = t;
training_info.w = w;
training_info.frames = frames;
training_info.pres_timing = pres_timing;
training_info.exe_timing = exe_timing;

cd(dirs.data_dir)
if ~exist(p.subj_ID,'dir')
    mkdir(p.subj_ID)
    disp(['/' p.subj_ID ' created.'])
end

save_filename = ['SD_SF_Apapt_Pilot_S' p.subj_ID '_Training_' p.display_setup '.mat'];
save([p.subj_ID '/' save_filename],'training_info','-mat','-v7.3');

cd(dirs.script_dir)

%% Turn off Kb and restore display

Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
if p.which_setup == 1, SetScreenDefault; end
Screen('CloseAll');
ShowCursor;