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

p.demo_run = 1; % halves the screen
p.simulate_response = 1;
p.training = 1;

p.which_setup = 3; % 0 = MacBook, 1 = 3329D, 2 = Scanner, 3 = S32D850

% Sync Test
if sum(p.which_setup == [0 3]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end

%% Set directories

p.subj_ID = '999';

dirs.script_dir = pwd;
dirs.init_dir = 'init'; addpath(dirs.init_dir);
dirs.modules_dir = 'script_modules'; addpath(dirs.modules_dir);
dirs.functions_dir = '../../functions'; addpath(dirs.functions_dir);
dirs.texture_dir = 'textures'; addpath(dirs.texture_dir);
dirs.data_dir = '../data'; addpath(dirs.data_dir);
dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/pclexp/Documents/Luis/MonitorCalibration'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window

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

    % Enter experiment
    run_info = run_experiment(p, w, dirs);

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
if p.which_setup == 1, SetScreenDefault; end
Screen('CloseAll');
ShowCursor;