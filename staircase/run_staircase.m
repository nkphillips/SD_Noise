%%% run_staircase

%% 
close all; clear all; clc

commandwindow;

% Grab date
t.the_date = datestr(now, 'yyyymmdd'); % Grab today's date
t.the_time = datestr(now,'HHMM'); % Grab current time

% Generate unique seed for random number generator (rng)
t.my_rng_seed = sum(100*clock);
rng(t.my_rng_seed);

%% Toggles

p.disp_on = 1;
p.half_screen = 1;
p.simulate_response = 0;
p.training = 0;

p.which_setup = 1; % 0 = MacBook, 1 = 3329C_ASUS, 2 = S32D850

% Sync Test
if sum(p.which_setup == [0 3]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
end 

save_staircase_data = 1;

%% Set directories

p.subj_ID = '999'; % dummy subj = 999

dirs.script_dir = pwd;
dirs.functions_dir = '../functions'; addpath(dirs.functions_dir);
dirs.data_dir = '../data'; addpath(dirs.data_dir);
dirs.init_dir = '../experiment/init'; addpath(dirs.init_dir);
dirs.modules_dir = '../experiment/script_modules'; addpath(dirs.modules_dir);
dirs.texture_dir = '../experiment/textures'; addpath(dirs.texture_dir);

dirs.logs_dir = [dirs.data_dir '/' p.subj_ID '/logs'];

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/serenceslabexp/Desktop/MonitorCalibration'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window
init_fixation

%% Loading screen

loading_text = 'Loading calibration...';

loading_text_boundary = Screen('TextBounds', w.window, loading_text);
loading_text_patch = CenterRectOnPoint(loading_text_boundary, w.centerX, w.centerY);

Screen('DrawText', w.window, loading_text, loading_text_patch(1),  loading_text_patch(2), w.white);
Screen('Flip', w.window);

%% Define stimulus

init_stimuli_params

%% Initialize staircase

init_staircase

%% Initialize (create) textures

init_textures

%% Make textures for drawing

make_textures

%% Define patches for drawing

make_patches

%% Define timing

init_timing

%% Define frames

frames = init_frames(t,p);

%% Loading complete

wait_for_trigger

%% Staircase 

% Run probe offset staircases
discrimination_staircase

%% Calculate final probe offsets for each condition & level

staircases.final_probe_offsets = nan(p.num_conds, p.num_levels);
for cond = 1:p.num_conds
    for lvl = 1:p.num_levels

        % Pre-allocate vector to store the mean probe offset for a staircase
        reversal_probe_offsets = nan(staircases.num_staircases_per_cond, 1);

        % Get the contrast levels at all the reversals from both staircases
        for n_sc = 1:staircases.num_staircases_per_cond

            % Get the number of reversals
            curr_num_reversals = nnz(staircases.reversal_indices(n_sc,:,lvl,cond) == 1);

            % Get the contrast levels at the reversals
            if curr_num_reversals > 0
                tmp_reversal_probe_offsets = staircases.probe_offsets(n_sc, staircases.reversal_indices(n_sc,:,lvl,cond) == 1, lvl,cond);
            end

            % Ignore the first reversals with respect to num_reversals_to_consider (counting back from the last reversal)
            tmp_reversal_probe_offsets = tmp_reversal_probe_offsets(end-staircases.num_reversals_to_consider+1:end);

            % Average all the reversals across a staircase
            reversal_probe_offsets(n_sc) = mean(tmp_reversal_probe_offsets, 'omitnan');

        end

        % Average the two levels across staircases
        staircases.final_probe_offsets(cond,lvl) = mean(reversal_probe_offsets);

    end
end

%% Turn off Kb and restore display

KbQueueStop(p.device_number);

Screen('LoadNormalizedGammaTable', w.window, w.DefaultCLUT);
if p.which_setup == 1, SetScreenDefault; end
Screen('CloseAll');
ShowCursor;

%% Save staircase data

staircases.contrast = stimuli.contrast;
staircases.bp_filter_width = stimuli.bp_filter_width;
staircases.p = p;
staircases.t = t;

if save_staircase_data

    save_filename = ['staircase_data_S' p.subj_ID '_' p.display_setup '_' t.the_date '_' t.the_time '.mat'];
    save([dirs.data_dir '/' p.subj_ID '/' save_filename], 'staircases');
    disp(['Staircase data saved to ' dirs.data_dir '/' p.subj_ID '/' save_filename]);
    
end

