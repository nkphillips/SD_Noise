% run_calibration

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

p.which_setup = 0; % 0 = MacBook, 1 = 3329C_ASUS, 2 = S32D850
p.disp_on = 0;
p.half_screen = 1;
p.simulate_response = 0;

% Sync Test
if sum(p.which_setup == [0 2]) > 0
    Screen('Preference', 'SkipSyncTests', 1); % Set to 1 if running on macOS
else
    Screen('Preference', 'SkipSyncTests', 0);
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

if p.which_setup == 1
    dirs.monitor_cal_dir = '/home/serenceslabexp/Desktop/MonitorCalibration/GammaTables'; addpath(dirs.monitor_cal_dir);
end

%% Set device and display; open window

init_device_input
init_display
open_window

%% Probe offset magnitudes

p.probe_offsets = round(linspace(0,15,7));

%% Initialize trackers

exit_session = 0;

%% Define stimuli parameters

p.num_levels = 10;
p.num_reps = 30;
p.num_trials_per_feature = p.num_levels * p.num_reps;

p.contrast_min = 0.05;
p.contrast_max = 0.9;
p.calibration_contrast_levels = round(logspace(log10(p.contrast_min), log10(p.contrast_max), p.num_levels),2);

p.filter_width_min = 2;
p.filter_width_max = 80;
p.calibration_filter_width_levels = round(logspace(log10(p.filter_width_min), log10(p.filter_width_max), p.num_levels),2);

%% Calibration loop
% datasample
% randi 
% repmat 
% nan

presentation_order = nan(p.num_trials_per_feature, 2); % matrix of indices. whatever approach you use, this var should store the order for both contrast and filter width sequences
responses = nan(p.num_reps, p.num_levels, 2); % page 1: contrast, page 2: filter width

for feature = 1:2 
    
    % whatever your solution, make sure you have p.num_reps of a given lvl
    datasample(p.calibration_contrast_levels, , 'Replace', false)

end


%% Estimate psychometric function



%% Extract constrast and filter width levels



%% Save data

