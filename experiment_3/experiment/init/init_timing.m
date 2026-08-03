%%% Initalize timing

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define event durations
% In seconds unless stated otherwise


t.test_dur = 0.5;
t.mask_dur = 0.5;
t.delay_dur = 2.0;
t.probe_dur = 0.5;
t.response_dur_est = 1;

t.trial_dur_est = sum([t.test_dur, t.mask_dur, t.delay_dur, t.response_dur_est]);

t.iti_min = 1;
t.iti_max = 1.5;

t.rest_dur = 10; % default = 10;



if toggles.training
    t.target_task_dur_mins = 5;
    t.max_task_dur_mins = 10;
    t.max_block_mins = 10;

   
elseif toggles.calibration
    t.target_task_dur_mins = 15;
    t.max_task_dur_mins = 20;
    t.max_block_mins = 10;

   
else % main exp
    t.target_task_dur_mins = 60;
    t.max_task_dur_mins = 70;
    t.max_block_mins = 10;
end

% FIND OUT HOW TO PRINT INFO

disp(['Target task duration: ', num2str(t.target_task_dur_mins), ' min'])
disp(['Maximum task duration: ', num2str(t.max_task_dur_mins), ' min'])
disp(['Maximum block duration: ', num2str(t.max_block_mins), ' min'])

%t.target_task_dur_mins = 60; % main exp = 60, training = 5, calibration = 15
%t.max_task_dur_mins = 70; % main exp = 70, training = 10, calibration = 20
%t.max_block_mins = 10; % main exp = 10, training  = 10, calibration = 10

%% Define noise sample update rate

t.noise_sample_update_rate = 10; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

%% Define frame rate and duration

[t.frame_dur, t.frame_rate] = get_framerate(t); % s, Hz

