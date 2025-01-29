%%% Initalize timing

%{ 

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define event durations

t.test_dur = 0.5; % s;
t.mask_dur = t.test_dur; % s
t.delay_dur = t.test_dur;

t.response_dur_est = 2; % s

t.trial_dur_est = sum([t.test_dur, t.mask_dur, t.delay_dur, t.response_dur_est]); % s
t.block_dur = 600; % s

t.exp_dur_est = [] / 60;

%% Define noise sample update rate

t.noise_sample_update_rate = 20; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

%% Define frame rate and duration

[t.frame_dur, t.frame_rate] = get_framerate(t); % s, Hz

