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

t.iti_min = 1;
t.iti_max = 1.5;

t.iti_dur = t.iti_min + (t.iti_max - t.iti_min) .* rand(p.num_trials, 1);

t.block_dur_est = ((t.trial_dur_est + mean(t.iti_dur(1:p.num_trials_per_block)))  * p.num_trials_per_block) / 60; % s
t.rest_dur = 30;

t.exp_dur_est = (p.num_trials * t.trial_dur_est + sum(t.iti_dur) + t.rest_dur * (p.num_blocks-1)) / 60;

%% Define noise sample update rate

t.noise_sample_update_rate = 20; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

%% Define frame rate and duration

[t.frame_dur, t.frame_rate] = get_framerate(t); % s, Hz

