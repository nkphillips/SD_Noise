%%% Initalize timing

%{ 

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define event durations

t.cue_dur = 2; % s

t.init_adaptation_dur = 10; % s; Altan et al. 2024
t.top_up_adaptation_dur = 0.5; % s; Altan et al. 2024 

t.test_dur = 0.5; % s; Altan et al. 2024

t.response_dur = 2; % s

t.rest_dur = t.init_adaptation_dur; % s; this is the minimum rest duration
t.optional_rest_dur = nan(1, p.num_blocks);

t.trial_dur = t.test_dur + t.response_dur; % s

t.exp_dur_est = ((p.num_blocks * p.num_trials_per_block) * (t.test_dur + t.response_dur) + ...
     ((p.num_blocks * p.num_trials_per_block) - p.num_blocks) * t.top_up_adaptation_dur + ...
    ((p.num_blocks-1) * t.rest_dur)) / 60;

%% Define noise sample update rate

t.noise_sample_update_rate = 20; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

%% Define frame rate and duration

[t.frame_dur, t.frame_rate] = get_framerate(t); % s, Hz

t.gap_after_adapt = t.frame_dur * 2;

