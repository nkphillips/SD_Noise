%%% Initalize timing

%{ 

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Define event durations
% In seconds unless stated otherwise

if p.longer_stim_dur
    t.test_dur = 3; 
else
    t.test_dur = 0.5; 
end

t.mask_dur = 0.5; 
t.delay_dur = 2; 
t.probe_dur = 0.5;
t.response_dur_est = 2; 

t.trial_dur_est = sum([t.test_dur, t.mask_dur, t.delay_dur, t.response_dur_est]); 

t.iti_min = 1;
t.iti_max = 1.5;

t.iti_dur = t.iti_min + (t.iti_max - t.iti_min) .* rand(p.num_trials-1, 1);

t.block_dur_est = round(((t.trial_dur_est + mean(t.iti_dur(1:p.num_trials_per_block)))  * p.num_trials_per_block) / 60, 2); 
t.rest_dur = 10; % default = 10;

t.exp_dur_est = round((p.num_trials * t.trial_dur_est + sum(t.iti_dur) + t.rest_dur * (p.num_blocks-1)) / 60, 2);

%% Define noise sample update rate

t.noise_sample_update_rate = 10; % Hz; default = 20
t.noise_sample_dur = 1/t.noise_sample_update_rate; % s

%% Define frame rate and duration

[t.frame_dur, t.frame_rate] = get_framerate(t); % s, Hz

