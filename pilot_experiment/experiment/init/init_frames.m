%%% init_frames

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

function frames = init_frames(t,p)

% Because rest periods can vary, each frame vector is a template that is used for every block

%% Calculate number of frames required for events

frames.noise_sample_update_frames_count = round(t.frame_rate/t.noise_sample_update_rate);

frames.test_frames_count = round(t.test_dur/t.frame_dur);
frames.response_frames_count = round(t.response_dur/t.frame_dur);

%% Test

frames.test_frame_deadlines = 0:t.frame_dur:t.test_dur-t.frame_dur;

frames.test_noise_sample_update = zeros(1, frames.test_frames_count);
frames.test_noise_sample_update(1:frames.noise_sample_update_frames_count:end) = 1;

% Generate sequence of sample frames
frames.test_noise_sample_update_seq = cell(1,p.num_blocks);
num_sample_updates = sum(frames.test_noise_sample_update);
for n_block = 1:p.num_blocks
    frames.test_noise_sample_update_seq{n_block} = gen_unique_seq([p.num_trials_per_block, num_sample_updates], 1:p.num_test_samples, p.num_test_samples/2);
end


%% Rest
% Minimum rest period in between blocks

frames.rest_frame_deadlines = 0:t.frame_dur:t.rest_dur-t.frame_dur;

end