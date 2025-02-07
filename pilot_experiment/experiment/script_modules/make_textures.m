%%% make_textures
% Make textures for drawing
% Window must be open for this to work (e.g., Screen('OpenWindow', ...)

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

disp('Making textures...')
tic

%% Make fixation space

fixation_space_made = Screen('MakeTexture', w.window, fixation_space);

%% Make wedge aperture

stimuli.aperture_made = Screen('MakeTexture', w.window, aperture_texture);

%% Make noise stimuli

stimuli.textures_made = nan(length(stimuli.contrast), length(stimuli.bp_filter_width), p.num_noise_samples);

for i = 1:size(stimuli.textures_made,1) % Contrasts
    for j = 1:size(stimuli.textures_made,2) % Filter width
        for k = 1:size(stimuli.textures_made,3) % noise sample
        
            stimuli.textures_made(i, j, k) = Screen('MakeTexture', w.window, stimuli.textures(:,:,i,j,k));
       
        end
    end
end
 
disp(['Elapsed time: ' num2str(toc) ' s'])
