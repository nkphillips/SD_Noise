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

stimuli.test_textures_made = nan(length(p.contrast), length(p.orientation_bp_filter_width), p.num_test_samples);
stimuli.mask_textures_made = nan(length(p.contrast), p.num_mask_samples);

for i = 1:size(stimuli.test_textures_made,1) % Contrasts
    for j = 1:size(stimuli.test_textures_made,2) % Filter width
        for k = 1:size(stimuli.test_textures_made,3) % noise sample
        
            stimuli.test_textures_made(i, j, k) = Screen('MakeTexture', w.window, stimuli.test_textures(:,:,i,j,k));

            if j == 1
                stimuli.mask_textures_made(i, k) = Screen('MakeTexture', w.window, stimuli.mask_textures(:,:,i,k));
            end
       
        end
    end
end
 
disp(['Elapsed time: ' num2str(toc) ' s'])
