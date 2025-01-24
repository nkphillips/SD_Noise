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

% Adaptor stimuli
stimuli.adaptor_sf_made = nan(sum(~isnan(p.adaptor_sfs)), p.num_adaptor_samples);

for i = 1:size(stimuli.adaptor_sf_made,1)
    for j = 1:size(stimuli.adaptor_sf_made,2)

        stimuli.adaptor_sf_made(i,j) = Screen('MakeTexture', w.window, stimuli.adaptor_sf_textures(:,:,i,j));

    end
end

% Reference stimuli
stimuli.reference_sf_made = nan(p.num_reference_samples,1);

for i = 1:p.num_reference_samples
    
    stimuli.reference_sf_made(i) = Screen('MakeTexture', w.window, stimuli.reference_sf_textures(:,:,i));

end

% Test stimuli
stimuli.test_sf_made = nan(p.num_test_sfs, p.num_test_samples);

for i = 1:size(stimuli.test_sf_made,1)
    for j = 1:size(stimuli.test_sf_made,2)
        
        stimuli.test_sf_made(i,j) = Screen('MakeTexture', w.window, stimuli.test_sf_textures(:,:,i,j));

    end
end

%% 

disp(['Elapsed time: ' num2str(toc) ' s'])
