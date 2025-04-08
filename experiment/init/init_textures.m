%%% init_textures

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Toggles

textures_filename = ['SD_Noise_textures_' p.display_setup '.mat'];

tic
generate_textures = 1; % this will be set to 0 if they don't need to be generated
if p.training
    save_textures = 0;
else
    save_textures = 1;
end

%% Aperture
% alpha level for aperture:
% 0 = completely transparent (the texture of the aperture is invisible)
% 255 = completely opaque (the texture of the aperture dominates)

aperture = create_circular_aperture(stimuli.aperture_width_px, stimuli.aperture_height_px, stimuli.aperture_radius_px); % texture size, radius of circle
% figure, subplot(1,2,1), imshow(aperture)

aperture = imgaussfilt(aperture, 0.1 * w.ppd);
% subplot(1,2,2), imshow(aperture)

% First layer of aperture is just a gray rectangle 
aperture_texture(:,:,1) = ones(size(aperture)) * w.gray;

% The second layer is the smoothened aperture made above
aperture_texture(:,:,2) = aperture * 255; 

% figure, imshow(aperture_texture(:,:,2), [0 255])

%% Generate textures 
% if needed

if generate_textures

    disp('Generating stimuli...')

    % Preallocate textures
    noise_textures = nan(stimuli.height_px, stimuli.width_px, length(stimuli.contrast), length(stimuli.bp_filter_width), p.num_test_samples);
    stimuli.test_textures = nan(stimuli.height_px, stimuli.width_px, length(stimuli.contrast), length(stimuli.bp_filter_width), p.num_test_samples);
    stimuli.mask_textures = nan(stimuli.height_px, stimuli.width_px, length(stimuli.contrast), p.num_mask_samples);

    for i = 1:size(noise_textures, 3) % Contrasts
        for j = 1:size(noise_textures, 4) % Filter widths
            for k = 1:size(noise_textures, 5) % Samples
                
                % Create base noise
                base_noise = create_noise_texture(stimuli.height_px, stimuli.width_px);
               
                % Store base noise as mask texture
                if j == 1
                    stimuli.mask_textures(:,:,i,k) = base_noise * w.gray * stimuli.contrast(i) + w.gray; % * stimuli.contrast(i)
                end
                
                % Make orientation-bandpass filtered noise
                noise_texture = make_orientation_bp_filtered_img(base_noise, 0, stimuli.bp_filter_width(j), w.ppd * 0.1);
                
                % Normalize noise texture
                noise_texture = normalize_array(noise_texture, 'min-max');
                
                % Ignore certain combos
                if i > 1 && j > 1
                    continue
                end

                % Convert to visible pixel values and scale by contrast
                stimuli.test_textures(:,:,i,j,k) = noise_texture * w.gray * stimuli.contrast(i) + w.gray; % 

            end
        end
    end

        
end

disp(['Elapsed time: ' num2str(toc) ' s'])
