%%% init_textures

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

tic

%% Toggles

textures_filename = ['SD_Noise_textures_' p.display_setup '.mat'];

if ~exist([dirs.texture_dir '/' textures_filename], 'file')
    generate_textures = 1;
else
    generate_textures = 0;
end

if p.training
    save_textures = 0;
else
    save_textures = 1;
end

%% Aperture
% alpha level for aperture:
% 0 = completely transparent (the texture of the aperture is invisible)
% 255 = completely opaque (the texture of the aperture dominates)

aperture = create_circular_aperture(p.aperture_width_px, p.aperture_height_px, p.aperture_radius_px); % texture size, radius of circle
% figure, subplot(1,2,1), imshow(aperture)

aperture = imgaussfilt(aperture, 0.1 * w.ppd);
% subplot(1,2,2), imshow(aperture)
aperture_texture(:,:,1) = ones(size(aperture)) * w.gray;
aperture_texture(:,:,2) = aperture * 255; 

% figure, imshow(aperture_texture(:,:,2), [0 255])

%% Generate textures 

if generate_textures

    disp('Generating stimuli...')

    % Preallocate textures
    noise_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_test_samples);
    p.test_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_test_samples);
    p.mask_textures = nan(p.height_px, p.width_px, length(p.contrast), p.num_mask_samples);

    for i = 1:size(noise_textures, 3) % Contrasts
        for j = 1:size(noise_textures, 4) % Filter widths
            for k = 1:size(noise_textures, 5) % Samples
                
                base_noise = create_noise_texture(p.height_px, p.width_px);
                base_noise = bandpassFilterImg(base_noise, [0, 180], [0.5 6], w.ppd * 0.1, w.f_Nyquist);
                base_noise = centerTextureContrast(base_noise, p.contrast(i), w.gray);
                
                if j == 1
                    stimuli.mask_textures(:,:,i,k) = base_noise;
                end
                
                % Make orientation- and spatial frequency-bandpass filtered noise
                % noise_texture = make_orientation_bp_filtered_img(base_noise, 0, p.orientation_bp_filter_width(j), w.ppd * 0.1);

                noise_texture = bandpassFilterImg(base_noise, [180 - p.orientation_bp_filter_width(j), 180 + p.orientation_bp_filter_width(j)], p.sf_bp_filter_cutoffs, w.ppd * 0.1, w.f_Nyquist);
                noise_texture = centerTextureContrast(noise_texture, p.contrast(i), w.gray);

                % Ignore certain combos
                if i > 1 && j > 1
                    continue
                end

                stimuli.test_textures(:,:,i,j,k) = noise_texture; % Convert to visible pixel values and scale by contrast

            end
        end
    end

else
    load([dirs.texture_dir '/' textures_filename]);
    stimuli = textures;
    clear textures;

end

disp(['Elapsed time: ' num2str(toc) ' s'])
