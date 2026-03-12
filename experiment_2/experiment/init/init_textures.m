%%% init_textures

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

function stimuli = init_textures(p, dirs, w)

tic

%% Toggles

    if p.training
        textures_filename = ['SD_Noise_Exp2_Training_textures_' p.display_setup '.mat'];
    elseif p.calibration
        textures_filename = ['SD_Noise_Exp2_Calibration_textures_' p.display_setup '.mat'];
    else
        textures_filename = ['SD_Noise_Exp2_S' p.subj_ID '_textures_' p.display_setup '.mat'];
    end

textures_path = [dirs.texture_dir '/' textures_filename];

if ~exist(textures_path, 'file')
    stimuli.generate_textures = 1;
    stimuli.save_textures = 1;
else
    stimuli.generate_textures = 0;
    stimuli.save_textures = 0;
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

stimuli.aperture_texture = aperture_texture;

% figure, imshow(aperture_texture(:,:,2), [0 255])

%% Generate textures

if stimuli.generate_textures

    disp('Generating stimuli...')

    % Preallocate textures
    noise_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_noise_samples);
    stimuli.test_textures = nan(p.height_px, p.width_px, length(p.contrast), length(p.orientation_bp_filter_width), p.num_noise_samples);
    stimuli.mask_textures = nan(p.height_px, p.width_px, length(p.contrast), p.num_noise_samples);

    for i = 1:size(noise_textures, 3) % Contrasts
        for j = 1:size(noise_textures, 4) % Orientation filter widths
            for k = 1:size(noise_textures, 5) % Samples

                % Ignore certain combos
                % We need all contrasts (i) at the easiest filter width (j=3, since 85% is the 3rd element)
                % We need all filter widths (j) at the easiest contrast (i=3, since 85% is the 3rd element)
                % So if we're not at the easiest filter width AND we're not at the easiest contrast, skip.

                if p.training
                    % In training, there is only 1 level for contrast and filter width,
                    % so we just generate the single combination.
                    if i > 1 || j > 1
                        continue
                    end
                elseif p.calibration
                    % In calibration, the easiest filter width is index 1 (2 degrees),
                    % and easiest contrast is index p.num_levels (90% contrast)
                    if i < length(p.contrast) && j > 1
                        continue
                    end
                else
                    % In main experiment, the easiest contrast is index p.num_levels (85% performance).
                    % The easiest filter width is index 1 (narrowest filter, sorted ascending).
                    if i ~= p.num_levels && j ~= 1
                        continue
                    end
                end

                % Create base noise (used as a mask)
                base_noise = create_noise_texture(p.height_px, p.width_px);
                base_noise = bandpassFilterImg(base_noise, [0, 180], [0.5 6], w.ppd * 0.1, w.f_Nyquist);
                base_noise = centerTextureContrast(base_noise, p.contrast(i), w.gray);

                % masks should be generated using the easiest baseline filter width (j=1, narrowest filter, sorted ascending)
                if j == 1
                    stimuli.mask_textures(:,:,i,k) = base_noise;
                end
                
                % Make orientation- and spatial frequency-bandpass filtered noise
                noise_texture = bandpassFilterImg(base_noise, [round(90 - p.orientation_bp_filter_width(j)/2), floor(90 + p.orientation_bp_filter_width(j)/2)], p.sf_bp_filter_cutoffs, w.ppd * 0.1, w.f_Nyquist);
                noise_texture = centerTextureContrast(noise_texture, p.contrast(i), w.gray);

                stimuli.test_textures(:,:,i,j,k) = noise_texture; % Convert to visible pixel values and scale by contrast

            end
        end
    end

    %% Save textures

    if stimuli.save_textures
        save(textures_path, 'stimuli', '-v7.3');
    end

else

    load(textures_path, 'stimuli');

end

%% Probe line

probe_line = ones(p.probe_length) * w.gray;

start_col = round(p.probe_length/2) - floor(p.probe_thickness/2);
end_col = round(p.probe_length/2) + floor(p.probe_thickness/2);
probe_line(:, start_col:end_col) = 0;

stimuli.probe_line = probe_line;

%%

disp(['Elapsed time: ' num2str(toc) ' s'])

end