%%% init_textures

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

%% Check if textures exist

textures_filename = ['SD_Noise_Pilot_S' p.subj_ID '_textures_' p.display_setup '_' num2str(p.run_num) '.mat'];

tic
generate_textures = 1; % this will be set to 0 if they don't need to be generated
save_textures = 1;

if p.run_num > 1

    cd(dirs.texture_dir)

    % Check if subject folder exists
    if ~exist(p.subj_ID,'dir')

        mkdir(p.subj_ID)

    elseif exist(p.subj_ID,'dir') ~= 0

        % Check if textures file exists
        cd(p.subj_ID)
        if exist(textures_filename,'file') ~= 0

            disp('Loading stimuli...')
            load(textures_filename)

            % Check compatibility
            struct_fieldnames = fieldnames(textures);

            incompat = zeros(1,numel(struct_fieldnames));

            for n_field = 1:numel(struct_fieldnames)
                
                incompat_test = stimuli.(struct_fieldnames{n_field}) == textures.(struct_fieldnames{n_field});
                
                if sum(incompat_test) < length(incompat_test)
                    incompat(n_field) = 1;
                end

            end

            % If the loaded textures are compatible, append them to the stimuli struct
            if sum(incompat) == 0

                disp('Loaded textures are compatible!')
                generate_textures = 0;
                save_textures = 0;

                struct_fieldnames = fieldnames(textures);
                struct_fieldnames(~contains(struct_fieldnames,'textures')) = [];

                for n_field = 1:numel(struct_fieldnames)
                
                    stimuli.(struct_fieldnames{n_field}) = textures.(struct_fieldnames{n_field});

                end

                clear textures

            end

        end
        
        cd(dirs.script_dir)

    end

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

    noise_textures = nan(stimuli.noise_height_px, stimuli.noise_width_px, length(stimuli.contrast), length(stimuli.bp_filter_width), p.num_noise_samples);

    for i = 1:size(noise_textures, 3) % Contrasts
        for j = 1:size(noise_textures, 4) % Filter widths
            for k = 1:size(noise_textures, 5) % Samples
                
                base_noise = create_noise_texture(stimuli.noise_height_px, stimuli.noise_width_px);
                
                noise_texture = make_orientation_bp_filtered_img(base_noise, p.orientations(n), stimuli);
                
                noise_texture = normalize_array(noise_texture, 'z-score');
                noise_texture(noise_texture < -2) = -2;
                noise_texture(noise_texture > 2) = 2;
                
                noise_textures(:,:,i,j,k) = noise_texture / std2(noise_texture);
                
                % Contrast normalization (jointly scale images linearly to visible range 0:255
                z_min = min(noise_textures(:));
                z_max = max(noise_textures(:));
                
                stimuli.textures(:,:,i,j,k) = ((noise_textures - z_min) / (z_max - z_min) ) * 255 * stimuli.test_sf_contrast;
                
            end
        end
    end

end

disp(['Elapsed time: ' num2str(toc) ' s'])
