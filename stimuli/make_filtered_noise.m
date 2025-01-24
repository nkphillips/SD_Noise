%%% make_filtered_noise
% Contrast normalization inspired by Perfetto, Wilder, & Walther 2020 Vision

%{

Written by Luis D. Ramirez
UCSD
lur003@ucsd.edu

%}

tic
%%

noise_made = nan(p.num_orientations,p.num_noise_samples);

for n_sample = 1:p.num_noise_samples

    %% Create base noise texture

    base_noise = create_noise_texture(stimuli);

    %% Create bandpass filtered versions and z-score images
    % Steps inspired by Perfetto, Wilder, & Walther 2020 Vision

    noise_images_z = nan(stimuli.diameter_px, stimuli.diameter_px, p.num_orientations);

    for n = 1:p.num_orientations 

        noise_image = make_orientation_bp_filtered_img(base_noise, p.orientations(n), stimuli);

        noise_image_z = normalize_array(noise_image,'z-score');

        noise_image_z(noise_image_z < -2) = -2;
        noise_image_z(noise_image_z > 2) = 2;

        noise_images_z(:,:,n) = noise_image_z ./ std2(noise_image_z); 

        %     imshow(noise_image_z,[-2 2])

    end

    %% Contrast normalization
    % Jointly scale images linearly to visible range [0,255]

    % Get joint min and max
    z_min = min(noise_images_z(:));
    z_max = max(noise_images_z(:));

    % min-max normalize images and scale to visible range
    noise_images_scaled = (noise_images_z - z_min) ./ (z_max - z_min) * 255;


    %{
for n = 1:p.num_orientations 

    imshow(noise_images_scaled(:,:,n),[0 255])

end
    %}

    %% Apply desired contrast

    stimuli.probe_noise = noise_images_scaled * stimuli.contrast;

    %{
for n = 1:p.num_orientations 

    imshow(probe_noise(:,:,n),[0 255])

end
    %}

    %% Make stimuli
    % Requires open window

    for n = 1:p.num_orientations 

        noise_made(n,n_sample) = Screen('MakeTexture', window, stimuli.probe_noise(:,:,n));

    end

end

disp(['Noise generation took ~' num2str(round(toc)) ' (s)'])
