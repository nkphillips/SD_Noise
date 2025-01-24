% Make SF bandpass-filtered noise mask texture

%{

Reads display setup (in 'w' struct), patch and filter parameters (in 'stimuli' struct) 
to create textures

%}

function mask_made = make_mask(img, stimuli, w, p, window)

%% Define lower and upper cutoff frequencies for bandpass filter

lower_cutoff = (p.min_sf - stimuli.bp_filter_width) / w.f_Nyquist;
upper_cutoff = (p.max_sf + stimuli.bp_filter_width) / w.f_Nyquist;

%% Create SF bandpass filter

bp_filter = Bandpass2(size(img), lower_cutoff, upper_cutoff);

bp_filter_smooth = imgaussfilt(bp_filter, stimuli.bp_filter_gauss_sd);

 % Normalize filter so values are between 0 and 1
bp_filter_norm = normalize_array(bp_filter_smooth, 'min-max');

% Set final filter 
% final_bp_filter = bp_filter_smooth;
final_bp_filter = bp_filter_norm;

%% FFT image

fft_img = fft2(double(img)); 

%% Store and remove DC value

% dc_comp = fft_img(1,1); % Store DC value
% 
% fft_img(1,1) = 0; % Eliminate DC value

%% Shift FFT

fft_shifted_img = fftshift(fft_img); % Moves lowest frequencies to the center

%% Store amplitude and phase info

amplitude = abs(fft_shifted_img);
phase = angle(fft_shifted_img);

%% Apply bandpass filter to amplitude of FFT image

filtered_amplitude = final_bp_filter .* amplitude;

%% Recombine filtered amplitude with original phase spectrum

filtered_fft_img = filtered_amplitude .* exp(1i*phase);

%% Inverse FFT shift

filtered_fft_img = ifftshift(filtered_fft_img); 

%% Re-introduce original DC value

% filtered_fft_img(1,1) = dc_comp; 

%% Inverse FFT to bring back to the spatial domain

filtered_img = ifft2(filtered_fft_img);

%% Store only the magnitude of the image

filtered_img = real(filtered_img);

%% z-score image

mask_img_z = normalize_array(filtered_img,'z-score');

mask_img_z(mask_img_z < -2) = -2;
mask_img_z(mask_img_z > 2) = 2;

mask_img_z = mask_img_z ./ std2(mask_img_z);

%% Contrast normalization
% Jointly scale images linearly to visible range [0,255]

mask_img = ( mask_img_z - min(mask_img_z(:)) ) ./ ( max(mask_img_z(:)) - min(mask_img_z(:)) ) * 255;

%% Apply desired contrast

mask_img = mask_img * stimuli.noise_contrast;

%% Make stimuli

mask_made = Screen('MakeTexture', window, mask_img);

end