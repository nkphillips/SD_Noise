function img = bandpassFilterImg(img, orientation_cutoffs, sf_cutoffs, smoothen_width, f_nyquist)

[rows, cols] = size(img);

%% Create the spatial frequency bandpass filter

spatial_frequency_bandpass_filter = Bandpass2(size(img), sf_cutoffs(1)/f_nyquist, sf_cutoffs(2)/f_nyquist);

% Smooth the SF filter to prevent spatial ringing (circularly symmetric, no angular distortion)
spatial_frequency_bandpass_filter = imgaussfilt(spatial_frequency_bandpass_filter, smoothen_width);

%% Create the orientation bandpass filter natively in polar coordinates
% This prevents the Cartesian blurring artifact at low spatial frequencies

u = (-cols/2 : cols/2-1);
v = (-rows/2 : rows/2-1);
[U, V] = meshgrid(u, v);
theta_2d = atan2d(V, U);
theta_2d(theta_2d < 0) = theta_2d(theta_2d < 0) + 180;

center_angle = mean(orientation_cutoffs);
half_width = abs(orientation_cutoffs(2) - orientation_cutoffs(1)) / 2;

ang_dist = abs(theta_2d - center_angle);
ang_dist = min(ang_dist, 180 - ang_dist); % Handle 180-deg wrapping

% Butterworth-style angular filter for a flat top with smooth roll-off
filter_order = 4; % Smooth but firm cutoff

if half_width == 0
    orientation_bandpass_filter = ones(rows, cols);
else
    orientation_bandpass_filter = 1 ./ (1 + (ang_dist ./ half_width).^(2 * filter_order));
end

%% Combine the filters

bandpass_filter = orientation_bandpass_filter .* spatial_frequency_bandpass_filter;
bandpass_filter = normalize_array(bandpass_filter, 'min-max');   % Normalize the bandpass filter

%% FFT the image

img_fft = fft2(img);

%% Shift the FFT

img_fft = fftshift(img_fft);

%% Apply the bandpass filter

img_fft = img_fft .* bandpass_filter;

%% Inverse FFT shift

img_fft = ifftshift(img_fft);

%% Inverse FFT to bring back to spatial domain

img = ifft2(img_fft);
img = real(img);

end