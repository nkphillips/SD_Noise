

function img = make_orientation_bp_filtered_img(img, orientation, filter_width, smoothen_width)

    %% Create the bandpass filter

    % Define the bandpass filter bounds
    lower_bound = orientation - filter_width/2;
    upper_bound = orientation + filter_width/2;

    % Define the bandpass filter
    bandpass_filter = OrientationBandpass(size(img), lower_bound, upper_bound);

    % Smooth the bandpass filter
    bandpass_filter = imgaussfilt(bandpass_filter, smoothen_width);

    % Normalize the bandpass filter
    bandpass_filter = normalize_array(bandpass_filter, 'min-max');

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

    % Store only the real part of the image
    img = real(img);

end