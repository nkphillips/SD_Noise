function img = bandpassFilterImg(img, orientation_cutoffs, sf_cutoffs, smoothen_width, f_nyquist)

    %% Create the orientationbandpass filter

    orientation_bandpass_filter = OrientationBandpass(size(img), orientation_cutoffs(1), orientation_cutoffs(2));

    %% Create the spatial frequency bandpass filter

    spatial_frequency_bandpass_filter = Bandpass2(size(img), sf_cutoffs(1)/f_nyquist, sf_cutoffs(2)/f_nyquist);

    %% Combine the filters

    bandpass_filter = orientation_bandpass_filter .* spatial_frequency_bandpass_filter;
    bandpass_filter = imgaussfilt(bandpass_filter, smoothen_width);  % Smooth the bandpass filter
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