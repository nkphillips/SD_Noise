function img = bandpassFilterImgVonMises(img, orientation_cutoffs, sf_cutoffs, smoothen_width, f_nyquist)
    %% Create the orientation bandpass filter using von Mises function
    
    % Get image dimensions
    [M, N] = size(img);
    
    % Create meshgrid for frequency coordinates
    [u, v] = meshgrid(-floor(N/2):ceil(N/2)-1, -floor(M/2):ceil(M/2)-1);
    
    % Convert to polar coordinates
    [theta, rho] = cart2pol(u, v);
    
    % Convert to degrees
    theta = theta * 180/pi;
    
    % Adjust theta to be in [0, 180] range (orientation space)
    theta = mod(theta, 180);
    
    % Mean orientation (midpoint of the cutoffs)
    mu = mean(orientation_cutoffs);
    
    % Calculate concentration parameter (kappa) based on orientation bandwidth
    % Higher kappa = narrower bandwidth
    bandwidth = abs(orientation_cutoffs(2) - orientation_cutoffs(1));
    kappa = 180 / bandwidth; % Simple heuristic, can be adjusted
    
    % Apply von Mises function for orientation
    % Converting back to radians for the von Mises formula
    theta_rad = theta * pi/180;
    mu_rad = mu * pi/180;
    orientation_bandpass_filter = exp(kappa * cos(2 * (theta_rad - mu_rad)));
    
    % Normalize the orientation filter
    orientation_bandpass_filter = normalize_array(orientation_bandpass_filter, 'min-max');

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