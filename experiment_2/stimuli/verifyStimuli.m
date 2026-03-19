% verifyStimuli
%   Verifies the generated textures for Experiment 2 by computing their
%   orientation power spectrum. Plots power spectra and saves the figures.
%   If textures don't exist, it can trigger their generation.
%
%   This ensures that the orientation bandpass filtering produced the
%   expected energy distribution.
%
% Syntax
%   verifyStimuli()
%

function verifyStimuli()

%% Toggles

toggles.save_texture_figures = true;
toggles.save_energy_figures = true;

%% Set directories

script_dir = fileparts(mfilename('fullpath'));
experiment_dir = fileparts(script_dir);
project_dir = fileparts(experiment_dir);

task_dir = fullfile(experiment_dir, 'experiment');
texture_dir = fullfile(task_dir, 'textures');
functions_dir = fullfile(task_dir, 'functions');

texture_figures_dir = fullfile(script_dir, 'figures', 'textures');
energy_figures_dir = fullfile(script_dir, 'figures', 'energy');

% Add directories to path
addpath(functions_dir);

% Create texture figures directory
if ~exist(texture_figures_dir, 'dir')
    mkdir(texture_figures_dir);
end

% Create energy figures directory
if ~exist(energy_figures_dir, 'dir')
    mkdir(energy_figures_dir);
end

%% Screen parameters (Match your default setup)

w.view_distance = 57; % cm
w.screen_width = 30; % cm
w.screen_width_px = 1512; % Typical Macbook Pro 14
w.screen_height_px = 982;
w.centerX = round(w.screen_width_px/2);
w.centerY = round(w.screen_height_px/2);

w.visual_angle = 2 * atan2d(w.screen_width/2,  w.view_distance);
w.ppd = w.screen_width_px/w.visual_angle;
w.px_size = w.screen_width/w.screen_width_px;
w.f_Nyquist = 1/(2*w.px_size);
w.gray = 127;

%% Stimuli parameters (Matches init_stimuli_params.m for Calibration)
p.calibration = 1;
p.training = 0;
p.subj_ID = '999';
p.display_setup = 'Macbook';

p.aperture_height_px = w.screen_height_px;
p.aperture_width_px = w.screen_width_px;
p.aperture_radius_px = round(w.ppd * 3);

p.height_px = round(p.aperture_radius_px*2 + w.ppd);
p.width_px = round(p.aperture_radius_px*2 + w.ppd);

p.num_levels = 7;
p.contrast_min = 0.05;
p.contrast_max = 0.9;
p.contrast = round(logspace(log10(p.contrast_min), log10(p.contrast_max), p.num_levels),2);

p.filter_width_min = 2;
p.filter_width_max = 80;
p.orientation_bp_filter_width = round(logspace(log10(p.filter_width_min), log10(p.filter_width_max), p.num_levels),2);

p.sf_bp_filter_cutoffs = [1 4]; % unit: cycles/degree

p.num_noise_samples = 5; % Only need a few to verify, 20 is heavy

%% Load textures if they exist, otherwise generate them

texture_filename = ['SD_Noise_Exp2_Calibration_textures_' p.display_setup '.mat'];
texture_filepath = fullfile(texture_dir, texture_filename);

fprintf('Target texture file to load: %s\n', texture_filename);

if ~exist(texture_filepath, 'file')
    fprintf('Texture file not found. Generating stimuli on the fly for verification...\n');

    n_filter_widths = length(p.orientation_bp_filter_width);
    textures = nan(p.height_px, p.width_px, p.num_levels, n_filter_widths, p.num_noise_samples);

    target_contrast_idx = p.num_levels; % Max contrast
    target_contrast_val = p.contrast(target_contrast_idx);

    % Only generate what we need to verify: 1 sample at max contrast for all filter widths
    for fw_idx = 1:n_filter_widths
        base_noise = create_noise_texture(p.height_px, p.width_px);
        noise_texture = bandpassFilterImg(base_noise, [round(0 - p.orientation_bp_filter_width(fw_idx)/2), floor(0 + p.orientation_bp_filter_width(fw_idx)/2)], p.sf_bp_filter_cutoffs, w.ppd * 0.1, w.f_Nyquist);
        noise_texture = centerTextureContrast(noise_texture, target_contrast_val, w.gray);
        textures(:,:,target_contrast_idx, fw_idx, 1) = noise_texture;
    end
else
    load(texture_filepath, 'stimuli');
    disp('Textures loaded.')
    textures = stimuli.test_textures;
end

%% Prepare figure handle

if toggles.save_texture_figures || toggles.save_energy_figures
    fg = figure('Visible','on','Color',[1 1 1], 'Position', [100 100 1200 400]);
    set(0,'CurrentFigure',fg);
end

%% Run texture energy verification (Focusing on Orientation, not SF)

fprintf('Starting texture verification...\n');

n_filter_widths = length(p.orientation_bp_filter_width);
target_contrast_idx = p.num_levels; % Max contrast
sample_idx = 1; % Pick the first exemplar

for fw_idx = 1:n_filter_widths

    fprintf('Processing filter width %d of %d (%.2f deg)... \n', fw_idx, n_filter_widths, p.orientation_bp_filter_width(fw_idx));

    %% Pull texture

    texture = textures(:,:,target_contrast_idx, fw_idx, sample_idx);

    % Save texture figure
    if toggles.save_texture_figures
        subplot(1,3,1);
        imagesc(texture), colormap(gray), axis square, box off, axis off;
        title(sprintf('Filter Width: %.2f\\circ', p.orientation_bp_filter_width(fw_idx)));
    end

    %% Compute 2D Power Spectrum

    % Subtract mean to remove DC component
    texture_zero_mean = texture - mean(texture(:));

    % Apply a 2D Hann window to smoothly taper the edges to zero.
    % This completely eliminates the artificial energy spikes at 0 and 90 degrees
    % caused by the sharp square boundaries of the image matrix.
    [tr, tc] = size(texture_zero_mean);
    win_r = 0.5 * (1 - cos(2*pi*(0:tr-1)'/(tr-1)));
    win_c = 0.5 * (1 - cos(2*pi*(0:tc-1)/(tc-1)));
    window_2d = win_r * win_c;
    
    texture_windowed = texture_zero_mean .* window_2d;

    % Pad to high power of 2 for better angular resolution
    padding = 2^nextpow2(max(size(texture_windowed)) * 2);

    fft_texture = fftshift(fft2(texture_windowed, padding, padding));
    power_spectrum = abs(fft_texture).^2;

    if toggles.save_energy_figures
        subplot(1,3,2);
        % Plotting amplitude/power linearly rather than log10 avoids exposing 
        % the low-energy sinc-ringing (grids) from the square edges, 
        % resulting in a clean white-on-black representation.
        imagesc(abs(fft_texture));
        colormap(gray);
        axis square; box off; axis off;
        title('2D Frequency Spectrum');
    end

    %% Create coordinate space for angular integration

    [rows, cols] = size(fft_texture);

    u = (-cols/2:cols/2-1) / cols;
    v = (-rows/2:rows/2-1) / rows;
    [U, V] = meshgrid(u, v);

    % Compute polar angle (in degrees, -180 to 180)
    theta = atan2d(V, U);

    % Wrap negative angles to 0-180 since power spectrum is symmetric
    theta(theta < 0) = theta(theta < 0) + 180;

    %% Angular average power spectrum

    % We want to see the energy profile as a function of orientation (0-180)
    step_sz = 1; % 1 degree steps
    bin_edges = 0:step_sz:180;

    [~, ~, bin_idx] = histcounts(theta(:), bin_edges);

    valid_mask = bin_idx > 0;
    bin_idx = bin_idx(valid_mask);
    power_vals = power_spectrum(valid_mask);

    num_bins = length(bin_edges) - 1;
    angular_power = accumarray(bin_idx, power_vals, [num_bins 1]);
    angular_counts = accumarray(bin_idx, 1, [num_bins 1]);

    % Compute mean angular energy profile
    energy_profile = angular_power ./ angular_counts;

    % Normalize
    energy_profile = energy_profile / max(energy_profile);

    % Generate bin centers
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    %% Plot power spectrum

    if toggles.save_energy_figures

        subplot(1,3,3);
        plot(bin_centers, energy_profile, 'k-', 'LineWidth', 1.5);
        hold on;

        % The target center is 0/180 in Fourier space (vertical spatial orientation).
        % Energy wraps around both edges of the 0-180 axis.
        target_width = p.orientation_bp_filter_width(fw_idx);

        % Plot expected bandwidth bounds at both edges
        xline(target_width/2, 'r--', 'LineWidth', 1.5);
        xline(180 - target_width/2, 'r--', 'LineWidth', 1.5);

        % Format figure
        title(sprintf('Orientation Energy (Target Width: %.2f\\circ)', target_width));
        xlabel('Orientation (\circ)');
        ylabel('Normalized power');
        xlim([0 180]);
        xticks(0:45:180);
        ylim([0 1.05]);
        box off; axis square;
        legend({'Energy', 'Filter Bounds'}, 'Location', 'northeast');

        % Save figure
        figure_name = sprintf('Filter_%02d_Orientation_Spectrum', fw_idx);
        figure_path = fullfile(energy_figures_dir, [figure_name '.pdf']);

        % Use exportgraphics to prevent clipping and ensure tight bounding box
        if exist('exportgraphics', 'file')
            exportgraphics(gcf, figure_path, 'Resolution', 300);
        else
            set(gcf, 'PaperPositionMode', 'auto');
            set(gcf, 'PaperOrientation', 'landscape');
            print(gcf, figure_path, '-dpdf', '-bestfit');
        end

        disp(['Saved ' figure_name '.pdf in ' energy_figures_dir]);
        clf;

    end

end

fprintf('Verification complete.\n');

if toggles.save_texture_figures || toggles.save_energy_figures
    close(fg)
end

end