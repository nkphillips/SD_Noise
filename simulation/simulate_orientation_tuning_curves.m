% --- Script to Visualize Orientation Tuning and Noise Effects ---
% This script generates two figures:
% 1. A bank of idealized orientation tuning curves and their response to a stimulus.
% 2. A comparison of the population response under different conditions:
%    - Baseline (high-contrast, no-noise)
%    - High Orientation Noise (imprecise stimulus)
%    - Low Contrast (unreliable stimulus)
% -----------------------------------------------------------------

%% --- Setup and Parameters ---
clear; clc; close all;

line_width = 1;

% Define the orientation space
orientations = 0:1:180;

% Neuron/Channel Parameters
preferred_orientations = 0:30:150; % 6 channels from 0 to 150 deg
amplitude = 1.0;                  % Max response (arbitrary units)
tuning_width_sigma = 10;          % Tuning width (degrees)

% Stimulus Parameters
stimulus_orientation = 90;        % Stimulus to present (degrees)

%% --- Figure 1: Tuning Curve Bank and Stimulus Response ---

figure(1);
set(gcf, 'Color', 'w');

% Subplot 1: Plot the bank of tuning curves
subplot(1, 2, 1);
hold on;
colors = lines(length(preferred_orientations)); % Assign a color to each curve

for i = 1:length(preferred_orientations)
   pref_ori = preferred_orientations(i);
   % Gaussian tuning curve function (handling circular wraparound for aesthetics)
   delta_ori = min(abs(orientations - pref_ori), 180 - abs(orientations - pref_ori));
   tuning_curve = amplitude * exp(-delta_ori.^2 / (2 * tuning_width_sigma^2));
   plot(orientations, tuning_curve, 'LineWidth', line_width, 'Color', colors(i,:));
end

set(gca,'TickDir','out');
title('Bank of Orientation Tuning Curves');
xlabel('Stimulus Orientation (\circ)');
ylabel('Response (a.u.)');
legend(arrayfun(@(x) [num2str(x), '^\circ'], preferred_orientations, 'UniformOutput', false), 'Location', 'northeast');
xlim([0 180]);
xticks(0:30:180);
ylim([0 amplitude*1.1]);
box off;
axis square;

% Subplot 2: Plot the response of the bank to a single stimulus
subplot(1, 2, 2);
hold on;

% Calculate response of each channel to the stimulus
delta_stim = min(abs(stimulus_orientation - preferred_orientations), 180 - abs(stimulus_orientation - preferred_orientations));
responses = amplitude * exp(-delta_stim.^2 / (2 * tuning_width_sigma^2));

% Plot as a curve and with markers for each channel
plot(preferred_orientations, responses, 'k--', 'LineWidth', line_width);
scatter(preferred_orientations, responses, 100, colors, 'filled');

set(gca,'TickDir','out');
title(['Response to a ', num2str(stimulus_orientation), '^\circ Stimulus']);
xlabel('Channel Preferred Orientation (\circ)');
ylabel('Response (a.u.)');
xlim([-10 190]);
ylim([0 amplitude*1.1]);
xticks(0:30:180);
yticks(0:0.2:1);
box off;
axis square;

%% --- Figure 2: Visualizing Noise vs. Contrast Effects ---
figure(2);
set(gcf, 'Color', 'w'); 

% --- Subplot 1: Effect of Orientation Noise ---
subplot(1, 2, 1);
hold on;

% Parameters for noise simulation
orientation_noise_sigmas = [1, 20, 40]; % decreasing precision (half the original values)
green_colors = [152, 251, 152; 46, 139, 87; 0, 100, 0] / 255; % Shades of green (light to dark)

% Baseline response (no added noise)
delta_base = min(abs(orientations - stimulus_orientation), 180 - abs(orientations - stimulus_orientation));
response_baseline = amplitude * exp(-delta_base.^2 / (2 * tuning_width_sigma^2));
plot(orientations, response_baseline, 'k', 'LineWidth', 1, 'DisplayName', 'Baseline');

% Plot responses for increasing noise levels
for i = 1:length(orientation_noise_sigmas)
    noise_sigma = orientation_noise_sigmas(i);
    
    % The effect of stimulus noise is to broaden the population response.
    % This can be modeled as a convolution, where variances add.
    effective_sigma = sqrt(tuning_width_sigma^2 + noise_sigma^2);
    
    % Amplitude is kept fixed, while width increases.
    response_noise = amplitude * exp(-delta_base.^2 / (2 * effective_sigma^2));
    
    plot(orientations, response_noise, 'Color', green_colors(i,:), 'LineWidth', line_width, ...
        'DisplayName', sprintf('Noise \\sigma = %d^\\circ', noise_sigma));
end

% Formatting for subplot 1
set(gca,'TickDir','out');
title('Effect of Orientation Noise (Precision)');
xlabel('Channel Preferred Orientation (\circ)');
ylabel('Response (a.u.)');
legend('show', 'Location', 'northeast');
xlim([0 180]);
ylim([0 amplitude*1.1]);
xticks(0:30:180);
yticks(0:0.2:1);
box off;
ax = gca;
ax.FontSize = 12;
axis square;

% --- Subplot 2: Effect of Contrast ---
subplot(1, 2, 2);
hold on;

% Parameters for contrast simulation
contrast_levels = [0.9, 0.5, 0.2]; % 90%, 50%, 20%
blue_colors = [173, 216, 230; 100, 149, 237; 0, 0, 139] / 255; % Shades of blue (light to dark)

% Baseline response (100% contrast)
plot(orientations, response_baseline, 'k', 'LineWidth', 1, 'DisplayName', '100% Contrast');

% Plot responses for decreasing contrast levels
for i = 1:length(contrast_levels)
    contrast = contrast_levels(i);
    
    % Low contrast multiplicatively scales down the response amplitude
    response_contrast = (amplitude * contrast) * exp(-delta_base.^2 / (2 * tuning_width_sigma^2));
    
   plot(orientations, response_contrast, 'Color', blue_colors(i,:), 'LineWidth', line_width, ...
        'DisplayName', sprintf('%d%% Contrast', contrast*100));
end

% Formatting for subplot 2
set(gca,'TickDir','out');
title('Effect of Stimulus Contrast');
xlabel('Channel Preferred Orientation (\circ)');
ylabel('Response (a.u.)');
legend('show', 'Location', 'northeast');
xlim([0 180]);
ylim([0 amplitude*1.1]);
xticks(0:30:180);
yticks(0:0.2:1);
box off;
ax = gca;
ax.FontSize = 12;
axis square;