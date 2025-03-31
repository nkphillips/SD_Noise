%%% simulate_serial_dependence


% Define possible orientations
orientation_range = [0 179];
probe_offset_range = [1 10];

%% Generate trials

num_trials = 1000;

stim_orientations = orientation_range(1) + (orientation_range(2) - orientation_range(1)) .* rand(num_trials, 1); 
probe_offsets = probe_offset_range(1) + (probe_offset_range(2) - probe_offset_range(1)) .* rand(num_trials, 1);
probe_orientations = calc_probe_orientation(stim_orientations, probe_offsets);

delta_theta = probe_orientations - stim_orientations;
Delta_theta = stim_orientations(1:end-1) - stim_orientations(2:end);

%% Simulate serial dependence

num_curves = 5;

amp = linspace(1, 5, num_curves);
w = linspace(1, 5, num_curves);

y = nan(num_curves, length(Delta_theta));

for n = 1:num_curves
    y(n,:) = gaussian_prime([amp(n), w(n)], Delta_theta);
end

%% Plot

figure('Color', 'w', 'Name', 'Serial dependence');

plot(Delta_theta, y);






