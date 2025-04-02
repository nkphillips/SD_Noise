%% Simulate data
num_trials = 1000;
offsets = -15:15;
sim_test_orientation = rand(1, num_trials) * 179;
rand_offsets = datasample(offsets, num_trials, 'Replace', true);
simulation_probe_orientations = sim_test_orientation + rand_offsets;

% Calculate probability of correct discrimination using cumulative normal
subj_mu = -3;
subj_sigma = 7;
threshold = 10;
slope = 5;
x = rand_offsets;
p_CW = 0.5 + 0.25 * normcdf(x, subj_mu, subj_sigma);
p_CW = p_CW / max(p_CW);

% Simulate response based on discrimination probability
responses = rand() < p_CW;

%% Set up optimization with iteration tracking
options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000, 'display','off');
init_params = [0, 1];
lower_bounds = [-15, 1];
upper_bounds = [15, 7];

% Create a structure to store iteration information
global iteration_data
iteration_data = struct('params', [], 'fval', []);

% Create a wrapper function that captures iteration data
function [f, g, H] = iteration_wrapper(params, fixed_params)
    global iteration_data
    iteration_data.params = [iteration_data.params; params];
    [f, g, H] = calc_fit(params, fixed_params);
    iteration_data.fval = [iteration_data.fval; f];
end

% Set up the optimization problem
fixed_params = [x', p_CW'];
response_model = @(free_params) iteration_wrapper(free_params, fixed_params);

% Run the optimization
[params_est, sse, exit_flag] = fmincon(response_model, init_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

%% Create animation
figure('Color', 'w', 'Position', [100 100 800 600]);
x = -20:0.01:20;

% Create video writer object
v = VideoWriter('optimization_animation.mp4', 'MPEG-4');
v.FrameRate = 5;  % 5 frames per second
open(v);

% Plot the actual data points
sorted_p_CW = sortrows([rand_offsets', p_CW']);
scatter(sorted_p_CW(:,1), sorted_p_CW(:,2), 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.3);
hold on;

% Plot the initial curve
init_mu = init_params(1);
init_sigma = init_params(2);
init_cdf = normcdf(x, init_mu, init_sigma);
init_pdf = normpdf(x, init_mu, init_sigma);
init_pdf = init_pdf / max(init_pdf);
h_cdf = plot(x, init_cdf, 'LineWidth', 2, 'Color', 'b');
h_pdf = plot(x, init_pdf, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
line([0,0], [0,1], 'LineWidth', 2, 'Color', 'k');

% Set up the plot
ylim([0, 1]);
xlim([-20, 20]);
box off
set(gca, 'TickDir', 'out');
xlabel('Probe Offset');
ylabel('P(Resp|CW)');

% Create title with parameter values
title_text = title(['\mu = ' num2str(round(init_mu, 2)) ', \sigma = ' num2str(round(init_sigma, 2))]);

% Animate through iterations
for i = 1:size(iteration_data.params, 1)
    % Update curves
    mu = iteration_data.params(i,1);
    sigma = iteration_data.params(i,2);
    
    cdf_response = normcdf(x, mu, sigma);
    pdf_response = normpdf(x, mu, sigma);
    pdf_response = pdf_response / max(pdf_response);
    
    set(h_cdf, 'YData', cdf_response);
    set(h_pdf, 'YData', pdf_response);
    set(title_text, 'String', ['\mu = ' num2str(round(mu, 2)) ', \sigma = ' num2str(round(sigma, 2))]);
    
    % Add frame to video
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    % Add a small pause to make the animation smoother
    pause(0.1);
end

% Close the video writer
close(v);

% Add final plot with optimal parameters
figure('Color', 'w');
scatter(sorted_p_CW(:,1), sorted_p_CW(:,2), 'k', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.3);
hold on;
plot(x, cdf_response, 'LineWidth', 2, 'Color', 'b');
plot(x, pdf_response, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'b');
line([0,0], [0,1], 'LineWidth', 2, 'Color', 'k');
title(['Final fit: \mu = ' num2str(round(params_est(1), 2)) ', \sigma = ' num2str(round(params_est(2), 2))]);
xlabel('Probe Offset');
ylabel('P(Resp|CW)');
ylim([0, 1]);
xlim([-20, 20]);
box off
set(gca, 'TickDir', 'out'); 