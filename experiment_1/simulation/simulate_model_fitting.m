
clear all; clc;
% close all; 

functions_dir = '../functions';
addpath(functions_dir);

%% Simulate data

num_trials = 1000;

% Random list of probe offsets
offsets = -15:0.01:15;

% Random list of test orientations from 0 - 179 deg
sim_test_orientation = rand(1, num_trials) * 179;

rand_offsets = datasample(offsets, num_trials, 'Replace', true);
simulation_probe_orientations = sim_test_orientation + rand_offsets;

% Calculate probability of correct discrimination using cumulative normal
subj_mu = 0;
subj_sigma = 3;
guess_rate = 0.25;

p_CW = calc_pCW(rand_offsets, subj_mu, subj_sigma, guess_rate);

%{
figure('Color', 'w');
subplot(1,4,1)
plot(sort(normcdf(rand_offsets, subj_mu, subj_sigma)))
ylim([0, 1])
title('CDF')
box off

subplot(1,4,2)
plot(sort((1 - guess_rate) * normcdf(rand_offsets, subj_mu, subj_sigma)))
title('Amplitude * CDF')
ylim([0, 1])
box off

subplot(1,4,3)
plot(sort(normcdf(rand_offsets, subj_mu, subj_sigma) + (0.5 * guess_rate)))
title('CDF + Baseline')
% ylim([0, 1])
box off

subplot(1,4,4)
plot(sort(p_CW))
title('Model')
ylim([0, 1])
box off
%}

% Simulate response based on discrimination probability
responses = rand() < p_CW;

%% Grid search for mu and sigma

unique_probe_offsets = sort([-linspace(3, 15, 7), linspace(3, 15, 7)]);
p_CW = [ 0         0         0    0.5000    0.4000    0.5455    0.2500    0.8333    0.7143    0.6667    1.0000    1.0000    0.5714    1.0000];

% [best_params, best_sse] = grid_search(rand_offsets, p_CW, guess_rate);
[best_params, best_sse] = grid_search(unique_probe_offsets, p_CW, guess_rate);


%% Estimate mu and sigma

options = optimoptions('fmincon');

init_params = best_params;
% init_params = [0,1];

lower_bounds = [-15, 1]; % lower bounds for [mu, sigma]
upper_bounds = [15, 20];

free_params = init_params;
% fixed_params{1} = [rand_offsets', p_CW'];
fixed_params{1} = [unique_probe_offsets', p_CW'];
fixed_params{2} = guess_rate;

response_model = @(free_params) calc_fit(free_params, fixed_params);

[params_est, sse, exit_flag] = fmincon(response_model, init_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

%% Generate estimated curves

mu = params_est(1);
sigma = params_est(2);

x = -20:0.01:20;
cdf_response = (1 - guess_rate) * normcdf(x, mu, sigma) + 0.5 * guess_rate;

x_exp = (1 - guess_rate) * normcdf(unique_probe_offsets, mu, sigma) + 0.5 * guess_rate;
r2 = 1 - (sum((p_CW - x_exp).^2) / sum((p_CW - mean(p_CW)).^2));

pdf_response = (1 - guess_rate) * normpdf(x, mu, sigma);
pdf_response = pdf_response / max(pdf_response); % normalize to peak at 1

sorted_p_CW = sortrows([unique_probe_offsets', p_CW']);

figure('Color', 'w','Name', '3 -> 1 contrast');
plot(x, cdf_response', 'LineWidth', 2,'Color', 'b');
hold on;
plot(x, pdf_response', 'LineWidth', 2,'LineStyle', '--','Color', 'b');
scatter(unique_probe_offsets, p_CW, 30, 'filled');
line([0,0], [0,1], 'LineWidth', 2, 'Color', 'k');
line([-20,20], [0.5,0.5], 'LineWidth', 2, 'Color', 'k');
title(['\mu = ' num2str(round(mu, 2)) ', \sigma = ' num2str(round(sigma, 2)) ' (R^2 = ' num2str(round(r2, 2)) ')']);
xlabel('Probe Offset');
ylabel('P(Resp|CW)');
ylim([0, 1]);
xlim([-20, 20]);
box off
set(gca, 'TickDir', 'out');

% figure('Color', 'w');
% plot(sorted_p_CW(:,1), sorted_p_CW(:,2), 'k.');
% ylim([0, 1]);
% xlim([-20, 20]);
% box off
% set(gca, 'TickDir', 'out');
