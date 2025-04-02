
clear all; close all; clc;

functions_dir = '../functions'; addpath(functions_dir);

%% Simulate data

num_trials = 1000;
offsets = -15:15;

rand_offsets = datasample(offsets, num_trials, 'Replace', true);

subj_mu = -3;
subj_sigma = 7;

p_CW = normcdf(sort(rand_offsets), subj_mu, subj_sigma);

%% Estimate mu and sigma

options = optimset('MaxFunEvals', 100000, 'MaxIter', 10000, 'display','off');

init_params = [0, 1];

lower_bounds = [-15, 1]; % lower bounds for [mu, sigma]
upper_bounds = [15, 7];

free_params = init_params;
fixed_params = [sort(rand_offsets)', p_CW'];

response_model = @(free_params) calc_fit(free_params, fixed_params);

[params_est, sse, exit_flag] = fmincon(response_model, init_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

%% Generate estimated curves

mu = params_est(1);
sigma = params_est(2);

x = -20:0.01:20;
cdf_response = normcdf(x, mu, sigma);

pdf_response = normpdf(x, mu, sigma);
pdf_response = pdf_response / max(pdf_response); % normalize to peak at 1

sorted_p_CW = sortrows([rand_offsets', p_CW']);

figure('Color', 'w');
plot(x, cdf_response', 'LineWidth', 2,'Color', 'b');
hold on;
plot(x, pdf_response', 'LineWidth', 2,'LineStyle', '--','Color', 'b');
line([0,0], [0,1], 'LineWidth', 2, 'Color', 'k');
title(['\mu = ' num2str(round(mu, 2)) ', \sigma = ' num2str(round(sigma, 2))]);
xlabel('Probe Offset');
ylabel('P(Resp|CW)');
ylim([0, 1]);
xlim([-20, 20]);
box off
set(gca, 'TickDir', 'out');

figure('Color', 'w');
scatter(sorted_p_CW(:,1), sorted_p_CW(:,2));
ylim([0, 1]);
xlim([-20, 20]);
box off
set(gca, 'TickDir', 'out');
