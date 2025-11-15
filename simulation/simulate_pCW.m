%%% simulate_pCW

%% Prepare workspace

clear all; close all; clc;

script_dir = pwd;
functions_dir =  '../analyses/functions'; addpath(functions_dir);

%%

guess_rate = 0;

max_offset = 20;
min_offset = -max_offset;
probe_offset = min_offset:0.1:max_offset;

%% No encoding bias

mu = 0;
sigma = 5;

p_cw = calc_pCW(probe_offset, mu, sigma, guess_rate);
norm_pdf = normpdf(probe_offset, mu, sigma);
norm_pdf = norm_pdf / max(norm_pdf);

figure('Color', 'white', 'Name', 'No encoding bias')
plot(probe_offset, p_cw, 'Color', 'black')
hold on
plot(probe_offset, norm_pdf, 'LineStyle', '--', 'Color', 'black')
xline(0)

ylim([0,1])
ylabel('P(RespCW)')
xlabel('\delta\theta')
box off;
set(gca, 'TickDir', 'out')

%% Counter clockwise encoding bias

mu = -5;
sigma = 5;

p_cw = calc_pCW(probe_offset, mu, sigma, guess_rate);
norm_pdf = normpdf(probe_offset, mu, sigma);
norm_pdf = norm_pdf / max(norm_pdf);

figure('Color', 'white', 'Name', 'Counter clockwise encoding bias')
plot(probe_offset, p_cw, 'Color', 'black')
hold on
plot(probe_offset, norm_pdf, 'LineStyle', '--', 'Color', 'black')
xline(0)

ylim([0,1])
ylabel('P(RespCW)')
xlabel('\delta\theta')
box off;
set(gca, 'TickDir', 'out')

%% Clockwise encoding respose

mu = 5;
sigma = 5;

