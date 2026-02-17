% run_staircase
% experiment 2

%% Prepare workspace



%% Toggles



%% Set directories

dirs.script_dir = pwd;
dirs.functions_dir = '../analyses/functions'; addpath(dirs.functions_dir);

%% Setup devices and display, open window



%% Loading screen


%% Define stimulus
% For precision (ie filter width), we could use a template noise sample.
%       staircase
% For contrast, we can just use template filtered noise texture ().

% max_contrast = 0.9; % Michelson contrast
% min_contrast = 0.05; % Michelson contrast
% max_filter_width = 85; % °
% min_filter_width = 2; % °

%%% Performance goals (orientation discrimination):
% lvl 1 = 90%
% lvl 2 = 75%
% lvl 3 = 65%
% (or we do linspace(lvl_1_perf, lvl_3_perf, num_lvls))

%%% Staircase goal:
% Each condition x level we use 2 staircases, a. 2-down 1-up, b. 3-down 1-up.
% Since we're assuming level 1 contrast and precision, the objective is to get contrast and precision levels for levels 2 and 3.

%%% Staircase mechanics
% contrast_init_step_size = 0.1;
% precision_init_step_size = 10;
% num_trials = 60;
% probe_offsets = round(linspace(0,15,7)); % ie  [0 3 5 8 10 13 15];
%
% Ignore first few reversals; after this ignoring window, start halving the step size when a reversal is detected
% When the max number of reversals are detected, stop the staircase.


%% assignment: simulate contrast staircase

% 1. simulate responses (1 = correct, 0 = wrong)
% check datasample() weights argument
%       you can use datasample()
% 2. simulate staircase (simulating stimulus values for every trial)
% staircase a: 2-down 1-up
% staircase b: 3-down 1-up
% the staircases are randomly interleaved.
%   init_contrast = 0.9;
%   max_contrast = 0.9;
%   init_step_size_contrast = 0.1;
% Have staircase end after certain amount of trials (eg., 60 trials, but usually are good with ~40 trials)
%
% 3. plot staircase: stimulus contrast as a function of trial

%% STEP 1: Simulate psychometric function (contrast)

contrast_vals = 0:0.01:1;
sim_psych_func_params = [];

free_params = sim_psych_func_params;
fixed_params = contrast_vals;
sim_psych_func = calcPsychometricFunction(free_params, fixed_params);

%% STEP 2: RESPONSE SIMULATION (contrast)

n_trials = 60;

probe_offsets = round(linspace(0,15,7));

offsets = datasample(probe_offsets, n_trials);
offsets = offsets .* datasample([-1 1], n_trials);
correct_ans = offsets > 0;

contrast = nan(n_trials,1);
resp     = nan(n_trials,1);   % 1 = correct, 0 = wrong
button_press = zeros(n_trials,1); % 1 = CW, 0 = CCW

min_contrast = 0.05;
max_contrast = 0.9;

contrast(1) = max_contrast;
step = 0.1;
min_step = 0.01;
max_step = 0.1;
step_tracker = nan;

% define observer bias
mu = 0;
sigma = 0.5;

num_correct = 0;


for t = 1:n_trials-1

    probe_offset = offsets(t);
    correct_resp = probe_offset > 0; % 1 = positive (CW), 0 = negative (CCW)

    p_CW = calc_pCW(probe_offset, mu, sigma, 0.1);
    button_press(t) = rand() < p_CW; % 1 = CW, 0 = CCW

    resp(t) = button_press(t) == correct_resp;

    if resp(t) == 1
        num_correct = num_correct + 1;
    elseif resp(t) == 0
        num_correct = 0;
        if ~isnan(step_tracker) && step_tracker == -1
            step = step/2;
            step = min(max(step,min_step),max_step);
        end
        step_tracker = 1;
        contrast(t+1) = contrast(t) + step;  % easier
    end

    % staircase rule (2-up-1-down)
    if num_correct == 2
        num_correct = 0;
        if ~isnan(step_tracker) && step_tracker == 1
            step = step/2;
            step = min(max(step,min_step),max_step);
        end
        step_tracker = -1;
        contrast(t+1) = contrast(t) - step;  % harder
    elseif num_correct == 1
        contrast(t+1) = contrast(t);
    end

    % staircase rule (3-up 1-down)

    % clamp
    contrast(t+1) = min(max(contrast(t+1),min_contrast),max_contrast);
end

figure;
yyaxis left
plot(1:n_trials,contrast,'LineWidth',1.5,'Marker','.')
ylabel('Contrast')
ylim([0 1])
yyaxis right
scatter(1:n_trials,resp,100,'LineWidth',1.5, 'Marker','.')
ylabel('Response (1=correct, 0=incorrect)')
xlabel('Trial')
title('Simulated contrast staircase')
box off;
set(gca, 'TickDir', 'out');

%% Estimate psychometric function

% Find all the unique contrast values, get the corresponding responses, then calculate performance for each value.
% This is what's fed into the psychometric function estimation.

unique_contrast_vals = unique(contrast);
performance = nan(length(unique_contrast_vals),1);
for i = 1:length(unique_contrast_vals)
    performance(i) = mean(resp(contrast == unique_contrast_vals(i)));
end

% Weibull psychometric function?
psych_func = @(free_params) calcPsychometricFunction(free_params, fixed_params);

init_params = [0.5, 1, 0.5, 0]; % this is what you're estimating (max performance, slope, C_50, baseline)
lower_bounds = [0, 0, 0, 0]; % (max performance, slope, C_50, baseline)
upper_bounds = [1, 10, 10, 1]; % (max performance, slope, C_50, baseline)

options = optimoptions('fmincon');
[params_est, sse, exit_flag] = fmincon(psych_func, init_params, [], [], [], [], lower_bounds, upper_bounds, [], options);

est_psych_func = psych_func(contrast_vals, params_est);


perc_correct = [0.65 0.75 0.90];

% get the x-values that correspond to the performance goals
x_vals = interp1(contrast_vals, est_psych_func, perc_correct); % look into this...

% Plot psychometric function and data
% Should be s-shaped curve that fits the data.
figure('Color', 'w');
plot(contrast_vals, est_psych_func, 'LineWidth', 1);
hold on;
scatter(contrast, resp, 100, 'LineWidth', 1, 'Marker','.');
xlabel('Contrast');
ylabel('Response probability');
title('Psychometric function');
box off;
set(gca, 'TickDir', 'out');
