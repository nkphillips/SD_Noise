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

% current challenge is the initial best case for contrast and filter width, 
% ie what is the easiest contrast level and filter width level to achieve 90% performance
% lvl 1 stimulus will the be same for 'contrast' and 'precision'

%%% Performance goals (orientation discrimination):
% lvl 1 = 90%
% lvl 2 = 75%
% lvl 3 = 65%
% (or we do linspace(lvl_1_perf, lvl_3_perf, num_lvls))

% a. lvl 1 precision staircase first at 90% contrast
% b. lvl 1 contrast staircase with lvl 1 precision threshold
%   b1. ideally, performance is still 90%. if so, we lock in lvl 1 precision and contrast.
%   b2. if performance is not 90%, keep staircasing contrast?
% c. Interleave contrast and precision staircase blocks
%   c1. example: contrast lvl 2 staircase (filter width fixed to precision lvl 1 threshold); 
%       precision lvl 3 staircase (contrast fixed to lvl 1 threshold); 
%       precision lvl 2 staircase (contrast fixed to lvl 1 threshold); 
%       contrast lvl 2 staircase (filter width fixed to precision lvl 1 threshold);

%%% Staircase goal:
% 3 precision levels that achieve lvl 1, 2, 3 performance
% 3 contrast levels that achieve lvl 1, 2, 3 performance

% Each condition x level we use 2 staircase, one starting at the min stimulus val, the other at the max

%%% Staircase mechanics 
% contrast_init_step_size = 0.1;
% precision_init_step_size = 10;
% num_reversals = 16;
% num_trials = 1000;
% probe_offsets = round(linspace(0,15,7)); % ie  [0 3 5 8 10 13 15];
%
% Ignore first few reversals; after this ignoring window, start halving the step size when a reversal is detected
% When the max number of reversals are detected, stop the staircase. 


%% assignment: simulate contrast staircase

% 1. simulate responses (1 = correct, 0 = wrong)
% check datasample() weights argument
%       you can use datasample()    
% 2. simulate staircase (simulating stimulus values for every trial)
%   init_contrast = 0.9;
%   max_contrast = 0.9;
%   init_step_size_contrast = 0.1;
% Have staircase end after certain amount of reversals
%   
% 3. plot staircase: stimulus contrast as a function of trial

%% STEP 1: RESPONSE SIMULATION (contrast)

n_trials = 50;

probe_offsets = round(linspace(0,15,7));

offsets = datasample(probe_offsets, n_trials);
offsets = offsets .* datasample([-1 1], n_trials);

contrast = nan(n_trials,1);
resp     = nan(n_trials,1);   % 1 = correct, 0 = wrong
button_press = zeros(n_trials,1); % 1 = CW, 0 = CCW

contrast(1) = 0.9;
step = 0.1;
min_step = 0.01;
max_step = 0.1;
step_tracker = nan;
min_contrast = 0.05;
max_contrast = 0.9;

mu = 0;
sigma = 0.5;
num_correct = 0;

perc_correct = 0.75;
% run 3 staircases for 3 levels? Start with a staircase to get general
% contrast range. 
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
grid on
