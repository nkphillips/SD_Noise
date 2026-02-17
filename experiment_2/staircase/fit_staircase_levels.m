% -------------------------------------------------------------------------
% PSYCHOMETRIC CURVE FITTING & INTERPOLATION SCRIPT
% -------------------------------------------------------------------------
% Purpose: Fit psychometric functions to pooled staircase data and
% interpolate stimulus levels for specific performance targets (65%, 80%, 90%).
%
% Models:
% 1. Contrast: Weibull Function
% 2. Filter Width: Equivalent Noise Model (Neuroscience-inspired)
% -------------------------------------------------------------------------

clear; clc; close all;

%% ------------------------------------------------------------------------
%  PART 1: SIMULATION / GROUND TRUTH DATA GENERATION
%  (Validating that we can recover known parameters from the functions)
% ------------------------------------------------------------------------

% --- Shared Constants ---
Signal_Tilt = 15; % Fixed signal strength (degrees) for the Width task
gamma = 0.5;      % Guess rate (2AFC)
lambda = 0.01;    % Lapse rate

% --- 1. Contrast Ground Truth (Weibull) ---
% Parameters: Alpha (Threshold), Beta (Slope)
True_Alpha_C = 0.12;
True_Beta_C  = 2.5;

% Generate Data
c_levels = [0.02, 0.05, 0.10, 0.20, 0.40, 0.80];
% Function: 0.5 + 0.49 * (1 - exp(-(x/alpha)^beta))
True_Weibull = @(b, x) gamma + (1 - gamma - lambda) .* (1 - exp(-(x ./ b(1)).^b(2)));
c_perf = True_Weibull([True_Alpha_C, True_Beta_C], c_levels);

% --- 2. Filter Width Ground Truth (Equivalent Noise) ---
% Parameters: Sigma_Int (Internal Noise)
True_SigmaInt = 8.5; % degrees

% Generate Data
w_levels = [2, 10, 20, 40, 60, 90];
% Function: normcdf( Signal / sqrt(2 * (sigma_ext^2 + sigma_int^2)) )
True_EqNoise = @(b, x) normcdf( Signal_Tilt ./ (sqrt(2) .* sqrt(x.^2 + b(1).^2)) );
w_perf = True_EqNoise(True_SigmaInt, w_levels);

% Targets we want to interpolate
targets = [0.65, 0.80, 0.90];

fprintf('--- SIMULATION SETTINGS ---\n');
fprintf('Contrast True Params: Alpha=%.2f, Beta=%.2f\n', True_Alpha_C, True_Beta_C);
fprintf('Width True Params:    Sigma_Int=%.2f\n', True_SigmaInt);
fprintf('---------------------------\n');

%% ------------------------------------------------------------------------
%  PART 2: CONTRAST ANALYSIS (Weibull Fit using fmincon)
% ------------------------------------------------------------------------
fprintf('\n--- Analyzing Contrast (Weibull) ---\n');

% 1. Define the Objective Function (Sum of Squared Errors)
% Same as True_Weibull defined above
Weibull_Func = True_Weibull;
Obj_Func_C = @(b) sum((c_perf - Weibull_Func(b, c_levels)).^2);

% 2. Fit the Data with fmincon
% Initial guess: Alpha=0.1, Beta=2.0
beta0_c = [0.1, 2.0];

% Constraints:
% Alpha (b1): Lower=0.001 (must be positive), Upper=1 (max contrast)
% Beta (b2):  Lower=0.1 (shallow slope), Upper=10 (steep slope)
lb_c = [0.001, 0.1];
ub_c = [1.0,   10.0];

options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
[beta_c, fval_c] = fmincon(Obj_Func_C, beta0_c, [], [], [], [], lb_c, ub_c, [], options);

fprintf('Fitted Alpha (Threshold): %.4f (True: %.4f)\n', beta_c(1), True_Alpha_C);
fprintf('Fitted Beta  (Slope):     %.4f (True: %.4f)\n', beta_c(2), True_Beta_C);

% 3. Inverse Calculation (Interpolation)
% Formula: x = alpha * [ -ln( (1-P-lambda)/(1-gamma-lambda) ) ]^(1/beta)
contrasts_needed = zeros(size(targets));

for i = 1:length(targets)
    P_t = targets(i);

    % Normalize target probability to 0-1 range (ignoring guess/lapse)
    K = (P_t - gamma) / (1 - gamma - lambda);

    % Solve for x
    term = -log(1 - K);
    contrasts_needed(i) = beta_c(1) * (term)^(1/beta_c(2));

    fprintf('Target: %.0f%% Correct -> Required Contrast: %.3f\n', ...
        P_t*100, contrasts_needed(i));
end

%% ------------------------------------------------------------------------
%  PART 3: FILTER WIDTH ANALYSIS (Equivalent Noise Model using fmincon)
% ------------------------------------------------------------------------
fprintf('\n--- Analyzing Filter Width (Eq. Noise Model) ---\n');

% 1. Define the Objective Function (Sum of Squared Errors)
% Same as True_EqNoise defined above
EqNoise_Func = True_EqNoise;
Obj_Func_W = @(b) sum((w_perf - EqNoise_Func(b, w_levels)).^2);

% 2. Fit the Data with fmincon
% Initial guess: Internal noise = 5 degrees
beta0_w = 5;

% Constraints:
% Sigma_Int (b1): Lower=0 (cannot be negative), Upper=100
lb_w = 0;
ub_w = 100;

[beta_w, fval_w] = fmincon(Obj_Func_W, beta0_w, [], [], [], [], lb_w, ub_w, [], options);
sigma_int = beta_w(1);

fprintf('Fitted Internal Noise (sigma_int): %.2f degrees (True: %.2f)\n', sigma_int, True_SigmaInt);

% 3. Inverse Calculation (Interpolation)
% Formula: x = sqrt( (Signal/d_prime)^2 - sigma_int^2 )
widths_needed = zeros(size(targets));

for i = 1:length(targets)
    P_t = targets(i);

    % A. Convert Target % to d' (sensitivity)
    d_prime_target = sqrt(2) * norminv(P_t);

    % B. Solve for External Noise (x)
    term_sq = (Signal_Tilt / d_prime_target)^2 - sigma_int^2;

    if term_sq < 0
        warning('Target %.2f%% is impossible (Internal noise too high)', P_t);
        widths_needed(i) = NaN;
    else
        widths_needed(i) = sqrt(term_sq);
    end

    fprintf('Target: %.0f%% Correct -> Required Width:    %.1f deg\n', ...
        P_t*100, widths_needed(i));
end

%% ------------------------------------------------------------------------
%  PART 4: VISUALIZATION
% ------------------------------------------------------------------------
figure('Color', 'w', 'Position', [100 100 1000 400]);

% --- Plot Contrast ---
subplot(1, 2, 1);
scatter(c_levels, c_perf, 60, 'k', 'filled'); hold on;
fplot(@(x) Weibull_Func(beta_c, x), [0 1], 'r-', 'LineWidth', 2);
plot(contrasts_needed, targets, 'b*', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Contrast (0-1)'); ylabel('% Correct');
title(sprintf('Contrast Fit (Weibull)\nTrue \\alpha=%.2f, \\beta=%.2f', True_Alpha_C, True_Beta_C));
legend('Sim Data', 'Fit', 'Interpolated Targets', 'Location', 'southeast');
axis square;

% --- Plot Filter Width ---
subplot(1, 2, 2);
scatter(w_levels, w_perf, 60, 'k', 'filled'); hold on;
fplot(@(x) EqNoise_Func(beta_w, x), [0 120], 'r-', 'LineWidth', 2);
plot(widths_needed, targets, 'b*', 'MarkerSize', 10, 'LineWidth', 1.5);
xlabel('Filter Width (deg)'); ylabel('% Correct');
title(sprintf('Width Fit (Eq. Noise)\nTrue \\sigma_{int} = %.1f^\\circ', True_SigmaInt));
legend('Sim Data', 'Fit', 'Interpolated Targets', 'Location', 'northeast');
axis square;
