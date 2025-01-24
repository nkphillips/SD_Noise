
function [params, SSE, J] = fit_DoG_split(delta, y_observed, initial_params)

% Initial parameters: [amplitude_neg, amplitude_pos, w_neg, w_pos]

% Define options (optional)
options = optimoptions('lsqcurvefit', 'Display', 'off');

% Define parameter bounds
lower_bounds = [-6, -6, 1/12, 1/12];
upper_bounds = [6, 6, 1/0.1, 1/0.1];

% Fit the model
[params, SSE, ~, ~, ~, ~, J] = lsqcurvefit(@gaussian_prime_split, initial_params, delta, y_observed, lower_bounds, upper_bounds, options);

end