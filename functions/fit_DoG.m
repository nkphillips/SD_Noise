% Estimate Derivative of Gaussian

function [params, SSE, J] = fit_DoG(delta, y_observed, initial_params)

% Initial parameters: [amplitude, w]

% Define options (optional)
options = optimoptions('lsqcurvefit', 'Display', 'off');

% Fit the model
[params, SSE, ~, ~, ~, ~, J] = lsqcurvefit(@gaussian_prime, initial_params, delta, y_observed, [], [], options);


end
