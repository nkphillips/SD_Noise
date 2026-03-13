function best_guess = gridSearch(nll_func, lb, ub, n_grid)
% GRIDSEARCH Performs a simple 2D grid search to find a good starting point
% for subsequent optimization (e.g., using fmincon).
%
% Inputs:
%   nll_func : Function handle computing the negative log-likelihood (or cost)
%              given a 1x2 parameter vector.
%   lb       : 1x2 vector of lower bounds.
%   ub       : 1x2 vector of upper bounds.
%   n_grid   : (Optional) Number of grid points per dimension. Default is 20.
%
% Output:
%   best_guess : 1x2 vector of the parameters that yielded the lowest cost
%                on the grid.

if nargin < 4
    n_grid = 20;
end

% Create grid for parameter 1
p1_grid = linspace(lb(1), ub(1), n_grid);

% Create grid for parameter 2
p2_grid = linspace(lb(2), ub(2), n_grid);

min_nll = inf;
best_guess = [lb(1), lb(2)]; % Default to lower bounds

% Evaluate cost function over the grid
for p1 = p1_grid
    for p2 = p2_grid
        curr_nll = nll_func([p1, p2]);
        if curr_nll < min_nll
            min_nll = curr_nll;
            best_guess = [p1, p2];
        end
    end
end
end
