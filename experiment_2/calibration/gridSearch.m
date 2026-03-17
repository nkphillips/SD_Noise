function best_guess = gridSearch(nll_func, lb, ub, n_grid)
% GRIDSEARCH N-dimensional grid search for optimization starting points.
%
% Inputs:
%   nll_func : Function handle computing the cost given a 1xN parameter vector.
%   lb       : 1xN vector of lower bounds.
%   ub       : 1xN vector of upper bounds.
%   n_grid   : (Optional) Points per dimension. Default scales with N to
%              keep total evaluations manageable.

n_dims = length(lb);

if nargin < 4
    % ~8000 total evaluations max; scale per-dim grid accordingly
    n_grid = max(5, floor(8000^(1/n_dims)));
end

grids = cell(1, n_dims);
for d = 1:n_dims
    grids{d} = linspace(lb(d), ub(d), n_grid);
end

% Build all combinations via ndgrid
grid_arrays = cell(1, n_dims);
[grid_arrays{:}] = ndgrid(grids{:});

n_total = numel(grid_arrays{1});
min_nll = inf;
best_guess = lb;

for idx = 1:n_total
    params = zeros(1, n_dims);
    for d = 1:n_dims
        params(d) = grid_arrays{d}(idx);
    end
    curr_nll = nll_func(params);
    if isfinite(curr_nll) && curr_nll < min_nll
        min_nll = curr_nll;
        best_guess = params;
    end
end
end
