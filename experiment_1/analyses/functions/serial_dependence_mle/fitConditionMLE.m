function [params, exit_flag] = fitConditionMLE(delta_theta, x_probe, response, opts)
% fitConditionMLE  Trial-level MLE for S&S DoG + baseline + Gaussian psychometric with 25% guess rate.
%
%   mu_t = dogIsolated(dtheta, A, w) + beta                                (S&S Eq. 1)
%   P_t  = (1 - lambda) * normcdf(x_probe, mu_t, sigma) + 0.5 * lambda     (lambda = 0.25)
%
% params = [A, w, sigma, beta] (1x4); w in 1/deg. NaN(1,4) if no trials.
% exit_flag is the fmincon exit flag (NaN on missing data, < 1 means failed/at-bound stop).

    if nargin < 4 || isempty(opts)
        opts = struct();
    end

    delta_theta = delta_theta(:);
    x_probe = x_probe(:);
    response = response(:);

    n = min([numel(delta_theta), numel(x_probe), numel(response)]);
    if n < 1
        params = nan(1, 4);
        exit_flag = NaN;
        return
    end

    delta_theta = delta_theta(1:n);
    x_probe = x_probe(1:n);
    response = double(response(1:n));
    response(response ~= 0 & response ~= 1) = NaN;
    valid = ~isnan(response);
    if ~any(valid)
        params = nan(1, 4);
        exit_flag = NaN;
        return
    end

    delta_theta = delta_theta(valid);
    x_probe = x_probe(valid);
    response = response(valid);

    defaults = unbinnedMLEFitDefaults();
    lb = defaults.lb;
    ub = defaults.ub;
    x0 = defaults.x0;

    if isfield(opts, 'lb') && ~isempty(opts.lb), lb = opts.lb(:); end
    if isfield(opts, 'ub') && ~isempty(opts.ub), ub = opts.ub(:); end
    if isfield(opts, 'x0') && ~isempty(opts.x0), x0 = opts.x0(:); end

    guess_rate = 0.25;
    if isfield(opts, 'guess_rate') && ~isempty(opts.guess_rate)
        guess_rate = opts.guess_rate;
    end

    alg = 'sqp';
    if isfield(opts, 'algorithm') && ~isempty(opts.algorithm)
        alg = opts.algorithm;
    end

    fmin_opts = optimoptions('fmincon', 'Algorithm', alg, 'Display', 'off');
    if isfield(opts, 'fmincon_options') && ~isempty(opts.fmincon_options)
        fmin_opts = opts.fmincon_options;
    end

    map_opts = resolveMapOptions(opts, lb, ub, x0);
    objective = @(theta) objectiveUnbinned(theta, delta_theta, x_probe, response, guess_rate, map_opts);

    try
        [theta_hat, ~, exit_flag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], fmin_opts);
        if exit_flag < 1
            params = nan(1, 4);
        else
            params = theta_hat(:)';
        end
    catch
        params = nan(1, 4);
        exit_flag = -99;
    end

end

function map_opts = resolveMapOptions(opts, lb, ub, x0)
    map_opts.use_map = false;
    if isfield(opts, 'use_map') && ~isempty(opts.use_map)
        map_opts.use_map = logical(opts.use_map);
    end

    map_opts.prior_means = x0(:)';
    if isfield(opts, 'prior_means') && ~isempty(opts.prior_means)
        map_opts.prior_means = opts.prior_means(:)';
    end

    map_opts.map_lambda = 0.5;
    if isfield(opts, 'map_lambda') && ~isempty(opts.map_lambda)
        map_opts.map_lambda = opts.map_lambda;
    end

    default_scales = ub(:)' - lb(:)';
    map_opts.param_scales = default_scales;
    if isfield(opts, 'param_scales') && ~isempty(opts.param_scales)
        map_opts.param_scales = opts.param_scales(:)';
    end

    if numel(map_opts.prior_means) ~= 4
        map_opts.prior_means = x0(:)';
    end
    if numel(map_opts.param_scales) ~= 4
        map_opts.param_scales = default_scales;
    end

    bad_scales = ~isfinite(map_opts.param_scales) | map_opts.param_scales <= 0;
    if any(bad_scales)
        fallback_scales = default_scales;
        fallback_scales(~isfinite(fallback_scales) | fallback_scales <= 0) = eps;
        map_opts.param_scales(bad_scales) = fallback_scales(bad_scales);
    end

    if ~isscalar(map_opts.map_lambda) || ~isfinite(map_opts.map_lambda)
        map_opts.map_lambda = 0.5;
    end
end

function cost = objectiveUnbinned(theta, delta_theta, x_probe, response, guess_rate, map_opts)
    nll = nllUnbinned(theta, delta_theta, x_probe, response, guess_rate);
    cost = nll;

    if map_opts.use_map
        z_scores = (theta(:)' - map_opts.prior_means) ./ map_opts.param_scales;
        penalty = map_opts.map_lambda * sum(z_scores .^ 2);
        cost = nll + penalty;
    end
end

function nll = nllUnbinned(theta, delta_theta, x_probe, response, guess_rate)
    A     = theta(1);
    w     = theta(2);
    sigma = theta(3);
    beta  = theta(4);

    mu = dogIsolated(delta_theta, A, w) + beta;
    p_psy = normcdf(x_probe, mu, sigma);
    p = (1 - guess_rate) .* p_psy + 0.5 .* guess_rate;        % S&S 2022 lapse model
    eps_log = eps;
    nll = -sum(response .* log(p + eps_log) + ...
               (1 - response) .* log(1 - p + eps_log), 'omitnan');
end
