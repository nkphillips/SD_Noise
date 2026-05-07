function [params, exit_flag] = fitConditionMLE(delta_theta, x_probe, response, opts)
% fitConditionMLE  Trial-level MLE for S&S DoG + baseline + Gaussian psychometric with 25% guess rate.
%
%   mu_t = dogIsolated(dtheta, A, w) + beta                                (S&S Eq. 1)
%   P_t  = (1 - lambda) * normcdf(x_probe, mu_t, sigma) + 0.5 * lambda     (lambda = 0.25)
%
% params = [A, w, sigma, beta] (1x4); w in 1/deg. NaN(1,4) if no trials.
% exit_flag is the fmincon exit flag (NaN on missing data, < 1 means failed/at-bound stop).

    if nargin < 4
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

    % Default bounds match legacy (S&S parameterization, w in 1/deg).
    % FWHM bounds 8 deg .. 120 deg under FWHM = 2*sqrt(log(2))/w.
    fwhm_min_deg = 8;
    fwhm_max_deg = 120;
    w_lb = (2 * sqrt(log(2))) / fwhm_max_deg;
    w_ub = (2 * sqrt(log(2))) / fwhm_min_deg;

    lb = [-30; w_lb;   1; -5];     % [A; w; sigma; beta] -- sigma_lb 1 deg (was 0.1: avoids lapse-sigma confound at the corner)
    ub = [ 30; w_ub;  90;  5];
    x0 = [  1;  0.1;   5;  0];

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

    objective = @(theta) nllUnbinned(theta, delta_theta, x_probe, response, guess_rate);

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
