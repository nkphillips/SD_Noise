function [mu_hat, sigma_hat, exit_flag] = fitSubjectBaselineBias(x_probe, response, opts)
% fitSubjectBaselineBias  Per-subject baseline psychometric fit.
%
% Maximum-likelihood fit of (mu, sigma) for a Gaussian psychometric with a
% Sheehan & Serences 2022 25% lapse rate, ignoring trial-by-trial Delta-theta
% structure. Used to estimate each subject's response-bias baseline before
% (optional) demeaning of x_probe in the unbinned SD pipeline.
%
%   P(CW | x) = (1 - lambda) * normcdf(x; mu, sigma) + 0.5 * lambda
%
% Inputs:
%   x_probe  - N x 1 probe offsets in degrees (signed: + = probe CW of stimulus)
%   response - N x 1 binary (1 = CW, 0 = CCW)
%   opts     - struct with optional fields .lb, .ub, .x0, .guess_rate,
%              .fmincon_options, .algorithm
%
% Outputs:
%   mu_hat, sigma_hat - estimated baseline parameters (deg)
%   exit_flag         - fmincon exit flag; NaN on insufficient data, < 1 on
%                       failed convergence

    if nargin < 3; opts = struct(); end

    x_probe = x_probe(:);
    response = double(response(:));
    n = min(numel(x_probe), numel(response));
    x_probe = x_probe(1:n);
    response = response(1:n);
    response(response ~= 0 & response ~= 1) = NaN;
    valid = ~isnan(response);
    x_probe = x_probe(valid);
    response = response(valid);

    if numel(x_probe) < 10
        mu_hat = NaN; sigma_hat = NaN; exit_flag = NaN;
        return
    end

    lb = [-15;  1];     % [mu; sigma]
    ub = [ 15; 90];
    x0 = [  0;  5];

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

    objective = @(theta) nllBaseline(theta, x_probe, response, guess_rate);

    try
        [theta_hat, ~, exit_flag] = fmincon(objective, x0, [], [], [], [], lb, ub, [], fmin_opts);
        if exit_flag < 1
            mu_hat = NaN; sigma_hat = NaN;
        else
            mu_hat = theta_hat(1);
            sigma_hat = theta_hat(2);
        end
    catch
        mu_hat = NaN; sigma_hat = NaN;
        exit_flag = -99;
    end
end

function nll = nllBaseline(theta, x_probe, response, guess_rate)
    mu = theta(1);
    sigma = theta(2);
    p_psy = normcdf(x_probe, mu, sigma);
    p = (1 - guess_rate) .* p_psy + 0.5 .* guess_rate;
    eps_log = eps;
    nll = -sum(response .* log(p + eps_log) + ...
               (1 - response) .* log(1 - p + eps_log), 'omitnan');
end
