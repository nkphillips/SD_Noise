function m = binaryFitR2Metrics(y, p_hat)
% binaryFitR2Metrics  Common R²-style summaries for binary outcomes vs predicted P(success).
%
%   y      : 0/1 (or logical)
%   p_hat  : model-predicted P(y=1), same length as y
%
%   .efron    — Efron (1978): 1 - sum((y-p̂)²)/sum((y-ȳ)²); same as calcR2.m on binary trials.
%              Often ~0.02–0.2 for trial-level psychometrics even when NLL fit is good.
%   .tjur     — Tjur (2009): mean(p̂|y=1) - mean(p̂|y=0); discrimination coefficient, typically > Efron.
%   .mcfadden — McFadden pseudo-R²: 1 - LL_model/LL_null (intercept-only); can be negative.

    y = double(y(:));
    p_hat = double(p_hat(:));
    n = min(numel(y), numel(p_hat));
    y = y(1:n);
    p_hat = p_hat(1:n);
    valid = ~(isnan(y) | isnan(p_hat));
    y = y(valid);
    p_hat = p_hat(valid);

    if isempty(y)
        m.efron = nan;
        m.tjur = nan;
        m.mcfadden = nan;
        return
    end

    m.efron = calcR2(y, p_hat);

    y1 = y == 1;
    y0 = y == 0;
    if any(y1) && any(y0)
        m.tjur = mean(p_hat(y1)) - mean(p_hat(y0));
    else
        m.tjur = nan;
    end

    eps_log = eps;
    pc = min(max(p_hat, eps_log), 1 - eps_log);
    LL_m = sum(y .* log(pc) + (1 - y) .* log(1 - pc));
    p0 = mean(y);
    p0 = min(max(p0, eps_log), 1 - eps_log);
    LL_0 = sum(y .* log(p0) + (1 - y) .* log(1 - p0));
    if LL_0 ~= 0 && isfinite(LL_m) && isfinite(LL_0)
        m.mcfadden = 1 - LL_m / LL_0;
    else
        m.mcfadden = nan;
    end

end
