function plotResponseBias(delta, mu, plt_opts, cond, mu_lo, mu_hi, baseline)

if nargin < 5
    mu_lo = [];
    mu_hi = [];
end
if nargin < 7
    baseline = [];
end

% Optionally subtract baseline (b) if toggle is on and baseline provided
    if isfield(plt_opts, 'rb_subtract_baseline') && plt_opts.rb_subtract_baseline && ~isempty(baseline)
    mu = mu - baseline;
    if ~isempty(mu_lo) && ~isempty(mu_hi)
        mu_lo = mu_lo - baseline;
        mu_hi = mu_hi - baseline;
    end
end

if nargin >= 6 && ~isempty(mu_lo) && ~isempty(mu_hi)
    delta = delta(:)';
    mu = mu(:)';
    mu_lo = mu_lo(:)';
    mu_hi = mu_hi(:)';
    idx = isfinite(mu) & isfinite(mu_lo) & isfinite(mu_hi);
    if any(idx)
        err_upper = mu_hi(idx) - mu(idx);
        err_lower = mu(idx) - mu_lo(idx);
        err = [err_upper; err_lower];
        if cond == 1
                lineProps = {'-','Color', plt_opts.colors.blue, 'LineWidth', plt_opts.line_width};
                lineColor = plt_opts.colors.blue;
        else
                lineProps = {'-','Color', plt_opts.colors.green, 'LineWidth', plt_opts.line_width};
                lineColor = plt_opts.colors.green;
        end
        shadedErrorBar(delta(idx), mu(idx), err, 'lineProps', lineProps, 'transparent', true, 'patchSaturation', 0.5);
        hold on;
            plot(delta, mu, 'LineWidth', plt_opts.line_width, 'Color', lineColor);
        hold on;
        return
    end
end

if cond == 1
        plot(delta, mu, 'LineWidth', plt_opts.line_width, 'Color', plt_opts.colors.blue);
elseif cond == 2
        plot(delta, mu, 'LineWidth', plt_opts.line_width, 'Color', plt_opts.colors.green);
end

hold on;
end