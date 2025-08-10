function result = processBootstrapWindowTask(task_data, p, bootstrap, toggles)
% processBootstrapWindowTask
% Bootstrap one window's trials to compute CIs for RB mu and performance.

    po = task_data.po(:);
    y  = task_data.y(:);
    B  = bootstrap.B;
    prc = bootstrap.ci;
    if ~isfield(bootstrap,'min_trials') || isempty(bootstrap.min_trials)
        bootstrap.min_trials = 0;
    end

    n = min(numel(po), numel(y));
    po = po(1:n); y = y(1:n);

    if n <= bootstrap.min_trials
        result.mu_lo = NaN; result.mu_hi = NaN;
        result.pc_lo = NaN; result.pc_hi = NaN;
        result.pccw_lo = NaN; result.pccw_hi = NaN;
        return
    end

    mu_b   = nan(B,1);
    pc_b   = nan(B,1);
    pccw_b = nan(B,1);

    for b = 1:B
        idx = randi(n, n, 1);
        po_b = po(idx);
        y_b  = y(idx);
        try
            [~, ~, params_est_b] = estimateResponseBias([po_b, y_b], p);
            mu_b(b) = params_est_b(1);
        catch ME %#ok<NASGU>
            mu_b(b) = NaN;
            if isfield(toggles,'disp_on') && toggles.disp_on && b == 1
                % warning suppressed to avoid flooding; uncomment to debug
                % warning('RB bootstrap fit failed: %s', ME.message);
            end
        end
        try
            is_correct = (y_b == (po_b > 0));
            pc_b(b) = 100 * mean(is_correct);
        catch
            pc_b(b) = NaN;
        end
        try
            pccw_b(b) = 100 * mean(1 - y_b);
        catch
            pccw_b(b) = NaN;
        end
    end

    mu_c = mu_b(~isnan(mu_b));
    if ~isempty(mu_c)
        q = prctile(mu_c, prc);
        result.mu_lo = q(1); result.mu_hi = q(2);
    else
        result.mu_lo = NaN; result.mu_hi = NaN;
    end

    pc_c = pc_b(~isnan(pc_b));
    if ~isempty(pc_c)
        q = prctile(pc_c, prc);
        result.pc_lo = q(1); result.pc_hi = q(2);
    else
        result.pc_lo = NaN; result.pc_hi = NaN;
    end

    pccw_c = pccw_b(~isnan(pccw_b));
    if ~isempty(pccw_c)
        q = prctile(pccw_c, prc);
        result.pccw_lo = q(1); result.pccw_hi = q(2);
    else
        result.pccw_lo = NaN; result.pccw_hi = NaN;
    end

end


