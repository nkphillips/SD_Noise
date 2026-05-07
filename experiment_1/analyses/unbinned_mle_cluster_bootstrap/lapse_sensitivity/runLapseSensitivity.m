function results = runLapseSensitivity(tbl_trials, varargin)
% runLapseSensitivity  Sweep the fixed lapse rate lambda over a grid; report how
% per-cell point estimates of (A, w, FWHM, sigma, beta) shift across lambda values.
% Bootstrap is intentionally NOT performed -- the goal is a sensitivity diagnostic
% on the directional conclusions of the manipulation contrasts, not CIs.
%
% Name-value pairs:
%   'lambdas'         -- lambda grid (default [0.02 0.05 0.10 0.25])
%   'fit_opts'        -- struct passed to fitConditionMLE; .lb, .ub, .x0 honored
%   'use_parallel'    -- logical (default true)
%   'num_workers'     -- [] (default) uses parpool default
%   'fig_subdir'      -- e.g. '1_back'; figures saved under
%                        unbinned_mle_cluster_bootstrap/lapse_sensitivity/figures/<subdir>/
%   'contrast_labels' -- {'90%','50%','25%'}
%   'precision_labels'-- {'2°','40°','80°'}

    ip = inputParser;
    addParameter(ip, 'lambdas', [0.02 0.05 0.10 0.25], @(x) isnumeric(x) && all(x >= 0 & x < 1));
    addParameter(ip, 'fit_opts', struct(), @isstruct);
    addParameter(ip, 'use_parallel', true, @(x) islogical(x) && isscalar(x));
    addParameter(ip, 'num_workers', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(ip, 'fig_subdir', '', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'contrast_labels', {'90%', '50%', '25%'}, @(x) iscell(x) && numel(x) == 3);
    addParameter(ip, 'precision_labels', {'2°', '40°', '80°'}, @(x) iscell(x) && numel(x) == 3);
    parse(ip, varargin{:});

    lambdas = ip.Results.lambdas(:)';
    n_lambda = numel(lambdas);
    use_parallel = ip.Results.use_parallel;
    num_workers = ip.Results.num_workers;
    fit_opts = ip.Results.fit_opts;

    this_dir = fileparts(mfilename('fullpath'));
    fig_subdir = char(ip.Results.fig_subdir);
    if isempty(fig_subdir)
        fig_dir = fullfile(this_dir, 'figures');
    else
        fig_dir = fullfile(this_dir, 'figures', fig_subdir);
    end
    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    addpath(fullfile(this_dir, '..'));        % access fitConditionMLE, dogIsolated, etc.
    addpath(fullfile(this_dir, '..', '..', 'functions'));

    num_conds = 18;

    params_grid = nan(n_lambda, num_conds, 4);   % [A, w, sigma, beta]
    fwhm_grid   = nan(n_lambda, num_conds);
    exit_grid   = nan(n_lambda, num_conds);

    if use_parallel
        pool = gcp('nocreate');
        if isempty(pool)
            if isempty(num_workers); parpool('local'); else; parpool('local', num_workers); end
        end
    end

    fprintf('Running lapse sensitivity sweep over %d values: %s\n', n_lambda, num2str(lambdas));
    t0 = tic;

    for il = 1:n_lambda
        lambda = lambdas(il);
        opts = fit_opts;
        opts.guess_rate = lambda;

        per_cell_params = nan(num_conds, 4);
        per_cell_exit = nan(num_conds, 1);

        if use_parallel
            parfor c = 1:num_conds
                [m, prev, curr] = conditionSubscriptsFromIndex(c);
                cm = tbl_trials.cond_manipulation;
                man = ones(height(tbl_trials), 1);
                man(cm == 'precision') = 2;
                mask = man == m & tbl_trials.cond_prev == prev & tbl_trials.cond_curr == curr;
                [pf, ef] = fitConditionMLE( ...
                    tbl_trials.delta_theta(mask), tbl_trials.x_probe(mask), tbl_trials.response(mask), opts);
                per_cell_params(c, :) = pf;
                per_cell_exit(c) = ef;
            end
        else
            for c = 1:num_conds
                [m, prev, curr] = conditionSubscriptsFromIndex(c);
                cm = tbl_trials.cond_manipulation;
                man = ones(height(tbl_trials), 1);
                man(cm == 'precision') = 2;
                mask = man == m & tbl_trials.cond_prev == prev & tbl_trials.cond_curr == curr;
                [pf, ef] = fitConditionMLE( ...
                    tbl_trials.delta_theta(mask), tbl_trials.x_probe(mask), tbl_trials.response(mask), opts);
                per_cell_params(c, :) = pf;
                per_cell_exit(c) = ef;
            end
        end

        params_grid(il, :, :) = reshape(per_cell_params, 1, num_conds, 4);
        exit_grid(il, :) = per_cell_exit(:)';
        for c = 1:num_conds
            wc = per_cell_params(c, 2);
            if isfinite(wc) && wc > 0
                fwhm_grid(il, c) = unbinnedWtoFwhm(wc);
            end
        end
        fprintf('  lambda = %.2f done\n', lambda);
    end

    elapsed = toc(t0);
    fprintf('Lapse sensitivity sweep complete in %.1f s.\n', elapsed);

    % --- Long-form table for export / inspection ---
    rows = cell(n_lambda * num_conds, 1);
    row_idx = 0;
    for il = 1:n_lambda
        for c = 1:num_conds
            [m, prev, curr] = conditionSubscriptsFromIndex(c);
            if m == 1; mlabel = 'contrast'; else; mlabel = 'precision'; end
            row_idx = row_idx + 1;
            rows{row_idx} = {lambdas(il), mlabel, prev, curr, ...
                params_grid(il, c, 1), params_grid(il, c, 2), fwhm_grid(il, c), ...
                params_grid(il, c, 3), params_grid(il, c, 4), exit_grid(il, c)};
        end
    end
    table_data = vertcat(rows{:});
    summary_table = cell2table(table_data, ...
        'VariableNames', {'lambda', 'cond_manipulation', 'cond_prev', 'cond_curr', ...
                          'A', 'w', 'fwhm', 'sigma', 'beta', 'exit_flag'});

    % --- Cross-lambda contrast directionality summary (referenced lambda = 0.25 by default) ---
    lambda_ref = 0.25;
    if any(lambdas == lambda_ref)
        ref_idx = find(lambdas == lambda_ref, 1);
    else
        ref_idx = n_lambda;
        lambda_ref = lambdas(ref_idx);
    end
    sign_match = nan(n_lambda, num_conds, 4);
    for il = 1:n_lambda
        for c = 1:num_conds
            for k = 1:4
                v_ref = params_grid(ref_idx, c, k);
                v_il = params_grid(il, c, k);
                if isfinite(v_ref) && isfinite(v_il)
                    sign_match(il, c, k) = double(sign(v_ref) == sign(v_il));
                end
            end
        end
    end

    results = struct();
    results.lambdas = lambdas;
    results.lambda_ref = lambda_ref;
    results.params_grid = params_grid;     % n_lambda x 18 x 4
    results.fwhm_grid = fwhm_grid;         % n_lambda x 18
    results.exit_grid = exit_grid;
    results.summary_table = summary_table;
    results.sign_match = sign_match;
    results.elapsed_sec = elapsed;
    results.fig_dir = fig_dir;

    % --- Save table to CSV ---
    csv_path = fullfile(fig_dir, 'lapse_sensitivity_summary.csv');
    writetable(summary_table, csv_path);
    fprintf('Wrote %s\n', csv_path);

    % --- Render figures ---
    plotLapseSensitivity(results, fig_dir, ip.Results.contrast_labels, ip.Results.precision_labels);

end
