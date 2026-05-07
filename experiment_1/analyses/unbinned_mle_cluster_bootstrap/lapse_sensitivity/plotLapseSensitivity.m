function plotLapseSensitivity(results, fig_dir, contrast_labels, precision_labels)
% plotLapseSensitivity  Render summary figures for the lapse-rate sensitivity sweep.
%
% Produces three PDFs in fig_dir:
%   lapse_sensitivity_A_FWHM_sigma_beta.pdf   -- 4 panels: each parameter vs lambda,
%                                                one trace per cell, colored by manipulation.
%   lapse_sensitivity_grid_contrast.pdf       -- 3x3 grid for contrast manip; one panel
%                                                per (prev, curr) showing all 4 params vs lambda.
%   lapse_sensitivity_grid_precision.pdf      -- same for precision manip.

    if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

    lambdas = results.lambdas(:);
    n_lambda = numel(lambdas);
    num_conds = 18;

    A_grid     = squeeze(results.params_grid(:, :, 1));   % n_lambda x 18
    fwhm_grid  = results.fwhm_grid;                       % n_lambda x 18
    sigma_grid = squeeze(results.params_grid(:, :, 3));
    beta_grid  = squeeze(results.params_grid(:, :, 4));

    manip = zeros(num_conds, 1);
    prev_lvl = zeros(num_conds, 1);
    curr_lvl = zeros(num_conds, 1);
    for c = 1:num_conds
        [manip(c), prev_lvl(c), curr_lvl(c)] = conditionSubscriptsFromIndex(c);
    end

    blue = [0.20 0.45 0.75];
    green = [0.30 0.65 0.30];

    %% ---- Figure 1: each param vs lambda, one trace per cell ----
    fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 900 700]);
    tl = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    panels = {A_grid, 'Amplitude A (deg)';
              fwhm_grid, 'FWHM (deg)';
              sigma_grid, '\sigma (deg)';
              beta_grid, '\beta (deg)'};
    for ip = 1:4
        nexttile(tl);
        hold on;
        for c = 1:num_conds
            color = blue * (manip(c) == 1) + green * (manip(c) == 2);
            plot(lambdas, panels{ip, 1}(:, c), '-o', 'Color', color, 'LineWidth', 0.8, ...
                 'MarkerSize', 3, 'MarkerFaceColor', color);
        end
        xlabel('\lambda (lapse rate)');
        ylabel(panels{ip, 2});
        title(panels{ip, 2}, 'Interpreter', 'tex');
        set(gca, 'TickDir', 'out', 'XTick', lambdas);
        if ip == 1
            yline(0, 'k:', 'HandleVisibility', 'off');
        end
        if ip == 4
            yline(0, 'k:', 'HandleVisibility', 'off');
        end
        % legend stub (only in panel 1)
        if ip == 1
            h_blue = plot(NaN, NaN, '-o', 'Color', blue, 'MarkerFaceColor', blue);
            h_green = plot(NaN, NaN, '-o', 'Color', green, 'MarkerFaceColor', green);
            legend([h_blue h_green], {'Contrast', 'Precision'}, 'Location', 'best', 'FontSize', 8);
        end
        box off;
    end

    title(tl, sprintf('Lapse-rate sensitivity (each line = one of 18 cells; ref \\lambda = %.2f)', ...
        results.lambda_ref), 'Interpreter', 'tex');

    out = fullfile(fig_dir, 'lapse_sensitivity_A_FWHM_sigma_beta.pdf');
    exportgraphics(fig, out, 'ContentType', 'vector');
    close(fig);

    %% ---- Figures 2 & 3: 3x3 grid per manipulation ----
    for m = 1:2
        if m == 1
            mname = 'contrast'; mtitle = 'Contrast'; labels = contrast_labels; mcolor = blue;
        else
            mname = 'precision'; mtitle = 'Precision'; labels = precision_labels; mcolor = green;
        end

        fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1100 900]);
        tl = tiledlayout(3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

        for prev = 1:3
            for curr = 1:3
                cidx = (m - 1) * 9 + (prev - 1) * 3 + curr;
                nexttile(tl);
                hold on;

                % Each param normalized to its lambda_ref value for visual comparability
                ref_idx = find(lambdas == results.lambda_ref, 1);
                if isempty(ref_idx); ref_idx = numel(lambdas); end

                A_trace     = A_grid(:, cidx);
                fwhm_trace  = fwhm_grid(:, cidx);
                sigma_trace = sigma_grid(:, cidx);
                beta_trace  = beta_grid(:, cidx);

                yyaxis left;
                p1 = plot(lambdas, A_trace, '-o', 'Color', mcolor, 'LineWidth', 1.2, 'MarkerFaceColor', mcolor);
                p4 = plot(lambdas, beta_trace, '-s', 'Color', [0.55 0.27 0.55], 'LineWidth', 1.0, 'MarkerFaceColor', [0.55 0.27 0.55]);
                ylabel('A, \beta (deg)');
                yline(0, 'k:', 'HandleVisibility', 'off');

                yyaxis right;
                p2 = plot(lambdas, fwhm_trace, '-^', 'Color', [0.85 0.55 0.20], 'LineWidth', 1.0, 'MarkerFaceColor', [0.85 0.55 0.20]);
                p3 = plot(lambdas, sigma_trace, '-d', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'MarkerFaceColor', [0.45 0.45 0.45]);
                ylabel('FWHM, \sigma (deg)');

                set(gca, 'TickDir', 'out', 'XTick', lambdas);
                xlabel('\lambda');
                title(sprintf('%s \\rightarrow %s', labels{prev}, labels{curr}), 'Interpreter', 'tex');
                if prev == 1 && curr == 1
                    legend([p1 p4 p2 p3], {'A', '\beta', 'FWHM', '\sigma'}, ...
                        'Location', 'best', 'FontSize', 7, 'Interpreter', 'tex');
                end
                box off;
            end
        end
        title(tl, sprintf('%s: parameter sensitivity to \\lambda', mtitle), 'Interpreter', 'tex');

        out = fullfile(fig_dir, sprintf('lapse_sensitivity_grid_%s.pdf', mname));
        exportgraphics(fig, out, 'ContentType', 'vector');
        close(fig);
    end

end
