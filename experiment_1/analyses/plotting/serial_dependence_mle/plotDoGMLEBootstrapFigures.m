function plotDoGMLEBootstrapFigures(ps, fig_dir, curve_boot, grid, params_boot, ci_pct, varargin)
% plotDoGMLEBootstrapFigures  3×3 grids per manipulation: binned μ MLE from trials + bootstrap μ(Δθ) + CI.
% curve_boot stores isolated DoG on grid; params_boot(:,c,4) adds β so ribbons match μ = DoG+β.
% Styling follows Super Subj Bias with DoG * (SD_Noise_Analyses_And_Figures.m): contrast=blue,
% precision=green; α, w, β on panels (σ only in table/results);
% R²_Δθ (binned visual fit) and McFadden R² (trial-level likelihood fit).

    if nargin < 6 || isempty(ci_pct)
        ci_pct = [2.5, 97.5];
    end

    ip = inputParser;
    ip.FunctionName = mfilename;
    addParameter(ip, 'contrast_labels', {'90%', '50%', '25%'}, @(x) iscell(x) && numel(x) == 3);
    addParameter(ip, 'precision_labels', {'2°', '40°', '80°'}, @(x) iscell(x) && numel(x) == 3);
    addParameter(ip, 'mu_bounds', [-20, 20], @(x) isnumeric(x) && numel(x) == 2);
    addParameter(ip, 'overlay', [], @(x) isempty(x) || isstruct(x));
    addParameter(ip, 'admitted', [], @(x) isempty(x) || (islogical(x) || isnumeric(x)));
    addParameter(ip, 'curve_lo', [], @(x) isempty(x) || isnumeric(x));   % num_conds x n_grid (S&S BCa or percentile)
    addParameter(ip, 'curve_hi', [], @(x) isempty(x) || isnumeric(x));
    addParameter(ip, 'ci_method', 'percentile', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'bootstrap_unit', 'subject', @(x) ischar(x) || isstring(x));
    addParameter(ip, 'folded_delta_theta', false, @(x) islogical(x) && isscalar(x));
    parse(ip, varargin{:});

    contrast_lbl = ip.Results.contrast_labels;
    precision_lbl = ip.Results.precision_labels;
    mu_b = ip.Results.mu_bounds(:)';
    rb_lo = mu_b(1);
    rb_hi = mu_b(2);
    overlay = ip.Results.overlay;
    admitted = ip.Results.admitted;                    % B x num_conds logical (or empty)
    curve_lo_in = ip.Results.curve_lo;                 % overrides in-figure percentile when supplied
    curve_hi_in = ip.Results.curve_hi;
    ci_method = lower(char(ip.Results.ci_method));
    bootstrap_unit = lower(char(ip.Results.bootstrap_unit));
    folded_delta_theta = ip.Results.folded_delta_theta;
    if strcmpi(bootstrap_unit, 'subject')
        ci_unit_label = 'subject-cluster';
    else
        ci_unit_label = 'trial-bootstrap';
    end
    ci_nominal_pct = ci_pct(2) - ci_pct(1);
    ci_label = sprintf('%.0f%% %s %s CI', ci_nominal_pct, local_ciLabel(ci_method), ci_unit_label);

    axes_fs = 14;
    tick_fs = 13;
    try
        axes_fs = ps.axes_label_font_size;
        tick_fs = ps.axes_tick_font_size;
    catch
    end

    num_levels = 3;
    n_boot = size(curve_boot, 1);
    n_grid = numel(grid);
    has_overlay = ~isempty(overlay) && isfield(overlay, 'params_point') && isfield(overlay, 'mu_bin_mle');

    for m = 1:2
        if m == 1
            labels = contrast_lbl;
            mname_file = 'contrast';
            cond_name_title = 'Contrast';
            cond_color = ps.colors.blue;
        else
            labels = precision_lbl;
            mname_file = 'precision';
            cond_name_title = 'Precision';
            cond_color = ps.colors.green;
        end

        dog_color = ps.colors.black;

        % --- Global y limits: bootstrap μ(Δθ)=DoG+β + binned μ MLE ---
        all_y = [];
        for prev = 1:num_levels
            for curr = 1:num_levels
                cidx = conditionIndexSerialDependence(m, prev, curr);
                slice_iso = reshape(curve_boot(:, cidx, :), [n_boot, n_grid]);
                beta_b = reshape(params_boot(:, cidx, 4), [n_boot, 1]);
                slice = slice_iso + beta_b;
                if ~isempty(admitted)
                    adm_c = logical(admitted(:, cidx));
                    if any(adm_c)
                        slice = slice(adm_c, :);
                    end
                end
                if size(slice, 1) <= 1
                    mu = slice;
                    lo = slice;
                    hi = slice;
                else
                    mu = meanColsOmitNan(slice);
                    if ~isempty(curve_lo_in) && ~isempty(curve_hi_in)
                        lo = curve_lo_in(cidx, :);
                        hi = curve_hi_in(cidx, :);
                    else
                        lo = prctileColsOmitNan(slice, ci_pct(1));
                        hi = prctileColsOmitNan(slice, ci_pct(2));
                    end
                end
                all_y = [all_y, mu(:)', lo(:)', hi(:)']; %#ok<AGROW>

                if has_overlay && cidx <= numel(overlay.mu_bin_mle) && ~isempty(overlay.mu_bin_mle{cidx})
                    mb = overlay.mu_bin_mle{cidx};
                    yy = mb.mu_deg(:)';
                    if folded_delta_theta
                        [~, yy] = local_mirrorFoldedMuBins(mb.delta_centers, yy);
                    end
                    yy = yy(isfinite(yy));
                    all_y = [all_y, yy(:)']; %#ok<AGROW>
                end
            end
        end
        all_y = all_y(isfinite(all_y));
        if isempty(all_y)
            y_lims = [rb_lo, rb_hi];
        else
            y_min = min(all_y);
            y_max = max(all_y);
            if y_min == y_max
                y_pad = max(abs(y_min), 1) * 0.1;
            else
                y_pad = 0.08 * (y_max - y_min);
            end
            y_lims = [min(rb_lo, y_min - y_pad), max(rb_hi, y_max + y_pad)];
        end

        fig = figure('Color', ps.figure_color, 'Visible', 'off', ...
            'Units', 'inches', 'Position', [1 1 11 11], ...
            'PaperUnits', 'inches', 'PaperSize', [11 11], ...
            'PaperPositionMode', 'manual', 'PaperPosition', [0 0 11 11]);
        tl = tiledlayout(num_levels, num_levels, 'Padding', 'compact', 'TileSpacing', 'compact');

        for prev = 1:num_levels
            for curr = 1:num_levels
                cidx = conditionIndexSerialDependence(m, prev, curr);
                ax = nexttile(tl);

                fg_title = sprintf('%s \\rightarrow %s', labels{prev}, labels{curr});

                hold(ax, 'on');

                % 1) Binned μ MLE from same trials as pooled fit (σ fixed to pooled σ̂); comparable to μ(Δθ)=DoG+β
                if has_overlay && cidx <= numel(overlay.mu_bin_mle) && ~isempty(overlay.mu_bin_mle{cidx})
                    mb = overlay.mu_bin_mle{cidx};
                    dx = mb.delta_centers;
                    ydat = mb.mu_deg;
                    if folded_delta_theta
                        [dx, ydat] = local_mirrorFoldedMuBins(dx, ydat);
                    end
                    idx = isfinite(ydat);
                    if any(idx)
                        plot(ax, dx(idx), ydat(idx), '-', 'LineWidth', ps.line_width, 'Color', cond_color);
                    end
                end

                slice_iso = reshape(curve_boot(:, cidx, :), [n_boot, n_grid]);
                beta_b = reshape(params_boot(:, cidx, 4), [n_boot, 1]);
                slice = slice_iso + beta_b;
                if ~isempty(admitted)
                    adm_c = logical(admitted(:, cidx));
                    if any(adm_c)
                        slice = slice(adm_c, :);
                    end
                end
                if size(slice, 1) <= 1
                    mu = slice;
                    lo = slice;
                    hi = slice;
                else
                    mu = meanColsOmitNan(slice);
                    if ~isempty(curve_lo_in) && ~isempty(curve_hi_in)
                        lo = curve_lo_in(cidx, :);
                        hi = curve_hi_in(cidx, :);
                    else
                        lo = prctileColsOmitNan(slice, ci_pct(1));
                        hi = prctileColsOmitNan(slice, ci_pct(2));
                    end
                end

                xf = [grid(:); flipud(grid(:))];
                yf = [lo(:); flipud(hi(:))];
                fill(ax, xf, yf, cond_color, 'FaceAlpha', 0.15, 'EdgeColor', 'none');

                plot(ax, grid, mu, '--', 'LineWidth', ps.line_width, 'Color', dog_color);

                title(ax, fg_title, 'Interpreter', 'tex', 'FontSize', tick_fs);

                axis(ax, 'square');
                box(ax, 'off');
                set(ax, 'TickDir', 'out', 'TickLength', [ps.tick_length, ps.tick_length], 'FontSize', tick_fs);

                xlim(ax, [-90, 90]);
                ylim(ax, y_lims);

                line(ax, [min(xlim(ax)), max(xlim(ax))], [0, 0], 'LineWidth', 1, 'Color', 'k');
                line(ax, [0, 0], y_lims, 'LineWidth', 1, 'Color', 'k');

                xticks(ax, -90:45:90);
                xtickangle(ax, 0);

                if prev == num_levels
                    xlabel(ax, '\Delta\theta (°)', 'FontSize', axes_fs, 'Interpreter', 'tex');
                else
                    set(ax, 'XTickLabel', []);
                end

                if curr == 1
                    ylabel(ax, 'Bias (°)', 'FontSize', axes_fs, 'Interpreter', 'none');
                end

                yt = floor(y_lims(1) / 5) * 5 : 5 : ceil(y_lims(2) / 5) * 5;
                if isempty(yt)
                    yticks(ax, 'auto');
                else
                    yticks(ax, yt);
                end

                % Annotations: binned Δθ R² + McFadden R² top-left; DoG + baseline bottom-right.
                if has_overlay && cidx <= size(overlay.params_point, 1)
                    pf = overlay.params_point(cidx, :);
                    r2_lines = {};
                    if isfield(overlay, 'r2_delta_bins') && ~isnan(overlay.r2_delta_bins(cidx))
                        r2_lines{end+1} = sprintf('R² Δθ = %.2f', overlay.r2_delta_bins(cidx)); %#ok<AGROW>
                    end
                    if isfield(overlay, 'r2_mcfadden') && ~isnan(overlay.r2_mcfadden(cidx))
                        r2_lines{end+1} = sprintf('R² MF = %.2f', overlay.r2_mcfadden(cidx)); %#ok<AGROW>
                    end
                    if ~isempty(r2_lines)
                        text(ax, 0.05, 0.90, strjoin(r2_lines, newline), ...
                            'Units', 'normalized', 'HorizontalAlignment', 'left', ...
                            'VerticalAlignment', 'top', 'FontWeight', 'normal', 'FontSize', 9, ...
                            'Color', ps.colors.black, 'Interpreter', 'none');
                    end
                    if ~all(isnan([pf(1), pf(2), pf(4)]))
                        if isfinite(pf(2)) && pf(2) > 0
                            fwhm_pf = serialDependenceWtoFwhm(pf(2));
                        else
                            fwhm_pf = NaN;
                        end
                        param_text = sprintf([ ...
                            'A = %.2f°' newline ...
                            'FWHM = %.1f°' newline ...
                            'β = %.2f°' ], pf(1), fwhm_pf, pf(4));
                        text(ax, 0.95, 0.10, param_text, ...
                            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                            'VerticalAlignment', 'bottom', 'FontWeight', 'normal', 'FontSize', 10, ...
                            'Color', ps.colors.black, 'Interpreter', 'none');
                    end
                end
            end
        end

        if folded_delta_theta
            title_prefix = 'Folded \Delta\theta Serial Dependence';
        else
            title_prefix = 'Serial Dependence';
        end
        title(tl, sprintf('%s: %s (%s)', title_prefix, cond_name_title, ci_label), 'FontSize', axes_fs + 1);

        out_pdf = fullfile(fig_dir, sprintf('serial_dependence_%s.pdf', mname_file));
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
        close(fig);
    end

end

function s = local_ciLabel(ci_method)
    ci_method = lower(char(ci_method));
    switch ci_method
        case 'bca'
            s = 'BCa';
        case 'percentile'
            s = 'percentile';
        otherwise
            s = char(ci_method);
    end
end

function [dx_out, mu_out] = local_mirrorFoldedMuBins(dx, mu)
    dx = dx(:);
    mu = mu(:);

    keep = isfinite(dx) & isfinite(mu) & dx >= 0;
    dx_pos = dx(keep);
    mu_pos = mu(keep);

    [dx_pos, ord] = sort(dx_pos);
    mu_pos = mu_pos(ord);

    is_zero = abs(dx_pos) < eps;
    dx_neg = -flipud(dx_pos(~is_zero));
    mu_neg = -flipud(mu_pos(~is_zero));

    dx_out = [dx_neg; dx_pos];
    mu_out = [mu_neg; mu_pos];
end
