function out = writeUnbinnedStatisticalResultsHtml(sd_noise, output_root, opts)
% writeUnbinnedStatisticalResultsHtml  Interactive HTML summary for saved analyses.

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'n_back_list') || isempty(opts.n_back_list)
        opts.n_back_list = sd_noise.config.n_back_list;
    end
    if ~exist(output_root, 'dir')
        mkdir(output_root);
    end

    all_tests = table();
    for i_nb = 1:numel(opts.n_back_list)
        n_back = opts.n_back_list(i_nb);
        key = sprintf('n%d', n_back);
        if ~isfield(sd_noise.results, key) || isempty(sd_noise.results.(key))
            continue
        end

        tests_i = computeTargetedDoGHypothesisTests(sd_noise.results.(key), n_back);
        if isempty(tests_i) || height(tests_i) == 0
            continue
        end

        super_dir = fullfile(output_root, sprintf('%d_back', n_back), 'super_subject');
        if ~exist(super_dir, 'dir')
            mkdir(super_dir);
        end
        writetable(tests_i, fullfile(super_dir, 'targeted_dog_hypothesis_tests.csv'));
        if isempty(all_tests) || width(all_tests) == 0
            all_tests = tests_i;
        else
            all_tests = [all_tests; tests_i]; %#ok<AGROW>
        end
    end

    if ~isempty(all_tests) && height(all_tests) > 0
        writetable(all_tests, fullfile(output_root, 'targeted_dog_hypothesis_tests_all_nback.csv'));
    end

    html_path = fullfile(output_root, 'statistical_results.html');
    local_writeHtml(html_path, sd_noise, opts, all_tests);

    out = struct('html_path', html_path, 'targeted_tests', all_tests);
end

function local_writeHtml(html_path, sd_noise, opts, all_tests)
    fid = fopen(html_path, 'w');
    if fid == -1
        error('writeUnbinnedStatisticalResultsHtml:openFailed', ...
            'Could not open HTML output for writing: %s', html_path);
    end
    cleanup_obj = onCleanup(@() fclose(fid));

    fprintf(fid, '<!doctype html>\n<html lang="en">\n<head>\n');
    fprintf(fid, '<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n');
    fprintf(fid, '<title>SD Noise Statistical Results</title>\n');
    fprintf(fid, '<style>\n');
    fprintf(fid, '%s\n', local_css());
    fprintf(fid, '</style>\n</head>\n<body>\n');

    fprintf(fid, '<main class="page">\n');
    fprintf(fid, '<header class="hero">\n');
    fprintf(fid, '<div><h1>Targeted Statistical Results</h1>\n');
    fprintf(fid, '<p>Analysis <code>%s</code>. Targeted bootstrap tests for whether contrast primarily changes DoG amplitude and orientation precision primarily changes DoG range. Rows are promoted only when the active bootstrap CI passes the diagnostic checks.</p></div>\n', ...
        local_escape(sd_noise.meta.analysis_datetime));
    fprintf(fid, '<a class="button" href="%s">All targeted tests CSV</a>\n', ...
        local_href('targeted_dog_hypothesis_tests_all_nback.csv'));
    fprintf(fid, '</header>\n');

    fprintf(fid, '<section class="controls" aria-label="Table filters">\n');
    fprintf(fid, '<label>n-back <select id="nbackFilter"><option value="all">All</option>');
    for i = 1:numel(opts.n_back_list)
        fprintf(fid, '<option value="%d">%d-back</option>', opts.n_back_list(i), opts.n_back_list(i));
    end
    fprintf(fid, '</select></label>\n');
    fprintf(fid, '<label>Parameter <select id="paramFilter"><option value="all">All</option><option value="A">Amplitude A</option><option value="FWHM">Range FWHM</option></select></label>\n');
    fprintf(fid, '<label>Family <select id="familyFilter"><option value="all">All</option><option value="within_manipulation">Within manipulation</option><option value="between_manipulation">Between manipulation</option></select></label>\n');
    fprintf(fid, '</section>\n');

    fprintf(fid, '<section class="grid">\n');
    fprintf(fid, '<article class="card"><h2>Primary Questions</h2><ul>');
    fprintf(fid, '<li><strong>Bias strength:</strong> does reducing contrast change DoG amplitude <code>A</code>?</li>');
    fprintf(fid, '<li><strong>Bias range:</strong> does reducing orientation precision change DoG range, reported as <code>FWHM</code>?</li>');
    fprintf(fid, '<li><strong>Dissociation:</strong> is the contrast effect larger for <code>A</code>, while the precision effect is larger for <code>FWHM</code>?</li>');
    fprintf(fid, '</ul></article>\n');
    fprintf(fid, '<article class="card"><h2>Interpretation Notes</h2><ul>');
    fprintf(fid, '<li>Effects are coded as lower information minus higher information: <code>L3 - L1</code>.</li>');
    fprintf(fid, '<li><code>FWHM</code> is used instead of raw <code>w</code> because it is directly interpretable as DoG range in degrees.</li>');
    fprintf(fid, '<li>The <code>average_prev_curr</code> rows average the previous-level and current-level spans; the other rows show each axis separately.</li>');
    fprintf(fid, '</ul></article>\n');
    fprintf(fid, '</section>\n');

    local_writeInterpretationSection(fid, all_tests);
    local_writeDiagnosticsSection(fid, all_tests);

    fprintf(fid, '<section class="card full"><h2>Targeted DoG Hypothesis Tests</h2>\n');
    fprintf(fid, '<div class="table-wrap"><table id="testsTable"><thead><tr>');
    headers = {'n-back', 'Question', 'Measure', 'Level Axis', 'Test Type', 'Tested Effect', ...
        'Estimate', 'Bootstrap median', '95% active CI', '95% percentile CI', ...
        'Admitted', 'Diagnostics', 'p', 'Holm p', 'Supports?'};
    for i = 1:numel(headers)
        fprintf(fid, '<th>%s</th>', headers{i});
    end
    fprintf(fid, '</tr></thead><tbody>\n');

    if ~isempty(all_tests) && height(all_tests) > 0
        for i = 1:height(all_tests)
            fprintf(fid, '%s\n', local_testRow(all_tests(i, :)));
        end
    end
    fprintf(fid, '</tbody></table></div></section>\n');

    fprintf(fid, '<section class="card full"><h2>Figure and CSV Links</h2>\n');
    fprintf(fid, '<div class="link-grid">\n');
    local_writeAcrossNBackLinkBlock(fid);
    for i = 1:numel(opts.n_back_list)
        local_writeLinkBlock(fid, opts.n_back_list(i));
    end
    fprintf(fid, '</div></section>\n');

    fprintf(fid, '</main>\n<script>\n%s\n</script>\n</body>\n</html>\n', local_js());
end

function local_writeInterpretationSection(fid, all_tests)
    fprintf(fid, '<section class="card full"><h2>Interpretation Guide</h2>\n');
    fprintf(fid, '<p class="section-copy">These takeaways use Holm-corrected <code>p &lt; .05</code> within the targeted test family and require the diagnostic validity flag. Treat axis-specific results as telling you where the modulation appears: previous stimulus level, current stimulus level, or their average.</p>\n');
    fprintf(fid, '<div class="interpretation-grid">\n');

    n_cards = 0;
    if ~isempty(all_tests) && height(all_tests) > 0
        keep = isfinite(all_tests.p_holm) & all_tests.p_holm < 0.05 & local_getLogicalColumn(all_tests, 'valid_for_inference', true);
        sig_tests = all_tests(keep, :);
        for i = 1:height(sig_tests)
            local_writeInterpretationCard(fid, sig_tests(i, :));
            n_cards = n_cards + 1;
        end
    end

    if n_cards == 0
        fprintf(fid, '<article class="takeaway neutral"><h3>No diagnostic-valid Holm-corrected targeted effects</h3><p>No targeted DoG hypothesis test both survived Holm correction at <code>p &lt; .05</code> and passed the bootstrap diagnostics. Inspect flagged rows and figures as descriptive evidence only.</p></article>\n');
    end

    fprintf(fid, '</div>\n');
    fprintf(fid, '<div class="interpretation-notes"><h3>How to read the results</h3><ul>');
    fprintf(fid, '<li>A significant <strong>precision effect on FWHM</strong> means serial-dependence range changes as orientation precision decreases.</li>');
    fprintf(fid, '<li>A significant <strong>precision effect &gt; contrast effect on FWHM</strong> supports a dissociation: precision modulates bias range more strongly than contrast over the sampled level range.</li>');
    fprintf(fid, '<li>A significant <strong>contrast effect on A</strong> means serial-dependence strength changes as contrast decreases.</li>');
    fprintf(fid, '<li>A significant <strong>contrast effect &gt; precision effect on A</strong> supports the complementary dissociation: contrast modulates bias strength more strongly than precision.</li>');
    fprintf(fid, '<li>These tests compare within-manipulation spans from high to low information. They do not claim that physical contrast and precision levels are numerically equivalent.</li>');
    fprintf(fid, '<li>Rows with low admitted bootstrap fraction, extreme <code>z0</code>, estimate/percentile-CI mismatch, or active/percentile disagreement are diagnostic warnings, not promoted evidence.</li>');
    fprintf(fid, '</ul></div></section>\n');
end

function local_writeDiagnosticsSection(fid, all_tests)
    fprintf(fid, '<section class="card full"><h2>Bootstrap Diagnostics</h2>\n');
    fprintf(fid, '<p class="section-copy">Rows below failed at least one bootstrap sanity check. They are retained for transparency but are not treated as hypothesis-supporting results.</p>\n');

    if isempty(all_tests) || height(all_tests) == 0 || ~ismember('valid_for_inference', all_tests.Properties.VariableNames)
        fprintf(fid, '<p>No diagnostic columns were available.</p></section>\n');
        return
    end

    flagged = all_tests(~all_tests.valid_for_inference, :);
    if isempty(flagged) || height(flagged) == 0
        fprintf(fid, '<p>No targeted rows failed the bootstrap diagnostics.</p></section>\n');
        return
    end

    fprintf(fid, '<div class="table-wrap"><table><thead><tr>');
    headers = {'n-back', 'Question', 'Measure', 'Level Axis', 'Estimate', '95% active CI', '95% percentile CI', 'Admitted', 'Diagnostics'};
    for i = 1:numel(headers)
        fprintf(fid, '<th>%s</th>', headers{i});
    end
    fprintf(fid, '</tr></thead><tbody>\n');

    for i = 1:height(flagged)
        row = flagged(i, :);
        fprintf(fid, ['<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td>' ...
            '<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n'], ...
            row.n_back(1), ...
            local_escape(char(row.hypothesis(1))), ...
            local_escape(char(row.parameter_label(1))), ...
            local_escape(local_axisLabel(char(row.axis(1)))), ...
            local_num(row.estimate(1)), ...
            local_ci(row.bca_lo(1), row.bca_hi(1)), ...
            local_ci(row.pc_lo(1), row.pc_hi(1)), ...
            local_escape(local_admitLabel(row)), ...
            local_escape(local_diagnosticLabel(row)));
    end
    fprintf(fid, '</tbody></table></div></section>\n');
end

function local_writeInterpretationCard(fid, row)
    support = logical(row.supports_hypothesis(1));
    if support
        cls = 'yes';
        tag = 'Supports target';
    else
        cls = 'warn';
        tag = 'Significant, check direction';
    end

    title_text = local_interpretationTitle(row);
    body_text = local_interpretationBody(row);
    fprintf(fid, '<article class="takeaway %s">', cls);
    fprintf(fid, '<div class="tag">%s</div>', local_escape(tag));
    fprintf(fid, '<h3>%s</h3>', local_escape(title_text));
    fprintf(fid, '<p>%s</p>', local_escape(body_text));
    fprintf(fid, '<p class="statline">Estimate %s, 95%% active CI %s, percentile CI %s, Holm p %s</p>', ...
        local_escape(local_num(row.estimate(1))), ...
        local_escape(local_ci(row.bca_lo(1), row.bca_hi(1))), ...
        local_escape(local_ci(row.pc_lo(1), row.pc_hi(1))), ...
        local_escape(char(row.p_holm_label(1))));
    fprintf(fid, '</article>\n');
end

function s = local_interpretationTitle(row)
    nback = row.n_back(1);
    question = char(row.hypothesis(1));
    axis_label = local_axisLabel(char(row.axis(1)));
    s = sprintf('%d-back: %s (%s)', nback, question, axis_label);
end

function s = local_interpretationBody(row)
    question = char(row.hypothesis(1));
    axis_label = lower(local_axisLabel(char(row.axis(1))));
    if contains(question, 'Precision effect > contrast effect on bias range')
        s = sprintf(['Orientation precision has a larger effect on DoG range than contrast for the %s axis. ' ...
            'This is evidence that degrading precision broadens or changes the spatial/orientation range of serial dependence more than reducing contrast.'], axis_label);
    elseif contains(question, 'Contrast effect > precision effect on bias strength')
        s = sprintf(['Contrast has a larger effect on DoG amplitude than precision for the %s axis. ' ...
            'This is evidence that reducing signal strength alters the magnitude of serial dependence more than degrading precision.'], axis_label);
    elseif contains(question, 'Precision effect on bias range')
        s = sprintf(['DoG range changes across precision levels for the %s axis. ' ...
            'This supports the idea that orientation uncertainty affects how broadly past stimuli influence current judgments.'], axis_label);
    elseif contains(question, 'Contrast effect on bias strength')
        s = sprintf(['DoG amplitude changes across contrast levels for the %s axis. ' ...
            'This supports the idea that signal strength affects the magnitude of attractive bias.'], axis_label);
    elseif contains(question, 'Precision effect on bias strength')
        s = sprintf(['DoG amplitude also changes across precision levels for the %s axis. ' ...
            'This suggests the strength of bias may not be exclusively contrast-dependent.'], axis_label);
    elseif contains(question, 'Contrast effect on bias range')
        s = sprintf(['DoG range also changes across contrast levels for the %s axis. ' ...
            'This suggests the range of bias may not be exclusively precision-dependent.'], axis_label);
    elseif contains(question, 'boundary interaction')
        s = sprintf(['The endpoint current-level effect differs between high- and low-information previous states. ' ...
            'This is evidence for a non-additive previous-current state combination, conditional on the bootstrap diagnostics.']);
    else
        s = char(row.comparison(1));
    end
end

function row_html = local_testRow(row)
    nback = row.n_back(1);
    param = char(row.parameter(1));
    family = char(row.test_family(1));
    support = 'No';
    support_class = 'no';
    if logical(row.supports_hypothesis(1))
        support = 'Yes';
        support_class = 'yes';
    end

    row_html = sprintf(['<tr data-nback="%d" data-param="%s" data-family="%s">' ...
        '<td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' ...
        '<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td>' ...
        '<td>%s</td><td>%s</td><td>%s</td><td>%s</td>' ...
        '<td><span class="pill %s">%s</span></td></tr>'], ...
        nback, local_escape(param), local_escape(family), ...
        nback, local_escape(char(row.hypothesis(1))), ...
        local_escape(char(row.parameter_label(1))), local_escape(local_axisLabel(char(row.axis(1)))), ...
        local_escape(local_familyLabel(family)), ...
        local_escape(char(row.comparison(1))), ...
        local_num(row.estimate(1)), ...
        local_num(local_getNumericScalar(row, 'boot_median')), ...
        local_ci(row.bca_lo(1), row.bca_hi(1)), ...
        local_ci(row.pc_lo(1), row.pc_hi(1)), ...
        local_escape(local_admitLabel(row)), ...
        local_escape(local_diagnosticLabel(row)), ...
        local_escape(char(row.p_bca_label(1))), ...
        local_escape(char(row.p_holm_label(1))), ...
        support_class, support);
end

function local_writeLinkBlock(fid, n_back)
    base = sprintf('%d_back/super_subject', n_back);
    fprintf(fid, '<article class="links"><h3>%d-back</h3>', n_back);
    local_link(fid, base, 'Targeted hypothesis tests CSV', 'targeted_dog_hypothesis_tests.csv');
    local_link(fid, base, 'Amplitude endpoint effects', 'targeted_dog_endpoint_effects_amplitude.pdf');
    local_link(fid, base, 'FWHM endpoint effects', 'targeted_dog_endpoint_effects_fwhm.pdf');
    local_link(fid, base, 'Amplitude scatter', 'unbinned_mle_super_sd_amplitude_scatter.pdf');
    local_link(fid, base, 'Range/FWHM scatter', 'unbinned_mle_super_sd_width_fwhm_scatter.pdf');
    local_link(fid, base, 'Amplitude pooled + subject points', 'unbinned_mle_super_sd_amplitude_pooled_subject_points.pdf');
    local_link(fid, base, 'FWHM pooled + subject points', 'unbinned_mle_super_sd_width_fwhm_pooled_subject_points.pdf');
    local_link(fid, base, 'DoG grid: contrast', 'unbinned_mle_isolated_dog_contrast.pdf');
    local_link(fid, base, 'DoG grid: precision', 'unbinned_mle_isolated_dog_precision.pdf');
    local_link(fid, base, 'All contrasts', 'contrasts.csv');
    local_link(fid, base, 'Close/far sigma delta tests', 'close_far_sigma_delta_tests.csv');
    fprintf(fid, '</article>');
end

function local_writeAcrossNBackLinkBlock(fid)
    base = 'across_n_back';
    fprintf(fid, '<article class="links"><h3>Across n-back</h3>');
    local_link(fid, base, 'Amplitude endpoint effects by n-back', 'targeted_dog_endpoint_effects_by_nback_amplitude.pdf');
    local_link(fid, base, 'FWHM endpoint effects by n-back', 'targeted_dog_endpoint_effects_by_nback_fwhm.pdf');
    local_link(fid, base, 'Amplitude by n-back, past level', 'unbinned_mle_sd_amplitude_by_nback_by_past_level.pdf');
    local_link(fid, base, 'FWHM by n-back, contrast', 'unbinned_mle_super_sd_width_fwhm_scatter_by_nback_contrast.pdf');
    fprintf(fid, '</article>');
end

function local_link(fid, base, label, fname)
    href = local_href([base '/' fname]);
    fprintf(fid, '<a href="%s">%s</a>', href, local_escape(label));
end

function s = local_css()
    s = [ ...
        "body{margin:0;background:#f7f7f5;color:#101010;font-family:Helvetica,Arial,sans-serif;}" newline ...
        ".page{max-width:1380px;margin:0 auto;padding:34px;}" newline ...
        ".hero{display:flex;justify-content:space-between;gap:24px;align-items:flex-start;margin-bottom:22px;}" newline ...
        "h1{margin:0 0 10px;font-size:38px;letter-spacing:-.03em;} h2{margin:0 0 14px;font-size:24px;} h3{margin:0 0 10px;font-size:18px;}" newline ...
        "p{font-size:18px;line-height:1.35;color:#444;margin:0;max-width:920px;}" newline ...
        ".button,.links a{display:inline-block;border:1.5px solid #101010;background:white;color:#101010;text-decoration:none;padding:9px 12px;margin:4px 6px 4px 0;font-size:14px;}" newline ...
        ".controls{display:flex;gap:16px;flex-wrap:wrap;background:white;border:1.5px solid #c7cbd1;padding:14px 16px;margin-bottom:18px;}" newline ...
        "label{font-weight:700;color:#6f3c00;} select{margin-left:8px;font-size:15px;padding:4px 8px;}" newline ...
        ".grid{display:grid;grid-template-columns:1fr 1fr;gap:18px;margin-bottom:18px;}" newline ...
        ".card{background:white;border:1.5px solid #c7cbd1;padding:20px;} .card.full{margin-bottom:18px;}" newline ...
        ".section-copy{margin-bottom:14px;}" newline ...
        ".interpretation-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:14px;margin:16px 0;}" newline ...
        ".takeaway{border:1.5px solid #c7cbd1;background:#fafafa;padding:16px;} .takeaway.yes{border-color:#5c9c5c;background:#f2f8f2;} .takeaway.warn{border-color:#b88a33;background:#fff8ea;} .takeaway.neutral{background:#f7f7f5;}" newline ...
        ".tag{display:inline-block;margin-bottom:8px;text-transform:uppercase;letter-spacing:.06em;font-size:11px;font-weight:700;color:#6f3c00;}" newline ...
        ".statline{margin-top:10px;font-size:14px;color:#555;}" newline ...
        ".interpretation-notes{margin-top:14px;border-top:1px solid #d9dce0;padding-top:12px;}" newline ...
        "li{font-size:16px;line-height:1.35;margin:7px 0;color:#333;}" newline ...
        ".table-wrap{overflow:auto;} table{width:100%;border-collapse:collapse;font-size:14px;} th,td{border-bottom:1px solid #d9dce0;text-align:left;vertical-align:top;padding:9px 8px;}" newline ...
        "th{position:sticky;top:0;background:#fff;color:#6f3c00;font-size:13px;text-transform:uppercase;letter-spacing:.04em;}" newline ...
        "code{font-family:Menlo,Consolas,monospace;background:#f1f1ef;padding:1px 4px;}" newline ...
        ".pill{display:inline-block;border-radius:999px;padding:3px 9px;font-weight:700;} .pill.yes{background:#e5f3e5;color:#146b14;} .pill.no{background:#f1f1ef;color:#555;} .pill.warn{background:#fff0c7;color:#704700;}" newline ...
        ".link-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(290px,1fr));gap:14px;} .links{border:1px solid #c7cbd1;padding:14px;background:#fafafa;}" newline ...
        "tr.hidden{display:none;}" ...
        ];
end

function s = local_js()
    s = [ ...
        "const filters=['nbackFilter','paramFilter','familyFilter'];" newline ...
        "function applyFilters(){const n=document.getElementById('nbackFilter').value;const p=document.getElementById('paramFilter').value;const f=document.getElementById('familyFilter').value;document.querySelectorAll('#testsTable tbody tr').forEach(r=>{const show=(n==='all'||r.dataset.nback===n)&&(p==='all'||r.dataset.param===p)&&(f==='all'||r.dataset.family===f);r.classList.toggle('hidden',!show);});}" newline ...
        "filters.forEach(id=>document.getElementById(id).addEventListener('change',applyFilters));" newline ...
        "applyFilters();" ...
        ];
end

function s = local_num(x)
    if ~isfinite(x)
        s = 'NA';
    else
        s = sprintf('%.3f', x);
    end
end

function s = local_ci(lo, hi)
    if ~isfinite(lo) || ~isfinite(hi)
        s = 'NA';
    else
        s = sprintf('[%.3f, %.3f]', lo, hi);
    end
end

function x = local_getNumericScalar(row, var_name)
    if ismember(var_name, row.Properties.VariableNames)
        x = row.(var_name)(1);
    else
        x = NaN;
    end
end

function tf = local_getLogicalColumn(tbl, var_name, default_value)
    if ismember(var_name, tbl.Properties.VariableNames)
        tf = logical(tbl.(var_name));
    else
        tf = repmat(default_value, height(tbl), 1);
    end
end

function s = local_admitLabel(row)
    n_admit = local_getNumericScalar(row, 'n_admit');
    admit_fraction = local_getNumericScalar(row, 'admit_fraction');
    if isfinite(n_admit) && isfinite(admit_fraction)
        s = sprintf('%d (%.1f%%)', round(n_admit), 100 * admit_fraction);
    elseif isfinite(n_admit)
        s = sprintf('%d', round(n_admit));
    else
        s = 'NA';
    end
end

function s = local_diagnosticLabel(row)
    if ismember('valid_for_inference', row.Properties.VariableNames) && logical(row.valid_for_inference(1))
        s = 'OK';
        return
    end

    flags = {};
    if local_getLogicalScalar(row, 'estimate_boot_mismatch')
        flags{end+1} = 'estimate outside percentile CI'; %#ok<AGROW>
    end
    if local_getLogicalScalar(row, 'bca_pc_conflict')
        flags{end+1} = 'active/percentile conflict'; %#ok<AGROW>
    end
    if local_getLogicalScalar(row, 'extreme_z0')
        flags{end+1} = 'extreme z0'; %#ok<AGROW>
    end
    if local_getLogicalScalar(row, 'low_admit_fraction')
        flags{end+1} = 'low admitted fraction'; %#ok<AGROW>
    end
    if isempty(flags)
        flags{end+1} = 'insufficient diagnostics';
    end
    s = strjoin(flags, '; ');
end

function tf = local_getLogicalScalar(row, var_name)
    if ismember(var_name, row.Properties.VariableNames)
        tf = logical(row.(var_name)(1));
    else
        tf = false;
    end
end

function s = local_axisLabel(axis_name)
    switch axis_name
        case 'average_prev_curr'
            s = 'Average of previous and current levels';
        case 'previous_level'
            s = 'Previous level';
        case 'current_level'
            s = 'Current level';
        case 'boundary_interaction'
            s = 'Boundary interaction';
        otherwise
            s = axis_name;
    end
end

function s = local_familyLabel(family_name)
    switch family_name
        case 'within_manipulation'
            s = 'Within manipulation';
        case 'between_manipulation'
            s = 'Between manipulations';
        otherwise
            s = local_titleCase(strrep(family_name, '_', ' '));
    end
end

function s = local_href(path_str)
    s = strrep(path_str, filesep, '/');
    s = strrep(s, ' ', '%20');
end

function s = local_escape(s)
    s = char(string(s));
    s = strrep(s, '&', '&amp;');
    s = strrep(s, '<', '&lt;');
    s = strrep(s, '>', '&gt;');
    s = strrep(s, '"', '&quot;');
end

function s = local_titleCase(s)
    if isempty(s)
        return
    end
    s(1) = upper(s(1));
end
